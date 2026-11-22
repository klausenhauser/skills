"""Helpers for surveying PDB crystallization precedent for a protein domain
and cross-referencing conditions against commercial screen formulations.
See SKILL.md for the full workflow this backs.

Python 3 standard library only. Import it from the skill's scripts folder:

    import sys; sys.path.insert(0, "<skill>/scripts"); import survey

Network access needed: rest.uniprot.org, data.rcsb.org, files.rcsb.org,
c6.csiro.au (screen formulations).
"""

import re


def fetch_uniprot_entry(accession):
    """Fetch a UniProt entry and return sequence, domain features, and PDB cross-refs.

    Returns dict: {accession, sequence, domains: [{type, start, end, description}],
                   pdb_xrefs: [{pdb_id, method, resolution, chains: [{ids, start, end}]}]}
    """
    import json
    import urllib.request

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    with urllib.request.urlopen(url, timeout=20) as r:
        d = json.loads(r.read())

    seq = d.get("sequence", {}).get("value", "")
    domains = []
    for f in d.get("features", []):
        if f.get("type") in ("Domain", "Region", "Repeat"):
            loc = f.get("location", {})
            domains.append(dict(
                type=f["type"],
                start=loc.get("start", {}).get("value"),
                end=loc.get("end", {}).get("value"),
                description=f.get("description", ""),
            ))

    pdb_xrefs = []
    for xr in d.get("uniProtKBCrossReferences", []):
        if xr.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in xr.get("properties", [])}
        chains = []
        for chain_spec in props.get("Chains", "").split(","):
            chain_spec = chain_spec.strip()
            m = re.match(r"([A-Za-z/]+)=(\d+)-(\d+)", chain_spec)
            if m:
                chains.append(dict(ids=m.group(1), start=int(m.group(2)), end=int(m.group(3))))
        pdb_xrefs.append(dict(
            pdb_id=xr.get("id"),
            method=props.get("Method"),
            resolution=props.get("Resolution"),
            chains=chains,
        ))

    return dict(accession=accession, sequence=seq, domains=domains, pdb_xrefs=pdb_xrefs)


def select_pdb_ids_for_domain(pdb_xrefs, domain_start, domain_end, min_overlap_frac=0.5):
    """Filter a pdb_xrefs list (from fetch_uniprot_entry) to entries whose chain
    range overlaps the given domain by at least min_overlap_frac of the domain length.
    Returns a sorted list of PDB ID strings.
    """
    domain_len = domain_end - domain_start + 1
    hits = []
    for entry in pdb_xrefs:
        for chain in entry.get("chains", []):
            overlap = min(chain["end"], domain_end) - max(chain["start"], domain_start) + 1
            if overlap > 0 and overlap / domain_len >= min_overlap_frac:
                hits.append(entry["pdb_id"])
                break
    return sorted(set(hits))


def canonical_range(entity_seq, canonical_seq, kmer=15):
    """Map a deposited entity sequence onto the canonical UniProt sequence.

    Depositors often number against another isoform while labelling the canonical
    accession, so the deposited residue numbers can't be trusted. This finds the first
    and last `kmer`-long stretches of the entity that occur in the canonical sequence.

    Returns dict: {start, end (1-based canonical numbering, None if no match),
                   n_extra (residues outside the matched stretch), n_term_extra,
                   c_term_extra (the extra sequences themselves)}
    """
    seq = entity_seq or ""
    hits = [i for i in range(len(seq) - kmer + 1) if seq[i:i + kmer] in canonical_seq]
    if not hits:
        return dict(start=None, end=None, n_extra=len(seq), n_term_extra="", c_term_extra="")
    a, b = hits[0], hits[-1] + kmer
    return dict(
        start=canonical_seq.find(seq[a:a + kmer]) + 1,
        end=canonical_seq.find(seq[b - kmer:b]) + kmer,
        n_extra=a + len(seq) - b,
        n_term_extra=seq[:a],
        c_term_extra=seq[b:],
    )


def tag_status(canonical):
    """Classify one entry's tag from canonical_range's result, end by end:
    'intact tag' if either end carries 8+ extra residues (His tag, uncleaved protease
    site...), 'scar' if an end carries 1-7 (a cleaved tag's remnant: GP, GPHM, SNA...),
    else 'none'. Confirm a family's call with get_construct_tag_diffs before quoting it.
    """
    longest = max(len(canonical.get("n_term_extra") or ""), len(canonical.get("c_term_extra") or ""))
    return "intact tag" if longest >= 8 else ("scar" if longest >= 1 else "none")


def construct_families(ranges, tol=3):
    """Cluster construct ranges that differ by at most `tol` residues at each end
    (1-313, 1-312 and 2-313 are one construct with or without the initiator Met).

    ranges: {pdb_id: (start, end)} in canonical numbering (from canonical_range).
    Returns families largest first: [{label: 'start–end' of the most common range,
    start, end, members: [pdb_ids], ranges: {(start, end): count}}]. Entries with no
    range (start None) form a family labelled 'unmapped'.
    """
    from collections import Counter

    counts = Counter(r for r in ranges.values() if r[0] is not None)
    fams = []
    for (s, e), _ in counts.most_common():
        for f in fams:
            if abs(f["start"] - s) <= tol and abs(f["end"] - e) <= tol:
                f["ranges"][(s, e)] = counts[(s, e)]
                break
        else:
            fams.append(dict(label=f"{s}–{e}", start=s, end=e, ranges={(s, e): counts[(s, e)]}))
    for f in fams:
        f["members"] = sorted(p for p, r in ranges.items() if r in f["ranges"])
    unmapped = sorted(p for p, r in ranges.items() if r[0] is None)
    if unmapped:
        fams.append(dict(label="unmapped", start=None, end=None, ranges={}, members=unmapped))
    return sorted(fams, key=lambda f: -len(f["members"]))


def fetch_rcsb_crystallization_batch(pdb_ids, sleep_s=0.05):
    """Batch-fetch crystallization conditions, resolution, and primary citation
    for a list of PDB IDs via the RCSB REST API (one polymer_entity/1 + entry call per ID).

    Returns dict: {pdb_id: {resolution, ph, temp_k, pdbx_details, citation_title,
                             citation_year, entity_seq}}
    """
    import json
    import time
    import urllib.request

    out = {}
    for pid in pdb_ids:
        entry = {}
        try:
            with urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=20) as r:
                d = json.loads(r.read())
            grow = (d.get("exptl_crystal_grow") or [{}])[0]
            cit = (d.get("rcsb_primary_citation") or {})
            entry["resolution"] = (d.get("rcsb_entry_info", {}) or {}).get("resolution_combined", [None])[0]
            entry["ph"] = grow.get("pH")
            entry["temp_k"] = grow.get("temp")
            entry["pdbx_details"] = grow.get("pdbx_details")
            entry["citation_title"] = cit.get("title")
            entry["citation_year"] = cit.get("year")
        except Exception as e:
            entry["error"] = str(e)
        try:
            with urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/1", timeout=20) as r:
                ed = json.loads(r.read())
            entry["entity_seq"] = ed.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        except Exception as e:
            entry["entity_seq_error"] = str(e)
        out[pid] = entry
        time.sleep(sleep_s)
    return out


def group_conditions(batch):
    """Group PDB IDs from fetch_rcsb_crystallization_batch by crystallization condition,
    ignoring case, whitespace and trailing punctuation, so 'Na Acetate, pH 6.5' and
    'NA ACETATE, pH 6.5' land together. Still reads the strings by eye afterwards:
    the same chemistry written differently stays in separate groups.

    Returns a list of {condition, pdb_ids}, largest group first. Entries with no
    recorded details are grouped under condition ''.
    """
    groups = {}
    for pid, v in batch.items():
        raw = " ".join((v.get("pdbx_details") or "").split())
        key = raw.lower().rstrip(" .,;")
        g = groups.setdefault(key, dict(condition=raw, pdb_ids=[]))
        g["pdb_ids"].append(pid)
    return sorted(groups.values(), key=lambda g: (-len(g["pdb_ids"]), g["condition"].lower()))


def get_construct_tag_diffs(pdb_id):
    """Report, per chain, how the deposited crystal sequence differs from its UniProt
    reference: the aligned accession and residue ranges, plus every residue flagged in
    `_struct_ref_seq_dif` (expression tag, conflict, engineered mutation).

    Parses the mmCIF directly because the RCSB JSON APIs do not expose
    `_struct_ref_seq_dif`. Handles BOTH mmCIF encodings: the `loop_` form used when an
    entry has several aligned chains, and the flat `_tag  value` form used when it has
    only one -- single-chain entries are the common case, so a parser that reads only
    the loop form silently returns [] for most entries.

    CAUTION on numbering: the residue numbers returned are as deposited. Depositors
    sometimes number against a splice isoform while still labelling the canonical
    accession, so a contiguous-looking db range here is NOT evidence that the construct
    is a contiguous stretch of the canonical sequence. Align the deposited entity
    sequence against the canonical sequence AND the annotated isoforms before quoting
    construct boundaries.
    """
    import re
    import shlex as _shlex
    import urllib.request

    url = "https://files.rcsb.org/download/{}.cif".format(pdb_id)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8", errors="ignore")

    def loop_rows(prefix):
        """Return (column_names, rows) for a `loop_` block, else (None, [])."""
        marker = "loop_"
        idx = 0
        while True:
            idx = text.find(marker, idx)
            if idx == -1:
                return None, []
            block = text[idx + len(marker):]
            lines = block.splitlines()
            cols, consumed = [], 0
            for line in lines:
                s = line.strip()
                consumed += 1
                if s.startswith(prefix):
                    cols.append(s.split(".", 1)[1].strip())
                elif cols:
                    break
                elif s and not s.startswith("_"):
                    break
                elif s.startswith("_"):
                    cols = []
                    break
            if not cols:
                idx += len(marker)
                continue
            rows = []
            for line in lines[consumed - 1:]:
                s = line.strip()
                if s.startswith("#") or s.startswith("loop_") or not s:
                    break
                try:
                    rows.append(_shlex.split(s))
                except ValueError:
                    rows.append(s.split())
            return cols, rows

    def flat_values(prefix):
        """Return a dict for the flat `_tag  value` encoding, else {}."""
        out = {}
        for line in text.splitlines():
            s = line.strip()
            if s.startswith(prefix):
                parts = s.split(None, 1)
                if len(parts) == 2:
                    out[parts[0].split(".", 1)[1]] = parts[1].strip()
        return out

    seq_cols, seq_rows = loop_rows("_struct_ref_seq.")
    if seq_cols:
        aligns = [dict(zip(seq_cols, r)) for r in seq_rows if len(r) == len(seq_cols)]
    else:
        flat = flat_values("_struct_ref_seq.")
        aligns = [flat] if flat else []

    dif_cols, dif_rows = loop_rows("_struct_ref_seq_dif.")
    if dif_cols:
        difs = [dict(zip(dif_cols, r)) for r in dif_rows if len(r) == len(dif_cols)]
    else:
        flat = flat_values("_struct_ref_seq_dif.")
        difs = [flat] if flat else []

    results = []
    for a in aligns:
        align_id = a.get("align_id")
        chain_difs = []
        for d in difs:
            if align_id is not None and d.get("align_id") not in (align_id, None):
                continue
            detail = d.get("details", "?")
            if detail in ("?", ".", None):
                continue
            detail = detail.strip("'\"")
            chain_difs.append({
                "resname": d.get("mon_id"),
                "auth_seq_num": d.get("pdbx_auth_seq_num"),
                "detail": detail,
            })
        results.append({
            "strand_id": a.get("pdbx_strand_id"),
            "ref_db_accession": a.get("pdbx_db_accession"),
            "db_align_beg": a.get("db_align_beg"),
            "db_align_end": a.get("db_align_end"),
            "auth_seq_align_beg": a.get("pdbx_auth_seq_align_beg"),
            "auth_seq_align_end": a.get("pdbx_auth_seq_align_end"),
            "diffs": chain_difs,
        })
    return results


def fetch_c6_screen(screen_id, screen_name=None):
    """Fetch a screen's formulation from CSIRO's C6 database as structured rows.

    Find the screen's numeric ID by searching its name on https://c6.csiro.au/screens
    (the ID is in the screen page's URL). C6 covers most commercial screens
    (Molecular Dimensions, Hampton, Jena...) and is preferred over parsing vendor
    PDFs: every well comes with ingredient, role, concentration, units and pH.

    Returns a list of dicts: {screen, well, raw, ph, components}, wells in plate order
    (A1, A2 ... H12 for 96-well screens). `raw` is a readable one-line formulation,
    the same shape parse_screen_pdf_text returns, so both feed
    search_screens_for_terms.
    """
    import json
    import urllib.request

    url = f"https://c6.csiro.au/export_screen/JSON-BCC?screen={screen_id}"
    with urllib.request.urlopen(url, timeout=30) as r:
        d = json.loads(r.read())
    name = screen_name or d.get("name") or f"C6 {screen_id}"
    wells = d.get("wells", [])
    cols = 12 if len(wells) > 48 else 6
    rows = []
    for k, comps in enumerate(wells):
        parts, ph = [], None
        for c in comps:
            conc = c.get("concentration")
            conc = f"{conc:g}" if isinstance(conc, (int, float)) else (conc or "")
            units = c.get("units") or ""
            units = f"% {units}" if units in ("v/v", "w/v", "w/w") else units
            part = f"{conc} {units} {c.get('ingredient', '')}".strip()
            if c.get("pH") is not None:
                part += f" pH {c['pH']:g}"
                ph = c["pH"]
            parts.append(part)
        rows.append(dict(screen=name, well=f"{'ABCDEFGH'[k // cols]}{k % cols + 1}",
                         raw="; ".join(parts), ph=ph, components=comps))
    return rows


def parse_screen_pdf_text(text, screen_name):
    """Fallback when a screen isn't in C6: parse text extracted from a formulation PDF
    into structured rows. Handles four layouts seen in practice:
      (a) leading condition number + well ID (JCSG+, JCSG Core I):
          '<num> <Well> <condition text>'
      (b) leading well ID only, no number (PACT premier):
          '<Well> <condition text>'
      (c) leading bare condition number, no letter-number well ID (Structure Screen 2 /
          classic sparse-matrix screens): '<num> <condition text starting with a
          concentration, e.g. "0.1 M ...">' -- well is reported as '#<num>'.
      (d) box-tube number (JCSG+ formulation sheet): '<box>-<tube> <condition text>'
          -- well is reported as '<box>-<tube>'.
    Returns a list of dicts: {screen, well, raw, ph}; ph is None because PDF tables
    don't label it reliably. Rows that don't match any pattern are dropped silently --
    inspect `text` directly for a new screen layout.
    """
    import re

    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        m = re.match(r"^(\d+)\s+([A-H]\d{1,2})\s+(.+)$", line)
        if m:
            rows.append(dict(screen=screen_name, well=m.group(2), raw=m.group(3), ph=None))
            continue
        m = re.match(r"^([A-H]\d{1,2})\s+(.+)$", line)
        if m:
            rows.append(dict(screen=screen_name, well=m.group(1), raw=m.group(2), ph=None))
            continue
        m = re.match(r"^(\d-\d{1,2})\s+(.+)$", line)
        if m:
            rows.append(dict(screen=screen_name, well=m.group(1), raw=m.group(2), ph=None))
            continue
        m = re.match(r"^(\d{1,2})\s+([\d.]+\s*M\s+.+)$", line)
        if m:
            rows.append(dict(screen=screen_name, well=f"#{m.group(1)}", raw=m.group(2), ph=None))
    return rows


def normalise_chemistry(text):
    """Lower-case a formulation or search term and spell PEGs one way, so
    'Polyethylene glycol 3,350', 'PEG3350', 'PEG 3.35K' and 'peg 3350' all read
    'peg 3350' and 'PEG MME 5,000' reads 'peg mme 5000'."""
    import re

    s = text.lower()
    s = re.sub(r"(?<=\d),(?=\d{3})", "", s)
    s = re.sub(r"poly\s?ethylene\s+glycol", "peg", s)
    s = re.sub(r"\bpeg\s*(monomethyl ether|mme)\s*", "peg mme ", s)
    s = re.sub(r"\bpeg\s*(mme\s+)?-?\s*(\d+(?:\.\d+)?)\s*k\b",
               lambda m: f"peg {m.group(1) or ''}{float(m.group(2)) * 1000:g}", s)
    s = re.sub(r"\bpeg(\d)", r"peg \1", s)
    return " ".join(s.split())


def search_screens_for_terms(screen_rows, terms, mode="all", ph=None, ph_tol=0.5):
    """Search parsed screen rows (from fetch_c6_screen or parse_screen_pdf_text, or
    several screens concatenated) for keyword matches.

    terms: list of case-insensitive words or phrases, e.g. ['sodium citrate', 'peg 3350'].
           Both terms and rows go through normalise_chemistry, and a term must match a
           whole word: 'tris' does not match 'bis-tris', 'mes' does not match 'times'.
    mode: 'all' requires every term present in the row; 'any' requires at least one.
    ph: if given, also require the well's pH within ph_tol. Rows with unknown pH
        (PDF rows) are kept and marked ph_checked=False.
    Returns the matching rows, each with an added 'matched' list of the terms found.
    """
    import re

    pats = [(t, re.compile(r"(?<![\w-])" + re.escape(normalise_chemistry(t)) + r"(?![\w-])"))
            for t in terms]
    out = []
    for row in screen_rows:
        raw = normalise_chemistry(row["raw"])
        matched = [t for t, p in pats if p.search(raw)]
        if not ((mode == "all" and len(matched) == len(pats)) or (mode == "any" and matched)):
            continue
        hit = dict(row, matched=matched, ph_checked=row.get("ph") is not None)
        if ph is not None and row.get("ph") is not None and abs(row["ph"] - ph) > ph_tol:
            continue
        out.append(hit)
    return out


def rank_screen_matches(screen_rows, terms, ph=None, ph_tol=0.5, top=3):
    """Rank screen wells by how much of one literature condition they reproduce.

    terms: the condition's components, e.g. ['lithium chloride', 'peg 6000', 'tris'].
    Each matched term scores 1 and a pH within ph_tol scores 1 more, so a well with
    every component at the right pH scores len(terms) + 1. Wells matching no term
    are dropped.
    Returns up to `top` rows, best first, each with 'score', 'max_score', 'matched',
    'missing' and 'ph_ok' (None when the well's pH is unknown).
    """
    out = []
    for row in search_screens_for_terms(screen_rows, terms, mode="any"):
        ph_ok = None
        if ph is not None and row.get("ph") is not None:
            ph_ok = abs(row["ph"] - ph) <= ph_tol
        score = len(row["matched"]) + (1 if ph_ok else 0)
        out.append(dict(row, score=score, max_score=len(terms) + (1 if ph is not None else 0),
                        missing=[t for t in terms if t not in row["matched"]], ph_ok=ph_ok))
    out.sort(key=lambda r: (-r["score"], r["well"]))
    return out[:top]


GROUP_HUES = ["#4e79a7", "#f28e2b", "#59a14f", "#b07aa1", "#76b7b2", "#edc948", "#ff9da7", "#9c755f"]
OTHER_HUE = "#8c8c8c"
STATUS_PILLS = {  # first word of a pill -> class; anything else is neutral
    "intact": "warn", "scar": "info", "cleaved": "info", "none": "muted",
    "exact": "good", "partial": "warn", "no": "bad", "✓": "good", "✗": "bad", "?": "muted",
}


def write_html_report(path, title, sections, sources=()):
    """Write the survey as one self-contained HTML page (no external files) that opens
    in any browser, follows light/dark mode and prints cleanly to PDF.

    Each section is a dict with `heading` and any of:
        text:    paragraphs; blank lines separate them, '- ' lines become bullets,
                 **bold** is kept
        table:   a flat list of row dicts
        groups:  [{key, label, note, rows}] -- one block per group: a tinted header
                 band (label, note) over its own table. Use it for condition groups
                 and screen hits per condition.
        pills:   column names rendered as status pills (tag status, match quality,
                 ✓/✗). A cell holding a list renders one pill per item.
        swatch:  column name that gets a coloured dot from each row's `_swatch` key
                 (e.g. its construct family)
        collapse: wrap each group in a <details> block, open if the group has more
                  than one row; all open when printed
    Keys starting with '_' are never shown as columns.

    Colours: every group `key` and row `_swatch` is a string like 'fam:29–313' or
    'cond:3'. Within one prefix, keys get hues in order of first appearance, the same
    hue everywhere in the report, so a condition's colour matches between the screen
    hits and the PDB list. After 8 hues, and for keys ending ':other', it is grey.
    sources: optional list of (label, url) pairs listed at the end.
    """
    import datetime
    import html

    hues = {}

    def hue(key):
        if key is None:
            return OTHER_HUE
        if key not in hues:
            prefix = key.split(":", 1)[0]
            used = sum(1 for k in hues if k.split(":", 1)[0] == prefix and hues[k] != OTHER_HUE)
            hues[key] = OTHER_HUE if key.endswith(":other") or used >= len(GROUP_HUES) else GROUP_HUES[used]
        return hues[key]

    def inline(s):
        return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", html.escape(s))

    def text_block(text):
        out = []
        for para in text.strip().split("\n\n"):
            lines = [l.strip() for l in para.strip().splitlines() if l.strip()]
            if lines and all(l.startswith("- ") for l in lines):
                out.append("<ul>" + "".join(f"<li>{inline(l[2:])}</li>" for l in lines) + "</ul>")
            elif lines:
                out.append(f"<p>{inline(' '.join(lines))}</p>")
        return "\n".join(out)

    def pill(v):
        word = str(v).split()[0].lower() if str(v).strip() else ""
        return f'<span class="pill {STATUS_PILLS.get(word, "muted")}">{html.escape(str(v))}</span>'

    def cell(row, col, s):
        v = row.get(col)
        if col in s.get("pills", ()):
            body = " ".join(pill(x) for x in (v if isinstance(v, (list, tuple)) else [v]) if x not in (None, ""))
        else:
            body = html.escape("" if v is None else ", ".join(map(str, v)) if isinstance(v, (list, tuple)) else str(v))
        if col == s.get("swatch") and row.get("_swatch"):
            body = f'<span class="dot" style="--g:{hue(row["_swatch"])}"></span>{body}'
        return f"<td>{body}</td>"

    def table(rows, s):
        cols = [c for c in rows[0] if not c.startswith("_")]
        head = "".join(f"<th>{html.escape(c.replace('_', ' '))}</th>" for c in cols)
        body = "".join("<tr>" + "".join(cell(r, c, s) for c in cols) + "</tr>" for r in rows)
        return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'

    def group(g, s):
        band = (f'<span class="dot" style="--g:{hue(g.get("key"))}"></span><strong>{inline(g["label"])}</strong>'
                + (f' <span class="note">{inline(g["note"])}</span>' if g.get("note") else ""))
        inner = table(g["rows"], s) if g.get("rows") else ""
        style = f'style="--g:{hue(g.get("key"))}"'
        if s.get("collapse"):
            opened = " open" if len(g.get("rows") or []) > 1 else ""
            return f'<details class="grp" {style}{opened}><summary class="band">{band}</summary>{inner}</details>'
        return f'<div class="grp" {style}><div class="band">{band}</div>{inner}</div>'

    parts = []
    for s in sections:
        parts.append(f"<h2>{html.escape(s['heading'])}</h2>")
        if s.get("text"):
            parts.append(text_block(s["text"]))
        if s.get("table"):
            parts.append(table(s["table"], s))
        for g in s.get("groups") or []:
            parts.append(group(g, s))
    if sources:
        parts.append("<h2>Sources</h2><ul>" + "".join(
            f'<li><a href="{html.escape(u)}">{html.escape(l)}</a></li>' for l, u in sources) + "</ul>")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{ --bg:#fff; --fg:#1d1d1f; --muted:#6e6e73; --line:#d9d9de; --head:#f4f4f6; --link:#0a58ca;
  --good:#1e7b34; --warn:#a15c00; --bad:#b42318; --info:#1f5fa8; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#16161a; --fg:#ececf1; --muted:#a1a1aa; --line:#33333b; --head:#222228;
  --link:#7cb3ff; --good:#6fcf86; --warn:#f0b35a; --bad:#ff8a80; --info:#8cbcff; }} }}
* {{ -webkit-print-color-adjust:exact; print-color-adjust:exact; }}
body {{ background:var(--bg); color:var(--fg); font:15px/1.5 -apple-system, "Segoe UI", Roboto, sans-serif; max-width:1100px; margin:0 auto; padding:24px 16px; }}
h1 {{ font-size:1.6em; margin:0 0 4px; }} h2 {{ font-size:1.2em; margin:32px 0 10px; border-bottom:1px solid var(--line); padding-bottom:4px; }}
.meta, .note {{ color:var(--muted); font-size:.9em; }} a {{ color:var(--link); }}
.scroll {{ overflow-x:auto; }} table {{ border-collapse:collapse; width:100%; font-size:.85em; }}
th, td {{ border:1px solid var(--line); padding:4px 6px; text-align:left; vertical-align:top; }}
th {{ background:var(--head); }}
.grp {{ margin:10px 0; border-left:6px solid var(--g); }}
.grp table {{ border-left:none; }}
.band {{ display:block; padding:6px 10px; background:color-mix(in srgb, var(--g) 18%, var(--bg)); cursor:default; }}
summary.band {{ cursor:pointer; }}
.dot {{ display:inline-block; width:.75em; height:.75em; border-radius:50%; background:var(--g); margin-right:6px; vertical-align:baseline; }}
.pill {{ display:inline-block; border:1px solid currentColor; border-radius:9px; padding:0 7px; margin:1px 2px 1px 0; font-size:.85em; white-space:nowrap; }}
.pill.good {{ color:var(--good); }} .pill.warn {{ color:var(--warn); }} .pill.bad {{ color:var(--bad); }}
.pill.info {{ color:var(--info); }} .pill.muted {{ color:var(--muted); }}
@media print {{ body {{ max-width:none; padding:0; font-size:11px; }} .scroll {{ overflow:visible; }} tr, .band {{ break-inside:avoid; }} }}
</style>
<script>addEventListener("beforeprint", () => document.querySelectorAll("details").forEach(d => d.open = true));</script>
</head><body>
<h1>{html.escape(title)}</h1>
<p class="meta">Crystallization precedent survey · generated {datetime.date.today().isoformat()}</p>
{chr(10).join(parts)}
</body></html>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path
