"""Helpers for deciding which host to express a protein in (E. coli / insect /
mammalian / yeast / cell-free) from UniProt, PDB and PubMed evidence.
See SKILL.md for the full workflow this backs.

Python 3 standard library only. Import it from the skill's scripts folder:

    import sys; sys.path.insert(0, "<skill>/scripts"); import expression

Network access needed: rest.uniprot.org, search.rcsb.org, data.rcsb.org,
eutils.ncbi.nlm.nih.gov, api.crossref.org.
"""

import html
import json
import re
import time
import urllib.parse
import urllib.request

TIMEOUT = 30


def _get_json(url, data=None, headers=None):
    req = urllib.request.Request(url, data=data, headers=headers or {})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return json.loads(r.read())


def _post_json(url, payload):
    return _get_json(url, data=json.dumps(payload).encode(),
                     headers={"Content-Type": "application/json"})


# ------------------------------------------------------------------ sequence

def fetch_uniprot_entry(accession):
    """Fetch a UniProt entry: sequence, the features that decide the host, and PDB
    cross-references.

    Returns dict: {accession, name, organism, sequence, length,
                   features: [{type, start, end, description}],
                   pdb_ids: [...], isoform_count}

    The features kept are the ones that drive the host choice: Signal and Propeptide
    (the mature-chain boundary), Chain, Domain, Disulfide bond, Glycosylation and
    Alternative sequence (splice isoforms -- see the numbering warning on
    canonical_range).
    """
    d = _get_json(f"https://rest.uniprot.org/uniprotkb/{accession}.json")
    keep = ("Signal", "Propeptide", "Transit peptide", "Chain", "Domain",
            "Glycosylation", "Disulfide bond", "Alternative sequence", "Region")
    feats = []
    for f in d.get("features", []):
        if f.get("type") in keep:
            loc = f.get("location", {})
            feats.append(dict(type=f["type"],
                              start=loc.get("start", {}).get("value"),
                              end=loc.get("end", {}).get("value"),
                              description=f.get("description", "")))
    pdb_ids = sorted(x["id"] for x in d.get("uniProtKBCrossReferences", [])
                     if x.get("database") == "PDB")
    names = d.get("proteinDescription", {}).get("recommendedName", {})
    return dict(
        accession=accession,
        name=names.get("fullName", {}).get("value", ""),
        organism=d.get("organism", {}).get("scientificName", ""),
        sequence=d.get("sequence", {}).get("value", ""),
        length=d.get("sequence", {}).get("length"),
        features=feats,
        pdb_ids=pdb_ids,
        isoform_count=sum(1 for c in d.get("comments", [])
                          if c.get("commentType") == "ALTERNATIVE PRODUCTS"
                          for _ in c.get("isoforms", [])),
    )


def fetch_isoform_sequences(accession, max_isoforms=12):
    """{isoform_id: sequence} for every splice isoform UniProt lists, canonical first.

    Feed these to tag_status as `isoform_seqs` whenever the entry has more than one
    isoform, or a construct solved from a longer isoform is miscalled as carrying a
    large N-terminal tag.
    """
    out = {}
    for n in range(1, max_isoforms + 1):
        iso = f"{accession}-{n}"
        try:
            req = urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{iso}.fasta")
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                txt = r.read().decode()
        except Exception:
            break
        seq = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
        if not seq:
            break
        out[iso] = seq
    return out


# Average residue masses (Da); water added once for the free chain.
_AA_MASS = {"A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
            "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
            "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
            "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326}
# Expasy/ProtParam pKa set: reproduces ProtParam's pI to ~0.01 unit on test chains.
_PKA = {"Cterm": 3.55, "Nterm": 7.50, "D": 4.05, "E": 4.45, "C": 9.00,
        "Y": 10.00, "H": 5.98, "K": 10.00, "R": 12.00}


def construct_props(seq):
    """Length / mass / pI / Cys count / epsilon-280 for the exact chain you will order
    (tag + linker + insert), not the database sequence.

    pI is a Henderson-Hasselbalch estimate on the Expasy/ProtParam pKa set, which it
    reproduces to about 0.01 unit on test chains -- enough to pick an IEX buffer, not
    a substitute for measurement.
    ext280_oxidised adds 125 M-1 cm-1 per cystine, assuming every Cys is paired.
    n_glyc_sequons lists N-X-S/T positions; a sequon is not an occupied site.
    """
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", (seq or "").upper())
    if not s:
        return dict(length=0)

    def charge(ph):
        q = 1 / (1 + 10 ** (ph - _PKA["Nterm"])) - 1 / (1 + 10 ** (_PKA["Cterm"] - ph))
        for aa in "KRH":
            q += s.count(aa) / (1 + 10 ** (ph - _PKA[aa]))
        for aa in "DECY":
            q -= s.count(aa) / (1 + 10 ** (_PKA[aa] - ph))
        return q

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    ext_red = 5500 * s.count("W") + 1490 * s.count("Y")
    return dict(length=len(s),
                mw_kda=round((sum(_AA_MASS[a] for a in s) + 18.0153) / 1000, 2),
                pI=round((lo + hi) / 2, 2),
                cys=s.count("C"),
                ext280_reduced=ext_red,
                ext280_oxidised=ext_red + 125 * (s.count("C") // 2),
                n_glyc_sequons=[m.start() + 1 for m in re.finditer(r"(?=N[^P][ST])", s)])


# ---------------------------------------------------------------- structures

def pdb_entries_for_uniprot(accession):
    """PDB IDs whose polymer entities reference this UniProt accession."""
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": ("rcsb_polymer_entity_container_identifiers"
                          ".reference_sequence_identifiers.database_accession"),
            "operator": "exact_match", "value": accession}},
         "return_type": "entry",
         "request_options": {"return_all_hits": True}}
    try:
        r = _post_json("https://search.rcsb.org/rcsbsearch/v2/query", q)
    except Exception:
        return []
    return sorted(h["identifier"] for h in r.get("result_set", []))


_PDB_GQL = """query($ids:[String!]!){entries(entry_ids:$ids){
 rcsb_id struct{title} exptl{method}
 rcsb_entry_info{resolution_combined}
 rcsb_accession_info{initial_release_date}
 citation{pdbx_database_id_DOI rcsb_authors year journal_abbrev title}
 polymer_entities{
   rcsb_polymer_entity{pdbx_description}
   entity_poly{pdbx_seq_one_letter_code_can}
   rcsb_entity_host_organism{ncbi_scientific_name}
   entity_src_gen{pdbx_host_org_scientific_name plasmid_name pdbx_host_org_vector_type}
   rcsb_polymer_entity_container_identifiers{auth_asym_ids
     reference_sequence_identifiers{database_accession database_name}}}}}"""


def pdb_construct_details(pdb_ids, chunk=50, sleep_s=0.1):
    """Per entry: method, resolution, primary citation, and for each polymer entity the
    expression host, plasmid, deposited sequence and UniProt cross-reference.

    Every deposition states the host and often the plasmid, so this recovers host and
    construct evidence from papers you cannot read. Queried in chunks of `chunk` IDs.

    Do NOT take construct boundaries from the deposited residue numbers: depositors
    often number against another isoform. Run canonical_range on each entity sequence
    instead.
    """
    out = []
    ids = list(pdb_ids)
    for i in range(0, len(ids), chunk):
        g = _post_json("https://data.rcsb.org/graphql",
                       {"query": _PDB_GQL, "variables": {"ids": ids[i:i + chunk]}})
        for e in (g.get("data", {}).get("entries") or []):
            if not e:
                continue
            cit = (e.get("citation") or [{}])[0]
            ents = []
            for pe in e["polymer_entities"]:
                src = (pe.get("entity_src_gen") or [{}])[0]
                hostorg = (pe.get("rcsb_entity_host_organism") or [{}])[0]
                ids_blk = pe["rcsb_polymer_entity_container_identifiers"]
                ents.append(dict(
                    description=pe["rcsb_polymer_entity"]["pdbx_description"],
                    chains=ids_blk["auth_asym_ids"],
                    host_raw=(hostorg.get("ncbi_scientific_name")
                              or src.get("pdbx_host_org_scientific_name")),
                    plasmid=src.get("plasmid_name"),
                    vector_type=src.get("pdbx_host_org_vector_type"),
                    sequence=pe["entity_poly"]["pdbx_seq_one_letter_code_can"],
                    xrefs=[(x["database_name"], x["database_accession"]) for x in
                           (ids_blk.get("reference_sequence_identifiers") or [])]))
            out.append(dict(
                pdb_id=e["rcsb_id"], title=e["struct"]["title"],
                method=(e.get("exptl") or [{}])[0].get("method"),
                resolution=((e.get("rcsb_entry_info") or {}).get("resolution_combined")
                            or [None])[0],
                released=(e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
                citation_doi=cit.get("pdbx_database_id_DOI"),
                citation_author=(cit.get("rcsb_authors") or [""])[0].split(",")[0],
                citation_year=cit.get("year"),
                citation_journal=cit.get("journal_abbrev"),
                citation_title=cit.get("title"),
                entities=ents))
        time.sleep(sleep_s)
    return out


# Host families, matched against the lower-cased host string, first hit wins.
HOST_FAMILIES = [
    ("cell-free", ("cell-free", "cell free", "in vitro", "wheat germ", "reticulocyte")),
    ("insect", ("spodoptera", "trichoplusia", "drosophila", "sf9", "sf21", "sf-9",
                "high five", "hi5", "baculovirus", "autographa")),
    ("mammalian", ("homo sapiens", "hek", "cricetulus", "cho", "mus musculus",
                   "chlorocebus", "cercopithecus", "cos-", "bos taurus")),
    ("yeast", ("saccharomyces", "pichia", "komagataella", "kluyveromyces",
               "schizosaccharomyces", "yarrowia")),
    ("E. coli", ("escherichia",)),
    ("other bacteria", ("bacillus", "lactococcus", "corynebacterium", "pseudomonas",
                        "mycobacterium", "brevibacillus")),
    ("plant", ("nicotiana", "arabidopsis", "oryza")),
]

# Vector-name fragments that imply a host family, used to cross-check the host field.
VECTOR_HINTS = {
    "insect": ("pfastbac", "fastbac", "pbev", "pvl13", "pacgp67", "pacsg", "pbac",
               "bac-to-bac", "baculo", "pomnibac"),
    "mammalian": ("pcdna", "phlsec", "ptt5", "pcep4", "pcaggs", "pef", "plvx",
                  "pcmv", "psecta", "pexpr"),
    "E. coli": ("pet", "pgex", "pmal", "plic", "pqe", "prset", "pcold", "psumo",
                "pbad", "ptrc", "ptwin", "pca24", "pdest", "pnic", "popin"),
    "yeast": ("ppic", "pyes", "pesc", "pgapz"),
    "cell-free": ("pcellfree", "pivex", "peu", "pf3a"),
}

NOT_STATED = "not stated"


def normalise_host(host_raw):
    """Collapse an RCSB expression-host string to a host family.

    RCSB records 'Escherichia coli', 'Escherichia coli BL21(DE3)' and
    'Escherichia coli BL21' as three distinct hosts; without this they count as
    three routes. Returns one of the HOST_FAMILIES labels, 'not stated', or the
    original string when nothing matches (so a genuinely new host is not hidden).
    """
    s = (host_raw or "").strip().lower()
    if not s or s in ("none", "null", "n/a", "unidentified"):
        return NOT_STATED
    for label, keys in HOST_FAMILIES:
        if any(k in s for k in keys):
            return label
    return host_raw.strip()


def host_from_vector(plasmid):
    """Host family implied by a plasmid/vector name, or None if it implies nothing."""
    s = (plasmid or "").strip().lower()
    if not s:
        return None
    for label, keys in VECTOR_HINTS.items():
        if any(s.startswith(k) or k in s for k in keys):
            return label
    return None


def target_entity(entry, accession=None):
    """The polymer entity of one pdb_construct_details record that is the target.

    An entry is often a complex -- a kinase with a substrate peptide, an antibody with
    its antigen -- and entity 1 is frequently the small partner, which carries no host
    at all. Picking entities[0] blindly therefore reports a large fraction of entries
    as 'host not stated'. Matches on the UniProt accession when given, else falls back
    to the first entity that records a host, else the first entity.
    """
    ents = entry.get("entities") or [{}]
    if accession:
        for e in ents:
            if any(a == accession for _, a in (e.get("xrefs") or [])):
                return e
    for e in ents:
        if e.get("host_raw"):
            return e
    return ents[0]


def summarise_hosts(details, accession=None):
    """Count host families across pdb_construct_details output and flag the metadata
    conflicts that are common in RCSB records.

    Pass the target `accession` so the host is read off the target's chain rather than
    a co-crystallised partner (see target_entity).

    Returns {counts: {family: n_entries}, by_family: {family: [pdb_ids]},
             conflicts: [{pdb_id, kind, host, plasmid}]}

    Three conflict kinds, all of which make a host assignment unverified:
      'vector_disagrees'  -- the plasmid implies a different host family
                             (a baculovirus vector recorded against E. coli)
      'vector_in_host'    -- the host-organism field holds a vector name
      'vector_no_host'    -- a plasmid is recorded but no host is

    Counted once per entry, on the target's entity.
    """
    counts, by_family, conflicts = {}, {}, []
    for d in details:
        ent = target_entity(d, accession)
        raw, plasmid = ent.get("host_raw"), ent.get("plasmid")
        fam = normalise_host(raw)
        implied = host_from_vector(plasmid)
        low = (raw or "").lower()
        if host_from_vector(raw) or "vector" in low or "cloning" in low:
            conflicts.append(dict(pdb_id=d["pdb_id"], kind="vector_in_host",
                                  host=raw, plasmid=plasmid))
        elif implied and fam not in (NOT_STATED, implied):
            conflicts.append(dict(pdb_id=d["pdb_id"], kind="vector_disagrees",
                                  host=raw, plasmid=plasmid))
        elif implied and fam == NOT_STATED:
            conflicts.append(dict(pdb_id=d["pdb_id"], kind="vector_no_host",
                                  host=raw, plasmid=plasmid))
        counts[fam] = counts.get(fam, 0) + 1
        by_family.setdefault(fam, []).append(d["pdb_id"])
    return dict(counts=dict(sorted(counts.items(), key=lambda kv: -kv[1])),
                by_family=by_family, conflicts=conflicts)


def canonical_range(entity_seq, canonical_seq, kmer=15):
    """Map a deposited entity sequence onto the canonical UniProt sequence.

    Depositors often number against another isoform while labelling the canonical
    accession, so neither the deposited residue numbers nor 'extra residues at the
    N terminus' can be read as a tag. Human PIM1 is the standard trap: most entries
    number against the 404-residue isoform while the canonical sequence is 313, and
    an entry carrying 96 extra N-terminal residues is that isoform, not a tag.

    Returns {start, end (1-based canonical numbering, None if no match),
             n_extra, n_term_extra, c_term_extra (the extra sequences themselves)}.
    """
    seq = entity_seq or ""
    hits = [i for i in range(len(seq) - kmer + 1) if seq[i:i + kmer] in canonical_seq]
    if not hits:
        return dict(start=None, end=None, n_extra=len(seq), n_term_extra="", c_term_extra="")
    a, b = hits[0], hits[-1] + kmer
    return dict(start=canonical_seq.find(seq[a:a + kmer]) + 1,
                end=canonical_seq.find(seq[b - kmer:b]) + kmer,
                n_extra=a + len(seq) - b,
                n_term_extra=seq[:a], c_term_extra=seq[b:])


# Tag fragments recognised in the extra residues canonical_range returns.
TAG_SIGNATURES = [("His", "HHHHHH"), ("3C/PreScission scar", "GP"),
                  ("TEV scar", "ENLYFQ"), ("thrombin scar", "LVPRGS"),
                  ("SUMO", "LQDSEV"), ("Strep", "WSHPQFEK"), ("FLAG", "DYKDDDDK")]


def tag_status(canonical, isoform_seqs=(), min_isoform=8):
    """Classify one entry's tag from canonical_range's result.

    Extra residues at an end are not automatically a tag. When the construct came
    from a longer splice isoform, the isoform's own extension sits outside the
    canonical sequence and reads as a large tag -- the human PIM1 trap, where entries
    solved from the 404-residue isoform carry ~91 'extra' N-terminal residues. Pass
    `isoform_seqs` (from fetch_isoform_sequences) and the longest end stretch that is
    really isoform sequence is removed before anything is called a tag, so a His-
    tagged isoform construct is reported as both.

    Returns {call, n_term_extra, c_term_extra, isoform_extra, tags}:
      call  'intact tag' (8+ residues left at an end, or a named tag such as His6
            matched outright), 'scar' (a cleaved tag's short remnant, such as GP from
            HRV-3C/PreScission), 'isoform' (everything extra was isoform sequence),
            or 'none'.
      tags  names matched against TAG_SIGNATURES, e.g. ['His', '3C/PreScission scar'].
    """
    n = canonical.get("n_term_extra") or ""
    c = canonical.get("c_term_extra") or ""
    iso_hit = ""
    for seq in isoform_seqs:
        for k in range(len(n), min_isoform - 1, -1):
            if len(n[-k:]) > len(iso_hit) and n[-k:] in seq:
                iso_hit = n[-k:]
                break
    n_rest = n[:len(n) - len(iso_hit)] if iso_hit else n
    longest = max(len(n_rest), len(c))
    found = [name for name, sig in TAG_SIGNATURES if sig in n_rest or sig in c]
    if longest >= 8 or [f for f in found if not f.endswith("scar")]:
        call = "intact tag"
    elif longest >= 1:
        call = "scar"
    else:
        call = "isoform" if iso_hit else "none"
    return dict(call=call, n_term_extra=n, c_term_extra=c,
                isoform_extra=iso_hit, tags=found)


def construct_families(ranges, tol=3):
    """Cluster construct ranges that differ by at most `tol` residues at each end
    (1-313, 1-312 and 2-313 are one construct with or without the initiator Met).

    ranges: {pdb_id: (start, end)} in canonical numbering, from canonical_range.
    Returns families largest first: [{label, start, end, members, ranges}].
    Entries with no mapping form a family labelled 'unmapped'.
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
            fams.append(dict(label=f"{s}–{e}", start=s, end=e,
                             ranges={(s, e): counts[(s, e)]}))
    for f in fams:
        f["members"] = sorted(p for p, r in ranges.items() if r in f["ranges"])
    unmapped = sorted(p for p, r in ranges.items() if r[0] is None)
    if unmapped:
        fams.append(dict(label="unmapped", start=None, end=None, ranges={}, members=unmapped))
    return sorted(fams, key=lambda f: -len(f["members"]))


# ---------------------------------------------------------------- literature

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def pubmed_search(term, retmax=15, sort="relevance", email=None, tool="expression-host-strategy"):
    """Field-restricted PubMed query -> [{pmid, year, first_author, journal, title, doi}].

    `email` is optional and OFF by default. NCBI asks for a contact address so it can
    warn you before blocking heavy use; passing one sends it to NCBI, so only pass it
    when the user has agreed to that. Unauthenticated use is rate-limited to about
    3 requests/second -- fine for a handful of themed queries.

    For a heavily studied target expect this to return disease and inhibitor papers
    and nothing about producing the protein; report that rather than padding.
    """
    p = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json",
         "sort": sort, "tool": tool}
    if email:
        p["email"] = email
    r = _get_json(_EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(p))
    ids = r["esearchresult"]["idlist"]
    if not ids:
        return []
    p2 = {"db": "pubmed", "id": ",".join(ids), "retmode": "json", "tool": tool}
    if email:
        p2["email"] = email
    res = _get_json(_EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(p2))["result"]
    out = []
    for i in ids:
        d = res.get(i, {})
        doi = [x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"]
        out.append(dict(pmid=i, year=d.get("pubdate", "")[:4],
                        first_author=(d["authors"][0]["name"] if d.get("authors") else ""),
                        journal=d.get("source", ""), title=d.get("title", ""),
                        doi=doi[0] if doi else ""))
    return out


def pubmed_abstracts(pmids, email=None, tool="expression-host-strategy"):
    """{pmid: {title, abstract}} via efetch. See pubmed_search on `email`."""
    import xml.etree.ElementTree as ET
    p = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
         "rettype": "abstract", "tool": tool}
    if email:
        p["email"] = email
    req = urllib.request.Request(_EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(p))
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        root = ET.fromstring(r.read())
    out = {}
    for art in root.findall(".//PubmedArticle"):
        out[art.find(".//PMID").text] = dict(
            title="".join(art.find(".//ArticleTitle").itertext()),
            abstract=" ".join("".join(a.itertext())
                              for a in art.findall(".//Abstract/AbstractText")))
    return out


# ----------------------------------------------------------------- citations

def _fold(s):
    """Lower-case, strip accents and punctuation, so 'Akue-Gedu' matches 'Akue-Gedu'
    written with its accents and a typographic hyphen."""
    import unicodedata
    d = unicodedata.normalize("NFKD", str(s or ""))
    return re.sub(r"[^a-z0-9]", "", "".join(c for c in d if not unicodedata.combining(c)).lower())


def check_dois(dois, claims=None, pause=0.15):
    """CrossRef check for every DOI before it goes in a deliverable.

    Returns {doi: {ok, first_author, year, title, journal, retracted, error?}}.

    A resolving DOI is NOT evidence that it supports the claim: a real paper about
    the right protein passes while being attributed to the wrong author or year.
    Pass `claims` as {doi: (first_author, year)} -- from pdb_citation_claims, or
    whatever you are about to write -- and each record also gets `claim_ok` and,
    when it fails, `claim_mismatch` naming what CrossRef says instead. Fix the
    deliverable from CrossRef, never the other way round.
    """
    out = {}
    for d in dois:
        if not d:
            continue
        try:
            m = _get_json("https://api.crossref.org/works/"
                          + urllib.parse.quote(d, safe=""))["message"]
            upd = " ".join(filter(None, (u.get("type") for u in (m.get("update-to") or []))))
            rec = dict(ok=True,
                       first_author=(m.get("author") or [{}])[0].get("family", ""),
                       year=(m.get("issued", {}).get("date-parts") or [[None]])[0][0],
                       title=html.unescape((m.get("title") or [""])[0]),
                       journal=html.unescape((m.get("container-title") or [""])[0]),
                       retracted=bool("retraction" in upd.lower()
                                      or m.get("relation", {}).get("is-retracted-by")))
        except Exception as exc:
            rec = dict(ok=False, error=str(exc)[:120])
        claim = (claims or {}).get(d)
        if claim and rec.get("ok"):
            author, year = claim
            a, b = _fold(author), _fold(rec["first_author"])
            a_ok = (not a) or (bool(b) and (a in b or b in a))
            y_ok = (not year) or abs(int(year) - int(rec["year"] or 0)) <= 1
            rec["claim_ok"] = bool(a_ok and y_ok)
            if not rec["claim_ok"]:
                rec["claim_mismatch"] = (f"claimed {author} {year}; "
                                         f"CrossRef says {rec['first_author']} {rec['year']}")
        out[d] = rec
        time.sleep(pause)
    return out


def pdb_citation_claims(details):
    """{doi: (first_author, year)} from pdb_construct_details, to feed check_dois.

    The PDB primary-citation record and CrossRef are independent, so disagreement
    between them is the cheapest catch there is for a misattributed reference.
    """
    return {d["citation_doi"]: (d.get("citation_author") or "", d.get("citation_year"))
            for d in details if d.get("citation_doi")}


# ------------------------------------------------------------------ graphics

HOST_HUES = {"E. coli": "#2a7f62", "insect": "#8a6bbe", "mammalian": "#d9731a",
             "yeast": "#b07aa1", "cell-free": "#4e79a7", "other bacteria": "#59a14f",
             "plant": "#7fa650", NOT_STATED: "#8c8c8c"}
SEGMENT_HUES = {"tag": "#3a3a3a", "linker": "#b8b8bd", "signal": "#d0d0d0",
                "propeptide": "#c7b299", "fusion": "#9fb8d6", "partner": "#e0a526",
                "target": "#1f5a96", "domain": "#4e79a7"}
_svg_uid = [0]


def track_svg(tracks, axis_len, axis_start=0, axis_label="residue", row_h=38,
              width=1000, left=190, title=None):
    """Draw construct/anatomy tracks as one inline SVG string for write_html_report.

    tracks: [{label, note, segments: [...], marks: [...]}]
      note is a few words; it sits in a fixed right gutter, so keep it short.
      A segment label is dropped when the segment is too narrow to hold it.
      segment: {start, end, label, color ('tag'/'target'/... or '#rrggbb'),
                hatch: True when the length is not stated -- never guess a tag length}
      mark:    {pos, label, filled: True for an occupied glycan / False for a sequon
                only, color}
    axis_len: the residue number the axis ends at, so every track shares one scale.
    axis_start: where the axis begins; pass a negative number to leave room for an
      N-terminal tag drawn before residue 1 (segments outside the axis are clipped).

    Scales to the page width and follows light/dark mode; it prints with the report.
    """
    _svg_uid[0] += 1
    uid = f"h{_svg_uid[0]}"
    # reserve a right gutter when any track carries a note, so bars never run under it
    x0 = left
    x1 = width - (128 if any(t.get("note") for t in tracks) else 24)
    lo, hi = axis_start, max(axis_start + 1, axis_len)
    span = hi - lo
    top = 34 if title else 14
    px = lambda p: x0 + ((max(lo, min(hi, p)) - lo) / span) * (x1 - x0)
    height = top + row_h * len(tracks) + 34
    hue = lambda c: SEGMENT_HUES.get(c, c or "#4e79a7")

    step = 10 ** max(0, len(str(int(span // 6))) - 1) or 1
    step *= next(m for m in (1, 2, 5, 10) if span / (step * m) <= 8)

    p = [f'<svg class="track" viewBox="0 0 {width} {height}" width="100%" '
         f'role="img" aria-label="{html.escape(title or axis_label)}">',
         f'<defs><pattern id="{uid}" width="6" height="6" patternUnits="userSpaceOnUse" '
         f'patternTransform="rotate(45)"><rect width="6" height="6" fill="var(--bg)"/>'
         f'<line x1="0" y1="0" x2="0" y2="6" stroke="currentColor" stroke-width="2.5" '
         f'opacity=".45"/></pattern></defs>']
    if title:
        p.append(f'<text x="8" y="18" class="t-title">{html.escape(title)}</text>')

    for i, tr in enumerate(tracks):
        y = top + i * row_h
        bar = y + row_h - 20
        p.append(f'<text x="8" y="{bar + 11}" class="t-lab">{html.escape(tr["label"])}</text>')
        if tr.get("note"):
            p.append(f'<text x="{x1 + 8}" y="{bar + 11}" class="t-note">'
                     f'{html.escape(tr["note"])}</text>')
        p.append(f'<line x1="{x0}" y1="{bar + 7}" x2="{x1}" y2="{bar + 7}" class="t-base"/>')
        for sg in tr.get("segments", []):
            a, b = px(sg["start"]), px(sg["end"])
            w = max(2.0, b - a)
            fill = hue(sg.get("color"))
            p.append(f'<rect x="{a:.1f}" y="{bar - 3}" width="{w:.1f}" height="20" rx="2.5" '
                     f'fill="{fill}" opacity=".92"/>')
            if sg.get("hatch"):
                p.append(f'<rect x="{a:.1f}" y="{bar - 3}" width="{w:.1f}" height="20" rx="2.5" '
                         f'fill="url(#{uid})"/>')
            # draw the label only when it actually fits, so it can never spill
            # over a neighbouring segment or the row label
            if sg.get("label") and w > 6.2 * len(str(sg["label"])) + 6:
                p.append(f'<text x="{(a + b) / 2:.1f}" y="{bar + 11}" class="t-seg" '
                         f'text-anchor="middle">{html.escape(str(sg["label"]))}</text>')
        for mk in tr.get("marks", []):
            x = px(mk["pos"])
            c = mk.get("color", "#b42318")
            p.append(f'<line x1="{x:.1f}" y1="{bar - 14}" x2="{x:.1f}" y2="{bar - 3}" '
                     f'stroke="{c}" stroke-width="1.4"/>'
                     f'<circle cx="{x:.1f}" cy="{bar - 16}" r="3.4" fill="'
                     + (c if mk.get("filled", True) else "var(--bg)")
                     + f'" stroke="{c}" stroke-width="1.4"/>')
            if mk.get("label"):
                p.append(f'<text x="{x:.1f}" y="{bar - 22}" class="t-mark" '
                         f'text-anchor="middle">{html.escape(str(mk["label"]))}</text>')

    axis_y = top + row_h * len(tracks) + 6
    p.append(f'<line x1="{x0}" y1="{axis_y}" x2="{x1}" y2="{axis_y}" class="t-axis"/>')
    first = int(step) * -(-lo // int(step))
    for t in range(first, hi + 1, int(step)):
        x = px(t)
        p.append(f'<line x1="{x:.1f}" y1="{axis_y}" x2="{x:.1f}" y2="{axis_y + 4}" class="t-axis"/>'
                 f'<text x="{x:.1f}" y="{axis_y + 16}" class="t-tick" '
                 f'text-anchor="middle">{t}</text>')
    p.append(f'<text x="{x1}" y="{axis_y + 28}" class="t-tick" text-anchor="end">'
             f'{html.escape(axis_label)}</text></svg>')
    return "\n".join(p)


def swatch_legend(items):
    """Inline legend: items = [(label, colour key or '#rrggbb')]. Pair with track_svg."""
    return '<p class="legend">' + "".join(
        f'<span class="key"><span class="dot" style="--g:{SEGMENT_HUES.get(c, c)}">'
        f'</span>{html.escape(l)}</span>' for l, c in items) + "</p>"


# -------------------------------------------------------------------- report

GROUP_HUES = ["#4e79a7", "#f28e2b", "#59a14f", "#b07aa1", "#76b7b2", "#edc948",
              "#ff9da7", "#9c755f"]
OTHER_HUE = "#8c8c8c"
STATUS_PILLS = {  # first word of a pill -> class; anything else is neutral
    "intact": "warn", "scar": "info", "isoform": "warn", "none": "muted",
    "verified": "good", "unverified": "warn", "conflict": "bad", "retracted": "bad",
    "mismatch": "bad", "n.d.": "muted", "yes": "good", "no": "bad", "✓": "good",
    "✗": "bad", "?": "muted", "recommended": "good", "back-up": "warn",
    "rejected": "bad",
}


def write_html_report(path, title, sections, sources=(), subtitle=None):
    """Write the strategy as one self-contained HTML page (no external files) that
    opens in any browser, follows light/dark mode and prints cleanly to PDF.

    Each section is a dict with `heading` and any of:
        text:    paragraphs; blank lines separate them, '- ' lines become bullets,
                 **bold** is kept, [label](url) becomes a link
        svg:     an inline SVG string from track_svg, or a list of them, optionally
                 with a legend from swatch_legend
        table:   a flat list of row dicts
        groups:  [{key, label, note, rows}] -- one block per group: a tinted header
                 band over its own table. Use it for one block per host family.
        pills:   column names rendered as status pills (tag call, host verdict,
                 verified/unverified). A cell holding a list renders one pill per item.
        swatch:  column name that gets a coloured dot from each row's `_swatch` key
        collapse: wrap each group in a <details> block, open if it has >1 row;
                  all open when printed
    Keys starting with '_' are never shown as columns.

    Colours: every group `key` and row `_swatch` is a string like 'host:E. coli' or
    'fam:29-313'. Within one prefix, keys take hues in order of first appearance and
    keep that hue everywhere in the report, so a host's colour matches between the
    summary and the entry list. A key of the form 'host:<family>' uses HOST_HUES.
    After 8 hues, and for keys ending ':other', it is grey.
    """
    import datetime

    hues = {}

    def hue(key):
        if key is None:
            return OTHER_HUE
        if key not in hues:
            prefix, _, rest = key.partition(":")
            if prefix == "host" and rest in HOST_HUES:
                hues[key] = HOST_HUES[rest]
            else:
                used = sum(1 for k in hues
                           if k.split(":", 1)[0] == prefix and hues[k] != OTHER_HUE)
                hues[key] = (OTHER_HUE if key.endswith(":other") or used >= len(GROUP_HUES)
                             else GROUP_HUES[used])
        return hues[key]

    def inline(s):
        s = html.escape(s)
        s = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2">\1</a>', s)
        return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)

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
            body = " ".join(pill(x) for x in (v if isinstance(v, (list, tuple)) else [v])
                            if x not in (None, ""))
        else:
            raw = ("" if v is None else ", ".join(map(str, v))
                   if isinstance(v, (list, tuple)) else str(v))
            body = inline(raw)
        if col == s.get("swatch") and row.get("_swatch"):
            body = f'<span class="dot" style="--g:{hue(row["_swatch"])}"></span>{body}'
        return f"<td>{body}</td>"

    def table(rows, s):
        cols = list(dict.fromkeys(c for r in rows for c in r if not c.startswith("_")))
        head = "".join(f"<th>{html.escape(c.replace('_', ' '))}</th>" for c in cols)
        body = "".join("<tr>" + "".join(cell(r, c, s) for c in cols) + "</tr>" for r in rows)
        return (f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
                f"<tbody>{body}</tbody></table></div>")

    def group(g, s):
        band = (f'<span class="dot" style="--g:{hue(g.get("key"))}"></span>'
                f'<strong>{inline(g["label"])}</strong>'
                + (f' <span class="note">{inline(g["note"])}</span>' if g.get("note") else ""))
        inner = table(g["rows"], s) if g.get("rows") else ""
        style = f'style="--g:{hue(g.get("key"))}"'
        if s.get("collapse"):
            opened = " open" if len(g.get("rows") or []) > 1 else ""
            return (f'<details class="grp" {style}{opened}><summary class="band">'
                    f"{band}</summary>{inner}</details>")
        return f'<div class="grp" {style}><div class="band">{band}</div>{inner}</div>'

    parts = []
    for s in sections:
        parts.append(f"<h2>{html.escape(s['heading'])}</h2>")
        if s.get("text"):
            parts.append(text_block(s["text"]))
        if s.get("svg"):
            svgs = s["svg"] if isinstance(s["svg"], (list, tuple)) else [s["svg"]]
            parts.extend(svgs)
        if s.get("table"):
            parts.append(table(s["table"], s))
        for g in s.get("groups") or []:
            parts.append(group(g, s))
    if sources:
        parts.append("<h2>Sources</h2><ul>" + "".join(
            f'<li><a href="{html.escape(u)}">{html.escape(l)}</a></li>'
            for l, u in sources) + "</ul>")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{ --bg:#fff; --fg:#1d1d1f; --muted:#6e6e73; --line:#d9d9de; --head:#f4f4f6; --link:#0a58ca;
  --good:#1e7b34; --warn:#a15c00; --bad:#b42318; --info:#1f5fa8; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#16161a; --fg:#ececf1; --muted:#a1a1aa; --line:#33333b;
  --head:#222228; --link:#7cb3ff; --good:#6fcf86; --warn:#f0b35a; --bad:#ff8a80; --info:#8cbcff; }} }}
* {{ -webkit-print-color-adjust:exact; print-color-adjust:exact; }}
body {{ background:var(--bg); color:var(--fg); font:15px/1.5 -apple-system, "Segoe UI", Roboto, sans-serif;
  max-width:1100px; margin:0 auto; padding:24px 16px; }}
h1 {{ font-size:1.6em; margin:0 0 4px; }}
h2 {{ font-size:1.2em; margin:32px 0 10px; border-bottom:1px solid var(--line); padding-bottom:4px; }}
.meta, .note {{ color:var(--muted); font-size:.9em; }} a {{ color:var(--link); }}
.scroll {{ overflow-x:auto; }} table {{ border-collapse:collapse; width:100%; font-size:.85em; }}
th, td {{ border:1px solid var(--line); padding:4px 6px; text-align:left; vertical-align:top; }}
th {{ background:var(--head); }}
.grp {{ margin:10px 0; border-left:6px solid var(--g); }}
.grp table {{ border-left:none; }}
.band {{ display:block; padding:6px 10px; background:color-mix(in srgb, var(--g) 18%, var(--bg)); cursor:default; }}
summary.band {{ cursor:pointer; }}
.dot {{ display:inline-block; width:.75em; height:.75em; border-radius:50%; background:var(--g);
  margin-right:6px; vertical-align:baseline; }}
.pill {{ display:inline-block; border:1px solid currentColor; border-radius:9px; padding:0 7px;
  margin:1px 2px 1px 0; font-size:.85em; white-space:nowrap; }}
.pill.good {{ color:var(--good); }} .pill.warn {{ color:var(--warn); }} .pill.bad {{ color:var(--bad); }}
.pill.info {{ color:var(--info); }} .pill.muted {{ color:var(--muted); }}
.legend {{ font-size:.85em; color:var(--muted); margin:2px 0 14px; }}
.key {{ margin-right:14px; white-space:nowrap; }}
svg.track {{ display:block; margin:6px 0 2px; max-width:100%; height:auto; color:var(--fg); }}
.t-title {{ font-size:14px; font-weight:600; fill:var(--fg); }}
.t-lab {{ font-size:12px; fill:var(--fg); }}
.t-note, .t-tick {{ font-size:11px; fill:var(--muted); }}
.t-mark {{ font-size:10px; fill:var(--muted); }}
.t-seg {{ font-size:11px; fill:#fff; font-weight:500; }}
.t-base {{ stroke:var(--line); stroke-width:1; }}
.t-axis {{ stroke:var(--muted); stroke-width:1; }}
@page {{ size:A4; margin:12mm; }}
@media print {{ body {{ max-width:none; padding:0; font-size:10.5px; }} .scroll {{ overflow:visible; }}
  tr, .band, svg.track, .grp {{ break-inside:avoid; }} h2 {{ break-after:avoid; }}
  summary {{ list-style:none; }} summary::-webkit-details-marker {{ display:none; }} }}
</style>
<!-- CSS cannot reveal a closed <details>, so this is what keeps collapsed groups
     from being dropped when the page is printed or converted to PDF. -->
<script>addEventListener("beforeprint", () => document.querySelectorAll("details").forEach(d => d.open = true));</script>
</head><body>
<h1>{html.escape(title)}</h1>
<p class="meta">{html.escape(subtitle or "Expression-host strategy")} · generated {datetime.date.today().isoformat()}</p>
{chr(10).join(parts)}
</body></html>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path


CHROMIUM_PATHS = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
    "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
    "microsoft-edge",
]


def find_chromium():
    """Path to an installed Chromium-family browser, or None."""
    import os
    import shutil
    for c in CHROMIUM_PATHS:
        if os.path.sep in c:
            if os.path.exists(c):
                return c
        else:
            found = shutil.which(c)
            if found:
                return found
    return None


def html_to_pdf(html_path, pdf_path=None, timeout=120):
    """Print a report written by write_html_report to PDF, A4, via headless Chrome.

    Deliver both files: the HTML is the one to read and click through, the PDF the
    one to forward or attach. No pandoc, LaTeX or Python PDF package is involved --
    the page already carries its own print stylesheet, and the print rules open every
    collapsed group so nothing is silently dropped from the PDF.

    Returns the PDF path, or None when no Chromium-family browser is installed; in
    that case tell the user to open the HTML and print it (Cmd/Ctrl-P, Save as PDF),
    which produces the same page.
    """
    import os
    import subprocess
    import tempfile

    exe = find_chromium()
    if not exe:
        return None
    html_path = os.path.abspath(html_path)
    pdf_path = os.path.abspath(pdf_path or os.path.splitext(html_path)[0] + ".pdf")
    with tempfile.TemporaryDirectory() as profile:
        subprocess.run([exe, "--headless=new", "--disable-gpu", "--no-sandbox",
                        f"--user-data-dir={profile}", "--no-pdf-header-footer",
                        f"--print-to-pdf={pdf_path}",
                        "file://" + urllib.parse.quote(html_path)],
                       check=True, timeout=timeout,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return pdf_path if os.path.exists(pdf_path) else None
