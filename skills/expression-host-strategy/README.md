# expression-host-strategy

Decide which host to make a protein in — *E. coli*, insect, mammalian HEK/CHO,
yeast or cell-free — and produce the decision as a self-contained HTML
page plus a matching A4 PDF.

## What it does

The question "which host should I express protein X in?" is usually treated
as a literature question. For any protein with deposited structures it is
mostly a **data** question: every PDB deposition records the expression host
and often the plasmid, so a few hundred entries answer it directly, and the
literature only has to fill in yields and the routes nobody deposited.

So the skill works in that order:

1. **UniProt** — sequence and the features that decide the host: signal
   peptide, propeptide, disulfides, glycosylation sequons, splice isoforms.
   `construct_props` gives length, mass, pI, Cys count and ε280 for the exact
   chain you will order.
2. **PDB** — one batched GraphQL query returns host, plasmid, deposited
   sequence and citation for every entry. `summarise_hosts` collapses host
   strings into families and counts them.
3. **Alignment, not residue numbers** — construct boundaries come from
   aligning each deposited sequence to the canonical one, never from the
   deposited numbering.
4. **PubMed** — themed, field-restricted queries for yields and the routes
   the PDB cannot show.
5. **CrossRef** — every DOI checked against the claim it is cited for.
6. **Report** — one HTML file with construct maps as inline SVG, host
   precedent, a route table with verdicts, the CRO brief and every entry
   grouped by host; then the same page printed to A4 PDF by `html_to_pdf`,
   so there is only ever one document to keep right.

## Why these pieces exist

Four failure modes drove the design, each found by running the skill on a real
target rather than imagined:

- **Host strings are not host families.** RCSB records `Escherichia coli`,
  `Escherichia coli BL21(DE3)` and `Escherichia coli BL21` as three different
  hosts. Counting raw strings splits one route into several.
- **Entity 1 is often the wrong chain.** For a kinase solved with a substrate
  peptide, the first polymer entity is the peptide, which carries no host at
  all. `target_entity` picks the chain that references the accession;
  without it a quarter of entries read as "host not stated".
- **Extra N-terminal residues are not automatically a tag.** Depositors
  frequently number against a longer splice isoform, and its extension then
  looks like a large affinity tag. Human PIM1 is the canonical trap: the
  canonical sequence is 313 residues, most entries number against the
  404-residue isoform, and entries carrying "91 extra residues" are that
  isoform. `tag_status` strips a matched isoform extension before calling
  anything a tag, and reports both when a construct has a tag *and* came from
  the longer isoform.
- **A resolving DOI is not a verified citation.** A real paper about the right
  protein passes a plain DOI check while being credited to the wrong author.
  `check_dois` takes a `claims` map and compares CrossRef's first author and
  year against it; `pdb_citation_claims` builds that map from the PDB's own
  citation records, which are an independent source.

It also reports, rather than silently resolving, the metadata conflicts that
are common in RCSB: baculovirus vectors recorded against an *E. coli* host, a
vector name sitting in the host-organism field, a plasmid with no host.

## Requirements

Python 3, standard library only. No `requests`, no Biopython, no matplotlib,
no pandoc, no LaTeX, no Python PDF package. Network access to
`rest.uniprot.org`, `search.rcsb.org`, `data.rcsb.org`,
`eutils.ncbi.nlm.nih.gov` and `api.crossref.org`.

The PDF is printed from the HTML by a headless Chromium-family browser
(Chrome, Edge, Chromium or Brave), which `html_to_pdf` locates itself. With
none installed it returns `None` and the HTML prints to the same PDF from the
browser's own Cmd/Ctrl-P → Save as PDF.

`pubmed_search` takes an optional `email` for NCBI's contact header. It is off
by default and only ever sends an address you pass explicitly.

## Files

| File | What it is |
|---|---|
| `SKILL.md` | The workflow Claude follows |
| `scripts/expression.py` | All helpers — fetching, normalising, alignment, SVG, HTML report |

## Using the helpers directly

```python
import sys; sys.path.insert(0, "skills/expression-host-strategy/scripts")
import expression as ex

entry = ex.fetch_uniprot_entry("P11309")
details = ex.pdb_construct_details(ex.pdb_entries_for_uniprot("P11309"))
hosts = ex.summarise_hosts(details, accession="P11309")
print(hosts["counts"])        # {'E. coli': 182, 'cell-free': 8, 'insect': 2}
print(len(hosts["conflicts"]))  # 9 entries whose host disagrees with their plasmid

iso = list(ex.fetch_isoform_sequences("P11309").values())
rng = {d["pdb_id"]: ex.canonical_range(ex.target_entity(d, "P11309")["sequence"],
                                       entry["sequence"]) for d in details}
tags = {p: ex.tag_status(c, isoform_seqs=iso) for p, c in rng.items()}

claims = ex.pdb_citation_claims(details)
checked = ex.check_dois(claims, claims=claims)   # claim_ok on every record
```

The 192-entry PIM1 survey above takes about four seconds, and

```python
ex.write_html_report(path, title, sections, sources)
ex.html_to_pdf(path)
```

produces the two deliverables.

## Known limitations

- One host is counted per PDB entry. For a complex whose chains were made in
  different hosts, read the per-entity records in `pdb_construct_details`.
- `normalise_host` matches a fixed keyword list; an unrecognised host is
  returned unchanged rather than hidden, so scan the counts and extend
  `HOST_FAMILIES` if something should have collapsed.
- `construct_props`' pI is a Henderson-Hasselbalch estimate on the
  Expasy/ProtParam pKa set. It tracks ProtParam closely but is not a
  measurement and ignores structural pKa shifts.
- A host in the PDB is a host the protein was *made* in, not evidence of
  yield. Depositions tell you what folds; only methods sections say how much.
