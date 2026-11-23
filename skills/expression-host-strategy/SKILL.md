---
name: expression-host-strategy
description: Decide which host to express a protein in (E. coli vs insect vs mammalian HEK/CHO vs yeast vs cell-free) and deliver the decision as a self-contained HTML report plus a matching A4 PDF, with construct maps, a host-precedent breakdown and a per-route reference table. Surveys UniProt, the PDB and PubMed for host-specific production precedent, recovers construct boundaries and tags from PDB depositions, and cross-checks every DOI against the claim it is cited for. Use when asked which host or expression system to make a protein in, for a production or purification strategy, or for a host comparison before designing constructs.
---

# Expression-host strategy survey

Answers "which host should I make protein X in, and what exactly do I order?"
The deliverable is two files, both named after the protein:

- `<protein>_expression_strategy.html` — one self-contained page, no external
  files, that opens in any browser; the one to read and click through.
- `<protein>_expression_strategy.pdf` — the same page as A4, for forwarding
  or attaching.

The PDF is printed from the HTML by `html_to_pdf`, so there is only ever one
document to keep right.

This is mostly a data-mining task. The PDB records the expression host for
every deposition, so for a protein with structures the host question is
usually already answered by machine-readable data; the literature fills in
yield and the routes nobody deposited. Go to papers to fill gaps, not to
start.

## Requirements

- Python 3, standard library only. The helpers are in `scripts/expression.py`
  (bundled with this skill); import it with
  `sys.path.insert(0, "<this skill>/scripts"); import expression`.
- Network access to rest.uniprot.org, search.rcsb.org, data.rcsb.org,
  eutils.ncbi.nlm.nih.gov and api.crossref.org. If this environment can't
  reach them, say so and stop rather than answering from memory.
- For the PDF only: any Chromium-family browser (Chrome, Edge, Chromium,
  Brave), which `html_to_pdf` finds itself. No pandoc, LaTeX or Python PDF
  package is involved. If none is installed it returns `None`; deliver the
  HTML and tell the user to print it with Cmd/Ctrl-P → Save as PDF, which
  produces the same page.

## Ask before assuming

This skill does not know the user's platform. Unless they have already said,
ask once, in one message, before step 5 — and carry the answers through the
report rather than inventing a house style:

- Which hosts they actually have (E. coli strains, insect, HEK/CHO, yeast,
  cell-free), and whether any are off the table.
- Whether refolding from inclusion bodies is acceptable.
- Preferred tag and protease site, and for a secreted protein which signal
  peptide — never choose one for them.
- Whether catalogue or GMP material would do instead of making it.

If they don't answer, use these defaults and say in the report that they are
defaults: soluble expression preferred over refolding; standard E. coli
strains (BL21(DE3), BL21-AI) count as routine while Rosetta-gami, SHuffle or
ArcticExpress make a result a caveat rather than a recommendation; insect is a
back-up arm unless it shows a real yield or quality advantage.

## Workflow

### 1. Resolve the molecule before searching

```python
import sys; sys.path.insert(0, "<this skill>/scripts"); import expression as ex
entry = ex.fetch_uniprot_entry("P11309")
ex.construct_props(entry["sequence"])
```

`fetch_uniprot_entry` returns the sequence, the features that decide the host
(Signal, Propeptide, Chain, Domain, Glycosylation, Disulfide bond, Alternative
sequence) and the PDB cross-references. Settle these first, because they
decide the host:

- **Mature-chain boundary.** Signal peptide *and* propeptide both come off.
  State the mature range in precursor numbering and give the offset, because
  the literature usually numbers the mature chain.
- **Disulfides and free cysteines** — how hard the folding problem is.
- **Glycosylation sequons, and which are actually occupied.** A sequon is not
  an occupied site. If no glycan is required for activity, a prokaryotic host
  stays on the table.
- **Obligate partner chains.** Some proteins do not fold or secrete alone;
  this can decide the whole strategy.
- `construct_props(chain)` on each chain you will actually order — length,
  mass, pI, Cys count, ε280 reduced and oxidised, and N-X-S/T sequons.

### 2. PDB depositions are host and construct data

```python
ids = ex.pdb_entries_for_uniprot("P11309")
details = ex.pdb_construct_details(ids)
hosts = ex.summarise_hosts(details)
```

Every deposition states the expression host and often the plasmid. Do this
**before** the literature sweep: it is the densest host evidence there is, and
it tells you what the literature still has to answer.

- `summarise_hosts` collapses host strings to families. Do not count raw
  strings: RCSB records `Escherichia coli`, `Escherichia coli BL21(DE3)` and
  `Escherichia coli BL21` as three different hosts.
- It also returns **`conflicts`**, which you must report rather than quietly
  resolve. RCSB host metadata disagrees with the plasmid more often than
  people expect — baculovirus vectors recorded against an *E. coli* host, a
  vector name sitting in the host-organism field, a plasmid with no host.
  Call those entries' host assignment unverified and say how many there are.
- If every structure was made in one host, the folding problem is solvable
  there even when the yield literature is thin. That is a real finding; say it
  plainly with the count.

### 3. Construct boundaries and tags

**Never take construct boundaries from the deposited residue numbers.**
Depositors often number against a different isoform while labelling the
canonical accession, and then "96 extra residues at the N terminus" is the
other isoform, not a His tag. Map the sequence instead:

```python
can = entry["sequence"]
rng  = {d["pdb_id"]: ex.canonical_range(d["entities"][0]["sequence"], can) for d in details}
tags = {p: ex.tag_status(c) for p, c in rng.items()}
fams = ex.construct_families({p: (c["start"], c["end"]) for p, c in rng.items()})
```

`tag_status` returns `intact tag` (8+ extra residues), `scar` (1–7, a cleaved
tag's remnant — GP is the HRV-3C/PreScission signature), `none`, or
`isoform` when the extra residues are themselves isoform sequence. Pass the
other isoforms' sequences as `isoform_seqs` whenever UniProt lists more than
one, or the isoform case is miscalled as a tag. `construct_families` clusters
ranges within 3 residues, which is usually 2–6 families across a protein's
whole deposition history.

Report tag status per family explicitly: entries that crystallised or were
assayed with the tag still on tell the user cleaving it is optional.

### 4. Literature, sized to the target

Themed, field-restricted PubMed queries, one per theme, not one broad query:
E. coli + inclusion body/refold; E. coli + soluble/SUMO/MBP/thioredoxin/
periplasm; HEK293/CHO/mammalian + purification/yield; baculovirus/Sf9/Hi5;
fusion architecture; glycosylation; signal peptide/secretion; the obligate
partner chain if there is one.

```python
hits = ex.pubmed_search("PIM1 AND (expression OR purification) AND Escherichia coli[tiab]")
ex.pubmed_abstracts([h["pmid"] for h in hits[:8]])
```

**Calibrate the expectation to how studied the target is.** For a heavily
studied drug target the sweeps return disease biology and inhibitor chemistry
and nothing about producing the protein — the PDB has already given the whole
answer. Report that the literature added nothing rather than padding the
report with papers that do not bear on the host choice. For a protein with few
or no structures the balance reverses and the sweeps carry the answer.

`pubmed_search` takes an optional `email`. NCBI asks for a contact address so
it can warn before blocking heavy use, but it is **off by default** because it
sends the user's address to NCBI; only pass one if the user has agreed.

For papers whose numbers you need and whose abstract doesn't carry them, fetch
the full text once by DOI (or with a literature skill if one is loaded). If
that returns only an abstract — common for closed-access journals — say so and
move on rather than retrying alternate routes.

### 5. Verify every DOI against the claim, not just the DOI

```python
checked = ex.check_dois(ex.pdb_citation_claims(details))
```

A resolving DOI proves nothing about attribution: a real paper on the right
protein will pass while being credited to the wrong author or year. Pass
`claims` as `{doi: (first_author, year)}` — `pdb_citation_claims` builds it
from the PDB's own citation records, which are independent of CrossRef — and
check `claim_ok` on every record before anything goes in the report. Also
check `retracted`. Fix the report from CrossRef, never the other way round.

Cite inline as `[Author Year](https://doi.org/...)`. Unpublished depositions
("to be published") are evidence of a working construct but are not a
citation; say which entries are unpublished rather than attributing them to a
paper that happens to share the protein.

### 6. Report

Write the page, then print it:

```python
ex.write_html_report("PIM1_expression_strategy.html", "...", sections, sources)
ex.html_to_pdf("PIM1_expression_strategy.html")   # -> ..._strategy.pdf, or None
```

Hand over both files. Say so if `html_to_pdf` returned `None` rather than
quietly delivering one file.

Sections, in this order:

- **Recommendation.** Host, construct, tag, in bold in the first sentence,
  then the two or three load-bearing numbers. Say what the recommendation
  rests on — precedent counts, a yield, or a single report.
- **Molecular anatomy.** One `track_svg` track over the precursor: signal and
  propeptide segments, the mature chain, domains, glycan lollipops (filled =
  occupied, open = sequon only) and free cysteines. Follow it with
  `swatch_legend` and a line of mass/pI/ε280 and the numbering convention.
- **Construct maps.** `track_svg` with one track per route on a shared residue
  axis, tags and linkers drawn to scale, `hatch=True` on any segment whose
  length is not stated — never guess a tag length. Colour by host with
  `HOST_HUES`.
- **Host precedent.** A table of host family → entries → representative PDB
  IDs → tag status, then the conflict list from `summarise_hosts` as its own
  block if it is not empty.
- **Routes.** One row per route: host, construct and tags, length, reported
  yield, reference/PDB, verdict pill (recommended / back-up / rejected). State
  **"no yield reported"** explicitly where that is the case; an empty cell
  reads like an oversight.
- **Rejected hosts.** One short paragraph each saying what specifically is
  missing — no precedent, no yield, cost not justified by the molecule.
- **For the CRO brief.** Numbering convention, signal peptide, glycan and
  stability liabilities, host-specific sequence gotchas, licensing flag on any
  named engineered variant, and the questions still open from "Ask before
  assuming".
- **All PDB entries**, as `groups` with `collapse=True`, one group per host
  family, keyed `host:<family>` so the colours match the host-precedent table.
- **Caveats and closed-access papers**, ranked by what each would change.

Save both in the user's project folder (or as downloadable files on
claude.ai).

Collapsed groups are expanded by a `beforeprint` handler in the page, so the
PDF carries every row even where the HTML folds them away. CSS cannot reveal
a closed `<details>`, so do not remove that script.

## Evidence conventions

- **Never fill a gap.** "not reported" in a cell, "not stated in the retrieved
  text" in prose. Vector, signal peptide and tag position are the fields most
  often missing from methods sections.
- **Titre is not yield.** Keep "mg/L of broth", "mg purified per litre" and
  "mg per N-litre fermentation" visibly distinct everywhere, and reconcile
  every yield against its stated volumes. A figure given as a concentration
  times one litre is not a yield until the step volumes are stated.
- **Single-report routes.** When a soluble or high-yield route rests on one
  paper, say "rests on a single report" rather than presenting it as
  established.
- **Flag closed-access papers** and list them at the end with what each would
  settle, ranked by what would actually change the recommendation.
- **Sanity-check with `construct_props`.** Compare a paper's stated extinction
  coefficient with the one computed from the sequence, and predicted with
  measured mass. A gap of 131 Da is a lost initiator Met, which changes the
  construct length you report.
- Catalogue and supplier claims come from the vendor datasheet; potencies are
  comparable only against the same standard and cell assay.
- **Identify attached PDFs by content, not filename.** Generic names such as
  `main.pdf` can be any paper, and the same paper can arrive twice with
  different hashes.

## Known limitations

- `summarise_hosts` counts one host per PDB entry, from the first polymer
  entity. For a complex whose chains were made in different hosts, read the
  per-entity records in `pdb_construct_details` instead.
- `normalise_host` matches a fixed keyword list. An unrecognised host is
  returned unchanged rather than hidden, so scan the counts for raw strings
  that should have been collapsed and extend `HOST_FAMILIES` if needed.
- `construct_props`' pI is a Henderson-Hasselbalch estimate. It tracks
  ProteinAnalysis closely but is not a measurement, and it ignores the pKa
  shifts a folded structure imposes.
- A host recorded in the PDB is a host the protein was *made* in, not
  evidence of yield. Depositions tell you what folds; only methods sections
  tell you how much.
