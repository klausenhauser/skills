# crystallization-precedent-survey

Before you design a crystallization construct for a protein that others have
already crystallized, this skill collects what worked for them: construct
boundaries, crystallization conditions and whether tags were left on. It then
finds the wells in commercial screens that come closest to those conditions,
so your first plate tests known precedent instead of a blind screen.

Claude reads [`SKILL.md`](SKILL.md) for the full instructions. This page is the
human overview.

## What you give it

A UniProt accession and the domain you want to crystallize, for example
"human PIM1, P11309, kinase domain". Optionally, the screens you have on the
shelf.

## What you get back

- **Table A:** every PDB entry for that domain, with its construct range in
  canonical UniProt numbering, tag status, resolution, crystallization
  condition and citation.
- **Table B:** each recurring condition, the screen wells that best reproduce
  it, and how good the match is.
- **A report** as one HTML page, with construct families and condition groups
  colour-coded the same way in every section: which constructs and conditions dominate,
  which conditions couldn't be recovered and why, which screens (or custom
  grid) to set up first, and both tables. It opens in any browser and prints
  to PDF.

## The logic

1. **Machine-readable sources first.** UniProt, RCSB and the CSIRO C6 screen
   database answer most of the question without reading papers. Papers are
   only for gaps.
2. **Never invent a missing condition.** Many depositions record only
   "vapour diffusion" or a pH. Those are reported as gaps, with the paper
   where the real value may be.
3. **Construct ranges in canonical numbering.** Depositors often number
   against another isoform, so ranges come from aligning each deposited
   sequence to the canonical UniProt sequence, not from the PDB's residue
   numbers.
4. **Tags from the deposition record.** The mmCIF lists every residue that
   differs from UniProt. A short scar means the tag was cleaved; a long run
   (His tag, protease site) means the crystal grew with the tag on.
5. **Honest screen matching.** Conditions are matched component by component
   (salt, buffer, pH, precipitant). When no well is close, it says so and
   suggests a custom grid rather than forcing a weak match.

## Steps

| Step | What happens |
|---|---|
| 1. Target | Fetches the UniProt entry: sequence, annotated domains and every linked PDB entry. Confirms the domain range with you if it's ambiguous. |
| 2. Entries | Keeps the PDB entries that cover the domain. |
| 3. Conditions | Fetches each entry's crystallization details, resolution and citation from RCSB, and groups entries that share a condition. |
| 4. Constructs and tags | Maps every construct onto the canonical sequence, groups them into families, and checks one entry per family for cleaved or intact tags. |
| 5. Screens | Fetches formulations from CSIRO C6 (vendor PDFs only as a fallback). |
| 6. Matching | Finds and ranks the wells closest to each recurring condition. |
| 7. Report | Writes the two tables as CSV and the HTML report. |

## Scripts

[`scripts/survey.py`](scripts/survey.py) holds the helpers Claude calls at
each step. Python 3 standard library only; no editing needed.

| Function | Step |
|---|---|
| `fetch_uniprot_entry`, `select_pdb_ids_for_domain` | 1–2 |
| `fetch_rcsb_crystallization_batch`, `group_conditions` | 3 |
| `canonical_range`, `tag_status`, `construct_families`, `get_construct_tag_diffs` | 4 |
| `fetch_c6_screen`, `parse_screen_pdf_text` | 5 |
| `search_screens_for_terms`, `rank_screen_matches`, `normalise_chemistry` | 6 |
| `write_html_report` | 7 |

## Requirements

- Python 3 and network access to UniProt, RCSB and c6.csiro.au.
- `pypdfium2` only if a screen isn't in C6 and has to be read from a PDF.
- Works in Claude Code, and on claude.ai when code execution with network
  access is enabled. A protein with ~200 PDB entries takes a few minutes.
