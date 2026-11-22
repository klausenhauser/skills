---
name: crystallization-precedent-survey
description: Survey every PDB entry that has crystallized a given protein domain, extract the construct boundaries and crystallization conditions each used, check whether affinity/purification tags were cleaved before crystallization, and cross-reference the recovered conditions against commercial screen formulations to prioritize which screens to test first. Use this before designing a new crystallization construct for a protein with prior crystallographic precedent (given a UniProt accession and a domain of interest), or when planning initial screening for a target that has been crystallized before.
---

# Crystallization precedent survey

Before designing a new construct to crystallize a domain that has already been
solved by others, mine the PDB depositions for that protein: what construct
boundaries worked, what buffer/salt/precipitant conditions produced
diffracting crystals, and whether the deposited chains still carry affinity
tags. Then map the recovered chemistry onto commercial sparse-matrix/grid
screens so the first plate poured tests the conditions closest to known
precedent, not a blind screen.

This is a data-mining task, not a literature-synthesis task — the answer
should be a table with citations, not prose. Most of the work below runs
directly against machine-readable sources (UniProt REST, RCSB REST/mmCIF,
CSIRO C6) with no LLM judgment calls. Go to the primary papers only to fill
gaps (step 3).

## Requirements

- Python 3, standard library only. The helpers are in `scripts/survey.py`
  (bundled with this skill); import it with
  `sys.path.insert(0, "<this skill>/scripts"); import survey`.
- Network access to rest.uniprot.org, data.rcsb.org, files.rcsb.org and
  c6.csiro.au. If this environment can't reach them (e.g. claude.ai without
  network access), say so and stop rather than answering from memory.
- Only for the PDF fallback in step 5: `pypdfium2` (`pip install pypdfium2`).

## Workflow

1. **Resolve the target and domain.** Call `fetch_uniprot_entry(accession)` to
   get the sequence, UniProt-annotated `domains` (PH, START, kinase, etc. with
   residue ranges), and the raw `pdb_xrefs` list (every PDB entry UniProt links
   to, with per-chain residue ranges). State the accession, organism and
   canonical length in the report. Confirm the domain of interest's range
   with the user if more than one domain-like feature is annotated near it —
   UniProt's `Domain`/`Region` boundaries are a starting point, not gospel;
   real constructs often extend a bit past them (see step 4).

2. **Select the relevant PDB entries.** Call
   `select_pdb_ids_for_domain(pdb_xrefs, domain_start, domain_end)` to filter
   to entries whose deposited chain overlaps the domain by at least 50% of its
   length (tune `min_overlap_frac` down if the domain is small and constructs
   are much larger, e.g. multi-domain complexes). Manually exclude any hits
   that are clearly a different domain of the same protein (e.g. an NMR
   structure of an N-terminal PH domain when you want the C-terminal START
   domain) — check the returned `chains` ranges against your target range
   before trusting the filter blindly.

3. **Pull crystallization conditions for every hit.**
   Call `fetch_rcsb_crystallization_batch(pdb_ids)` — two RCSB REST calls per
   ID (about 1 s per entry; run it in the background for more than ~50
   entries), returning resolution, pH, temperature, the free-text
   `pdbx_details` (reservoir composition, when deposited), the primary
   citation title/year and the deposited entity sequence. Then call
   `group_conditions(batch)`, which groups entries by condition ignoring case
   and whitespace (many depositions from the same paper reuse one condition),
   rather than treating each PDB ID as an independent condition. Read the
   groups by eye and merge the ones that are the same chemistry written
   differently. **Expect gaps**: some depositions — especially older ones
   (pre-~2010) — record only `"VAPOR DIFFUSION, HANGING DROP"`, only a pH, or
   only the protein buffer with no salt/buffer identity. Don't fabricate a
   plausible-sounding condition to fill the gap; report it as unrecovered and
   name the paper where the real value might live (methods or SI), then try
   fetching that paper's full text once (via its DOI, or a literature skill
   if one is available) — if that returns only an abstract (common for
   closed-access journals), say so and move on rather than retrying
   alternate fetch routes.

4. **Construct boundaries and tags.** Deposited residue numbers can't be
   trusted: depositors often number against another isoform while labelling
   the canonical accession (human PIM1, for example, is mostly numbered
   against its 404-residue long form while the canonical sequence is 313).
   So take construct boundaries from `canonical_range(entity_seq,
   canonical_sequence)` for every entry, which maps the deposited sequence
   onto the canonical one and returns any extra residues at each end, and
   `tag_status(...)` on that result for each entry's tag call (intact / scar /
   none). Group entries into construct families with
   `construct_families({pdb_id: (start, end)})`, which merges ranges within 3
   residues of each other (usually 2-6 families across a protein's full
   deposition history), and report all ranges in canonical numbering. If a family maps with an internal gap or not at all, check the
   UniProt isoforms — a construct may have been solved from a natural splice
   isoform lacking a loop the canonical sequence has; report this so the
   user's construct (presumably canonical) isn't surprised by an extra
   flexible loop none of the precedent structures had to deal with.

   Then call `get_construct_tag_diffs(representative_pdb_id)` on one entry per
   family. This parses `_struct_ref_seq`/`_struct_ref_seq_dif` from the mmCIF —
   the authoritative record of every residue where the deposited chain
   diverges from the UniProt reference — and separates real point differences
   (`'conflict'`, `'engineered mutation'`) from cloning/tag scars
   (`'expression tag'`), and confirms or corrects the `tag_status` calls for
   that family. A **2-6 residue** N- or C-terminal scar flagged
   `expression tag` or `conflict` means a tag WAS present and was cleaved
   before the crystal was grown (a Gly-Pro doublet is the classic PreScission/
   HRV-3C signature; Ser-Asn-Ser or similar short runs are generic
   cloning/thrombin scars). An intact tag shows up as 8+ extra residues
   matching a known tag (His6/His8, often with an uncleaved protease site,
   GST, MBP, SUMO). Intact tags are common — report each family's tag status
   explicitly, since crystals grown with the tag on tell the user that
   cleaving it is optional for this protein.

5. **Get the screen formulations.** Use CSIRO's C6 database: find each
   screen's numeric ID by searching its name on https://c6.csiro.au/screens
   (the ID is in the screen page's URL; JCSG+ is 867), then call
   `fetch_c6_screen(screen_id)`. It returns one row per well with ingredient,
   role, concentration and pH, in plate order. Ask the user which screens
   they have on the shelf; if they don't say, survey a few broad ones (e.g.
   JCSG+, PACT premier, Morpheus, an Index-type screen) and say which.
   Only if a screen isn't in C6, fall back to the vendor's formulation PDF:
   find it by web search (the vendor's URLs change, so don't hardcode one),
   extract its text with `pypdfium2` (`pdfium.PdfDocument(path)` →
   `page.get_textpage().get_text_range()` per page), then parse with
   `parse_screen_pdf_text(text, screen_name)`. This handles four layouts —
   if the result is empty or far short of the screen's well count, print a
   slice of the raw text and extend the regex rather than guessing the table
   structure blind, and spot-check a handful of parsed rows against the PDF.

6. **Cross-reference and rank.** For each distinct recovered literature
   condition, call `search_screens_for_terms(all_rows, terms, mode="all",
   ph=<condition pH>)` with the salt/buffer/precipitant names from that
   condition (e.g. `["sodium citrate", "peg 3350"]`), spelled the way the
   screen spells them — matching is by name, not chemistry, so check the
   screen's ingredient names for salts with several names (e.g. ammonium
   phosphate is "di-ammonium hydrogen phosphate" in JCSG+). Terms and wells are
   normalised (case, PEG spellings such as "PEG3350" vs "Polyethylene glycol
   3,350") and matched as whole words, so "tris" doesn't match "bis-tris";
   the pH filter applies ±0.5 to wells with a known pH. For near matches,
   call `rank_screen_matches(all_rows, terms, ph=<pH>)`, which scores every
   well by how many of the condition's components it reproduces plus one for
   the pH, and lists what each well is missing. Judge the ranking by the
   literature condition's components (salt identity, buffer identity, pH,
   precipitant identity, precipitant concentration) — an exact salt+buffer+pH+precipitant match is "test
   first"; a same-chemistry-family match with a different counter-ion or PEG
   size is a reasonable second choice; note explicitly when NO screen well is
   a good match (e.g. a detergent additive like octyl-glucoside rarely appears
   in standard screens, and high-salt conditions such as 1 M sodium acetate
   are often absent from PEG-based sparse matrices) rather than forcing a
   weak match into the top recommendation. When the most-used precedent
   conditions have no good match, recommend a custom grid around them.

7. **Report.** Produce two tables and a short recommendation: (a) PDB ID →
   canonical construct range → tag status → resolution → condition → citation
   (one row per PDB ID, grouped visually by shared condition), and (b)
   literature condition → best-matching screen wells → match quality. State
   plainly which conditions could not be recovered and why (paywalled paper,
   undeposited header) rather than omitting them silently — a gap is still
   useful information for construct/screening planning. Save both tables as
   CSV, then build the report with `write_html_report(path, title, sections,
   sources)`: one self-contained HTML page that opens in any browser and
   prints to PDF. Sections, in order:
   - **Target and domain:** accession, organism, canonical length, domain
     range, number of entries, the numbering caveat if it applied.
   - **Construct families:** one row per family from `construct_families`
     (`_swatch: "fam:<label>"`), with the entry count and a tag breakdown as
     pills (e.g. `["intact tag 20", "scar 30", "none 13"]`).
   - **Screen hits:** `groups`, one per condition group you merged in step 3
     (`key: "cond:<n>"`, label = the condition, note = its entry count),
     largest first; rows = the wells from `rank_screen_matches` with one ✓/✗
     column per component plus pH, the score and a quality pill (exact /
     partial / no match).
   - **Recommendation** and **unrecovered conditions**.
   - **All PDB entries:** `groups` with `collapse=True`, one per condition
     group with the same `cond:<n>` keys (so colours match the screen hits),
     the condition written once as the label and left out of the rows; then
     one block for single-entry conditions (keeping a condition column) and
     one for unrecorded conditions. Within a group, sort by family (largest
     first) then resolution; each row gets `_swatch: "fam:<label>"` on the construct
     column and a tag pill.
   - **Sources:** UniProt, RCSB, C6 and paper links.
   Save everything in the user's project folder (or as downloadable files on
   claude.ai).

## Known limitations

- RCSB's `pdbx_details` field is free text with no controlled vocabulary —
  don't try to auto-parse salt/buffer/precipitant identities out of it
  beyond simple keyword search; read it and transcribe by eye for the final
  table.
- `get_construct_tag_diffs`'s regex parsing of raw mmCIF text is a
  lighter-weight alternative to the RCSB `polymer_entity`/`entry` JSON APIs
  (which don't expose `struct_ref_seq_dif` at all) — if a future mmCIF
  formatting change breaks the regex, fall back to downloading the file and
  inspecting the `_struct_ref_seq_dif` loop by eye.
- C6 transcribes vendor formulations; for a final plate, check the wells you
  pick against the vendor's current sheet, since screens are revised.
