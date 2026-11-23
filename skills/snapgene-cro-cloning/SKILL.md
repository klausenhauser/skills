---
name: snapgene-cro-cloning
description: Design, annotate and verify in-silico cloning constructs from SnapGene .dna vector files - replace an ORF with a new insert while preserving reading frame, choose seamless (homology-arm) vs restriction-site junctions, census restriction sites for single-cutter uniqueness across the whole plasmid (both strands, circular-aware), screen Dam/Dcm methylation blockage, set CDS feature qualifiers (accession-linked /product, /note, residue numbering, cached translation), and emit a CRO synthesis order sheet as .xlsx. Reads and writes .dna directly with the Python standard library, so it needs neither SnapGene nor any pip package. Use when swapping an ORF in a SnapGene plasmid, preparing a gene-synthesis or cloning order for a contract research organization (CRO), annotating or QC-ing a .dna construct, or when SnapGene displays restriction sites at the wrong position after an edit.
---

# SnapGene ORF replacement and CRO order sheet

Take a SnapGene vector, cut out one ORF, splice in a new one while preserving
the reading frame of any upstream fusion partner, annotate the result properly,
verify it, and produce the paperwork a CRO needs to synthesise and clone it.

## Requirements

Python 3 standard library only. No Biopython, no openpyxl, no autosnapgene,
no snapgene-reader, and no SnapGene installation — the .dna reader and writer
are implemented against the file format itself.

```python
import sys; sys.path.insert(0, "<path to this skill>/scripts")
import snapgene as sg
```

## Ask before assuming

This skill has no house style of its own. The conventions below are defaults
to offer, not rules to apply silently. At the first checkpoint, ask for:

- **The vector.** Ask which plasmid the construct goes into and have them
  supply the `.dna` file. A host choice does not imply a vector. Do not parse a
  vector, pick cloning sites or draft a sequence before they have answered.
- **The CRO order-sheet template.** Ask for their own sheet as a template file.
  Read its headers with `read_order_sheet_columns()` so their columns drive the
  output. The fallback column list below is a description of one common
  template, not a substitute for theirs.
- **The name of every final plasmid — ask for each one, every time.**
  Construct names usually belong to a numbered series that continues from
  existing constructs, and only the person ordering knows which numbers are
  already taken. Vector filenames, sibling constructs and a visible series
  prefix are not an answer: they show the pattern, not the next free number.
  Ask for the name of each construct explicitly and do not fill one in as a
  placeholder.
- **The signal peptide, tags and linkers**, for any secreted construct. Ask
  whether to use the native signal peptide or a heterologous one, and which.
  Offer options with a recommendation; they decide.
- **Where and how to clone — ask for each construct, every time.** This covers
  the strategy (restriction-site or seamless, and the cost difference at their
  CRO), the specific enzyme pair or junction points, and how much is
  synthesised: the ORF alone, or a whole cassette including tags and linkers.
  Report the single-cutter census and a recommendation, then let them choose.
  The vector's own annotations are evidence, never the answer — a feature
  marking a previously synthesised region, a site named in the filename, or
  the pair a sibling construct used all describe how the backbone was built
  last time, not what this order should do. Present what you found and ask.

An answer covers the construct it was given for. Two constructs in one batch
get two questions, and a convention confirmed for an earlier order is a default
to re-offer, not consent to reuse.

Take every tag, linker and signal-peptide sequence from a source record — a
UniProt feature range, a vendor sequence, or a file they supply — and cite the
accession and residue range. Never type one from recall. After building,
compare the translated leader in the construct against the record.

## Stages and human checkpoints

Six stages. **Every stage ends with a checkpoint: stop, report, and wait for an
explicit reply before starting the next.** Silence, or approval of something
else, is not approval. If they say stop or ask for a change, do not continue
past that point.

| Stage | Steps | Ends with |
|---|---|---|
| Brief | host, tags, construct boundaries (carried in by the human) | Checkpoint 0 |
| Verify | 1 parse vector and ORF, 2 parse inserts | Checkpoints 1, 2 |
| Build | 3 draft ORF replacement | Checkpoint 3 |
| Design | 4 junction strategy, then the final fragment | Checkpoint 4 |
| Check | 5 annotate, write, verify | Checkpoint 5 |
| Package | 6 CRO order sheet | Checkpoint 6 |

Step 3 builds a draft because the site census in step 4 must run against the
assembled construct. The fragment actually submitted for synthesis is built
inside step 4, once the strategy is approved.

Checkpoint 4 is where the cloning strategy, the enzyme pair or junction points
and the extent of the synthesised fragment are decided, and Checkpoint 5 is
where each final plasmid is named. Neither is answerable from the vector, so
neither may be inferred — see **Ask before assuming**. Where a run has been
told to report everything at the end instead of stopping at each checkpoint,
these two still have to be asked up front, because every later step is built on
them and a wrong guess invalidates the construct, the order sheet and the name
on both.

Each checkpoint reports, in this shape: **Done** (2–4 lines), **Findings**
(measurements, and anything unexpected — report only what a check actually
returned, and say "not assessed" for anything not run), **Decisions needed**
(open choices, with a recommendation), **Ask** ("Continue to <next stage>?
yes / change <item> / stop").

## Workflow

### 1. Verify — parse the vector, locate the ORF

```python
vec = sg.read_dna(vector_path)        # Dna(5784 bp, circular, 23 features)
orf = vec.feature("GeneName", ftype="CDS")
```

Feature coordinates come back 0-based half-open and slice the sequence
directly: `vec.seq[orf.start:orf.end]` is the annotated CDS. Do not subtract 1
to "correct" for 1-based indexing — the conversion is already done, and doing
it twice shifts the frame and yields a plausible-looking but wrong translation.

Confirm frame integrity rather than trusting the annotation:

```python
sg.translate_and_validate(sg.feature_span(vec, orf))
```

Check that the length is a multiple of 3, that the upstream fusion partner
reads in frame right up to `orf.start`, and where the real first stop codon
falls — a downstream stop cassette or C-terminal tag can sit a few bases past
the annotated feature end, so translate a wide window past `orf.end` instead of
trusting the boundary.

**Checkpoint 1.** Report the vector name and length, the ORF coordinates and
length, whether it is a multiple of 3, the upstream fusion partner and whether
it reads in frame, and where the first real stop codon falls. Ask whether this
is the right ORF.

### 2. Verify — parse the insert

```python
pairs = sg.parse_protein_dna_fasta_pairs(path)   # {name: {'protein', 'dna'}}
sg.translate_and_validate(dna, expected_protein=protein)
```

Always prefer a supplied DNA sequence over back-translating from protein, and
validate that it translates cleanly to the stated protein before using it.

**Protein supplied but no DNA.** This skill has no codon optimiser and no
reverse translator. Do not back-translate silently. Stop and offer exactly two
options: the human provides codon-optimised DNA (from the CRO, from a codon
optimisation tool, or from an existing construct), or — if browser tools are
available in this session — open the CRO's codon-optimisation page, paste the
protein, select the host they named and list any restriction sites to avoid.
Many such tools email the result, so if the page does not show the DNA, ask
them to paste it back. Do not retry a failing browser action more than two or
three times.

Only with explicit approval may a self-made back-translation be used, and only
as a clearly labelled **placeholder** — never as the order sequence — so
junctions and sites can be checked while the real optimised DNA is awaited.

**Checkpoint 2.** Report each insert's source, length, `translate_and_validate`
result and the host it was optimised for.

### 3. Build — draft the replacement

```python
new_seq = sg.replace_orf(vec.seq, orf.start, orf.end, insert_dna)
delta = len(insert_dna) - orf.length(len(vec.seq))
```

To carry the annotation over, copy the features and shift everything
downstream of the edit:

```python
import copy
con = sg.Dna(seq=new_seq, circular=vec.circular, dam=vec.dam, dcm=vec.dcm,
             features=copy.deepcopy(vec.features), notes=dict(vec.notes))
spanning = sg.shift_features(con, orf.end, delta)
```

`shift_features` returns the segments that *span* the edit point instead of
guessing at them — they cannot be shifted correctly without knowing what the
edit meant, so resolve each one deliberately. Then replace the target
feature's own segment with the insert's extent.

**If the human has already made the edit in SnapGene and supplied a new
`.dna`,** treat that file as authoritative: extract the junction sequences from
it, do not regenerate or overwrite it.

**Checkpoint 3.** Report the draft length, the change against the original ORF,
and that the draft is an input to the census, not a deliverable.

### 4. Design — seamless or restriction sites

Apply their stated preference. **Never mix the two on one construct.** Each
construct gets one strategy end to end — do not cut one side and
homology-arm the other. If the topology seems to force a hybrid (an N-terminal
in-frame fusion with no site at that junction, a real site at the other), that
is a signal to extend the *restriction* strategy back to a site further
upstream, not to call the near junction seamless. In a vector where the insert
fuses in frame right after a leader with no site at that exact junction,
extending the fragment back to a unique site upstream of the translation start
lets the whole insert ship as one restriction-cloned fragment.

Census candidate enzymes **across the full assembled construct**, not the
insert. The requirement is a single cutter in the whole plasmid: an insert that
is internally clean still breaks the strategy if the backbone carries a second
copy.

```python
sg.check_unique_and_clean(insert_dna, {"EcoRI": "GAATTC"})      # insert screen
sg.check_sites_in_construct(new_seq, {"NdeI": "CATATG", "XhoI": "CTCGAG"})
# -> {name: {'positions': [...], 'count': n, 'unique': bool}}; want unique=True
```

`find_sites` takes a literal IUPAC recognition sequence, **not a regex**, and
by default scans **both strands** and **wraps the circular junction**. Both
matter: every Type IIS enzyme (BsaI `GGTCTC`, BsmBI `CGTCTC`, BbsI `GAAGAC`) is
non-palindromic, so a forward-only scan reports a contaminated insert as clean;
and a site straddling the origin is invisible to a linear scan. Palindromic
sites de-duplicate, so nothing double-counts. Pass `circular=False` only for a
linear fragment. `sg.ENZYME_SITES` ships a panel of common enzymes to census
against.

**When the intended enzyme is blocked** — by an internal copy, or by
methylation — work down this ladder and stop at the first rung that works:

1. **Silent synonymous codon change** removing the site from the insert.
   Prefer a codon that is also host-preferred. `remove_site_synonymously()`
   proposes the changes rather than applying them; confirm the protein is
   unchanged with `translate_and_validate`, then re-census the **full**
   construct.
2. **Extend the synthesis boundary** a few hundred bases further into the
   backbone to reach an alternate site already unique there. Costs extra
   synthesis; leaves the insert untouched.
3. **Seamless / homology-arm cloning** — last resort, only when no synonymous
   fix and no alternate site exist.

Record which rung was used and why in the construct's `/note` and the
file-level description: a reviewer cannot otherwise tell a deliberate silent
mutation from an error.

Screen methylation mechanically, not from memory:

```python
sg.screen_methylation_context(new_seq, site_position, site_len=6, enzyme="XbaI")
# -> {'dam': [...], 'dcm': [...], 'context': str, 'blocked_candidate': bool}
```

Dam (`GATC`) and Dcm (`CCWGG`) methylation silently blocks a site that looks
perfectly unique. Treat a hit as a candidate to confirm against REBASE, not a
verdict — whether overlap actually blocks a given enzyme is enzyme-specific. If
a blocked site is the only option, shift the junction a few bases rather than
reaching for a methylation-deficient strain.

Junctions:

```python
five, three = sg.seamless_junction_arms(vec.seq, orf.start, orf.end,
                                        insert_dna, arm_len=18)
fragment = sg.restriction_synthesis_fragment(vec.seq, five_cut_start, orf.start,
                                             orf.end, insert_dna, three_cut_end)
```

The synthesis fragment **must include the flanking restriction sites
themselves**. Submitting the coding sequence alone and leaving the enzyme sites
out of the synthesised DNA is the standard mistake — they have to be physically
present for the cut to work.

**Checkpoint 4.** Report the strategy (restriction or seamless, never mixed),
each enzyme's whole-plasmid census, any Dam/Dcm hits, which rung of the ladder
was used and why, and the fragment boundaries. A silent codon change or an
extended boundary needs explicit approval here.

### 5. Check — annotate, write, verify

```python
sg.apply_feature_metadata(feat, con.seq,
    description="<protein description>", accession="<ACCESSION>", db="UNIPROT",
    gene="<GENE>", note="Mature chain (residues N-M) of UniProt <ACCESSION>",
    numbering_start=N)
sg.set_construct_description(con, "<html><body>…this construct…</body></html>")
sg.write_dna(con, out_path)
```

Then verify the file **as written to disk**, and read the result:

```python
report = sg.verify_construct(sg.read_dna(out_path), insert_feature="<name>",
                             expected_protein=mature_protein,
                             expect_unique=["NdeI", "XhoI"])
```

`verify_construct` measures frame integrity, internal stops, the translation
against `expected_protein`, the cached `translation` qualifier against the
actual sequence, the numbering attribute, whether `/product` carries an
accession link, a non-empty file-level description, and single-cutter status
for each named enzyme. It returns measurements and a `problems` list — it
asserts nothing on your behalf.

**Report only what it returns.** Do not describe a construct as verified for a
property the report does not cover, and never state the outcome of a check that
was not run — say "not assessed". Sequence length and feature ranges being
correct tells you nothing about the qualifier-level defects, which is exactly
how an empty `translation` ships undetected.

`sg.compute_protein_mw(protein)` cross-checks mass against a value quoted in
the source record. `sg.check_construct_name(name)` checks the name against the
length cap.

**Checkpoint 5.** Show the report as returned: `all_ok`, the `problems` list,
the measurements. If there are problems, propose fixes and wait — do not fix
and continue on your own. Say explicitly that opening the file in SnapGene is
the one check this skill cannot perform, and ask them to confirm it opens and
draws correctly.

### 6. Package — the CRO order sheet

```python
columns = sg.read_order_sheet_columns(template_path)   # their own headers
sg.write_order_sheet(rows, "CRO_order_sheet.xlsx", columns=columns)
```

Map each value to the template's own headers and keep its column order. Do not
invent columns. `write_order_sheet` reports any row key that matches no column
in `unknown_keys` rather than dropping it silently; if a needed value has no
column, stop and ask. It raises rather than writing a file a spreadsheet
application will refuse if a cell exceeds the 32,767-character limit — a long
vector sequence is the usual cause.

Only when they say there is no template, fall back to: `No`, `Final plasmid
Name`, `Final plasmid sequence`, `Insert sequence`, `Insert Source`, `Cloning
sites (5' and 3')`, `5'homologous sequence (~30 bp)`, `3' hom (~30 bp)`,
`Vector`, `Vector source`, `Vector sequence`, `Shippable`, `Maxiprep`.

Filling that sheet:

- **Insert Source** is `"to be synthesized"` once the insert is being ordered
  as a synthesis job, not a database accession.
- **Vector** is the plain vector name; a CRO tracking code belongs in
  **Vector source**.
- **Cloning sites and the homology columns are mutually exclusive.** For
  restriction cloning, name both enzymes in `Cloning sites` and leave the
  homology columns blank. For seamless, leave `Cloning sites` **entirely
  blank** and put the real junction sequence in the homology columns — never a
  restriction-site label inherited from the vector's own name, which describes
  how the backbone was built, not where the insert joins.

**Checkpoint 6.** Show the rows — name with its length against the cap, cloning
sites or homology arms, vector, insert source — and ask them to confirm, and to
open the sheet, before it is treated as final. This skill never submits an
order; sending it is the human's step.

## Conventions worth offering as defaults

Offer these; do not apply them without asking.

- **Naming**: `<SERIES><n>_<vector family>_<tags/cleavage sites>_<ORF
  name>[-FL or a residue range]`, with a length cap — 50 characters is common.
  Use the insert's original source-file identifier verbatim as the ORF-name
  component rather than a prettified form. C-terminal tags go at the end of the
  name when they are C-terminal in the construct, not grouped with the
  N-terminal ones. Check the cap with `check_construct_name()`; it is easy to
  breach once tags and a range suffix are appended.
- **Feature annotation**: keep the insert's original dataset or file ID as the
  CDS feature name rather than renaming to a bare gene symbol, and populate
  `/gene`, `/product` and `/note` from the source record. When the insert is a
  truncated or mature chain relative to the full reference, the `/note` must
  cite both the accession **and** the residue range it corresponds to.
- **`/product` must carry the accession visibly.** A plain-text `/product` does
  not surface a visible database link in SnapGene's feature inspector. Build it
  with `make_product_html(description, accession, db=...)`, in addition to —
  not instead of — a descriptive `/note`.
- **Replace the file-level description every time.** See the pitfall below.
- **Species shorthand**: the common short form in construct names and metadata
  prose (`Human`, `Mouse`), with the full binomial reserved for the
  provenance note citing the source database.

## Pitfalls

- **The cached `translation` qualifier goes stale and nothing warns you.**
  SnapGene caches each CDS feature's peptide as a plain string. An ORF
  replacement rewrites the sequence and shifts the feature range but does not
  refresh that string, so a construct can ship with an empty or wrong
  translation while its sequence, feature range, reading frame and enzyme
  positions are all correct. Always set it with `apply_feature_metadata` and
  confirm with `verify_construct`.
- **Two translation conventions are easy to get wrong**, and both were
  established by reading real files. A terminal stop codon is **kept** as a
  trailing `*`; stripping it makes every stop-terminated CDS look stale. And a
  **multi-segment feature is translated per segment and joined with a comma**,
  so a two-segment CDS reads `MSIQ…VFA,HPETL…W*`. Use
  `feature_translation()` for the stored form and `protein_of()` for a plain
  protein to compare against a database record.
- **A feature can span the origin** of a circular plasmid, in which case
  `end <= start` and the span is `seq[start:] + seq[:end]`, not the empty slice
  Python gives. Use `segment.subseq(seq)` and `feature.length(len(seq))` rather
  than slicing or subtracting by hand; `spans_origin` tells you which case you
  are in.
- **Don't leave the vector's legacy file-level description.** Every vector
  carries a description in its notes block left over from whatever it was last
  used for, and an ORF replacement copies it through verbatim — it is not a
  feature, so nothing about editing the ORF touches it, and it is wrong for the
  new construct. Replace it every time with
  `set_construct_description(dna, html)`, not only when it is reported stale:
  `verify_construct` can see that it is empty but cannot tell a stale one from
  a correct one.
- **`/product`, `/note` and `/gene` are HTML source, not plain text.** SnapGene
  stores them wrapped as `<html><body>…</body></html>` and renders the markup
  in the feature inspector. Pass real markup — `make_product_html()` builds it —
  and let the writer escape it; pre-escaping it yourself, or escaping it twice,
  makes SnapGene display literal `<html>` tags. A read-write-read test cannot
  catch that, because a matching error in the reader cancels it. **Open the
  written file and look at the inspector** — a qualifier that renders as tags
  is the only symptom.
- **Don't reuse the vector's own flanking sites** — the ones often embedded in
  its name — as the cloning sites for a seamless junction. They describe how
  the backbone was built, not where the insert joins.
- **Don't forget the restriction sites in the synthesised fragment** when the
  strategy is restriction cloning.
- **After any session reset, re-derive coordinates and insert data from the
  source files** rather than assuming prior in-memory state. Files on disk
  persist; variables do not.
- **When naming sources conflict, the human's stated convention wins.** Prefer,
  in order: an explicit instruction in the current conversation, then a stated
  convention, then a correction already applied to a delivered file, then older
  plan text. Reverting a name they already corrected is worse than leaving it
  alone; if the sources genuinely conflict, ask.

## Known limitations

- **Opening the written file in SnapGene is not verifiable from here.** Every
  programmatic check can pass on a file the application still refuses. Ask the
  human to open it, at Checkpoint 5, every time.
- **The writer keeps only what it models**: sequence, topology and the
  strandedness and methylation flags, all features with their segments,
  colours, directionality, translation and numbering attributes and qualifiers,
  and the file-level notes. It drops the cached restriction-site list, the
  position overlay, edit history, primers, alignments and enzyme visibility.
  Of those only **primers** are user data — if the vector carries primers that
  matter, say so and keep the original file alongside, or re-add them in
  SnapGene. That trade is deliberate: see `sg.WRITER_NOTE`.
- No codon optimiser, no reverse translator, no cut-site offsets or overhangs
  (recognition sequences only), no star activity, no isoschizomer expansion,
  and no Gibson/Golden-Gate assembly simulation.
- `read_order_sheet_columns` reads a template's headers, not its formatting.
  The sheet is written fresh with those headers rather than edited in place.
