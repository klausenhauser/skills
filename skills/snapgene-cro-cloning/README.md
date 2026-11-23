# snapgene-cro-cloning

Swap an ORF in a SnapGene vector, pick the cloning strategy, annotate the
construct properly, verify it, and produce the order sheet a CRO needs.

## What it does

Six stages, each ending at a human checkpoint:

1. **Verify the vector** — parse the `.dna`, locate the ORF, confirm the frame
   and find where the real first stop codon falls.
2. **Verify the insert** — parse protein/DNA pairs and check the DNA translates
   to the stated protein.
3. **Build a draft** — splice the insert in and shift the downstream features.
4. **Design the junctions** — restriction sites or seamless assembly, with a
   whole-plasmid single-cutter census and a methylation screen.
5. **Annotate and verify** — set the CDS qualifiers, write the file, then
   re-read it from disk and measure it.
6. **Package** — fill the CRO order sheet as `.xlsx`.

The deliverables are a `.dna` construct file and an `.xlsx` order sheet.

## Why it is built this way

### It reads and writes `.dna` itself

A `.dna` file is a flat sequence of typed blocks: one byte of type, four bytes
of big-endian length, then the payload. The sequence, the features and the
notes are three of those blocks, and the features and notes are plain XML. So
the skill parses and writes the format directly with `struct` and
`xml.etree`, and needs neither SnapGene nor any third-party package.

That is not only about dependencies. It fixes a real bug class.

SnapGene writes several blocks that cache absolute nucleotide positions
computed against the sequence as it stood when the file was saved — a packed
position overlay and a chain of edit-history blocks. A library that loads a
file, edits the sequence and writes it back copies every block it does not
model through as raw bytes, stale caches included. The result is a file whose
sequence, features and translations are all correct while SnapGene desktop
still draws a restriction site tens to hundreds of bases from its true
position, because it reads one of those caches when rendering.

This writer emits only the blocks it models, so there is no stale state to
strip, and no "rebuild clean before writing" step to remember. The size
difference is the clearest evidence of how much of a working file is cache: on
the vectors this was tested against, a 356 KB `.dna` rewrites to 19 KB with
every base, feature, qualifier and note preserved.

The trade is explicit. Primers, alignments, the cached restriction-site list,
edit history and enzyme visibility are dropped. Of those only primers are user
data, so the skill says to keep the original file alongside when the vector
carries primers that matter.

### What testing against real files changed

The module was checked against twenty-two real vectors, which corrected seven
things a reasonable first implementation gets wrong:

- **The cached translation keeps its terminal stop.** SnapGene stores a
  stop-terminated CDS with a trailing `*`. Stripping it makes every such
  feature look like a stale cache — 44 of the 135 cached translations in the
  test set.
- **A multi-segment feature is translated per segment, joined with a comma.**
  A two-segment CDS is stored as `MSIQ…VFA,HPETL…W*`. Translating the spliced
  span as one string gives the same residues without the separator, and
  comparing that against the stored value reports a false mismatch.
- **A feature can span the origin.** On a circular plasmid a segment may have
  `end <= start`, and its span is `seq[start:] + seq[:end]` — not the empty
  slice Python returns. One real vector has such a feature, and a naive reader
  raises or silently annotates nothing.
- **The numbering attribute is `consecutiveNumberingStartsFrom`.** The
  companion attribute `consecutiveTranslationNumbering="1"` sits on every
  feature in every file — all 553 of them — so it is a default, not the switch
  that turns custom numbering on.
- **A qualifier is HTML escaped exactly once.** SnapGene stores `/product`,
  `/note` and `/gene` as HTML source inside an XML attribute, so the file holds
  `&lt;html&gt;&lt;body&gt;…`. Escaping on write while also unescaping on read
  is self-consistent and round-trips perfectly, but it is one escape too many:
  SnapGene decodes once, so the feature inspector shows literal `<html><body>`
  tags instead of formatted text. A read-write-read test cannot see this,
  because the extra unescape on the way in cancels the extra escape on the way
  out. Only comparing the written bytes against SnapGene's own catches it.
- **Line breaks in an attribute must be character references.** An XML parser
  applies attribute-value normalisation and silently rewrites a literal newline
  inside an attribute to a space, so a feature label containing a line break —
  SnapGene writes those — loses one character per break on re-read, invisibly.
  `quoteattr()` handles this; hand-rolled quoting has to do it explicitly.
- **Not every file has a notes block.** Defaulting a `Type` into an empty one
  looks harmless and is how a new construct gets written, but applied
  unconditionally it invents metadata the source never carried, so a file
  without notes does not come back the way it went in.

With those seven fixed, all 135 translation qualifiers across the twenty-two
files reproduce exactly, all twenty-two round-trip through the writer with no
loss of sequence, topology, features, segments, colours, qualifiers, numbering
or notes, and all 647 text qualifiers are written byte-for-byte as SnapGene
itself wrote them.

The HTML escaping bug was found only because a construct was opened in SnapGene
and the annotation looked wrong. That is the argument for the final human
checkpoint: every automated check the module ran on itself passed while the
file was being written in a form SnapGene would not render. A round-trip test
proves the writer agrees with the reader, not that either agrees with
SnapGene — an escape too many on the way out and a decode too many on the way
in cancel exactly. `scripts/selftest.py` therefore compares the written bytes
against a SnapGene-authored file directly, which is also what turned up the
last two findings.

### The checkpoints are the point

Cloning design is cheap to redo on screen and expensive to redo after
synthesis. Every stage stops and reports rather than running to a finished
answer, because the decisions that matter — which vector, which series, which
signal peptide, whether a silent codon change is acceptable — belong to the
person ordering the DNA. The skill asks instead of inferring, and it treats
approval of one thing as approval of that thing only.

### Verification reports, it does not assert

`verify_construct` returns measurements and a `problems` list rather than
raising. This matters because the failure that motivated it was not a wrong
sequence: it was a construct whose sequence, feature range, reading frame and
enzyme positions were all correct, and whose cached translation qualifier was
empty. Every sequence-level check passed. So the report covers the
qualifier-level properties too, and the skill is explicit that a property the
report does not cover must be described as "not assessed" rather than verified.

The one check it cannot make is whether SnapGene opens the file. That is a
human checkpoint, every time.

## Layout

```
SKILL.md              the workflow, checkpoints and pitfalls
scripts/snapgene.py   .dna reader/writer, sequence and site tools, xlsx writer
scripts/selftest.py   checks the reader/writer against real .dna files
```

Run the self-test against your own vectors after changing anything in
`snapgene.py` — no test data ships with it, and the more varied the folder the
more it covers:

```
python3 scripts/selftest.py ~/path/to/dna/files
```

`scripts/snapgene.py` is importable on its own:

```python
import sys; sys.path.insert(0, "skills/snapgene-cro-cloning/scripts")
import snapgene as sg

vec = sg.read_dna("vector.dna")
orf = vec.feature("GeneName", ftype="CDS")
census = sg.check_sites_in_construct(vec.seq)       # both strands, circular
```

Python 3 standard library only.
