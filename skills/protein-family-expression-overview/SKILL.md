---
name: protein-family-expression-overview
description: Build an expression and purification overview of a protein or gene family as two matching deliverables - a self-contained tabbed HTML dashboard and a styled Word report - from one researched data file, with members grouped by what drives the purification strategy, a four-tier difficulty rating per member, and every row traceable to a database accession. Renders both outputs in one pass from one dataset so their counts cannot diverge, and reports inconsistencies rather than leaving them to the eye. Use when asked for an expression or purification overview of a protein family, a family or target reference, a comparison table across related proteins, or a dashboard plus a matching Word document.
---

# Protein family expression overview: dashboard + Word report

Two deliverables, one dataset, one run. Written for someone who already knows
how to purify proteins: lead with the non-obvious and leave out the basics.

## Requirements

Python 3 standard library only. `scripts/reference.py` renders both outputs and
`scripts/docx.py` writes the `.docx`; nothing is installed and no Node package
or converter is needed.

## Ask before assuming

- **The family and its members.** Never infer the member list. Ask for it, or
  for the accession set to derive it from.
- **The axes.** Which one or two categorical axes matter for *this* family
  (class, compartment, signalling mode, chain type). Do not carry over axes
  from a previous family; see the warning below.
- **Who it is for and what they will do with it.** The difficulty ratings and
  "key issues" differ completely between someone choosing a host and someone
  designing a construct.
- **Whether a difficulty rating may be stated at all.** It is judgement, not
  data. If they want only what is in the records, leave it out rather than
  inventing a ramp.

## Workflow

1. **Research first.** Gather facts before building anything. Cite only sources
   actually retrieved, and verify every reference — authors, journal, year —
   against the record rather than from recall. A reference with neither a DOI
   nor a URL is reported as unverifiable.
2. **Write the data file.** One JSON, described below. Every member carries its
   accession so each row is traceable.
3. **Build.** `python3 scripts/reference.py reference.json --out-dir .`
   Both files come out of one pass.
4. **Read the `PROBLEM:` lines** before looking at the output; they are listed
   below.
5. **Open both and look.** Check the filters work, the dashboard is legible in
   dark mode, and the Word tables are not clipped. This is the step no check
   replaces.
6. **Report in two or three lines** what each file contains. No recap of steps.

## The data file

```json
{
  "title": "...",
  "subtitle": "...",
  "db": "UNIPROT",
  "heuristic": "the one thing to remember, shown as an infobox",
  "axes": {
    "class":        {"label": "Class",      "values": ["...", "..."]},
    "localisation": {"label": "Compartment","values": ["...", "..."]}
  },
  "members": [
    {"name": "...", "accession": "P12345",
     "class": "...", "localisation": "...",
     "difficulty": "easy|moderate|hard|special",
     "system": "...", "issues": ["..."], "capture": "..."}
  ],
  "supergroups": [
    {"name": "...", "difficulty": "hard",
     "members": ["..."], "rationale": ["..."]}
  ],
  "structure_notes": [{"heading": "Conserved", "points": ["..."]}],
  "secondary": {"label": "Records", "columns": ["..."], "rows": [["..."]]},
  "matrix": {
    "columns": ["..."],
    "states": {"yes": {"label": "...", "colour": "#27500A",
                       "short": "Yes", "fill": "EAF3DE", "text": "27500A"}},
    "cells": {"MemberName": {"ColumnName": "yes"}},
    "note": "what absence means"
  },
  "references": [{"text": "...", "doi": "..."}]
}
```

**Declare `values` for every axis.** Without them the values are read off the
members, so a misspelt one validates against itself and appears as its own
category with its own colour. The generator reports an axis with no declared
values for exactly this reason.

**Never hard-code an axis into the skill.** The predecessor of this skill fixed
its second axis to one family's signalling classes and fixed a matrix cell to
one family's co-receptor. For any other family that renders a column which does
not apply and a matrix expecting a co-factor that does not exist. Axes and
matrix states are declared per family; a family with one axis declares one.

## Conventions

- **Difficulty is a fixed four-tier scale**, used exactly:
  Easy `●●○○` / Moderate `●●●○` / Hard `●●●●` / Special `★★★○`.
  Special means *a different workflow*, not *harder* — an intracellular target
  among secreted ones, say. A fifth circle could not say that.
- **Supergroups are defined by what drives the purification strategy** — the
  capture handle, stability, disulfides, glycosylation, oligomerisation — not
  by phylogeny. This is the part of the document that is actually worth
  reading; sequence similarity is already in the database.
- **"Key issues" are concrete and specific**: a construct boundary that is not
  residue 1, a required cofactor a chelating buffer will strip, an isoform or
  start-codon ambiguity, a monomer–dimer split on SEC. Not "use a His-tag".
- Dense bullets, not prose. One infobox, carrying the master heuristic.
- State plainly which ratings are judgement rather than data.

## What the generator reports

Each as a `PROBLEM:` line, all of them in one run rather than failing on the
first:

- a difficulty outside the four-tier scale, on a member or a group
- a member with no accession, so the row cannot be traced to a record
- an axis value that is not declared, and an axis that declares none
- a supergroup naming a member that is not in the member list
- a member that is in no supergroup
- a matrix row, column or cell state that is not declared
- a reference with neither DOI nor URL
- a duplicated or unnamed member
- anything `validate_docx` finds in the written Word file

## Why both outputs come from one run

The obvious design is two generators reading the same file, which then needs a
QA step that diffs one output against the other for drifted counts and
ratings. Here both renderers read one in-memory object built once, so there is
nothing to diff and no drift to check for. The member count in the dashboard
header and the count in the Word overview are the same expression.

## Pitfalls

- **A dashboard table will crush its widest column** rather than scroll, if the
  table has no minimum width: the browser satisfies `width: 100%` by wrapping
  the longest text column to about one character per line. The table carries a
  `min-width` and the wrapper scrolls, with the member column pinned.
- **Do not put a rating in the document that you would not defend out loud.**
  The difficulty column is the most quoted part of a reference like this and
  the least evidenced. Where it is a guess, say so in the subtitle.
- **Absence in the matrix is not absence in nature.** It means the record does
  not annotate it. Say which in the matrix `note`, every time.
- **Check the glyphs survive.** `● ○ ★ β α` render in Arial, but confirm after
  any conversion to PDF rather than assuming.
- **The visual pass is not optional.** Every check here is structural; none of
  them can see an unreadable colour pairing, a table clipped at the page edge
  or a filter that returns nothing.
