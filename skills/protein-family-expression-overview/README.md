# protein-family-expression-overview

Build an expression and purification overview of a protein or gene family as a
tabbed HTML dashboard and a matching styled Word report, from one data file.

## What it does

1. **Research first**, citing only sources actually retrieved.
2. **Write one JSON** describing the family: members with accessions, one or two
   declared categorical axes, purification groups, structural notes, an optional
   secondary table and an optional matrix.
3. **Run one command.** Both deliverables come out of a single pass:
   `<stem>_Dashboard.html` and `<stem>_Reference.docx`.
4. **Read the reported problems**, then open both and look at them.

```
python3 scripts/reference.py reference.json --out-dir .
```

## Why it is built this way

### Both outputs come from one dataset, in one run

The obvious design is two generators reading the same file, which then needs a
QA step diffing one output against the other for counts and ratings that have
drifted apart. Here both renderers read one in-memory object built once, so
there is nothing to diff: the member count in the dashboard header and the one
in the Word overview are the same expression, not two expressions that agree.
A whole class of checking disappears instead of being automated.

### The axes are declared, not built in

This replaces an earlier skill that hard-coded its second axis to one protein
family's signalling classes, and hard-coded a matrix cell to that family's
co-receptor. Those were invisible while it was used for that family and wrong
for every other one: another family gets a column that does not apply and a
matrix expecting a co-factor that does not exist.

So axes, their values and the matrix states are all declared per family. The
module contains no family. A family with one axis declares one; a family with
three declares three, each getting its own palette, filter and legend.

Declaring the axis *values* also makes them checkable. If they are left out
they are read off the members, which means a misspelt value validates against
itself and shows up as its own category with its own colour — so the generator
reports an axis that declares none, rather than appearing to check something it
cannot.

### It reports what the checklist used to ask a human to notice

The predecessor ended with a QA checklist: same counts in both outputs, filters
work, every claim traceable, uncertain ratings flagged. Most of that is now
either structural or reported — a difficulty outside the scale, a supergroup
naming someone who is not a member, a member in no supergroup, a matrix cell in
an undeclared state, a reference with neither DOI nor URL, a member with no
accession. All problems come back from one run rather than the first one
raising.

What is left on the checklist is the part a program genuinely cannot do: open
both files and look at them. No check here can see an unreadable colour
pairing, a table clipped at the page edge, or a filter combination that returns
nothing.

### Difficulty is four tiers and says so

Easy `●●○○`, Moderate `●●●○`, Hard `●●●●`, Special `★★★○`. The ramp reads as
ordinal at a glance; the star tier means *a different workflow* rather than
*harder*, which a fifth circle could not express. It is also the most quoted
and least evidenced column in a document like this, which is why the skill
insists on stating where a rating is judgement rather than data.

### What testing it changed

- `docx.py` duck-typed any non-string cell as a sequence, so an integer — a
  length, a count, a year — failed deep inside the writer with `'int' object is
  not iterable`. Scalars are now coerced with `str()`.
- The dashboard table crushed its widest column to about one character per line
  instead of scrolling, because `width: 100%` with no `min-width` is satisfied
  by wrapping rather than overflowing. It now has a minimum width, scrolls
  inside its card, and pins the member column.
- Checking against a deliberately corrupted data file was what showed that the
  axis check was vacuous without declared values.

The Word output was verified by an independent reader: pandoc parses the
generated `.docx` and recovers the headings, the infobox, the tables and the
numbered references.

## Layout

```
SKILL.md               the workflow, data-file spec, conventions and pitfalls
scripts/reference.py   reference.json -> dashboard HTML + Word report
scripts/docx.py        .docx reader, writer and structural validator
```

`docx.py` is also carried by another skill. Skills in this marketplace are
self-contained — a skill may not reference a file outside its own folder, and
there is no shared folder — so the duplication is deliberate. A fix to one copy
has to be applied to the other, which is the price of each skill working on its
own.

## Requirements

Python 3 standard library only.
