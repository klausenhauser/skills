#!/usr/bin/env python3
"""Build a family reference as an HTML dashboard and a matching Word report.

    python3 reference.py reference.json --out-dir .

Writes `<stem>_Dashboard.html` and `<stem>_Reference.docx` from one data file.

Both deliverables are rendered from the same in-memory object in a single run,
so their member counts, difficulty ratings and group membership cannot drift
apart. That is the point: the obvious design is two generators reading the same
file, which then needs a QA step that diffs one output against the other. Here
there is nothing to diff.

Standard library only; the Word half uses the sibling `docx.py`.

The data file
-------------
Everything family-specific lives in the JSON, including the axis names. A
reference for a secreted-ligand family and one for a kinase family differ only
in their data, never in this module:

    {
      "title": "...", "subtitle": "...",
      "axes": {
        "class": {"label": "Subfamily"},
        "mode":  {"label": "Signalling mode"}
      },
      "members": [
        {"name": "...", "accession": "P12345", "class": "...", "mode": "...",
         "difficulty": "moderate", "system": "E. coli",
         "issues": ["..."], "capture": "..."}
      ],
      "supergroups": [
        {"name": "...", "difficulty": "hard",
         "members": ["..."], "rationale": ["..."]}
      ],
      "structure_notes": [{"heading": "...", "points": ["..."]}],
      "secondary": {"label": "...", "columns": [...], "rows": [[...]]},
      "matrix": {"columns": [...], "states": {...}, "cells": {...}},
      "references": [{"text": "...", "doi": "..."}]
    }

`axes` is declarative on purpose. Hard-coding the second axis -- the original
of this skill fixed it to one family's signalling classes, and fixed a matrix
cell to one family's co-receptor -- makes the skill silently wrong for any
other family: the dashboard renders a column that does not apply and a matrix
expecting a co-factor that does not exist. Declared axes each get a palette, a
filter and a legend, and a family with one axis simply declares one.

What it checks
--------------
`check()` reports what the original left to a human eye: a difficulty outside
the scale, a supergroup naming a member that is not in the member list, a
matrix row or column that is not a member, a member with no accession, a
reference with neither DOI nor URL, and any member missing from every
supergroup. It returns problems rather than raising, so one run lists them all.
"""

import argparse
import html
import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import docx as D  # noqa: E402

# ----------------------------------------------------------- vocabulary

# Four tiers, fixed. The glyphs are deliberate: a filled-circle ramp reads as
# ordinal at a glance, and the star tier marks "different workflow" rather than
# "harder", which a fifth circle could not say.
DIFFICULTY = {
    "easy":     {"label": "Easy",     "glyph": "●●○○",
                 "bg": "EAF3DE", "fg": "27500A", "rank": 1},
    "moderate": {"label": "Moderate", "glyph": "●●●○",
                 "bg": "FAEEDA", "fg": "633806", "rank": 2},
    "hard":     {"label": "Hard",     "glyph": "●●●●",
                 "bg": "FCEBEB", "fg": "791F1F", "rank": 3},
    "special":  {"label": "Special",  "glyph": "★★★○",
                 "bg": "FBEAF0", "fg": "72243E", "rank": 4},
}

# Cycled for any declared categorical axis, so no axis is privileged.
AXIS_PALETTE = [
    ("EEEDFE", "3C3489"), ("E6F1FB", "0C447C"), ("EAF3DE", "27500A"),
    ("FAEEDA", "633806"), ("FAECE7", "993C1D"), ("FBEAF0", "72243E"),
    ("F1EFE8", "444441"),
]

NAVY, BLUE, MUTED, BORDER, ALTROW = "1B3A6B", "2563A8", "555550", "CCCCCC", "F5F8FF"
INFO_BG, INFO_FG = "E6F1FB", "0C447C"

DB_URL = {
    "UNIPROT": "https://www.uniprot.org/uniprotkb/{acc}",
    "PDB": "https://www.rcsb.org/structure/{acc}",
    "GENBANK": "https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
}


def _esc(value):
    return html.escape(str(value if value is not None else ""))


def _plural(n, word):
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def _slug(value):
    return re.sub(r"[^a-z0-9]+", "-", str(value).lower()).strip("-") or "x"


# -------------------------------------------------------------- the data

class Reference:
    """One family reference. Both renderers read this and nothing else."""

    def __init__(self, raw):
        self.raw = raw
        self.title = raw.get("title") or "Family reference"
        self.subtitle = raw.get("subtitle") or ""
        self.db = (raw.get("db") or "UNIPROT").upper()
        self.members = list(raw.get("members") or [])
        self.supergroups = list(raw.get("supergroups") or [])
        self.structure_notes = list(raw.get("structure_notes") or [])
        self.secondary = raw.get("secondary") or None
        self.matrix = raw.get("matrix") or None
        self.references = list(raw.get("references") or [])
        self.heuristic = raw.get("heuristic") or ""
        self.axes = self._axes(raw.get("axes") or {})

    def _axes(self, declared):
        """Resolve declared axes to {key: {label, values, colours}}."""
        out = {}
        for key, spec in declared.items():
            spec = spec or {}
            values = spec.get("values") or self._observed(key)
            colours = {}
            for i, v in enumerate(values):
                colours[v] = AXIS_PALETTE[i % len(AXIS_PALETTE)]
            out[key] = {"label": spec.get("label") or key.title(),
                        "values": values, "colours": colours,
                        # Whether the value list came from the data file or was
                        # inferred. An inferred list cannot catch a typo,
                        # because the typo is in it; check() says so.
                        "declared": bool(spec.get("values"))}
        return out

    def _observed(self, key):
        seen = []
        for m in self.members:
            v = (m.get(key) or "").strip()
            if v and v not in seen:
                seen.append(v)
        return seen

    def difficulty(self, member):
        return DIFFICULTY.get(str(member.get("difficulty", "")).lower())

    def url(self, accession):
        if not accession:
            return None
        return DB_URL.get(self.db, DB_URL["UNIPROT"]).format(acc=accession)

    def counts(self):
        """Per-tier and per-axis counts, computed once and used by both."""
        tiers = {k: 0 for k in DIFFICULTY}
        for m in self.members:
            d = str(m.get("difficulty", "")).lower()
            if d in tiers:
                tiers[d] += 1
        axes = {}
        for key, axis in self.axes.items():
            axes[key] = {v: sum(1 for m in self.members
                                if (m.get(key) or "").strip() == v)
                         for v in axis["values"]}
        return {"members": len(self.members), "difficulty": tiers,
                "axes": axes, "supergroups": len(self.supergroups)}


def load(path):
    return Reference(json.loads(pathlib.Path(path).read_text(encoding="utf-8")))


# ------------------------------------------------------------- checking

def check(ref):
    """Report inconsistencies the eye would have to catch. Never raises."""
    problems = []
    names = [m.get("name", "") for m in ref.members]
    seen = set()
    for n in names:
        if not n:
            problems.append("a member has no name")
        elif n in seen:
            problems.append(f"member {n!r} appears more than once")
        seen.add(n)

    for m in ref.members:
        n = m.get("name", "?")
        if not ref.difficulty(m):
            problems.append(
                f"{n}: difficulty {m.get('difficulty')!r} is not one of "
                f"{', '.join(sorted(DIFFICULTY))}")
        if not m.get("accession"):
            problems.append(f"{n}: no accession, so the entry cannot be "
                            f"traced to a database record")
        for key in ref.axes:
            v = (m.get(key) or "").strip()
            if v and v not in ref.axes[key]["values"]:
                problems.append(f"{n}: {key}={v!r} is not a declared value "
                                f"of that axis")

    for key, axis in ref.axes.items():
        if not axis["declared"]:
            problems.append(
                f"axis {key!r} declares no `values`, so they were read off the "
                f"members and a misspelt one validates against itself; list "
                f"them in the data file to have it checked")

    grouped = set()
    for g in ref.supergroups:
        gn = g.get("name", "?")
        if not ref.difficulty(g) and g.get("difficulty"):
            problems.append(f"supergroup {gn}: difficulty "
                            f"{g.get('difficulty')!r} is not in the scale")
        for member in g.get("members") or []:
            if member not in seen:
                problems.append(f"supergroup {gn} names {member!r}, which is "
                                f"not in the member list")
            grouped.add(member)
    if ref.supergroups:
        for n in names:
            if n and n not in grouped:
                problems.append(f"{n} is in no supergroup")

    if ref.matrix:
        cols = ref.matrix.get("columns") or []
        states = ref.matrix.get("states") or {}
        for row, cells in (ref.matrix.get("cells") or {}).items():
            if row not in seen:
                problems.append(f"matrix row {row!r} is not a member")
            for col, state in (cells or {}).items():
                if col not in cols:
                    problems.append(f"matrix column {col!r} of {row} is not "
                                    f"in the declared columns")
                if state not in states:
                    problems.append(f"matrix cell {row}/{col} is {state!r}, "
                                    f"which is not a declared state")
        if not states:
            problems.append("matrix has no `states`, so its cells cannot be "
                            "drawn or explained in a legend")

    for i, r in enumerate(ref.references, 1):
        if not (r.get("doi") or r.get("url")):
            problems.append(f"reference {i} has neither a DOI nor a URL, so "
                            f"it cannot be verified")
    return problems


# ------------------------------------------------------------- dashboard

_CSS = """
:root {
  --page: #eeecea; --card: #ffffff; --ink: #1a1a18; --muted: #6a6a64;
  --line: rgba(0,0,0,.12); --accent: #2563A8;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --page: #1c1c1a; --card: #262624; --ink: #ecebe6; --muted: #a3a29b;
    --line: rgba(255,255,255,.14); --accent: #7fb0e6;
  }
}
:root[data-theme="dark"] {
  --page: #1c1c1a; --card: #262624; --ink: #ecebe6; --muted: #a3a29b;
  --line: rgba(255,255,255,.14); --accent: #7fb0e6;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
  Helvetica, Arial, sans-serif; }
.wrap { max-width: 1400px; margin: 0 auto; padding: 24px 16px 64px; }
header { display: flex; justify-content: space-between; align-items: flex-end;
  gap: 24px; flex-wrap: wrap; margin-bottom: 18px; }
h1 { font-size: 24px; margin: 0 0 4px; letter-spacing: -.2px; }
.sub { color: var(--muted); max-width: 70ch; }
.meta { color: var(--muted); font-size: 12.5px; text-align: right; }
.tabs { display: flex; gap: 4px; flex-wrap: wrap; border-bottom: 1px solid var(--line);
  margin-bottom: 18px; }
.tabs button { appearance: none; border: 0; background: transparent;
  color: var(--muted); font: inherit; font-size: 13px; padding: 9px 13px;
  border-radius: 8px 8px 0 0; cursor: pointer; }
.tabs button:hover { color: var(--ink); }
.tabs button.active { color: var(--ink); background: var(--card);
  box-shadow: inset 0 -2px 0 var(--accent); }
.panel { display: none; }
.panel.active { display: block; }
.card { background: var(--card); border: .5px solid var(--line);
  border-radius: 12px; padding: 16px; margin-bottom: 14px; }
.filters { display: flex; gap: 10px; flex-wrap: wrap; align-items: center;
  margin-bottom: 12px; }
.filters label { font-size: 12.5px; color: var(--muted); }
.filters select, .filters input { font: inherit; font-size: 13px;
  padding: 5px 8px; border-radius: 8px; border: .5px solid var(--line);
  background: var(--card); color: var(--ink); }
.n-of { margin-left: auto; font-size: 12.5px; color: var(--muted);
  font-variant-numeric: tabular-nums; }
.scroll { max-height: 70vh; overflow: auto; border-radius: 8px; }
/* Without a minimum width the browser satisfies width:100% by crushing the
   widest text column to roughly one character per line rather than letting the
   wrapper scroll. Narrow viewports are the common case for a reference someone
   opens beside a protocol, so the table scrolls horizontally instead. */
table { border-collapse: collapse; width: 100%; min-width: 900px;
  font-size: 12.5px; }
td:first-child, th:first-child { position: sticky; left: 0;
  background: var(--card); }
thead th { position: sticky; top: 0; background: var(--card); text-align: left;
  font-weight: 600; padding: 9px 10px; border-bottom: 1px solid var(--line);
  white-space: nowrap; z-index: 1; }
tbody td { padding: 9px 10px; border-bottom: .5px solid var(--line);
  vertical-align: top; }
tbody tr:hover td { background: rgba(127,127,127,.06); }
.pill { display: inline-block; border-radius: 100px; padding: 2px 9px;
  font-size: 10.5px; font-weight: 500; white-space: nowrap; }
.issues { margin: 0; padding-left: 15px; }
.issues li { margin: 1px 0; }
.groups { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  gap: 14px; }
.group { background: var(--card); border: .5px solid var(--line);
  border-left: 3px solid var(--g, var(--accent)); border-radius: 8px;
  padding: 14px 16px; }
.group h3 { margin: 0 0 4px; font-size: 15px; }
.chips { display: flex; gap: 5px; flex-wrap: wrap; margin: 8px 0 10px; }
.chip { font-size: 11px; padding: 2px 8px; border-radius: 100px;
  background: rgba(127,127,127,.12); }
.lab { font-size: 10.5px; letter-spacing: .6px; text-transform: uppercase;
  color: var(--muted); margin-bottom: 5px; }
.dash { margin: 0; padding: 0; list-style: none; }
.dash li { position: relative; padding-left: 14px; margin: 3px 0;
  font-size: 12.5px; }
.dash li:before { content: "–"; position: absolute; left: 0;
  color: var(--muted); }
.info { background: var(--i-bg, #E6F1FB); color: var(--i-fg, #0C447C);
  border: 1px solid var(--accent); border-radius: 10px; padding: 12px 14px;
  margin-bottom: 14px; }
.dot { display: inline-block; width: 11px; height: 11px; border-radius: 50%; }
.legend { display: flex; gap: 18px; flex-wrap: wrap; margin-top: 10px;
  font-size: 12.5px; color: var(--muted); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
footer { margin-top: 26px; color: var(--muted); font-size: 12px; }
footer a { color: var(--accent); }
a { color: var(--accent); }
.vis-hidden { display: none !important; }
@media (max-width: 760px) {
  .wrap { padding: 16px; }
  header { flex-direction: column; align-items: flex-start; }
  .meta { text-align: left; }
}
"""

_JS = """
function switchTab(id, btn) {
  document.querySelectorAll('.panel').forEach(function (p) {
    p.classList.toggle('active', p.id === id);
  });
  document.querySelectorAll('.tabs button').forEach(function (b) {
    b.classList.toggle('active', b === btn);
  });
}
function applyFilters() {
  var sels = Array.prototype.slice.call(
    document.querySelectorAll('[data-filter]'));
  var q = (document.getElementById('q') || {}).value || '';
  q = q.trim().toLowerCase();
  var rows = document.querySelectorAll('#master tbody tr');
  var shown = 0;
  rows.forEach(function (tr) {
    var ok = sels.every(function (s) {
      return !s.value || tr.getAttribute('data-' + s.dataset.filter) === s.value;
    });
    if (ok && q) { ok = tr.textContent.toLowerCase().indexOf(q) !== -1; }
    tr.classList.toggle('vis-hidden', !ok);
    if (ok) { shown++; }
  });
  var c = document.getElementById('nof');
  if (c) { c.textContent = shown + ' of ' + rows.length; }
}
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('[data-filter], #q').forEach(function (el) {
    el.addEventListener(el.tagName === 'INPUT' ? 'input' : 'change', applyFilters);
  });
  var clear = document.getElementById('clear');
  if (clear) {
    clear.addEventListener('click', function () {
      document.querySelectorAll('[data-filter]').forEach(function (s) {
        s.value = '';
      });
      var q = document.getElementById('q');
      if (q) { q.value = ''; }
      applyFilters();
    });
  }
  applyFilters();
});
"""


def _pill(label, bg, fg):
    return (f'<span class="pill" style="background:#{bg};color:#{fg}">'
            f'{_esc(label)}</span>')


def _difficulty_pill(spec):
    if not spec:
        return '<span class="pill" style="background:#eee;color:#555">?</span>'
    return _pill(f'{spec["label"]} {spec["glyph"]}', spec["bg"], spec["fg"])


def render_dashboard(ref):
    c = ref.counts()
    axis_keys = list(ref.axes)

    # --- master table -------------------------------------------------
    head = ["Member"] + [ref.axes[k]["label"] for k in axis_keys] + [
        "Difficulty", "Preferred system", "Key issues", "Capture step"]
    body = []
    for m in ref.members:
        spec = ref.difficulty(m)
        attrs = "".join(
            f' data-{k}="{_esc((m.get(k) or "").strip())}"' for k in axis_keys)
        attrs += f' data-difficulty="{_esc(str(m.get("difficulty","")).lower())}"'
        acc = m.get("accession")
        name = (f'<a href="{_esc(ref.url(acc))}">{_esc(m.get("name"))}</a>'
                if acc else _esc(m.get("name")))
        cells = [f"<strong>{name}</strong>"
                 + (f'<br><span style="color:var(--muted);font-size:11px">'
                    f'{_esc(acc)}</span>' if acc else "")]
        for k in axis_keys:
            v = (m.get(k) or "").strip()
            if v:
                bg, fg = ref.axes[k]["colours"].get(v, ("F1EFE8", "444441"))
                cells.append(_pill(v, bg, fg))
            else:
                cells.append("")
        cells.append(_difficulty_pill(spec))
        cells.append(_esc(m.get("system")))
        issues = m.get("issues") or []
        cells.append("<ul class='issues'>"
                     + "".join(f"<li>{_esc(i)}</li>" for i in issues)
                     + "</ul>" if issues else "")
        cells.append(_esc(m.get("capture")))
        body.append(f"<tr{attrs}>"
                    + "".join(f"<td>{v}</td>" for v in cells) + "</tr>")

    filters = []
    for k in axis_keys:
        opts = "".join(
            f'<option value="{_esc(v)}">{_esc(v)} '
            f'({c["axes"][k].get(v, 0)})</option>'
            for v in ref.axes[k]["values"])
        filters.append(
            f'<label>{_esc(ref.axes[k]["label"])} '
            f'<select data-filter="{_esc(k)}">'
            f'<option value="">All</option>{opts}</select></label>')
    dopts = "".join(
        f'<option value="{k}">{v["label"]} ({c["difficulty"][k]})</option>'
        for k, v in sorted(DIFFICULTY.items(), key=lambda kv: kv[1]["rank"])
        if c["difficulty"][k])
    filters.append('<label>Difficulty <select data-filter="difficulty">'
                   f'<option value="">All</option>{dopts}</select></label>')
    filters.append('<label>Search <input id="q" type="search" '
                   'placeholder="any text"></label>')
    filters.append('<button id="clear" class="chip" '
                   'style="cursor:pointer;border:0">Clear filters</button>')

    master = (
        '<div class="card"><div class="filters">' + "".join(filters)
        + '<span class="n-of" id="nof"></span></div>'
        '<div class="scroll"><table id="master"><thead><tr>'
        + "".join(f"<th>{_esc(h)}</th>" for h in head)
        + "</tr></thead><tbody>" + "".join(body) + "</tbody></table></div></div>")

    # --- supergroups ---------------------------------------------------
    cards = []
    for g in ref.supergroups:
        spec = ref.difficulty(g)
        colour = f'#{spec["fg"]}' if spec else "var(--accent)"
        chips = "".join(f'<span class="chip">{_esc(x)}</span>'
                        for x in g.get("members") or [])
        points = "".join(f"<li>{_esc(p)}</li>" for p in g.get("rationale") or [])
        cards.append(
            f'<div class="group" style="--g:{colour}">'
            f'<h3>{_esc(g.get("name"))}</h3>'
            + (_difficulty_pill(spec) if spec else "")
            + f'<div class="chips">{chips}</div>'
            f'<div class="lab">Why they group</div>'
            f'<ul class="dash">{points}</ul></div>')
    groups_panel = (
        (f'<div class="info">{_esc(ref.heuristic)}</div>' if ref.heuristic else "")
        + f'<div class="groups">{"".join(cards)}</div>')

    # --- structure notes -----------------------------------------------
    notes = []
    for n in ref.structure_notes:
        pts = "".join(f"<li>{_esc(p)}</li>" for p in n.get("points") or [])
        notes.append(f'<div class="card"><h3 style="margin-top:0">'
                     f'{_esc(n.get("heading"))}</h3>'
                     f'<ul class="dash">{pts}</ul></div>')

    # --- secondary table ------------------------------------------------
    secondary = ""
    if ref.secondary:
        cols = ref.secondary.get("columns") or []
        rows = ref.secondary.get("rows") or []
        secondary = (
            '<div class="card"><div class="scroll"><table><thead><tr>'
            + "".join(f"<th>{_esc(x)}</th>" for x in cols)
            + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join(f"<td>{_esc(v)}</td>" for v in r)
                      + "</tr>" for r in rows)
            + "</tbody></table></div></div>")

    # --- matrix ----------------------------------------------------------
    matrix = ""
    if ref.matrix:
        cols = ref.matrix.get("columns") or []
        states = ref.matrix.get("states") or {}
        trs = []
        for m in ref.members:
            name = m.get("name")
            cells = (ref.matrix.get("cells") or {}).get(name) or {}
            tds = []
            for col in cols:
                st = states.get(cells.get(col))
                tds.append(
                    f'<td style="text-align:center">'
                    f'<span class="dot" style="background:'
                    f'{_esc((st or {}).get("colour", "#cccccc"))}" '
                    f'title="{_esc((st or {}).get("label", "not stated"))}">'
                    f"</span></td>")
            trs.append(f"<tr><td><strong>{_esc(name)}</strong></td>"
                       + "".join(tds) + "</tr>")
        legend = "".join(
            f'<span><span class="dot" style="background:'
            f'{_esc(v.get("colour", "#ccc"))}"></span>{_esc(v.get("label", k))}'
            f"</span>" for k, v in states.items())
        note = ref.matrix.get("note") or ""
        matrix = ('<div class="card"><div class="scroll"><table><thead><tr>'
                  '<th>Member</th>'
                  + "".join(f'<th style="text-align:center">{_esc(x)}</th>'
                            for x in cols)
                  + "</tr></thead><tbody>" + "".join(trs)
                  + f'</tbody></table></div><div class="legend">{legend}</div>'
                  + (f'<p style="color:var(--muted);font-size:12px">'
                     f'{_esc(note)}</p>' if note else "")
                  + "</div>")

    tabs = [("master", "Master table", master),
            ("groups", "Purification groups", groups_panel)]
    if secondary:
        tabs.append(("secondary",
                     ref.secondary.get("label") or "Secondary", secondary))
    if notes:
        tabs.append(("structure", "Structural features", "".join(notes)))
    if matrix:
        tabs.append(("matrix", "Interaction matrix", matrix))

    bar = "".join(
        f'<button class="{"active" if i == 0 else ""}" '
        f"onclick=\"switchTab('{t[0]}', this)\">{_esc(t[1])}</button>"
        for i, t in enumerate(tabs))
    panels = "".join(
        f'<section class="panel {"active" if i == 0 else ""}" '
        f'id="{t[0]}">{t[2]}</section>' for i, t in enumerate(tabs))

    refs = "".join(
        f'<li>{_esc(r.get("text"))}'
        + (f' <a href="https://doi.org/{_esc(r["doi"])}">'
           f'doi:{_esc(r["doi"])}</a>' if r.get("doi")
           else (f' <a href="{_esc(r["url"])}">link</a>' if r.get("url") else ""))
        + "</li>" for r in ref.references)

    return (
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{_esc(ref.title)}</title><style>{_CSS}</style></head><body>"
        f'<div class="wrap"><header><div><h1>{_esc(ref.title)}</h1>'
        f'<div class="sub">{_esc(ref.subtitle)}</div></div>'
        f'<div class="meta">{c["members"]} members · '
        f'{_plural(c["supergroups"], "group")} · '
        f'{_plural(len(ref.references), "reference")}</div>'
        f'</header><nav class="tabs">{bar}</nav>{panels}'
        + (f"<footer><strong>References</strong><ol>{refs}</ol></footer>"
           if refs else "")
        + f"</div><script>{_JS}</script></body></html>")


# ----------------------------------------------------------------- report

def render_report(ref, path):
    c = ref.counts()
    doc = D.Document(page="LETTER", margins=(1080, 1080, 1080, 1080),
                     font="Arial", size_pt=10, title=ref.title)
    doc.heading_style(1, size_pt=16, colour=NAVY)
    doc.heading_style(2, size_pt=13, colour=BLUE)
    doc.heading_style(3, size_pt=11, colour=BLUE)
    doc.header(ref.title, right="Family reference")
    doc.footer_page_x_of_y(left=f"{c['members']} members")

    doc.para(ref.title, bold=True, size_pt=22, colour=NAVY, space_after=D.pt(2))
    if ref.subtitle:
        doc.para(ref.subtitle, size_pt=11, colour=MUTED, space_after=D.pt(10))
    doc.rule(size=8, colour=NAVY)

    if ref.heuristic:
        t = doc.table([doc.content_width], borders={
            s: (6, BLUE) for s in ("top", "left", "bottom", "right")},
            margins=(120, 140, 120, 140))
        t.row([D.Cell(ref.heuristic, fill=INFO_BG, colour=INFO_FG,
                      size_pt=10)])
        doc.spacer(8)

    doc.heading("1 Overview", level=1)
    axis_summary = "; ".join(
        f"{ref.axes[k]['label']}: "
        + ", ".join(f"{v} ({n})" for v, n in c["axes"][k].items() if n)
        for k in ref.axes)
    doc.para(f"{c['members']} members in {c['supergroups']} purification "
             f"groups. " + (axis_summary or ""))
    doc.bullets([f"{DIFFICULTY[k]['label']} {DIFFICULTY[k]['glyph']}: {n}"
                 for k, n in sorted(c["difficulty"].items(),
                                    key=lambda kv: DIFFICULTY[kv[0]]["rank"])
                 if n])

    if ref.supergroups:
        doc.heading("2 Purification groups", level=1)
        w = doc.columns(1.1, 0.5, 2.4)
        t = doc.table(w, header=["Group", "Difficulty", "Why they group"],
                      header_fill=NAVY, header_colour="FFFFFF")
        for i, g in enumerate(ref.supergroups):
            spec = ref.difficulty(g)
            blocks = [D.Paragraph(
                D._runs("– " + p, size_pt=9, font="Arial"),
                space_after=D.pt(1), space_before=0)
                for p in g.get("rationale") or []]
            members = ", ".join(g.get("members") or [])
            t.row([
                D.Cell(blocks=[D.Paragraph(
                    D._runs(g.get("name"), bold=True, size_pt=9, font="Arial"),
                    space_after=D.pt(1), space_before=0),
                    D.Paragraph(D._runs(members, size_pt=8, colour=MUTED,
                                        font="Arial"),
                                space_after=0, space_before=0)]),
                D.Cell(f'{spec["label"]}\n{spec["glyph"]}' if spec else "—",
                       fill=spec["bg"] if spec else None,
                       colour=spec["fg"] if spec else None, size_pt=8,
                       align="center"),
                D.Cell(blocks=blocks or [D.Paragraph([])]),
            ], fill=ALTROW if i % 2 else None, size_pt=9)

    doc.heading("3 Master table", level=1)
    axis_keys = list(ref.axes)
    fracs = [1.1] + [0.8] * len(axis_keys) + [0.55, 0.9, 2.0, 1.0]
    w = doc.columns(*fracs)
    head = ["Member"] + [ref.axes[k]["label"] for k in axis_keys] + [
        "Diff.", "System", "Key issues", "Capture"]
    t = doc.table(w, header=head, header_fill=NAVY, header_colour="FFFFFF")
    for i, m in enumerate(ref.members):
        spec = ref.difficulty(m)
        issues = [D.Paragraph(D._runs("– " + x, size_pt=8, font="Arial"),
                              space_after=D.pt(1), space_before=0)
                  for x in m.get("issues") or []]
        cells = [D.Cell(blocks=[
            D.Paragraph(D._runs(m.get("name"), bold=True, size_pt=9,
                                font="Arial"), space_after=0, space_before=0),
            D.Paragraph(D._runs(m.get("accession") or "", size_pt=7.5,
                                colour=MUTED, font="Arial"),
                        space_after=0, space_before=0)])]
        for k in axis_keys:
            cells.append(D.Cell((m.get(k) or ""), size_pt=8))
        cells.append(D.Cell(spec["glyph"] if spec else "?",
                            fill=spec["bg"] if spec else None,
                            colour=spec["fg"] if spec else None,
                            size_pt=8, align="center"))
        cells.append(D.Cell(m.get("system") or "", size_pt=8))
        cells.append(D.Cell(blocks=issues or [D.Paragraph([])]))
        cells.append(D.Cell(m.get("capture") or "", size_pt=8))
        t.row(cells, fill=ALTROW if i % 2 else None)

    n = 4
    if ref.structure_notes:
        doc.heading(f"{n} Structural features", level=1)
        n += 1
        for note in ref.structure_notes:
            doc.heading(note.get("heading") or "", level=2)
            doc.bullets(note.get("points") or [], size_pt=9)

    if ref.secondary:
        doc.heading(f"{n} {ref.secondary.get('label') or 'Secondary'}", level=1)
        n += 1
        cols = ref.secondary.get("columns") or []
        if cols:
            w = doc.columns(*([1] * len(cols)))
            t = doc.table(w, header=cols, header_fill=NAVY,
                          header_colour="FFFFFF")
            for i, row in enumerate(ref.secondary.get("rows") or []):
                t.row([D.Cell(v, size_pt=8) for v in row],
                      fill=ALTROW if i % 2 else None)

    if ref.matrix:
        doc.heading(f"{n} Interaction matrix", level=1)
        n += 1
        cols = ref.matrix.get("columns") or []
        states = ref.matrix.get("states") or {}
        w = doc.columns(*([1.4] + [1] * len(cols)))
        t = doc.table(w, header=["Member"] + cols, header_fill=NAVY,
                      header_colour="FFFFFF")
        for i, m in enumerate(ref.members):
            cells = [D.Cell(m.get("name"), bold=True, size_pt=8)]
            row = (ref.matrix.get("cells") or {}).get(m.get("name")) or {}
            for col in cols:
                st = states.get(row.get(col)) or {}
                cells.append(D.Cell(st.get("short", "–"),
                                    fill=st.get("fill"),
                                    colour=st.get("text"),
                                    size_pt=8, align="center"))
            t.row(cells, fill=ALTROW if i % 2 else None)
        doc.spacer(4)
        doc.para("  ".join(f'{v.get("short", k)} = {v.get("label", k)}'
                           for k, v in states.items()),
                 size_pt=8, colour=MUTED)
        if ref.matrix.get("note"):
            doc.para(ref.matrix["note"], size_pt=8, colour=MUTED)

    if ref.references:
        doc.heading(f"{n} References", level=1)
        doc.numbered([
            r.get("text", "") + (f'  doi:{r["doi"]}' if r.get("doi")
                                 else (f'  {r["url"]}' if r.get("url") else ""))
            for r in ref.references], size_pt=9)

    doc.save(path)
    return path


# -------------------------------------------------------------------- main

def build(data_path, out_dir=".", stem=None):
    ref = load(data_path)
    problems = check(ref)
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = stem or _slug(ref.title)

    html_path = out / f"{stem}_Dashboard.html"
    html_path.write_text(render_dashboard(ref), encoding="utf-8")
    docx_path = out / f"{stem}_Reference.docx"
    render_report(ref, docx_path)

    report = D.validate_docx(docx_path)
    problems.extend(f"docx: {p}" for p in report["problems"])
    return {"problems": problems, "written": [html_path, docx_path],
            "counts": ref.counts(), "docx_checks": report["checks"]}


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--stem")
    a = ap.parse_args(argv[1:])
    r = build(a.data, a.out_dir, a.stem)
    c = r["counts"]
    print(f"  {c['members']} members, {c['supergroups']} groups")
    print("  difficulty: " + ", ".join(
        f"{DIFFICULTY[k]['label']} {v}" for k, v in c["difficulty"].items() if v))
    for k, v in r["docx_checks"].items():
        print(f"  docx {k}: {v}")
    for p in r["written"]:
        print(f"  wrote {p}")
    for p in r["problems"]:
        print(f"  PROBLEM: {p}")
    return 1 if r["problems"] else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
