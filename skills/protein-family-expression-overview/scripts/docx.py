#!/usr/bin/env python3
"""Read and write Word .docx files with the Python standard library.

A .docx is a zip of OOXML parts. Writing one needs no library: the parts are
XML, and this module emits exactly the ones Word and LibreOffice require. It
replaces two external dependencies at once -- docx-js (Node) for writing and
pandoc for reading -- so a skill that uses it runs on a bare Python install.

Writing
-------
    doc = Document(page="A4", margin_cm=2.0, font="Arial", size_pt=10)
    doc.header("COURSE · Worksheet · v1-2026")
    doc.footer_page_x_of_y()
    doc.heading("Purification of a fluorescent protein", level=1)
    doc.para("Read the whole protocol before you start.")
    doc.callout("Goal", "Purify and quantify one His-tagged protein.")
    t = doc.table([2400, 1800, 5438], header=["Step", "Volume", "Tube"])
    t.row(["Load", "5 mL", "A1"])
    doc.bullets(["Check the calculation.", "Then pipette."])
    doc.save("handout.docx")

Reading
-------
    for p in read_docx("input.docx"):
        print(p.style, p.text)
    print(docx_text("input.docx"))

Checking
--------
    report = validate_docx("handout.docx")   # {"ok": bool, "problems": [...]}

Units. Word measures in DXA (twentieths of a point, 1440 per inch) for page
and table geometry, half-points for font size, and eighths of a point for
border widths. `cm()` and `pt()` convert. A4 is 11906 x 16838 DXA; at 2 cm
margins the content width is 9638, and every table's column widths must sum to
exactly that or Word silently rescales the table.

Why the spacing helpers always set lineRule. Passing `w:line` without
`w:lineRule="auto"` makes LibreOffice treat the value as an exact height
rather than a multiple, and lines overlap. Every spacing call here sets both.
"""

import re
import zipfile
from xml.sax.saxutils import escape as _xml_escape

# --------------------------------------------------------------- units

DXA_PER_CM = 566.929133858
DXA_PER_INCH = 1440

PAGES = {                      # (width, height) in DXA
    "A4": (11906, 16838),
    "LETTER": (12240, 15840),
    "A4_LANDSCAPE": (16838, 11906),
    "LETTER_LANDSCAPE": (15840, 12240),
}


def cm(value):
    """Centimetres to DXA."""
    return int(round(value * DXA_PER_CM))


def pt(value):
    """Points to DXA."""
    return int(round(value * 20))


def _attr(value):
    """Quote an XML attribute value, escaping exactly once.

    quoteattr() would also do this, but it switches to single quotes when the
    value contains a double quote, and whitespace characters must become
    character references or attribute-value normalisation rewrites them.
    """
    return '"' + _xml_escape(str(value), {'"': "&quot;", "\n": "&#10;",
                                          "\r": "&#13;", "\t": "&#9;"}) + '"'


def _text(value):
    """Escape element text exactly once."""
    return _xml_escape(str(value))


def _hex_colour(value, default="auto"):
    """Normalise a colour to six hex digits without '#'."""
    if value is None:
        return default
    v = str(value).lstrip("#").upper()
    if not re.fullmatch(r"[0-9A-F]{6}", v):
        raise ValueError(f"colour {value!r} is not six hex digits")
    return v


# --------------------------------------------------------------- runs

class Run:
    """One stretch of text with uniform formatting."""

    def __init__(self, text, bold=False, italic=False, size_pt=None,
                 colour=None, font=None, underline=False, strike=False,
                 superscript=False, subscript=False, highlight=None):
        self.text = text
        self.bold = bold
        self.italic = italic
        self.size_pt = size_pt
        self.colour = colour
        self.font = font
        self.underline = underline
        self.strike = strike
        self.superscript = superscript
        self.subscript = subscript
        self.highlight = highlight

    def xml(self):
        p = []
        if self.font:
            p.append(f'<w:rFonts w:ascii={_attr(self.font)} '
                     f'w:hAnsi={_attr(self.font)} w:cs={_attr(self.font)}/>')
        if self.bold:
            p.append("<w:b/><w:bCs/>")
        if self.italic:
            p.append("<w:i/><w:iCs/>")
        if self.strike:
            p.append("<w:strike/>")
        if self.underline:
            p.append('<w:u w:val="single"/>')
        if self.colour:
            p.append(f'<w:color w:val={_attr(_hex_colour(self.colour))}/>')
        if self.size_pt:
            half = int(round(self.size_pt * 2))
            p.append(f'<w:sz w:val="{half}"/><w:szCs w:val="{half}"/>')
        if self.highlight:
            p.append(f'<w:highlight w:val={_attr(self.highlight)}/>')
        if self.superscript:
            p.append('<w:vertAlign w:val="superscript"/>')
        if self.subscript:
            p.append('<w:vertAlign w:val="subscript"/>')
        rpr = f"<w:rPr>{''.join(p)}</w:rPr>" if p else ""
        # A run never contains a newline: Word has no such character inside
        # <w:t>. Callers split paragraphs instead; a stray newline becomes an
        # explicit <w:br/> so nothing is silently dropped.
        body = []
        for i, chunk in enumerate(str(self.text).split("\n")):
            if i:
                body.append("<w:br/>")
            if chunk:
                body.append(f'<w:t xml:space="preserve">{_text(chunk)}</w:t>')
        return f"<w:r>{rpr}{''.join(body)}</w:r>"


def _runs(content, **kw):
    """Accept a string, a Run, a list of either, or any scalar.

    Scalars are coerced with str() rather than duck-typed as iterable. A table
    cell holding a number is ordinary -- a length, a count, a year -- and
    treating a non-string as a sequence turns that into `TypeError: 'int'
    object is not iterable` from deep inside the writer, which reads like a
    bug in the caller's data rather than a missing conversion here.
    """
    if content is None:
        return []
    if isinstance(content, Run):
        return [content]
    if isinstance(content, (list, tuple)):
        out = []
        for item in content:
            out.extend(_runs(item, **kw))
        return out
    return [Run(content if isinstance(content, str) else str(content), **kw)]


# ---------------------------------------------------------- paragraphs

_ALIGN = {"left": "left", "right": "right", "center": "center",
          "centre": "center", "justify": "both", "both": "both"}


class Paragraph:
    def __init__(self, runs=(), style=None, align=None, space_before=None,
                 space_after=None, line=None, indent=None, hanging=None,
                 keep_next=False, keep_lines=False, page_break_before=False,
                 num_id=None, num_level=0, borders=None, shading=None,
                 tabs=()):
        self.runs = list(runs)
        self.style = style
        self.align = align
        self.space_before = space_before
        self.space_after = space_after
        self.line = line
        self.indent = indent
        self.hanging = hanging
        self.keep_next = keep_next
        self.keep_lines = keep_lines
        self.page_break_before = page_break_before
        self.num_id = num_id
        self.num_level = num_level
        self.borders = borders or {}
        self.shading = shading
        self.tabs = list(tabs)

    def ppr(self):
        p = []
        if self.style:
            p.append(f'<w:pStyle w:val={_attr(self.style)}/>')
        if self.keep_next:
            p.append("<w:keepNext/>")
        if self.keep_lines:
            p.append("<w:keepLines/>")
        if self.page_break_before:
            p.append("<w:pageBreakBefore/>")
        if self.num_id is not None:
            p.append(f'<w:numPr><w:ilvl w:val="{self.num_level}"/>'
                     f'<w:numId w:val="{self.num_id}"/></w:numPr>')
        if self.borders:
            bits = []
            for side in ("top", "left", "bottom", "right"):
                b = self.borders.get(side)
                if not b:
                    continue
                size, colour = b
                bits.append(f'<w:{side} w:val="single" w:sz="{int(size)}" '
                            f'w:space="4" w:color={_attr(_hex_colour(colour))}/>')
            if bits:
                p.append(f"<w:pBdr>{''.join(bits)}</w:pBdr>")
        if self.shading:
            p.append(f'<w:shd w:val="clear" w:color="auto" '
                     f'w:fill={_attr(_hex_colour(self.shading))}/>')
        if self.tabs:
            t = "".join(
                f'<w:tab w:val={_attr(kind)} w:pos="{int(pos)}"/>'
                for pos, kind in self.tabs)
            p.append(f"<w:tabs>{t}</w:tabs>")
        sp = {}
        if self.space_before is not None:
            sp["w:before"] = int(self.space_before)
        if self.space_after is not None:
            sp["w:after"] = int(self.space_after)
        if self.line is not None:
            # lineRule must accompany line, or LibreOffice overlaps the lines.
            sp["w:line"] = int(self.line)
            sp["w:lineRule"] = "auto"
        if sp:
            p.append("<w:spacing " +
                     " ".join(f"{k}={_attr(v)}" for k, v in sp.items()) + "/>")
        ind = {}
        if self.indent is not None:
            ind["w:left"] = int(self.indent)
        if self.hanging is not None:
            ind["w:hanging"] = int(self.hanging)
        if ind:
            p.append("<w:ind " +
                     " ".join(f"{k}={_attr(v)}" for k, v in ind.items()) + "/>")
        if self.align:
            if self.align not in _ALIGN:
                raise ValueError(f"unknown alignment {self.align!r}")
            p.append(f'<w:jc w:val="{_ALIGN[self.align]}"/>')
        return f"<w:pPr>{''.join(p)}</w:pPr>" if p else ""

    def xml(self):
        return f"<w:p>{self.ppr()}{''.join(r.xml() for r in self.runs)}</w:p>"


class _RawBlock:
    def __init__(self, xml):
        self._xml = xml

    def xml(self):
        return self._xml


# -------------------------------------------------------------- tables

class Cell:
    def __init__(self, content=None, width=None, fill=None, span=1,
                 valign=None, align=None, bold=False, size_pt=None,
                 colour=None, font=None, margins=None, borders=None,
                 blocks=None):
        self.content = content
        self.width = width
        self.fill = fill
        self.span = span
        self.valign = valign
        self.align = align
        self.bold = bold
        self.size_pt = size_pt
        self.colour = colour
        self.font = font
        self.margins = margins
        self.borders = borders
        self.blocks = blocks          # pre-built blocks, overrides content

    def xml(self, default_width, default_margins, default_borders):
        p = [f'<w:tcW w:w="{int(self.width or default_width)}" w:type="dxa"/>']
        if self.span > 1:
            p.append(f'<w:gridSpan w:val="{int(self.span)}"/>')
        borders = self.borders if self.borders is not None else default_borders
        if borders:
            bits = []
            for side in ("top", "left", "bottom", "right"):
                b = borders.get(side)
                if not b:
                    continue
                size, colour = b
                bits.append(f'<w:{side} w:val="single" w:sz="{int(size)}" '
                            f'w:space="0" w:color={_attr(_hex_colour(colour))}/>')
            if bits:
                p.append(f"<w:tcBorders>{''.join(bits)}</w:tcBorders>")
        if self.fill:
            p.append(f'<w:shd w:val="clear" w:color="auto" '
                     f'w:fill={_attr(_hex_colour(self.fill))}/>')
        m = self.margins if self.margins is not None else default_margins
        if m:
            top, right, bottom, left = m
            p.append("<w:tcMar>"
                     f'<w:top w:w="{int(top)}" w:type="dxa"/>'
                     f'<w:left w:w="{int(left)}" w:type="dxa"/>'
                     f'<w:bottom w:w="{int(bottom)}" w:type="dxa"/>'
                     f'<w:right w:w="{int(right)}" w:type="dxa"/>'
                     "</w:tcMar>")
        if self.valign:
            p.append(f'<w:vAlign w:val={_attr(self.valign)}/>')
        tcpr = f"<w:tcPr>{''.join(p)}</w:tcPr>"

        if self.blocks is not None:
            body = "".join(b.xml() for b in self.blocks)
        else:
            runs = _runs(self.content, bold=self.bold, size_pt=self.size_pt,
                         colour=self.colour, font=self.font)
            body = Paragraph(runs, align=self.align, space_after=0,
                             space_before=0).xml()
        # A table cell must end with a paragraph or the file is corrupt.
        if "<w:p>" not in body and "<w:p " not in body:
            body += "<w:p/>"
        return f"<w:tc>{tcpr}{body}</w:tc>"


class Table:
    """A table whose column widths must sum to the table width.

    Word accepts a mismatch and then silently rescales every column, which is
    why both the grid and each cell carry an explicit DXA width here.
    """

    def __init__(self, col_widths, borders=None, margins=(80, 120, 80, 120),
                 indent=0, layout="fixed"):
        self.col_widths = [int(w) for w in col_widths]
        self.borders = borders if borders is not None else {
            side: (4, "CCCCCC") for side in ("top", "left", "bottom", "right")}
        self.margins = margins
        self.indent = int(indent)
        self.layout = layout
        self.rows = []

    @property
    def width(self):
        return sum(self.col_widths)

    def row(self, cells, header=False, cant_split=True, height=None,
            fill=None, bold=None, size_pt=None, colour=None, align=None,
            font=None, borders=None):
        """Append a row. `cells` may be strings, Cells, or lists of blocks."""
        built, col = [], 0
        for item in cells:
            if isinstance(item, Cell):
                c = item
            else:
                c = Cell(item)
            if c.width is None:
                span = max(1, c.span)
                c.width = sum(self.col_widths[col:col + span])
            if fill is not None and c.fill is None:
                c.fill = fill
            if bold is not None and not c.bold:
                c.bold = bold
            if size_pt is not None and c.size_pt is None:
                c.size_pt = size_pt
            if colour is not None and c.colour is None:
                c.colour = colour
            if align is not None and c.align is None:
                c.align = align
            if font is not None and c.font is None:
                c.font = font
            if borders is not None and c.borders is None:
                c.borders = borders
            built.append(c)
            col += max(1, c.span)
        if col != len(self.col_widths):
            raise ValueError(
                f"row spans {col} columns, table has {len(self.col_widths)}")
        self.rows.append((built, header, cant_split, height))
        return self

    def xml(self):
        bits = []
        for side in ("top", "left", "bottom", "right", "insideH", "insideV"):
            b = self.borders.get(side)
            if not b:
                continue
            size, colour = b
            bits.append(f'<w:{side} w:val="single" w:sz="{int(size)}" '
                        f'w:space="0" w:color={_attr(_hex_colour(colour))}/>')
        tblpr = [
            f'<w:tblW w:w="{self.width}" w:type="dxa"/>',
            f'<w:tblLayout w:type={_attr(self.layout)}/>',
        ]
        if self.indent:
            tblpr.append(f'<w:tblInd w:w="{self.indent}" w:type="dxa"/>')
        if bits:
            tblpr.append(f"<w:tblBorders>{''.join(bits)}</w:tblBorders>")
        grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in self.col_widths)

        rows = []
        for cells, header, cant_split, height in self.rows:
            trpr = []
            if cant_split:
                trpr.append("<w:cantSplit/>")
            if header:
                trpr.append("<w:tblHeader/>")
            if height:
                trpr.append(f'<w:trHeight w:val="{int(height)}" '
                            f'w:hRule="atLeast"/>')
            pr = f"<w:trPr>{''.join(trpr)}</w:trPr>" if trpr else ""
            body = "".join(
                c.xml(self.col_widths[0], self.margins, self.borders)
                for c in cells)
            rows.append(f"<w:tr>{pr}{body}</w:tr>")

        return (f"<w:tbl><w:tblPr>{''.join(tblpr)}</w:tblPr>"
                f"<w:tblGrid>{grid}</w:tblGrid>{''.join(rows)}</w:tbl>")


# ------------------------------------------------------------ document

_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/>
<Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
{extra}</Types>"""

_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>"""

_SETTINGS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:settings xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:compat><w:compatSetting w:name="compatibilityMode" w:uri="http://schemas.microsoft.com/office/word" w:val="15"/></w:compat>
<w:updateFields w:val="true"/>
</w:settings>"""

_W_NS = ('xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
         'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"')


class Document:
    """Builds a .docx. Blocks are appended in order and written on save()."""

    def __init__(self, page="A4", margin_cm=2.0, margins=None, font="Arial",
                 size_pt=10, title=None, landscape=False):
        key = page.upper() + ("_LANDSCAPE" if landscape else "")
        if key not in PAGES:
            raise ValueError(f"unknown page size {page!r}; "
                             f"choose from {sorted(PAGES)}")
        self.page_w, self.page_h = PAGES[key]
        if margins is None:
            m = cm(margin_cm)
            self.margins = (m, m, m, m)       # top, right, bottom, left
        else:
            self.margins = tuple(int(x) for x in margins)
        self.font = font
        self.size_pt = size_pt
        self.title = title
        self.blocks = []
        self._header = None
        self._footer = None
        self._heading_sizes = {1: 14, 2: 11, 3: 10}
        self._heading_colour = {1: None, 2: None, 3: None}
        self._num_ids = {}

    # -- geometry -----------------------------------------------------

    @property
    def content_width(self):
        """Usable width in DXA: the number table columns must sum to."""
        return self.page_w - self.margins[1] - self.margins[3]

    def columns(self, *fractions, gap=0):
        """Split the content width into DXA widths, summing exactly.

        Rounding is absorbed by the last column so the total always matches
        content_width, which is the thing Word checks.
        """
        total = self.content_width - gap * (len(fractions) - 1)
        scale = sum(fractions)
        out = [int(total * f / scale) for f in fractions[:-1]]
        out.append(total - sum(out))
        return out

    # -- styling ------------------------------------------------------

    def heading_style(self, level, size_pt=None, colour=None):
        if size_pt is not None:
            self._heading_sizes[level] = size_pt
        if colour is not None:
            self._heading_colour[level] = colour
        return self

    # -- blocks -------------------------------------------------------

    def add(self, block):
        self.blocks.append(block)
        return block

    def para(self, content=None, **kw):
        run_kw = {k: kw.pop(k) for k in
                  ("bold", "italic", "size_pt", "colour", "font", "underline")
                  if k in kw}
        return self.add(Paragraph(_runs(content, **run_kw), **kw))

    def heading(self, content, level=1, **kw):
        kw.setdefault("keep_next", True)
        kw.setdefault("space_before", pt(10) if level == 1 else pt(8))
        kw.setdefault("space_after", pt(4))
        colour = kw.pop("colour", self._heading_colour.get(level))
        size = kw.pop("size_pt", self._heading_sizes.get(level, 12))
        runs = _runs(content, bold=True, size_pt=size, colour=colour,
                     font=self.font)
        return self.add(Paragraph(runs, style=f"Heading{level}", **kw))

    def spacer(self, height_pt=6):
        return self.add(Paragraph([], space_after=pt(height_pt),
                                  space_before=0))

    def page_break(self):
        return self.add(_RawBlock('<w:p><w:r><w:br w:type="page"/></w:r></w:p>'))

    def rule(self, size=8, colour="B3203A", space_after=None):
        """A horizontal rule drawn as a paragraph border, never as a table."""
        return self.add(Paragraph(
            [], borders={"bottom": (size, colour)},
            space_after=space_after if space_after is not None else pt(6),
            space_before=0))

    def bullets(self, items, level=0, size_pt=None, space_after=None,
                bold=False):
        """Bulleted list via numbering definitions, never a literal glyph."""
        num_id = self._numbering("bullet")
        out = []
        for item in items:
            out.append(self.add(Paragraph(
                _runs(item, size_pt=size_pt, font=self.font, bold=bold),
                num_id=num_id, num_level=level,
                space_after=space_after if space_after is not None else pt(2),
                space_before=0, indent=pt(18) * (level + 1),
                hanging=pt(9))))
        return out

    def numbered(self, items, level=0, size_pt=None, bold=False):
        num_id = self._numbering("decimal")
        out = []
        for item in items:
            out.append(self.add(Paragraph(
                _runs(item, size_pt=size_pt, font=self.font, bold=bold),
                num_id=num_id, num_level=level, space_after=pt(2),
                space_before=0, indent=pt(18) * (level + 1), hanging=pt(9))))
        return out

    def _numbering(self, kind):
        if kind not in self._num_ids:
            self._num_ids[kind] = len(self._num_ids) + 1
        return self._num_ids[kind]

    def table(self, col_widths=None, header=None, **kw):
        widths = col_widths or [self.content_width]
        header_kw = {k: kw.pop(k) for k in ("header_fill", "header_colour",
                                            "header_bold", "header_size")
                     if k in kw}
        t = Table(widths, **kw)
        if header:
            t.row(header, header=True,
                  fill=header_kw.get("header_fill", "F4F4F4"),
                  colour=header_kw.get("header_colour"),
                  bold=header_kw.get("header_bold", True),
                  size_pt=header_kw.get("header_size", self.size_pt - 1),
                  font=self.font)
        self.add(t)
        return t

    def callout(self, title, body=None, fill="FBEFF1", bar="B3203A",
                bar_width=18, size_pt=None, items=None):
        """A tinted single-cell table with a coloured left bar.

        Built as a one-cell table because a paragraph border cannot carry a
        thick coloured edge on one side only with a fill behind it.
        """
        blocks = []
        if title:
            blocks.append(Paragraph(
                _runs(title, bold=True, size_pt=size_pt or self.size_pt,
                      font=self.font), space_after=pt(2), space_before=0))
        if body:
            blocks.append(Paragraph(
                _runs(body, size_pt=size_pt or self.size_pt, font=self.font),
                space_after=0, space_before=0, line=240))
        for item in items or []:
            blocks.append(Paragraph(
                _runs("– " + str(item), size_pt=size_pt or self.size_pt,
                      font=self.font),
                space_after=pt(1), space_before=0, indent=pt(8)))
        t = Table([self.content_width],
                  borders={"left": (bar_width, bar)},
                  margins=(100, 140, 100, 140))
        t.row([Cell(blocks=blocks, fill=fill)], cant_split=True)
        self.add(t)
        return t

    def sequence(self, seq, per_block=10, per_line=60, fill="F4F4F4",
                 highlight=None, font="Courier New", size_pt=9):
        """A nucleotide or peptide sequence in numbered blocks."""
        seq = re.sub(r"\s+", "", str(seq))
        lines = []
        for start in range(0, len(seq), per_line):
            chunk = seq[start:start + per_line]
            blocks = " ".join(chunk[i:i + per_block]
                              for i in range(0, len(chunk), per_block))
            lines.append(f"{start + 1:>6}  {blocks}")
        blocks = [Paragraph(_runs(line, font=font, size_pt=size_pt),
                            space_after=0, space_before=0)
                  for line in lines]
        t = Table([self.content_width], borders={}, margins=(100, 120, 100, 120))
        t.row([Cell(blocks=blocks, fill=fill)], cant_split=False)
        self.add(t)
        return t

    # -- header and footer --------------------------------------------

    def header(self, left=None, right=None, size_pt=None, colour="6B6B6B"):
        runs = _runs(left or "", size_pt=size_pt or self.size_pt - 2,
                     colour=colour, font=self.font)
        tabs = [(self.content_width, "right")]
        if right:
            runs.append(Run("\t", size_pt=size_pt or self.size_pt - 2))
            runs.extend(_runs(right, size_pt=size_pt or self.size_pt - 2,
                              colour=colour, font=self.font))
        align = None if right else "right"
        self._header = Paragraph(runs, tabs=tabs, align=align, space_after=0)
        return self._header

    def footer_page_x_of_y(self, left=None, size_pt=None, colour="6B6B6B"):
        """Footer with 'Page X of Y' as real fields, right-aligned."""
        size = size_pt or self.size_pt - 2
        runs = _runs(left or "", size_pt=size, colour=colour, font=self.font)
        runs.append(Run("\t", size_pt=size))
        runs.extend(_runs("Page ", size_pt=size, colour=colour, font=self.font))
        self._footer = Paragraph(runs, tabs=[(self.content_width, "right")],
                                 space_after=0)
        self._footer_fields = True
        return self._footer

    def footer(self, left=None, right=None, size_pt=None, colour="6B6B6B"):
        size = size_pt or self.size_pt - 2
        runs = _runs(left or "", size_pt=size, colour=colour, font=self.font)
        if right:
            runs.append(Run("\t", size_pt=size))
            runs.extend(_runs(right, size_pt=size, colour=colour,
                              font=self.font))
        self._footer = Paragraph(runs, tabs=[(self.content_width, "right")],
                                 space_after=0)
        self._footer_fields = False
        return self._footer

    def _field(self, instr, size_pt, colour):
        return (f'<w:fldSimple w:instr={_attr(" " + instr + " ")}>'
                f'<w:r><w:rPr>'
                f'<w:rFonts w:ascii={_attr(self.font)} w:hAnsi={_attr(self.font)}/>'
                f'<w:color w:val={_attr(_hex_colour(colour))}/>'
                f'<w:sz w:val="{int(size_pt * 2)}"/></w:rPr>'
                f'<w:t>1</w:t></w:r></w:fldSimple>')

    # -- parts ---------------------------------------------------------

    def _styles_xml(self):
        half = int(self.size_pt * 2)
        out = [f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               f"<w:styles {_W_NS}>",
               "<w:docDefaults><w:rPrDefault><w:rPr>"
               f'<w:rFonts w:ascii={_attr(self.font)} w:hAnsi={_attr(self.font)} '
               f'w:cs={_attr(self.font)}/>'
               f'<w:sz w:val="{half}"/><w:szCs w:val="{half}"/>'
               "</w:rPr></w:rPrDefault>"
               "<w:pPrDefault><w:pPr>"
               '<w:spacing w:after="120" w:line="252" w:lineRule="auto"/>'
               "</w:pPr></w:pPrDefault></w:docDefaults>",
               '<w:style w:type="paragraph" w:default="1" w:styleId="Normal">'
               '<w:name w:val="Normal"/></w:style>']
        for level in (1, 2, 3):
            size = int(self._heading_sizes.get(level, 12) * 2)
            colour = self._heading_colour.get(level)
            col = (f'<w:color w:val={_attr(_hex_colour(colour))}/>'
                   if colour else "")
            out.append(
                f'<w:style w:type="paragraph" w:styleId="Heading{level}">'
                f'<w:name w:val="heading {level}"/>'
                f'<w:basedOn w:val="Normal"/>'
                f'<w:next w:val="Normal"/>'
                f'<w:pPr><w:keepNext/><w:outlineLvl w:val="{level - 1}"/>'
                f'<w:spacing w:before="200" w:after="80" w:line="240" '
                f'w:lineRule="auto"/></w:pPr>'
                f'<w:rPr><w:b/><w:bCs/>{col}'
                f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/></w:rPr>'
                f"</w:style>")
        out.append("</w:styles>")
        return "".join(out)

    def _numbering_xml(self):
        out = [f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               f"<w:numbering {_W_NS}>"]
        kinds = self._num_ids or {}
        for kind, num_id in sorted(kinds.items(), key=lambda kv: kv[1]):
            levels = []
            for lvl in range(3):
                if kind == "bullet":
                    fmt, text, font = "bullet", "•", "Symbol"
                else:
                    fmt, text, font = "decimal", f"%{lvl + 1}.", self.font
                levels.append(
                    f'<w:lvl w:ilvl="{lvl}">'
                    f'<w:start w:val="1"/>'
                    f'<w:numFmt w:val="{fmt}"/>'
                    f'<w:lvlText w:val={_attr(text)}/>'
                    f'<w:lvlJc w:val="left"/>'
                    f'<w:pPr><w:ind w:left="{pt(18) * (lvl + 1)}" '
                    f'w:hanging="{pt(9)}"/></w:pPr>'
                    f'<w:rPr><w:rFonts w:ascii={_attr(font)} '
                    f'w:hAnsi={_attr(font)} w:hint="default"/></w:rPr>'
                    f"</w:lvl>")
            out.append(f'<w:abstractNum w:abstractNumId="{num_id}">'
                       f'<w:multiLevelType w:val="hybridMultilevel"/>'
                       f'{"".join(levels)}</w:abstractNum>')
        for kind, num_id in sorted(kinds.items(), key=lambda kv: kv[1]):
            out.append(f'<w:num w:numId="{num_id}">'
                       f'<w:abstractNumId w:val="{num_id}"/></w:num>')
        out.append("</w:numbering>")
        return "".join(out)

    def _document_xml(self, header_rel, footer_rel):
        body = "".join(b.xml() for b in self.blocks)
        refs = ""
        if header_rel:
            refs += f'<w:headerReference w:type="default" r:id={_attr(header_rel)}/>'
        if footer_rel:
            refs += f'<w:footerReference w:type="default" r:id={_attr(footer_rel)}/>'
        top, right, bottom, left = self.margins
        sect = (f"<w:sectPr>{refs}"
                f'<w:pgSz w:w="{self.page_w}" w:h="{self.page_h}"/>'
                f'<w:pgMar w:top="{top}" w:right="{right}" '
                f'w:bottom="{bottom}" w:left="{left}" '
                f'w:header="{pt(36)}" w:footer="{pt(36)}" w:gutter="0"/>'
                f"</w:sectPr>")
        return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                f"<w:document {_W_NS}><w:body>{body}{sect}</w:body></w:document>")

    def _hf_xml(self, tag, paragraph, with_fields=False):
        xml = paragraph.xml()
        if with_fields:
            size = self.size_pt - 2
            fields = (self._field("PAGE", size, "6B6B6B")
                      + f'<w:r><w:rPr><w:sz w:val="{int(size * 2)}"/>'
                        f'<w:color w:val="6B6B6B"/></w:rPr>'
                        f'<w:t xml:space="preserve"> of </w:t></w:r>'
                      + self._field("NUMPAGES", size, "6B6B6B"))
            xml = xml.replace("</w:p>", fields + "</w:p>")
        return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                f"<w:{tag} {_W_NS}>{xml}</w:{tag}>")

    def _core_xml(self):
        t = _text(self.title or "")
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<cp:coreProperties '
                'xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
                'xmlns:dc="http://purl.org/dc/elements/1.1/">'
                f"<dc:title>{t}</dc:title></cp:coreProperties>")

    def save(self, path):
        """Write the .docx. Returns the path."""
        has_header = self._header is not None
        has_footer = self._footer is not None
        header_rel = "rId10" if has_header else None
        footer_rel = "rId11" if has_footer else None

        rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering" Target="numbering.xml"/>'
                '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/settings" Target="settings.xml"/>']
        extra_types = []
        if has_header:
            rels.append('<Relationship Id="rId10" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/header" Target="header1.xml"/>')
            extra_types.append('<Override PartName="/word/header1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.header+xml"/>')
        if has_footer:
            rels.append('<Relationship Id="rId11" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Target="footer1.xml"/>')
            extra_types.append('<Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>')
        rels.append("</Relationships>")

        parts = {
            "[Content_Types].xml": _CONTENT_TYPES.format(
                extra="".join(extra_types)),
            "_rels/.rels": _RELS,
            "docProps/core.xml": self._core_xml(),
            "docProps/app.xml": ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                                 '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
                                 "<Application>docx.py</Application></Properties>"),
            "word/styles.xml": self._styles_xml(),
            "word/numbering.xml": self._numbering_xml(),
            "word/settings.xml": _SETTINGS,
            "word/_rels/document.xml.rels": "".join(rels),
            "word/document.xml": self._document_xml(header_rel, footer_rel),
        }
        if has_header:
            parts["word/header1.xml"] = self._hf_xml("hdr", self._header)
        if has_footer:
            parts["word/footer1.xml"] = self._hf_xml(
                "ftr", self._footer, getattr(self, "_footer_fields", False))

        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            # [Content_Types].xml must be the first entry in the archive.
            z.writestr("[Content_Types].xml", parts.pop("[Content_Types].xml"))
            for name, data in parts.items():
                z.writestr(name, data)
        return path


# -------------------------------------------------------------- reading

_W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


class ReadParagraph:
    def __init__(self, text, style=None, in_table=False):
        self.text = text
        self.style = style
        self.in_table = in_table

    def __repr__(self):
        return f"ReadParagraph({self.text[:40]!r}, style={self.style!r})"


def read_docx(path):
    """Paragraphs of a .docx, in document order, table cells included.

    Replaces `pandoc file.docx -t markdown` for the common case of pulling the
    wording out of an input document. Returns ReadParagraph objects carrying
    the paragraph style name, so headings can be told from body text.
    """
    import xml.etree.ElementTree as ET
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")
    root = ET.fromstring(xml)
    body = root.find(f"{_W}body")
    out = []
    if body is None:
        return out

    def paragraphs(node, in_table):
        for child in node:
            if child.tag == f"{_W}p":
                style = None
                ppr = child.find(f"{_W}pPr")
                if ppr is not None:
                    st = ppr.find(f"{_W}pStyle")
                    if st is not None:
                        style = st.get(f"{_W}val")
                text = []
                for t in child.iter():
                    if t.tag == f"{_W}t":
                        text.append(t.text or "")
                    elif t.tag in (f"{_W}br", f"{_W}cr"):
                        text.append("\n")
                    elif t.tag == f"{_W}tab":
                        text.append("\t")
                out.append(ReadParagraph("".join(text), style, in_table))
            elif child.tag == f"{_W}tbl":
                for row in child.findall(f"{_W}tr"):
                    for cellnode in row.findall(f"{_W}tc"):
                        paragraphs(cellnode, True)
            elif child.tag in (f"{_W}sdt", f"{_W}sdtContent"):
                paragraphs(child, in_table)

    paragraphs(body, False)
    return out


def docx_text(path, keep_empty=False):
    """The whole document as plain text, one paragraph per line."""
    lines = [p.text for p in read_docx(path)]
    if not keep_empty:
        lines = [l for l in lines if l.strip()]
    return "\n".join(lines)


# ------------------------------------------------------------ checking

_REQUIRED = ["[Content_Types].xml", "_rels/.rels", "word/document.xml",
             "word/_rels/document.xml.rels", "word/styles.xml"]


def validate_docx(path):
    """Structural check of a written .docx.

    Reports rather than raising, so a build script can print every problem at
    once. This is not a substitute for looking at the rendered pages -- it
    cannot see a clipped table or an orphaned callout -- but it does catch the
    corruption that stops the file opening at all.
    """
    import xml.etree.ElementTree as ET
    problems = []
    checks = {}
    try:
        with zipfile.ZipFile(path) as z:
            bad = z.testzip()
            if bad:
                problems.append(f"corrupt archive member: {bad}")
            names = z.namelist()
            checks["parts"] = len(names)
            if names and names[0] != "[Content_Types].xml":
                problems.append("[Content_Types].xml is not the first entry")
            for need in _REQUIRED:
                if need not in names:
                    problems.append(f"missing part: {need}")
            for name in names:
                if not name.endswith((".xml", ".rels")):
                    continue
                try:
                    ET.fromstring(z.read(name))
                except ET.ParseError as e:
                    problems.append(f"{name}: not well-formed XML ({e})")
            # every r:id used in document.xml must exist in its rels part
            if ("word/document.xml" in names
                    and "word/_rels/document.xml.rels" in names):
                doc = z.read("word/document.xml").decode("utf-8", "replace")
                relxml = ET.fromstring(z.read("word/_rels/document.xml.rels"))
                have = {r.get("Id") for r in relxml}
                used = set(re.findall(r'r:id="([^"]+)"', doc))
                checks["relationships_used"] = len(used)
                for rid in sorted(used - have):
                    problems.append(f"document.xml uses {rid}, "
                                    f"which is not in document.xml.rels")
                for target in (r.get("Target") for r in relxml):
                    if target and not target.startswith(("http", "mailto")) \
                            and f"word/{target}" not in names:
                        problems.append(f"relationship target missing: {target}")
            # table geometry: grid total must equal the declared table width
            if "word/document.xml" in names:
                doc = z.read("word/document.xml").decode("utf-8", "replace")
                tables = re.findall(
                    r'<w:tbl>.*?<w:tblW w:w="(\d+)".*?<w:tblGrid>(.*?)</w:tblGrid>',
                    doc, re.S)
                checks["tables"] = len(tables)
                for i, (declared, grid) in enumerate(tables, 1):
                    cols = [int(x) for x in
                            re.findall(r'<w:gridCol w:w="(\d+)"/>', grid)]
                    if cols and sum(cols) != int(declared):
                        problems.append(
                            f"table {i}: columns sum to {sum(cols)} but "
                            f"tblW is {declared}; Word will rescale it")
                empty = re.findall(r"<w:tc>(?:(?!</w:tc>).)*?</w:tc>", doc, re.S)
                checks["cells"] = len(empty)
                for cellxml in empty:
                    if "<w:p" not in cellxml:
                        problems.append("a table cell has no paragraph, "
                                        "which Word treats as corrupt")
                        break
            if "word/document.xml" in names:
                doc = z.read("word/document.xml").decode("utf-8", "replace")
                for m in re.finditer(r'<w:spacing([^/]*)/>', doc):
                    if "w:line=" in m.group(1) and "w:lineRule=" not in m.group(1):
                        problems.append("a w:spacing sets w:line without "
                                        "w:lineRule; LibreOffice will overlap "
                                        "the lines")
                        break
    except zipfile.BadZipFile as e:
        problems.append(f"not a zip file: {e}")
    except FileNotFoundError:
        problems.append("file does not exist")
    return {"ok": not problems, "problems": problems, "checks": checks}


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    target = sys.argv[1]
    report = validate_docx(target)
    for k, v in report["checks"].items():
        print(f"  {k}: {v}")
    for p in report["problems"]:
        print(f"  PROBLEM: {p}")
    print("ok" if report["ok"] else "FAILED")
    raise SystemExit(0 if report["ok"] else 1)
