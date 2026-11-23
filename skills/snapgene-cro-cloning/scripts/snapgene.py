"""SnapGene .dna cloning helpers: read, edit, annotate, verify, and order-sheet export.

Python standard library only. No Biopython, no openpyxl, no autosnapgene,
no snapgene-reader. Import it with::

    import sys; sys.path.insert(0, "<skill>/scripts")
    import snapgene as sg

The .dna reader and writer are implemented directly against the file format
(1-byte block type, 4-byte big-endian length, payload). The writer emits only
the blocks it models -- header, sequence, features, notes -- so a file it
writes carries no stale cached state from the vector it was derived from.
That is deliberate: see WRITER_NOTE.
"""

import re
import struct
import xml.etree.ElementTree as ET
import zipfile
from xml.sax.saxutils import escape

WRITER_NOTE = """\
Why this writer rebuilds rather than patches.

A .dna file is a flat sequence of typed blocks. SnapGene itself writes
several blocks that cache absolute nucleotide positions computed against the
sequence as it stood when the file was saved: block 2 (a packed position
overlay) and the block-11 HistoryNode chain are the two that matter. A
library that loads a file, edits the sequence and writes it back copies every
block it does not model through as raw bytes -- stale caches included. The
sequence, the features and the translations are then all correct, while
SnapGene desktop can still draw a restriction site tens to hundreds of bases
away from its true position, because it reads one of those caches when
rendering.

This module never carries an unmodelled block across. write_dna() emits a
fresh header, the sequence, the features and the notes, and nothing else --
the same thing SnapGene produces when you paste a sequence into a new file,
plus the annotation. There is no stale state to strip because none is kept,
so no "rebuild clean" step is needed before writing.

What is preserved: sequence, topology and the strandedness/methylation
flags, every feature with its segments, colours, directionality, translation
and numbering attributes and all qualifiers, and the file-level notes
including the description shown in the info pane.

What is dropped: the cached restriction-site list, the position overlay, edit
history, primers, alignments, enzyme visibility and custom enzyme sets. Of
these only primers are user data; if the vector carries primers you need,
copy them in SnapGene after the construct is built, or keep the original file
alongside.
"""

# ---------------------------------------------------------------- constants

BLOCK_HEADER = 9
BLOCK_SEQUENCE = 0
BLOCK_FEATURES = 10
BLOCK_PRIMERS = 5
BLOCK_NOTES = 6

# Sequence-block flag bits (byte 0 of the sequence block payload).
FLAG_CIRCULAR = 0x01
FLAG_DOUBLE_STRANDED = 0x02
FLAG_DAM_METHYLATED = 0x04
FLAG_DCM_METHYLATED = 0x08
FLAG_ECOKI_METHYLATED = 0x10

IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "R": "Y", "Y": "R",
    "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}

# A census panel of common cloning enzymes. Recognition sequences only --
# this module does not model cut offsets, overhangs or star activity.
ENZYME_SITES = {
    "AgeI": "ACCGGT", "AscI": "GGCGCGCC", "AvrII": "CCTAGG",
    "BamHI": "GGATCC", "BbsI": "GAAGAC", "BclI": "TGATCA",
    "BglII": "AGATCT", "BmtI": "GCTAGC", "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC", "BspEI": "TCCGGA", "BsrGI": "TGTACA",
    "BstBI": "TTCGAA", "ClaI": "ATCGAT", "DraIII": "CACNNNGTG",
    "EagI": "CGGCCG", "EcoRI": "GAATTC", "EcoRV": "GATATC",
    "FseI": "GGCCGGCC", "HindIII": "AAGCTT", "KpnI": "GGTACC",
    "MfeI": "CAATTG", "MluI": "ACGCGT", "NcoI": "CCATGG",
    "NdeI": "CATATG", "NheI": "GCTAGC", "NotI": "GCGGCCGC",
    "NsiI": "ATGCAT", "PacI": "TTAATTAA", "PciI": "ACATGT",
    "PmeI": "GTTTAAAC", "PmlI": "CACGTG", "PstI": "CTGCAG",
    "PvuI": "CGATCG", "SacI": "GAGCTC", "SacII": "CCGCGG",
    "SalI": "GTCGAC", "SbfI": "CCTGCAGG", "ScaI": "AGTACT",
    "SmaI": "CCCGGG", "SpeI": "ACTAGT", "SphI": "GCATGC",
    "StuI": "AGGCCT", "SwaI": "ATTTAAAT", "XbaI": "TCTAGA",
    "XhoI": "CTCGAG", "XmaI": "CCCGGG",
}

# Enzymes whose cutting is blocked or impaired by overlapping Dam (GATC) or
# Dcm (CCWGG) methylation in a standard cloning strain. The real answer is
# enzyme- and context-specific: treat a hit as a candidate to confirm against
# REBASE, not as a verdict.
DAM_SENSITIVE = ("BclI", "ClaI", "MboI", "XbaI", "TaqI", "HphI", "MboII", "NruI", "NarI")
DCM_SENSITIVE = ("StuI", "EaeI", "AlwI", "ApaI", "BsaBI", "NaeI", "ScrFI")

DB_URL_TEMPLATES = {
    "UNIPROT": "https://www.uniprot.org/uniprotkb/{acc}",
    "GENBANK": "https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
    "PDB": "https://www.rcsb.org/structure/{acc}",
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Average isotopic residue masses (Da), for a mass cross-check against a
# source record header. Not monoisotopic.
_AA_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
    "U": 150.0379, "O": 237.3018,
}
_WATER = 18.0153

# Codons ranked by E. coli usage, most used first, for suggesting a
# synonymous swap that also does not worsen codon usage.
ECOLI_PREFERRED = {
    "A": ["GCG", "GCC", "GCA", "GCT"], "R": ["CGC", "CGT", "CGG", "AGA", "AGG", "CGA"],
    "N": ["AAC", "AAT"], "D": ["GAT", "GAC"], "C": ["TGC", "TGT"],
    "Q": ["CAG", "CAA"], "E": ["GAA", "GAG"], "G": ["GGC", "GGT", "GGG", "GGA"],
    "H": ["CAT", "CAC"], "I": ["ATC", "ATT", "ATA"],
    "L": ["CTG", "TTA", "TTG", "CTC", "CTT", "CTA"], "K": ["AAA", "AAG"],
    "M": ["ATG"], "F": ["TTT", "TTC"], "P": ["CCG", "CCA", "CCT", "CCC"],
    "S": ["AGC", "TCT", "TCC", "AGT", "TCG", "TCA"],
    "T": ["ACC", "ACG", "ACT", "ACA"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTG", "GTT", "GTC", "GTA"], "*": ["TAA", "TGA", "TAG"],
}

# Fallback order-sheet columns. A template the human supplies overrides this.
ORDER_SHEET_COLUMNS = [
    "No", "Final plasmid Name", "Final plasmid sequence", "Insert sequence",
    "Insert Source", "Cloning sites (5' and 3')", "5'homologous sequence (~30 bp)",
    "3' hom (~30 bp)", "Vector", "Vector source", "Vector sequence",
    "Shippable", "Maxiprep",
]

ORDER_SHEET_RED_COLUMNS = (
    "Cloning sites (5' and 3')", "5'homologous sequence (~30 bp)", "3' hom (~30 bp)",
)

XLSX_CELL_LIMIT = 32767


# ------------------------------------------------------------ block codec

def read_blocks(path):
    """Yield (block_type, payload_bytes) for every block in a .dna file."""
    with open(path, "rb") as fh:
        data = fh.read()
    out, i = [], 0
    while i < len(data):
        if i + 5 > len(data):
            raise ValueError(f"truncated block header at byte {i}")
        btype = data[i]
        (length,) = struct.unpack(">I", data[i + 1:i + 5])
        payload = data[i + 5:i + 5 + length]
        if len(payload) != length:
            raise ValueError(f"truncated block {btype} at byte {i}: "
                             f"want {length} bytes, got {len(payload)}")
        out.append((btype, payload))
        i += 5 + length
    return out


def _pack_block(btype, payload):
    return bytes([btype]) + struct.pack(">I", len(payload)) + payload


# -------------------------------------------------------------- data model

class Segment:
    """One drawn span of a feature. `start`/`end` are 0-based half-open and
    slice the sequence directly: seq[start:end] is the annotated span.

    A segment on a circular plasmid may span the origin, in which case
    `end <= start` and `spans_origin` is True -- the span is then
    seq[start:] + seq[:end], not the empty slice Python would give. Use
    `subseq(seq)` rather than slicing by hand.
    """

    def __init__(self, start, end, color="#99ccff", seg_type="standard",
                 translated=False, name=None):
        self.start = int(start)
        self.end = int(end)
        self.color = color
        self.type = seg_type
        self.translated = bool(translated)
        self.name = name

    @property
    def spans_origin(self):
        return self.end <= self.start

    def subseq(self, seq):
        """The nucleotides this segment covers, origin wrap included."""
        if self.spans_origin:
            return seq[self.start:] + seq[:self.end]
        return seq[self.start:self.end]

    def length(self, seq_len):
        """Length in bases; needs the sequence length to resolve a wrap."""
        if self.spans_origin:
            return seq_len - self.start + self.end
        return self.end - self.start

    def __len__(self):
        if self.spans_origin:
            raise ValueError(
                f"segment {self.start}-{self.end} spans the origin; "
                f"use .length(len(seq)) instead of len()")
        return self.end - self.start

    def __repr__(self):
        wrap = ", spans_origin" if self.spans_origin else ""
        return f"Segment({self.start}-{self.end}, translated={self.translated}{wrap})"


class Feature:
    """A SnapGene feature.

    `qualifiers` maps name -> value. Values are plain Python: ints stay ints
    (codon_start, transl_table) and strings stay strings, and the writer
    serialises each to the matching XML value type. SnapGene stores most
    descriptive qualifier text as an HTML fragment; this module hands it back
    as stored, and make_product_html() builds the accession-linked form.
    """

    def __init__(self, name, ftype="misc_feature", segments=None,
                 directionality=None, qualifiers=None, attrs=None):
        self.name = name
        self.type = ftype
        self.segments = list(segments or [])
        self.directionality = directionality
        self.qualifiers = dict(qualifiers or {})
        # Attributes carried through verbatim so a round-trip is lossless.
        self.attrs = dict(attrs or {})

    @property
    def spans_origin(self):
        return any(s.spans_origin for s in self.segments)

    @property
    def start(self):
        """Lowest 0-based coordinate. For a feature spanning the origin this
        is the start of the wrapping segment, and `start < end` does not
        hold -- test `spans_origin` before treating the pair as a slice."""
        if self.spans_origin:
            return min(s.start for s in self.segments if s.spans_origin)
        return min(s.start for s in self.segments)

    @property
    def end(self):
        if self.spans_origin:
            return max(s.end for s in self.segments if s.spans_origin)
        return max(s.end for s in self.segments)

    @property
    def numbering_start(self):
        """The residue number SnapGene labels the first translated residue
        with, or None. The file attribute is consecutiveNumberingStartsFrom."""
        v = self.attrs.get("consecutiveNumberingStartsFrom")
        return int(v) if v is not None else None

    @numbering_start.setter
    def numbering_start(self, value):
        if value is None:
            self.attrs.pop("consecutiveNumberingStartsFrom", None)
        else:
            self.attrs["consecutiveNumberingStartsFrom"] = str(int(value))
            # Present on every feature in real files; set it so a feature
            # built from scratch matches what SnapGene writes.
            self.attrs.setdefault("consecutiveTranslationNumbering", "1")

    @property
    def translation(self):
        return self.qualifiers.get("translation")

    def length(self, seq_len):
        """Total annotated length in bases, origin wrap included."""
        return sum(s.length(seq_len) for s in self.segments)

    def __repr__(self):
        return (f"Feature({self.name!r}, {self.type}, "
                f"{self.start}-{self.end}, {len(self.segments)} segment(s))")


class Dna:
    """A parsed .dna file."""

    def __init__(self, seq="", circular=True, double_stranded=True,
                 dam=False, dcm=False, ecoki=False, features=None,
                 notes=None, export_version=15, import_version=19):
        self.seq = seq.upper()
        self.circular = circular
        self.double_stranded = double_stranded
        self.dam = dam
        self.dcm = dcm
        self.ecoki = ecoki
        self.features = list(features or [])
        # Notes block children, name -> text. "Description" is the info-pane
        # summary; "Type", "UUID", "Created" and others also live here.
        self.notes = dict(notes or {})
        self.export_version = export_version
        self.import_version = import_version

    @property
    def description(self):
        """The file-level description shown in SnapGene's info pane."""
        return self.notes.get("Description")

    @description.setter
    def description(self, value):
        self.notes["Description"] = value

    def feature(self, name, ftype=None, exact=True):
        """Return the first feature matching `name` (and optionally type).
        With exact=False, matches case-insensitively on a substring."""
        for f in self.features:
            hit = (f.name == name) if exact else (name.lower() in (f.name or "").lower())
            if hit and (ftype is None or f.type == ftype):
                return f
        raise KeyError(f"feature {name!r} (type={ftype}) not found; "
                       f"have {[f.name for f in self.features]}")

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        topo = "circular" if self.circular else "linear"
        return f"Dna({len(self.seq)} bp, {topo}, {len(self.features)} features)"


# ----------------------------------------------------------------- reading

def _parse_features_xml(xml_text):
    root = ET.fromstring(xml_text)
    features = []
    for fe in root.findall("Feature"):
        attrs = dict(fe.attrib)
        name = attrs.pop("name", "")
        ftype = attrs.pop("type", "misc_feature")
        direction = attrs.pop("directionality", None)
        segments = []
        for se in fe.findall("Segment"):
            rng = se.get("range", "")
            m = re.match(r"^(\d+)-(\d+)$", rng.strip())
            if not m:
                continue
            # SnapGene ranges are 1-based inclusive.
            segments.append(Segment(
                start=int(m.group(1)) - 1,
                end=int(m.group(2)),
                color=se.get("color", "#99ccff"),
                seg_type=se.get("type", "standard"),
                translated=se.get("translated") == "1",
                name=se.get("name"),
            ))
        qualifiers = {}
        for qe in fe.findall("Q"):
            qname = qe.get("name")
            ve = qe.find("V")
            if ve is None or qname is None:
                continue
            if ve.get("int") is not None:
                qualifiers[qname] = int(ve.get("int"))
            elif ve.get("text") is not None:
                # ElementTree has already decoded the XML entities, and that is
                # the only decode there is: SnapGene stores a qualifier as HTML
                # source escaped once into an XML attribute, so the decoded
                # value is HTML markup ("<html><body>..."), not display text.
                # Decoding a second time here would turn a literal "&amp;" into
                # a bare "&" and, paired with a second escape on write, hide a
                # double-escaping bug from any read-write-read test.
                qualifiers[qname] = ve.get("text")
            elif ve.get("predef") is not None:
                qualifiers[qname] = ve.get("predef")
        features.append(Feature(name, ftype, segments, direction, qualifiers, attrs))
    return features


def _parse_notes_xml(xml_text):
    root = ET.fromstring(xml_text)
    notes = {}
    for child in root:
        notes[child.tag] = (child.text or "").strip()
    return notes


def read_dna(path):
    """Parse a SnapGene .dna file into a Dna object.

    Feature coordinates come back 0-based half-open and slice the sequence
    directly -- seq[f.start:f.end] reproduces the annotated span. Do not
    subtract 1 to "correct" for 1-based indexing; the conversion is already
    done here, and doing it twice shifts the reading frame and yields a
    plausible-looking but wrong translation.
    """
    dna = Dna(seq="", features=[], notes={})
    seen_sequence = False
    for btype, payload in read_blocks(path):
        if btype == BLOCK_HEADER:
            if not payload.startswith(b"SnapGene"):
                raise ValueError(f"{path}: not a SnapGene file "
                                 f"(header block starts {payload[:8]!r})")
            if len(payload) >= 14:
                ftype, exp, imp = struct.unpack(">HHH", payload[8:14])
                if ftype != 1:
                    raise ValueError(f"{path}: file type {ftype}, expected 1 (DNA)")
                dna.export_version, dna.import_version = exp, imp
        elif btype == BLOCK_SEQUENCE:
            flags = payload[0]
            dna.seq = payload[1:].decode("ascii").upper()
            dna.circular = bool(flags & FLAG_CIRCULAR)
            dna.double_stranded = bool(flags & FLAG_DOUBLE_STRANDED)
            dna.dam = bool(flags & FLAG_DAM_METHYLATED)
            dna.dcm = bool(flags & FLAG_DCM_METHYLATED)
            dna.ecoki = bool(flags & FLAG_ECOKI_METHYLATED)
            seen_sequence = True
        elif btype == BLOCK_FEATURES:
            dna.features = _parse_features_xml(payload.decode("utf-8"))
        elif btype == BLOCK_NOTES:
            dna.notes = _parse_notes_xml(payload.decode("utf-8"))
    if not seen_sequence:
        raise ValueError(f"{path}: no sequence block")
    return dna


# ----------------------------------------------------------------- writing

def _attr(value):
    """Quote an XML attribute value the way SnapGene does.

    SnapGene always delimits with double quotes and escapes an embedded one as
    &quot;. quoteattr() instead switches the delimiter to single quotes, which
    is equivalent XML but does not match SnapGene's own bytes -- and SnapGene
    parses this XML with its own reader, not a conformant one, so matching its
    convention is cheaper than assuming the alternative is accepted.

    Newline, carriage return and tab must be written as character references,
    which is the other thing quoteattr() does and the reason it cannot simply
    be replaced by quoting. An XML parser applies attribute-value
    normalisation and silently rewrites a literal newline inside an attribute
    to a space, so a feature label carrying a line break -- SnapGene writes
    those, and they survive in real files -- comes back one space per break
    short. The loss is invisible in the written file and only shows up on a
    re-read.
    """
    return '"' + escape(str(value), {'"': "&quot;", "\n": "&#10;",
                                     "\r": "&#13;", "\t": "&#9;"}) + '"'


def _features_xml(features):
    parts = ['<?xml version="1.0"?>',
             f'<Features nextValidID="{len(features)}">']
    for i, f in enumerate(features):
        attrs = dict(f.attrs)
        attrs.setdefault("allowSegmentOverlaps", "0")
        attrs.setdefault("consecutiveTranslationNumbering", "1")
        bits = [f'recentID="{i}"', f"name={_attr(f.name or '')}",
                f"type={_attr(f.type or 'misc_feature')}"]
        if f.directionality is not None:
            bits.append(f"directionality={_attr(str(f.directionality))}")
        for k in sorted(attrs):
            if k in ("recentID", "name", "type", "directionality"):
                continue
            bits.append(f"{k}={_attr(str(attrs[k]))}")
        parts.append("<Feature " + " ".join(bits) + ">")
        for s in f.segments:
            sbits = [f'range="{s.start + 1}-{s.end}"',
                     f"color={_attr(s.color or '#99ccff')}",
                     f"type={_attr(s.type or 'standard')}"]
            if s.translated:
                sbits.append('translated="1"')
            if s.name:
                sbits.append(f"name={_attr(s.name)}")
            parts.append("<Segment " + " ".join(sbits) + "/>")
        for qname in sorted(f.qualifiers):
            value = f.qualifiers[qname]
            parts.append(f"<Q name={_attr(qname)}>")
            if isinstance(value, bool):
                raise TypeError(f"qualifier {qname!r}: bool is not a SnapGene "
                                f"value type; use int 1/0 or a string")
            if isinstance(value, int):
                parts.append(f'<V int="{value}"/>')
            else:
                # _attr() escapes as well as quotes, so it is the whole
                # XML-level escape. Wrapping it around escape() would encode
                # twice: SnapGene decodes once, so "<html>" would reach the
                # feature inspector as literal "&lt;html&gt;" text and render
                # as visible tags instead of formatted HTML.
                parts.append(f"<V text={_attr(str(value))}/>")
            parts.append("</Q>")
        parts.append("</Feature>")
    parts.append("</Features>")
    return "".join(parts)


def _notes_xml(notes):
    parts = ["<Notes>"]
    for key in notes:
        value = notes[key]
        if value is None:
            continue
        parts.append(f"<{key}>{escape(str(value))}</{key}>")
    parts.append("</Notes>")
    return "".join(parts)


def write_dna(dna, path):
    """Write a Dna object to a .dna file SnapGene can open.

    Emits the header, sequence, features and notes blocks and nothing else.
    See WRITER_NOTE for why that matters: no unmodelled block is carried
    across, so no stale position cache can survive a sequence edit, and no
    separate "rebuild clean" pass is needed before writing.
    """
    if not dna.seq:
        raise ValueError("refusing to write an empty sequence")
    bad = set(dna.seq) - set("ACGTRYSWKMBDHVN")
    if bad:
        raise ValueError(f"sequence contains non-IUPAC characters: {sorted(bad)}")
    for f in dna.features:
        for s in f.segments:
            if not (0 <= s.start <= len(dna.seq) and 0 <= s.end <= len(dna.seq)):
                raise ValueError(f"feature {f.name!r} segment {s.start}-{s.end} "
                                 f"is outside the {len(dna.seq)} bp sequence")
            if s.spans_origin and not dna.circular:
                raise ValueError(f"feature {f.name!r} segment {s.start}-{s.end} "
                                 f"spans the origin of a linear sequence")

    flags = 0
    if dna.circular:
        flags |= FLAG_CIRCULAR
    if dna.double_stranded:
        flags |= FLAG_DOUBLE_STRANDED
    if dna.dam:
        flags |= FLAG_DAM_METHYLATED
    if dna.dcm:
        flags |= FLAG_DCM_METHYLATED
    if dna.ecoki:
        flags |= FLAG_ECOKI_METHYLATED

    header = b"SnapGene" + struct.pack(">HHH", 1, dna.export_version,
                                      dna.import_version)
    chunks = [
        _pack_block(BLOCK_HEADER, header),
        _pack_block(BLOCK_SEQUENCE, bytes([flags]) + dna.seq.encode("ascii")),
        _pack_block(BLOCK_FEATURES, _features_xml(dna.features).encode("utf-8")),
    ]
    # Only write a notes block when there is something to put in it. Not every
    # SnapGene file has one, so defaulting a Type into an empty one invents
    # metadata the source never carried and stops the file round-tripping.
    if dna.notes:
        notes = dict(dna.notes)
        notes.setdefault("Type", "Synthetic")
        chunks.append(_pack_block(BLOCK_NOTES,
                                  _notes_xml(notes).encode("utf-8")))
    with open(path, "wb") as fh:
        fh.write(b"".join(chunks))
    return path


# --------------------------------------------------------- sequence basics

def revcomp(seq):
    """Reverse complement, IUPAC-aware."""
    return "".join(IUPAC_COMPLEMENT.get(c, "N") for c in reversed(seq.upper()))


def translate(dna_seq, to_stop=False):
    """Translate a DNA string in frame 1. Incomplete trailing codons are
    dropped; unknown codons become 'X'."""
    seq = dna_seq.upper().replace("U", "T")
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and to_stop:
            break
        out.append(aa)
    return "".join(out)


def translate_and_validate(dna_seq, expected_protein=None, allow_trailing_stop=True):
    """Translate and report frame integrity against an expected protein.

    Returns a dict of measurements -- it asserts nothing. `match` is None
    when no expected_protein was given, so `if result["match"]` and
    `if result["match"] is not False` mean different things; test explicitly.
    """
    seq = dna_seq.upper().replace("U", "T")
    protein = translate(seq)
    trailing_stop = protein.endswith("*")
    body = protein[:-1] if trailing_stop else protein
    internal_stops = [i + 1 for i, a in enumerate(body) if a == "*"]

    result = {
        "protein": body,
        "length_nt": len(seq),
        "length_aa": len(body),
        "in_frame": len(seq) % 3 == 0,
        "trailing_stop": trailing_stop,
        "internal_stops": internal_stops,
        "starts_with_met": body.startswith("M"),
        "match": None,
        "first_mismatch": None,
    }
    if not allow_trailing_stop and trailing_stop:
        result["internal_stops"] = internal_stops + [len(body) + 1]

    if expected_protein is not None:
        exp = expected_protein.upper().rstrip("*")
        result["match"] = (body == exp)
        if body != exp:
            for i, (a, b) in enumerate(zip(body, exp)):
                if a != b:
                    result["first_mismatch"] = {
                        "residue": i + 1, "got": a, "expected": b,
                    }
                    break
            else:
                result["first_mismatch"] = {
                    "residue": min(len(body), len(exp)) + 1,
                    "got": f"length {len(body)}",
                    "expected": f"length {len(exp)}",
                }
    return result


def feature_segments_in_order(feature):
    """A feature's segments in biological 5'->3' order along the feature's own
    direction: ascending for a forward feature, descending for a reverse one."""
    reverse = str(feature.directionality) == "2"
    return sorted(feature.segments, key=lambda s: s.start, reverse=reverse)


def feature_span(dna_or_seq, feature):
    """The feature's nucleotide sequence, segments spliced in biological order
    and reverse-complemented for a reverse feature."""
    seq = dna_or_seq.seq if hasattr(dna_or_seq, "seq") else dna_or_seq
    reverse = str(feature.directionality) == "2"
    parts = [s.subseq(seq) for s in feature_segments_in_order(feature)]
    return "".join(revcomp(p) if reverse else p for p in parts)


def feature_translation(dna_or_seq, feature):
    """The feature's peptide in SnapGene's own cached-qualifier convention.

    Two conventions matter and neither is obvious; both were established by
    reading real files rather than from documentation:

    - A terminal stop codon is KEPT as a trailing '*'. Stripping it makes
      every stop-terminated CDS look like a stale cache.
    - A multi-segment feature is translated one segment at a time and the
      pieces are joined with a comma, so a two-segment CDS reads
      'MSIQ...VFA,HPETL...W*'. Translating the spliced span as one string
      gives the same residues without the separator, and comparing that
      against the cached value reports a false mismatch.
    """
    seq = dna_or_seq.seq if hasattr(dna_or_seq, "seq") else dna_or_seq
    reverse = str(feature.directionality) == "2"
    pieces = []
    for s in feature_segments_in_order(feature):
        part = s.subseq(seq)
        pieces.append(translate(revcomp(part) if reverse else part))
    return ",".join(pieces)


def protein_of(dna_or_seq, feature):
    """The feature's protein as a plain sequence: no segment separators and no
    trailing stop, for comparison against a database record."""
    return feature_translation(dna_or_seq, feature).replace(",", "").rstrip("*")


def compute_protein_mw(protein_seq):
    """Average molecular mass in Da, for a cross-check against a value quoted
    in a source record. Residue masses are average, not monoisotopic."""
    seq = protein_seq.upper().rstrip("*")
    return round(sum(_AA_MASS.get(a, 0.0) for a in seq) + _WATER, 2)


def parse_protein_dna_fasta_pairs(path):
    """Parse a file of alternating protein and DNA FASTA records into
    {name: {'protein': str, 'dna': str}}.

    A record is classified by composition, not by order, so a file that
    lists DNA first still parses. Names are taken from the protein record
    where both are present.
    """
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.replace(" ", ""))
    if name is not None:
        records.append((name, "".join(chunks)))

    out = {}
    for rec_name, seq in records:
        seq = seq.upper()
        letters = set(seq) - {"*", "-"}
        is_dna = letters and letters <= set("ACGTUN")
        key = re.sub(r"[_\s-]*(protein|aa|dna|nt|nucleotide)$", "", rec_name,
                     flags=re.I).strip() or rec_name
        slot = out.setdefault(key, {"protein": None, "dna": None})
        slot["dna" if is_dna else "protein"] = seq
    return out


# ------------------------------------------------------------ construct ops

def replace_orf(vector_seq, orf_start, orf_end, insert_dna):
    """Splice insert_dna in place of vector_seq[orf_start:orf_end]."""
    if not 0 <= orf_start <= orf_end <= len(vector_seq):
        raise ValueError(f"ORF {orf_start}-{orf_end} outside "
                         f"{len(vector_seq)} bp vector")
    return vector_seq[:orf_start] + insert_dna.upper() + vector_seq[orf_end:]


def shift_features(dna, edit_start, delta):
    """Shift every feature segment that lies downstream of edit_start by
    delta bases, in place. Segments spanning the edit are returned in a
    list of (feature, segment) pairs for the caller to resolve -- they
    cannot be shifted correctly without knowing what the edit meant."""
    spanning = []
    for f in dna.features:
        for s in f.segments:
            if s.start >= edit_start:
                s.start += delta
                s.end += delta
            elif s.end > edit_start:
                spanning.append((f, s))
    return spanning


def seamless_junction_arms(vector_seq, orf_start, orf_end, insert_dna, arm_len=18):
    """Homology arms for seamless assembly: (5' arm, 3' arm), each arm_len
    bases of vector flanking the insertion point, oriented 5'->3'."""
    five = vector_seq[max(0, orf_start - arm_len):orf_start]
    three = vector_seq[orf_end:orf_end + arm_len]
    return five, three


def restriction_synthesis_fragment(vector_seq, five_prime_start, orf_start,
                                   orf_end, insert_dna, three_prime_end):
    """The fragment to submit for synthesis, running from five_prime_start to
    three_prime_end in vector coordinates with the insert spliced in.

    The returned fragment INCLUDES the flanking restriction sites: the common
    mistake is to submit the coding sequence alone and leave the enzyme sites
    out of the synthesised DNA, where they must physically be for the cut to
    work.
    """
    if not five_prime_start <= orf_start <= orf_end <= three_prime_end:
        raise ValueError(f"boundaries out of order: {five_prime_start} <= "
                         f"{orf_start} <= {orf_end} <= {three_prime_end}")
    return (vector_seq[five_prime_start:orf_start] + insert_dna.upper()
            + vector_seq[orf_end:three_prime_end])


# -------------------------------------------------------------- site census

def iupac_to_regex(site):
    """Turn a literal IUPAC recognition sequence into a regex. The input is a
    recognition sequence, NOT a regex -- 'GGTCTC', not 'GG(T|C)CTC'."""
    try:
        return "".join(IUPAC_REGEX[c] for c in site.upper())
    except KeyError as exc:
        raise ValueError(f"{site!r} is not an IUPAC sequence: {exc}") from None


def find_sites(seq, site, circular=True, both_strands=True):
    """0-based positions of every occurrence of a recognition sequence.

    Scans both strands and wraps the circular junction by default, and both
    matter. Every Type IIS enzyme (BsaI GGTCTC, BsmBI CGTCTC, BbsI GAAGAC) is
    non-palindromic, so a forward-only scan reports a contaminated insert as
    clean; and a site straddling the origin is invisible to a linear scan.
    Palindromic sites de-duplicate, so nothing is double-counted.

    Pass circular=False for a linear fragment.
    """
    seq = seq.upper()
    pattern = iupac_to_regex(site)
    targets = [pattern]
    if both_strands:
        rc = iupac_to_regex(revcomp(site))
        if rc != pattern:
            targets.append(rc)

    scan = seq
    if circular and len(seq) > len(site):
        scan = seq + seq[:len(site) - 1]

    hits = set()
    for pat in targets:
        for m in re.finditer(f"(?={pat})", scan):
            hits.add(m.start() % len(seq))
    return sorted(hits)


def check_sites_in_construct(construct_seq, enzymes=None, circular=True):
    """Census enzymes across a whole construct.

    Returns {name: {'positions': [...], 'count': n, 'unique': bool}}. For a
    cloning enzyme you want unique=True: a single cutter in the WHOLE plasmid.
    An insert that is internally clean still breaks the strategy if the
    backbone carries a second copy, which is why this takes the assembled
    construct rather than the insert.
    """
    enzymes = enzymes if enzymes is not None else ENZYME_SITES
    out = {}
    for name, site in enzymes.items():
        pos = find_sites(construct_seq, site, circular=circular)
        out[name] = {"positions": pos, "count": len(pos),
                     "unique": len(pos) == 1, "site": site}
    return out


def check_unique_and_clean(insert_dna, enzymes):
    """Insert-level screen: which candidate enzymes do not cut the insert.
    A pass here is necessary but not sufficient -- confirm with
    check_sites_in_construct on the assembled plasmid."""
    out = {}
    for name, site in enzymes.items():
        pos = find_sites(insert_dna, site, circular=False)
        out[name] = {"positions": pos, "count": len(pos), "clean": not pos,
                     "site": site}
    return out


def screen_methylation_context(seq, position, site_len=6, circular=True,
                               enzyme=None):
    """Dam/Dcm overlap around a site at `position`.

    Returns the overlapping/abutting GATC and CCWGG positions and the local
    context string. A hit is a candidate to confirm against REBASE, not a
    verdict: whether overlap actually blocks a given enzyme is
    enzyme-specific. `blocked_candidate` is True only when a hit coincides
    with an enzyme this module lists as sensitive.
    """
    seq = seq.upper()
    n = len(seq)
    lo, hi = position - 6, position + site_len + 6
    if circular:
        window = "".join(seq[i % n] for i in range(lo, hi))
    else:
        lo, hi = max(0, lo), min(n, hi)
        window = seq[lo:hi]

    dam, dcm = [], []
    for m in re.finditer("(?=GATC)", window):
        dam.append(lo + m.start())
    for m in re.finditer("(?=CC[AT]GG)", window):
        dcm.append(lo + m.start())

    site_lo, site_hi = position, position + site_len
    def overlaps(p, plen):
        return not (p + plen < site_lo - 1 or p > site_hi)

    dam_hits = [p % n for p in dam if overlaps(p, 4)]
    dcm_hits = [p % n for p in dcm if overlaps(p, 5)]

    blocked = None
    if enzyme is not None:
        blocked = bool((dam_hits and enzyme in DAM_SENSITIVE)
                       or (dcm_hits and enzyme in DCM_SENSITIVE))

    return {"dam": dam_hits, "dcm": dcm_hits, "context": window,
            "context_start": lo % n, "enzyme": enzyme,
            "blocked_candidate": blocked}


def synonymous_alternatives(dna_seq, codon_index, host="ecoli"):
    """Synonymous codons for the codon at codon_index, ranked host-preferred
    first, for removing a conflicting site without changing the protein."""
    codon = dna_seq[codon_index * 3:codon_index * 3 + 3].upper()
    aa = CODON_TABLE.get(codon)
    if aa is None:
        raise ValueError(f"codon {codon!r} at index {codon_index} is not valid")
    ranked = ECOLI_PREFERRED.get(aa, [])
    others = sorted(c for c, a in CODON_TABLE.items() if a == aa and c not in ranked)
    return {"codon": codon, "aa": aa,
            "alternatives": [c for c in ranked + others if c != codon]}


def remove_site_synonymously(insert_dna, site, frame_offset=0, host="ecoli"):
    """Propose silent codon changes that remove every copy of `site` from an
    insert, preferring host-optimal synonymous codons.

    `frame_offset` is the index within insert_dna where the reading frame
    starts. Returns one proposal per site occurrence; each gives the codon to
    change, the replacement and the resulting sequence. Returns the proposals
    rather than applying them -- a silent change needs explicit approval, and
    the full construct must be re-censused afterwards.
    """
    proposals = []
    for pos in find_sites(insert_dna, site, circular=False):
        for offset in range(len(site)):
            abs_pos = pos + offset
            if abs_pos < frame_offset:
                continue
            ci = (abs_pos - frame_offset) // 3
            try:
                alts = synonymous_alternatives(insert_dna[frame_offset:], ci, host)
            except ValueError:
                continue
            for alt in alts["alternatives"]:
                lo = frame_offset + ci * 3
                candidate = insert_dna[:lo] + alt + insert_dna[lo + 3:]
                if not find_sites(candidate, site, circular=False):
                    proposals.append({
                        "site_position": pos,
                        "codon_index": ci,
                        "codon_position": lo,
                        "from_codon": alts["codon"],
                        "to_codon": alt,
                        "aa": alts["aa"],
                        "host_rank": (ECOLI_PREFERRED.get(alts["aa"], []).index(alt)
                                      if alt in ECOLI_PREFERRED.get(alts["aa"], [])
                                      else None),
                        "sequence": candidate,
                    })
                    break
            if proposals and proposals[-1]["site_position"] == pos:
                break
    return proposals


# ------------------------------------------------------------- annotation

def make_product_html(description, accession=None, db="UNIPROT", url=None):
    """Build the accession-linked /product HTML SnapGene renders as a visible
    database link in the feature inspector. A plain-text /product does not
    surface the identifier the same way.
    """
    if accession is None:
        return f"<html><body>{escape(description)}</body></html>"
    db = db.upper()
    if url is None:
        template = DB_URL_TEMPLATES.get(db)
        if template is None:
            raise ValueError(f"no URL template for db={db!r}; pass url=")
        url = template.format(acc=accession)
    return (f"<html><body>{escape(description)}<br/>{escape(db)} "
            f'<u><a href="{escape(url)}">{escape(accession)}</a></u>'
            f"</body></html>")


def apply_feature_metadata(feature, construct_seq, description=None,
                           accession=None, db="UNIPROT", gene=None, note=None,
                           numbering_start=None, url=None):
    """Set a CDS feature's qualifiers and refresh its cached translation.

    The translation qualifier is the trap: SnapGene caches each CDS feature's
    peptide as a plain string, and an ORF replacement that rewrites the
    sequence and shifts the feature range does NOT refresh it. A construct can
    ship with an empty or wrong translation while its sequence, feature range,
    reading frame and enzyme positions are all correct. This recomputes it
    from construct_seq every time.

    codon_start and transl_table are written as ints, not strings: the two
    serialise to different XML value types.
    """
    feature.qualifiers["translation"] = feature_translation(construct_seq, feature)
    feature.qualifiers["codon_start"] = 1
    feature.qualifiers["transl_table"] = 1
    if gene is not None:
        feature.qualifiers["gene"] = gene
    if description is not None:
        feature.qualifiers["product"] = make_product_html(
            description, accession=accession, db=db, url=url)
    if note is not None:
        feature.qualifiers["note"] = f"<html><body>{escape(note)}</body></html>"
    if numbering_start is not None:
        feature.numbering_start = numbering_start
    for s in feature.segments:
        s.translated = True
    return feature


def set_construct_description(dna, description_html):
    """Replace the file-level info-pane description.

    Every vector carries a description left over from whatever it was last
    used for, and an ORF replacement copies it through verbatim -- it lives in
    the notes block, not in any feature. It is wrong for the new construct, so
    replace it every time rather than only when it looks stale; verify_construct
    can see that it is empty but cannot tell a stale one from a correct one.
    """
    dna.notes["Description"] = description_html
    return dna


# ---------------------------------------------------------- verification

def verify_construct(dna, insert_feature=None, expected_protein=None,
                     expect_unique=None, expect_absent=None):
    """Measure a construct and return what was measured plus a problems list.

    It asserts nothing on the caller's behalf. Report only what this returns:
    do not describe a construct as verified for a property the report does not
    cover, and say "not assessed" for a check that did not run. Sequence
    length and feature ranges being right tells you nothing about the
    qualifier-level defects, which is how an empty translation ships unnoticed.
    """
    problems, checks = [], {}
    checks["length"] = len(dna.seq)
    checks["topology"] = "circular" if dna.circular else "linear"
    checks["feature_count"] = len(dna.features)

    desc = dna.description
    checks["description_present"] = bool(desc and desc.strip())
    if not checks["description_present"]:
        problems.append("file-level description is empty")

    if insert_feature is not None:
        try:
            feat = dna.feature(insert_feature)
        except KeyError as exc:
            problems.append(str(exc))
            feat = None
        if feat is not None:
            span = feature_span(dna, feat)
            tr = translate_and_validate(span, expected_protein=expected_protein)
            if feat.spans_origin:
                problems.append(f"{feat.name}: spans the origin; check the "
                                f"junction by hand before trusting the frame")
            checks["insert"] = {
                "name": feat.name, "type": feat.type,
                "start": feat.start, "end": feat.end, "length_nt": len(span),
                "in_frame": tr["in_frame"],
                "internal_stops": tr["internal_stops"],
                "protein_length": tr["length_aa"],
                "protein_matches_expected": tr["match"],
                "first_mismatch": tr["first_mismatch"],
                "numbering_start": feat.numbering_start,
            }
            if not tr["in_frame"]:
                problems.append(f"{feat.name}: length {len(span)} is not a "
                                f"multiple of 3")
            if tr["internal_stops"]:
                problems.append(f"{feat.name}: internal stop codon(s) at "
                                f"residue {tr['internal_stops']}")
            if expected_protein is not None and not tr["match"]:
                problems.append(f"{feat.name}: translation does not match the "
                                f"expected protein ({tr['first_mismatch']})")

            cached = feat.qualifiers.get("translation")
            actual = feature_translation(dna, feat)
            checks["cached_translation_present"] = bool(cached)
            checks["cached_translation_matches"] = (cached == actual) if cached else False
            if not cached:
                problems.append(f"{feat.name}: cached translation qualifier is "
                                f"empty")
            elif cached != actual:
                problems.append(f"{feat.name}: cached translation qualifier is "
                                f"stale ({len(cached)} aa cached vs "
                                f"{len(actual)} aa actual)")

            product = feat.qualifiers.get("product", "")
            checks["product_has_link"] = "<a href=" in str(product)
            if not checks["product_has_link"]:
                problems.append(f"{feat.name}: /product carries no accession link")
            if not feat.qualifiers.get("note"):
                problems.append(f"{feat.name}: /note is empty")

    if expect_unique:
        census = check_sites_in_construct(
            dna.seq, {n: ENZYME_SITES[n] for n in expect_unique
                      if n in ENZYME_SITES}, circular=dna.circular)
        missing = [n for n in expect_unique if n not in ENZYME_SITES]
        if missing:
            problems.append(f"no recognition sequence on file for {missing}; "
                            f"pass them via check_sites_in_construct(enzymes=...)")
        checks["unique_sites"] = census
        for name, info in census.items():
            if not info["unique"]:
                problems.append(f"{name}: {info['count']} site(s) in the "
                                f"construct at {info['positions']}, wanted "
                                f"exactly 1")

    if expect_absent:
        census = check_sites_in_construct(
            dna.seq, {n: ENZYME_SITES[n] for n in expect_absent
                      if n in ENZYME_SITES}, circular=dna.circular)
        checks["absent_sites"] = census
        for name, info in census.items():
            if info["count"]:
                problems.append(f"{name}: expected absent, found "
                                f"{info['count']} at {info['positions']}")

    return {"all_ok": not problems, "problems": problems, "checks": checks}


def check_construct_name(name, max_len=50):
    """Measure a construct name against the length cap. The cap is easy to
    breach once tags and a range suffix are appended, so check it rather than
    eyeballing it."""
    return {"name": name, "length": len(name), "max": max_len,
            "ok": len(name) <= max_len, "over_by": max(0, len(name) - max_len)}


# ------------------------------------------------------- order sheet (xlsx)

def read_order_sheet_columns(path, sheet_index=0):
    """Read the header row of an .xlsx template, so the human's own template
    drives the column list instead of a hard-coded fallback.

    Returns the header texts in order. Handles both shared-string and inline
    cells. Formatting is not read: filling a template in place is not possible
    from the standard library, so write_order_sheet() reproduces the headers
    rather than editing the file.
    """
    with zipfile.ZipFile(path) as zf:
        shared = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
            for si in root.findall(f"{ns}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{ns}t")))
        sheets = sorted(n for n in zf.namelist()
                        if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
        if not sheets:
            raise ValueError(f"{path}: no worksheet found")
        root = ET.fromstring(zf.read(sheets[sheet_index]))

    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rows = root.iter(f"{ns}row")
    first = next(rows, None)
    if first is None:
        return []
    out = []
    for c in first.findall(f"{ns}c"):
        ctype = c.get("t")
        v = c.find(f"{ns}v")
        if ctype == "s" and v is not None:
            out.append(shared[int(v.text)])
        elif ctype == "inlineStr":
            is_el = c.find(f"{ns}is")
            out.append("".join(t.text or "" for t in is_el.iter(f"{ns}t"))
                       if is_el is not None else "")
        elif v is not None:
            out.append(v.text or "")
        else:
            out.append("")
    return out


def _col_letter(index):
    """1-based column index to spreadsheet letters."""
    letters = ""
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


_XLSX_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>"""

_XLSX_ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""

_XLSX_WORKBOOK_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>"""

# cellXfs: 0 = default, 1 = bold on grey, 2 = bold red on grey, 3 = wrapped body
_XLSX_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="3">
<font><sz val="11"/><name val="Calibri"/></font>
<font><b/><sz val="11"/><name val="Calibri"/></font>
<font><b/><color rgb="FFFF0000"/><sz val="11"/><name val="Calibri"/></font>
</fonts>
<fills count="3">
<fill><patternFill patternType="none"/></fill>
<fill><patternFill patternType="gray125"/></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFD9D9D9"/><bgColor indexed="64"/></patternFill></fill>
</fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="4">
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment vertical="center" wrapText="1"/></xf>
<xf numFmtId="0" fontId="2" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment vertical="center" wrapText="1"/></xf>
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top"/></xf>
</cellXfs>
</styleSheet>"""

_XLSX_WORKBOOK = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets>
</workbook>"""


def write_order_sheet(rows, out_path, columns=None, red_font_columns=None,
                      widths=None):
    """Write a CRO synthesis order sheet as .xlsx.

    `rows` is a list of dicts keyed by column name; a missing key becomes a
    blank cell and a key not in `columns` is reported, not silently dropped.
    Pass `columns` from read_order_sheet_columns() on the human's own template
    so their headers and order drive the sheet.

    Raises if any cell exceeds the 32767-character spreadsheet limit rather
    than writing a file the spreadsheet application will refuse -- a long
    vector sequence is the usual cause.
    """
    columns = list(columns) if columns else list(ORDER_SHEET_COLUMNS)
    red = set(red_font_columns if red_font_columns is not None
              else ORDER_SHEET_RED_COLUMNS)

    unknown = sorted({k for row in rows for k in row} - set(columns))
    too_long = []
    for i, row in enumerate(rows, start=2):
        for col in columns:
            v = row.get(col)
            if v is not None and len(str(v)) > XLSX_CELL_LIMIT:
                too_long.append((i, col, len(str(v))))
    if too_long:
        detail = ", ".join(f"row {r} {c!r} is {n} chars" for r, c, n in too_long)
        raise ValueError(f"cell(s) over the {XLSX_CELL_LIMIT}-character "
                         f"spreadsheet limit: {detail}")

    default_widths = {
        "No": 5, "Final plasmid Name": 42, "Final plasmid sequence": 30,
        "Insert sequence": 30, "Insert Source": 18,
        "Cloning sites (5' and 3')": 24, "5'homologous sequence (~30 bp)": 30,
        "3' hom (~30 bp)": 30, "Vector": 42, "Vector source": 20,
        "Vector sequence": 20, "Shippable": 11, "Maxiprep": 11,
    }
    default_widths.update(widths or {})

    # Element order is fixed by the OOXML schema: sheetViews, then cols,
    # then sheetData. A file with them in another order opens as repaired
    # or not at all.
    parts = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
             '<worksheet xmlns="http://schemas.openxmlformats.org/'
             'spreadsheetml/2006/main">',
             f'<dimension ref="A1:{_col_letter(len(columns))}'
             f'{len(rows) + 1}"/>',
             '<sheetViews><sheetView tabSelected="1" workbookViewId="0">'
             '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" '
             'state="frozen"/></sheetView></sheetViews>',
             '<sheetFormatPr defaultRowHeight="15"/>',
             "<cols>"]
    for i, col in enumerate(columns, start=1):
        parts.append(f'<col min="{i}" max="{i}" '
                     f'width="{default_widths.get(col, 20)}" customWidth="1"/>')
    parts.append("</cols><sheetData>")

    parts.append(f'<row r="1" ht="30" customHeight="1">')
    for i, col in enumerate(columns, start=1):
        style = 2 if col in red else 1
        parts.append(f'<c r="{_col_letter(i)}1" s="{style}" t="inlineStr">'
                     f"<is><t>{escape(str(col))}</t></is></c>")
    parts.append("</row>")

    for r, row in enumerate(rows, start=2):
        parts.append(f'<row r="{r}">')
        for i, col in enumerate(columns, start=1):
            value = row.get(col)
            if value is None or value == "":
                continue
            ref = f"{_col_letter(i)}{r}"
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                parts.append(f'<c r="{ref}" s="3"><v>{value}</v></c>')
            else:
                text = escape(str(value))
                space = ' xml:space="preserve"' if str(value) != str(value).strip() else ""
                parts.append(f'<c r="{ref}" s="3" t="inlineStr"><is>'
                             f"<t{space}>{text}</t></is></c>")
        parts.append("</row>")

    parts.append("</sheetData></worksheet>")
    sheet = "".join(parts)

    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", _XLSX_CONTENT_TYPES)
        zf.writestr("_rels/.rels", _XLSX_ROOT_RELS)
        zf.writestr("xl/workbook.xml", _XLSX_WORKBOOK)
        zf.writestr("xl/_rels/workbook.xml.rels", _XLSX_WORKBOOK_RELS)
        zf.writestr("xl/styles.xml", _XLSX_STYLES)
        zf.writestr("xl/worksheets/sheet1.xml", sheet)

    return {"path": out_path, "rows": len(rows), "columns": columns,
            "unknown_keys": unknown}
