#!/usr/bin/env python3
"""Check this module's reader and writer against real SnapGene files.

    python3 selftest.py /path/to/a/folder/of/dna/files

The folder is searched recursively for *.dna. Nothing is written into it; all
output goes to a temporary directory. No files ship with this test -- point it
at your own vectors, the more varied the better.

Three checks, and the third is the one that matters most:

  translations  Every cached `translation` qualifier is recomputed from the
                sequence and compared to the stored string. This pins the two
                conventions that are easy to get wrong: a terminal stop is kept
                as a trailing `*`, and a multi-segment feature is translated per
                segment and comma-joined.

  round-trip    Read, write, read again, and compare everything modelled.
                Catches loss of sequence, topology, features, segments,
                colours, directionality, qualifiers, numbering and notes.

  encoding      Compare the <V text=...> bytes this writer produces against the
                bytes SnapGene itself wrote for the same qualifier.

The third exists because the second cannot replace it. A round-trip proves the
writer agrees with the reader, not that either agrees with SnapGene: a value
escaped once too often on the way out and decoded once too often on the way in
round-trips perfectly while SnapGene renders the qualifier as literal
`<html><body>` tags. Only a comparison against SnapGene's own bytes sees that.
"""

import pathlib
import re
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import snapgene as sg  # noqa: E402

_QUALIFIER = re.compile(
    r'<Q name="([^"]*)"><V text=(".*?"|\'.*?\')/></Q>', re.S)
_FEATURE = re.compile(
    r'<Feature\b[^>]*name="([^"]*)"[^>]*>(.*?)</Feature>', re.S)


def _text_qualifier_bytes(path):
    """Map (feature name, qualifier) -> the raw attribute value as written,
    quotes and entity references intact, read straight out of the file."""
    for btype, payload in sg.read_blocks(str(path)):
        if btype != sg.BLOCK_FEATURES:
            continue
        xml = payload.decode("utf-8", "replace")
        found = {}
        for feat in _FEATURE.finditer(xml):
            for qual in _QUALIFIER.finditer(feat.group(2)):
                found[(feat.group(1), qual.group(1))] = qual.group(2)
        return found
    return {}


def check_translations(files):
    total = bad = 0
    problems = []
    for path in files:
        dna = sg.read_dna(str(path))
        for feature in dna.features:
            cached = feature.qualifiers.get("translation")
            if not cached:
                continue
            total += 1
            if sg.feature_translation(dna, feature) != cached:
                bad += 1
                problems.append(f"{path.name}: {feature.name}")
    return total, bad, problems


def check_round_trip(files, workdir):
    total = bad = 0
    problems = []
    for path in files:
        total += 1
        before = sg.read_dna(str(path))
        copy = workdir / ("rt_" + path.name)
        sg.write_dna(before, str(copy))
        after = sg.read_dna(str(copy))

        issues = []
        if before.seq != after.seq:
            issues.append("sequence")
        for attr in ("circular", "double_stranded", "dam", "dcm", "ecoki"):
            if getattr(before, attr) != getattr(after, attr):
                issues.append(attr)
        if before.notes != after.notes:
            issues.append("notes")
        if len(before.features) != len(after.features):
            issues.append(f"feature count {len(before.features)}"
                          f"->{len(after.features)}")
        else:
            for a, b in zip(before.features, after.features):
                if (a.name, a.type, str(a.directionality)) != \
                   (b.name, b.type, str(b.directionality)):
                    issues.append(f"identity of {a.name!r}")
                if [(s.start, s.end, s.translated, s.color, s.type)
                        for s in a.segments] != \
                   [(s.start, s.end, s.translated, s.color, s.type)
                        for s in b.segments]:
                    issues.append(f"segments of {a.name!r}")
                if a.qualifiers != b.qualifiers:
                    differing = sorted(
                        k for k in set(a.qualifiers) | set(b.qualifiers)
                        if a.qualifiers.get(k) != b.qualifiers.get(k))
                    issues.append(f"qualifiers {differing} of {a.name!r}")
                if a.numbering_start != b.numbering_start:
                    issues.append(f"numbering of {a.name!r}")
        if issues:
            bad += 1
            problems.append(f"{path.name}: " + "; ".join(issues))
    return total, bad, problems


def check_encoding(files, workdir):
    """Compare written qualifier bytes against SnapGene's own."""
    total = bad = html_valued = 0
    problems = []
    for path in files:
        theirs = _text_qualifier_bytes(path)
        if not theirs:
            continue
        copy = workdir / ("enc_" + path.name)
        sg.write_dna(sg.read_dna(str(path)), str(copy))
        ours = _text_qualifier_bytes(copy)
        for key, their_bytes in theirs.items():
            if key not in ours:
                continue
            total += 1
            if "&lt;html&gt;" in their_bytes:
                html_valued += 1
            if their_bytes != ours[key]:
                bad += 1
                if len(problems) < 5:
                    problems.append(
                        f"{path.name} {key[0]}/{key[1]}\n"
                        f"      snapgene: {their_bytes[:100]}\n"
                        f"      ours:     {ours[key][:100]}")
    return total, bad, html_valued, problems


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    root = pathlib.Path(argv[1]).expanduser()
    files = sorted(root.rglob("*.dna"))
    if not files:
        print(f"no .dna files under {root}")
        return 2

    failed = False
    with tempfile.TemporaryDirectory(prefix="snapgene-selftest-") as tmp:
        workdir = pathlib.Path(tmp)
        print(f"{len(files)} file(s) under {root}\n")

        total, bad, problems = check_translations(files)
        print(f"translations  {total - bad}/{total} reproduce exactly")
        for p in problems[:5]:
            print(f"    {p}")
        failed |= bool(bad)

        total, bad, problems = check_round_trip(files, workdir)
        print(f"round-trip    {total - bad}/{total} lossless")
        for p in problems[:5]:
            print(f"    {p}")
        failed |= bool(bad)

        total, bad, html_valued, problems = check_encoding(files, workdir)
        print(f"encoding      {total - bad}/{total} byte-identical to "
              f"SnapGene's own ({html_valued} carry HTML)")
        for p in problems:
            print(f"    {p}")
        failed |= bool(bad)

    print("\n" + ("FAILED" if failed else "all checks passed"))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
