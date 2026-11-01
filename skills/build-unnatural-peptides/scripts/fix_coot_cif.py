#!/usr/bin/env python3
"""
fix_coot_cif.py -- repair mmCIF files written by Coot so Phenix can read them.

THE PROBLEM
    Coot writes residues added or rebuilt in the current session with invalid
    mmCIF label_* fields, while the auth_* fields stay correct:

        _atom_site.label_atom_id   .          (a nameless atom)
        _atom_site.label_comp_id   ---
        _atom_site.label_asym_id   .
        _atom_site.label_seq_id    -32767

    Phenix reads label_*, not auth_*, so it cannot type those atoms:

        Sorry: Fatal problems interpreting model file:
          Unknown nonbonded energy type

    A nameless atom (label_atom_id = ".") is the usual trigger, and it is
    easy to miss because the file looks fine in Coot and in ChimeraX, both
    of which read auth_*.

THE FIX
    A gemmi round trip. gemmi reads the model preferring auth_*, so it
    recovers the real chemistry; clearing label_seq, then setup_entities()
    and assign_label_seq_id(force=True), rebuilds the entity and subchain
    bookkeeping and regenerates every label_* field from scratch.
    Coordinates never move.

    label_seq_id must be cleared explicitly -- gemmi otherwise carries
    Coot's -32767 sentinel straight through to the output.

    The round trip also drops the comp_*/link_*/mod_* restraint blocks that
    Phenix embeds in its output models and Coot passes along. Those silently
    override link dictionaries given on the command line, which is its own
    hard-to-find bug, so removing them is deliberate. --keep-restraints
    keeps them.

WHAT THIS DOES NOT FIX
    "Unknown nonbonded energy type" has a second, unrelated cause: a ligand
    or non-standard residue whose restraint CIF was not supplied to Phenix
    at all. If this script reports a clean file and Phenix still fails,
    check that every component listed at the end of the run is covered by a
    restraint file.

ALSO FIXED: undefined ANISOU records
    Some Coot sessions leave a *phantom* anisotropic-ADP record on a
    rebuilt residue: the atom is listed in the ANISOU/atom_site_anisotrop
    loop (so Phenix treats it as anisotropic), but all six u[i][j] values
    are "?" -- never actually computed. gemmi parses "?" as NaN, and the
    first structure-factor calculation in phenix.refine hits a
    non-positive-definite U tensor and dies immediately, before any
    refinement runs, with:

        RuntimeError: cctbx Internal Error: .../sampling_base.h(159):
        CCTBX_ASSERT(d > 0) failure.

    This has nothing to do with label_* fields and isn't caught by the fix
    described above -- a file with clean labels can still carry this and
    still crash phenix.refine. This script now also detects atoms with a
    present-but-non-finite ANISOU record and clears them back to isotropic
    (matching how Coot correctly left its *other* rebuilt residues in the
    same session, which have no ANISOU record at all). The existing B_iso
    is kept; coordinates are never touched.

USAGE
    python3 fix_coot_cif.py coot-0.cif              -> coot-0_fixed.cif
    python3 fix_coot_cif.py coot-0.cif -o ready.cif
    python3 fix_coot_cif.py *.cif                   # batch
    python3 fix_coot_cif.py coot-0.cif --in-place    # overwrites, keeps .bak
    python3 fix_coot_cif.py coot-0.cif --keep-restraints

    If gemmi is not in your default python, use one you already have:
        phenix.python fix_coot_cif.py coot-0.cif
        ccp4-python  fix_coot_cif.py coot-0.cif

VERIFICATION
    Reports atom count before and after, the maximum coordinate deviation,
    and re-checks the output for invalid label_* fields. Max deviation must
    read 0.000000. Exit status is nonzero if any check fails, so this is
    safe to chain ahead of a refinement in a shell script.

REQUIRES
    gemmi   (pip install gemmi, or the copy bundled with Phenix / CCP4)
"""

import argparse
import math
import shutil
import sys
from pathlib import Path

try:
    import gemmi
except ImportError:
    sys.exit("error: gemmi not found.  Try:  pip install gemmi\n"
             "       or run this with  phenix.python  /  ccp4-python")

# label_comp_id and label_atom_id must name something; "." is not acceptable
# for either. For label_seq_id, "." IS correct for a non-polymer residue, so
# only Coot's sentinel and other non-positive values are wrong.
BAD_NAME = {"---", "--", "?", ".", ""}


def model_block(doc):
    """The data block holding the coordinates.

    Not sole_block(): a real Coot or Phenix mmCIF usually carries comp_*
    and link_* restraint blocks next to the model, and sole_block() raises
    on those.
    """
    for block in doc:
        if block.find_values("_atom_site.Cartn_x"):
            return block
    return doc[0] if len(doc) else None


def _col(block, tag):
    """Values of _atom_site.<tag>, or None if that tag is absent.

    Deliberately one tag at a time. Asking cif.find() for several tags at
    once returns an EMPTY table if any single one is missing, which reads
    as "nothing wrong here" and hides exactly the problem this script
    exists to find.
    """
    values = block.find_values("_atom_site." + tag)
    if not values:
        return None
    return [v.strip().strip("'\"") for v in values]


def _aniso_components(atom):
    u = atom.aniso
    return (u.u11, u.u22, u.u33, u.u12, u.u13, u.u23)


def find_bad_aniso(path):
    """Atoms carrying an ANISOU record whose u[i][j] values are undefined
    (Coot writes "?" for all six; gemmi parses that as NaN).

    Returns a list of (chain, comp, seq, atom_name) using auth_* identifiers
    -- the same identifiers diagnose() reports for label_* problems, so the
    two checks can share the same "affected residues" listing.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    bad = []
    for model in st:
        for chain in model:
            for res in chain:
                for a in res:
                    if a.aniso.nonzero() and not all(
                            math.isfinite(v) for v in _aniso_components(a)):
                        bad.append((chain.name, res.name, res.seqid.num, a.name))
        break  # first model only; mmCIF from Coot/Phenix has one model
    return bad


def diagnose(path):
    """Return (problems, n_rows, per-row descriptions, affected residues)."""
    block = model_block(gemmi.cif.read(str(path)))
    if block is None:
        return ["file contains no data block"], 0, [], []

    comp = _col(block, "label_comp_id")
    atom = _col(block, "label_atom_id")
    asym = _col(block, "label_asym_id")
    seq = _col(block, "label_seq_id")
    a_comp = _col(block, "auth_comp_id")
    a_seq = _col(block, "auth_seq_id")
    a_asym = _col(block, "auth_asym_id")

    n = max((len(c) for c in (comp, atom, asym, seq) if c), default=0)
    if n == 0:
        return ["no _atom_site loop found"], 0, [], []

    problems = []
    for tag, col in (("label_comp_id", comp), ("label_atom_id", atom)):
        if col is None:
            problems.append(f"_atom_site.{tag} is missing entirely")

    counts = {"label_atom_id": 0, "label_comp_id": 0,
              "label_asym_id": 0, "label_seq_id": 0}
    residues = {}

    for i in range(n):
        issues = []
        if atom and atom[i] in BAD_NAME:
            counts["label_atom_id"] += 1
            issues.append("nameless atom")
        if comp and comp[i] in BAD_NAME:
            counts["label_comp_id"] += 1
            issues.append(f"comp_id={comp[i] or 'empty'}")
        if asym and asym[i] in BAD_NAME:
            counts["label_asym_id"] += 1
            issues.append("no asym_id")
        if seq:
            v = seq[i]
            if v not in (".", "?"):          # "." is legal for non-polymers
                try:
                    if int(v) < 1:
                        counts["label_seq_id"] += 1
                        issues.append(f"seq_id={v}")
                except ValueError:
                    counts["label_seq_id"] += 1
                    issues.append(f"seq_id={v}")
        if issues:
            key = (a_comp[i] if a_comp else "?",
                   a_asym[i] if a_asym else "?",
                   a_seq[i] if a_seq else "?")
            residues.setdefault(key, set()).update(issues)

    for tag, c in counts.items():
        if c:
            problems.append(f"{c} row(s) with invalid {tag}")

    bad_aniso = find_bad_aniso(path)
    if bad_aniso:
        problems.append(f"{len(bad_aniso)} atom(s) with undefined ANISOU (u=?)")
        for ch, comp, seq, atom_name in bad_aniso:
            key = (comp, ch, str(seq))
            residues.setdefault(key, set()).add("undefined ANISOU")

    rows = [f"{c} {ch} {s}: {', '.join(sorted(iss))}"
            for (c, ch, s), iss in residues.items()]
    return problems, n, rows, list(residues)


def strip_restraints(doc):
    removed = []
    for i in range(len(doc) - 1, -1, -1):
        if doc[i].name.startswith(("comp_", "link_", "mod_")):
            removed.append(doc[i].name)
            del doc[i]
    return removed


def atom_snapshot(path):
    """(n_atoms, [(x, y, z), ...]) in file order."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    pos = [(a.pos.x, a.pos.y, a.pos.z)
           for model in st for chain in model for res in chain for a in res]
    return len(pos), pos


def fix(path, out_path, keep_restraints=False, in_place=False):
    path = Path(path)
    print(f"\n=== {path.name} ===")

    problems, n_rows, rows, _ = diagnose(path)
    if problems:
        for p in problems:
            print(f"  found: {p}")
        print(f"  affected residues ({len(rows)}):")
        for r in rows[:15]:
            print(f"    {r}")
        if len(rows) > 15:
            print(f"    ... and {len(rows) - 15} more")
    else:
        print(f"  {n_rows} atom rows, no invalid label_* fields "
              "(round trip will still normalise the file)")

    n_before, pos_before = atom_snapshot(path)
    if n_before == 0:
        print("  ** gemmi read 0 atoms -- this is not a model gemmi can "
              "parse (truncated file, or an incomplete _atom_site loop). "
              "Nothing written. **")
        return False

    # ---- the repair ------------------------------------------------------
    st = gemmi.read_structure(str(path))

    # Clear phantom ANISOU records: present but with undefined (u=?, NaN
    # once parsed) values. cctbx trips CCTBX_ASSERT(d > 0) in
    # sampling_base.h on these the moment it computes structure factors --
    # before refinement even starts -- so this has to be fixed regardless
    # of whether label_* is otherwise clean. Clearing it drops the atom
    # back to isotropic, which is how Coot correctly left its *other*
    # rebuilt residues in the same session. B_iso and coordinates untouched.
    cleared_aniso = 0
    for model in st:
        for chain in model:
            for res in chain:
                for a in res:
                    if a.aniso.nonzero() and not all(
                            math.isfinite(v) for v in _aniso_components(a)):
                        a.aniso = gemmi.SMat33f(0, 0, 0, 0, 0, 0)
                        cleared_aniso += 1
    if cleared_aniso:
        print(f"  cleared {cleared_aniso} undefined ANISOU record(s) "
              "(u=? -> isotropic)")

    # Clear Coot's -32767 sentinel first, or gemmi writes it back out.
    for model in st:
        for chain in model:
            for res in chain:
                if res.label_seq is not None and res.label_seq < 1:
                    res.label_seq = None

    st.setup_entities()                  # rebuild entities and subchains
    st.assign_label_seq_id(True)         # regenerate label_seq_id

    # assign_label_seq_id() will not renumber a polymer for which no reference
    # sequence exists -- it keeps whatever the file had. Coot's own numbering
    # already SKIPS the residues it added, so those residues land in the
    # polymer subchain with label_seq_id = "." , which is inconsistent mmCIF.
    # Renumber 1..N over each polymer, but only for entities that declare no
    # reference sequence: where _entity_poly_seq exists, label_seq_id refers
    # to the full construct and must not be touched.
    no_reference = {e.name for e in st.entities if len(e.full_sequence) == 0}
    renumbered = 0
    for model in st:
        for chain in model:
            poly = chain.get_polymer()
            if len(poly) == 0:
                continue
            ent = st.get_entity_of(poly)
            if ent is not None and ent.name not in no_reference:
                continue
            for i, res in enumerate(poly, start=1):
                if res.label_seq != i:
                    res.label_seq = i
                    renumbered += 1
    if renumbered:
        print(f"  renumbered label_seq_id for {renumbered} polymer residue(s) "
              "(no reference sequence in file)")

    doc = st.make_mmcif_document()       # regenerate every label_* field

    if not keep_restraints:
        removed = strip_restraints(doc)
        if removed:
            print(f"  removed {len(removed)} restraint block(s): "
                  + ", ".join(removed[:5])
                  + (f" ... +{len(removed) - 5}" if len(removed) > 5 else ""))

    if in_place:
        backup = path.with_suffix(path.suffix + ".bak")
        shutil.copy2(path, backup)
        print(f"  backup: {backup.name}")
        out_path = path

    doc.write_file(str(out_path))

    # ---- verify ----------------------------------------------------------
    n_after, pos_after = atom_snapshot(out_path)
    if n_after == 0:
        print("  ** output has 0 atoms -- do not use this file **")
        return False
    if n_before != n_after:
        print(f"  ** atom count changed {n_before} -> {n_after} -- "
              "do not use this file **")
        return False

    max_dev = max((max(abs(a - b) for a, b in zip(p, q))
                   for p, q in zip(pos_before, pos_after)), default=0.0)
    print(f"  atoms: {n_before} -> {n_after}   "
          f"max coordinate deviation: {max_dev:.6f} A")
    if max_dev > 1e-4:
        print("  ** coordinates moved -- inspect before refining **")
        return False

    still, _, still_rows, _ = diagnose(out_path)
    if still:
        print(f"  ** output still has problems: {'; '.join(still)} **")
        for r in still_rows[:5]:
            print(f"    {r}")
        return False

    print(f"  wrote {Path(out_path).name}   label_* regenerated, "
          "coordinates identical")
    return True


def main():
    ap = argparse.ArgumentParser(
        description="Repair Coot-written mmCIF files so Phenix can read them.")
    ap.add_argument("inputs", nargs="+", help="mmCIF file(s) from Coot")
    ap.add_argument("-o", "--output", help="output path (single input only)")
    ap.add_argument("--in-place", action="store_true",
                    help="overwrite the input, keeping a .bak copy")
    ap.add_argument("--keep-restraints", action="store_true",
                    help="keep embedded comp_*/link_* blocks "
                         "(default is to strip them)")
    args = ap.parse_args()

    if args.output and len(args.inputs) > 1:
        sys.exit("error: -o takes a single input file")
    if args.output and args.in_place:
        sys.exit("error: -o and --in-place are mutually exclusive")

    ok = True
    for item in args.inputs:
        p = Path(item)
        if not p.is_file():
            print(f"skipping {item}: not a file")
            ok = False
            continue
        out = (Path(args.output) if args.output
               else p.with_name(p.stem + "_fixed" + p.suffix))
        ok &= fix(p, out, keep_restraints=args.keep_restraints,
                  in_place=args.in_place)

    print()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
