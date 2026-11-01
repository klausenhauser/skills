#!/usr/bin/env python3
"""
fix_grade2_cif.py -- post-process Grade2 restraint CIFs for use as peptide
residues in Phenix / CCP4, and for PDB deposition.

Does three things, driven by a JSON config so it is reusable across projects:

  1. renames atom IDs consistently across EVERY loop that references atoms
     (_chem_comp_atom, _chem_comp_bond, _chem_comp_angle, _chem_comp_tor,
      _chem_comp_chir, _chem_comp_plane_atom, ...)
  2. sets _chem_comp.group and _chem_comp.type so the component is recognised
     as a peptide-linking residue (or left as a non-polymer cap)
  3. sets a human-readable _chem_comp.name (Grade2 writes the InChIKey there)

Usage:
    python3 fix_grade2_cif.py config.json

Config format:
{
  "outdir": ".",
  "components": {
    "RES1": {
      "input":  "../RES1/RES1.restraints.cif",
      "output": "RES1.restraints.fixed.cif",
      "rename": {"C5": "C", "O2": "O", "O1": "OXT", "H1": "HXT"},
      "group":  "non-polymer",
      "type":   "NON-POLYMER",
      "name":   "example N-terminal cap carboxylic acid"
    }
  }
}

Any of rename / group / type / name may be omitted.
"""
import json
import sys
from pathlib import Path

import gemmi


def _columns(loop):
    """Return loop contents as a list of columns (gemmi .values is a copy)."""
    w, n = loop.width(), loop.length()
    vals = loop.values
    return [[vals[r * w + c] for r in range(n)] for c in range(w)]


def rename_atoms(block, mapping):
    """Rename atom IDs in every tag containing 'atom_id' in this block."""
    if not mapping:
        return 0
    changed = 0
    for item in block:
        loop = item.loop
        if loop is None:
            continue
        cols = [i for i, t in enumerate(loop.tags) if "atom_id" in t]
        if not cols:
            continue
        columns = _columns(loop)
        for c in cols:
            for r, v in enumerate(columns[c]):
                raw = gemmi.cif.as_string(v)
                if raw in mapping:
                    columns[c][r] = gemmi.cif.quote(mapping[raw])
                    changed += 1
        loop.set_all_values(columns)
    return changed


def set_comp_fields(doc, comp_id, fields):
    """Update _chem_comp.{group,type,name} for comp_id in data_comp_list."""
    done = {}
    for block in doc:
        for item in block:
            loop = item.loop
            if loop is None or "_chem_comp.id" not in loop.tags:
                continue
            tags = list(loop.tags)
            columns = _columns(loop)
            i_id = tags.index("_chem_comp.id")
            touched = False
            for row, v in enumerate(columns[i_id]):
                if gemmi.cif.as_string(v) != comp_id:
                    continue
                for key, value in fields.items():
                    tag = "_chem_comp." + key
                    if tag not in tags:
                        continue
                    columns[tags.index(tag)][row] = gemmi.cif.quote(value)
                    done[key] = value
                    touched = True
            if touched:
                loop.set_all_values(columns)
    return done


def verify(path):
    """Sanity check: every atom referenced in a loop exists in _chem_comp_atom."""
    doc = gemmi.cif.read_file(str(path))
    problems = []
    for block in doc:
        names = set()
        for item in block:
            loop = item.loop
            if loop is None or "_chem_comp_atom.atom_id" not in loop.tags:
                continue
            c = loop.tags.index("_chem_comp_atom.atom_id")
            names = {gemmi.cif.as_string(v) for v in _columns(loop)[c]}
        if not names:
            continue
        for item in block:
            loop = item.loop
            if loop is None:
                continue
            columns = _columns(loop)
            for i, t in enumerate(loop.tags):
                if "atom_id" not in t or t.startswith("_chem_comp_atom."):
                    continue
                for v in columns[i]:
                    raw = gemmi.cif.as_string(v)
                    if raw and raw not in names:
                        problems.append(f"{t}: '{raw}' not in atom list")
    return sorted(set(problems))


def main(config_path):
    cfg = json.loads(Path(config_path).read_text())
    root = Path(config_path).parent
    outdir = (root / cfg.get("outdir", ".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    for comp_id, spec in cfg["components"].items():
        src = (root / spec["input"]).resolve()
        dst = outdir / spec.get("output", f"{comp_id}.restraints.fixed.cif")
        doc = gemmi.cif.read_file(str(src))

        n = 0
        for block in doc:
            n += rename_atoms(block, spec.get("rename", {}))

        fields = {k: spec[k] for k in ("group", "type", "name") if k in spec}
        applied = set_comp_fields(doc, comp_id, fields)

        doc.write_file(str(dst))

        print(f"{comp_id}:")
        print(f"  {src.name} -> {dst.name}")
        if spec.get("rename"):
            pairs = ", ".join(f"{a}->{b}" for a, b in spec["rename"].items())
            print(f"  renamed {n} references ({pairs})")
        for k, v in applied.items():
            print(f"  _chem_comp.{k} = {v}")
        bad = verify(dst)
        print("  verify: OK" if not bad else "  verify: FAILED")
        for p in bad:
            print(f"    ! {p}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
