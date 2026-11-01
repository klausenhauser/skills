#!/usr/bin/env python3
"""
check_embedded_restraints.py  <model.cif>  [more.cif ...]

Warns if a Phenix (or Coot) output mmCIF still carries embedded restraint
dictionaries (comp_*/link_*/mod_* blocks) before it gets reused as the input
model for another refinement round.

Phenix writes every comp_* and link_* dictionary it used into its own output
model. Feed that file back in as the next round's input and those embedded
definitions are silently re-read, overriding any link CIF you supply
separately on the command line -- with no error. Symptom: you edit a link
CIF (tighten a sigma, change a target) and the .geo shows the OLD values
regardless, because the model's own embedded copy won. See
strip_phenix_restraints.py for the fix.

Run this on any Phenix/Coot output before pointing phenix.refine at it as
the input model:

    python3 check_embedded_restraints.py refine_031.cif

Exit status: 0 if clean (safe to use as-is), 1 if embedded dictionary blocks
were found (usable as a gate in a shell pipeline: `... && phenix.refine ...`).

Confirmed incident this caught: one round's output carried an embedded
data_link_* block at the original loose sigma (0.050 A). It was fed straight
into the next round without stripping, so a freshly tightened external link
CIF (sigma 0.020 A) was silently ignored -- that round's own .geo still showed
sigma 0.050, and the linked bond barely moved.
"""
import sys
from pathlib import Path
import gemmi


def check_one(path):
    path = Path(path)
    doc = gemmi.cif.read_file(str(path))
    embedded = [b.name for b in doc if not b.find_loop("_atom_site.id")]
    if not embedded:
        print(f"{path.name}: clean, no embedded comp_/link_/mod_ blocks. Safe to use as input.")
        return 0
    print(f"{path.name}: EMBEDDED RESTRAINT BLOCKS FOUND -- {len(embedded)} block(s):")
    for name in embedded:
        print(f"  data_{name}")
    print(
        "  Phenix wrote these into its own output. Using this file as the input\n"
        "  model for another run will silently re-supply these definitions and\n"
        "  override any link CIF you pass separately -- e.g. a retightened sigma\n"
        "  in your external link file will have NO effect, with no error.\n"
        f"  Fix: python3 strip_phenix_restraints.py {path} <clean_out.cif>"
    )
    return 1


def main(paths):
    status = 0
    for p in paths:
        status |= check_one(p)
    return status


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1:]))
