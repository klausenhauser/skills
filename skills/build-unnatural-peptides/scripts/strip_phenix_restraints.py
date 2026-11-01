#!/usr/bin/env python3
"""
strip_phenix_restraints.py  <in.cif>  [out.cif]

Phenix writes its refined model as an mmCIF that ALSO contains every
restraint dictionary used in that run -- comp_* blocks and, critically,
link_* blocks. Feed that file back in as the model for the next round and
those definitions are re-read, silently overriding the link CIFs you
supply on the command line. Symptom: the .geo shows sigmas you did not ask
for, and no amount of swapping link files changes anything.

This keeps only the coordinate block and drops the rest.

    python3 strip_phenix_restraints.py refine_020.cif refine_020_clean.cif
"""
import sys
from pathlib import Path
import gemmi

def main(src, dst=None):
    src = Path(src)
    dst = Path(dst) if dst else src.with_name(src.stem + "_clean.cif")
    doc = gemmi.cif.read_file(str(src))
    keep, drop = [], []
    for block in doc:
        (keep if block.find_loop("_atom_site.id") else drop).append(block.name)
    if not keep:
        sys.exit(f"{src.name}: no block with _atom_site found")
    for i in range(len(doc) - 1, -1, -1):
        if doc[i].name in drop:
            del doc[i]
    doc.write_file(str(dst))
    print(f"{src.name} -> {dst.name}")
    print(f"  kept:    {', '.join(keep)}")
    print(f"  removed: {', '.join(drop) if drop else '(nothing)'}")

if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    main(*sys.argv[1:])
