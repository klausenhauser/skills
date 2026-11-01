# build-unnatural-peptides

Takes a crystal structure with an unnatural peptide or a covalently bound
ligand from "fitted into density" to "refined with proper restraints and ready
for PDB deposition".

It is for ligands that Phenix and the PDB cannot treat as a single standard
residue:

- peptides with non-standard residues, caps, staples or macrocycles, which must
  be deposited as one chemical component per residue, linked together;
- ligands covalently bonded to a protein side chain through any warhead (SuFEx,
  sulfonyl fluoride, acrylamide, epoxide, boronate…).

Claude reads [`SKILL.md`](SKILL.md) for the full instructions. This page is the
human overview.

## How the work is shared

The work alternates between Claude and you. **Claude** prepares the restraint
dictionaries, the split ligand, link definitions, params files and deposition
files. **You** run Grade2, build in Coot and refine in Phenix. Claude stops at
each hand-over and waits for your build or refinement.

## The logic

1. **Fit the whole ligand first, then split it.** The ligand is fitted into the
   density as a single residue. Only then is it cut into per-residue components,
   from those fitted coordinates. Atoms never move, so every junction comes out
   at a real bond length. Building from separate monomers gives junctions about
   1 Å too long, which nothing will link.
2. **Prove the split, don't guess it.** `split_ligand.py` matches each
   component's atom graph onto the fitted ligand and requires an *exact cover*:
   every heavy atom claimed exactly once. A unique cover proves the assignment.
   If there is no cover, the component list is wrong and the script says so.
3. **Restrain links softly, then tighten.** A covalent link to the protein
   starts loose and is tightened only once the geometry is close. Otherwise one
   strained bond distorts the model.
4. **Measure from coordinates, not from Phenix's `.geo`.** The `.geo` file
   reports geometry from before refinement ran.
5. **Deposit with the links recorded.** OneDep splits a multi-residue ligand
   into separate entities. The `struct_conn` covalent-bond rows are the only
   record that they were one linked molecule.

## Steps

| Step | Who | What happens |
|---|---|---|
| 0. Overview | Claude asks, you confirm | Target, ligand, resolution, the residue that reacts, what already exists, the deliverable. **The protein chain and the ligand chain are defined first** and used in every later step. |
| 1. Full ligand | You, with Claude preparing | Grade2 dictionary for the whole ligand. Claude decides how to split it into components and gives you the SMILES of each. You run Grade2 on them and Claude fixes the output (`fix_grade2_cif.py`). You build the full ligand in Coot and refine in Phenix. |
| 2. Split | Claude, then you | Claude splits the fitted ligand (`split_ligand.py`), checks every junction, and writes the peptide link dictionary and params. You refine. |
| 2a. Covalent link | Claude, then you | Reactive ligands only, after the peptide links work. Claude writes the warhead link with literature geometry targets and stages its restraint from soft to tight. You refine, possibly over several rounds. |
| Each round | Claude | Checks for embedded dictionaries before reusing an output (`check_embedded_restraints.py`, `strip_phenix_restraints.py`), repairs Coot saves (`fix_coot_cif.py`), and keeps a running record of R-work, R-free, bond lengths and clashscore. |
| 3. Deposition | Claude | Cleans the final model for OneDep (`make_deposition_cif.py`), records the covalent bonds in `struct_conn`, verifies nothing moved, and writes deposition notes for the annotators. |

## Scripts

All in [`scripts/`](scripts/). See [`scripts/README.md`](scripts/README.md) for
usage. Claude copies a script into your project and adapts it there; the
originals here stay untouched.

| Script | Use | Edit before running? |
|---|---|---|
| `fix_grade2_cif.py` | Fix Grade2 output: atom names, residue type, chemical name | JSON config |
| `fix_coot_cif.py` | Repair Coot-saved mmCIF so Phenix can read it; run on every Coot save | No |
| `split_ligand.py` | Split the fitted ligand into components | Yes: the CONFIGURE block (model, chains, components, anchor) |
| `check_embedded_restraints.py` | Warn if a Phenix output still carries its own dictionaries | No |
| `strip_phenix_restraints.py` | Remove those dictionaries | No |
| `make_deposition_cif.py` | Clean the final model for deposition and add `struct_conn` rows | Yes: the `LINKS` table |

## Requirements

- **Phenix** and **Coot**, installed by you.
- **Grade2:** needs a Global Phasing licence. Without a local install, use the
  free [Grade Web Server](http://grade.globalphasing.org).
- **Python with gemmi and RDKit** to run the scripts. `phenix.python` includes
  both.

Works best in **Claude Code** on the machine where Phenix lives. On claude.ai,
Claude can prepare and check files you upload, but you run Grade2, Coot and
Phenix yourself.

## Project folder layout

Claude writes its files beside Phenix's job folders, never inside them:

```
<project>--phenix/
  Refine_01/ ... Refine_NN/   Phenix's own
  Split-Ligand/               ligand_split.cif, adapted split script, logs
  Restraints/                 component CIFs, combined link CIF, params files
  Deposition/                 final CIF and deposition notes
```
