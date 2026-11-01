# build-unnatural-peptides — scripts

Companion to the `build-unnatural-peptides` skill, which explains when and why to
use each of these. They ship inside the skill so they are available wherever the
skill is installed.

    pip install gemmi rdkit --break-system-packages

Phenix's own Python bundles both, so `phenix.python <script>` works too.

| Script | Use | Reusability |
|---|---|---|
| `split_ligand.py` | Split an already-fitted single-residue ligand into per-residue components with correct atom names and junction geometry. | **reference implementation — adapt the CONFIGURE block** |
| `fix_grade2_cif.py` | Post-process Grade2 output: rename atoms across every loop, set `_chem_comp.group` / `.type` / `.name`. | drop-in, config-driven |
| `fix_coot_cif.py` | Repair a Coot-saved mmCIF whose new residues have invalid `label_*` fields. Run on every Coot save. | drop-in |
| `check_embedded_restraints.py` | Exit 1 if a model still carries embedded `comp_*`/`link_*` blocks, before you point phenix.refine at it as input. | drop-in |
| `strip_phenix_restraints.py` | Remove those blocks. The fix when the check fires. | drop-in |
| `make_deposition_cif.py` | Add `struct_conn` rows with distances measured from the coordinates, preserving rows already present. | edit the `LINKS` table at the top |

## Typical order

    python3 fix_grade2_cif.py config.json          # all Grade2 output, full ligand included
    # ... user builds the FULL ligand in Coot, refines in Phenix ...
    python3 fix_coot_cif.py coot-0.cif             # every Coot save, always
    python3 split_ligand.py                        # after adapting it
    # ... link CIF + params, then refine ...
    python3 check_embedded_restraints.py refine_NN.cif   # before reusing as input
    python3 strip_phenix_restraints.py refine_NN.cif     # if the check fires
    python3 make_deposition_cif.py refine_final.cif

## The two traps that cost the most time

**1. A Phenix output carries its own dictionaries.** Every `comp_*` and `link_*`
block Phenix used is written into its output mmCIF. Reuse that file as the next
round's input and those definitions are re-read and silently override the link
CIF you supply — no error. Symptom: you tighten a sigma and the `.geo` keeps
showing the old value. This has happened twice, on two different projects, and
cost a refinement round each time. Run `check_embedded_restraints.py` before
every reuse, not only when something looks wrong.

**2. `.geo`'s printed bond values are stale.** They are computed once, right
after `pdb_interpretation` builds the model from the *input* coordinates, before
the macrocycles run. Use `.geo` for topology only — how many restraints at a
junction, which mods and `Link IDs` applied. For any restraint that started far
from its target, measure the final geometry from the output coordinate file
(gemmi `read_structure`, `Position.dist`). A bond has been seen reading an
identical `.geo` value across three consecutive rounds while really closing from
2.9 to 2.0 Å.

## split_ligand.py is a worked example, not a finished tool

It was written for one ligand and carries that ligand's assumptions: hardcoded
input paths, a `COMPONENTS` list naming specific components and their leaving
atoms, a `NUM` map of component to residue number, `manual()` builds for
standard CCD components with no Grade2 output to read, and a fixed output chain.

**Change all of that.** What is general, and worth keeping unchanged, is the
matching logic in the middle:

1. `comp_from_cif()` builds each component's heavy-atom graph from its restraint
   CIF, dropping the leaving atoms so the graph matches the residue as it exists
   *inside* the ligand.
2. Connectivity on the fitted ligand is perceived from coordinates with
   `rdDetermineBonds.DetermineConnectivity` — no bond orders, which avoids every
   aromatic-perception problem.
3. `generic()` strips bond orders from both query and target so matching is pure
   element-plus-topology.
4. The `solve()` backtracker finds placements that are mutually disjoint and
   together cover every heavy atom.

That last step is the safeguard. If an exact cover exists and each component
matched uniquely, the assignment is correct — not a guess. If no cover is found,
something is wrong with the component list or the leaving atoms, and the script
says so rather than producing a plausible-looking wrong answer.

Always verify the junction distances afterwards. They should be real bond
lengths, because they are the bonds from the fitted model.

## Where output goes

Alongside the Phenix `Refine_NN/` folders, never inside them:

    <project>--phenix/
      Refine_01/ ... Refine_NN/     Phenix's own; read only
      Split-Ligand/                 ligand_split.cif, adapted split script, logs
      Restraints/                   component CIFs, the combined link CIF, both params files
      Deposition/                   the finished CIF and the deposition notes

An input that already lives in a subfolder gets its output in that subfolder —
a Grade2 run folder keeps its own fixed CIF beside the raw output, and the copy
Phenix is pointed at goes in `Restraints/`. Sharing components between projects
means copying the component CIFs; splits and link restraints are per-structure
and stay separate.

Projects started under the older single-folder convention (`components/` or
`Claude-Unnaturals` holding everything) keep it — renaming would break paths
already referenced by adapted scripts and running Phenix jobs.
