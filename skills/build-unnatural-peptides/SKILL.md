---
name: "build-unnatural-peptides"
description: "Split an unnatural peptide or covalent ligand into per-residue PDB chemical components, generate and fix Grade2 restraints, link residues to each other or to a protein side chain for Phenix refinement, and assemble a deposition-ready mmCIF. Use when refining or depositing a crystal structure that contains a peptide with non-standard residues, a stapled or macrocyclic peptide, or a ligand covalently bonded to the protein, and Phenix or Coot need custom restraints or link definitions."
---

# Linking unnatural residues and covalent ligands for refinement and deposition

For peptidic ligands that must become per-residue chemical components, and for
ligands covalently bonded to a protein side chain (SuFEx, acrylamide, epoxide,
sulfonyl fluoride, boronate, any warhead). The method is the same whether the
link is residue-to-residue or residue-to-ligand.

**The work alternates.** Claude prepares components, restraints, link
dictionaries, splits and deposition files. The user builds in Coot and runs
Phenix. Each step ends where the user takes over — do not race ahead of a build
or a refinement that has not happened yet.

## Three rules that override everything below

### 1. Use the scripts bundled with this skill

The reference scripts ship inside this skill's own `scripts/` folder.
**Read its README first, every time. Use these scripts — do not write ad-hoc
replacements**, not even quick inline ones.

| Script | Use |
|---|---|
| `split_ligand.py` | Split the fitted full ligand into components. Reference implementation: adapt the CONFIGURE block, leave the matching logic alone. |
| `fix_grade2_cif.py` | Post-process Grade2 CIFs from a JSON config. |
| `fix_coot_cif.py` | Repair Coot-written mmCIF. Run on every Coot save. |
| `check_embedded_restraints.py` | Exit 1 if a model still carries embedded `comp_*`/`link_*` blocks. |
| `strip_phenix_restraints.py` | Remove those blocks. |
| `make_deposition_cif.py` | Add `struct_conn`, preserving rows already present. Edit its `LINKS` table. |

If the `scripts/` folder is not reachable, stop and say so before starting. Never edit the
reference copies; adapted copies go in the project. Run them wherever gemmi and
RDKit are available (`phenix.python` bundles both) and report their own output.

There is no script for link dictionaries or params — write those by hand
following Step 2a.

### 2. Full ligand built into place first, then split

**Split only from a CIF model in which the complete, current ligand has been
fitted as one residue into the density**, using the dictionary that will be
deposited. The split carves components out of that fitted ligand, so every piece
lands where the ligand was built.

Do not split from: a ligand fitted with an older or different definition of the
compound (even in the same dataset), dictionary or SMILES coordinates, a ligand
transplanted from another dataset, or hand-placed monomers. Hand-placed fragments
give junctions ~1 A too long, which nothing will link, and they fly apart on
refinement.

**If that model does not exist, stop and ask for it.** Step 1 can proceed while
waiting.

### 3. File layout

Alongside Phenix's `Refine_NN/` job folders, never inside them:

```
<project>--phenix/
  Refine_01/ ... Refine_NN/     Phenix's own; read only
  Split-Ligand/                 ligand_split.cif, adapted split script, logs
  Restraints/                   component CIFs, the combined link CIF, both params files
  Deposition/                   the finished CIF and the deposition notes
```

**Component dictionaries the user generates in Grade2 sit according to their
scope.** A set shared by two or more crystals lives **one level above the project
folders**; a set specific to one crystal lives inside that project. Either way a
Grade2 run folder keeps its own fixed CIF beside the raw output, and the copy
Phenix is pointed at is copied into that project's `Restraints/`. **So when
several datasets share a ligand, look one level up before asking for SMILES or
new Grade2 runs** — the set is probably already there.

**Everything Claude generates is per-project, always** — the split, the combined
link CIF, both params files. They carry coordinates, chains, residue numbers and
sigmas, none of which transfer between structures even for an identical ligand.
Copying a link CIF in as a starting point is fine, but it becomes that structure's
own file the moment a sigma is tightened, because two datasets reach "close enough
to tighten" at different rounds.

Repaired Coot models and stripped Phenix outputs go where they feed, never back
into a `Refine_NN/`. Use per ligand-and-dataset subfolders when a project has more
than one.

State where each file went, in plain words, as it is written.

Projects already running under the older single-folder convention (`components/`
or `Claude-Unnaturals` holding everything) keep it; renaming would break paths in
adapted scripts and live Phenix jobs. Likewise, a link file already wired into a
running job may stay put for that job.

## Step 0 — project overview

Lay out what is being done and confirm it: target, ligand, resolution, which
residue the warhead reacts with, what already exists (fitted ligand? Grade2
output, in this project or one level up? refinement history?), where the folders
are, and what the deliverable is. A wrong assumption here is cheap; the same one
found at Step 3 is not.

**First, define the chains.** Before anything else, ask which protein chain the
ligand binds and which chain the ligand is (or will be) in. With several copies
in the asymmetric unit, residue names and numbers repeat across chains, so every
later step — the split, link definitions, distance checks and `struct_conn`
rows — must name the chain explicitly. Record both and use them throughout:
`LIGAND`, `OUT_CHAIN` and `ANCHOR` in `split_ligand.py`, the chain selections in
link params, and the chain columns of `LINKS` in `make_deposition_cif.py`.

## Step 1 — establish the full ligand

1. **Grade2 on the whole ligand**, so there is a dictionary to build and refine
   the intact molecule with. Confirm it exists before going on.
2. **Decide the decomposition now**, say why, and **give the user the SMILES for
   each building block** so they can run Grade2 on them in parallel.
3. **Run `fix_grade2_cif.py` on everything Grade2 returns**, whole ligand
   included.
4. **The user builds the full ligand in Coot and refines in Phenix.** Wait.

**Decomposing.** One component per backbone — anything with N-CA-C-O becomes a
residue, and side-chain decorations stay inside their residue however large (a
warhead arm is not a residue). Caps are components too, as ACE and NH2 are in the
CCD; if a CCD component already combines a residue and its cap, use it. Generate
restraints from the complete free molecule, then omit leaving atoms from the
coordinates — Phenix silently drops restraints referencing absent atoms.

Use 5-character placeholder IDs; wwPDB reassigns at OneDep. **This commits the
whole pipeline to mmCIF** — legacy PDB truncates to three. Warn the user.

**Grade2.** Grade2 needs a licence from Global Phasing. Unless the user has
installed it locally with their own licence, they run it on the free Grade Web
Server (http://grade.globalphasing.org) and upload the results. Ask which before
starting, and do not try to run `grade2` yourself without a local install. The
command-line options below are for a local install; the web server may not offer
them, so with server output check the points below in the returned CIF instead.

One run per component, `--no_charging --aa_loose`, submitting the
neutral free molecule (COOH, NH2) with any leaving group present. Verify from the
log that `--aa_loose` picked the alpha-amino N and not a side-chain amide N; it
should not fire on acyl caps, and if it does rerun that one with
`--no_aa_labels`. **Read the `CHECK:` line** — Grade2 tests the CCD by InChIKey,
and if a component already exists you must use the existing ID or deposition
rejects the duplicate.

**The CIP trap.** Do not infer L/D from the CIP letter. When the CB substituent
outranks the carboxyl — a ring sulfur, as in thioproline or cysteine — the
descriptor inverts, so L is (R) and (S) means D. A ligand whose centres are all
*S* can still contain one D residue. Check each centre against its own reference
compound and against the synthesis building block. Re-derive from the user's
original SMILES rather than from what you typed, and show the comparison. Getting
this wrong bakes a wrong chiral volume into everything downstream.

**Fixing the CIFs** (`fix_grade2_cif.py`, not by hand). Set `_chem_comp.group`
and `_chem_comp.type` to `L-peptide`/`L-PEPTIDE LINKING`,
`D-peptide`/`D-PEPTIDE LINKING`, or leave acyl caps `non-polymer`/`NON-POLYMER`.
Set `_chem_comp.name` — Grade2 writes the InChIKey there and OneDep reads that
field. Rename only a cap's acyl atoms to `C`, `O`, `OXT`, `HXT` so the link uses
peptide-bond names; **leave side-chain names as Grade2 wrote them**, since wwPDB
remediates naming and every rename risks desynchronising the bond, angle, torsion
and plane loops. Check the script's `verify: OK` line.

## Step 2 — split the fitted ligand

Prerequisite: rule 2, repaired with `fix_coot_cif.py` if it came from Coot.

Run the adapted `split_ligand.py`. It perceives connectivity from the fitted
coordinates, builds each component's heavy-atom graph from its fixed CIF minus
leaving atoms, matches with bond orders ignored, and solves an exact cover so
every heavy atom is claimed once. **A unique cover is proof the assignment is
right, not a guess**; if no cover is found the component list or the leaving atoms
are wrong, and the script says so rather than producing a plausible wrong answer.
Coordinates never move. Output `ligand_split.cif` to `Split-Ligand/`.

**Validate before refining:** every junction 1.30-1.40 A, atom count unchanged,
maximum coordinate deviation 0.000000, each residue's chirality checked against
the input. Handle every copy of the ligand to be modelled.

Write the peptide link dictionary and params (format in Step 2a) into
`Restraints/`, and drop stale whole-ligand dictionaries from the run.

**If the ligand is covalently reactive, do not add the warhead link yet** — get
the peptide links in and refined first, then Step 2a.

**Link detection is distance-based and happens once**, at cutoffs near
1.9-2.2 A. Residues further apart never link, and refinement will not bring them
back.

The user builds in Coot and refines. Wait.

## Step 2a — the covalent link (reactive ligands only)

**Block structure is where this goes wrong.** A CCP4 link dictionary needs
`data_link_list` (the `_chem_link` summary loop), `data_mod_list` (the
`_chem_mod` summary loop **only**), **one `data_mod_<id>` block per
modification** holding its `_chem_mod_atom` loop, and `data_link_<id>` with the
bond, angle and plane loops. Putting `_chem_mod_atom` inside `data_mod_list`
gives `KeyError: '<mod id>'`. Use short hyphen-free mod ids.

**One combined `.cif` for every link, never per-stage files.** Which restraint
files reach a run is set by hand each time, and a file left out is dropped with
no error — only a quiet "Not linked" in the log. This has happened: a staged
warhead link was added as its own file, the run carried only that file, and the
previously-working peptide links went unrestrained for a round. Merge a new link
into the existing combined set, and when tightening a sigma edit that one file in
place.

**Two params files, one per program**, each with a header saying which:
`apply_links.params` with the inner `pdb_interpretation { }` only, for
`phenix.pdb_interpretation`; `apply_links_refine.params` with the same content
inside `refinement { }`, for `phenix.refine`. The wrong one gives "PHIL
parameters are not recognized". Always the nested block form, one
`apply_cif_link` sub-block per link:

```
refinement {
  pdb_interpretation {
    apply_cif_link {
      data_link = LINKID1
      residue_selection_1 = chain C and resname AAAAA and resid 1 and name C
      residue_selection_2 = chain C and resname BBB and resid 2 and name N
    }
    apply_cif_link { ... }
  }
}
```

not repeated one-line `refinement.pdb_interpretation.apply_cif_link { }` scopes.
A params file is needed once there are two links regardless — on the command line
the second overwrites the first.

**Mods.** They delete leaving atoms at link time. Keep the mod that deletes a
**hydrogen on a heteroatom that becomes substituted** even when that atom is
already absent from the coordinates — riding-H addition otherwise re-creates it
and over-coordinates the atom. With heavy-atom-only split coordinates Phenix logs
the delete as "skipped"; that is expected, since the mod is what tells the
residue's topology to expect one fewer substituent. If the model *does* still
carry such a hydrogen, delete it from the coordinates too or Phenix cannot type
it.

**Geometry targets.** For the ligand's own inter-residue junctions, the
whole-ligand Grade2 dictionary already holds Mogul (CSD) targets for exactly
those bonds, angles and planes — map its atom names to component names by
coordinates and reuse them. For a warhead link, use published small-molecule
structures of the right linkage, not the PDB: protein precedents are few and
their distances track each depositor's dictionary rather than chemistry. Match
the linkage type properly — a Tyr-OH sulfonate ester is an aryl-aryl ester, not
an alkyl one. State the targets and their provenance. Leave genuinely soft
torsions unrestrained (planes keep amides flat without choosing cis or trans).

**Standard residue pairs may already be linked for you.** Two consecutive
standard-type peptide-linking residues (both recognised CCD/GeoStd components)
are usually linked automatically by Phenix as a TRANS bond at its own default
(~1.329 A). Define the explicit link as a fallback but **leave its
`apply_cif_link` commented out** — applying it alongside the automatic bond
restrains the same bond twice. Confirmed: with only the flanking custom links
applied, the `.geo` showed exactly one bond at that junction and `Link IDs`
included `TRANS`. With Phenix's GUI automatic linking on, a more specific
supplied definition wins over the generic default — also confirmed, one restraint
at the custom Mogul target with the `apply_cif_link` still commented out. Check
the `.geo` for the specific case rather than assuming either way.

**Do not judge automatic linking, or an external CCD/GeoStd dictionary, from an
isolated test.** A minimal reproduction on the ligand chain alone, without the
full protein and monomer library, can show two standard residues unlinked when
the real full-model run links them fine — the automatic pass needs that context.
And a dictionary from GeoStd on GitHub may not match what phenix.refine loads
from its own bundled library (one component's backbone amine H is `HN` in one and
`H` in the other). Trust the real run's `.log` and `.geo`.

**Stage the restraint.** Compute `(delta/sigma)^2` first — if one restraint
carries a large fraction of total bond strain it will win against everything and
distort the model. Refine soft (sigma ~0.10) until the geometry is close, then
tighten. Ordinary macrocycles do close gaps of ~1 A over a few runs on their own;
simulated annealing helps further on a far-from-target restraint but is not
required to see movement. Once tightened, confirm the tightened numbers actually
reached the run before judging whether they helped.

**Before forcing a bond that will not form, check the approach direction.**
Compute the vacant coordination direction at the electrophile as the negative sum
of unit vectors to its existing substituents, and measure the angle to the
nucleophile. Near 180 degrees means the reactive face is turned away and no
restraint will help — rotate the group about its attachment bond first and report
the best achievable alignment. A bond that does close tends to start within ~15
degrees of aligned.

Validate, then wait for the user's build and refinement.

## Following the refinement

Keep a running record across rounds — round number, what changed, R-work, R-free,
the restrained bond length **measured from coordinates**, clashscore — and carry
it forward. Thirty rounds in, "which run introduced this" cannot be answered
retrospectively.

**Two traps, both silent, both confirmed to have cost runs:**

**1. `.geo`'s printed values are stale.** The `.eff` records what was requested;
the `.geo` shows the restraint actually built — but its per-atom model
bond/angle/torsion values are computed once, right after `pdb_interpretation`
builds the model from the *input* coordinates, before the macrocycles run.
Timestamps confirm it: `.geo` is written minutes before the `.log` finishes. For a
junction that started near its target this is harmless. For one that started far
away it is actively misleading — a warhead bond read an identical `.geo` value
across three consecutive rounds while really closing from 2.9 to 2.0 A, which
produced a wrong "it isn't moving" conclusion. **Use `.geo` for topology only**
(exactly one bond restraint per junction — zero means it failed to link, two means
double-restrained — plus which mods and `Link IDs` applied), and **measure any
far-from-target restraint from the output coordinates** (gemmi `read_structure`,
`Position.dist`).

**2. A Phenix output carries its own dictionaries.** Every `comp_*` and `link_*`
block Phenix used is written into its output mmCIF. Reuse that file as the next
round's input and those definitions are re-read and override the link CIF given
on the command line, with no error — you tighten a sigma and the `.geo` keeps
showing the old one. **Run `check_embedded_restraints.py` before pointing
phenix.refine at any Phenix or Coot output as input, every round**, and
`strip_phenix_restraints.py` when it fires. This has happened on two different
projects. A Coot round trip strips them as a side effect; do not rely on that.

## Coot

Coot draws inter-residue bonds on distance alone — **a drawn bond is not a
restrained one**. It will not link 5-character components it has no link
definition for, and its Refine can pull them apart. Use Edit Chi Angles and the
torsion tools, which change coordinates without building restraints, and leave
refinement to Phenix. For edits Coot cannot make on non-standard residues, such
as changing an omega, tell the user instead of writing a script for it.

**Coot writes newly added or rebuilt residues to mmCIF with invalid `label_*`
fields** while `auth_*` stays correct: `label_atom_id` `.` (a nameless atom),
`label_comp_id` `---`, `label_asym_id` `.`, `label_seq_id` `-32767`. Phenix reads
`label_*`, so it aborts with **"unknown nonbonded energy type"**. A nameless atom
is the usual trigger and easy to miss, because the file looks fine in Coot and
ChimeraX — both read `auth_*`. **Run `fix_coot_cif.py` on every Coot save.** It is
safe on an already clean file, so run it routinely instead of diagnosing. Note
that Coot's own `label_seq_id` numbering skips the residues it added, so the
repair renumbers the polymer — but only where the file declares no reference
sequence, since `_entity_poly_seq` means `label_seq_id` refers to the full
construct.

That Phenix error has a second, unrelated cause: a component whose restraint CIF
was never supplied. If the model is clean and Phenix still fails, check coverage.

## Step 3 — deposition

When the final refinement is done, the user may need `make_deposition_cif.py`
to clean up that model for deposition: it strips the restraint blocks Phenix
embeds and adds the covalent-bond `struct_conn` rows. Edit its `LINKS` table
for this ligand first.

Make `Deposition/` and copy every file the deposition needs into it.

**`struct_conn`: read what is there, keep every row, append only what is
missing.** Phenix *does* write `struct_conn` — both the links you applied and
bonds it found itself, notably crystallographic-symmetry disulfides. Never
rebuild the loop: a script using `init_loop` once replaced it wholesale and
silently dropped a symmetry disulfide from a finished deposition file.

Query it one tag at a time. Asking gemmi's `find()` for several tags at once
returns an **empty table if any single one is missing**, which reads as "no rows"
and is not — Phenix omits `pdbx_dist_value`, and that false negative is what
caused the mistake above. `make_deposition_cif.py` does this correctly and also
skips a link Phenix already recorded, so no duplicate row is added.

**Phenix writes no `ptnr2_symmetry`**, so a symmetry bond arrives looking like a
chemically impossible self-bond — same chain, residue and atom on both sides.
Either add the operator in wwPDB notation or flag it in the notes.

Include the protein-ligand covale row always, and inter-residue covale rows for a
multi-residue ligand: OneDep gives each distinct chemical component its own
entity regardless of connectivity, so **a multi-residue ligand is deposited as
several entities, never merged into one**, whatever the dictionary says. The
covale rows are the only record that the fragments were one linked construct. Set
`pdbx_leaving_atom_flag`; measure distances from the coordinates.

**Keep the leaving group in the dictionary, out of the coordinates.** The
component describes the compound as soaked, the bound instance omits the leaving
atoms, the covale record explains why. Precedent: DU1 in 5UEN, KBH in 6RNU.

**Strip the embedded restraint blocks**, then **verify**: atom count unchanged,
maximum coordinate deviation 0.000000, single data block,
`check_embedded_restraints.py` clean, and every row that was there before still
there.

**Write deposition notes** in `Deposition/` beside the CIF, covering: that the
ligand's residues arrive as separate entities despite being covalently bonded,
with the covale records as the evidence; any D-configured residue and why its
descriptor looks inverted; the covalent bond and its independent confirmation;
any symmetry bond lacking an operator; and any validation outlier with the
evidence that justifies it. Annotators read these.

## Judging the result honestly

When a restrained bond refuses to reach its target, **first re-measure it from
the output coordinates, not `.geo`** — it may be moving and only look static. If
it genuinely is not, consider in order: the restraint is still soft; a stale
embedded dictionary is overriding the link CIF; the geometry is wrong-faced; the
chemistry did not happen. Say which, and say when the evidence is thin. **Do not
conclude a reaction did not occur while an unexplained soft or overridden
restraint is still in play** — that conclusion has been reached wrongly, argued at
length, and contradicted by mass spectrometry.

**R-free is the wrong tool for local questions.** Refinement absorbs misplaced
atoms by inflating their B-factors, rescuing global R while leaving the local fit
wrong. For a single bond or amide use real-space CC per residue, B-factors
relative to neighbours, and clashscore.

To test an alternative conformation properly, build it, **re-fit dependent atoms
independently** (flipping an omega drags the whole side chain), refine to
convergence, and compare local metrics. Report which model wins and on what.

If independent evidence exists — mass spectrometry for a covalent adduct — weigh
it above the crystallography at moderate resolution, and say so.

A residual gap between the final geometry and the literature target (~0.1 A on a
warhead bond) is not automatically a failure — say so plainly alongside the
number, and note whether the surrounding junctions refined normally. Whether to
keep iterating or accept is the user's call.

## Tools

The Skill-folder scripts first (rule 1), outputs in the project folders (rule 3).
Underneath: gemmi for CIF and coordinate work, RDKit for graph matching and
stereochemistry. For any file operation, verify atom count and maximum coordinate
deviation before and after, and report both.