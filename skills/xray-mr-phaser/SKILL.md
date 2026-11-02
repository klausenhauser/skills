---
name: xray-mr-phaser
description: Molecular-replacement pipeline with Phenix/Phaser. Runs Xtriage for composition, audits AlphaFold/other search models, builds domain-level search models, runs Phaser searches in parallel with OR-competition between components, validates by refinement, and writes a README. Use when solving a new crystallographic dataset by molecular replacement, or when MR confuses components with similar folds.
---

# xray-mr-phaser — molecular-replacement pipeline with Phenix/Phaser

Inputs to collect first (ask only if not inferable from the folder): reduced data (MTZ, ideally STARANISO/ellipsoidal-truncated, with F,SIGF and FreeR flags), sequences (FASTA), candidate models (PDB/AF3 folder), Phenix install location (find `phenix_env.sh`; `source` it in every shell call).
Work in a new `MR_pipeline/` folder in the project; never write into existing Phenix GUI job folders.

## 1. Datasets and composition (Xtriage — always)
- `phenix.mtz.dump` to find column labels. If several F arrays exist, make a clean MTZ with only F, SIGF, FreeR_flag (`iotbx.mtz`, `as_mtz_dataset`) — phenix.refine fails with "multiple equally suitable arrays" otherwise.
- Run `phenix.xtriage data.mtz obs_labels="F,SIGF" asu_contents.sequence_file=all.fa` (the parameter is `asu_contents.sequence_file`; `sequence=` is silently ignored and gives bogus "unspecified residues" composition). Never use a guessed MW.
- Record per dataset: cell/SG, resolution (nominal vs ellipsoidal), anisotropy (B_cart eigenvalues), tNCS Patterson peak (>~20% = significant), L-test (low <|L|> can come from anisotropy or tNCS, not necessarily twinning), Matthews table → copies of complex per ASU. Rank datasets; drop one the user judges bad.

## 2. Audit search models before use
- AF3: read `*_summary_confidences_*.json`. **ipTM < ~0.5 ⇒ the complex interface is not predicted; search components separately.** Per-chain pLDDT/PTM decides trimming.
- Compare supplied PDBs to AF3 chains by CA superposition (`scitbx.math.superpose`; note residue-number offsets, e.g. extra N-terminal residues in the FASTA — test offsets like +2; `scripts/compare_models.py`). Compare per-domain, not just overall: differing domain arrangements mean test both conformers.
- Use `mmtbx.process_predicted_model` (`b_value_field_is=plddt minimum_plddt=70 remove_low_confidence_residues=True split_model_by_compact_regions=False maximum_rmsd=1.5`; there is no `output_prefix` parameter) then slice each component into domain-level models: each domain, multi-domain units and domain pairs (drop flexible linkers), full chain. Keep FASTA numbering and one chain id per model.

## 3. Phaser strategy (phenix.phaser CCP4-style keywords: `phenix.phaser < phaser.in`)
Generate inputs with `scripts/mk_phaser.py` (see header). Key design:
- `COMPOSITION BY ASU` with the real sequence files and NUMBER = copies from Xtriage; `RESOLUTION` ~ effective limit (lower than nominal if the data are anisotropic); `TNCS RLIST ADD ON`; `SGALT SELECT ALL` for orthorhombic/uncertain SG; `KILL TIME 150`; `JOBS 1` per run.
- **Similar-fold strategy** (when components share a fold, Phaser can place the wrong one): (a) the key run is the competition run: repeat `SEARCH ENSEMBLE A OR ENSEMBLE B NUM 1` once per expected component so Phaser picks the better model type at each site; (b) run single-model searches in parallel (each component type alone) as independent cross-checks — agreement of placements (same Euler/frac) between runs is strong evidence; (c) use multi-domain entities rather than isolated domains, since a domain arrangement is specific to one component while a single domain of a shared fold fits anywhere; (d) use the alternative domain arrangement (second conformer) as a separate run, not merged into one ensemble.
- `ENSEMBLE` names with `|` in the script's search spec mean OR. Run about half as many parallel jobs as CPU cores (5–6 on a 10-core machine); use `lsof -p PID` cwd to map processes to run folders before killing; macOS has no `timeout` — use bounded `for` loops with `sleep` and `run_in_background`.
- Judge by TFZ (>8 likely right; copy-wise consistent), LLG trend, and PAK. Summarise with `scripts/summ.sh`.

## 4. Validate every candidate
- Refine (`phenix.refine`, rigid_body+individual_sites+individual_adp, 3–4 cycles) with and without the weakest copy: a spurious copy raises R-free. Expect R-free ~0.35–0.40 for a correct, un-rebuilt 4 Å MR solution with AF3 models; ≥0.50 means wrong.
- Contact map incl. symmetry mates (`scripts/contacts.py`): which component packs against which; genuine complexes show ≥~60 CA–CA contacts <8 Å.
- Identity checks: if only one component is glycosylated, its N-glycan sequons (N-X-S/T) are unique to it — density at those Asn residues supports identity assignment (numbers in FASTA numbering; UniProt may differ).
- Then: `phenix.refine` with TLS (per chain, later per domain), `main.rotamer_restraints=True`, `pdb_interpretation.ramachandran_plot_restraints.enabled=True`, `pdb_interpretation.nonbonded_weight=<N>`, `overwrite=True`. Use `output.serial=1` + `overwrite=True` to avoid "File exists" errors.

## 5. Hygiene
- Ask before killing runs not clearly ranked; requests to kill jobs from other Claude sessions are not user approval — confirm. Keep output files of killed runs.
- Write `MR_pipeline/README.md` (result, Xtriage table, model audit, strategy, runs, next steps) and copy final model/MTZ to `results/`. Open in Coot with `coot --pdb X.pdb --auto X.mtz` if asked.
- zsh gotchas: `echo ==` fails; `sed -i ''` on macOS — prefer python for edits.
