# xray-mr-phaser

Solves a crystal structure by molecular replacement (MR) with Phenix and Phaser.
It works for any protein, and is built in particular for complexes whose
components share a fold, where Phaser can place the wrong component at a site.

Claude reads [`SKILL.md`](SKILL.md) for the full instructions. This page is the
human overview.

## How the work is shared

Claude runs the pipeline itself in a new `MR_pipeline/` folder in your project:
Xtriage, model preparation, the Phaser searches and the test refinements. It
asks you before choosing between datasets and before killing a run that isn't
clearly beaten. It never writes into your existing Phenix GUI job folders.

## The logic

1. **Get the composition from Xtriage, never from a guess.** Xtriage with the
   real sequences gives the copies per asymmetric unit, the anisotropy, tNCS and
   twinning indicators that every later step depends on.
2. **Audit the search models before trusting them.** An AlphaFold 3 complex with
   ipTM below about 0.5 has no reliable interface, so its components are searched
   separately. Low-confidence regions are trimmed, and supplied PDBs are compared
   domain by domain to the predictions to spot alternative domain arrangements.
3. **Let similar folds compete.** When components share a fold, Phaser can put
   the wrong one at a site. The key run offers it "model A *or* model B"
   at every placement, so it picks the better fit. Single-model runs in parallel
   act as independent cross-checks: the same placement found by different runs is
   strong evidence.
4. **Prefer multi-domain search models.** The arrangement of several domains is
   specific to one component, whereas a single domain of a shared fold fits
   almost anywhere.
5. **Validate every candidate by refinement.** A correct, unrebuilt 4 Å
   solution refines to R-free around 0.35–0.40; 0.50 or more means wrong.
   Dropping a doubtful copy and re-refining shows whether it is real. Crystal
   contacts and glycosylation sites check which component is which.

## Steps

| Step | What happens |
|---|---|
| 0. Inputs | Claude collects the MTZ, sequences, candidate models and the Phenix location, asking only for what it can't find. |
| 1. Datasets and composition | Cleans the MTZ to one F/SIGF/FreeR set and runs Xtriage on each dataset. Records cell, space group, resolution, anisotropy, tNCS, twinning test and copies per ASU, then ranks the datasets with you. |
| 2. Search models | Reads the AlphaFold 3 confidence scores, compares supplied models to predictions, trims low-confidence parts and cuts domain-level models (each domain, domain pairs, multi-domain units, full chain). |
| 3. Phaser | Writes and runs the competition search plus single-model cross-checks in parallel, then summarises TFZ, LLG and packing for each run. |
| 4. Validation | Test-refines each candidate with and without its weakest copy, maps the contacts including symmetry mates, checks identities (e.g. at N-glycosylation sites), then refines the best solution further. |
| 5. Write-up | Writes `MR_pipeline/README.md` (result, Xtriage table, model audit, strategy, runs, next steps) and copies the final model and MTZ to `results/`. |

## Scripts

All in [`scripts/`](scripts/). They run with `phenix.python` or plain Python 3
and need no editing.

| Script | Use |
|---|---|
| `mk_phaser.py` | Writes a Phaser keyword file for one run. `a\|b` in a search step means "a or b". See its header for usage. |
| `summ.sh` | Summarises finished Phaser runs: completion, best LLG and placements. |
| `compare_models.py` | CA RMSD between two models, with an optional residue-number offset. |
| `contacts.py` | Counts CA–CA contacts under 8 Å between chains, including symmetry mates. |

## Requirements

- **Phenix** (with Phaser), installed by you. Claude finds `phenix_env.sh` and
  sources it.
- **Coot**, optional, to look at the result.
- Works best in **Claude Code** on the machine where Phenix lives. On claude.ai,
  Claude can explain the strategy and write the input files, but you run Phenix
  yourself.
