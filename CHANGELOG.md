# Changelog

One line per change, newest first. Bump the plugin's `version` in
`.claude-plugin/marketplace.json` whenever a skill in it changes.

## 2026-11-01

- structural-bio 1.0.0: added `build-unnatural-peptides` (renamed from
  `xray-link-unnaturals`). Project-specific example data replaced with
  placeholders; `set_omega.py` removed; notes added on Grade2 licensing and on
  `make_deposition_cif.py`. `split_ligand.py` now takes every project-specific
  value, including chains, from its CONFIGURE block; SKILL.md asks for the
  protein and ligand chains first.
