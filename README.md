# skills

Structural biology and lab skills for Claude, usable on claude.ai and in
Claude Code.

| Plugin | Skills |
|---|---|
| `structural-bio` | [`build-unnatural-peptides`](skills/build-unnatural-peptides/README.md) |

Each skill's own `README.md` (in `skills/<name>/`) explains what it does and
its steps.

## Install

**Claude Code** (no GitHub account needed; this repo is public):

```
/plugin marketplace add klausenhauser/skills
/plugin install <plugin>@klausenhauser-skills
```

Get updates later with `/plugin marketplace update klausenhauser-skills`.

**claude.ai:** download the skill's ZIP from the
[latest release](https://github.com/klausenhauser/skills/releases/latest) and
upload it in your claude.ai settings under Skills. Members of our Team account
get the skills from the team admin and don't need to do anything.

---

## Maintainer recipes

All commands run from the repo folder.

### Add a skill

1. `python3 tools/new_skill.py <name>`
2. Edit `skills/<name>/SKILL.md`. The description must say what the skill does
   **and when to use it**.
3. Add `"./skills/<name>"` to one plugin's `skills` list in
   `.claude-plugin/marketplace.json`.
4. `python3 tools/validate.py`
5. Commit.

### Update a skill

1. Edit the files in `skills/<name>/`.
2. Bump that plugin's `version` in `.claude-plugin/marketplace.json`
   (e.g. `1.0.0` → `1.0.1`). Claude Code users only get the change when this
   number changes.
3. Add a line to `CHANGELOG.md`.
4. `python3 tools/validate.py`, then commit and push.
5. Deploy to claude.ai (below).

### Deploy to the team on claude.ai

1. Get the ZIPs, either way:
   - push a version tag (`git tag v1.0.1 && git push --tags`); GitHub builds
     the ZIPs and attaches them to a Release, or
   - build them locally: `python3 tools/build_zips.py` (or
     `python3 tools/build_zips.py <name>` for one); they land in `dist/`.
2. In claude.ai → Organization settings → Skills, upload each **changed** ZIP.
3. If a skill was renamed, delete its old version there too.

### My personal account

Upload the same ZIP in your personal claude.ai settings if you want it there.

### Outside contributors

Point them to [CONTRIBUTING.md](CONTRIBUTING.md).

## Layout

```
skills/<name>/          one folder per skill (SKILL.md, scripts/, reference/)
.claude-plugin/         marketplace.json: groups skills into plugins
template/               starting point for new skills
tools/                  validate.py, build_zips.py, new_skill.py
dist/                   built ZIPs (not committed)
inbox/                  unsorted imports (not committed)
```

## Licence

[MIT](LICENSE)
