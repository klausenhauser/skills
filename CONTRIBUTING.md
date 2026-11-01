# Contributing

Suggestions and fixes are welcome.

- **Found a problem?** Open an issue using the bug report template. Say which
  skill, what you asked Claude, and what went wrong.
- **Proposing a change?** Fork the repo, edit the skill, run
  `python3 tools/validate.py`, and open a pull request.
- **Adding a skill?** Run `python3 tools/new_skill.py <name>` to start from the
  template, then add it to a plugin in `.claude-plugin/marketplace.json`.

Rules for every skill:

- Everything the skill needs lives inside its own folder. No references to
  other skills or to files outside the folder.
- The description says what the skill does **and when to use it**.
- No secrets, personal details, internal URLs or unpublished data. This repo is
  public.
- Tools use the Python standard library only.
