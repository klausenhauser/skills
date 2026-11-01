#!/usr/bin/env python3
"""
new_skill.py -- start a new skill from template/.

    python3 tools/new_skill.py <name>

Copies template/ to skills/<name>/ and fills in the name. Standard library only.
"""
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main(args):
    if len(args) != 1:
        print(__doc__.strip())
        return 1
    name = args[0]
    if not re.fullmatch(r"[a-z0-9-]{1,64}", name) or re.search(r"claude|anthropic", name):
        print("Name must be lowercase letters, digits and hyphens, at most 64 characters,")
        print("and must not contain 'claude' or 'anthropic'.")
        return 1
    dest = ROOT / "skills" / name
    if dest.exists():
        print(f"skills/{name} already exists.")
        return 1
    shutil.copytree(ROOT / "template", dest)
    skill_md = dest / "SKILL.md"
    skill_md.write_text(skill_md.read_text().replace("SKILL-NAME", name))
    print(f"Created skills/{name}/")
    print("Next:")
    print(f"  1. Edit skills/{name}/SKILL.md (description: what it does AND when to use it).")
    print(f"  2. Add \"./skills/{name}\" to one plugin's \"skills\" list in .claude-plugin/marketplace.json.")
    print("  3. Run: python3 tools/validate.py")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
