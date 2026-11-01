#!/usr/bin/env python3
"""
build_zips.py -- build claude.ai upload ZIPs into dist/.

    python3 tools/build_zips.py            # every skill
    python3 tools/build_zips.py <name>     # one skill

Runs validate.py first and refuses to build if it fails. Each ZIP holds
<name>/SKILL.md and the rest of the skill folder at its root.
Standard library only.
"""
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate import ROOT, SKILLS, validate  # noqa: E402

DIST = ROOT / "dist"
EXCLUDE_NAMES = {".DS_Store"}
EXCLUDE_PARTS = {"__pycache__"}


def wanted(path):
    if path.name in EXCLUDE_NAMES or path.name.startswith(".git"):
        return False
    return not (EXCLUDE_PARTS & set(path.parts))


def build(folder):
    DIST.mkdir(exist_ok=True)
    out = DIST / f"{folder.name}.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(folder.rglob("*")):
            if f.is_file() and wanted(f.relative_to(folder)):
                z.write(f, Path(folder.name) / f.relative_to(folder))
    return out


def human(n):
    return f"{n / 1024:.1f} KB" if n < 1024 * 1024 else f"{n / 1024 / 1024:.1f} MB"


def main(args):
    if not validate():
        print("\nValidation failed -- nothing built.")
        return 1
    if args:
        folders = [SKILLS / a for a in args]
        missing = [f.name for f in folders if not f.is_dir()]
        if missing:
            print(f"No such skill: {', '.join(missing)}")
            return 1
    else:
        folders = sorted(p for p in SKILLS.iterdir() if p.is_dir())
    print()
    for folder in folders:
        out = build(folder)
        print(f"built {out.relative_to(ROOT)}  ({human(out.stat().st_size)})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
