#!/usr/bin/env python3
"""
validate.py -- check every skill in skills/ and print a pass/fail table.

    python3 tools/validate.py

Exits 1 if any skill fails, so it can gate build_zips.py and CI.
Standard library only.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKILLS = ROOT / "skills"
MARKETPLACE = ROOT / ".claude-plugin" / "marketplace.json"

MAX_NAME = 64
MAX_DESC = 1024
BIG_FILE = 10 * 1024 * 1024
TEXT_SUFFIXES = {".md", ".py", ".txt", ".json", ".yaml", ".yml", ".sh", ".cfg",
                 ".toml", ".ini", ".html", ".css", ".js", ".ts", ".csv", ".tsv",
                 ".cif", ".eff", ".params", ".xml"}
WHEN_PHRASES = re.compile(r"\buse (this )?(when|for|whenever|before|after)\b|\bwhen (asked|the user)\b", re.I)
SECRETS = [
    ("Anthropic API key", re.compile(r"sk-ant-[A-Za-z0-9_-]{20,}")),
    ("OpenAI-style key", re.compile(r"\bsk-[A-Za-z0-9]{32,}")),
    ("GitHub token", re.compile(r"\b(ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}")),
    ("AWS access key", re.compile(r"\bAKIA[0-9A-Z]{16}\b")),
    ("Google API key", re.compile(r"\bAIza[0-9A-Za-z_-]{35}\b")),
    ("Slack token", re.compile(r"\bxox[abprs]-[A-Za-z0-9-]{10,}")),
    ("private key", re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----")),
    ("hard-coded password", re.compile(r"(password|passwd|api_key|apikey|secret|token)\s*[:=]\s*['\"][^'\"\s]{6,}['\"]", re.I)),
]
HOME_PATH = re.compile(r"/Users/(?!<)[A-Za-z0-9._-]+|/home/(?!<)[A-Za-z0-9._-]+")


def parse_frontmatter(text):
    """Minimal YAML frontmatter reader: top-level `key: value` pairs, quoted or
    plain, plus `>`/`|` block scalars. Enough for SKILL.md; no PyYAML needed."""
    if not text.startswith("---"):
        return None
    lines = text.splitlines()
    try:
        end = next(i for i, l in enumerate(lines[1:], 1) if l.strip() == "---")
    except StopIteration:
        return None
    data, key, block = {}, None, None
    for line in lines[1:end]:
        if block is not None and (line.startswith((" ", "\t")) or not line.strip()):
            block.append(line.strip())
            data[key] = (" " if style == ">" else "\n").join(b for b in block).strip()
            continue
        block = None
        m = re.match(r"^([A-Za-z_][\w-]*):\s*(.*)$", line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if val in (">", "|", ">-", "|-"):
            block, style = [], val[0]
            data[key] = ""
        elif len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
            data[key] = val[1:-1]
        else:
            data[key] = val
    return data


def marketplace_paths():
    """Return {skill folder name: [plugin names]} and a list of errors."""
    if not MARKETPLACE.exists():
        return {}, ["missing .claude-plugin/marketplace.json"]
    try:
        mp = json.loads(MARKETPLACE.read_text())
    except json.JSONDecodeError as e:
        return {}, [f"marketplace.json is not valid JSON: {e}"]
    listed, errors = {}, []
    for plugin in mp.get("plugins", []):
        pname = plugin.get("name", "?")
        for path in plugin.get("skills", []):
            p = (ROOT / plugin.get("source", "./") / path).resolve()
            if not (p / "SKILL.md").exists():
                errors.append(f"plugin '{pname}' lists {path}, which has no SKILL.md")
            if p.parent == SKILLS:
                listed.setdefault(p.name, []).append(pname)
    return listed, errors


def check_skill(folder, listed):
    errors, warnings = [], []
    skill_md = folder / "SKILL.md"
    if not skill_md.exists():
        return None, ["no SKILL.md"], warnings

    fm = parse_frontmatter(skill_md.read_text())
    if fm is None:
        return None, ["SKILL.md has no YAML frontmatter"], warnings
    name, desc = fm.get("name", ""), fm.get("description", "")

    if not name:
        errors.append("frontmatter has no name")
    else:
        if name != folder.name:
            errors.append(f"name '{name}' does not match folder '{folder.name}'")
        if not re.fullmatch(r"[a-z0-9-]+", name):
            errors.append("name must be lowercase letters, digits and hyphens")
        if len(name) > MAX_NAME:
            errors.append(f"name longer than {MAX_NAME} characters")
        if re.search(r"claude|anthropic", name):
            errors.append("name must not contain 'claude' or 'anthropic'")

    if not desc:
        errors.append("frontmatter has no description")
    else:
        if len(desc) > MAX_DESC:
            errors.append(f"description is {len(desc)} characters (max {MAX_DESC})")
        if "\n" in desc:
            errors.append("description must be one paragraph")
        if not WHEN_PHRASES.search(desc):
            warnings.append("description has no 'Use when...' phrase")

    plugins = listed.get(folder.name, [])
    if len(plugins) == 0:
        errors.append("not listed in any plugin in marketplace.json")
    elif len(plugins) > 1:
        errors.append(f"listed in more than one plugin: {', '.join(plugins)}")

    siblings = {p.name for p in SKILLS.iterdir() if p.is_dir()} - {folder.name}
    for f in sorted(folder.rglob("*")):
        if not f.is_file() or "__pycache__" in f.parts:
            continue
        rel = f.relative_to(folder)
        if f.name == ".env" or f.suffix == ".env":
            errors.append(f"{rel}: .env file")
        if f.stat().st_size > BIG_FILE:
            warnings.append(f"{rel}: larger than 10 MB")
        if f.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = f.read_text(errors="replace")
        for label, pat in SECRETS:
            if pat.search(text):
                errors.append(f"{rel}: possible {label}")
        for m in re.finditer(r"\.\./([A-Za-z0-9_.-]+)", text):
            if m.group(1) in siblings or m.group(1) in ("..", "_shared", "skills", "tools"):
                errors.append(f"{rel}: refers outside the skill folder ({m.group(0)})")
        for m in set(HOME_PATH.findall(text)):
            warnings.append(f"{rel}: home-directory path {m}")
    return name, errors, warnings


def validate(quiet=False):
    """Return True if every skill passes. Prints a table unless quiet."""
    listed, mp_errors = marketplace_paths()
    folders = sorted(p for p in SKILLS.iterdir() if p.is_dir()) if SKILLS.exists() else []
    rows, seen, ok = [], {}, not mp_errors

    for folder in folders:
        name, errors, warnings = check_skill(folder, listed)
        if name:
            if name in seen:
                errors.append(f"duplicate name, also used by {seen[name]}")
            seen[name] = folder.name
        ok &= not errors
        rows.append((folder.name, "FAIL" if errors else ("warn" if warnings else "pass"),
                     errors + [f"warning: {w}" for w in warnings]))

    for orphan in sorted(set(listed) - {f.name for f in folders}):
        mp_errors.append(f"marketplace.json lists skills/{orphan}, which does not exist")
        ok = False

    if not quiet:
        width = max([len(r[0]) for r in rows] + [5])
        print(f"{'skill':<{width}}  result  notes")
        print(f"{'-' * width}  ------  -----")
        for skill, status, notes in rows:
            print(f"{skill:<{width}}  {status:<6}  {notes[0] if notes else ''}")
            for n in notes[1:]:
                print(f"{'':<{width}}          {n}")
        for e in mp_errors:
            print(f"marketplace: {e}")
        print(f"\n{len(rows)} skill(s): {'all passed' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if validate() else 1)
