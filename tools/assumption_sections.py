#!/usr/bin/env python3
"""List the topic write-ups (topics/**/README.md) that have no "Assumptions" heading yet.

Convention (CLAUDE.md): each topic write-up says what its results rest on, under a heading
that contains the word "Assumptions": the hypotheses its theorems carry, any modeling
assumptions (which data, which equation describes the system), and which claims are only
checked numerically. This script only reports; it never fails the build.

Usage: python3 tools/assumption_sections.py
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HEAD = re.compile(r"^#{1,6}\s.*assumption", re.I | re.M)

files = sorted(p for p in (ROOT / "topics").rglob("README.md") if "lean" not in p.parent.parts)
missing = [p for p in files if not HEAD.search(p.read_text(errors="replace"))]
print(f"{len(files) - len(missing)} of {len(files)} topic write-ups have an Assumptions section.")
for p in missing:
    print(f"  missing: {p.parent.relative_to(ROOT)}")
