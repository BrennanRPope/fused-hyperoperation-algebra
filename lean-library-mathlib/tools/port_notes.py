#!/usr/bin/env python3
"""Check the name-mapping tables in port-notes/ and print a per-module summary.

Usage:  python3 tools/port_notes.py            # coverage report
        python3 tools/port_notes.py --table    # Markdown summary table for PORT-NOTES.md

Every theorem listed in ../lean-library/theorems.tsv must appear in the table of its module.
"""
import collections, glob, os, sys

here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
orig = collections.defaultdict(list)
for line in open(os.path.join(here, "..", "lean-library", "theorems.tsv"), encoding="utf-8"):
    mod, name = line.rstrip("\n").split("\t")[:2]
    orig[mod].append(name)

STATUSES = ("ported", "mathlib", "dropped", "not-ported")
rows, problems = {}, []
for mod in sorted(orig):
    path = os.path.join(here, "port-notes", mod + ".tsv")
    if not os.path.exists(path):
        problems.append(f"{mod}: no table")
        continue
    table = []
    for i, line in enumerate(open(path, encoding="utf-8")):
        if i == 0:
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3 or parts[1] not in STATUSES:
            problems.append(f"{mod}: bad row {i + 1}: {line.strip()[:60]}")
            continue
        table.append(parts + [""] * (4 - len(parts)))
    rows[mod] = table
    have = {r[0] for r in table}
    # theorems.tsv gives fully qualified names; tables use the name without its namespace
    for full in orig[mod]:
        if not any(full == h or full.endswith("." + h) for h in have):
            problems.append(f"{mod}: {full} missing from table")

if "--table" in sys.argv:
    print("| Module | Original theorems | Rows | Ported | Mathlib | Dropped | Not ported |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    tot = collections.Counter()
    for mod in sorted(orig):
        t = rows.get(mod, [])
        c = collections.Counter(r[1] for r in t)
        tot.update(c); tot["orig"] += len(orig[mod]); tot["rows"] += len(t)
        print(f"| `{mod.removeprefix('FusedAlgebra.')}` | {len(orig[mod])} | {len(t)} | "
              + " | ".join(str(c[s]) for s in STATUSES) + " |")
    print(f"| **Total** | {tot['orig']} | {tot['rows']} | "
          + " | ".join(str(tot[s]) for s in STATUSES) + " |")
else:
    print(f"{len(rows)}/{len(orig)} modules have tables")
    for p in problems:
        print("PROBLEM", p)
    sys.exit(1 if problems else 0)
