#!/usr/bin/env python3
"""Compile every topic Lean file that imports Mathlib, against the pinned Mathlib (report only).

tools/topic_lean.py skips files that import Mathlib, because core Lean can't build them. This script
builds the ones it can, one file at a time, inside the Lake environment of topics/mathlib-map (same
Lean 4.37.0-rc1 and Mathlib v4.37.0-rc1 as lean-library-mathlib, plus the core library FusedAlgebra).

A file is built here when every import is Mathlib.* or FusedAlgebra.*. Files that import their own
sibling modules (multi-file drafts such as proof-complexity or flt-tracking) are listed as skipped:
they need their own Lake project. topics/mathlib-map itself is built by its own CI steps.

Run from the repository root after `lake exe cache get` in topics/mathlib-map:
    python3 tools/mathlib_topic_files.py [files...]
Prints OK, FAIL (with the compiler's errors) or SKIP per file and a summary. Exit code is 1 if any
file failed, so CI can show the step as failed; the CI step is continue-on-error and does not block.
"""
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LAKE_DIR = ROOT / "topics" / "mathlib-map"
TIMEOUT = 1800


def imports(path):
    return re.findall(r"^\s*import\s+(\S+)", path.read_text(errors="replace"), re.M)


def candidates(args):
    if args:
        return [Path(a).resolve() for a in args]
    out = []
    for f in sorted((ROOT / "topics").rglob("*.lean")):
        if LAKE_DIR in f.parents or ".lake" in f.parts:
            continue
        if any(m.split(".")[0] == "Mathlib" for m in imports(f)):
            out.append(f)
    return out


def main():
    files = candidates(sys.argv[1:])
    ok, bad, skipped = [], [], []
    for f in files:
        rel = f.relative_to(ROOT)
        other = [m for m in imports(f) if m.split(".")[0] not in ("Mathlib", "FusedAlgebra")]
        if other:
            skipped.append(rel)
            print(f"SKIP   {rel}: imports {', '.join(other)} (needs its own Lake project)", flush=True)
            continue
        t0 = time.time()
        try:
            r = subprocess.run(["lake", "env", "lean", str(f)], cwd=LAKE_DIR,
                               capture_output=True, text=True, timeout=TIMEOUT)
            out = (r.stdout + r.stderr).strip()
            failed = r.returncode != 0 or "declaration uses 'sorry'" in out
        except subprocess.TimeoutExpired:
            out, failed = f"timed out after {TIMEOUT} s", True
        dt = time.time() - t0
        if failed:
            bad.append(rel)
            print(f"FAIL   {rel}  ({dt:.0f}s)", flush=True)
            for line in out.splitlines():
                print(f"       {line}")
        else:
            ok.append(rel)
            note = " (warnings)" if "warning" in out else ""
            print(f"OK     {rel}  ({dt:.0f}s){note}", flush=True)
    print()
    print(f"{len(ok)} of {len(ok) + len(bad)} Mathlib topic files build; {len(bad)} fail; "
          f"{len(skipped)} skipped (sibling imports).")
    if bad:
        print("Failing: " + ", ".join(str(b) for b in bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
