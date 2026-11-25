#!/usr/bin/env python3
"""Compile every Lean file under topics/ on its own and report which ones still build.

Most topic folders keep the Lean file a result was first proved in. The shared library in
lean-library/ is where the proofs live and what CI checks for sorry and axioms; these copies are
the historical record. This script makes sure each one is still checked: it compiles every
topics/**/*.lean with core Lean, with the library's build/ on LEAN_PATH (so `import FusedAlgebra...`
works) and with sibling modules compiled first (so `import Fermat3` or `import Axioms.Rederive`
works). Files that need Mathlib are skipped here: they are built in lean-library-mathlib/.

A file that fails is an error unless tools/topic_lean_known.tsv lists it with a reason (a copy
that broke when the library twin was repaired, a sketch, and so on). A listed file that now builds
is reported so the list can be trimmed.

Run after `bash lean-library/build.sh`.
Usage: python3 tools/topic_lean.py [--jobs N]
"""
import argparse, concurrent.futures as cf, os, re, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIB_BUILD = ROOT / "lean-library" / "build"
KNOWN = ROOT / "tools" / "topic_lean_known.tsv"
IMPORT = re.compile(r"^\s*(?:public\s+)?import\s+([A-Za-z0-9_.]+)", re.M)


def find_local(f, mod):
    """The file and source root for a module that lives beside f (not in the library)."""
    rel = Path(*mod.split(".")).with_suffix(".lean")
    d = f.parent
    for _ in range(4):
        cand = d / rel
        if cand.exists():
            return cand, d
        if d == ROOT / "topics":
            break
        d = d.parent
    return None, None


def compile_file(f, timeout):
    """Compile f, compiling its local imports first into a scratch build dir. Returns (ok, msg)."""
    with tempfile.TemporaryDirectory() as tmp:
        lean_path = f"{LIB_BUILD}{os.pathsep}{tmp}/_mods"
        env = dict(os.environ, LEAN_PATH=lean_path)
        done = set()
        work = Path(tmp) / "_work"
        (work / "build").mkdir(parents=True)

        def build(path, root, as_module):
            text = path.read_text(errors="replace")
            for mod in IMPORT.findall(text):
                if mod.split(".")[0] in ("FusedAlgebra", "Lean", "Init", "Std", "Lake"):
                    continue
                if mod.split(".")[0] == "Mathlib":
                    return None, "needs Mathlib"
                dep, droot = find_local(path, mod)
                if dep is None:
                    return False, f"import {mod} not found"
                if dep in done:
                    continue
                ok, msg = build(dep, droot, True)
                if ok is None:
                    return None, msg
                if not ok:
                    return False, f"{dep.relative_to(ROOT)}: {msg}"
                done.add(dep)
            args = ["lean"]
            if as_module:
                rel = path.relative_to(root).with_suffix(".olean")
                out = Path(tmp) / "_mods" / rel
                out.parent.mkdir(parents=True, exist_ok=True)
                args += [f"--root={root}", "-o", str(out)]
            args.append(str(path))
            try:
                # Some files write their theorem lists to build/ when compiled; keep that in scratch.
                p = subprocess.run(args, cwd=work, env=env, capture_output=True,
                                   text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                return False, f"timed out after {timeout}s"
            if p.returncode != 0:
                err = next((l for l in (p.stdout + p.stderr).splitlines() if "error" in l),
                           (p.stdout + p.stderr).strip()[:200])
                return False, err.replace(str(ROOT) + "/", "").strip()[:200]
            return True, ""

        return build(f, f.parent, False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()
    if not (LIB_BUILD / "FusedAlgebra.olean").exists():
        sys.exit("lean-library/build is missing: run bash lean-library/build.sh first")
    known = {}
    if KNOWN.exists():
        for line in KNOWN.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                path, _, reason = line.partition("\t")
                known[path.strip()] = reason.strip()
    files, skipped = [], []
    for f in sorted((ROOT / "topics").rglob("*.lean")):
        if ".lake" in f.parts:
            continue
        if re.search(r"^\s*import\s+Mathlib", f.read_text(errors="replace"), re.M):
            skipped.append(f)
        else:
            files.append(f)
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        results = list(ex.map(lambda f: (f, *compile_file(f, args.timeout)), files))
    bad, fixed, ok = [], [], 0
    for f, good, msg in list(results):
        rel = str(f.relative_to(ROOT))
        if good is None:
            skipped.append(f)
            results.remove((f, good, msg))
            continue
        if good:
            ok += 1
            if rel in known:
                fixed.append(rel)
        elif rel in known:
            print(f"KNOWN  {rel}: {known[rel]}")
        else:
            bad.append(rel)
            print(f"FAIL   {rel}: {msg}")
    for rel in fixed:
        print(f"NOW BUILDS (remove from {KNOWN.name})  {rel}")
    for f in skipped:
        print(f"SKIP   {f.relative_to(ROOT)}: needs Mathlib")
    listed_missing = [k for k in known if not (ROOT / k).exists()]
    for k in listed_missing:
        print(f"STALE  {k} is listed in {KNOWN.name} but does not exist")
    print(f"\n{ok} of {len(results)} topic Lean files build; {len(results) - ok - len(bad)} known "
          f"failures; {len(bad)} new failures; {len(skipped)} skipped (Mathlib).")
    sys.exit(1 if bad or listed_missing else 0)


if __name__ == "__main__":
    main()
