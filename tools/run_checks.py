#!/usr/bin/env python3
"""Run every topic's numeric checks and compare them with the saved output.

For each topics/**/checks*.py this runs the script in its own folder and fails if:
  - the script exits with an error or times out,
  - any output line reports a failure (a line starting with FAIL, or "N failed" with N > 0),
  - the saved checks*_output.txt states a pass count ("29/29", "18 of 18", "49 passed",
    "101 PASS") and the new run does not reach the same count.
Scripts with no saved output still have to run cleanly with no FAIL line.
Numeric checks are evidence, not proofs; the Lean library is where proofs live.

Usage: python3 tools/run_checks.py [--jobs N] [--timeout SECONDS] [paths...]
"""
import argparse, concurrent.futures as cf, os, re, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COUNT_PATTERNS = [
    re.compile(r"(\d+)\s*/\s*(\d+)\s+(?:checks?\s+)?(?:pass|passed|ok)", re.I),
    re.compile(r"(\d+)\s+of\s+(\d+)\s+(?:checks?\s+)?(?:pass|passed)", re.I),
    re.compile(r"(\d+)\s*/\s*(\d+)\s*$"),
    re.compile(r"(?:all\s+)?(\d+)\s+(?:checks?\s+)?(?:pass|passed|PASS)\b", re.I),
]
FAIL_LINE = re.compile(r"^\s*(FAIL|FAILED)\b")
FAILED_COUNT = re.compile(r"\b([1-9]\d*)\s+(?:checks?\s+)?(?:fail|failed|FAIL)\b", re.I)


def pass_count(text):
    """The last pass count stated in the text, as an int, or None."""
    best = None
    for line in text.splitlines():
        for pat in COUNT_PATTERNS:
            m = pat.search(line)
            if m:
                best = int(m.group(1))
                break
    return best


def run_one(script, timeout):
    expected_file = script.with_name(script.stem + "_output.txt")
    expected = expected_file.read_text(errors="replace") if expected_file.exists() else None
    t0 = time.time()
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONHASHSEED="0")
    try:
        p = subprocess.run([sys.executable, script.name], cwd=script.parent, capture_output=True,
                           text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return script, False, f"timed out after {timeout}s", time.time() - t0
    out = p.stdout + p.stderr
    dt = time.time() - t0
    problems = []
    if p.returncode != 0:
        tail = " | ".join(out.strip().splitlines()[-3:])
        problems.append(f"exit code {p.returncode}: {tail[:300]}")
    fails = [l.strip() for l in p.stdout.splitlines() if FAIL_LINE.match(l) or FAILED_COUNT.search(l)]
    if fails:
        problems.append(f"{len(fails)} failing line(s), first: {fails[0][:200]}")
    if expected is not None:
        want, got = pass_count(expected), pass_count(p.stdout)
        if want is not None and (got is None or got < want):
            problems.append(f"pass count {got} is below the saved {want}")
    note = "; ".join(problems) if problems else (
        f"ok ({pass_count(p.stdout)} pass)" if pass_count(p.stdout) is not None else "ok")
    return script, not problems, note, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("paths", nargs="*")
    a = ap.parse_args()
    scripts = [Path(p).resolve() for p in a.paths] if a.paths else sorted((ROOT / "topics").rglob("checks*.py"))
    results = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for script, ok, note, dt in ex.map(lambda s: run_one(s, a.timeout), scripts):
            rel = script.relative_to(ROOT)
            print(f"{'PASS' if ok else 'FAIL'}  {rel}  ({dt:.0f}s)  {note}", flush=True)
            results.append(ok)
    bad = results.count(False)
    print(f"\n{len(results) - bad} of {len(results)} check scripts pass.")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
