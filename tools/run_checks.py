#!/usr/bin/env python3
"""Run every topic's numeric checks and compare them with the saved output.

For each topics/**/checks*.py, and each script listed in tools/extra_checks.txt, this runs the
script in its own folder and fails if:
  - the script exits with an error or times out,
  - any output line reports a failure (a line starting with FAIL, or "N failed" with N > 0),
  - the saved checks*_output.txt states a pass count ("29/29", "18 of 18", "49 passed",
    "101 PASS") and the new run does not reach the same count.
Scripts with no saved output still have to run cleanly with no FAIL line.
Some scripts name their working folder as /mnt/project-files/math-theory/<topic> (to import a
helper such as topics/gr-bridges/ga.py, or to write their output). Those paths are rewritten to
the repo's topics/<topic> before the script runs, so the repo copy is what gets checked.
Numeric checks are evidence, not proofs; the Lean library is where proofs live.
A script listed in tools/checks_needing_packages.txt with a Python module CI does not install is
reported as SKIP (not a failure) while that module can't be imported; it runs as soon as it can.
A check line listed in tools/checks_known_stale.txt (it tests a snapshot of a file that a later round
changed, such as an older theorem count in a topic's Lean file) is reported as STALE, not counted as a
failure, and lowers the saved pass count by one. Every other line still has to pass.

Usage: python3 tools/run_checks.py [--jobs N] [--timeout SECONDS] [paths...]
"""
import argparse, concurrent.futures as cf, importlib.util, os, re, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COUNT_PATTERNS = [
    re.compile(r"(\d+)\s*/\s*(\d+)\s+(?:checks?\s+)?(?:pass|passed|ok)", re.I),
    re.compile(r"(\d+)\s+of\s+(\d+)\s+(?:checks?\s+)?(?:pass|passed)", re.I),
    re.compile(r"(\d+)\s*/\s*(\d+)\s*$"),
    re.compile(r"(?:all\s+)?(\d+)\s+(?:checks?\s+)?(?:pass|passed|PASS)\b", re.I),
]
FAIL_LINE = re.compile(r"^\s*(FAIL|FAILED)\b")
# A line that itself reports a pass ("PASS ... 12 fail Turán") is not a failure count.
PASS_LINE = re.compile(r"^\s*(PASS|OK)\b")
WORK_PREFIX = "/mnt/project-files/math-theory/"
RUNNER = (
    "import sys\n"
    "path, prefix, topics = sys.argv[1], sys.argv[2], sys.argv[3]\n"
    "src = open(path, encoding='utf-8').read().replace(prefix, topics)\n"
    "sys.argv = [path]\n"
    "exec(compile(src, path, 'exec'), {'__name__': '__main__', '__file__': path})\n"
)
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


def run_one(script, timeout, stale_ids=()):
    expected_file = script.with_name(script.stem + "_output.txt")
    expected = expected_file.read_text(errors="replace") if expected_file.exists() else None
    t0 = time.time()
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONHASHSEED="0")
    # Checks run side by side (--jobs), so each one gets a single math-library
    # thread; otherwise numpy/scipy threads oversubscribe the cores and some
    # checks slow from seconds to past the timeout. Results are unchanged.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env.setdefault(var, "1")
    try:
        cmd = [sys.executable, script.name]
        if WORK_PREFIX in script.read_text(errors="replace"):
            cmd = [sys.executable, "-c", RUNNER, str(script), WORK_PREFIX, str(ROOT / "topics") + "/"]
        p = subprocess.run(cmd, cwd=script.parent, capture_output=True,
                           text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return script, False, f"timed out after {timeout}s", time.time() - t0
    out = p.stdout + p.stderr
    dt = time.time() - t0
    problems = []
    if p.returncode != 0:
        tail = " | ".join(out.strip().splitlines()[-3:])
        problems.append(f"exit code {p.returncode}: {tail[:300]}")
    fails = [l.strip() for l in p.stdout.splitlines()
             if FAIL_LINE.match(l) or (FAILED_COUNT.search(l) and not PASS_LINE.match(l))]
    stale = [l for l in fails if any(re.match(rf"(FAIL|FAILED)\s+{re.escape(i)}\b", l) for i in stale_ids)]
    fails = [l for l in fails if l not in stale]
    if fails:
        problems.append(f"{len(fails)} failing line(s), first: {fails[0][:200]}")
    if expected is not None:
        want, got = pass_count(expected), pass_count(p.stdout)
        if want is not None:
            want -= len(stale)
        if want is not None and (got is None or got < want):
            problems.append(f"pass count {got} is below the saved {want}")
    note = "; ".join(problems) if problems else (
        f"ok ({pass_count(p.stdout)} pass)" if pass_count(p.stdout) is not None else "ok")
    if stale:
        note += f"; STALE (known, see tools/checks_known_stale.txt): {'; '.join(l[:120] for l in stale)}"
    return script, not problems, note, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("paths", nargs="*")
    a = ap.parse_args()
    if a.paths:
        scripts = [Path(p).resolve() for p in a.paths]
    else:
        scripts = sorted((ROOT / "topics").rglob("checks*.py"))
        extra = ROOT / "tools" / "extra_checks.txt"
        if extra.exists():
            for line in extra.read_text().splitlines():
                if line.strip() and not line.startswith("#"):
                    p = ROOT / line.strip()
                    if not p.exists():
                        sys.exit(f"tools/extra_checks.txt lists {line.strip()}, which does not exist")
                    if p not in scripts:
                        scripts.append(p)
    needs = {}
    needs_file = ROOT / "tools" / "checks_needing_packages.txt"
    if needs_file.exists():
        for line in needs_file.read_text().splitlines():
            parts = line.split()
            if parts and not line.startswith("#"):
                needs[(ROOT / parts[0]).resolve()] = parts[1]
    kept = []
    for s in scripts:
        mod = needs.get(s.resolve())
        if mod and importlib.util.find_spec(mod) is None:
            print(f"SKIP  {s.relative_to(ROOT)}  needs the Python module {mod}, which is not installed", flush=True)
        else:
            kept.append(s)
    scripts = kept
    stale_ids = {}
    stale_file = ROOT / "tools" / "checks_known_stale.txt"
    if stale_file.exists():
        for line in stale_file.read_text().splitlines():
            parts = line.split()
            if parts and not line.startswith("#"):
                stale_ids.setdefault((ROOT / parts[0]).resolve(), []).append(parts[1])
    # Scripts that need most of the runner's memory (tools/checks_run_alone.txt) run after the
    # pool, one at a time, so they never share the machine with other checks.
    alone = set()
    alone_file = ROOT / "tools" / "checks_run_alone.txt"
    if alone_file.exists():
        for line in alone_file.read_text().splitlines():
            parts = line.split()
            if parts and not line.startswith("#"):
                alone.add((ROOT / parts[0]).resolve())
    pooled = [s for s in scripts if s.resolve() not in alone]
    later = [s for s in scripts if s.resolve() in alone]
    results = []

    def report(script, ok, note, dt):
        rel = script.relative_to(ROOT)
        print(f"{'PASS' if ok else 'FAIL'}  {rel}  ({dt:.0f}s)  {note}", flush=True)
        results.append(ok)

    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(lambda s: run_one(s, a.timeout, stale_ids.get(s.resolve(), ())), pooled):
            report(*r)
    for s in later:
        report(*run_one(s, a.timeout, stale_ids.get(s.resolve(), ())))
    bad = results.count(False)
    print(f"\n{len(results) - bad} of {len(results)} check scripts pass.")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
