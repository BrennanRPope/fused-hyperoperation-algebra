"""Constant census across every result page, tested against a fair null.

Reads every *output*.txt and *results*.json under math-theory (not the Lean
library, not this folder), pulls out every number printed with at least 6
significant digits, and asks which named constants it matches (the constant,
its negative, its reciprocal, or one less / one more, since mult-add shifts
by 1). Then it builds a null: the same numbers with every digit after the
second significant one replaced at random, 2000 times, and counts matches the
same way. A constant is a real repeat only if it turns up well above its null.

Writes census_sig<MIN_SIG>.json (counts, pages, null mean/sd, p-values) and prints a table.
"""
import json, math, os, random, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
phi = (1 + 5 ** 0.5) / 2

CONSTANTS = {
    # golden family (the thing asked about)
    "golden ratio φ": phi,
    "φ² (Beraha B5)": phi * phi,
    # other quadratic irrationals of similar simplicity (controls for φ)
    "√2": 2 ** 0.5,
    "silver 1+√2": 1 + 2 ** 0.5,
    "√3": 3 ** 0.5,
    "√5": 5 ** 0.5,
    "2+√3": 2 + 3 ** 0.5,
    # cubic / cyclotomic neighbours
    "2cos(π/7)": 2 * math.cos(math.pi / 7),
    "plastic ρ": 1.324717957244746,
    # transcendental classics
    "π": math.pi,
    "π/2": math.pi / 2,
    "e": math.e,
    "1/e": 1 / math.e,
    "ln 2": math.log(2),
    "e^(1/e)": math.e ** (1 / math.e),
    "omega W(1)": 0.5671432904097838,
    "Euler γ": 0.5772156649015329,
    "Catalan G": 0.915965594177219,
    "ζ(3)": 1.2020569031595942,
    "ln 3": math.log(3),
    "2/π": 2 / math.pi,
}

TRANSFORMS = {
    "x": lambda c: c, "-x": lambda c: -c, "1/x": lambda c: 1 / c,
    "x-1": lambda c: c - 1, "x+1": lambda c: c + 1,
}
REL_TOL = float(os.environ.get("REL_TOL", "3e-6"))
MIN_SIG = int(os.environ.get("MIN_SIG", "6"))
NUM = re.compile(r"(?<![\w.])-?\d+\.\d+(?:[eE][-+]?\d+)?")


def sig_digits(s):
    m = s.lstrip("-").split("e")[0].split("E")[0].replace(".", "").lstrip("0")
    return len(m)


def gather():
    items = []  # (page, string)
    for dp, dn, fn in os.walk(ROOT):
        if "lean-library" in dp or os.sep + "patterns" in dp:
            continue
        for f in fn:
            if not (("output" in f and f.endswith(".txt")) or ("results" in f and f.endswith(".json"))):
                continue
            page = os.path.relpath(dp, ROOT).split(os.sep)[0]
            try:
                txt = open(os.path.join(dp, f), errors="ignore").read()
            except OSError:
                continue
            for s in NUM.findall(txt):
                if sig_digits(s) >= MIN_SIG:
                    items.append((page, s))
    return items


TARGETS = [(name, tname, tf(c)) for name, c in CONSTANTS.items() for tname, tf in TRANSFORMS.items()]


def matches(x):
    out = set()
    if x == 0 or not math.isfinite(x):
        return out
    for name, tname, t in TARGETS:
        if t != 0 and abs(x - t) <= REL_TOL * abs(t):
            out.add(name)
    return out


def scramble(s, rng):
    neg = s.startswith("-")
    body = s.lstrip("-")
    mant, _, exp = body.replace("E", "e").partition("e")
    chars = list(mant)
    seen = 0
    for i, ch in enumerate(chars):
        if ch.isdigit():
            if seen == 0 and ch == "0":
                continue
            seen += 1
            if seen > 2:
                chars[i] = str(rng.randrange(10))
    out = ("-" if neg else "") + "".join(chars) + (("e" + exp) if exp else "")
    return float(out)


def main(reps=2000):
    items = gather()
    # de-duplicate the same printed number on the same page (repeated lines)
    uniq = sorted(set(items))
    obs = {n: 0 for n in CONSTANTS}
    pages = {n: {} for n in CONSTANTS}
    for page, s in uniq:
        for n in matches(float(s)):
            obs[n] += 1
            pages[n][page] = pages[n].get(page, 0) + 1
    rng = random.Random(20261125)
    null = {n: [] for n in CONSTANTS}
    for _ in range(reps):
        cnt = {n: 0 for n in CONSTANTS}
        for page, s in uniq:
            for n in matches(scramble(s, rng)):
                cnt[n] += 1
        for n in CONSTANTS:
            null[n].append(cnt[n])
    res = []
    for n in CONSTANTS:
        arr = null[n]
        mu = sum(arr) / reps
        sd = (sum((a - mu) ** 2 for a in arr) / reps) ** 0.5
        p = (1 + sum(1 for a in arr if a >= obs[n])) / (reps + 1)
        res.append({"constant": n, "observed": obs[n], "null_mean": mu, "null_sd": sd,
                    "p": p, "pages": pages[n], "n_pages": len(pages[n])})
    res.sort(key=lambda r: (-r["n_pages"], -r["observed"]))
    out = {"numbers_scanned": len(uniq), "pages_scanned": len({p for p, _ in uniq}),
           "rel_tol": REL_TOL, "reps": reps, "results": res}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"census_sig{MIN_SIG}.json"), "w"), indent=1)
    print(f"{len(uniq)} distinct numbers (>={MIN_SIG} sig. digits) from {out['pages_scanned']} pages")
    print(f"{'constant':18s} {'obs':>4s} {'pages':>5s} {'null':>6s} {'p':>7s}  pages")
    for r in res:
        print(f"{r['constant']:18s} {r['observed']:4d} {r['n_pages']:5d} {r['null_mean']:6.2f} {r['p']:7.4f}  {dict(sorted(r['pages'].items()))}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2000)
