"""Corner-weight signs of every 5-region holographic entropy inequality (DRAFT, Nov 26 2026).
Source of the 8 orbit representatives: S. Hernandez-Cuenca, 'The holographic entropy cone for five regions',
Phys. Rev. D 100, 026004 (2019), arXiv:1903.09148, Table 1 (orbit lengths 15, 20, 45, 72, 10, 60, 60, 90).
Each instance is written in the 5 regions A..E (the purifier O removed by S(X) = S(complement)), then in
corner weights y = Mobius inverse of S, so S(I) = sum over corners J inside I of y_J."""
import itertools, json, re
from collections import Counter
SYM = "ABCDEO"; ALL = frozenset(SYM)
REPS = {
 1: ("SA + SB", "SAB"),
 2: ("SAB + SAC + SBC", "SA + SB + SC + SABC"),
 3: ("SABC + SADE + SBCDE", "SA + SBC + SDE + SABCDE"),
 4: ("SABC + SABD + SACE + SBDE + SCDE", "SAB + SAC + SBD + SCE + SDE + SABCDE"),
 5: ("SABC + SABD + SABE + SACD + SACE + SADE + SBCE + SBDE + SCDE", "SAB + SAC + SAD + SBE + SCE + SDE + SBCD + SABCE + SABDE + SACDE"),
 6: ("3SABC + 3SABD + SABE + SACD + 3SACE + SADE + SBCD + SBCE + SBDE + SCDE", "2SAB + 2SAC + SAD + SAE + SBC + 2SBD + 2SCE + SDE + 2SABCD + 2SABCE + SABDE + SACDE"),
 7: ("2SABC + SABD + SABE + SACD + SADE + SBCE + SBDE", "SAB + SAC + SAD + SBC + SBE + SDE + SABCD + SABCE + SABDE"),
 8: ("SAD + SBC + SABE + SACE + SADE + SBDE + SCDE", "SA + SB + SC + SD + SAE + SDE + SBCE + SABDE + SACDE"),
}
ORBIT = {1: 15, 2: 20, 3: 45, 4: 72, 5: 10, 6: 60, 7: 60, 8: 90}
NAME = {1: "SA", 2: "MMI", 3: "MMI lifted", 4: "cyclic", 5: "new (10)", 6: "new (60a)", 7: "new (60b)", 8: "new (90)"}
def parse(side, sign):
    out = Counter()
    for term in side.split("+"):
        m = re.fullmatch(r"\s*(\d*)S([A-E]+)\s*", term)
        out[frozenset(m.group(2))] += sign*(int(m.group(1)) if m.group(1) else 1)
    return out
def canon(ineq):
    out = Counter()
    for I, c in ineq.items():
        if "O" in I: I = ALL - I
        if I and I != ALL - {"O"} | set() or I:
            out[I] += c
    return frozenset((I, c) for I, c in out.items() if c and I)
def corner_coeffs(inst):
    regs = "ABCDE"; c = {}
    for r in range(1, 6):
        for J in itertools.combinations(regs, r):
            v = sum(cI for I, cI in inst if set(J) <= I)
            if v: c["".join(J)] = v
    return c
def minkowski_coeffs(inst):
    regs = "ABCDE"; c = {}
    for r in range(1, 6):
        for J in itertools.combinations(regs, r):
            v = sum(cI for I, cI in inst if set(J) & I)
            if v: c["".join(J)] = v
    return c
report = {}; allok = True
for k, (L, R) in REPS.items():
    base = parse(L, 1) + Counter(); base.update(parse(R, -1))
    insts = set()
    for p in itertools.permutations(SYM):
        mp = dict(zip(SYM, p))
        insts.add(canon(Counter({frozenset(mp[x] for x in I): c for I, c in base.items()})))
    stats = Counter(); examples = {}
    for inst in insts:
        cc = corner_coeffs(inst); mk = minkowski_coeffs(inst)
        grades = sorted({len(J) for J in cc}); signs = {v > 0 for v in cc.values()}
        key = ("one sign" if len(signs) == 1 else "mixed signs") + ", lowest corner grade " + str(grades[0])
        stats[key] += 1; examples.setdefault(key, cc)
        stats["minkowski weights equal corner weights"] += (mk == cc)
    report[k] = {"name": NAME[k], "instances": len(insts), "paper_orbit": ORBIT[k], "stats": dict(stats),
                 "example": examples}
    ok = len(insts) == ORBIT[k]; allok &= ok
    print(f"[{k}] {NAME[k]:11s} instances {len(insts):3d} (paper {ORBIT[k]:3d}) {'OK' if ok else 'MISMATCH'}")
    for key, v in stats.items(): print(f"      {key}: {v}")
    for key, ex in examples.items(): print(f"      e.g. ({key}): {ex}")
tot = sum(r["instances"] for r in report.values())
print(f"\ntotal facets {tot} (paper 372); orbit sizes all match: {allok}")
json.dump(report, open("hec5_signs.json", "w"), indent=1, default=str)
