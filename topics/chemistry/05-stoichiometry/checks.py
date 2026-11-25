"""Stoichiometry, definite and multiple proportions through the circle operation.
1. Yields multiply, so losses the circle operation at k = −1; a synthesis tree's yield is set by its depth.
2. Isotope patterns: the lift (p_light + p_heavy·x)^n; "at least one heavy atom" is the k = −1 power.
3. Multiple proportions: small whole-number ratios; mixed oxides are mediants of their parents.
4. Balancing: an equation is a null vector of the element matrix (Haber, respiration, photosynthesis).
Run: python3 checks.py"""
import json, math
from functools import reduce
import numpy as np
import sympy as sp
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b

# 1. yields
ys = [0.92, 0.85, 0.78, 0.95, 0.88, 0.90, 0.81]
loss = reduce(lambda a, b: w(-1, a, b), [1-y for y in ys])
check("seven steps: overall loss = fold of step losses under k = −1 = 1 − Π yields", abs(loss - (1 - np.prod(ys))) < 1e-15,
      f"overall yield {100*np.prod(ys):.1f}%")
import itertools
orders = {round(reduce(lambda a, b: w(-1, a, b), [1-ys[i] for i in p]), 14) for p in itertools.permutations(range(5))}
check("the k = −1 fold is order-free (any order of the same steps gives one loss)", len(orders) == 1)
lin = 0.8**7; conv = 0.8**3
R.update(linear=lin, convergent=conv)
check("eight fragments at 80% per coupling: linear 0.8⁷ = 21.0%, convergent (three levels) 0.8³ = 51.2%", abs(lin-0.2097) < 1e-4 and abs(conv-0.512) < 1e-12)
ten = 0.9**10
check("ten steps at 90% give 34.9%: k = −1 power wpow(−1, 0.1, 10) = 1 − 0.9¹⁰ = 65.1% lost", abs((1 - 0.9**10) - 0.6513) < 1e-4)

# 2. isotope patterns (IUPAC representative abundances)
ab = {"C": [0.9893, 0.0107], "H": [0.999885, 0.000115], "Cl": [0.7576, 0.0, 0.2424], "Br": [0.5069, 0.0, 0.4931],
      "N": [0.99636, 0.00364], "O": [0.99757, 0.00038, 0.00205]}
def pattern(formula):
    p = np.array([1.0])
    for el, n in formula.items():
        for _ in range(n): p = np.convolve(p, ab[el])
    return p
cl2 = pattern({"Cl": 2}); br2 = pattern({"Br": 2})
r_cl = cl2[[0, 2, 4]]/cl2[0]; r_br = br2[[0, 2, 4]]/br2[2]
R.update(Cl2=list(cl2[[0, 2, 4]]), Br2=list(br2[[0, 2, 4]]))
check("Cl₂ peaks M : M+2 : M+4 = 1 : 0.640 : 0.102 (the 9 : 6 : 1 of (3 + x)²)", np.allclose(r_cl, [1, 2*0.2424/0.7576, (0.2424/0.7576)**2]),
      f"{r_cl.round(3)}")
check("Br₂ ≈ 1 : 2 : 1 ((1 + x)²): measured-abundance ratio 0.514 : 1 : 0.486", np.allclose(br2[[0, 2, 4]]/br2[2], [0.5069/(2*0.4931), 1, 0.4931/(2*0.5069)]),
      f"{r_br.round(3)}")
for n in (6, 20, 60):
    p = pattern({"C": n})
    assert abs(p[1]/p[0] - n*0.0107/0.9893) < 1e-12
    assert abs((1 - p[0]) - (1 - (1-0.0107)**n)) < 1e-12
c60 = pattern({"C": 60})
R["C60"] = dict(M1=c60[1]/c60[0], any13=1-c60[0])
check("Cₙ: M+1/M = n·p₁₃/p₁₂ exactly; the chance of at least one ¹³C is the k = −1 power 1 − (1 − p)ⁿ",
      True, f"C60: M+1 is {100*c60[1]/c60[0]:.1f}% of M; {100*(1-c60[0]):.1f}% of molecules carry ¹³C")
# generating function: the pattern of a compound is the product of its atoms' polynomials (lift)
a = pattern({"C": 2, "H": 6}); b = pattern({"O": 1}); c = pattern({"C": 2, "H": 6, "O": 1})
check("patterns multiply: ethanol's = ethane's × oxygen's (generating functions are the lift)", np.allclose(np.convolve(a, b), c))

# 3. multiple and definite proportions
mN, mO, mFe, mPb, mMn = 14.007, 15.999, 55.845, 207.2, 54.938
oxides = {"N₂O": (2, 1), "NO": (1, 1), "N₂O₃": (2, 3), "NO₂": (1, 2), "N₂O₅": (2, 5)}
perg = {k: o*mO/(n*mN) for k, (n, o) in oxides.items()}
ratios = [v/perg["N₂O"] for v in perg.values()]
R["NOx_per_gN"] = perg
check("Dalton: grams of O per gram of N in N₂O, NO, N₂O₃, NO₂, N₂O₅ stand 1 : 2 : 3 : 4 : 5",
      np.allclose(ratios, [1, 2, 3, 4, 5]), ", ".join(f"{k} {v:.3f}" for k, v in perg.items()))
def frac(nm, mm, no):  # mass fraction of metal
    return nm*mm/(nm*mm + no*mO)
mixes = [("Fe₃O₄ = FeO + Fe₂O₃", mFe, (1, 1), (2, 3), (3, 4)), ("Pb₃O₄ = 2 PbO + PbO₂", mPb, (2, 2), (1, 2), (3, 4)),
         ("Mn₃O₄ = MnO + Mn₂O₃", mMn, (1, 1), (2, 3), (3, 4))]
good = True; info = []
for name, mm, A, B, C in mixes:
    assert (A[0]+B[0], A[1]+B[1]) == C           # atom counts add: the formula is a mediant
    fA, fB, fC = frac(A[0], mm, A[1]), frac(B[0], mm, B[1]), frac(C[0], mm, C[1])
    good &= min(fA, fB) < fC < max(fA, fB)
    info.append(f"{name}: metal {100*fC:.2f}% between {100*fB:.2f}% and {100*fA:.2f}%")
check("mixed-valence oxides: metal:oxygen 3/4 is the mediant of the parents' ratios, and the mass fraction lies between", good, "; ".join(info))
R["mediants"] = info

# 4. balancing as a null vector
def balance(reac, prod):
    species = reac + prod
    els = sorted({e for s in species for e in s})
    M = sp.Matrix([[(s.get(e, 0) if i < len(reac) else -s.get(e, 0)) for i, s in enumerate(species)] for e in els])
    ns = M.nullspace()
    assert len(ns) == 1
    v = ns[0]; v = v * sp.ilcm(*[x.q for x in v]); v = v / sp.igcd(*[int(x) for x in v])
    return [abs(int(x)) for x in v]
eqs = {
 "Haber: N₂ + 3 H₂ → 2 NH₃": ([{"N": 2}, {"H": 2}], [{"N": 1, "H": 3}], [1, 3, 2]),
 "respiration: C₆H₁₂O₆ + 6 O₂ → 6 CO₂ + 6 H₂O": ([{"C": 6, "H": 12, "O": 6}, {"O": 2}], [{"C": 1, "O": 2}, {"H": 2, "O": 1}], [1, 6, 6, 6]),
 "photosynthesis: 6 CO₂ + 6 H₂O → C₆H₁₂O₆ + 6 O₂": ([{"C": 1, "O": 2}, {"H": 2, "O": 1}], [{"C": 6, "H": 12, "O": 6}, {"O": 2}], [6, 6, 1, 6]),
 "nitrogenase: N₂ + 8 H⁺… as N₂ + 8 H → 2 NH₃ + H₂": ([{"N": 2}, {"H": 1}], [{"N": 1, "H": 3}, {"H": 2}], None),
}
for name, (re_, pr, want) in eqs.items():
    if want is None: continue
    got = balance(re_, pr)
    check(f"balanced as the null vector of the element matrix — {name}", got == want, str(got))
# nitrogenase: N2 + 8H -> 2NH3 + H2 has a 2-dim family? check: element matrix 2 rows, 4 species -> nullity 2
M = sp.Matrix([[2, 0, -1, 0], [0, 1, -3, -2]])
check("N₂ + H → NH₃ + H₂ has a two-parameter family: stoichiometry alone cannot fix nitrogenase's one H₂ per N₂", len(M.nullspace()) == 2)

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
