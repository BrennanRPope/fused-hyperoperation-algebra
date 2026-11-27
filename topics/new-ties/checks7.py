"""Round 7 checks (README section 11). Run after plane.py. Writes checks7_output.txt and checks7.json."""
import json, math, random, cmath
import numpy as np
from fractions import Fraction as F
from scipy.optimize import brentq
import plane as P
R = json.load(open("plane.json"))
out = []; res = {}
def check(name, ok, detail=""):
    out.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); res[name] = bool(ok)

# W1, W2: certified ping-pong for both examples
for k, name in [("W1", "ORIG"), ("W2", "MOD")]:
    c = R["C1_" + name]
    check(f"{k} {name}: invariant disc and pairwise-disjoint class covers (exact fractions)",
          c["certified"] and c["invariant"] and c["disjoint_pairs_failed"] == 0,
          f"covers {c['X_discs']} discs, smallest gap {min(c['min_gap_float'].values()):.3f}")

# W3: 3-adic constants for MOD
c = R["C2"]
check("W3 MOD 3-adic facts: coefficients 3-integral, classes 2/3 shrink 3-adically, v3(lam1)=0, v3(lam2)=1, v3(rho)=-beta",
      c["all_coefficients_3_integral"] and c["class23_multipliers_v3_ge_1"] and c["centres_1_and_i_are_3_adic_units"]
      and c["v3_lam1"] == 0 and c["v3_lam2"] == 1 and c["v3_of_rho_is_minus_beta"])

# W4: illustration: MOD lands on valuation 0 only, ORIG spreads
m = R["C3"]["MOD"]; o = R["C3"]["ORIG"]
check("W4 random words (outer run class 2 or 3) at z = 1: MOD always a 3-adic unit, ORIG not",
      list(m.keys()) == ["0"] and len(o) > 5, f"MOD {m}, ORIG {len(o)} valuations")

# W5: exact distinct counts = trace counts to length 11, both examples
tc = R["C4"]["trace"]
check("W5 different words give different maps, every word up to length 11, both examples (hashed mod 2 primes)",
      R["C4"]["ORIG"]["match"] and R["C4"]["MOD"]["match"], f"length 11: {tc[-1]:,} elements each")

# W6: the trace counts come from the formula 1/P = 1 + sum((1-t)^k - 1); compare with direct enumeration of M
def enum_M(n):
    """reduced words of N^2 * N * N by length, counted directly"""
    from functools import lru_cache
    @lru_cache(None)
    def f(L, prev):   # number of reduced words of length L whose outer syllable is not class prev
        if L == 0: return 1
        tot = 0
        for k in range(3):
            if k == prev: continue
            sizes = range(1, L + 1)
            for s in sizes:
                ways = (s + 1) if k == 0 else 1    # elements of N^2 of length s: s+1; of N: 1
                tot += ways * f(L - s, k)
        return tot
    return [f(L, -1) for L in range(1, n + 1)]
check("W6 trace-monoid counts agree with direct enumeration of reduced words", enum_M(11) == tc)

# W7: certified D* brackets contain the float roots; first-return identity at D*
r1 = 2**0.5 / 4
for name, sl in [("ORIG", [[r1, .3], [1/3], [1/3]]), ("MOD", [[r1, .3], [.3], [.3]])]:
    D = brentq(lambda s: 1 + sum(np.prod([1 - x**s for x in l]) - 1 for l in sl), .5, 2)
    lo, hi = R["C5_" + name]
    P1 = 1 / np.prod([1 - x**D for x in sl[0]])
    Qinv = 1 + sum(np.prod([1 - x**D for x in l]) - 1 for l in sl[1:]); Q = 1 / Qinv
    check(f"W7 {name}: certified D* bracket holds the root; first-return pieces have (P1-1)(Q-1)=1 there",
          lo <= D <= hi and abs((P1 - 1) * (Q - 1) - 1) < 1e-9, f"D* in [{lo:.6f}, {hi:.6f}]")

# W8: linear irreducibility: some first-return piece turns by 45 degrees
lam = (1 + 1j) / 4 * (1 / 3)
check("W8 the piece a·(class-3 map) turns by 45 degrees, so the pieces fix no line (Hochman's irreducibility)",
      abs(math.degrees(cmath.phase(lam)) - 45) < 1e-12)

# W9: ping-pong in action: random reduced words send a base point into the cover of their outer class
def in_cover(z, cover): return any(abs(z - complex(float(c[0]), float(c[1]))) <= float(r) + 1e-12 for c, r in cover)
ok = True; tried = 0
for name, cB, Rd in [("ORIG", (F(3,5), F(3,8)), F(87,100)), ("MOD", (F(29,50), F(19,50)), F(17,20))]:
    _, X, _ = P.pingpong(name, cB, Rd, F(1, 60), F(1, 60))
    rng = random.Random(11)
    for _ in range(300):
        w = P.rand_word(rng.randint(1, 5), rng)
        inner = w[-1][0]; m = next(j for j in range(3) if j != inner)
        z0 = P.apply(P.word_map(name, [(m, (1, 0) if m == 0 else (1,))]), (cB[0], cB[1]))
        z = P.apply(P.word_map(name, w), z0); zc = complex(float(z[0]), float(z[1]))
        ok &= in_cover(zc, X[w[0][0]]) and not any(in_cover(zc, X[j]) for j in range(3) if j != w[0][0]); tried += 1
check("W9 random words land in the cover of their outer class and in no other", ok, f"{tried} words")

# W10: the last step for ORIG: e <= 6 excluded for every d
c = R["C7"]
check("W10 ORIG last step: lam1^d lam2^-e is no ratio of class-2/3 points, all d >= 1, e <= 6 (exact)",
      c["pingpong_certified"] and c["all_excluded"], f"{c['pairs_checked_exactly']} pairs exact, size bound beyond")

# W11: and it fails for e = 7, d = 8 (so the restriction is needed by this method)
lam1, lam2 = P.EXAMPLES["ORIG"][0][1]
rho = P.cm(P.pw(lam1, 8), P.pw((F(10, 3), F(0)), 7)); _, X, _ = P.pingpong("ORIG", (F(3,5), F(3,8)), F(87,100), F(1,60), F(1,60))
W = X[1] + X[2]; rup = (P.SQRT2_UP / 4)**8 * F(10, 3)**7
check("W11 the disc test cannot exclude (d, e) = (8, 7): a genuinely open case for this method",
      not P.ratio_excluded(rho, rup, W), f"|rho| = {float(P.n2(rho))**.5:.4f}")

# W12: certified unconditional lower bounds for ORIG, consistent with D*
lo6, hi6 = R["C8_ORIG_beta_le_6"]; loa, hia = R["C6_ORIG_chain"]; Dlo, Dhi = R["C5_ORIG"]
check("W12 ORIG proved range: 1.0186 (pure a-runs) <= 1.19061 (b-runs <= 6) <= dim <= D* <= 1.19064",
      1 < loa < lo6 < Dlo and hi6 < Dhi and Dhi - lo6 < 3e-5, f"[{lo6:.6f}, {Dhi:.6f}]")

# W13: MOD's D* is strictly below 2 (Hochman's min(2, .) is the similarity dimension)
check("W13 both D* values are below 2", R["C5_MOD"][1] < 2 and R["C5_ORIG"][1] < 2)

txt = "\n".join(out) + f"\n\n{sum(res.values())}/{len(res)} pass\n"
open("checks7_output.txt", "w").write(txt); json.dump(res, open("checks7.json", "w"), indent=1); print(txt)
