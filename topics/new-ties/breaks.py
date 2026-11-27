"""Round 9 (DRAFT): which parameters break the formula, and by how much (README section 13).

Theorem P10: if two different words give the same map, dim K <= s_c < D*.
Corollary P11: so (for algebraic parameters, or any real parameters on the line given Kittle-Kogler) with D* <= d,
               dim K = D*  exactly when no two different words collide. The breaking set E is a countable union of
               proper algebraic subsets of parameter space (proved here by an algebraic-independence argument).

This script
  B1  a concrete break inside the separated region, on the line: centres 0, 1, 2; class 1 = {19/200, 1/20} about 0,
      class 2 = {1/19} about 1, class 3 = {1/10} about 2. The two pieces (19/200 o 1/19) and (1/20 o 1/10) are the
      same map z -> (18 + z)/200 (exact check). Certified D* and the certified upper bound s_Y (root of
      f(s) - (1/200)^s = 1, f = (P1 - 1)(Q - 1)); certified lower bounds from the collision-free subsystems
      with one class-1 map (pattern 1,1,1: test (B) is empty, (A) holds by Lemma A').
  B2  the counts: distinct maps vs trace-monoid elements by length (hashed mod two primes); the deficit appears at
      length 2 and grows.
  B1' the collision spreads to a whole family of letters, (a lam).(nu w) = (a mu).(kappa w); removing that family gives
      the sharper certified upper bound (root of f(s) - P1(s) Q(s) (mu kappa)^s = 1)
  B2' the family model's counts match the distinct-map counts exactly to length 9
  B3  the bounds side by side
  B4  the drop over the whole collision surface lam = mu (2 - kappa), nu = kappa mu / lam (inside the box
      |.| < 1/9): D* minus the family bound on a grid (for the figure).
  B5  no collisions at algebraically independent parameters: the polynomial argument, checked by evaluating random
      word pairs at random real (float) small parameters and confirming the maps differ (sanity check only; the
      proof is in the README).
Writes breaks.json."""
import json, math, random, itertools
from fractions import Fraction as F
import numpy as np
from mpmath import iv, mp
import plane as P
import general as G
iv.dps = 30

C = [F(0), F(1), F(2)]
LAM, MU, NU, KAP = F(19, 200), F(1, 20), F(1, 19), F(1, 10)
EX = [((F(0), F(0)), [(LAM, F(0)), (MU, F(0))]), ((F(1), F(0)), [(NU, F(0))]), ((F(2), F(0)), [(KAP, F(0))])]

def ivq(q):
    q = F(q); return iv.mpf(q.numerator) / q.denominator

def hom(c, l):  # z -> c + l (z - c) as (slope, intercept)
    return (l, c * (1 - l))
def comp(f, g):  # f o g
    return (f[0] * g[0], f[0] * g[1] + f[1])

def fseries(s, cl1, others):
    """f(s) = (P1 - 1)(Q - 1): sum over first-return pieces u.v of |lam|^s (interval arithmetic)"""
    P1 = iv.mpf(1)
    for l in cl1: P1 /= (1 - ivq(l)**s)
    invQ = iv.mpf(0)
    for ls in others:
        Pj = iv.mpf(1)
        for l in ls: Pj /= (1 - ivq(l)**s)
        invQ += 1 / Pj
    invQ -= (len(others) - 1)
    Q = 1 / invQ
    return (P1 - 1) * (Q - 1)

def Dstar_f(s, classes):
    tot = iv.mpf(1)
    for ls in classes:
        pr = iv.mpf(1)
        for l in ls: pr *= 1 - ivq(l)**s
        tot += pr - 1
    return tot

def certify_decreasing(g, lo, hi):
    """certified root of a decreasing g with g(lo) > 0 > g(hi)"""
    for _ in range(90):
        mid = (lo + hi) / 2; v = g(iv.mpf(mid))
        if v.a > 0: lo = mid
        elif v.b < 0: hi = mid
        else: break
        if hi - lo < 1e-9: break
    assert g(iv.mpf(lo)).a > 0 and g(iv.mpf(hi)).b < 0
    return lo, hi

if __name__ == "__main__":
    OUT = {}
    # B1: the exact collision
    x = comp(hom(C[0], LAM), hom(C[1], NU)); y = comp(hom(C[0], MU), hom(C[2], KAP))
    eps = max(LAM, MU, NU, KAP)
    OUT["B1_collision"] = dict(x=[str(v) for v in x], y=[str(v) for v in y], equal=(x == y),
                               lemmaA=eps * (4 * 2 + 1) < 1, eps=str(eps))
    print("B1 collision", x, y, x == y, "Lemma A' holds:", eps * 9 < 1)
    cl = [[LAM, MU], [NU], [KAP]]
    # D*: the formula's root (increasing-to-decreasing: Dstar_f is increasing in s from negative to positive?)
    g = lambda s: -Dstar_f(s, cl)           # Dstar_f increases from 1 - 4 < 0 at s = 0 to 1 at infinity
    Dlo, Dhi = certify_decreasing(g, 0.05, 1.5)
    # check: f(D*) = 1 bracket
    fD = fseries(iv.mpf(Dlo), [LAM, MU], [[NU], [KAP]]), fseries(iv.mpf(Dhi), [LAM, MU], [[NU], [KAP]])
    # s_Y: root of f(s) - (lam*nu)^s - 1 (remove the duplicate letter)
    lamY = LAM * NU
    gy = lambda s: fseries(s, [LAM, MU], [[NU], [KAP]]) - ivq(lamY)**s - 1
    Ylo, Yhi = certify_decreasing(gy, 0.05, 1.5)
    # the collision spreads: (a lam).(nu w) = (a mu).(kappa w) for every class-1 a and every word w of classes 2, 3,
    # so a whole family of letters, of total weight P1(s) Q(s) (mu kappa)^s, duplicates kept letters
    def P1Q(s):
        P1 = 1 / ((1 - ivq(LAM)**s) * (1 - ivq(MU)**s)); Q = 1 / (1 - ivq(NU)**s - ivq(KAP)**s)
        return P1 * Q
    gfam = lambda s: fseries(s, [LAM, MU], [[NU], [KAP]]) - P1Q(s) * ivq(lamY)**s - 1
    Flo, Fhi = certify_decreasing(gfam, 0.05, 1.5)
    # lower bounds: collision-free subsystems with a single class-1 map (Theorem P8, (B) empty, Lemma A')
    lb = {}
    for name, c1 in [("lam only", [LAM]), ("mu only", [MU])]:
        gl = lambda s, c1=c1: -Dstar_f(s, [c1, [NU], [KAP]])
        lb[name] = certify_decreasing(gl, 0.05, 1.5)
    OUT["B1_numbers"] = dict(Dstar=[Dlo, Dhi], f_at_Dstar=[str(fD[0]), str(fD[1])], s_Y_upper=[Ylo, Yhi], s_family_upper=[Flo, Fhi], drop_family_at_least=Dlo - Fhi,
                             drop_at_least=Dlo - Yhi, lower_bounds=lb)
    print("B1 D* in", (Dlo, Dhi), " s_Y in", (Ylo, Yhi), " family bound", (Flo, Fhi), " drop >=", Dlo - Fhi, " lower bounds", lb)

    # B2: counts
    n = 9
    dc = G.distinct_counts(EX, n); tc = P.trace_counts([2, 1, 1], n)
    import sympy as sp
    t = sp.symbols("t"); P1t = 1 / (1 - t)**2; Qt = 1 / (1 - 2 * t); St = (P1t - 1) * (Qt - 1)
    model = sp.Poly(sp.series(Qt * P1t / (1 - (St - t**2 * P1t * Qt)), t, 0, n + 1).removeO(), t).all_coeffs()[::-1][1:]
    model = [int(v) for v in model]
    OUT["B2_counts"] = dict(lengths=list(range(1, n + 1)), distinct=dc, trace=tc, family_model=model,
                            model_matches=(model == dc), deficit=[t_ - d for t_, d in zip(tc, dc)])
    print("   model   ", model, "matches:", model == dc)
    print("B2 distinct", dc); print("   trace   ", tc)

    # B3: the three upper bounds side by side (Theorem P10's general block bound with n = 1, the single letter, the family)
    OUT["B3_bounds"] = dict(Dstar=[Dlo, Dhi], single_letter=[Ylo, Yhi], family=[Flo, Fhi], lower=max(lb.values()))
    # B4: drop over the collision surface (floats)
    mp.dps = 20
    def droot(fun, lo=0.02, hi=1.5):
        from scipy.optimize import brentq
        return brentq(fun, lo, hi)
    grid = []
    for mu in np.linspace(0.005, 0.11, 43):
        for kap in np.linspace(0.005, 0.11, 43):
            lam = mu * (2 - kap); nu = kap * mu / lam
            if max(lam, mu, nu, kap) >= 1 / 9: continue
            sQ = droot(lambda s: 1 - nu**s - kap**s, 1e-6, 1.5)
            fl = lambda s: (1/((1-lam**s)*(1-mu**s)) - 1) * (1/(1 - nu**s - kap**s) - 1)
            D = droot(lambda s: 1 + (1-lam**s)*(1-mu**s) - 1 + (1-nu**s) - 1 + (1-kap**s) - 1)
            Y = droot(lambda s: fl(s) - (lam*nu)**s / ((1-lam**s)*(1-mu**s)*(1 - nu**s - kap**s)) - 1, sQ + 1e-9, 1.5)
            grid.append([float(mu), float(kap), D, Y])
    OUT["B4_surface"] = grid
    drops = [g_[2] - g_[3] for g_ in grid]
    print("B4 grid points", len(grid), "drop range", min(drops), max(drops))

    # B5: sanity check, random small real parameters: random distinct trace elements give different maps
    rng = random.Random(5); bad = 0; tried = 0
    gens = lambda pr: [(0.0, pr[0]), (0.0, pr[1]), (1.0, pr[2]), (2.0, pr[3])]
    letters_class = [0, 0, 1, 2]
    def normal(word):  # canonical form: within each maximal class-0 run, sort letters (they commute); same-class
        out = []; run = []
        for a in word:
            if letters_class[a] == 0: run.append(a)
            else:
                out += sorted(run); run = []; out.append(a)
        return tuple(out + sorted(run))
    for t in range(3000):
        pr = [rng.uniform(0.01, 0.11) for _ in range(4)]
        w1 = [rng.randrange(4) for _ in range(rng.randint(1, 7))]; w2 = [rng.randrange(4) for _ in range(rng.randint(1, 7))]
        if normal(w1) == normal(w2): continue
        def ev(w):
            s, b = 1.0, 0.0
            for a in reversed(w):          # apply innermost first: word read outer->inner
                c, l = gens(pr)[a]; s, b = l * s, l * b + c * (1 - l)
            return s, b
        tried += 1
        (s1, b1), (s2, b2) = ev(w1), ev(w2)
        if abs(s1 - s2) < 1e-14 and abs(b1 - b2) < 1e-14: bad += 1
    OUT["B5_random_real"] = dict(pairs=tried, equal_maps=bad)
    print("B5 random real parameters:", tried, "pairs of different words,", bad, "equal maps")
    json.dump(OUT, open("breaks.json", "w"), indent=1, default=str)
