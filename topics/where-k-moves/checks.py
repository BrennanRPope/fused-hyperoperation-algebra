"""Checks for where-k-moves/README.md (lead R4): a forgotten hidden state, frozen or running.

Plane law (the plane dial): F(x, y) = (x + y + k x y) / (1 - d x y).  Its invariant differential is
dx / (1 + k x + d x^2), so a law conjugate to addition, F(x, y) = phi(psi(x) + psi(y)), sits on the plane
dial exactly when q(x) := dF/dy(x, 0) is a quadratic 1 + k x + d x^2 (a Riccati equation for phi).
"""
import json
import math

import numpy as np
import sympy as sp
import mpmath as mp

results = []
store = {}


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("   " + detail) if detail else ""))


x, y, z, k, d = sp.symbols("x y z k d")


def plane(a, b, kk, dd):
    return (a + b + kk * a * b) / (1 - dd * a * b)


# ---------------------------------------------------------------------------------------------- A. frozen bath,
# chances proportional to the hidden level h: event i hits with chance a_i * h, independently given h.
a1, a2, a3 = sp.symbols("a1 a2 a3", positive=True)
hs = sp.symbols("h0:3", nonnegative=True)
ws = sp.symbols("w0:3", positive=True)


def union_A(doses, levels, weights):
    tot = 0
    for h, w in zip(levels, weights):
        miss = 1
        for a in doses:
            miss *= (1 - a * h)
        tot += w * (1 - miss)
    return sp.expand(tot)


def moments(levels, weights, n):
    return sum(w * h ** n for h, w in zip(levels, weights))


# A1: pairwise k is flat (any doses) and equals -E[h^2]/E[h]^2, for any number of levels
lv, wt = hs, ws
m1, m2, m3 = (moments(lv, wt, n) for n in (1, 2, 3))
p1, p2, p3 = a1 * m1, a2 * m1, a3 * m1
u12 = union_A([a1, a2], lv, wt)
kpair = sp.simplify((u12 - p1 - p2) / (p1 * p2))
check("A1 frozen bath, any number of levels: two events glue by p + q + k p q with k = -E[h^2]/E[h]^2, "
      "the same for every pair of doses (k is flat)", sp.simplify(kpair + m2 / m1 ** 2) == 0,
      "k = " + str(sp.simplify(kpair)))

# A2: three events: the law closes (associative, same k) iff E[h^3] E[h] = E[h^2]^2
u123 = union_A([a1, a2, a3], lv, wt)
kk_ = -m2 / m1 ** 2
glued = sp.expand(plane(plane(p1, p2, kk_, 0), p3, kk_, 0))
defect = sp.simplify((u123 - glued) / (p1 * p2 * p3))
target = (m3 * m1 - m2 ** 2) / m1 ** 4
check("A2 three events: (true union) - (p1 o p2 o p3) = p1 p2 p3 (E[h^3]E[h] - E[h^2]^2)/E[h]^4",
      sp.simplify(defect - target) == 0)

# A3: Cauchy-Schwarz: E[h^2]^2 <= E[h] E[h^3] with equality iff h takes one nonzero value (plus 0)
rng = np.random.default_rng(0)
worst_gap_multi, cs_ok = 1e9, True
for _ in range(20000):
    nlev = rng.integers(2, 5)
    hv = rng.uniform(0.05, 1, nlev); wv = rng.dirichlet(np.ones(nlev))
    M1, M2, M3 = (wv * hv).sum(), (wv * hv ** 2).sum(), (wv * hv ** 3).sum()
    gap = M1 * M3 - M2 ** 2
    cs_ok &= gap >= -1e-15
    worst_gap_multi = min(worst_gap_multi, gap / M1 ** 4)
cure = union_A([a1, a2, a3], [0, hs[1]], [1 - ws[1], ws[1]])
pi_ = ws[1]
q1, q2, q3 = a1 * pi_ * hs[1], a2 * pi_ * hs[1], a3 * pi_ * hs[1]
glued_cure = sp.expand(plane(plane(q1, q2, -1 / pi_, 0), q3, -1 / pi_, 0))
check("A3 Cauchy-Schwarz: the defect is >= 0 always and zero only for one 'on' level plus an 'off' level "
      "(the mixture cure model, k = -1/pi); random 2-4 level baths with all levels on never close",
      cs_ok and worst_gap_multi > 0 and sp.simplify(cure - glued_cure) == 0,
      f"smallest defect over 20000 random baths with every level on: {worst_gap_multi:.2e} > 0")
store["A_defect_examples"] = []
for hv, wv in [([0, 1], [0.6, 0.4]), ([0.5, 1], [0.5, 0.5]), ([0, 0.5, 1], [0.4, 0.3, 0.3]),
               ([0.2, 0.6, 1], [1 / 3] * 3)]:
    hv, wv = np.array(hv), np.array(wv)
    M1, M2, M3 = (wv * hv).sum(), (wv * hv ** 2).sum(), (wv * hv ** 3).sum()
    store["A_defect_examples"].append({"levels": hv.tolist(), "weights": wv.tolist(),
                                       "k": -M2 / M1 ** 2, "defect": (M1 * M3 - M2 ** 2) / M1 ** 4})

# ---------------------------------------------------------------------------------------------- B. frozen bath,
# exposures add hazard: P(no event | h) = exp(-h * Lambda), and two exposures add their Lambda.
s = sp.symbols("s", nonnegative=True)


def law_from_laplace(L, Linv):
    """p = 1 - L(Lambda); F(p1, p2) = 1 - L(Linv(1-p1) + Linv(1-p2))."""
    return lambda P1, P2: 1 - L(Linv(1 - P1) + Linv(1 - P2))


def q_of_laplace(Lexpr):
    """q(p) = dF/dy(p, 0) = L'(Lambda)/L'(0) written as a function of p = 1 - L(Lambda)."""
    return sp.diff(Lexpr, s) / sp.diff(Lexpr, s).subs(s, 0)


# B1: any frailty: associative (numeric, 3-level frailty)
hv = np.array([0.3, 1.0, 2.5]); wv = np.array([0.2, 0.5, 0.3])
Ln = lambda t: float((wv * np.exp(-hv * t)).sum())
Linvn = lambda u: float(mp.findroot(lambda t: (wv * np.exp(-hv * float(t))).sum() - u, 0.5))
Fn = law_from_laplace(Ln, Linvn)
worst = 0
for _ in range(200):
    P = rng.uniform(0.01, 0.6, 3)
    worst = max(worst, abs(Fn(Fn(P[0], P[1]), P[2]) - Fn(P[0], Fn(P[1], P[2]))))
check("B1 frozen bath where exposures add hazard: a 3-level bath gives an associative law "
      "(it is conjugate to adding exposures)", worst < 1e-10, f"worst associativity gap {worst:.1e}")

# B2: geometric frailty (levels 1, 2, 3, ... with weights (1-t) t^(n-1)) = the Ali-Mikhail-Haq law = the plane dial
t, pz = sp.symbols("theta pi0", positive=True)
Lgeo = (1 - t) * sp.exp(-s) / (1 - t * sp.exp(-s))
u = sp.symbols("u", positive=True)
# S = L(Lambda), e^{-Lambda} = S / (1 - t + t S)
Sgeo = lambda P1, P2: None
P1, P2 = sp.symbols("P1 P2", positive=True)
e1 = (1 - P1) / (1 - t + t * (1 - P1)); e2 = (1 - P2) / (1 - t + t * (1 - P2))
F_geo = sp.simplify(1 - (1 - t) * e1 * e2 / (1 - t * e1 * e2))
ok_geo = sp.simplify(F_geo - plane(P1, P2, -(1 + t), t)) == 0
check("B2 geometric frailty (levels 1, 2, 3, ... weights (1-theta) theta^(n-1)) gives exactly the plane law "
      "with k = -(1 + theta), d = theta (the Ali-Mikhail-Haq copula)", ok_geo)

# B3: zero-modified geometric (cure share 1 - pi0 on top): still the plane dial
Lzm = (1 - pz) + pz * Lgeo
qzm = sp.simplify(q_of_laplace(Lzm))
pp = sp.symbols("p")
# express q in p: e^{-s} = E solves 1 - p = (1 - pz) + pz (1-t)E/(1-tE)
E = sp.symbols("E", positive=True)
Esol = sp.solve(sp.Eq(1 - pp, (1 - pz) + pz * (1 - t) * E / (1 - t * E)), E)[0]
q_in_p = sp.simplify(qzm.subs(sp.exp(-s), E).subs(E, Esol))
poly = sp.Poly(sp.simplify(sp.together(q_in_p)).as_numer_denom()[0], pp)
num, den = sp.fraction(sp.simplify(sp.together(q_in_p)))
is_quad = sp.degree(sp.Poly(sp.expand(num), pp)) <= 2 and sp.degree(sp.Poly(sp.expand(den), pp)) == 0
qq = sp.expand(sp.simplify(q_in_p))
kz, dz = sp.simplify(qq.coeff(pp, 1)), sp.simplify(qq.coeff(pp, 2))
check("B3 add a cured share (zero-modified geometric frailty): still a plane law", is_quad and sp.simplify(qq.coeff(pp, 0) - 1) == 0,
      f"k = {kz}, d = {dz}")
store["B3"] = {"k": str(kz), "d": str(dz)}

# B4: gamma frailty (Clayton): plane only at theta = 1 (exponential frailty), where it touches the parabola
th = sp.symbols("th", positive=True)
Lgam = (1 + th * s) ** (-1 / th)
qg = sp.simplify(q_of_laplace(Lgam))   # = L^(1+th) = (1-p)^(1+th)
qg_p = sp.simplify(qg.subs(s, ((1 - pp) ** (-th) - 1) / th))
ok_b4 = True
for thv in (0.5, 1.0, 2.0):
    for pv in (0.1, 0.4, 0.8):
        sv = ((1 - pv) ** (-thv) - 1) / thv
        ok_b4 &= abs(float(qg.subs({th: thv, s: sv})) - (1 - pv) ** (1 + thv)) < 1e-12
check("B4 gamma frailty (Clayton copula): q(p) = (1 - p)^(1 + theta), a quadratic only at theta = 1 "
      "(exponential frailty: k = -2, d = 1, on the parabola k^2 = 4d)",
      ok_b4 and sp.expand((1 - pp) ** 2) == 1 - 2 * pp + pp ** 2)
check("B5 the geometric line k + d = -1 touches the parabola k^2 = 4d exactly at theta = 1",
      sp.factor((1 + t) ** 2 - 4 * t) == (t - 1) ** 2)

# B6: the bilinear term of a hazard-exposure law is also -E[h^2]/E[h]^2; for the geometric ladder that is -(1 + theta)
Eh1, Eh2 = 1 / (1 - t), (1 + t) / (1 - t) ** 2
n_ = sp.Symbol("n", positive=True, integer=True)
ok_m = all(abs(sum(nn * (1 - tv) * tv ** (nn - 1) for nn in range(1, 4000)) - 1 / (1 - tv)) < 1e-9 and
           abs(sum(nn * nn * (1 - tv) * tv ** (nn - 1) for nn in range(1, 4000)) - (1 + tv) / (1 - tv) ** 2) < 1e-7
           for tv in (0.2, 0.5, 0.8))
check("B6 for hazard exposures the bilinear term is again -E[h^2]/E[h]^2: for the geometric ladder -(1 + theta), "
      "matching B2", ok_m and sp.simplify(-Eh2 / Eh1 ** 2 + (1 + t)) == 0)

# ---------------------------------------------------------------------------------------------- C. running hidden
# state: crossover interference.  Map distance adds along the chromosome, so r(D1 + D2) = F(r1, r2) is always a law.
Dm = sp.symbols("D", nonnegative=True)
rH = (1 - sp.exp(-2 * Dm)) / 2
rK = sp.tanh(2 * Dm) / 2
rC = (1 - sp.exp(-4 * Dm) * (1 + 2 * Dm)) / 2      # chi-square model m = 1 (every second chiasma)
check("C1 Haldane map (no interference): r' = 1 - 2r, the plane law k = -2, d = 0",
      sp.simplify(sp.diff(rH, Dm) - (1 - 2 * rH)) == 0)
check("C2 Kosambi map: r' = 1 - 4 r^2, the plane law k = 0, d = -4: a bend, not a new law",
      sp.simplify(sp.diff(rK, Dm) - (1 - 4 * rK ** 2)) == 0)
# C3: Felsenstein family: coincidence linear in r <=> Riccati <=> plane dial on the line d = -2k - 4
K_ = sp.symbols("K")
kF, dF = -2 * K_, 4 * (K_ - 1)
rr = sp.symbols("r")
check("C3 every plane-dial map function with r -> 1/2 lies on the line d = -2k - 4; its coincidence "
      "c(r) = (1 - r')/(2r) = -(k + d r)/2 is linear (Felsenstein 1979: c = K - 2(K - 1) r)",
      sp.simplify((1 + kF * rr + dF * rr ** 2).subs(rr, sp.Rational(1, 2))) == 0
      and sp.simplify(-(kF + dF * rr) / 2 - (K_ - 2 * (K_ - 1) * rr)) == 0)
# C4: the chi-square m = 1 map is not on the plane dial: r' as a function of r is not a quadratic
Dg = np.linspace(0.0, 3.0, 3001)
rc = 0.5 * (1 - np.exp(-4 * Dg) * (1 + 2 * Dg))
qc = np.exp(-4 * Dg) * (1 + 4 * Dg)
coef = np.polyfit(rc, qc, 2)
resid = np.max(np.abs(np.polyval(coef, rc) - qc))
c_chi = np.where(rc > 1e-9, (1 - qc) / (2 * np.maximum(rc, 1e-12)), np.nan)
check("C4 chi-square m = 1 (the running two-phase hidden state): r' is not a quadratic in r, so it is off the "
      "plane dial", resid > 1e-3, f"best quadratic misses by {resid:.3f} (r' runs 0 to 1)")
store["C4"] = {"best_quadratic": coef.tolist(), "max_miss": float(resid)}
# C5: its coincidence starts at 0, so it is close to Kosambi at short range; compare the k-drift
kk = lambda rf, D: (rf(2 * D) - 2 * rf(D)) / rf(D) ** 2
rfC = lambda D: 0.5 * (1 - np.exp(-4 * D) * (1 + 2 * D))
rfK = lambda D: 0.5 * np.tanh(2 * D)
Dq = np.linspace(0.02, 1.5, 200)
xs = rfK(Dq)
plane_drift = (0 + 2 * (-4) * xs) / (1 + 4 * xs ** 2)
check("C5 a plane law also 'drifts' when read through one-dial k: k_eff(x) = (k + 2dx)/(1 - dx^2); "
      "Kosambi runs from 0 to -2 exactly like the chi-square model, so drift alone does not mean leaving the family",
      np.max(np.abs(kk(rfK, Dq) - plane_drift)) < 1e-9 and abs(plane_drift[-1] + 2) < 0.01)
store["C5"] = {"D": Dq.tolist(), "k_chi": kk(rfC, Dq).tolist(), "k_kos": kk(rfK, Dq).tolist(),
               "x_chi": rfC(Dq).tolist(), "x_kos": xs.tolist()}
store["C_coincidence"] = {"r": rc[1::30].tolist(), "c_chi": c_chi[1::30].tolist()}

# ---------------------------------------------------------------------------------------------- D. Haldane / Wu
# one-state weight Z_g(x) = 1 + 1/w with w^g (1 + w)^(1 - g) = 1/x (Wu 1994).  Frozen mixtures of honest
# elements have counts >= 0; check the counts of Z_g.
mp.mp.dps = 40


def Zg_series(g, N=24):
    # series of Z in x via mpmath taylor of Z(x) = 1 + 1/w(x), w solved by Newton near x -> 0 (w ~ 1/x)
    def Z(xv):
        if xv == 0:
            return mp.mpf(1)
        f = lambda wv: g * mp.log(wv) + (1 - g) * mp.log(1 + wv) + mp.log(xv)
        w0 = 1 / xv
        return 1 + 1 / mp.findroot(f, w0)
    return mp.taylor(Z, 0, N, method="step", h=mp.mpf("1e-4")) if False else None


# do it exactly instead: put v = 1/w; then v^(-g) (1 + 1/v)^(1-g) = 1/x  <=>  x = v (1 + v)^(g - 1)
# so v = x (1 + v)^(1 - g): Lagrange inversion gives [x^n] v = (1/n) [v^(n-1)] (1+v)^(n(1-g)) = binom(n(1-g), n-1)/n
def Zg_coeffs(g, N):
    cs = [mp.mpf(1)]
    for n in range(1, N + 1):
        cs.append(mp.binomial(n * (1 - g), n - 1) / n)
    return cs


negs = {}
for g in [0, mp.mpf(1) / 4, mp.mpf(1) / 3, mp.mpf(1) / 2, mp.mpf(2) / 3, mp.mpf(3) / 4, 1]:
    cs = Zg_coeffs(g, 40)
    first = next((n for n, c in enumerate(cs) if c < -mp.mpf("1e-30")), None)
    negs[str(float(g))] = {"first_negative_count": first, "coeffs": [float(c) for c in cs[:12]]}
check("D1 Wu's one-state weight: Bose (g = 0) and Fermi (g = 1) have honest counts; every g strictly "
      "between (1/4, 1/3, 1/2, 2/3, 3/4) has a negative count, so it is not a frozen mixture of honest states",
      negs["0.0"]["first_negative_count"] is None and negs["1.0"]["first_negative_count"] is None
      and all(negs[str(float(g))]["first_negative_count"] is not None for g in [1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4]),
      "first negative count at x^n: " + ", ".join(f"g={kk_}: n={v['first_negative_count']}" for kk_, v in negs.items()))
store["D"] = negs
# D2: semion: Z = (x/2 + sqrt(1 + x^2/4))^2 = exp(2 asinh(x/2)); (ln Z)' = 1/sqrt(1 + x^2/4), a square root of a
# quadratic, not a quadratic: so semions are not even on the bent (plane) family of one-state weights
Zs = (x / 2 + sp.sqrt(1 + x ** 2 / 4)) ** 2
vv = sp.symbols("v", positive=True)
check("D2 semions (g = 1/2): Z = (x/2 + sqrt(1 + x^2/4))^2, and (ln Z)' = 1/sqrt(1 + x^2/4)",
      sp.simplify(sp.diff(sp.log(Zs), x) - 1 / sp.sqrt(1 + x ** 2 / 4)) == 0
      and all(abs(float((Zs - 1).subs(x, vn * (1 + vn) ** -0.5)) - vn) < 1e-12 for vn in (0.3, 1.0, 2.5)))
# D3: the one-dial k of Wu's statistics runs from 2g - 1 (sparse) to g (crowded)
gg = 0.3
w_sparse, w_crowd = 1e8, 1e-10
kw = lambda w: w + gg - w ** gg * (1 + w) ** (1 - gg)
check("D3 Wu's k = 1/n - 1/x runs from 2g - 1 when sparse to g when crowded (g = 0.3)",
      abs(kw(w_sparse) - (2 * gg - 1)) < 1e-6 and abs(kw(w_crowd) - gg) < 1e-3)

print("ALL PASS" if all(ok for _, ok in results) else "SOME FAIL", f"({sum(ok for _, ok in results)} of {len(results)})")
json.dump(store, open("checks.json", "w"), indent=1, default=str)
