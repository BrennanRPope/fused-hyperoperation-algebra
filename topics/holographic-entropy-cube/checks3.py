"""checks3.py: third pass (DRAFT). Chasing the fuzzy lead, and our own proof of the gluing law.
K  fuzzy (soft) min cuts: OURS (proof), technique from Ellis-Monroe-Newman 1976 (four copies + rotation).
   Soft entropy S_h(T) = -h log(Z_T/Z_0), Z_T = sum over bulk 0/1 labels of exp(-cut/h)  (= free energy cost of
   pinning region T to 1 in a ferromagnetic Ising model; = annealed Renyi-2 of a random tensor network).
   Theorem (ours): S_h obeys monogamy (MMI) for every graph with positive weights and every h > 0.
L  the gluing law in AdS3: OURS (proof), the ingredients (cross-ratio, RT, hyperbolic geometry) are standard.
Run: python3 checks3.py   (needs numpy, scipy, sympy)"""
import itertools, math, json, sys
import numpy as np, sympy as sp
from scipy.optimize import minimize
from scipy.special import logsumexp
rng = np.random.default_rng(20261125)
OUT = {}; NPASS = [0, 0]
def check(name, ok, detail=""):
    NPASS[1] += 1; NPASS[0] += bool(ok)
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""), flush=True)

# ------------------------------------------------------------------ K
print("== K. fuzzy min cuts (ours): monogamy holds at every softness")
P4 = list(itertools.product((1, -1), repeat=4))           # sign vectors
v = np.array([[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]])   # the 'odd frame'
Had = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])  # columns used below = rows (symmetric)
check("K1 the odd frame v1..v4 is orthogonal, |v_i|^2 = 4", (v @ v.T == 4 * np.eye(4)).all())
alpha = {s: tuple((v @ np.array(s)) // 2) for s in P4}
Sprime = set(alpha.values())
even_pm1 = {s for s in P4 if np.prod(s) == 1}
axes = {tuple(int(2 * sg * (i == j)) for j in range(4)) for i in range(4) for sg in (1, -1)}
check("K2 rotated single-site set S' = {even sign vectors} U {+-2 e_i} (16 points)", Sprime == (even_pm1 | axes) and len(Sprime) == 16)
moms = [sum(np.prod([a[i] ** n[i] for i in range(4)]) for a in Sprime) for n in itertools.product(range(7), repeat=4)]
check("K3 every moment E[a1^n1 a2^n2 a3^n3 a4^n4] over S' is >= 0 (all n_i <= 6): Griffiths' condition", min(moms) >= 0, f"{len(moms)} moments, min {min(moms)}")
flip4 = Had.copy(); flip4[:, 3] *= -1        # negate copy 4 in each pinned column
check("K4 negating copy 4 sends the Hadamard frame (even patterns) to (v1, v2, v3, -v4)",
      (flip4 == np.vstack([v[0], v[1], v[2], -v[3]])).all())
check("K5 pinned terminal j sits at alpha = 2 e_j (it feeds only coordinate j)", all(alpha[tuple(v[j])] == tuple(2 * (i == j) for i in range(4)) for j in range(4)))

def graph(rng, nb, dens=0.7, scale=1.0):
    N = 4 + nb; W = rng.exponential(scale, (N, N)) * (rng.random((N, N)) < dens); W = np.triu(W, 1); return W + W.T
def Z4(W, nb, s):
    """sum over bulk spins of exp(sum_{i<j} J_ij s_i s_j), J = W/2 (h = 1); terminals 0..3 pinned to s."""
    N = 4 + nb; J = W / 2; tot = []
    for b in itertools.product((1, -1), repeat=nb):
        x = np.array(list(s) + list(b)); tot.append(0.5 * x @ J @ x)
    return logsumexp(tot)
def softS(W, nb):   # S(T) for T subset of {a,b,c}, terminal 3 = purifier pinned 0, h = 1
    N = 4 + nb; F = {}
    for T in itertools.product((0, 1), repeat=3):
        lab0 = list(T) + [0]; vals = []
        for b in itertools.product((0, 1), repeat=nb):
            lab = np.array(lab0 + list(b)); vals.append(-(W * (lab[:, None] != lab[None, :])).sum() / 2)
        F[T] = -logsumexp(vals)
    return {T: F[T] - F[(0, 0, 0)] for T in F}
def I3(S): return sum((-1) ** (sum(T) + 1) * S[T] for T in S if sum(T))
errs = []; Js = []
for _ in range(300):
    nb = int(rng.integers(1, 5)); W = graph(rng, nb)
    lz = {s: Z4(W, nb, s) for s in P4}
    J = sum(np.prod(s) * lz[s] for s in P4) / 16
    errs.append(abs(I3(softS(W, nb)) - 8 * J)); Js.append(J)
check("K6 soft I3 = 8 h J_abco, J_abco = 4-body coupling of the terminals' marginal (exact reformulation)", max(errs) < 1e-9, f"max err {max(errs):.1e}")
check("K7 J_abco <= 0 on 300 random graphs (so I3 <= 0)", max(Js) < 1e-12, f"max J {max(Js):.2e}")
# the four-copy identity: prod over even patterns = Z_4copy(v1,v2,v3,-v4), odd = Z_4copy(v1..v4); compute 4-copy directly
W = graph(rng, 2); nb = 2
def Z4copy(W, nb, cols):
    tot = 0.0
    for copy in range(4):
        s = tuple(int(cols[i][copy]) for i in range(4)); tot += Z4(W, nb, s)
    return tot
odd = [tuple(int(x) for x in r) for r in (v.T)]  # rows of the odd patterns: copy c gives config (v1[c],..,v4[c])
lhs = Z4copy(W, nb, [v[0], v[1], v[2], v[3]]); rhs = Z4copy(W, nb, [v[0], v[1], v[2], -v[3]])
pe = sum(Z4(W, nb, tuple(r)) for r in Had.T); po = sum(Z4(W, nb, tuple(r)) for r in v.T)
check("K8 four copies: log prod_even Z = log Z4copy(v1,v2,v3,-v4), log prod_odd Z = log Z4copy(v1,v2,v3,v4)",
      abs(rhs - pe) < 1e-9 and abs(lhs - po) < 1e-9)
# rotated model: bulk variables in S', pair energy J * (a.a'), terminal j field 2 J alpha_j -> same partition function
def Zrot(W, nb, sign4):
    J = W / 2; pts = sorted(Sprime); tot = []
    pins = [np.array([2 * (i == j) for i in range(4)]) for j in range(4)]; pins[3] = sign4 * pins[3]
    for bs in itertools.product(pts, repeat=nb):
        X = pins + [np.array(b) for b in bs]
        tot.append(sum(J[i, k] * X[i] @ X[k] for i in range(4 + nb) for k in range(i + 1, 4 + nb)))
    return logsumexp(tot)
check("K9 the rotation keeps the energy: rotated model gives the same Z4copy for both pins",
      abs(Zrot(W, nb, 1) - lhs) < 1e-9 and abs(Zrot(W, nb, -1) - rhs) < 1e-9)
# adversarial search (J <= 0 claim): optimise weights to push I3 up
best = -1
for trial in range(60):
    nb = int(rng.integers(1, 5)); N = 4 + nb; iu = np.triu_indices(N, 1)
    labs = {T: np.array([list(T) + [0] + list(b) for b in itertools.product((0, 1), repeat=nb)]) for T in itertools.product((0, 1), repeat=3)}
    cutm = {T: (L[:, iu[0]] != L[:, iu[1]]).astype(float) for T, L in labs.items()}
    def f(th):
        w = np.exp(th); F = {T: -logsumexp(-(c @ w)) for T, c in cutm.items()}
        S = {T: F[T] - F[(0, 0, 0)] for T in F}; return I3(S) / (sum(abs(x) for x in S.values()) + 1e-9)
    r = minimize(lambda t: -f(t), rng.normal(0, 2, len(iu[0])), method="L-BFGS-B", bounds=[(-8, 5)] * len(iu[0]))
    best = max(best, -r.fun)
check("K10 adversarial search (60 optimised graphs, up to 4 bulk vertices) never makes soft I3 positive", best < 1e-12, f"best normalized I3 {best:.1e}")
x = sp.symbols("x", positive=True)
check("K11 star with 4 legs (Kadanoff decimation): K4 = (1/8)[log cosh 4K - 4 log cosh 2K] <= 0 since cosh^4 2K - cosh 4K = (cosh^2 2K - 1)^2",
      sp.simplify(sp.cosh(2 * x) ** 4 - sp.cosh(4 * x) - (sp.cosh(2 * x) ** 2 - 1) ** 2) == 0)
# Potts: the Z2 symmetry matters. star with q colours: Z_T = sum_c exp(-sum_i w_i [c != colour_i])
def potts_star_I3(w, q):
    out = {}
    for T in itertools.product((0, 1), repeat=3):
        cols = list(T) + [0]; out[T] = sum(math.exp(-sum(wi * (c != ci) for wi, ci in zip(w, cols))) for c in range(q))
    S = {T: -math.log(out[T] / out[(0, 0, 0)]) for T in out}; return I3(S)
bestq = {}
for q in (2, 3, 8, 20):
    bb = -9
    for _ in range(40):
        r = minimize(lambda t: -potts_star_I3(np.exp(t), q), rng.normal(0, 1.5, 4), method="Nelder-Mead"); bb = max(bb, -r.fun)
    bestq[q] = bb
check("K12 Potts (q colours) star: q = 2 never breaks MMI, q = 8 and q = 20 do (the proof needs the two-state symmetry)",
      bestq[2] < 1e-9 and bestq[8] > 1e-7 and bestq[20] > 1e-4, ", ".join(f"q={q}: {b:.1e}" for q, b in bestq.items()))

# an explicit 3-state counterexample (two bulk vertices); terminals a, b, c, o = 0..3, bulk = 4, 5
EDG = list(itertools.combinations(range(6), 2))
W3 = [0.585243, 0.181209, 0.368999, 2.640805, 0.0, 0.62655, 1.726188, 0.0, 0.306167, 1.030605, 0.0, 1.125454, 0.493235, 0.0, 0.112386]
def potts_S(q):
    Z = {}
    for T in itertools.product((0, 1), repeat=3):
        Z[T] = sum(math.exp(-sum(w * (col[u] != col[vv]) for w, (u, vv) in zip(W3, EDG)))
                   for col in (list(T) + [0] + list(b) for b in itertools.product(range(q), repeat=2)))
    return {T: -math.log(Z[T] / Z[(0, 0, 0)]) for T in Z}
i3_q3, i3_q2 = I3(potts_S(3)), I3(potts_S(2))
check("K12b explicit graph (positive weights, 2 bulk vertices): 3-state Potts breaks MMI, the same graph with 2 states does not",
      i3_q3 > 1e-6 and i3_q2 < 0, f"I3: 3 states {i3_q3:.2e}, 2 states {i3_q2:.2e}")
OUT["K12b"] = {"weights": W3, "edges": EDG, "I3_q3": i3_q3, "I3_q2": i3_q2}
OUT["K"] = {"potts_star_best_I3": bestq, "J_samples": Js}

# five-party facets (Hernandez-Cuenca 2019 table items 3-7; item 8 not retrieved): adversarial soft-min search
import re
PP = "ABCDE"
def parse(s): return [(int(m.group(1) or 1), tuple(1 if p in m.group(2) else 0 for p in PP)) for m in (re.match(r"(\d*)(\w+)", t) for t in s.split())]
FAC = {"I3 (MMI uplift)": ("ABC ADE BCDE", "A BC DE ABCDE"), "I4 cyclic": ("ABC ABD ACE BDE CDE", "AB AC BD CE DE ABCDE"),
       "I5": ("ABC ABD ABE ACD ACE ADE BCE BDE CDE", "AB AC AD BE CE DE BCD ABCE ABDE ACDE"),
       "I6": ("3ABC 3ABD ABE ACD 3ACE ADE BCD BCE BDE CDE", "2AB 2AC AD AE BC 2BD 2CE DE 2ABCD 2ABCE ABDE ACDE"),
       "I7": ("2ABC ABD ABE ACD ADE BCE BDE", "AB AC AD BC BE DE ABCD ABCE ABDE")}
FAC = {k: (parse(a), parse(b)) for k, (a, b) in FAC.items()}
worst5 = {}
for name, (L, Rr) in FAC.items():
    bb = -9
    for trial in range(12):
        nb = int(rng.integers(1, 4)); N = 6 + nb; iu = np.triu_indices(N, 1)
        Ts = list(itertools.product((0, 1), repeat=5)); idx = {T: i for i, T in enumerate(Ts)}
        bulk = np.array(list(itertools.product((0, 1), repeat=nb)))
        cut = np.array([((c := np.hstack([np.tile(np.array(T + (0,)), (len(bulk), 1)), bulk]))[:, iu[0]] != c[:, iu[1]]).astype(float) for T in Ts])
        def g(th):
            F = -logsumexp(-(cut @ np.exp(th)), axis=1); S = F - F[0]
            return (sum(k * S[idx[T]] for k, T in L) - sum(k * S[idx[T]] for k, T in Rr)) / (np.abs(S).sum() + 1e-9)
        r = minimize(g, rng.normal(0, 2, len(iu[0])), method="L-BFGS-B", bounds=[(-8, 5)] * len(iu[0])); bb = max(bb, -r.fun)
    worst5[name] = bb
check("K13 five-party holographic facets: adversarial search finds no soft-min violation (evidence, not proof)",
      max(worst5.values()) < 1e-12, ", ".join(f"{k}: {b:.0e}" for k, b in worst5.items()))
OUT["K"]["five_party_worst"] = worst5

# ------------------------------------------------------------------ L
print("== L. the gluing law, proved by shadows (ours)")
a1, a2, b1, m, b2, p = sp.symbols("a1 a2 b1 m b2 p", real=True)
xr = lambda A1, A2, B1, B2: (A2 - A1) * (B2 - B1) / ((B1 - A1) * (B2 - A2))
check("L1 Ptolemy on a line: (a2-a1)(b2-b1) + (b1-a2)(b2-a1) = (b1-a1)(b2-a2)",
      sp.expand((a2 - a1) * (b2 - b1) + (b1 - a2) * (b2 - a1) - (b1 - a1) * (b2 - a2)) == 0)
f = lambda q: (q - a1) / (q - a2)          # the Moebius map sending a1 -> 0, a2 -> infinity
check("L2 shadow lemma: 1 - x(A,B) = f(b2)/f(b1), f the Moebius map a1 -> 0, a2 -> oo", sp.simplify(1 - xr(a1, a2, b1, b2) - f(b2) / f(b1)) == 0)
check("L3 shadows add (telescoping): f(b2)/f(b1) = [f(m)/f(b1)] [f(b2)/f(m)], so 1-x multiplies and x glues by o_{-1}",
      sp.simplify(xr(a1, a2, b1, b2) - (lambda u, w: u + w - u * w)(xr(a1, a2, b1, m), xr(a1, a2, m, b2))) == 0)
def hdist(z, w): return math.acosh(1 + abs(z - w) ** 2 / (2 * z.imag * w.imag))
def geo(A1, A2, u):
    th = math.pi / (1 + math.exp(-max(-30.0, min(30.0, u))))   # u in R -> th in (0, pi)
    c, r = (A1 + A2) / 2, (A2 - A1) / 2; return complex(c + r * math.cos(th), r * math.sin(th))
def foot(A1, A2, q, eps=1e-7):   # nearest point on gamma_A to q + i eps (independent of the Moebius map)
    r = minimize(lambda th: hdist(geo(A1, A2, th[0]), complex(q, eps)), [0.0], method="Nelder-Mead", options={"xatol": 1e-13, "fatol": 1e-15})
    best = r; 
    for s in (-3.0, -1.0, 1.0, 3.0):
        rr = minimize(lambda th: hdist(geo(A1, A2, th[0]), complex(q, eps)), [s], method="Nelder-Mead", options={"xatol": 1e-13, "fatol": 1e-15})
        if rr.fun < best.fun: best = rr
    return geo(A1, A2, best.x[0])
e_sh = []; e_add = []; e_del = []; e_mi = []
for _ in range(60):
    A1, A2, B1, M, B2 = np.sort(rng.random(5) * 10)
    t = lambda u, w: hdist(foot(A1, A2, u), foot(A1, A2, w))
    X0 = (A2 - A1) * (B2 - B1) / ((B1 - A1) * (B2 - A2))
    e_sh.append(abs(math.exp(-t(B1, B2)) - (1 - X0)))
    e_add.append(abs(t(B1, B2) - t(B1, M) - t(M, B2)))
    dl = min(minimize(lambda z: hdist(geo(A1, A2, z[0]), geo(B1, B2, z[1])), [s, u], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14}).fun for s in (-2.0, 0.0, 2.0) for u in (-2.0, 0.0, 2.0))
    e_del.append(abs((1 - X0) - math.tanh(dl / 2) ** 2))
    # RT mutual information on the line (Poincare patch): S(interval) = (c/3) log(len/eps); two intervals: min over pairings
    eps = 1e-3; Sl = lambda u, w: math.log(abs(w - u) / eps)
    SA, SB = Sl(A1, A2), Sl(B1, B2); SAB = min(SA + SB, Sl(A2, B1) + Sl(A1, B2))
    e_mi.append(abs((SA + SB - SAB) - max(0.0, -2 * math.log(math.sinh(dl / 2)))))
check("L4 the shadow length t of B on A's geodesic (feet found by direct minimisation) obeys e^-t = 1 - x", max(e_sh) < 1e-5, f"max err {max(e_sh):.1e}")
check("L5 shadows of neighbouring B1, B2 add: t(B1 u B2) = t(B1) + t(B2)", max(e_add) < 1e-5, f"max err {max(e_add):.1e}")
check("L6 with delta = distance between the two geodesics: x = sech^2(delta/2), 1 - x = tanh^2(delta/2)", max(e_del) < 1e-8, f"max err {max(e_del):.1e}")
check("L7 RT mutual information = (c/3) max(0, log X), X = odds = csch^2(delta/2): connected iff delta < 2 arcsinh 1", max(e_mi) < 1e-6, f"max err {max(e_mi):.1e}")
phi = lambda t: max(0.0, math.log(math.expm1(t))) if t > 0 else 0.0
bad = 0
for _ in range(20000):
    t1, t2 = rng.exponential(1.5, 2); bad += phi(t1 + t2) < phi(t1) + phi(t2) - 1e-12
check("L8 monogamy for A next to B1, B2: phi(t) = max(0, log(e^t - 1)) is superadditive (I3 <= 0), 20000 pairs", bad == 0)
xx, yy = sp.symbols("xx yy", positive=True)
check("L9 the step behind L8: (e^s - 1)(e^u - 1) <= e^(s+u) - 1 because the gap is (e^s - 1) + (e^u - 1) = X o_1 Y - XY",
      sp.simplify((xx * yy + xx + yy) - xx * yy - (xx + yy)) == 0)
# conformal factors cancel (covers the circle and the BTZ line)
w_, r_ = sp.Function("w"), sp.Function("rho")
d = lambda u, s: (w_(s) - w_(u)) * r_(u) * r_(s)
check("L10 conformal factors cancel: with d(u,v) = (w(v)-w(u)) rho(u) rho(v) the cross-ratio only sees w (circle: w = tan(theta/2); BTZ: w = e^(2 pi u/beta))",
      sp.simplify(d(a1, a2) * d(b1, b2) / (d(a1, b1) * d(a2, b2)) - xr(w_(a1), w_(a2), w_(b1), w_(b2))) == 0)

print(f"\n{NPASS[0]}/{NPASS[1]} checks pass")
json.dump(OUT, open("results3.json", "w"), default=float)
