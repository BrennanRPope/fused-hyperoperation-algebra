"""DRAFT checks: how quantum path sums aggregate into classical laws, and path sums as
semiring aggregations, read through the circle operation a o_k b = a + b + k*a*b.

Run:  python3 checks.py  > checks_output.txt   (writes results.json for figures.py)
Needs numpy, scipy, sympy, mpmath.
"""
import json, math, cmath, itertools, random
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import logsumexp

R = {}
passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
    else:
        failed += 1
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))


def circ(a, b, k):
    return a + b + k * a * b


def lift(a, k):
    return 1 + k * a


rng = np.random.default_rng(20261125)

# ---------------------------------------------------------------------------
print("\n== A. Along one path: composing phases by o_k is phase addition ==")
# A1: with weights w = exp(i S / hbar) and dial k = i/hbar, a = (w - 1)/k; chaining by o_k
#     reproduces exp(i * sum S / hbar). The lift is an isomorphism, so nothing new is added.
hbar = 0.37
k = 1j / hbar
S = rng.normal(size=12)
a = (np.exp(1j * S / hbar) - 1) / k
acc = 0
for aj in a:
    acc = circ(acc, aj, k)
err = abs(lift(acc, k) - np.exp(1j * S.sum() / hbar))
check("A1 chaining 12 segment phases by o_k at k=i/hbar = exp(i*sum(S)/hbar)", err < 1e-12, f"err={err:.1e}")
# A1b: in the action coordinate s, a = E(k s)/k with E(x)=e^x-1 (the height-flow map) is the
#      lift's inverse; o_k becomes plain + on actions.
s1, s2 = 0.3, -1.1
lhs = circ(np.expm1(k * s1) / k, np.expm1(k * s2) / k, k)
rhs = np.expm1(k * (s1 + s2)) / k
check("A1b a = E(k s)/k turns o_k into + on actions (s1 o s2 -> s1 + s2)", abs(lhs - rhs) < 1e-13)
# A2: the tan law. o_k at imaginary k, seen through the slide x -> x/(1 - (k/2) x) (Cayley chart),
#     is t1 (+) t2 = (t1+t2)/(1 - t1 t2 * hbar^-2/4)... check with hbar = 1: t = 2 tan(theta/2).
th1, th2 = 0.7, 1.9
t = lambda th: 2 * math.tan(th / 2)
tanlaw = (t(th1) + t(th2)) / (1 - t(th1) * t(th2) / 4)
check("A2 phases in the chart t=2tan(theta/2) compose by the tan law (a+b)/(1-ab/4)",
      abs(tanlaw - t(th1 + th2)) < 1e-12)
# the same law from o_k at k=i via the slide c = k/2: x -> x/(1 + c x) maps (p,q)=(k,0) to (0, 1/4)
kk = 1j
x1 = (cmath.exp(1j * th1) - 1) / kk
x2 = (cmath.exp(1j * th2) - 1) / kk
slide = lambda x: x / (1 + (kk / 2) * x)
check("A2b that chart is o_k at k=i after the slide x/(1+(k/2)x) (Cayley); slid values are real",
      abs(slide(circ(x1, x2, kk)) - t(th1 + th2)) < 1e-12 and abs(slide(x1).imag) < 1e-12)
# A3: real k never gives the circle class: discriminant k^2 - 4d = k^2 > 0 for real k (d=0).
check("A3 real k (any sign) is the line class, discriminant k^2>0; the tan/circle class needs k^2<0",
      all((kv ** 2) > 0 for kv in [-3, -1, -0.2, 0.5, 1, 4]) and (1j) ** 2 == -1)

# ---------------------------------------------------------------------------
print("\n== B. Summing over paths: stationary phase picks Newton's path (quartic, non-Gaussian) ==")
lam = 0.5
V = lambda x: 0.5 * x ** 2 + lam * x ** 4
dV = lambda x: x + 4 * lam * x ** 3
xa, xb, T = -0.6, 0.9, 0.8


def shoot(p0):
    return solve_ivp(lambda t_, y: [y[1], -dV(y[0])], [0, T], [xa, p0],
                     rtol=1e-11, atol=1e-12, dense_output=True)


ps = np.linspace(-15, 15, 600)
vals = [shoot(p).y[0, -1] - xb for p in ps]
roots = [brentq(lambda p: shoot(p).y[0, -1] - xb, ps[i], ps[i + 1])
         for i in range(len(ps) - 1) if vals[i] * vals[i + 1] < 0]
p0 = min(roots, key=abs)
sol = shoot(p0)
xmid = float(sol.sol(T / 2)[0])
print(f"   Newton paths from {xa} to {xb} in T={T}: initial momenta {np.round(roots,4)}; using p0={p0:.6f}, midpoint {xmid:.6f}")

L, N = 10.0, 16384
X = np.linspace(-L / 2, L / 2, N, endpoint=False)
dx = X[1] - X[0]
P = 2 * np.pi * np.fft.fftfreq(N, dx)


def evolve(psi, t_, hb, dt=5e-4, sign=1, record=None):
    """Split-step chain: each step multiplies phases e^{-i V dt/2hb} e^{-i hb p^2 dt/2} e^{-i V dt/2hb}.
    That is the o_k chain at k = -i/hb on the weights (why-i section 5); the sum over paths is the FFT."""
    n = int(round(t_ / dt)); dt = t_ / n
    eV = np.exp(-1j * sign * V(X) * dt / (2 * hb))
    eK = np.exp(-1j * sign * hb * P ** 2 * dt / 2)
    out = []
    for i in range(n):
        psi = eV * psi
        psi = np.fft.ifft(eK * np.fft.fft(psi))
        psi = eV * psi
        if record is not None and (i + 1) % record == 0:
            out.append(psi.copy())
    return (psi, out) if record is not None else psi


Bres = {"xmid_newton": xmid, "roots": roots, "hbars": [], "stat_pt": [], "curv_hbar": [], "window_frac": [],
        "spiral": {}, "integrand": {}}
for hb in [0.08, 0.04, 0.02, 0.01]:
    w = 0.2 * hb  # narrow start: momentum spread hb/(2w) = 2.5, so many paths are summed, not just one
    g = lambda x0: np.exp(-(X - x0) ** 2 / (4 * w ** 2)).astype(complex)
    pf = evolve(g(xa), T / 2, hb)
    pb = evolve(g(xb), T / 2, hb, sign=-1)
    f = np.conj(pb) * pf  # contribution of all paths through midpoint X at time T/2
    df = np.gradient(f, dx)
    kloc = np.imag(df / (f + 1e-300))  # (p_forward - p_backward)/hbar
    near = np.where(np.abs(X - xmid) < 0.6)[0]
    sgn = np.sign(kloc[near])
    zc = near[np.where(np.diff(sgn) != 0)[0]]
    zs = X[zc] + dx * np.abs(kloc[zc]) / (np.abs(kloc[zc]) + np.abs(kloc[zc + 1]))
    xs = zs[np.argmin(np.abs(zs - xmid))]
    i0 = np.argmin(np.abs(X - xs))
    curv = (kloc[i0 + 3] - kloc[i0 - 3]) / (6 * dx)  # d^2 phase / dx^2
    amp = np.sum(f) * dx
    win = np.abs(X - xmid) < 3 * np.sqrt(hb)
    Bres["hbars"].append(hb); Bres["stat_pt"].append(float(xs)); Bres["curv_hbar"].append(float(curv * hb))
    Bres["window_frac"].append(float(abs(np.sum(f[win]) * dx / amp)))
    if hb == 0.01:
        regs = {}
        for lo, hi in [(-1.0, -0.3), (-0.3, 0.65), (0.65, 1.5)]:
            m_ = (X >= lo) & (X < hi)
            regs[f"{lo}..{hi}"] = (float(abs(np.sum(f[m_]) * dx / amp)), float(np.abs(f[m_]).max() / np.abs(f).max()))
        Bres["regions"] = regs
    sel = slice(np.argmin(np.abs(X + 1.2)), np.argmin(np.abs(X - 1.5)))
    norm = np.abs(f).max()
    Bres["integrand"][str(hb)] = {"x": X[sel][::2].tolist(), "re": (f[sel].real / norm)[::2].tolist()}
    cs = np.cumsum(f[sel]) * dx / abs(amp)
    Bres["spiral"][str(hb)] = {"re": cs.real[::2].tolist(), "im": cs.imag[::2].tolist(),
                               "x": X[sel][::2].tolist()}
    print(f"   hbar={hb}: stationary point {xs:.5f} (Newton {xmid:.5f}); curvature*hbar={curv*hb:.4f}; "
          f"window +-3sqrt(hbar) gives {Bres['window_frac'][-1]:.3f} of the amplitude")
check("B1 the phase of the summed contributions is stationary at Newton's midpoint (all four hbar)",
      max(abs(s_ - xmid) for s_ in Bres["stat_pt"]) < 2 * dx,
      f"max offset {max(abs(s_ - xmid) for s_ in Bres['stat_pt']):.1e}, grid {dx:.1e}")
cv = np.array(Bres["curv_hbar"])
check("B2 phase curvature grows like 1/hbar (curvature*hbar constant to 3%), so the useful window shrinks like sqrt(hbar)",
      np.ptp(cv) / abs(cv.mean()) < 0.03, f"curv*hbar = {np.round(cv,4)}")
check("B3 a window of +-3 sqrt(hbar) around Newton's midpoint gives the whole amplitude to within 10% at every hbar",
      max(abs(w_ - 1) for w_ in Bres["window_frac"]) < 0.1, f"{np.round(Bres['window_frac'],3)}")
rg = Bres["regions"]
check("B5 cancelling: at hbar=0.01 the strips either side carry contributions ~90% as big as at Newton's point, yet each nets <5% of the amplitude; the middle strip nets ~100%",
      all(v[0] < 0.05 and v[1] > 0.8 for k_, v in rg.items() if k_ != "-0.3..0.65") and abs(rg["-0.3..0.65"][0] - 1) < 0.08,
      str({k_: (round(v[0], 3), round(v[1], 3)) for k_, v in rg.items()}))
# B4: Ehrenfest: a coherent packet's centre follows Newton; the gap shrinks in proportion to hbar.
Tlong = 2.0
solE = solve_ivp(lambda t_, y: [y[1], -dV(y[0])], [0, Tlong], [xa, p0], rtol=1e-11, atol=1e-12, dense_output=True)
Bres["ehrenfest"] = {}
devs = []
for hb in [0.04, 0.02, 0.01, 0.005]:
    s0 = math.sqrt(hb / 2)
    psi = (np.exp(-(X - xa) ** 2 / (4 * s0 ** 2)) * np.exp(1j * p0 * X / hb)).astype(complex)
    psi /= np.sqrt(np.sum(abs(psi) ** 2) * dx)
    _, snaps = evolve(psi, Tlong, hb, dt=1e-3, record=20)
    ts = np.arange(1, len(snaps) + 1) * 20 * 1e-3
    xm = np.array([np.sum(X * abs(s_) ** 2) * dx for s_ in snaps])
    xn = solE.sol(ts)[0]
    devs.append(float(np.max(abs(xm - xn))))
    Bres["ehrenfest"][str(hb)] = {"t": ts.tolist(), "xq": xm.tolist(), "xn": xn.tolist()}
slope = np.polyfit(np.log([0.04, 0.02, 0.01, 0.005]), np.log(devs), 1)[0]
Bres["ehrenfest_dev"] = devs; Bres["ehrenfest_slope"] = float(slope)
check("B4 Ehrenfest: packet centre minus Newton's x(t) over t<=2 shrinks like hbar^1 (slope 0.9-1.1)",
      0.9 < slope < 1.1, f"max gaps {np.round(devs,5)}, slope {slope:.3f}")
R["B"] = Bres

# ---------------------------------------------------------------------------
print("\n== C. The dial: classical is |k| -> infinity; k -> 0 is the opposite end (counting) ==")
# A tiny route network: paths with actions S_p.  Z(k) = sum_p exp(k S_p).
Sp = np.array([1.00, 1.35, 1.6, 2.2, 2.9, 3.0])
Z = lambda kv: np.sum(np.exp(kv * Sp))
check("C1 k=0: every path weighs 1, the path sum just counts paths (6)", abs(Z(0) - 6) < 1e-15)
soft = lambda kv: logsumexp(kv * Sp) / kv
gaps = [abs(soft(-b) - Sp.min()) for b in [1, 10, 100, 1000]]
check("C2 k -> -infinity (real): (1/k) log Z -> min action (least action); gap <= log(6)/|k|",
      all(g <= math.log(6) / b + 1e-12 for g, b in zip(gaps, [1, 10, 100, 1000])) and gaps[-1] < 1e-3,
      f"gaps {np.round(gaps,5)}")
check("C3 k -> +infinity: -> max (Viterbi / most-likely path after a sign flip)",
      abs(soft(1000) - Sp.max()) < 1e-3)
# C4: in the relativity reading the composed quantity is a velocity with k = -1/c^2, and k -> 0 is
#     Galilean.  In the path reading the composed quantity is a weight with k = i/hbar, and the
#     classical limit is hbar -> 0, i.e. |k| -> infinity.  Opposite ends of the dial.
v1, v2 = 0.6, 0.7
c = 1.0
rap = math.tanh(math.atanh(v1) + math.atanh(v2))
check("C4 velocities compose by (v1+v2)/(1+v1v2/c^2) (rapidities add); c->inf (dial->0) gives Galilean v1+v2",
      abs((v1 + v2) / (1 + v1 * v2) - rap) < 1e-14 and abs((v1 + v2) / (1 + v1 * v2 / 1e12) - (v1 + v2)) < 1e-11)
# C5: hbar -> infinity (k -> 0) on the quartic: amplitude -> counting measure (all phases equal); the
#     stationary-phase window grows like sqrt(hbar), i.e. paths stop cancelling.
cvh = Bres["curv_hbar"][0]
widths = [math.sqrt(hb / abs(cvh)) for hb in [0.01, 1, 100]]
check("C5 window width sqrt(hbar/curv) grows without bound as hbar->inf (k->0): no cancelling, no classical path",
      widths[2] > 50 * widths[0], f"widths {np.round(widths,3)}")
R["C"] = {"Sp": Sp.tolist(), "gaps": gaps}

# ---------------------------------------------------------------------------
print("\n== D. Path sums as semirings: which aggregations are genuine, and is o_k a dial between them ==")


def is_semiring(add, mul, zero, one, sample, tol=1e-9):
    bad = []
    for a_, b_, c_ in itertools.product(sample, repeat=3):
        tests = {
            "add assoc": (add(add(a_, b_), c_), add(a_, add(b_, c_))),
            "add comm": (add(a_, b_), add(b_, a_)),
            "mul assoc": (mul(mul(a_, b_), c_), mul(a_, mul(b_, c_))),
            "left distrib": (mul(a_, add(b_, c_)), add(mul(a_, b_), mul(a_, c_))),
            "right distrib": (mul(add(b_, c_), a_), add(mul(b_, a_), mul(c_, a_))),
            "zero": (add(a_, zero), a_),
            "one": (mul(a_, one), a_),
            "zero absorbs": (mul(a_, zero), zero),
        }
        for nm, (u, v) in tests.items():
            if isinstance(u, float) and math.isinf(u) and u == v:
                continue
            if abs(u - v) > tol * (1 + abs(u) + abs(v)):
                bad.append(nm)
    return sorted(set(bad))


# D1: the o_k semiring.  oplus_k: a + b + 1/k (lift of +), otimes = o_k (lift of x). zero = wall -1/k, one = 0.
D = {"semiring_ok": {}}
for kv in [-2, -1, -0.3, 0.7, 1, 1j, 2j, 1 + 1j]:
    smp = [complex(z) for z in rng.normal(size=4) + 1j * rng.normal(size=4)] + [-1 / kv, 0]
    bad = is_semiring(lambda a_, b_, kv=kv: a_ + b_ + 1 / kv, lambda a_, b_, kv=kv: circ(a_, b_, kv),
                      -1 / kv, 0, smp)
    D["semiring_ok"][str(kv)] = not bad
check("D1 for every k tried (real and complex) (oplus = a+b+1/k, otimes = o_k) is a commutative semiring; zero is the wall -1/k, one is 0",
      all(D["semiring_ok"].values()), str(D["semiring_ok"]))
check("D2 the lift 1+k a carries oplus_k to + and o_k to x (so every finite k is (C,+,x) re-charted)",
      abs(lift(0.3 + 0.8 + 1 / 2.5, 2.5) - (lift(0.3, 2.5) + lift(0.8, 2.5))) < 1e-14)
# D3: o_k as the "or" (oplus) with ordinary x as "and": distributes only at k = 0 or on {0, 1}.
fails = {}
for kv in [-1, -0.5, 0.5, 1]:
    a_, b_, c_ = 0.6, 0.3, 0.45
    fails[kv] = abs(a_ * circ(b_, c_, kv) - circ(a_ * b_, a_ * c_, kv))
bool_ok = all(a_ * circ(b_, c_, -1) == circ(a_ * b_, a_ * c_, -1) for a_, b_, c_ in itertools.product([0, 1], repeat=3))
check("D3 the union rule o_{-1} as oplus with x as otimes fails distributivity off {0,1} (error = |k|a(1-a)bc), but on {0,1} it is OR/AND, the Boolean semiring",
      all(abs(fails[kv] - abs(kv) * 0.6 * 0.4 * 0.3 * 0.45) < 1e-14 for kv in fails) and bool_ok, f"{fails}")
D["union_fail"] = {str(kv): v for kv, v in fails.items()}
# D4: the classical catalogue (all genuine semirings).
INF = float("inf")
cat = {
    "counting / probability (R>=0, +, x)": (lambda a_, b_: a_ + b_, lambda a_, b_: a_ * b_, 0.0, 1.0, [0.0, 0.5, 1.3, 2.0]),
    "tropical min-plus (least action, shortest path)": (min, lambda a_, b_: a_ + b_, INF, 0.0, [0.0, 0.4, 1.7, 3.2, INF]),
    "max-plus": (max, lambda a_, b_: a_ + b_, -INF, 0.0, [0.0, -0.4, 1.7, 3.2, -INF]),
    "Viterbi max-times (most likely path)": (max, lambda a_, b_: a_ * b_, 0.0, 1.0, [0.0, 0.2, 0.7, 1.0]),
    "bottleneck max-min (widest path)": (max, min, -INF, INF, [-INF, 0.3, 1.0, 2.5, INF]),
    "Boolean OR-AND (reachability)": (lambda a_, b_: float(bool(a_) or bool(b_)), lambda a_, b_: float(bool(a_) and bool(b_)), 0.0, 1.0, [0.0, 1.0]),
    "log semiring (log-sum-exp, +) at h=0.3": (lambda a_, b_: -0.3 * np.logaddexp(-a_ / 0.3, -b_ / 0.3), lambda a_, b_: a_ + b_, INF, 0.0, [0.0, 0.4, 1.7, 3.0]),
}
D["catalogue"] = {}
for nm, (ad, mu, z0, o1, smp) in cat.items():
    bad = is_semiring(ad, mu, z0, o1, smp)
    D["catalogue"][nm] = bad
check("D4 counting, min-plus, max-plus, Viterbi, bottleneck, Boolean and the log semiring all pass every semiring axiom",
      all(not b for b in D["catalogue"].values()), str({k_: v for k_, v in D["catalogue"].items() if v}))
# D5: softening the second operation too (to reach bottleneck max-min) breaks distributivity at every finite softness.
errs = []
for tsoft in [1, 4, 16, 64]:
    oplus = lambda a_, b_: np.logaddexp(tsoft * a_, tsoft * b_) / tsoft          # soft max
    otimes = lambda a_, b_: -np.logaddexp(-tsoft * a_, -tsoft * b_) / tsoft      # soft min
    a_, b_, c_ = 0.9, 0.62, 0.6
    errs.append(abs(otimes(a_, oplus(b_, c_)) - oplus(otimes(a_, b_), otimes(a_, c_))))
check("D5 soft-max with soft-min (one dial toward bottleneck) is not distributive at finite softness t; the gap only closes as t -> infinity",
      all(e > 1e-9 for e in errs[:3]) and errs[0] > errs[1] > errs[2], f"errors {np.round(errs,6)}")
D["bottleneck_soft_err"] = errs
# D6: Maslov: in action coordinates the o_k semiring becomes (log-sum-exp_k, +); its limits are min-plus and max-plus.
sa, sb = 0.4, 1.3
for kv in [-1, -5, -20]:
    via_circ = math.log(lift(np.expm1(kv * sa) / kv + np.expm1(kv * sb) / kv + 1 / kv, kv)) / kv
    lse = np.logaddexp(kv * sa, kv * sb) / kv
    assert abs(via_circ - lse) < 1e-12
check("D6 pulled back to actions by a = E(ks)/k, oplus_k is log-sum-exp at temperature 1/|k| (Maslov's family); k->-inf gives min",
      abs(np.logaddexp(-1000 * sa, -1000 * sb) / -1000 - min(sa, sb)) < 1e-12)
# D7: provenance (Green, Karvounarakis, Tannen 2007): compute the path polynomial once in N[X], then any
#     aggregation is an evaluation of it.  Small graph s->a->t, s->b->t, s->t, a->b.
e = sp.symbols("e_sa e_at e_sb e_bt e_st e_ab")
esa, eat, esb, ebt, est, eab = e
paths = [[esa, eat], [esb, ebt], [est], [esa, eab, ebt]]
poly = sum(sp.Mul(*p_) for p_ in paths)
wS = {esa: 1.0, eat: 2.0, esb: 0.5, ebt: 2.0, est: 3.4, eab: 0.3}
cap = {esa: 5, eat: 2, esb: 3, ebt: 4, est: 1, eab: 6}
mn = min(sum(wS[x_] for x_ in p_) for p_ in paths)
wd = max(min(cap[x_] for x_ in p_) for p_ in paths)
cnt = poly.subs({x_: 1 for x_ in e})
check("D7 one provenance polynomial, many aggregates: count=4, shortest=" + str(mn) + ", widest=" + str(wd),
      cnt == 4 and abs(mn - 2.5) < 1e-12 and wd == 4)
D["provenance_poly"] = str(poly)
# D8: one dial on one graph: weights e^{k S}, read (1/k) log Z.  k=-inf shortest path, k=0 counts paths,
#     k imaginary gives the quantum amplitude.  Record the sweep.
Spath = np.array([sum(wS[x_] for x_ in p_) for p_ in paths])
ks = -np.logspace(-2, 3, 60)
D["sweep"] = {"k": ks.tolist(), "soft": [float(logsumexp(kv * Spath) / kv) for kv in ks],
              "Spath": Spath.tolist(),
              "quantum_abs": [float(abs(np.sum(np.exp(1j * (-kv) * Spath))) / len(Spath)) for kv in ks]}
check("D8 on one graph the dial slides from counting (|k| small: Z~4) to shortest path (|k| large: 2.5)",
      abs(np.sum(np.exp(-0.001 * Spath)) - 4) < 0.02 and abs(D["sweep"]["soft"][-1] - 2.5) < 1e-3)
R["D"] = D

# ---------------------------------------------------------------------------
print("\n== E. Zooming out: the height flow is a critical branching population ==")
# E1: E(x) = f(1+x) - 1 with f(y) = e^{y-1}, the Poisson(1) offspring generating function (critical: mean 1).
#     So E^n(-1) = f_n(0) - 1 = -(probability a Galton-Watson family line survives n generations).
q = 0.0
surv = []
xh = -1.0
for n in range(1, 2001):
    q = math.exp(q - 1)
    xh = math.expm1(xh)
    surv.append(1 - q)
    assert abs((q - 1) - xh) < 1e-12
check("E1 E^n(-1) = -(survival to generation n) of a critical Poisson(1) Galton-Watson tree, n=1..2000", True)
check("E2 Kolmogorov: survival * n -> 2/sigma^2 = 2 (sigma^2=1): the height flow's 'fades like 2/height'",
      abs(surv[1999] * 2000 - 2) < 0.01, f"n*P(survive) at n=2000: {surv[1999]*2000:.4f}")
# Monte Carlo
pois = np.random.default_rng(7)
trials, gens = 40000, 30
alive = np.ones(trials, dtype=np.int64)
mc = []
for g_ in range(gens):
    alive = np.array([pois.poisson(1.0, size=m).sum() if m > 0 else 0 for m in alive]) if g_ < 12 else \
        pois.poisson(alive.astype(float))
    mc.append(float(np.mean(alive > 0)))
check("E3 simulation of 40,000 family trees matches -E^n(-1) for n=1..30 within 3 sigma",
      all(abs(mc[i] - surv[i]) < 3 * math.sqrt(surv[i] * (1 - surv[i]) / trials) + 1e-3 for i in range(gens)))
R["E"] = {"surv": surv[:200], "mc": mc}
# E4: Bell numbers (already on the combinatorics page): E twice = set of nonempty sets of a set.
xs = sp.symbols("x")
E2 = sp.series(sp.exp(sp.exp(xs) - 1) - 1, xs, 0, 8).removeO()
bells = [sp.factorial(n) * E2.coeff(xs, n) for n in range(1, 8)]
check("E4 E twice counts groupings of a set (Bell numbers 1,2,5,15,52,203,877): 'a set of sets' as an aggregate",
      [int(b) for b in bells] == [1, 2, 5, 15, 52, 203, 877])
# E5: Ising chain decimation.  Bonds in series multiply t = tanh K; bonds in parallel add K, so t's combine by
#     velocity addition (a+b)/(1+ab).  Check against exact partition functions.
def Z_chain(Ks):
    Tm = [np.array([[math.exp(K_), math.exp(-K_)], [math.exp(-K_), math.exp(K_)]]) for K_ in Ks]
    M = np.eye(2)
    for t_ in Tm:
        M = M @ t_
    return M  # two-end transfer matrix
K1, K2, K3 = 0.4, 0.9, 0.25
M = Z_chain([K1, K2])
t_eff = (M[0, 0] - M[0, 1]) / (M[0, 0] + M[0, 1])
check("E5 series bonds: t_eff = tanh K1 * tanh K2 (the lift's multiplication)",
      abs(t_eff - math.tanh(K1) * math.tanh(K2)) < 1e-14)
check("E6 parallel bonds: t = (t1+t2)/(1+t1 t2), relativistic velocity addition, exactly",
      abs(math.tanh(K1 + K2) - (math.tanh(K1) + math.tanh(K2)) / (1 + math.tanh(K1) * math.tanh(K2))) < 1e-15)
tt = [math.tanh(K_) for K_ in (K1, K2, K3)]
vadd = lambda a_, b_: (a_ + b_) / (1 + a_ * b_)
derr = abs(tt[0] * vadd(tt[1], tt[2]) - vadd(tt[0] * tt[1], tt[0] * tt[2]))
check("E7 series-over-parallel does not distribute (different networks): series/parallel bonds are not a semiring",
      derr > 1e-3, f"gap {derr:.4f}")
# E8: the 1-D chain's own RG flow t -> t^b is linear in u=-ln t (relevant), not E's marginal 2/height fade.
check("E8 1-D Ising RG: u=-ln tanh K doubles per 2x zoom (u' = 2u): linear, not marginal like E",
      abs(-math.log(math.tanh(0.7) ** 2) - 2 * (-math.log(math.tanh(0.7)))) < 1e-14)
# E9: fractional heights: E^(1/2)'s formal series (Baker 1958: radius 0) has coefficients of both signs.
cH = [sp.Rational(1, 1)]
# coefficients of formal half-iterate h(h(x)) = e^x - 1, h = x + sum c_n x^n
cs_ = sp.symbols("c2:9")
h = xs + sum(cs_[i] * xs ** (i + 2) for i in range(7))
eq = sp.expand(sp.series(h.subs(xs, h) - (sp.exp(xs) - 1), xs, 0, 9).removeO())
solh = {}
for n in range(2, 9):
    cn = eq.coeff(xs, n).subs(solh)
    solh[cs_[n - 2]] = sp.solve(cn, cs_[n - 2])[0]
half = [solh[c_] for c_ in cs_]
check("E9 the half-height flow E^(1/2) has coefficients of both signs (1/4, 1/48, 0, 1/3840, -7/92160, ...)",
      any(v < 0 for v in half) and half[0] == sp.Rational(1, 4), str(half[:6]))
R["E"]["half"] = [str(v) for v in half]

# ---------------------------------------------------------------------------
print("\n== F. Decoherence: which-path records compose by the union rule ==")
# Two paths, M environment qubits; path 1 leaves qubit i in |0>, path 2 rotates it by angle th_i.
# Overlap c_i = cos th_i.  Fringe visibility V = prod |c_i|.  Distinguishability 1 - V^2 = o_{-1} of p_i = sin^2 th_i.
th = rng.uniform(0.05, 0.6, size=6)
pvals = np.sin(th) ** 2
un = 0.0
for p_ in pvals:
    un = circ(un, p_, -1)
Vis = np.prod(np.cos(th))
check("F1 1 - visibility^2 = p1 o_{-1} p2 o_{-1} ... (union rule over the environment), exact", abs((1 - Vis ** 2) - un) < 1e-14)
# brute force: full state of particle (2 paths) x 6 qubits, reduced density matrix
def kron_all(vs):
    out = np.array([1.0 + 0j])
    for v_ in vs:
        out = np.kron(out, v_)
    return out
env1 = kron_all([np.array([1, 0], complex)] * 6)
env2 = kron_all([np.array([math.cos(t_), math.sin(t_)], complex) for t_ in th])
rho01 = 0.5 * np.vdot(env2, env1)
check("F2 brute-force reduced density matrix: 2|rho_12| equals prod cos(theta_i) (64-dim environment)",
      abs(2 * abs(rho01) - Vis) < 1e-14)
# F3: law of large numbers: many weak records kill the fringes: V ~ exp(-sum th^2 /2)
Ms = [0, 2, 8, 32, 128]
thw = 0.15
vis_M = [math.cos(thw) ** m for m in Ms]
check("F3 many weak records: visibility = cos(0.15)^M ~ exp(-M*0.0113) -> 0; at M=128 it is <0.25",
      vis_M[-1] < 0.25 and abs(vis_M[-1] - math.exp(-128 * thw ** 2 / 2)) < 0.02)
xsc = np.linspace(-3, 3, 400)
R["F"] = {"M": Ms, "vis": vis_M, "x": xsc.tolist(),
          "fringes": {str(m): (1 + v_ * np.cos(6 * xsc)).tolist() for m, v_ in zip(Ms, vis_M)}}
a1, a2 = 0.5 * cmath.exp(0.3j), 0.4 * cmath.exp(2.1j)
Pdec = lambda v_: abs(a1) ** 2 + abs(a2) ** 2 + 2 * v_ * (a1 * a2.conjugate()).real
check("F4 P = |a1|^2+|a2|^2+2V Re(a1 a2*): V=1 is |a1+a2|^2 (amplitudes add), V=0 is P1+P2 (probabilities add)",
      abs(Pdec(1) - abs(a1 + a2) ** 2) < 1e-15 and abs(Pdec(0) - (abs(a1) ** 2 + abs(a2) ** 2)) < 1e-15)

print(f"\n{passed} passed, {failed} failed")
R["summary"] = {"passed": passed, "failed": failed}
json.dump(R, open("results.json", "w"), default=lambda o: float(o) if isinstance(o, (np.floating,)) else str(o))
