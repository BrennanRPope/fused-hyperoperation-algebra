"""DRAFT round 23 checks: a lump of matter in the law's de Sitter spacetime.

Round 20 built de Sitter spacetime from the law (one dial q, a dot product with one minus sign); round 22 put light and free
matter into it, with no pull back on the geometry. Here matter is allowed to curve the law back: a point mass gives
Schwarzschild-de Sitter, and a ball of dust is the law's own slices with matter in them, glued to that outside. The question is
how much of the bent geometry is still the law's.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, time
import numpy as np
import sympy as spy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

# ------------------------------------------------------------------ the law with a Lorentzian dot product (rounds 20-22, unchanged)
ETA = np.array([-1.0, 1, 1, 1])
def dot4(a, b):
    return (a * b * ETA).sum(-1, keepdims=True)
def place(x, qq=1.0):             # round 20's placement into 5D (W, T, X, Y, Z), on W^2 - T^2 + X^2 + Y^2 + Z^2 = 1/q
    e = dot4(x, x); return np.concatenate([(1 - qq * e), 2 * math.sqrt(qq) * x], -1) / (math.sqrt(qq) * (1 + qq * e))

# ------------------------------------------------------------------ curvature of a diagonal metric (symbolic), used by E1, E4, E5
def curvature(gdiag, X):
    n = len(X); g = spy.diag(*gdiag); gi = spy.diag(*[1 / v for v in gdiag])
    Gam = [[[spy.simplify(sum(gi[a, d] * (spy.diff(g[d, b], X[c]) + spy.diff(g[d, c], X[b]) - spy.diff(g[b, c], X[d])) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    R = [[[[spy.simplify(spy.diff(Gam[a][b][d], X[c]) - spy.diff(Gam[a][b][c], X[d])
                         + sum(Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c] for e in range(n)))
            for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]          # R^a_{bcd}
    Ric = spy.Matrix(n, n, lambda b, d: spy.simplify(sum(R[a][b][a][d] for a in range(n))))
    Rs = spy.simplify(sum(gi[b, b] * Ric[b, b] for b in range(n)))
    Rl = [[[[g[a, a] * R[a][b][c][d] for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]   # R_{abcd}
    return g, gi, Ric, Rs, Rl
def weyl_sq(g, gi, Ric, Rs, Rl, n=4):
    def C(a, b, c, d):
        return (Rl[a][b][c][d] - (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c] - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c]) / (n - 2)
                + Rs * (g[a, c] * g[b, d] - g[a, d] * g[b, c]) / ((n - 1) * (n - 2)))
    W2 = 0; K2 = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    f_ = gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d]
                    W2 += f_ * C(a, b, c, d)**2; K2 += f_ * Rl[a][b][c][d]**2
    return spy.simplify(W2), spy.simplify(K2)

# ================================================================== E. A lump of matter in the law's de Sitter
say("== E. Matter curving the law's spacetime back")
t, r, th, ph, M, q = spy.symbols('t r theta phi M q', positive=True)
f_s = 1 - 2 * M / r - q * r**2
gS, giS, RicS, RsS, RlS = curvature([-f_s, 1 / f_s, r**2, r**2 * spy.sin(th)**2], [t, r, th, ph])
ein_ok = all(spy.simplify(RicS[i, j] - 3 * q * gS[i, j]) == 0 for i in range(4) for j in range(4))
W2S, K2S = weyl_sq(gS, giS, RicS, RsS, RlS)
k_ok = spy.simplify(K2S - (24 * q**2 + 48 * M**2 / r**6)) == 0
w_ok = spy.simplify(W2S - 48 * M**2 / r**6) == 0
check("E1 (exact, symbolic; textbook) a point mass M in the law's de Sitter: f(r) = 1 - 2M/r - q r^2, metric -f dt^2 + dr^2/f + r^2 dOmega^2, has Ricci = 3q g for every M. So outside the mass it solves Einstein's equations with the law's Lambda = 3q and no matter: Schwarzschild-de Sitter (Kottler 1918). Its total curvature (Kretschmann) is 24 q^2 + 48 M^2/r^6: the 24 q^2 is the law's de Sitter, unchanged by the mass, and the whole of the mass's part is Weyl (tidal) curvature, 48 M^2/r^6",
      ein_ok and k_ok and w_ok, f"Ricci - 3q g = 0: {ein_ok}; Kretschmann = 24q^2 + 48M^2/r^6: {k_ok}; Weyl^2 = 48M^2/r^6: {w_ok}; scalar curvature {spy.simplify(RsS)}")
results["kottler"] = dict(einstein_lambda_3q=bool(ein_ok), kretschmann=str(spy.simplify(K2S)), weyl_sq=str(spy.simplify(W2S)), scalar=str(spy.simplify(RsS)))

# E2: at M = 0 this is the law's de Sitter, by an explicit change of coordinates from the law's chart
rg = np.random.default_rng(23)
def static_coords(x, qq):          # from the law's chart: place in 5D, then static time t_s and radius r (W > |T|: the pole's static patch)
    P = place(x, qq); W, T = P[:, 0], P[:, 1]; Xs = P[:, 2:]
    return np.arctanh(T / W) / math.sqrt(qq), np.linalg.norm(Xs, axis=1), Xs
def law_metric(x, qq):
    e = dot4(x, x)[:, 0]; return 4 / (1 + qq * e)**2
rowsE2 = []
for qq in (1.0, 0.37, 2.5):
    worst = 0.0; n_ok = 0
    while n_ok < 200:
        x = rg.normal(size=(1, 4)) * 0.3 / math.sqrt(qq)
        P = place(x, qq)
        if not (P[0, 0] > abs(P[0, 1]) + 0.05 and abs(1 + qq * dot4(x, x)[0, 0]) > 0.3):
            continue
        n_ok += 1; h = 1e-6
        def Y(xx):                                  # static chart (t_s, cartesian X) so the static metric is -f dt^2 + spatial part
            ts, rr, Xs = static_coords(xx, qq); return np.concatenate([ts[:, None], Xs], 1)
        J = np.array([(Y(x + h * np.eye(4)[k]) - Y(x - h * np.eye(4)[k]))[0] / (2 * h) for k in range(4)]).T   # dY/dx
        ts, rr, Xs = static_coords(x, qq); fr = 1 - qq * rr[0]**2; nX = Xs[0] / rr[0]
        G = np.zeros((4, 4)); G[0, 0] = -fr; G[1:, 1:] = np.eye(3) + np.outer(nX, nX) * (1 / fr - 1)   # dr^2/f + r^2 dOmega^2
        pull = J.T @ G @ J; law = law_metric(x, qq)[0] * np.diag(ETA)
        worst = max(worst, np.abs(pull - law).max() / np.abs(law).max())
    rowsE2.append(dict(q=qq, worst=float(worst)))
okE2 = all(rw["worst"] < 1e-7 for rw in rowsE2)
check("E2 (exact; textbook coordinates) with no mass, Schwarzschild-de Sitter is the law's own spacetime: placing the law's chart in 5D (round 20) and reading static time t_s = artanh(T/W)/sqrt(q) and radius r = |(X, Y, Z)| turns the law's metric 4 eta(dx, dx)/(1 + q eta)^2 into -(1 - q r^2) dt_s^2 + dr^2/(1 - q r^2) + r^2 dOmega^2. Its horizon, r = 1/sqrt(q), is the pole's view of round 22's horizon. So the mass M is the one new number; everything else is the law",
      okE2, "; ".join(f"q = {rw['q']}: 200 points of the pole's static patch, pulled-back metric vs the law's {rw['worst']:.0e} (relative)" for rw in rowsE2))
results["static_patch"] = rowsE2

# E3: horizons: a black hole inside the law's cosmic horizon, up to the Nariai limit
def horizons(Mv, qq):
    rr = np.roots([-qq, 0, 1, -2 * Mv]); rr = np.sort(rr[np.abs(rr.imag) < 1e-12].real); return rr[rr > 0]
rowsE3 = []
for qq in (1.0, 0.25):
    Mn = 1 / math.sqrt(27 * qq)
    for frac in (0.0001, 0.3, 0.9, 0.999, 1.001, 1.5):
        hs = horizons(frac * Mn, qq); rowsE3.append(dict(q=qq, M_over_Mnariai=frac, horizons=hs.tolist()))
ok3 = True
for rw in rowsE3:
    two = len(rw["horizons"]) == 2
    ok3 &= (two if rw["M_over_Mnariai"] < 1 else not two)
    if rw["M_over_Mnariai"] == 0.0001:
        ok3 &= abs(rw["horizons"][0] - 2 * 0.0001 / math.sqrt(27 * rw["q"])) < 1e-6 and abs(rw["horizons"][1] - 1 / math.sqrt(rw["q"])) < 1e-4
    if rw["M_over_Mnariai"] == 0.999:
        ok3 &= abs(rw["horizons"][0] - 1 / math.sqrt(3 * rw["q"])) < 0.05 / math.sqrt(rw["q"]) and abs(rw["horizons"][1] - 1 / math.sqrt(3 * rw["q"])) < 0.05 / math.sqrt(rw["q"])
check("E3 (exact; textbook) a mass makes two horizons: a black hole's (about 2M for small M) inside the law's cosmic one (about 1/sqrt(q)). They exist together only while 27 q M^2 < 1 and merge at r = 1/sqrt(3q) at the Nariai limit M = 1/sqrt(27 q) (Nariai 1950). So the law's dial sets the largest black hole its spacetime can hold",
      ok3, "; ".join(f"q {rw['q']}, M/M_N {rw['M_over_Mnariai']}: horizons {[round(h, 4) for h in rw['horizons']]}" for rw in rowsE3 if rw["q"] == 1.0))
results["horizons"] = rowsE3

# E4: a ball of dust: inside, the law's slices with matter; outside, Schwarzschild-de Sitter; the edge falls freely
tt_, chi_, a_ = spy.symbols('t chi a', positive=True); A = spy.Function('a')(tt_); rho = spy.Function('rho')(tt_)
gF, giF, RicF, RsF, RlF = curvature([-1, A**2, A**2 * spy.sin(chi_)**2, A**2 * spy.sin(chi_)**2 * spy.sin(th)**2], [tt_, chi_, th, ph])
Gtt = spy.simplify(RicF[0, 0] - RsF * gF[0, 0] / 2); Gcc = spy.simplify(RicF[1, 1] - RsF * gF[1, 1] / 2)
fried = spy.simplify(Gtt - 3 * q - 8 * spy.pi * rho)                         # G_tt + Lambda g_tt = 8 pi rho, with g_tt = -1
fried_target = spy.simplify(3 * (spy.diff(A, tt_)**2 + 1) / A**2 - 3 * q - 8 * spy.pi * rho)
f1_ok = spy.simplify(fried - fried_target) == 0
# the space-space equation with rho a^3 = const (dust) follows from the time derivative of the Friedmann equation: check on solutions
W2F, K2F = weyl_sq(gF, giF, RicF, RsF, RlF)
wF_ok = spy.simplify(W2F) == 0
say(f"   Friedmann from Einstein on the law's slices with dust: G_tt + Lambda g_tt - 8 pi rho = 3(a'^2 + 1)/a^2 - 3q - 8 pi rho: {f1_ok}; Weyl^2 inside = {spy.simplify(W2F)}")
rowsE4 = []
for name, qq, a0, chi0 in (("collapses to a black hole", 1.0, 0.45, 0.6), ("expands with the law", 1.0, 0.95, 0.6)):
    C = a0 * (1 - qq * a0**2)                            # turning point at a0: a'^2 = C/a + q a^2 - 1 = 0
    rho0 = 3 * C / (8 * math.pi * a0**3); q_in = 1 / a0**2
    Mball = C * math.sin(chi0)**3 / 2; E = math.cos(chi0)
    # inside: second-order Friedmann a'' = -C/(2a^2) + q a, from rest at a0, integrated in the dust's proper time
    def rin(s, y): return [y[1], -C / (2 * y[0]**2) + qq * y[0]]
    def rout(s, y): return [y[1], -(2 * Mball / y[0]**2 - 2 * qq * y[0]) / 2]          # radial geodesic in Kottler: R'' = -f'(R)/2
    stop = lambda s, y: y[0] - 1e-3; stop.terminal = True
    big = lambda s, y: y[0] - 30; big.terminal = True
    s_end = 6.0
    si = solve_ivp(rin, (0, s_end), [a0, 0.0], rtol=1e-11, atol=1e-13, events=[stop, big], dense_output=True, method="DOP853")
    so = solve_ivp(rout, (0, s_end), [a0 * math.sin(chi0), 0.0], rtol=1e-11, atol=1e-13, events=[stop, big], dense_output=True, method="DOP853")
    s_common = min(si.t[-1], so.t[-1]); ss = np.linspace(0, s_common, 400)
    Rin = si.sol(ss)[0] * math.sin(chi0); Rout = so.sol(ss)[0]; Rdot = so.sol(ss)[1]
    match = float(np.abs(Rin - Rout).max() / np.abs(Rout).max())
    energy = float(np.abs(Rdot**2 + (1 - 2 * Mball / Rout - qq * Rout**2) - E**2).max())
    hs = horizons(Mball, qq); ham = abs(q_in - (qq + 8 * math.pi * rho0 / 3))
    inside_bh = bool(len(hs) == 2 and Rout[-1] < hs[0])
    rowsE4.append(dict(case=name, q=qq, a0=a0, chi0=chi0, C=C, rho0=rho0, q_inside=q_in, M=Mball, E=E, horizons=hs.tolist(), match=match, energy=energy,
                       hamiltonian=ham, end_R=float(Rout[-1]), end_s=float(s_common), crossed_bh_horizon=inside_bh,
                       s=ss[::8].tolist(), R=Rout[::8].tolist(), collapse=bool(si.t[-1] < s_end and si.y[0, -1] < 0.01)))
okE4 = f1_ok and wF_ok and all(rw["match"] < 1e-8 and rw["energy"] < 1e-8 and rw["hamiltonian"] < 1e-12 for rw in rowsE4) and rowsE4[0]["crossed_bh_horizon"] and rowsE4[0]["collapse"] and rowsE4[1]["end_R"] > 10 * rowsE4[1]["horizons"][1]
check("E4 (exact; textbook) a ball of dust: inside, the slices stay the law's round 3-spheres, with Einstein giving a'^2 = C/a + q a^2 - 1 (closed Friedmann with dust and Lambda = 3q) and no Weyl curvature. At its moment of rest each inside slice is the law's 3-sphere at a larger dial, q_in = q + 8 pi rho/3. The ball's edge, at angle chi0, falls exactly as a free particle of energy cos(chi0) in Schwarzschild-de Sitter of mass M = C sin^3(chi0)/2 (Oppenheimer-Snyder with Lambda). A dense ball collapses through the black-hole horizon; a light one is carried off by the law's expansion",
      okE4, "; ".join(f"{rw['case']} (a0 {rw['a0']}, chi0 {rw['chi0']}): q_in {rw['q_inside']:.3f} = q + 8 pi rho/3 to {rw['hamiltonian']:.0e}, M {rw['M']:.4f}, horizons {[round(h, 4) for h in rw['horizons']]}, edge inside vs Kottler geodesic {rw['match']:.0e}, energy {rw['energy']:.0e}, R at s = {rw['end_s']:.2f}: {rw['end_R']:.3f}" for rw in rowsE4))
results["dust_ball"] = dict(friedmann=bool(f1_ok), weyl_inside=str(spy.simplify(W2F)), cases=rowsE4)

# E5: what the law cannot carry: tidal curvature. The law's spacetime and its slices have none; the mass's outside is all tidal
x0, x1, x2, x3 = spy.symbols('x0 x1 x2 x3', real=True)
Om2 = 4 / (1 + q * (-x0**2 + x1**2 + x2**2 + x3**2))**2
gL, giL, RicL, RsL, RlL = curvature([-Om2, Om2, Om2, Om2], [x0, x1, x2, x3])
W2L, K2L = weyl_sq(gL, giL, RicL, RsL, RlL)
wL_ok = spy.simplify(W2L) == 0
# time-symmetric slice (t = const) of Kottler: h = dr^2/f + r^2 dOmega^2; its two sectional curvatures
Kt = spy.simplify((1 - f_s) / r**2); Krt = spy.simplify(-spy.diff(f_s, r) / (2 * r))
scal3 = spy.simplify(2 * (Kt + 2 * Krt))
sec_ok = spy.simplify(Kt - (q + 2 * M / r**3)) == 0 and spy.simplify(Krt - (q - M / r**3)) == 0 and spy.simplify(scal3 - 6 * q) == 0
# the law's 3D slice at the dial q (round 13's placement) has every sectional curvature = q: check K for its metric 4/(1+q r^2)^2 (dr^2 + r^2 dOmega^2)
psi = 4 / (1 + q * r**2)**2; areal = spy.sqrt(psi) * r              # areal radius R(r); slice is dR^2/F(R) + R^2 dOmega^2 with F = (dR/dr)^2/psi
F_law = spy.simplify(spy.diff(areal, r)**2 / psi)
Kt_law = spy.simplify((1 - F_law) / areal**2); Krt_law = spy.simplify(-spy.diff(F_law, r) / spy.diff(areal, r) / (2 * areal))
law_sec_ok = spy.simplify(Kt_law - q) == 0 and spy.simplify(Krt_law - q) == 0
check("E5 (exact, symbolic) where the law stops: its spacetime has no Weyl (tidal) curvature (Weyl^2 = 0 for 4 eta/(1 + q eta)^2; conformally flat, textbook for de Sitter), and neither do the dust ball's insides (E4). On a slice at rest, the law's 3-sphere bends the same way in every plane (sectional curvature q), while the slice outside a mass bends q + 2M/r^3 across the sphere and q - M/r^3 along the radius. The average is still the law's q (scalar curvature 6q, the constraint), but the split, which falls as M/r^3 and is what pulls an object apart and in, is not the law's. The law gives the Ricci part (matter and Lambda); the pull of a lump on the space around it is the part it lacks",
      wL_ok and sec_ok and law_sec_ok, f"Weyl^2 of the law's spacetime = {spy.simplify(W2L)}; Kottler slice: across {Kt}, along {Krt}, scalar {scal3}; the law's slice: across {Kt_law}, along {Krt_law}")
results["tidal"] = dict(weyl_law=str(spy.simplify(W2L)), kottler_slice_across=str(Kt), kottler_slice_along=str(Krt), kottler_slice_scalar=str(scal3),
                        law_slice_across=str(Kt_law), law_slice_along=str(Krt_law))

say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed, seconds=round(time.time() - T_START, 1))
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
