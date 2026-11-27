"""DRAFT round 20 checks: de Sitter spacetime from the law, with space and time in one law.

Rounds 18-19 found the law's two signs are the two groups: at q > 0 the 3-sphere's rotations, at q < 0 the Lorentz boosts.
Round 16 needed two dials of opposite sign (velocity and space, k^2 = -q_v q_s). Here the law is kept exactly as it is, with one
dial q, and only its dot product is made Lorentzian, a.b = -a0 b0 + a1 b1 + a2 b2 + a3 b3, on 4D points x = (t, x, y, z).
Run: python3 checks.py  (numpy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, time
import numpy as np
import sympy as spy

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

ETA = np.array([-1.0, 1, 1, 1])
def dot(a, b, eta=ETA):
    return (a * b * eta).sum(-1, keepdims=True)
def oplus(a, b, qq, eta=ETA):     # the law of rounds 12-19, word for word, with its dot product as a parameter
    ab = dot(a, b, eta); aa = dot(a, a, eta); bb = dot(b, b, eta)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def place(x, qq=1.0):             # inverse stereographic placement (rounds 13, 18, 19) into 5D: (W, T, X, Y, Z)
    e = dot(x, x); return np.concatenate([(1 - qq * e), 2 * math.sqrt(abs(qq)) * x], -1) / (math.sqrt(abs(qq)) * (1 + qq * e))
G5 = np.diag([1.0, -1, 1, 1, 1])  # ambient form W^2 - T^2 + X^2 + Y^2 + Z^2
amb = lambda P: (P * P * np.diag(G5)).sum(-1)
rg = np.random.default_rng(20)
def pts(n, s=0.4):                # random chart points kept away from the chart's edge (where 1 +/- eta(x, x) vanishes)
    out = np.zeros((0, 4))
    while len(out) < n:
        c = rg.normal(size=(4 * n, 4)) * s; e = dot(c, c)[:, 0]
        out = np.vstack([out, c[(np.abs(e) < 0.5)]])
    return out[:n]
def fit_linear(A, B):
    L = np.linalg.lstsq(A, B, rcond=None)[0].T; return L, np.abs(A @ L.T - B).max()

# ================================================================== D. de Sitter from the law
say("== D. The law with a Lorentzian dot product")
Ys = pts(400); onds = np.abs(amb(place(Ys)) - 1).max()
rows1 = []
for x in pts(30):
    L, res = fit_linear(place(Ys), place(oplus(np.tile(x, (400, 1)), Ys, 1.0)))
    rows1.append((res, np.abs(L.T @ G5 @ L - G5).max(), abs(np.linalg.det(L) - 1), np.abs(L[:, 0] - place(x[None])[0]).max()))
r1 = np.array(rows1).max(0)
check("D1 (exact) with a Lorentzian dot product and q = +1, the law's points lie on 4D de Sitter space (W^2 - T^2 + X^2 + Y^2 + Z^2 = 1), and adding x is a linear map of the 5D space that keeps that form, has determinant 1 and carries rest to x: a de Sitter isometry in SO(1,4). So the law is de Sitter's own move from rest to x, as at q = +1 with a plain dot it is the 3-sphere's (round 18)",
      onds < 1e-12 and r1[0] < 1e-10 and r1[1] < 1e-10 and r1[2] < 1e-10 and r1[3] < 1e-12,
      f"400 points on de Sitter to {onds:.0e}; 30 random x: worst linear fit {r1[0]:.0e}, form kept {r1[1]:.0e}, det - 1 {r1[2]:.0e}, rest -> x {r1[3]:.0e}")
results["isometry"] = dict(on_ds=float(onds), linear=float(r1[0]), form=float(r1[1]), det=float(r1[2]), rest_to_x=float(r1[3]))

eta4 = np.diag(ETA); rows2 = []; boost_part = 0.0
Zs = pts(300)
for x, y in zip(pts(20), pts(20)):
    xy = oplus(x[None], y[None], 1.0)
    lhs = oplus(np.tile(x, (300, 1)), oplus(np.tile(y, (300, 1)), Zs, 1.0), 1.0)
    # find the 4x4 turn U with x (+) (y (+) z) = (x (+) y) (+) U z : solve by undoing (x (+) y) on the left, then fit
    Uz = oplus(-np.tile(xy[0], (300, 1)), lhs, 1.0)
    U, res = fit_linear(Zs, Uz)
    rows2.append((res, np.abs(U.T @ eta4 @ U - eta4).max(), abs(np.linalg.det(U) - 1)))
    boost_part = max(boost_part, np.abs(U[0, 1:]).max())
r2 = np.array(rows2).max(0)
check("D2 (exact) the law's twist is now a Lorentz transformation of the local frame: the turn U by which the law fails to associate is linear, keeps -t^2 + x^2 + y^2 + z^2 and has determinant 1, and in general it mixes time with space (a boost), not only space with space. So the local frames of the law's de Sitter space carry the Lorentz group of round 19, as twists of the same law",
      r2[0] < 1e-10 and r2[1] < 1e-10 and r2[2] < 1e-10 and boost_part > 0.01,
      f"20 random pairs, 300 test points each: worst linear fit {r2[0]:.0e}, Minkowski form kept {r2[1]:.0e}, det - 1 {r2[2]:.0e}; largest time-space entry of U {boost_part:.3f}")
results["twist"] = dict(linear=float(r2[0]), form=float(r2[1]), det=float(r2[2]), boost_part=float(boost_part))

# D3: the three kinds of direction
X3, Y3 = rg.normal(size=(2000, 3)) * 0.8, rg.normal(size=(2000, 3)) * 0.8
pad = lambda v: np.concatenate([np.zeros((len(v), 1)), v], 1)
e_space = np.abs(oplus(pad(X3), pad(Y3), 1.0)[:, 1:] - oplus(X3, Y3, 1.0, np.ones(3))).max()
e_space_t = np.abs(oplus(pad(X3), pad(Y3), 1.0)[:, 0]).max()
a, b = rg.uniform(-0.9, 0.9, 2000), rg.uniform(-0.9, 0.9, 2000)
tline = oplus(np.stack([a, 0 * a, 0 * a, 0 * a], 1), np.stack([b, 0 * b, 0 * b, 0 * b], 1), 1.0)[:, 0]
e_time = np.abs(tline - (a + b) / (1 + a * b)).max()
e_time_qm = np.abs(tline - oplus(a[:, None], b[:, None], -1.0, np.ones(1))[:, 0]).max()
n = rg.normal(size=(2000, 3)); n /= np.linalg.norm(n, axis=1)[:, None]; null = np.concatenate([np.ones((2000, 1)), n], 1)
s1, s2 = rg.uniform(-3, 3, (2000, 1)), rg.uniform(-3, 3, (2000, 1))
e_null = np.abs(oplus(s1 * null, s2 * null, 1.0) - (s1 + s2) * null).max()
check("D3 (exact) one dial, three behaviours. Along space (t = 0) the law is round 12's 3D law at +q. Along time it is the law at -q: (a + b)/(1 + ab), rapidity addition, the Lorentz side of round 19. Along a light ray (a null direction) it is plain addition, untouched by the dial. So round 16's two dials of opposite sign come from one dial and the signature of the dot product",
      e_space < 1e-12 and e_space_t == 0 and e_time < 1e-12 and e_time_qm < 1e-12 and e_null < 1e-12,
      f"2000 pairs each: space vs round 12's law {e_space:.0e} (stays in space: {e_space_t:.0e}); time vs (a+b)/(1+ab) {e_time:.0e}, vs the law at q = -1 {e_time_qm:.0e}; light rays vs plain addition {e_null:.0e}")
results["directions"] = dict(space=float(e_space), time=float(e_time), time_vs_q_minus=float(e_time_qm), null=float(e_null))

# D4: curvature of the law's metric, symbolically.  Nearby points: the law's distance gives ds^2 = 4 eta(dx,dx)/(1 + q eta(x,x))^2.
qs_ = spy.symbols('q', positive=True); xs = spy.symbols('t x y z', real=True); et = [-1, 1, 1, 1]
E = sum(e * v**2 for e, v in zip(et, xs)); f = 4 / (1 + qs_ * E)**2
g = spy.diag(*[f * e for e in et]); gi = g.inv()
Gam = [[[sum(gi[a_, d] * (spy.diff(g[d, b_], xs[c]) + spy.diff(g[d, c], xs[b_]) - spy.diff(g[b_, c], xs[d])) for d in range(4)) / 2
         for c in range(4)] for b_ in range(4)] for a_ in range(4)]
def Riem(a_, b_, c, d):
    return (spy.diff(Gam[a_][b_][d], xs[c]) - spy.diff(Gam[a_][b_][c], xs[d])
            + sum(Gam[a_][c][e] * Gam[e][b_][d] - Gam[a_][d][e] * Gam[e][b_][c] for e in range(4)))
Ric = spy.Matrix(4, 4, lambda b_, d: spy.simplify(sum(Riem(a_, b_, a_, d) for a_ in range(4))))
ric_ok = spy.simplify(Ric - 3 * qs_ * g) == spy.zeros(4, 4)
Rs = spy.simplify(sum(gi[i, j] * Ric[i, j] for i in range(4) for j in range(4)))
Gt = spy.simplify(Ric - Rs / 2 * g)
Lam = spy.simplify(-Gt[1, 1] / g[1, 1])
ein_ok = spy.simplify(Gt + Lam * g) == spy.zeros(4, 4)
# numerical check that the law's own distance gives that metric, and that the 5D placement pulls back to it
x0 = pts(1, 0.3)[0]; hstep = 1e-6; dvs = rg.normal(size=(6, 4)) * hstep
law_d2 = [float((dot(oplus(-x0[None], (x0 + dv)[None], 1.0), oplus(-x0[None], (x0 + dv)[None], 1.0))[0, 0])) * 4 for dv in dvs]
met_d2 = [float(4 * dot(dv[None], dv[None])[0, 0] / (1 + dot(x0[None], x0[None])[0, 0])**2) for dv in dvs]
pull = [float(((place((x0 + dv)[None]) - place(x0[None]))**2 * np.diag(G5)).sum()) for dv in dvs]
e_met = max(abs(l - m) / abs(m) for l, m in zip(law_d2, met_d2)); e_pull = max(abs(p - m) / abs(m) for p, m in zip(pull, met_d2))
check("D4 (exact, symbolic) the law's own metric ds^2 = 4 eta(dx, dx)/(1 + q eta(x, x))^2 has Ricci = 3q g, so it solves Einstein's equations with no matter and Lambda = 3q: de Sitter space of radius 1/sqrt(q). This is the same Lambda = 3q that round 16 read from the glued slabs at a moment of rest. The law's distance for nearby points and the 5D placement both give that metric",
      ric_ok and ein_ok and spy.simplify(Lam - 3 * qs_) == 0 and e_met < 1e-4 and e_pull < 1e-4,
      f"Ricci - 3q g = 0: {ric_ok}; Einstein tensor + Lambda g = 0: {ein_ok} with Lambda = {Lam}; scalar curvature {Rs}; "
      f"six small steps: law's distance^2 vs the metric {e_met:.0e}, placement vs the metric {e_pull:.0e} (relative, step 1e-6)")
results["einstein"] = dict(ricci_is_3qg=bool(ric_ok), einstein=bool(ein_ok), Lambda=str(Lam), scalar=str(Rs), law_metric=e_met, placement=e_pull)

# D5: the dial's sign picks de Sitter or anti-de Sitter; any other q is a rescaling
qa = place(Ys, -1.0); e_ads = np.abs(-qa[:, 0]**2 - qa[:, 1]**2 + (qa[:, 2:]**2).sum(1) + 1).max()
qs_list = [-2.0, -0.5, 0.3, 1.7]
e_scale = max(np.abs(oplus(Ys * 0.4, Ys[::-1] * 0.4, qq) - oplus(math.sqrt(abs(qq)) * Ys * 0.4, math.sqrt(abs(qq)) * Ys[::-1] * 0.4,
                                                              math.copysign(1.0, qq)) / math.sqrt(abs(qq))).max() for qq in qs_list)
check("D5 (exact) the dial's sign picks the spacetime: q > 0 gives de Sitter (Lambda = 3q > 0, an expanding universe), q < 0 gives anti-de Sitter (Lambda = 3q < 0, -W^2 - T^2 + X^2 + Y^2 + Z^2 = -1), and any other q is one of these rescaled. The symbolic Einstein check in D4 holds for either sign",
      e_ads < 1e-12 and e_scale < 1e-12,
      f"400 points at q = -1 on anti-de Sitter to {e_ads:.0e}; rescaling for q = {qs_list} {e_scale:.0e}")
results["sign"] = dict(ads=float(e_ads), rescale=float(e_scale))

# D6: the expansion.  De Sitter's global time slices: the event reached from rest by a time kick of rapidity tau has T = sinh(tau);
# the space slice T = sinh(tau) is a 3-sphere of radius cosh(tau) (round 14's a(t) = cosh(kt)/k). Read it from the law's own moves.
taus = np.array([0.0, 0.5, 1.0, 1.5])
rows6 = []
for tau in taus:
    a0 = math.tanh(tau / 2)                                       # chart position of rapidity tau along time (2 artanh a0 = tau)
    P0 = place(np.array([[a0, 0, 0, 0]]))[0]
    # events on the slice T = sinh(tau): chart points with 2t/(1 + eta) = sinh(tau); sample, then measure the sphere's radius in 5D
    ang = rg.normal(size=(500, 3)); ang /= np.linalg.norm(ang, axis=1)[:, None]
    th = rg.uniform(0, 0.9 * math.pi, (500, 1))      # the chart leaves out the far side, W = -1
    Pv = np.concatenate([np.cosh(tau) * np.cos(th), np.full((500, 1), math.sinh(tau)), np.cosh(tau) * np.sin(th) * ang], 1)
    on = np.abs(amb(Pv) - 1).max()
    # pull those 5D points back to the law's chart and check the law's spatial kicks keep them on de Sitter
    W, Tt, V = Pv[:, :1], Pv[:, 1:2], Pv[:, 2:]
    xc = np.concatenate([Tt, V], 1) / (1 + W)
    back = np.abs(place(xc) - Pv).max()
    rows6.append(dict(tau=float(tau), T=float(P0[1]), sinh=float(math.sinh(tau)), radius=float(P0[0]), sampled_radius=float(np.linalg.norm(Pv[:, [0, 2, 3, 4]], axis=1).mean()),
                      cosh=float(math.cosh(tau)), on_ds=float(on), chart_round_trip=float(back)))
e6 = max(max(abs(r["T"] - r["sinh"]), abs(r["radius"] - r["cosh"]), abs(r["sampled_radius"] - r["cosh"]), r["on_ds"], r["chart_round_trip"]) for r in rows6)
check("D6 (exact) the law's time kick of rapidity tau from rest reaches the event (W, T) = (cosh tau, sinh tau). The space slice of constant T through it is (textbook) the 3-sphere W^2 + X^2 + Y^2 + Z^2 = cosh^2 tau, so its radius is the event's W, cosh tau, and the law's chart covers that slice except its far side (W = -1). These are the slices of rounds 14-16, radius cosh(kt)/k with k = sqrt(q): the expansion they assumed is the law's own spacetime",
      e6 < 1e-10,
      "; ".join(f"tau {r['tau']}: T {r['T']:.4f} = sinh {r['sinh']:.4f}, slice radius {r['radius']:.4f} = cosh {r['cosh']:.4f}" for r in rows6) + f"; worst {e6:.0e}")
results["slices"] = rows6

# for the picture: 2D de Sitter (t, x) from the law: a grid of the law's kicks in space and time, and two light rays
grid = []
for i in range(-4, 5):
    for j in range(-7, 8):
        p = np.zeros((1, 4))
        for _ in range(abs(j)):
            p = oplus(np.array([[0, math.copysign(math.tan(math.pi / 16), j), 0, 0]]), p, 1.0)
        for _ in range(abs(i)):
            p = oplus(np.array([[math.copysign(math.tanh(0.25 / 2), i), 0, 0, 0]]), p, 1.0)
        P_ = place(p)[0]; grid.append(dict(i=i, j=j, W=float(P_[0]), T=float(P_[1]), X=float(P_[2])))
results["grid"] = grid

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
