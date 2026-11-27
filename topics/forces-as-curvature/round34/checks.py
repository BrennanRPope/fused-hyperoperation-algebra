"""DRAFT round 34: gravitational waves on spacetime glued from flat pieces (Regge's method, round 11), in 4 dimensions.
Round 11 showed that on glued flat pieces curvature is exactly the law's loop twist (the deficit angle at each joint), and that
Regge's action, joint size times deficit, tends to Einstein's (1/2) integral of R. Round 30 found that the law leaves 2 functions
per point free on a slice: the gravitational-wave data. This round asks what the glued pieces do with those 2 functions: do ripples
of the pieces' edge lengths move, at what speed, with how many polarizations, and what does the gluing scale change?
Lattice: the Freudenthal (Kuhn) triangulation of the hypercubic lattice Z^4, 24 flat 4-simplices per unit cell, 15 edges per
vertex (directions in {0,1}^4 minus 0), 50 triangle (joint) types. Variables: squared edge lengths. Second-order Regge action,
by Schlafli's identity: H_ef = sum over joints t of dA_t/dx_e * d(deficit_t)/dx_f. In Fourier form a 15 x 15 matrix M(k).
Checks: (G1) the flat lattice and the symmetries of M(k); (G2) long wavelengths: M(k) becomes Einstein's linearized action;
(G3) real frequency: exactly two extra flat directions appear on a light cone, the two polarizations; (G4) the light cone itself:
its exact shape, its long-wavelength speed, and what the gluing scale would have to be to show up.
Numbers: numpy finite differences (step 1e-6) for the angle derivatives, sympy for the continuum comparison."""
import json, io, time, itertools
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ---------------------------------------------------------------- the lattice
D = 4
def simplex(p):
    v = [np.zeros(D, int)]
    for a in p:
        w = v[-1].copy(); w[a] += 1; v.append(w)
    return [tuple(x) for x in v]
SIMP = [simplex(p) for p in itertools.permutations(range(D))]            # 24 simplex types per unit cell
DIRS = [u for u in itertools.product((0, 1), repeat=D) if any(u)]       # 15 edge directions
DIDX = {u: i for i, u in enumerate(DIRS)}
def edge_key(a, b):
    a, b = np.array(a), np.array(b); d = b - a
    if np.all(d >= 0): return tuple(a), DIDX[tuple(d)]
    return tuple(b), DIDX[tuple(-d)]
def canon_tri(t):
    t = sorted(t); base = np.array(t[0]); return tuple(tuple(np.array(x) - base) for x in t)
TRIS = sorted({canon_tri(c) for s in SIMP for c in itertools.combinations(s, 3)})   # 50 joint types
INST = [tuple(tuple(np.array(x) + np.array(n)) for x in s) for s in SIMP for n in itertools.product(range(-2, 2), repeat=D)]
def simplices_with(tri):
    st = set(tri); return [s for s in INST if st <= set(s)]
def coords_from_sq(sq):
    G = np.zeros((4, 4))
    for i in range(1, 5):
        for j in range(1, 5):
            G[i-1, j-1] = (sq[(0, i)] + sq[(0, j)] - (sq[(min(i, j), max(i, j))] if i != j else 0)) / 2
    return np.vstack([np.zeros(4), np.linalg.cholesky(G)])
def dihedral(X, hinge):
    others = [i for i in range(5) if i not in hinge]
    W = np.linalg.inv((X[1:] - X[0]).T)
    nrm = lambda i: W[i-1] if i > 0 else -W.sum(axis=0)
    n1, n2 = nrm(others[0]), nrm(others[1])
    return np.pi - np.arccos(np.clip(n1 @ n2 / np.linalg.norm(n1) / np.linalg.norm(n2), -1, 1))
def area_grad(x1, x2, x3):
    A = np.sqrt(2 * (x1*x2 + x2*x3 + x3*x1) - (x1**2 + x2**2 + x3**2)) / 4
    return A, np.array([2*(x2 + x3 - x1), 2*(x1 + x3 - x2), 2*(x1 + x2 - x3)]) / (32 * A)
PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
deficits = []; ST = []; hstep = 1e-6
for t in TRIS:
    es = [edge_key(t[i], t[j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    A, g = area_grad(*[float(np.dot(DIRS[e[1]], DIRS[e[1]])) for e in es])
    dd = {}; angle_sum = 0.0
    for s in simplices_with(t):
        lab = {v: i for i, v in enumerate(s)}; hinge = [lab[v] for v in t]
        sq0 = {(i, j): float(np.sum((np.array(s[j]) - np.array(s[i]))**2)) for i, j in PAIRS}
        angle_sum += dihedral(coords_from_sq(sq0), hinge)
        for (i, j) in PAIRS:
            sp_, sm_ = dict(sq0), dict(sq0); sp_[(i, j)] += hstep; sm_[(i, j)] -= hstep
            dth = (dihedral(coords_from_sq(sp_), hinge) - dihedral(coords_from_sq(sm_), hinge)) / (2 * hstep)
            ek = edge_key(s[i], s[j]); dd[ek] = dd.get(ek, 0.0) - dth
    deficits.append(2 * np.pi - angle_sum)
    ST.append((A, list(zip(es, g)), list(dd.items())))
def Mk(k):
    k = np.asarray(k, complex); M = np.zeros((15, 15), complex)
    for A, dA, dd in ST:
        al = np.zeros(15, complex); be = np.zeros(15, complex)
        for (p, u), v in dA: al[u] += v * np.exp(-1j * k @ np.array(p))
        for (p, u), v in dd: be[u] += v * np.exp(1j * k @ np.array(p))
        M += np.outer(al, be)
    return M
def gauge_vec(k, mu):          # moving every vertex by e_mu exp(i k.p): dx_u = 2 u_mu (exp(i k.u) - 1)
    return np.array([2 * u[mu] * (np.exp(1j * np.dot(k, u)) - 1) for u in DIRS])
sv = lambda k: np.linalg.svd(Mk(k), compute_uv=False)
NULL = 1e-8

# ---------------------------------------------------------------- G1
say("== G1. The flat lattice and the symmetries of the second-order Regge action")
worst_def = float(np.max(np.abs(deficits)))
rng = np.random.default_rng(0)
herm, gauge_res, nnull = [], [], []
for _ in range(5):
    k = rng.uniform(-np.pi, np.pi, 4); M = Mk(k)
    herm.append(np.linalg.norm(M - M.conj().T) / np.linalg.norm(M))
    gauge_res.append(max(np.linalg.norm(M @ gauge_vec(k, mu)) / np.linalg.norm(gauge_vec(k, mu)) for mu in range(4)) / np.linalg.norm(M, 2))
    nnull.append(int(np.sum(sv(k) < NULL)))
M0 = Mk(np.zeros(4)).real; M0 = (M0 + M0.T) / 2
w0, V0 = np.linalg.eigh(M0)
stiff = sorted(w0[np.abs(w0) > NULL].round(9).tolist()); nsoft = int(np.sum(np.abs(w0) < NULL))
idx = [(m, n) for m in range(4) for n in range(m, 4)]
P = np.array([[u[m] * u[n] * (1 if m == n else 2) for (m, n) in idx] for u in DIRS], float)   # constant metric change h -> edge changes
Nsp = V0[:, np.abs(w0) < NULL]; Msp = V0[:, np.abs(w0) > NULL]
Qp, _ = np.linalg.qr(P); ext = Nsp - Qp @ (Qp.T @ Nsp); Xm = np.linalg.svd(ext)[0][:, :1]   # the soft mode that is not a metric change
perm_ok = True
for pm in itertools.permutations(range(4)):
    Pm = np.zeros((15, 15))
    for a, u in enumerate(DIRS): Pm[DIDX[tuple(u[pm[i]] for i in range(4))], a] = 1
    perm_ok &= min(np.linalg.norm(Pm @ Xm - Xm), np.linalg.norm(Pm @ Xm + Xm)) < 1e-8
check("G1 (numbers) spacetime glued from the Freudenthal triangulation of the 4D hypercubic lattice is flat: every joint's deficit, 2 pi minus the angles meeting there, is zero. The second-order Regge action M(k) is Hermitian (Schlafli's identity). Moving the vertices costs nothing (4 exact zero directions at every wavevector, as diffeomorphisms should), and there is exactly one more flat direction at every wavevector. At k = 0 that extra one is the lattice's own flex: a change of edge lengths, the same for every choice of axes, that keeps every piece flat to first order but is not a change of metric. The other 4 directions at k = 0 are stiff",
      worst_def < 1e-12 and max(herm) < 1e-9 and max(gauge_res) < 1e-8 and set(nnull) == {5} and nsoft == 11 and len(stiff) == 4 and perm_ok,
      f"joint types {len(TRIS)}, worst deficit {worst_def:.1e}; |M - M^H|/|M| <= {max(herm):.1e}; vertex moves |M g|/(|M||g|) <= {max(gauge_res):.1e}; zero directions at 5 random k: {nnull}; at k = 0: {nsoft} soft (10 metric + 1 flex), stiff eigenvalues {stiff}; flex unchanged by all 24 axis permutations: {perm_ok}")
results["G1"] = dict(joints=len(TRIS), worst_deficit=worst_def, herm=max(herm), gauge=max(gauge_res), nnull=nnull, nsoft=nsoft, stiff=stiff, flex_symmetric=perm_ok)

# ---------------------------------------------------------------- G2
say("== G2. Long wavelengths: the glued pieces' action is Einstein's linearized action")
B = np.hstack([P, Xm, Msp]); Binv = np.linalg.inv(B)
def Eh(k):
    """metric block: the 4 stiff directions integrated out, the flex held at 0 (it is an exact flat direction, so fixing it is a
    choice of how to describe the same flat shape, like fixing coordinates)"""
    Mb = B.T @ Mk(k) @ B; a, b, c, d = Mb[:11, :11], Mb[:11, 11:], Mb[11:, :11], Mb[11:, 11:]
    return (a - b @ np.linalg.solve(d, c))[:10, :10]
def FP(k):   # Fierz-Pauli symbol: k^2 h.h - 2 (k h)^2 + 2 (k h k) tr h - k^2 (tr h)^2, on idx coordinates (off-diagonals h_mn)
    k = np.asarray(k, complex); k2 = k @ k; F = np.zeros((10, 10), complex)
    def H(b):
        h = np.zeros((4, 4)); m, n = idx[b]; h[m, n] = h[n, m] = 1; return h
    for i in range(10):
        for j in range(10):
            hi, hj = H(i), H(j)
            F[i, j] = k2 * np.sum(hi * hj) - 2 * (k @ hi) @ (k @ hj) + (k @ hi @ k) * np.trace(hj) + (k @ hj @ k) * np.trace(hi) - k2 * np.trace(hi) * np.trace(hj)
    return F
dirs4 = [np.array(v, float) / np.linalg.norm(v) for v in ([1, 0, 0, 0], [1, 1, 0, 0], [1, 2, -1, 1], [0.3, -1, 2, 0.7])]
G2rows = []
for kd in dirs4:
    for e in (0.04, 0.02, 0.01, 0.005):
        k = e * kd; E = Eh(k); F = FP(k)
        G2rows.append(dict(dir=kd.round(3).tolist(), k=e, rel=float(np.linalg.norm(E + F / 4) / np.linalg.norm(F / 4))))
slopes = []
for i in range(0, len(G2rows), 4):
    r = G2rows[i:i + 4]; slopes.append(float(np.polyfit(np.log([x["k"] for x in r]), np.log([x["rel"] for x in r]), 1)[0]))
# continuum: average of the second-order part of sqrt(g) R for g = 1 + eps H cos(k.x), exact in sympy
X = sp.symbols('x0:4'); ep = sp.symbols('epsilon')
kv = [1, 2, -1, 1]; Hm = sp.Matrix([[1, 2, 0, -1], [2, -1, 1, 0], [0, 1, 3, 2], [-1, 0, 2, 1]])
ph = sum(a * b for a, b in zip(kv, X)); g = sp.eye(4) + ep * Hm * sp.cos(ph)
gi = sp.eye(4) - ep * Hm * sp.cos(ph) + ep**2 * (Hm * Hm) * sp.cos(ph)**2
tr2 = lambda e: sp.series(sp.expand(e), ep, 0, 3).removeO()
Gam = [[[tr2(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4)) / 2) for c in range(4)] for b in range(4)] for a in range(4)]
Ric = lambda b, c: sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c]) + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(4)) for a in range(4))
Rs = tr2(sum(gi[b, c] * Ric(b, c) for b in range(4) for c in range(4)))
sg = tr2(1 + ep * sp.cos(ph) * Hm.trace() / 2 + ep**2 * sp.cos(ph)**2 * (Hm.trace()**2 / 8 - (Hm * Hm).trace() / 4))
L2 = sp.expand(tr2(sg * Rs)).coeff(ep, 2); s_ = sp.Symbol('s')
avg = sp.nsimplify(sp.integrate(L2.subs({X[1]: 0, X[2]: 0, X[3]: 0}).subs(X[0], s_), (s_, 0, 2 * sp.pi)) / (2 * sp.pi))
Fnum = float(np.real(FP(np.array(kv, float)) @ np.array([float(Hm[m, n]) for (m, n) in idx]) @ np.array([float(Hm[m, n]) for (m, n) in idx])))
ratio = sp.nsimplify(avg / sp.nsimplify(Fnum))
check("G2 (numbers plus exact sympy; textbook Fierz-Pauli and Roček-Williams 1981) at long wavelengths the glued pieces' second-order action, with the stiff directions integrated out and the flex held at 0, is -F(k)/4 on the metric, where F is the Fierz-Pauli form of linearized gravity. The mismatch falls in proportion to k in every direction tested. In the continuum, the average second-order part of sqrt(g) R for a ripple is exactly -F/8. With one lattice vertex per unit 4-volume and <cos^2> = 1/2, both give -F/16 per unit volume: Regge's sum of joint size times deficit is exactly Einstein's (1/2) integral of sqrt(g) R for ripples",
      max(r["rel"] for r in G2rows if r["k"] == 0.005) < 0.02 and all(0.9 < s < 1.1 for s in slopes) and ratio == sp.Rational(-1, 8),
      f"relative mismatch |E + F/4|/|F/4| at k = 0.04, 0.02, 0.01, 0.005: " + "; ".join(f"dir {G2rows[i]['dir']}: " + ", ".join(f"{x['rel']:.1e}" for x in G2rows[i:i+4]) for i in range(0, len(G2rows), 4)) + f"; slopes {[round(s, 3) for s in slopes]}; continuum <(sqrt g R)_2> = {avg} with F = {Fnum:g}, ratio {ratio}")
results["G2"] = dict(rows=G2rows, slopes=slopes, continuum_avg=str(avg), F=Fnum, ratio=str(ratio))

# ---------------------------------------------------------------- G3
say("== G3. Real frequency: two extra flat directions, the two polarizations, appear exactly on a light cone")
def cone(ks): return 2 * np.arcsinh(np.sqrt(np.sum(np.sin(np.asarray(ks) / 2)**2)))   # plain lattice wave equation's frequency
def kvec(ks, w): return np.concatenate([ks, [1j * w]])                                    # k_4 -> i omega (Euclidean continued)
# (a) generic complex points on sum_mu sin^2(k_mu/2) = 0, and just off it
gen = []
for _ in range(6):
    k = rng.normal(size=4) * 0.8 + 1j * rng.normal(size=4) * 0.3
    k[3] = 2 * np.arcsin(np.sqrt(-np.sum(np.sin(k[:3] / 2)**2) + 0j))
    k_off = k.copy(); k_off[3] *= 1.01
    gen.append((int(np.sum(sv(k) < NULL)), int(np.sum(sv(k_off) < NULL))))
# (b) scan in frequency at one spatial wavevector, for the picture
eps = 0.1; ks_axis = np.array([eps, 0, 0]); wl = cone(ks_axis)
scan_w = np.linspace(0.985, 1.01, 121) * eps
scan = [sv(kvec(ks_axis, w))[-7:].tolist() for w in scan_w]
# (c) polarizations: which metric shapes the two new flat directions carry, against the transverse-traceless + and x
def vec(h3):
    h = np.zeros((4, 4)); h[:3, :3] = h3; return np.array([h[m, n] for (m, n) in idx])
pol = []
for e in (0.1, 0.05, 0.025):
    ks = np.array([e, 0, 0]); k = kvec(ks, cone(ks))
    Vh = np.linalg.svd(Mk(k))[2]; N = Vh.conj().T[:, -7:]; Hn = (Binv @ N)[:10]; Qn, _ = np.linalg.qr(Hn)
    a_, c_ = np.array([0, 1., 0]), np.array([0, 0, 1.])
    ep_, ex_ = vec(np.outer(a_, a_) - np.outer(c_, c_)), vec(np.outer(a_, c_) + np.outer(c_, a_)) / 2
    htt = np.zeros(10); htt[idx.index((3, 3))] = 1
    res = lambda v: float(np.linalg.norm(v - Qn @ (Qn.conj().T @ v)) / np.linalg.norm(v))
    pol.append(dict(k=e, plus=res(ep_), cross=res(ex_), h44=res(htt)))
xslope = float(np.polyfit(np.log([p["k"] for p in pol]), np.log([p["cross"] for p in pol]), 1)[0])
check("G3 (numbers) continuing the fourth direction to real frequency (k_4 -> i omega), the glued pieces gain exactly 2 more flat directions, and only on the surface sum_mu sin^2(k_mu/2) = 0: at 6 random complex points on it there are 7 (4 vertex moves, the flex, and 2 new), and 1% off it there are 5. Along a real spatial wavevector both appear at one frequency, so the two polarizations move together. Their metric shapes are the transverse-traceless + and x of a gravitational wave: + exactly, x up to a lattice correction that falls in proportion to k (this lattice has no mirror y -> -y, so the x wave picks up a small admixture at finite k). A time-time ripple, which would be a third polarization, is not among them",
      all(a == 7 and b == 5 for a, b in gen) and max(p["plus"] for p in pol) < 1e-5 and 0.9 < xslope < 1.1 and min(p["h44"] for p in pol) > 0.2,
      f"flat directions on / off the cone at 6 random complex points: {gen}; residuals of + , x , h_44 in the new directions' metric span at k = 0.1, 0.05, 0.025: " + "; ".join(f"{p['plus']:.1e}, {p['cross']:.1e}, {p['h44']:.2f}" for p in pol) + f"; x residual slope in k: {xslope:.3f}")
results["G3"] = dict(generic=gen, scan_w=(scan_w / eps).tolist(), scan=scan, cone_w=wl / eps, pol=pol, cross_slope=xslope)

# ---------------------------------------------------------------- G4
say("== G4. The light cone: its exact shape, its long-wavelength speed, and the gluing scale it would take to see it")
dirs3 = {"axis (1,0,0)": [1, 0, 0], "face (1,1,0)": [1, 1, 0], "body (1,1,1)": [1, 1, 1], "generic (2,-1,1/2)": [2, -1, 0.5]}
G4rows = []
for name, v in dirs3.items():
    kd = np.array(v, float) / np.linalg.norm(v)
    for e in (0.1, 0.3, 0.6, 1.0):
        ks = e * kd; wl = cone(ks)
        f6 = lambda r: sv(kvec(ks, r * e))[-6]; f7 = lambda r: sv(kvec(ks, r * e))[-7]
        r6 = minimize_scalar(f6, bounds=(wl / e - 0.03, wl / e + 0.03), method='bounded', options={'xatol': 1e-13}).x
        r7 = minimize_scalar(f7, bounds=(wl / e - 0.03, wl / e + 0.03), method='bounded', options={'xatol': 1e-13}).x
        G4rows.append(dict(dir=name, k=e, speed=r6, speed2=r7, cone=wl / e, beta_pred=-(1 + np.sum(kd**4)) / 24,
                           beta_meas=(r6 - 1) / e**2 if e == 0.1 else None))
worst_cone = max(max(abs(r["speed"] - r["cone"]), abs(r["speed2"] - r["cone"])) for r in G4rows)
beta_err = max(abs(r["beta_meas"] - r["beta_pred"]) / abs(r["beta_pred"]) for r in G4rows if r["k"] == 0.1)
# what it would take to see it: GW170817 / GRB 170817A, -3e-15 <= (v_gw - c)/c <= 7e-16 (Abbott et al. 2017)
c, f_gw, bound = 2.998e8, 100.0, 3e-15
kphys = 2 * np.pi * f_gw / c; worst_beta = 1 / 12
a_max = np.sqrt(bound / worst_beta) / kphys
lP = 1.616e-35; shift_planck = worst_beta * (kphys * lP)**2
check("G4 (numbers; textbook lattice dispersion, GW170817 bound from Abbott et al. 2017) the cone is exactly the plain lattice wave equation's, 4 sinh^2(omega/2) = sum_i 4 sin^2(k_i/2), at every wavevector tested up to one gluing length, not just at long wavelengths. So the waves move at light speed at long wavelengths and lag by beta (ka)^2 with beta = -(1 + sum k_i^4/k^4)/24: 1/12 along an axis, 1/18 along a body diagonal. Both polarizations lag the same, so there is no splitting. GW170817 limits the lag to 3 x 10^-15 at about 100 Hz, which only requires the gluing length a to be under about 9 cm; at the Planck length the lag is about 10^-83",
      worst_cone < 1e-6 and beta_err < 0.01 and 0.05 < a_max < 0.2,
      "speed omega/k found (two flat directions) vs cone: " + "; ".join(f"{r['dir']} ka={r['k']}: {r['speed']:.8f}, {r['speed2']:.8f} vs {r['cone']:.8f}" for r in G4rows) + f"; worst gap {worst_cone:.1e}; beta at ka = 0.1 vs -(1 + sum khat^4)/24: " + ", ".join(f"{r['dir']} {r['beta_meas']:.5f} vs {r['beta_pred']:.5f}" for r in G4rows if r['k'] == 0.1) + f" (worst {beta_err:.1e}); a_max at 100 Hz = {a_max:.3f} m; lag at Planck length {shift_planck:.1e}")
results["G4"] = dict(rows=G4rows, worst_cone_gap=worst_cone, beta_rel_err=beta_err, a_max_m=a_max, f_gw=f_gw, bound=bound, planck_shift=shift_planck)

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
