"""DRAFT round 35: the energy gravitational waves carry, and how that energy bends spacetime.
Round 34 built gravitational waves on spacetime glued from flat pieces, with the law's loop twist (the deficit angle at each joint,
round 11) as curvature: 2 polarizations at light speed. This round follows their energy. (E1) the energy a wave carries on the
pieces: the pieces' action for the two wave shapes, and the speed its energy travels at. (E2) how that energy bends spacetime:
averaged over a wavelength, a wave's second-order curvature is exactly that of a beam of light-like matter with the wave's energy,
and on a slice it raises the law's local dial by the same rule as matter (rounds 23, 25, 30). (E3) the energy an orbiting pair
sends off, by integrating the wave's energy over a sphere, for any eccentricity. (E4) the Hulse-Taylor binary pulsar's shrinking
orbit. The lattice block below is copied verbatim from round 34's checks.py. E2 is exact sympy; E3 numbers (spectral derivatives,
quadrature); E4 uses published pulsar parameters."""
import json, io, time, itertools
import numpy as np
import sympy as sp

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

# metric coordinates and the reduced (stiff integrated out, flex held at 0) metric block, as in round 34's G2
idx = [(m, n) for m in range(4) for n in range(m, 4)]
P = np.array([[u[m] * u[n] * (1 if m == n else 2) for (m, n) in idx] for u in DIRS], float)
M0 = Mk(np.zeros(4)).real; M0 = (M0 + M0.T) / 2; w0, V0 = np.linalg.eigh(M0)
Nsp = V0[:, np.abs(w0) < NULL]; Msp = V0[:, np.abs(w0) > NULL]
Qp, _ = np.linalg.qr(P); Xm = np.linalg.svd(Nsp - Qp @ (Qp.T @ Nsp))[0][:, :1]
B = np.hstack([P, Xm, Msp])
def Eh(k):
    Mb = B.T @ Mk(k) @ B; a, b, c, d = Mb[:11, :11], Mb[:11, 11:], Mb[11:, :11], Mb[11:, 11:]
    return (a - b @ np.linalg.solve(d, c))[:10, :10]
def hvec(h):  return np.array([h[m, n] for (m, n) in idx])     # idx coordinates: off-diagonal coordinate = h_mn
def cone(ks): return 2 * np.arcsinh(np.sqrt(np.sum(np.sin(np.asarray(ks) / 2)**2)))

# ---------------------------------------------------------------- E1
say("== E1. The energy a wave carries on the pieces")
E1rows = []
for kd in ([1, 0, 0, 0], [1, 1, 0, 0], [1, 2, -1, 1]):
    kd = np.array(kd, float) / np.linalg.norm(kd)
    basis = np.linalg.svd(kd[None, :])[2][1:]                       # 3 directions across the wave
    a_, c_ = basis[0], basis[1]
    ep, ex = np.outer(a_, a_) - np.outer(c_, c_), np.outer(a_, c_) + np.outer(c_, a_)
    for e in (0.02, 0.01, 0.005):
        E = Eh(e * kd); V = np.array([hvec(ep), hvec(ex)]).T
        blk = V.T @ E @ V                                           # 2x2 action on (+, x)
        target = -e**2 / 4 * np.sum(ep * ep)                         # -(k^2/4) h.h, h.h = 2 for both
        E1rows.append(dict(dir=kd.round(3).tolist(), k=e, plus=float(blk[0, 0].real / target), cross=float(blk[1, 1].real / target),
                           mix=float(abs(blk[0, 1]) / abs(target))))
# speed of the energy: group velocity from the cone, along an axis and a body diagonal
vg = []
for name, v in (("axis", [1, 0, 0]), ("body", [1, 1, 1])):
    d = np.array(v, float) / np.linalg.norm(v)
    for ka in (0.1, 0.3, 0.6, 1.0):
        hh = 1e-6; g = (cone((ka + hh) * d) - cone((ka - hh) * d)) / (2 * hh)
        beta = -(1 + np.sum(d**4)) / 24
        vg.append(dict(dir=name, ka=ka, group=float(g), phase=float(cone(ka * d) / ka), lowk=1 + 3 * beta * ka**2))
worst = max(max(abs(r["plus"] - 1), abs(r["cross"] - 1)) for r in E1rows if r["k"] == 0.005)
vg_ok = all(abs(r["group"] - r["lowk"]) < 3e-5 for r in vg if r["ka"] == 0.1) and all(r["group"] < r["phase"] < 1 for r in vg)
check("E1 (numbers; textbook wave energy and group velocity) on the pieces, the action for the two wave shapes + and x is -(k^2/4) h.h, the same for both and with no mixing, up to a lattice correction that falls in proportion to k. With Regge's action normalised as (1/2) integral of sqrt(g) R (that is, G = 1/8 pi), this is linearized Einstein gravity's energy density (1/32 pi G) <hdot_ij hdot_ij>: both shapes carry energy alike. The energy moves at the cone's group velocity, which is light speed at long wavelengths and lags 3 times as much as the wave crests do, (ka)^2/4 along an axis",
      worst < 0.01 and max(r["mix"] for r in E1rows) < 0.01 and vg_ok,
      "(+, x) stiffness / -(k^2/4) h.h at k = 0.02, 0.01, 0.005: " + "; ".join(f"dir {E1rows[i]['dir']}: " + ", ".join(f"{r['plus']:.4f}/{r['cross']:.4f}" for r in E1rows[i:i+3]) for i in range(0, len(E1rows), 3)) + f"; worst at 0.005 {worst:.1e}, mixing <= {max(r['mix'] for r in E1rows):.1e}; group velocity: " + ", ".join(f"{r['dir']} ka={r['ka']}: {r['group']:.6f} (crests {r['phase']:.6f}, 1+3 beta (ka)^2 = {r['lowk']:.6f})" for r in vg))
results["E1"] = dict(rows=E1rows, group=vg)

# ---------------------------------------------------------------- E2
say("== E2. How a wave's energy bends spacetime (exact, second order, averaged over a wavelength; G = 1)")
t, x, y, z = X = sp.symbols('t x y z', real=True); A_, B_, ph = sp.symbols('A B phi', real=True); w = sp.symbols('omega', positive=True); e_ = sp.symbols('epsilon')
u = t - z
h = sp.zeros(4); h[1, 1] = A_ * sp.cos(w * u); h[2, 2] = -h[1, 1]; h[1, 2] = h[2, 1] = B_ * sp.cos(w * u + ph)   # + and x, moving along z
eta = sp.diag(-1, 1, 1, 1); g = eta + e_ * h; gi = eta - e_ * eta * h * eta + e_**2 * eta * h * eta * h * eta
tr2 = lambda E: sp.series(sp.expand(E), e_, 0, 3).removeO()
def christoffel(g, gi, Xs, n):
    return [[[tr2(sum(gi[a, d] * (sp.diff(g[d, b], Xs[c]) + sp.diff(g[d, c], Xs[b]) - sp.diff(g[b, c], Xs[d])) for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]
def ricci(Gm, Xs, n, b, c):
    return tr2(sum(sp.diff(Gm[a][b][c], Xs[a]) - sp.diff(Gm[a][b][a], Xs[c]) + sum(Gm[a][a][d] * Gm[d][b][c] - Gm[a][c][d] * Gm[d][b][a] for d in range(n)) for a in range(n)))
s_ = sp.Symbol('s')
def avg(E): return sp.simplify(sp.integrate(sp.expand(E).subs({t: s_, x: 0, y: 0, z: 0}), (s_, 0, 2 * sp.pi / w)) * w / (2 * sp.pi))
Gm = christoffel(g, gi, X, 4)
R1 = sp.zeros(4); R2 = sp.zeros(4)
for b in range(4):
    for c in range(4):
        r = ricci(Gm, X, 4, b, c); R1[b, c] = sp.simplify(r.coeff(e_, 1)); R2[b, c] = avg(r.coeff(e_, 2))
dh = lambda m: [sp.diff(h[i, j], X[m]) for i in range(1, 4) for j in range(1, 4)]
tI = sp.Matrix(4, 4, lambda m, n: avg(sum(p * q for p, q in zip(dh(m), dh(n)))) / (32 * sp.pi))      # Isaacson's wave energy tensor
einstein_ok = sp.simplify(R2 + 8 * sp.pi * tI) == sp.zeros(4)
R2scalar = sp.simplify(sum(eta[a, a] * R2[a, a] for a in range(4)))
null_ok = sp.simplify(tI[0, 0] - tI[3, 3]) == 0 and sp.simplify(tI[0, 0] + tI[0, 3]) == 0
# on the slice t = 0: 3-metric delta + h_ij, extrinsic curvature K_ij = -(1/2) dh_ij/dt
h3 = h[1:, 1:]; g3 = sp.eye(3) + e_ * h3; g3i = sp.eye(3) - e_ * h3 + e_**2 * h3 * h3; Y = [x, y, z]
G3 = christoffel(g3, g3i, Y, 3)
R3 = tr2(sum(g3i[b, c] * ricci(G3, Y, 3, b, c) for b in range(3) for c in range(3)))
K = -e_ * sp.diff(h3, t) / 2
KK = tr2(sum(K[i, j] * K[k, l] * g3i[i, k] * g3i[j, l] for i in range(3) for j in range(3) for k in range(3) for l in range(3)))
trK = tr2(sum(g3i[i, j] * K[i, j] for i in range(3) for j in range(3)))
R3_2, KK2, K2 = avg(R3.coeff(e_, 2)), avg(KK.coeff(e_, 2)), avg((trK**2).coeff(e_, 2))
rho_gw = sp.simplify((KK2 - K2 - R3_2) / (16 * sp.pi))
dial_shift = sp.simplify(16 * sp.pi * rho_gw / 6)                    # law's local dial = mean bending of the slice = 3R/6
rule_ok = sp.simplify(rho_gw - tI[0, 0]) == 0 and sp.simplify(dial_shift - 8 * sp.pi * tI[0, 0] / 3) == 0
check("E2 (exact sympy; textbook Isaacson 1968) a + and x wave solves Einstein's equations at first order. At second order, averaged over a wavelength, its curvature is exactly -8 pi times Isaacson's energy tensor t = <d h d h>/32 pi, which is the curvature a beam of light-like matter with the wave's energy would make: energy density equal to its flux, no pressure across it, Ricci scalar 0. On a slice, the wave's own ripples and its motion each supply half of that energy, so the smooth slice underneath must bend by 16 pi rho_gw, and the law's local dial (the mean bending, round 30) rises by 8 pi rho_gw / 3: exactly the rule rounds 23 and 25 found for matter. So in the law's terms a wave's energy feeds the dial like any other matter",
      R1 == sp.zeros(4) and einstein_ok and R2scalar == 0 and null_ok and rule_ok and sp.simplify(KK2 + R3_2) == 0,
      f"first-order Ricci zero: {R1 == sp.zeros(4)}; <R2_tt> = {R2[0, 0]}, <R2_tz> = {R2[0, 3]}, <R2_zz> = {R2[3, 3]}, others 0; t_tt = {tI[0, 0]}; <R2> + 8 pi t = 0: {einstein_ok}; <R2 scalar> = {R2scalar}; light-like (t_tt = t_zz = -t_tz): {null_ok}; slice: <K_ij K^ij> = {KK2}, <K^2> = {K2}, <3R2> = {R3_2}; rho_gw = {rho_gw}; dial shift = {dial_shift} = 8 pi rho_gw/3: {rule_ok}")
results["E2"] = dict(R2_tt=str(R2[0, 0]), R2_tz=str(R2[0, 3]), t_tt=str(tI[0, 0]), KK=str(KK2), K2=str(K2), R3_2=str(R3_2), rho_gw=str(rho_gw), dial_shift=str(dial_shift), einstein=einstein_ok, null=null_ok, rule=rule_ok)

# ---------------------------------------------------------------- E3
say("== E3. The energy an orbiting pair sends off (G = c = 1)")
# (a) integrating the wave energy over a sphere: P = (1/8 pi) * integral dOmega Lambda_ijkl Q'''_ij Q'''_kl, against (1/5) Q'''.Q'''
nth, nph = 200, 400
th = (np.arange(nth) + 0.5) * np.pi / nth; phs = np.arange(nph) * 2 * np.pi / nph
TH, PH = np.meshgrid(th, phs, indexing="ij"); wts = np.sin(TH) * (np.pi / nth) * (2 * np.pi / nph)
nv = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
Pp = np.eye(3)[None, None] - nv[..., :, None] * nv[..., None, :]
def sphere_power(Q):            # Q: traceless symmetric 3x3 standing for Q'''
    PQP = np.einsum('abik,kl,ablj->abij', Pp, Q, Pp); tr = np.einsum('abii->ab', PQP)
    tt = PQP - 0.5 * Pp * tr[..., None, None]
    return float(np.sum(wts * np.einsum('abij,abij->ab', tt, tt)) / (8 * np.pi))
ang = []; rng = np.random.default_rng(0)
for _ in range(4):
    Qr = rng.normal(size=(3, 3)); Qr = Qr + Qr.T; Qr -= np.trace(Qr) / 3 * np.eye(3)
    ang.append(sphere_power(Qr) / (np.sum(Qr * Qr) / 5))
# (b) a Kepler orbit (total mass 1, reduced mass mu, semi-major axis 1), Q_ij = mu (x_i x_j - r^2 delta/3), third time derivative spectrally
def kepler_power(ecc, N=8192, mu=0.25):
    Mm = 2 * np.pi * np.arange(N) / N; Ea = Mm.copy()
    for _ in range(60): Ea = Ea - (Ea - ecc * np.sin(Ea) - Mm) / (1 - ecc * np.cos(Ea))
    xs = np.cos(Ea) - ecc; ys = np.sqrt(1 - ecc**2) * np.sin(Ea)          # period 2 pi with G M = a = 1
    pos = np.stack([xs, ys, 0 * xs]); r2 = xs**2 + ys**2
    Q = mu * (pos[:, None] * pos[None] - r2 / 3 * np.eye(3)[:, :, None])
    f = np.fft.fftfreq(N, d=1.0 / N)
    Q3 = np.real(np.fft.ifft((1j * f)**3 * np.fft.fft(Q, axis=-1), axis=-1))
    return float(np.mean(np.einsum('ijt,ijt->t', Q3, Q3)) / 5), 32 / 5 * mu**2 * (1 + 73 / 24 * ecc**2 + 37 / 96 * ecc**4) / (1 - ecc**2)**3.5
kep = []
for ecc in (0.0, 0.3, 0.617134, 0.8):
    num, pm = kepler_power(ecc, N=32768 if ecc > 0.7 else 8192); kep.append(dict(e=ecc, numeric=num, peters_mathews=pm, ratio=num / pm))
check("E3 (numbers; textbook quadrupole formula, Peters and Mathews 1963) adding up E1's energy flux over a sphere around a source gives the quadrupole formula, power = (G/5c^5) <Q'''_ij Q'''_ij>: the sphere integral is exactly 1/5 for any shape of source. For a pair on an orbit of any eccentricity this gives Peters and Mathews' power, (32/5) G^4 mu^2 M^3 / (c^5 a^5) times (1 + 73 e^2/24 + 37 e^4/96)/(1 - e^2)^(7/2): eccentric orbits radiate far more, about 12 times at the Hulse-Taylor pulsar's e = 0.617",
      all(abs(r - 1) < 1e-4 for r in ang) and all(abs(r["ratio"] - 1) < 1e-8 for r in kep),
      f"sphere integral / (1/5): {[round(r, 7) for r in ang]}; orbit power numeric / Peters-Mathews: " + ", ".join(f"e={r['e']}: {r['ratio']:.10f} (enhancement {r['peters_mathews'] / kep[0]['peters_mathews']:.3f})" for r in kep))
results["E3"] = dict(sphere=ang, kepler=kep)

# ---------------------------------------------------------------- E4
say("== E4. The Hulse-Taylor binary pulsar's shrinking orbit")
# Weisberg and Huang 2016 (ApJ 829, 55): Pb = 0.322997448918 d, e = 0.6171340, m_p = 1.438, m_c = 1.390 solar masses;
# GR prediction for dPb/dt = -2.40263e-12; measured, after the galactic correction, intrinsic / GR = 0.9983 +- 0.0016.
Tsun = 4.925490947e-6                     # G M_sun / c^3 in seconds
Pb, ecc, mp, mc = 0.322997448918 * 86400, 0.6171340, 1.438, 1.390
fe = (1 + 73 / 24 * ecc**2 + 37 / 96 * ecc**4) / (1 - ecc**2)**3.5
Pdot = -(192 * np.pi / 5) * Tsun**(5 / 3) * (Pb / (2 * np.pi))**(-5 / 3) * fe * mp * mc * (mp + mc)**(-1 / 3)
G_, c_, Msun = 6.67430e-11, 299792458.0, 1.32712440018e20 / 6.67430e-11
Mt, mu = (mp + mc) * Msun, mp * mc / (mp + mc) * Msun
a = (G_ * Mt * (Pb / (2 * np.pi))**2)**(1 / 3)
Lgw = 32 / 5 * G_**4 / c_**5 * mu**2 * Mt**3 / a**5 * fe
years = np.linspace(0, 42, 85); shift = 0.5 * (Pdot / Pb) * (years * 3.15576e7)**2
ok4 = abs(Pdot / -2.40263e-12 - 1) < 2e-3
check("E4 (numbers; published parameters, Weisberg and Huang 2016) with the pulsar's measured orbit and masses, E3's power predicts the orbital period shrinks by 2.40 x 10^-12 seconds per second, as the published general-relativity value. The pulsar, after correcting for its motion in the Galaxy, shows 0.9983 +- 0.0016 of that. The pair radiates about 8 x 10^24 watts in gravitational waves, and after 40 years its periastron arrives about 70 seconds earlier than it would on an orbit that did not radiate",
      ok4, f"f(e) = {fe:.4f}; dPb/dt = {Pdot:.5e} (published GR value -2.40263e-12, ratio {Pdot / -2.40263e-12:.5f}); semi-major axis {a / 1e3:.4e} km; gravitational-wave power {Lgw:.3e} W; cumulative periastron shift after 40 years {0.5 * (Pdot / Pb) * (40 * 3.15576e7)**2:.1f} s")
results["E4"] = dict(fe=fe, Pdot=Pdot, ratio_to_published=Pdot / -2.40263e-12, measured_ratio=[0.9983, 0.0016], a_m=a, L_W=Lgw, years=years.tolist(), shift_s=shift.tolist())

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
