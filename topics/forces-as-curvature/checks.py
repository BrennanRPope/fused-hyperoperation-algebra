"""DRAFT checks: every force as the curvature of a chain of circle operations.

Law: a o_k b = a + b + k*a*b, here with a, b square MATRICES (Jacobson's circle
operation on a matrix ring, scaled by k).  Lift: 1 + k*a.  Way back: abar = -a(1+ka)^-1.

Run: python3 checks.py  (numpy, scipy, sympy, mpmath).  Writes results.json.
"""
import json, math, itertools
import numpy as np
from scipy.linalg import expm
import sympy as sp
import mpmath as mp

rng = np.random.default_rng(20261127)
results, passed, failed = {}, 0, 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
    else:
        failed += 1
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))


def circ(a, b, k):
    return a + b + k * a @ b


def way_back(a, k):
    n = a.shape[0]
    return -a @ np.linalg.inv(np.eye(n) + k * a)


def chain(vals, k):
    """Compose link values left to right: v1 o v2 o ... (lift = L1 L2 ...)."""
    out = np.zeros_like(vals[0])
    for v in vals:
        out = circ(out, v, k)
    return out

# ---------------------------------------------------------------- Part A: the algebra
print("\n== A. The circle operation on matrices")
n = 4
for k in (1.0, -0.7, 0.3j):
    a = rng.normal(size=(n, n)) * 0.4 + 0j
    b = rng.normal(size=(n, n)) * 0.4 + 0j
    c = rng.normal(size=(n, n)) * 0.4 + 0j
    I = np.eye(n)
    lift_ok = np.abs(I + k * circ(a, b, k) - (I + k * a) @ (I + k * b)).max()
    assoc = np.abs(circ(circ(a, b, k), c, k) - circ(a, circ(b, c, k), k)).max()
    wb = np.abs(circ(a, way_back(a, k), k)).max()
    ab, bb = way_back(a, k), way_back(b, k)
    loop = chain([a, b, ab, bb], k)
    exact = k * (a @ b - b @ a) @ (I + k * ab) @ (I + k * bb)
    ident = np.abs(loop - exact).max()
    check(f"A1 lift turns o_k into the matrix product (k={k})", lift_ok < 1e-12, f"{lift_ok:.1e}")
    check(f"A2 associative on matrices (k={k})", assoc < 1e-12, f"{assoc:.1e}")
    check(f"A3 way back undoes a value (k={k})", wb < 1e-12, f"{wb:.1e}")
    check(f"A4 loop identity a o b o abar o bbar = k[a,b](1+k abar)(1+k bbar) (k={k})",
          ident < 1e-12, f"{ident:.1e}")
# commuting values: the loop closes
a = np.diag(rng.normal(size=n)); b = np.diag(rng.normal(size=n))
lp = chain([a, b, way_back(a, 1.0), way_back(b, 1.0)], 1.0)
check("A5 values that commute: the four-step loop leaves nothing", np.abs(lp).max() < 1e-13)
# small loop: leftover / eps^2 -> k[A,B]
A = rng.normal(size=(n, n)); B = rng.normal(size=(n, n))
errs = []
for eps in (1e-1, 1e-2, 1e-3):
    k = 0.8
    lp = chain([eps * A, eps * B, way_back(eps * A, k), way_back(eps * B, k)], k)
    errs.append(np.abs(lp / eps**2 - k * (A @ B - B @ A)).max())
check("A6 small loop: leftover/eps^2 -> k[A,B], error falls like eps", errs[2] < 1e-2 and errs[1] / errs[2] > 8,
      ", ".join(f"{e:.1e}" for e in errs))

# ---------------------------------------------------------------- Part B: gauge plaquettes
print("\n== B. Plaquettes for U(1), SU(2), SU(3): the k*ab term is the self-interaction")
def su_basis(N):
    if N == 1:
        return [np.array([[1.0 + 0j]])]
    T = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; T.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; T.append(m)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(l):
            m[i, i] = 1
        m[l, l] = -l
        T.append(m / math.sqrt(2 * l * (l + 1)))
    return T

def make_field(N, seed):
    r = np.random.default_rng(seed)
    T = su_basis(N)
    modes = [(r.normal(size=2), r.normal(size=len(T)), r.normal(size=len(T)), r.uniform(0, 6.3)) for _ in range(4)]
    def A(mu, x):
        out = np.zeros((N, N), complex)
        for kv, c0, c1, ph in modes:
            s = math.sin(kv @ x + ph + mu)
            cc = c0 if mu == 0 else c1
            for t, coef in zip(T, cc):
                out += coef * s * t
        return out
    return A, T

def deriv(A, mu, nu, x, h=1e-5):
    e = np.zeros(2); e[nu] = h
    return (A(mu, x + e) - A(mu, x - e)) / (2 * h)

def plaquette(A, x0, eps, k, sub=8):
    """Chain of link values around the square centred at x0 (bottom, right, top back, left back)."""
    vals = []
    def edge(mu, start, sign):
        d = eps / sub
        e = np.zeros(2); e[mu] = 1
        for s in range(sub):
            mid = start + e * d * (s + 0.5) * sign
            L = expm(k * sign * d * A(mu, mid))
            vals.append((L - np.eye(L.shape[0])) / k)
    h = eps / 2
    edge(0, x0 + np.array([-h, -h]), +1)
    edge(1, x0 + np.array([h, -h]), +1)
    edge(0, x0 + np.array([h, h]), -1)
    edge(1, x0 + np.array([-h, h]), -1)
    return chain(vals, k)

gauge = {}
for N, name in ((1, "U(1) electromagnetism"), (2, "SU(2) weak"), (3, "SU(3) strong")):
    A, T = make_field(N, 7 + N)
    k = 1j * 0.9           # k = i g : imaginary dial, unitary links
    x0 = np.array([0.3, -0.2])
    dA = deriv(A, 1, 0, x0) - deriv(A, 0, 1, x0)
    comm = k * (A(0, x0) @ A(1, x0) - A(1, x0) @ A(0, x0))
    F_full = dA + comm
    row = []
    for eps in (0.08, 0.04, 0.02, 0.01):
        # the loop is based at a corner, so it reads F there, turned by the corner's frame: an error ~ eps.
        # Two sizes (Richardson) remove it and leave eps^2.
        est = 2 * plaquette(A, x0, eps / 2, k) / (eps / 2)**2 - plaquette(A, x0, eps, k) / eps**2
        raw = plaquette(A, x0, eps, k) / eps**2
        row.append((eps, float(np.abs(est - F_full).max()), float(np.abs(est - dA).max()), float(np.abs(raw - F_full).max())))
    gauge[name] = {"rows": row, "comm_size": float(np.abs(comm).max()), "F_size": float(np.abs(F_full).max())}
    e_last, err_full, err_nocomm, _ = row[-1]
    order = math.log(row[0][1] / row[-1][1]) / math.log(row[0][0] / row[-1][0])
    check(f"B1 {name}: plaquette/eps^2 -> dA - dA + k[A,A], two loop sizes (error order {order:.2f})",
          err_full < 2e-3 and order > 1.8, f"err {err_full:.1e} at eps {e_last}")
    if N == 1:
        check(f"B2 {name}: no self-interaction term (values commute)", gauge[name]["comm_size"] < 1e-12)
    else:
        check(f"B2 {name}: leaving out k[A,A] misses by an order-one amount", err_nocomm > 0.05,
              f"miss {err_nocomm:.2f} vs F size {gauge[name]['F_size']:.2f}")
    # unitary leftover: 1 + k p is in the group
    p = plaquette(A, x0, 0.04, k)
    U = np.eye(N) + k * p
    check(f"B3 {name}: loop lift stays unitary (det 1 for SU)", np.abs(U.conj().T @ U - np.eye(N)).max() < 1e-12
          and (N == 1 or abs(np.linalg.det(U) - 1) < 1e-12))
results["gauge"] = gauge

# colour factors from the commutator (size of the k*ab term)
print("\n== B4. Size of the self-interaction: adjoint Casimir C_A from commutators")
CA = {}
for N in (1, 2, 3):
    T = su_basis(N)
    if N == 1:
        CA[N] = 0.0; continue
    D = len(T)
    f = np.zeros((D, D, D))
    for a_, b_, c_ in itertools.product(range(D), repeat=3):
        f[a_, b_, c_] = (-2j * np.trace((T[a_] @ T[b_] - T[b_] @ T[a_]) @ T[c_])).real
    M = np.einsum("acd,bcd->ab", f, f)
    CA[N] = float(M[0, 0])
    check(f"B4 SU({N}): f^acd f^bcd = C_A delta^ab with C_A = {N}", np.abs(M - N * np.eye(D)).max() < 1e-12)
CF3 = float(sum((t @ t) for t in su_basis(3))[0, 0].real)
check("B5 SU(3): quark charge C_F = 4/3", abs(CF3 - 4 / 3) < 1e-12)
meas = {"ALEPH": (2.93, math.hypot(0.14, 0.58), 1.35, math.hypot(0.07, 0.26)),
        "OPAL": (3.02, math.hypot(0.25, 0.49), 1.34, math.hypot(0.13, 0.22))}
for exp, (ca, dca, cf, dcf) in meas.items():
    check(f"B6 {exp} four-jet C_A = {ca} +- {dca:.2f}: within 1 sigma of 3, {ca/dca:.1f} sigma from 0 (abelian)",
          abs(ca - 3) < dca and ca / dca > 4)
results["casimir"] = {"CA": CA, "CF_SU3": CF3, "measured": meas}

# ---------------------------------------------------------------- Part C: gravity
print("\n== C. Gravity: the same loop on the Lorentz (tangent-frame) ring, Schwarzschild")
t, r, th, ph, Ms = sp.symbols("t r theta phi M", positive=True)
X = [t, r, th, ph]
f_ = 1 - 2 * Ms / r
g = sp.diag(-f_, 1 / f_, r**2, r**2 * sp.sin(th)**2)
gi = g.inv()
Gam = [[[sp.simplify(sum(gi[a_, d] * (sp.diff(g[d, b_], X[c_]) + sp.diff(g[d, c_], X[b_]) - sp.diff(g[b_, c_], X[d]))
                         for d in range(4)) / 2) for c_ in range(4)] for b_ in range(4)] for a_ in range(4)]
# Riemann R^a_{b c d} = d_c G^a_{db} - d_d G^a_{cb} + G^a_{ce} G^e_{db} - G^a_{de} G^e_{cb}
Riem = [[[[sp.simplify(sp.diff(Gam[a_][d][b_], X[c_]) - sp.diff(Gam[a_][c_][b_], X[d])
                       + sum(Gam[a_][c_][e] * Gam[e][d][b_] - Gam[a_][d][e] * Gam[e][c_][b_] for e in range(4)))
           for d in range(4)] for c_ in range(4)] for b_ in range(4)] for a_ in range(4)]
Gam_f = sp.lambdify((t, r, th, ph, Ms), [[[Gam[a_][b_][c_] for c_ in range(4)] for b_ in range(4)] for a_ in range(4)])
Riem_f = sp.lambdify((t, r, th, ph, Ms), Riem)
g_f = sp.lambdify((t, r, th, ph, Ms), g)

def conn(x, mu, M=1.0):
    G = np.array(Gam_f(*x, M), dtype=float)   # G[a][b][c] = Gamma^a_{bc}
    return G[:, mu, :]                        # (Gamma_mu)^a_c

def holonomy(x0, mu, nu, eps, M=1.0, sub=12):
    """Parallel-transport chain around a square centred at x0 in the (mu,nu) plane.
    Link value of a step = expm(-Gamma_mu d) - 1, k = 1; later steps sit to the LEFT."""
    vals = []
    def edge(dirn, start, sign, size):
        d = size / sub
        for s in range(sub):
            mid = start.copy(); mid[dirn] += sign * d * (s + 0.5)
            vals.append(expm(-conn(mid, dirn, M) * sign * d) - np.eye(4))
    em, en = eps[mu], eps[nu]
    c = np.array(x0, float)
    p0 = c.copy(); p0[mu] -= em / 2; p0[nu] -= en / 2
    p1 = p0.copy(); p1[mu] += em
    p2 = p1.copy(); p2[nu] += en
    p3 = p2.copy(); p3[mu] -= em
    edge(mu, p0, +1, em); edge(nu, p1, +1, en); edge(mu, p2, -1, em); edge(nu, p3, -1, en)
    out = np.zeros((4, 4))
    for v in vals:              # later step on the left: out = v o out
        out = v + out + v @ out
    return out

x0 = (0.0, 6.0, 1.0, 0.3)
R0 = np.array(Riem_f(*x0, 1.0), dtype=float)
scale = (2e-3, 2e-3, 2e-3, 2e-3)
worst, sgn = 0.0, None
holo_store = {}
for mu, nu in itertools.combinations(range(4), 2):
    H = 2 * holonomy(x0, mu, nu, [e / 2 for e in scale]) * 4 - holonomy(x0, mu, nu, scale)
    H = H / 1.0
    pred = R0[:, :, mu, nu] * scale[mu] * scale[nu]
    s = -1.0 if np.abs(H + pred).max() < np.abs(H - pred).max() else 1.0
    sgn = s if sgn is None else sgn
    worst = max(worst, np.abs(H - s * pred).max() / (scale[mu] * scale[nu]))
    holo_store[(mu, nu)] = H
check("C1 all 6 planes x 16 entries: loop leftover/area = -R^a_{b mu nu} (one sign for all)",
      worst < 1e-5 and sgn == -1.0, f"max error {worst:.1e} (Riemann entries up to {np.abs(R0).max():.2e})")
iso = 0.0
for (mu, nu) in holo_store:
    Hraw = holonomy(x0, mu, nu, scale)
    pb = np.array(x0, float); pb[mu] -= scale[mu] / 2; pb[nu] -= scale[nu] / 2   # base corner
    gb = np.array(g_f(*pb, 1.0), dtype=float)
    iso = max(iso, np.abs((np.eye(4) + Hraw).T @ gb @ (np.eye(4) + Hraw) - gb).max())
check("C2 every loop lift keeps the metric at the base point (it is a Lorentz transformation)", iso < 1e-10, f"{iso:.1e}")
Htr = np.eye(4) + holo_store[(0, 1)]; Hthph = np.eye(4) + holo_store[(2, 3)]
boost_like = abs(Htr[0, 1]) > 0 and abs(Htr[2, 3]) < 1e-14
rot_like = abs(Hthph[2, 3]) > 0 and abs(Hthph[0, 1]) < 1e-14
check("C3 a loop in the t-r plane returns a boost; in the theta-phi plane, a rotation", boost_like and rot_like)

# sectional curvatures from loops vs exact, over r
names = {(0, 1): "t-r", (0, 2): "t-theta", (0, 3): "t-phi", (1, 2): "r-theta", (1, 3): "r-phi", (2, 3): "theta-phi"}
sec = {v: [] for v in names.values()}
rs = [2.6, 3.5, 5.0, 7.0, 10.0, 14.0, 20.0]
for rr in rs:
    xx = (0.0, rr, 1.0, 0.3)
    gg = np.array(g_f(*xx, 1.0), dtype=float)
    for (mu, nu), nm in names.items():
        ee = [2e-3 * max(1.0, rr / 6)] * 4
        H = 2 * 4 * holonomy(xx, mu, nu, [e / 2 for e in ee]) - holonomy(xx, mu, nu, ee)
        Rup = -H[mu, nu] / (ee[mu] * ee[nu])          # R^mu_{nu mu nu}
        K = gg[mu, mu] * Rup / (gg[mu, mu] * gg[nu, nu])
        Kex = float(sp.N(sum(g[mu, a_] * Riem[a_][nu][mu][nu] for a_ in range(4)) / (g[mu, mu] * g[nu, nu])
                         ).subs({t: 0, r: rr, th: 1.0, ph: 0.3, Ms: 1.0}))
        sec[nm].append((rr, K, Kex))
werr = max(abs(K - Kex) * rr**3 for v in sec.values() for rr, K, Kex in v)
check("C4 six sectional curvatures from loops match exact, 2.6 <= r <= 20", werr < 1e-4, f"max error x r^3 = {werr:.1e}")
exact_forms = {"t-r": 2, "theta-phi": 2, "t-theta": -1, "t-phi": -1, "r-theta": -1, "r-phi": -1}
ok = all(abs(Kex * rr**3 / Mfac - 1) < 1e-9 for nm, Mfac in exact_forms.items() for rr, K, Kex in sec[nm])
check("C5 the loops read the mass directly: K(t-r) = K(theta-phi) = +2M/r^3, the other four planes -M/r^3",
      ok, "K = R_{mu nu mu nu}/(g_mu mu g_nu nu), signature -+++")
ric = max(abs(sum(dict((rr2, K) for rr2, K, _ in sec[nm])[rr] for nm in planes))
          for rr in rs for planes in (("t-r", "t-theta", "t-phi"), ("t-r", "r-theta", "r-phi")))
check("C6 summing the planes through one direction gives zero (Ricci = 0 outside the mass)", ric < 1e-5, f"{ric:.1e}")
results["schwarzschild_sectional"] = {k_: [[a, b, c] for a, b, c in v] for k_, v in sec.items()}
results["riemann_count"] = {"2D": [1, 1, 0], "3D": [6, 6, 0], "4D": [20, 10, 10]}
check("C7 components: Riemann = Ricci + Weyl gives 1=1+0 (2D), 6=6+0 (3D), 20=10+10 (4D)",
      all(n_ ** 2 * (n_ ** 2 - 1) // 12 == R + W for n_, (_, R, W) in zip((2, 3, 4), results["riemann_count"].values()))
      and all(n_ * (n_ + 1) // 2 == R for n_, (_, R, W) in zip((2, 3, 4), [(0, 1, 0), (0, 6, 0), (0, 10, 0)]) if n_ > 2) )

# ---------------------------------------------------------------- Part D: Gravity Probe B by one 4D chain
print("\n== D. Gravity Probe B: one chain of identical link values around the orbit")
mp.mp.dps = 50
GM, c_ = mp.mpf("3.986004418e14"), mp.mpf(299792458)
Mg = GM / c_**2
rg = mp.mpf("7027.4e3")
Om = mp.sqrt(Mg / rg**3)                 # d phi / d(ct)
fv = 1 - 2 * Mg / rg
Gs = [[[sp.N(Gam[a_][b_][c_].subs({r: sp.Float(str(rg), 50), th: sp.pi / 2, Ms: sp.Float(str(Mg), 50)}), 50)
        for c_ in range(4)] for b_ in range(4)] for a_ in range(4)]
Amat = mp.matrix(4, 4)
for a_ in range(4):
    for c2 in range(4):
        Amat[a_, c2] = mp.mpf(str(Gs[a_][0][c2])) / Om + mp.mpf(str(Gs[a_][3][c2]))   # Gamma_t dt/dphi + Gamma_phi
Hexact = mp.expm(-2 * mp.pi * Amat)
def angle_after(Hm):
    S0 = mp.matrix([0, 1 / mp.sqrt(1 / fv), 0, 0])          # unit radial, orthogonal to u
    S = Hm * S0
    gD = [-fv, 1 / fv, rg**2, rg**2]
    u = mp.matrix([1, 0, 0, Om]); un = mp.sqrt(-(gD[0] * u[0]**2 + gD[3] * u[3]**2)); u = u / un
    # unit vector along phi, orthogonal to u
    w = mp.matrix([Om * gD[3] / (-gD[0]) * 1, 0, 0, 1])     # w.u = 0
    wn = mp.sqrt(gD[0] * w[0]**2 + gD[3] * w[3]**2); w = w / wn
    dot = lambda a, b: sum(gD[i] * a[i] * b[i] for i in range(4))
    return mp.atan2(dot(S, w), dot(S, S0))
ang = angle_after(Hexact)
yr, mas = mp.mpf(365.25 * 86400), 180 / mp.pi * 3600 * 1000
period = 2 * mp.pi / (Om * c_)
rate = ang / period * yr * mas
closed = (1 - mp.sqrt(1 - 3 * Mg / rg)) * 2 * mp.pi / period * yr * mas
check("D1 chain of n identical steps -> expm: angle per orbit = 2 pi (1 - sqrt(1-3M/r)) exactly",
      abs(rate - closed) < mp.mpf("1e-20"), f"{float(rate):.2f} mas/yr")
conv = []
# Finite chains: n identical Cayley link values a = (1 - A d/2)^-1 (A d) ... lift (1 - A d/2)^-1 (1 + A d/2), d = -2 pi / n.
# (A first-order link 1 - A d is not a Lorentz step, and its error swamps a 6e-9 rad effect.)
for m in (10, 14, 18, 22, 26):
    nsteps = 2**m
    d = -2 * mp.pi / nsteps
    Hn = mp.inverse(mp.eye(4) - Amat * (d / 2)) * (mp.eye(4) + Amat * (d / 2))
    for _ in range(m):
        Hn = Hn * Hn                     # composing the chain with itself: 2^m identical steps
    conv.append((nsteps, float(angle_after(Hn) / period * yr * mas)))
errs_d = [abs(b - float(rate)) for a, b in conv]
check("D2 finite chains of 2^m identical Lorentz steps converge to it like 1/n^2",
      errs_d[-1] < 0.05 and errs_d[0] / errs_d[1] > 200, ", ".join(f"n=2^{int(math.log2(a))}: {b:.1f}" for a, b in conv))
gpb_meas, gpb_err = 6601.8, 18.3
check(f"D3 against Gravity Probe B: {float(rate):.1f} vs {gpb_meas} +- {gpb_err} mas/yr",
      abs(float(rate) - gpb_meas) < gpb_err, f"{(float(rate)-gpb_meas)/gpb_err:+.2f} sigma")
results["gpb"] = {"rate": float(rate), "convergence": conv, "measured": [gpb_meas, gpb_err], "GR_with_oblateness": 6606.1}

# ---------------------------------------------------------------- Part E: Aharonov-Bohm, a U(1) chain round a solenoid
print("\n== E. Aharonov-Bohm: the U(1) chain round a solenoid, field zero on the path")
def ab_phase(flux, path, q_over_hbar=1.0, nstep=4000):
    """A = flux/(2 pi rho) phi-hat outside; link value (e^{i q A.dl/hbar} - 1)/k with k = i."""
    k = 1j
    total = 0j
    ts = np.linspace(0, 1, nstep + 1)
    pts = np.array([path(s) for s in ts])
    for p0, p1 in zip(pts[:-1], pts[1:]):
        mid = (p0 + p1) / 2; dl = p1 - p0
        rho2 = mid @ mid
        Avec = flux / (2 * math.pi * rho2) * np.array([-mid[1], mid[0]])
        L = np.exp(k * q_over_hbar * (Avec @ dl))
        a = (L - 1) / k
        total = total + a + k * total * a
    return np.angle(1 + k * total)
circle = lambda s: 2.0 * np.array([math.cos(2 * math.pi * s), math.sin(2 * math.pi * s)])
ellipse = lambda s: np.array([3.0 * math.cos(2 * math.pi * s) + 0.5, 1.5 * math.sin(2 * math.pi * s)])
phis = []
for fl in (0.5, 1.3, 2.9):
    p1, p2 = ab_phase(fl, circle), ab_phase(fl, ellipse)
    phis.append((fl, p1, p2))
okab = all(abs(np.angle(np.exp(1j * (p1 - fl)))) < 1e-5 and abs(np.angle(np.exp(1j * (p2 - fl)))) < 1e-5 for fl, p1, p2 in phis)
check("E1 loop phase = q*flux/hbar, same for a circle and an off-centre ellipse", okab)
outside = lambda s: np.array([5.0 + 0.8 * math.cos(2 * math.pi * s), 0.8 * math.sin(2 * math.pi * s)])
check("E2 a loop that does not enclose the solenoid gets no phase", abs(ab_phase(2.0, outside)) < 1e-6)
# electrons (q = e) round flux h/2e: phase pi.  Tonomura et al. 1986 saw pi (and 0 when the trapped flux was even).
check("E3 electrons round one superconducting flux quantum h/2e: phase = pi (Tonomura 1986 saw a half-fringe shift)",
      abs(abs(ab_phase(math.pi, circle)) - math.pi) < 1e-6)
results["ab"] = phis

# ---------------------------------------------------------------- Part F: twist (torsion) next to curvature
print("\n== F. Twist: one 5x5 chain (turn + shift) gives curvature AND torsion together")
# Weitzenbock (teleparallel) connection from the static tetrad e^0 = sqrt(f) dt, e^1 = dr/sqrt(f), e^2 = r dtheta, e^3 = r sin dphi
E = sp.diag(sp.sqrt(f_), 1 / sp.sqrt(f_), r, r * sp.sin(th))      # E[a, nu] = e^a_nu
Einv = E.inv()                                                      # Einv[rho, a] = e_a^rho
GamW = [[[sp.simplify(sum(Einv[a_, q] * sp.diff(E[q, c_], X[b_]) for q in range(4)))
          for c_ in range(4)] for b_ in range(4)] for a_ in range(4)]   # Gamma^a_{b c}, derivative index b
GamW_f = sp.lambdify((t, r, th, ph, Ms), GamW)
Tors = [[[sp.simplify(GamW[a_][b_][c_] - GamW[a_][c_][b_]) for c_ in range(4)] for b_ in range(4)] for a_ in range(4)]
Tors_f = sp.lambdify((t, r, th, ph, Ms), Tors)

def conn5(x, mu, which, M=1.0):
    G = np.array((Gam_f if which == "LC" else GamW_f)(*x, M), dtype=float)
    A = np.zeros((5, 5)); A[:4, :4] = G[:, mu, :]; A[mu, 4] = -1.0    # turn block + solder (shift) column
    return A

def holonomy5(x0, mu, nu, eps, which, M=1.0, sub=12):
    vals = []
    def edge(dirn, start, sign, size):
        d = size / sub
        for s_ in range(sub):
            mid = start.copy(); mid[dirn] += sign * d * (s_ + 0.5)
            vals.append(expm(-conn5(mid, dirn, which, M) * sign * d) - np.eye(5))
    em, en = eps, eps
    p0 = np.array(x0, float); p0[mu] -= em / 2; p0[nu] -= en / 2
    p1 = p0.copy(); p1[mu] += em; p2 = p1.copy(); p2[nu] += en; p3 = p2.copy(); p3[mu] -= em
    edge(mu, p0, +1, em); edge(nu, p1, +1, en); edge(mu, p2, -1, em); edge(nu, p3, -1, en)
    out = np.zeros((5, 5))
    for v in vals:
        out = v + out + v @ out
    return out

def rich5(x0, mu, nu, eps, which, M=1.0):
    return 8 * holonomy5(x0, mu, nu, eps / 2, which, M) - holonomy5(x0, mu, nu, eps, which, M)

xF = (0.0, 6.0, 1.0, 0.3); epsF = 2e-3
R0F = np.array(Riem_f(*xF, 1.0), dtype=float); T0 = np.array(Tors_f(*xF, 1.0), dtype=float)
lc_turn = lc_shift = w_turn = w_shift_err = 0.0
twist_store = {}
for mu, nu in itertools.combinations(range(4), 2):
    Hl = rich5(xF, mu, nu, epsF, "LC") / epsF**2
    Hw = rich5(xF, mu, nu, epsF, "W") / epsF**2
    lc_turn = max(lc_turn, np.abs(Hl[:4, :4] + R0F[:, :, mu, nu]).max())
    lc_shift = max(lc_shift, np.abs(Hl[:4, 4]).max())
    w_turn = max(w_turn, np.abs(Hw[:4, :4]).max())
    w_shift_err = max(w_shift_err, np.abs(Hw[:4, 4] - T0[:, mu, nu]).max())
    twist_store[names[(mu, nu)]] = {"LC_turn": float(np.abs(Hl[:4, :4]).max()), "LC_shift": float(np.abs(Hl[:4, 4]).max()),
                                    "W_turn": float(np.abs(Hw[:4, :4]).max()), "W_shift": float(np.abs(Hw[:4, 4]).max())}
check("F1 Levi-Civita (GR): turn block = -R, shift column = 0 (no twist)", lc_turn < 1e-5 and lc_shift < 1e-6,
      f"turn err {lc_turn:.1e}, shift {lc_shift:.1e}")
check("F2 teleparallel (Weitzenbock): turn block = 0 (no curvature), shift column = +T (all twist)",
      w_turn < 1e-6 and w_shift_err < 1e-5, f"turn {w_turn:.1e}, shift err {w_shift_err:.1e}")
# the mass sits in the t-component of the twist: T^t_{r t} = hover acceleration / redshift
rr_ = sp.symbols("rr", positive=True)
Ttrt = sp.simplify(Tors[0][1][0])
check("F3 the twist that carries the mass: T^t_{rt} = M/(r^2 (1-2M/r)) = hovering acceleration M/(r^2 sqrt f) over sqrt f",
      sp.simplify(Ttrt - Ms / (r**2 * f_)) == 0, f"T^t_rt = {Ttrt}")
# same metric, same physics: the contorsion turns one connection into the other
subF = {t: 0, r: 6, th: 1.0, ph: 0.3, Ms: 1}
Tn = np.array([[[float(Tors[a_][b_][c_].subs(subF)) for c_ in range(4)] for b_ in range(4)] for a_ in range(4)])
Wn = np.array([[[float(GamW[a_][b_][c_].subs(subF)) for c_ in range(4)] for b_ in range(4)] for a_ in range(4)])
Ln = np.array([[[float(Gam[a_][b_][c_].subs(subF)) for c_ in range(4)] for b_ in range(4)] for a_ in range(4)])
gn = np.array(g.subs(subF), dtype=float)
Tl = np.einsum("la,abc->lbc", gn, Tn)                    # T_{a b c}, b = derivative slot
Kl = 0.5 * (Tl + np.transpose(Tl, (1, 0, 2)) + np.transpose(Tl, (1, 2, 0)))   # K_{abc} = (T_abc + T_bac + T_cab)/2
recon = np.abs(np.einsum("la,abc->lbc", np.linalg.inv(gn), Kl) - (Wn - Ln)).max()
check("F4 trade: Weitzenbock = Levi-Civita + contorsion K_abc = (T_abc + T_bac + T_cab)/2, built from the twist alone", recon < 1e-12, f"{recon:.1e}")
# torsion vector and the scalar: TEGR's Lagrangian T equals -R up to a total derivative; R = 0 here, so T = -2 div(T_vector)
Tvec = [sum(Tors[a_][a_][m_] for a_ in range(4)) for m_ in range(4)]   # T^a_{a mu}
check("F5 twist with M = 0 in these round frames is not zero (frame artefact): T^theta_{r theta} = 1/r",
      sp.simplify(Tors[2][1][2].subs(Ms, 0) - 1 / r) == 0)
abel = max(abs(sp.simplify(sum(E[a_, q] * Tors[q][b_][c_] for q in range(4)) - (sp.diff(E[a_, c_], X[b_]) - sp.diff(E[a_, b_], X[c_]))))
           for a_ in range(4) for b_ in range(4) for c_ in range(4))
check("F6 in frame components the twist is T^a = d e^a, linear like electromagnetism's F = dA (no k*ab term)", abel == 0)
results["twist"] = twist_store

# ---------------------------------------------------------------- Part G: frame dragging as the twist of the time leg
print("\n== G. Frame dragging = twist of the static observers' time leg (de^0), Gravity Probe B")
Gc, cc = 6.674e-11, 299792458.0
I_E, w_E = 8.034e37, 7.2921159e-5            # Earth's moment of inertia (kg m^2), spin rate (rad/s)
J = I_E * w_E * np.array([0, 0, 1.0])
GMe = 3.986004418e14
def g0i(x):                                   # weak-field rotating Earth, x0 = ct: g_0i = -(2G/c^3) (J x r)_i / r^3
    rr_ = np.linalg.norm(x)
    return -(2 * Gc / cc**3) * np.cross(J, x) / rr_**3
def g00(x):
    return -(1 - 2 * GMe / (cc**2 * np.linalg.norm(x)))
def e0(x):                                    # time leg of static observers: e^0_i = g_0i / sqrt(-g_00)
    return g0i(x) / math.sqrt(-g00(x))
def twist_time_leg(x, h=1.0):
    """T^0_{jk} = d_j e^0_k - d_k e^0_j, by central differences (step h metres); returned as an axial vector."""
    D = np.zeros((3, 3))
    for j in range(3):
        dx = np.zeros(3); dx[j] = h
        D[j] = (e0(x + dx) - e0(x - dx)) / (2 * h)       # D[j,k] = d_j e^0_k
    T = D - D.T
    return np.array([T[1, 2], T[2, 0], T[0, 1]])            # (T_yz, T_zx, T_xy)
# check against the closed form at one point
xp = np.array([3.1e6, -2.2e6, 5.4e6])
Omega_twist = -(cc / 2) * twist_time_leg(xp)
nvec = xp / np.linalg.norm(xp)
Omega_LT = Gc / (cc**2 * np.linalg.norm(xp)**3) * (3 * (J @ nvec) * nvec - J)
rel = np.abs(Omega_twist - Omega_LT).max() / np.abs(Omega_LT).max()
check("G1 gyroscope drift = -(c/2) x twist of the time leg = Lense-Thirring (G/c^2 r^3)(3(J.n)n - J)", rel < 1e-5, f"rel err {rel:.1e}")
# average round GP-B's polar orbit
rgp = 7027.4e3
ths = np.linspace(0, 2 * np.pi, 721)[:-1]
Om_avg = np.mean([-(cc / 2) * twist_time_leg(rgp * np.array([math.sin(a), 0, math.cos(a)]), h=2.0) for a in ths], axis=0)
masyr = 365.25 * 86400 * 180 / math.pi * 3600e3
fd = Om_avg[2] * masyr
closed_fd = Gc * np.linalg.norm(J) / (2 * cc**2 * rgp**3) * masyr
check("G2 polar-orbit average of the time-leg twist = G J / (2 c^2 r^3)", abs(fd - closed_fd) / closed_fd < 1e-4,
      f"{fd:.2f} mas/yr (closed form {closed_fd:.2f})")
check(f"G3 against Gravity Probe B frame dragging: {fd:.1f} vs 37.2 +- 7.2 mas/yr (GP-B's own GR figure 39.2)",
      abs(fd - 37.2) < 7.2, f"{(fd-37.2)/7.2:+.2f} sigma")
# the other twist pieces carry no frame dragging: drop g_0i and the time-leg twist is zero
J_saved = J.copy(); J = np.zeros(3)
zero_tw = np.abs(twist_time_leg(xp)).max()
J = J_saved
check("G4 without the Earth's spin the time leg has no twist (Newton's pull sits in T^0_0i, not here)", zero_tw < 1e-30, f"{zero_tw:.1e}")
results["frame_drag_twist"] = {"mas_per_yr": fd, "closed": closed_fd, "measured": [37.2, 7.2], "GPB_GR": 39.2,
                               "orbit_profile": [[float(a), float((-(cc / 2) * twist_time_leg(rgp * np.array([math.sin(a), 0, math.cos(a)]), h=2.0))[0] * masyr),
                                                  float((-(cc / 2) * twist_time_leg(rgp * np.array([math.sin(a), 0, math.cos(a)]), h=2.0))[2] * masyr)] for a in ths[::10]]}

# ---------------------------------------------------------------- Part H: one chain for gravity + light (Kaluza-Klein, 5D)
print("\n== H. One 5D chain carries both gravity and light (Kaluza 1921, Klein 1926)")
T5, X5, Y5, Z5, W5, Bf, kap = sp.symbols("t x y z w B kappa", real=True)
XX = [T5, X5, Y5, Z5, W5]
Avec4 = [0, -Bf * Y5 / 2, Bf * X5 / 2, 0]                    # uniform magnetic field B along z
eta = sp.diag(-1, 1, 1, 1)
g5 = sp.zeros(5, 5)
for m_ in range(4):
    for n_ in range(4):
        g5[m_, n_] = eta[m_, n_] + kap**2 * Avec4[m_] * Avec4[n_]
    g5[m_, 4] = g5[4, m_] = kap * Avec4[m_]
g5[4, 4] = 1
g5i = sp.simplify(g5.inv())
G5 = [[[sp.simplify(sum(g5i[a_, d] * (sp.diff(g5[d, b_], XX[c_]) + sp.diff(g5[d, c_], XX[b_]) - sp.diff(g5[b_, c_], XX[d]))
                        for d in range(5)) / 2) for c_ in range(5)] for b_ in range(5)] for a_ in range(5)]
R5 = [[[[sp.simplify(sp.diff(G5[a_][d][b_], XX[c_]) - sp.diff(G5[a_][c_][b_], XX[d])
                     + sum(G5[a_][c_][e] * G5[e][d][b_] - G5[a_][d][e] * G5[e][c_][b_] for e in range(5)))
          for d in range(5)] for c_ in range(5)] for b_ in range(5)] for a_ in range(5)]
Ric5 = sp.simplify(sum(g5i[b_, d] * R5[a_][b_][a_][d] for a_ in range(5) for b_ in range(5) for d in range(5)))
F2 = 2 * Bf**2                                                 # F_{mu nu} F^{mu nu} for a pure magnetic field
check("H1 5D curvature scalar = -(kappa^2/4) F^2 with flat 4D: light appears as curvature of the fifth direction",
      sp.simplify(Ric5 + kap**2 * F2 / 4) == 0, f"R5 = {Ric5}")
G5_f = sp.lambdify(XX + [Bf, kap], G5); R5_f = sp.lambdify(XX + [Bf, kap], R5)
def holonomyN(x0, mu, nu, eps, Gf, pars, N, sub=12):
    def conn_(x, d):
        return np.array(Gf(*x, *pars), dtype=float)[:, d, :]
    vals = []
    def edge(dirn, start, sign, size):
        dd = size / sub
        for s_ in range(sub):
            mid = start.copy(); mid[dirn] += sign * dd * (s_ + 0.5)
            vals.append(expm(-conn_(mid, dirn) * sign * dd) - np.eye(N))
    p0 = np.array(x0, float); p0[mu] -= eps / 2; p0[nu] -= eps / 2
    p1 = p0.copy(); p1[mu] += eps; p2 = p1.copy(); p2[nu] += eps; p3 = p2.copy(); p3[mu] -= eps
    edge(mu, p0, +1, eps); edge(nu, p1, +1, eps); edge(mu, p2, -1, eps); edge(nu, p3, -1, eps)
    out = np.zeros((N, N))
    for v in vals:
        out = v + out + v @ out
    return out
xK = (0.0, 0.4, -0.7, 0.2, 0.0); parsK = (0.6, 1.0); eK = 2e-3
RK = np.array(R5_f(*xK, *parsK), dtype=float)
worstK = 0.0; kk_store = {}
lab5 = ["t", "x", "y", "z", "w"]
for mu, nu in itertools.combinations(range(5), 2):
    H = (8 * holonomyN(xK, mu, nu, eK / 2, G5_f, parsK, 5) - holonomyN(xK, mu, nu, eK, G5_f, parsK, 5)) / eK**2
    worstK = max(worstK, np.abs(H + RK[:, :, mu, nu]).max())
    kk_store[lab5[mu] + "-" + lab5[nu]] = float(np.abs(H).max())
check("H2 all 10 planes x 25 entries of the 5D curvature come out of loops of one 5x5 chain", worstK < 1e-5, f"max error {worstK:.1e}")
g5n = np.array(g5.subs({X5: xK[1], Y5: xK[2], Bf: parsK[0], kap: parsK[1]}), dtype=float)
# contract from the loop-measured tensor (rebuild it once)
Rloop = np.zeros((5, 5, 5, 5))
for mu, nu in itertools.combinations(range(5), 2):
    H = (8 * holonomyN(xK, mu, nu, eK / 2, G5_f, parsK, 5) - holonomyN(xK, mu, nu, eK, G5_f, parsK, 5)) / eK**2
    Rloop[:, :, mu, nu] = -H; Rloop[:, :, nu, mu] = H
g5in = np.linalg.inv(g5n)
R5scal_loops = np.einsum("bd,abad->", g5in, Rloop)
check("H3 contracting the loop-measured curvature gives R5 = -B^2/2: the field's energy, read off loops",
      abs(R5scal_loops + parsK[0]**2 / 2) < 1e-5, f"{R5scal_loops:.7f} vs {-parsK[0]**2/2}")
# a geodesic in 5D is a charged particle in 4D: cyclotron circle, omega = kappa * p_w * B
Gfun = lambda x: np.array(G5_f(*x, *parsK), dtype=float)
def geo(y, h):
    x, u = y[:5], y[5:]
    a = -np.einsum("abc,b,c->a", Gfun(x), u, u)
    return np.concatenate([u, a])
pw = 0.8                                                      # conserved momentum along the fifth direction
x_ = np.array([0, 0.5, 0, 0, 0.0]); uxy = np.array([0, 0.3])
u_ = np.array([0, 0, 0.3, 0, 0.0]); u_[4] = pw - kap.subs(kap, parsK[1]) * float(Avec4[2].subs({Bf: parsK[0], X5: 0.5})) * u_[2]
u_[0] = math.sqrt(1 + 0.3**2)
y = np.concatenate([x_, u_]); h = 1e-3; traj = [y[:5].copy()]
for i in range(int(2 * math.pi / (parsK[1] * pw * parsK[0]) / h * 1.0)):
    k1 = geo(y, h); k2 = geo(y + h / 2 * k1, h); k3 = geo(y + h / 2 * k2, h); k4 = geo(y + h * k3, h)
    y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); traj.append(y[:5].copy())
traj = np.array(traj)
closure = np.linalg.norm(traj[-1, 1:3] - traj[0, 1:3])
radius = 0.3 / (parsK[1] * pw * parsK[0])
rfit = (traj[:, 1].max() - traj[:, 1].min()) / 2
check("H4 a straight line in 5D is a cyclotron circle in 4D: closes after 2 pi/(kappa p_w B), radius v/(kappa p_w B)",
      closure < 1e-3 and abs(rfit - radius) / radius < 1e-3, f"gap {closure:.1e}, radius {rfit:.4f} vs {radius:.4f}")
results["kk"] = {"planes": kk_store, "R5_loops": float(R5scal_loops), "traj": traj[::40, 1:3].tolist(), "radius": radius}

# ---------------------------------------------------------------- Part I: the weak force's group from a hidden 3-sphere
print("\n== I. Non-abelian Kaluza-Klein: SU(2) from a hidden 3-sphere; the self-interaction is the sphere's twist")
# Everything numeric: S^3 = SU(2) in Euler angles y = (phi, theta, psi), g = e^{i phi s3/2} e^{i theta s2/2} e^{i psi s3/2}.
P1 = np.array([[0, 1], [1, 0]], complex); P2 = np.array([[0, -1j], [1j, 0]]); P3 = np.array([[1, 0], [0, -1]], complex)
def Ug_(y):
    e = lambda sig, a: math.cos(a / 2) * np.eye(2) + 1j * math.sin(a / 2) * sig
    return e(P3, y[0]) @ e(P2, y[1]) @ e(P3, y[2])
_ph, _th, _ps = sp.symbols("phi1 theta1 psi1", real=True)
# right-invariant forms of SU(2) in these Euler angles: dg g^-1 = (i/2) sigma_a rho^a (derived with sympy, then hard-coded)
RHO = sp.Matrix([[0, sp.sin(_ph), -sp.cos(_ph) * sp.sin(_th)],
                 [0, sp.cos(_ph), sp.sin(_ph) * sp.sin(_th)],
                 [1, 0, sp.cos(_th)]])
RHO_f = sp.lambdify((_ph, _th, _ps), RHO)
def rhoF_(y, h=None):
    return np.array(RHO_f(*y), dtype=float)
# check the hard-coded forms against dg g^-1 numerically
_yy = np.array([0.4, 1.1, 0.7]); _hh = 1e-6; _err = 0.0
for _m in range(3):
    _dy = np.zeros(3); _dy[_m] = _hh
    _M = (Ug_(_yy + _dy) - Ug_(_yy - _dy)) / (2 * _hh) @ np.linalg.inv(Ug_(_yy))
    for _a, _sg in enumerate((P1, P2, P3)):
        _err = max(_err, abs((-1j * np.trace(_sg @ _M)).real - rhoF_(_yy)[_a, _m]))
check("I0 hard-coded right-invariant forms match dg g^-1", _err < 1e-8, f"{_err:.1e}")
def Rvec_(y):
    return np.linalg.inv(rhoF_(y))                    # columns a: right-invariant vector field R_a^m
def gamS3(y, rho_):
    rf = rhoF_(y); return (rho_**2 / 4) * rf.T @ rf
y0 = np.array([0.4, 1.1, 0.7]); rhov = 1.3
KKn = Rvec_(y0).T @ gamS3(y0, rhov) @ Rvec_(y0)
check("I1 the sphere's symmetry fields have constant products R_a.R_b = (rho^2/4) delta_ab",
      np.abs(KKn - rhov**2 / 4 * np.eye(3)).max() < 1e-8, f"{np.abs(KKn - rhov**2/4*np.eye(3)).max():.1e}")
def lie_(a_, b_, y, h=1e-4):
    Va, Vb = Rvec_(y)[:, a_], Rvec_(y)[:, b_]
    dVb = np.zeros((3, 3)); dVa = np.zeros((3, 3))
    for n_ in range(3):
        dy = np.zeros(3); dy[n_] = h
        dVb[:, n_] = (Rvec_(y + dy)[:, b_] - Rvec_(y - dy)[:, b_]) / (2 * h)
        dVa[:, n_] = (Rvec_(y + dy)[:, a_] - Rvec_(y - dy)[:, a_]) / (2 * h)
    return dVb @ Va - dVa @ Vb
br = lie_(0, 1, y0); R3_ = Rvec_(y0)[:, 2]
cbr = float(br @ R3_ / (R3_ @ R3_))
check("I2 the sphere's symmetries don't commute: [R_1, R_2] = c R_3 with c = +-1 (SU(2)'s bracket)",
      np.abs(br - cbr * R3_).max() < 1e-6 and abs(abs(cbr) - 1) < 1e-6, f"c = {cbr:+.6f}")
c1v, c2v = 0.7, 0.5
Acon = np.array([[0, c1v, 0, 0], [0, 0, c2v, 0], [0, 0, 0, 0.0]])   # A^a_mu: constant, no derivatives at all
etan = np.diag([-1.0, 1, 1, 1])
def metric7(x, hidden):
    y = np.array(x[4:])
    if hidden == "sphere":
        Rv, gm = Rvec_(y), gamS3(y, rhov)
    else:
        Rv, gm = np.eye(3), np.eye(3) * rhov**2 / 4               # flat 3-torus: commuting symmetries
    Nm = Acon.T @ Rv.T                                            # N[mu, m] = A^a_mu R_a^m
    G = np.zeros((7, 7))
    G[:4, :4] = etan + Nm @ gm @ Nm.T
    G[:4, 4:] = -Nm @ gm; G[4:, :4] = G[:4, 4:].T
    G[4:, 4:] = gm
    return G
def christoffel_num(hidden, h=None):
    Xs = sp.symbols("x0:7", real=True)
    if hidden == "sphere":
        Rh = RHO.subs({_ph: Xs[4], _th: Xs[5], _ps: Xs[6]})
        Rv = sp.simplify(Rh.inv()); gm = (sp.Rational(rhov**2).limit_denominator(10**6) if False else rhov**2) / 4 * (Rh.T * Rh)
    else:
        Rv, gm = sp.eye(3), sp.eye(3) * rhov**2 / 4
    Am = sp.Matrix(Acon)
    Nm = Am.T * Rv.T
    G = sp.zeros(7, 7)
    G[:4, :4] = sp.Matrix(etan) + Nm * gm * Nm.T
    G[:4, 4:] = -Nm * gm; G[4:, :4] = G[:4, 4:].T; G[4:, 4:] = gm
    Gf = sp.lambdify(Xs, G)
    dGf = sp.lambdify(Xs, [[[sp.diff(G[i_, j_], Xs[l_]) for l_ in range(7)] for j_ in range(7)] for i_ in range(7)])
    def Gam(*x):
        Gm = np.array(Gf(*x), dtype=float); Gi = np.linalg.inv(Gm)
        dG = np.array(dGf(*x), dtype=float)                  # dG[i, j, l] = d_l G_ij
        lower = 0.5 * (dG + np.transpose(dG, (0, 2, 1)) - np.transpose(dG, (2, 0, 1)))   # [d, b, c]
        return np.einsum("ad,dbc->abc", Gi, lower)
    return Gam
x7 = (0.0, 0.2, -0.1, 0.3, 0.4, 1.1, 0.7)
def scalar_from_loops(hidden, x0, eps=2e-3, N=7):
    Gf = christoffel_num(hidden)
    Rl = np.zeros((N, N, N, N)); planes = {}
    for mu, nu in itertools.combinations(range(N), 2):
        H = (8 * holonomyN(x0, mu, nu, eps / 2, Gf, (), N, sub=10) - holonomyN(x0, mu, nu, eps, Gf, (), N, sub=10)) / eps**2
        Rl[:, :, mu, nu] = -H; Rl[:, :, nu, mu] = H; planes[(mu, nu)] = float(np.abs(H).max())
    gin = np.linalg.inv(metric7(np.array(x0), hidden))
    return float(np.einsum("bd,abad->", gin, Rl)), planes
R7_loops, planes7 = scalar_from_loops("sphere", x7)
F3 = c1v * c2v                                                     # |F^3_xy| = c1 c2: only the [A,A] term
pred7 = 6 / rhov**2 - (rhov**2 / 4) * (2 * F3**2) / 4
check("I3 7D curvature from loops = R(S^3) - (rho^2/16) F^a F^a, with F^3_xy = c1 c2 made ONLY of the [A,A] term",
      abs(R7_loops - pred7) < 1e-5, f"{R7_loops:.5f} vs {pred7:.5f}")
R7t, planes7t = scalar_from_loops("torus", x7)
maxleft = max(planes7t.values())
check("I4 same constant field, hidden torus (commuting symmetries): every loop closes, no field at all",
      maxleft < 1e-5, f"largest leftover/area {maxleft:.1e}")
sweep = []
for c2s in (0.0, 0.3, 0.6, 0.9, 1.2):
    Acon[1, 2] = c2s
    rs_, _ = scalar_from_loops("sphere", x7)
    rt_, _ = scalar_from_loops("torus", x7)
    sweep.append((c1v * c2s, rs_ - 6 / rhov**2, rt_))
Acon[1, 2] = c2v
okS = all(abs(a - (-(rhov**2 / 8) * f**2)) < 1e-5 for f, a, b in sweep) and all(abs(b) < 1e-6 for f, a, b in sweep)
check("I5 sweep: sphere's field term follows -(rho^2/8)(c1 c2)^2; torus stays flat", okS)
results["nonabelian_kk"] = {"R7_loops": R7_loops, "pred": pred7, "R_S3": 6 / rhov**2, "F3": F3, "sweep": sweep,
                            "planes_sphere": {f"{a}-{b}": v for (a, b), v in planes7.items()},
                            "planes_torus": {f"{a}-{b}": v for (a, b), v in planes7t.items()}}

print(f"\n{passed} passed, {failed} failed")
results["summary"] = {"passed": passed, "failed": failed}
json.dump(results, open("results.json", "w"), indent=1, default=lambda o: str(o))
