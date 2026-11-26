"""Polish (Nov 26): refine the numbers that search.py stopped at on a grid or at a loose tolerance.

Runs with numpy only (this machine has no scipy), so the pieces of quintessence.py and sm_rge.py that are needed
are copied here with scipy's solve_ivp / brentq / minimize_scalar replaced by a fixed-step RK4 and a small Brent.
Each section first re-computes the stored number at the stored point (to show the copy agrees), then polishes.

  1. exponential slope profile: search.py sampled chi2(lambda) at 17 points 0.1 apart and read the best slope and
     the 1-sigma edges off a straight-line interpolation.  Polish: Brent on the profile, root-finds for the edges.
  2. hilltop and cosine "free" fits: search.py picked the shape parameter from a 12-point grid.  Polish: Brent in
     the parameter (Om refitted at each point), full data; see the note in section 2 about the shooting bracket.
  3. Standard Model running: the scale where lambda turns negative and where beta_lambda is smallest were read off
     a 2000-point grid.  Polish: root-find lambda(t) = 0 and beta_lambda(t) = 0 on the RK4 solution.
"""
import os, math, time, json
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # keep CPU = wall
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DS = os.path.join(HERE, "..", "dark-sector", "data")
T0 = time.time()
OUT = []
def say(s=""):
    print(s, flush=True); OUT.append(s)

# ------------------------------------------------------------------ small numerical tools (numpy only)
def brent_min(f, a, b, tol=1e-10, maxit=200):
    """Brent's minimiser on [a, b] (golden section + parabolic steps). Returns (x, f(x), n_evals)."""
    g = 0.3819660112501051
    x = w = v = a + g * (b - a); fx = fw = fv = f(x); n = 1; d = e = 0.0
    for _ in range(maxit):
        m = 0.5 * (a + b); tol1 = tol * abs(x) + 1e-15; tol2 = 2 * tol1
        if abs(x - m) <= tol2 - 0.5 * (b - a):
            break
        if abs(e) > tol1:
            r = (x - w) * (fx - fv); q = (x - v) * (fx - fw); p = (x - v) * q - (x - w) * r; q = 2 * (q - r)
            if q > 0: p = -p
            q = abs(q); et = e; e = d
            if abs(p) >= abs(0.5 * q * et) or p <= q * (a - x) or p >= q * (b - x):
                e = (a - x) if x >= m else (b - x); d = g * e
            else:
                d = p / q; u = x + d
                if u - a < tol2 or b - u < tol2: d = tol1 if m >= x else -tol1
        else:
            e = (a - x) if x >= m else (b - x); d = g * e
        u = x + (d if abs(d) >= tol1 else (tol1 if d > 0 else -tol1))
        fu = f(u); n += 1
        if fu <= fx:
            if u >= x: a = x
            else: b = x
            v, w, x, fv, fw, fx = w, x, u, fw, fx, fu
        else:
            if u < x: a = u
            else: b = u
            if fu <= fw or w == x: v, w, fv, fw = w, u, fw, fu
            elif fu <= fv or v == x or v == w: v, fv = u, fu
    return x, fx, n

def root(f, a, b, xtol=1e-12, maxit=200):
    """Bracketed root (Illinois false position, falls back to bisection). Returns (x, n_evals)."""
    fa, fb = f(a), f(b); n = 2
    if fa * fb > 0: raise ValueError("no sign change")
    side = 0
    for _ in range(maxit):
        c = (a * fb - b * fa) / (fb - fa)
        if not (min(a, b) < c < max(a, b)): c = 0.5 * (a + b)
        fc = f(c); n += 1
        if fc * fb > 0:
            b, fb = c, fc
            if side == -1: fa *= 0.5
            side = -1
        else:
            a, fa = c, fc
            if side == 1: fb *= 0.5
            side = 1
        if abs(b - a) < xtol or fc == 0: break
    return c, n

# ------------------------------------------------------------------ data and chi2, copied from quintessence.py
P = os.path.join(DS, "pantheon")
tab = np.genfromtxt(os.path.join(P, "PantheonPlusSH0ES.dat"), names=True, dtype=None, encoding=None)
zhd, zhel, mb = tab["zHD"], tab["zHEL"], tab["m_b_corr"]
raw = np.loadtxt(os.path.join(P, "PantheonPlusSH0ES_STAT+SYS.cov"))
n_ = int(raw[0]); cov = raw[1:].reshape(n_, n_)
mask = zhd > 0.01
zhd, zhel, mb = zhd[mask], zhel[mask], mb[mask]
Cinv = np.linalg.inv(cov[np.ix_(mask, mask)])
ones = np.ones(len(mb)); S11 = ones @ Cinv @ ones
B = os.path.join(DS, "desi")
bz, bv, bq = [], [], []
for l in open(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_mean.txt")):
    if l.startswith("#") or not l.strip(): continue
    a_ = l.split(); bz.append(float(a_[0])); bv.append(float(a_[1])); bq.append(a_[2])
bz = np.array(bz); bv = np.array(bv)
bCinv = np.linalg.inv(np.loadtxt(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_cov.txt")))
ZG = np.concatenate([[0.0], np.logspace(-4, np.log10(3.0), 1500)])
Z_START = 10.0
NOUT = np.append(-np.log(1 + ZG[::-1])[:-1], 0.0)          # output N grid, increasing (same as quintessence.py)
NI = -math.log(1 + Z_START)
say(f"data loaded in {time.time() - T0:.1f} s")

def cumtrapz(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])

def chi2s(E):
    chi = cumtrapz(1.0 / E, ZG)
    dm = np.interp(zhd, ZG, chi); mu = 5 * np.log10((1 + zhel) * dm); d = mb - mu
    sn = float(d @ Cinv @ d - (ones @ Cinv @ d) ** 2 / S11)
    DM = np.interp(bz, ZG, chi); Eb = np.interp(bz, ZG, E); DH = 1 / Eb; DV = (bz * DM ** 2 * DH) ** (1 / 3)
    t = np.array([{"DM_over_rs": DM[i], "DH_over_rs": DH[i], "DV_over_rs": DV[i]}[q] for i, q in enumerate(bq)])
    hrd = (t @ bCinv @ bv) / (t @ bCinv @ t); r = bv - hrd * t
    return sn, float(r @ bCinv @ r)

def lam_of(shape, par, phi):
    if shape == "exp": return par
    if shape == "hilltop":
        den = 1 - phi * phi / (2 * par * par)
        if den < 1e-3: raise ValueError("rolled past the zero of V")
        return (phi / par ** 2) / den
    den = 1 + math.cos(phi / par)
    if den < 1e-4: raise ValueError("at the bottom of the circle")
    return math.sin(phi / par) / (par * den)

S15, S6 = math.sqrt(1.5), math.sqrt(6)
def rhs(x, y, phi, shape, par):
    lam = lam_of(shape, par, phi)
    if abs(x) > 1.2 or abs(y) > 1.2 or not math.isfinite(lam) or abs(lam) > 1e3: raise ValueError("field ran away")
    g = 1 + x * x - y * y
    return -3 * x + S15 * lam * y * y + 1.5 * x * g, -S15 * lam * x * y + 1.5 * y * g, S6 * x

H_STEP = 0.004
def integrate(yi, shape, par, phi0, dense=False):
    """RK4 from N = -ln 11 to 0 with step <= H_STEP; if dense, also the state at every output N (Hermite)."""
    x, y, p = 0.0, yi, phi0
    nk = np.concatenate([np.linspace(NI, NOUT[0], int(math.ceil((NOUT[0] - NI) / H_STEP)) + 1),
                         np.linspace(NOUT[0], 0.0, int(math.ceil(-NOUT[0] / H_STEP)) + 1)[1:]])
    S = np.empty((len(nk), 3)); D = np.empty((len(nk), 3))
    S[0] = (x, y, p); D[0] = rhs(x, y, p, shape, par)
    for i in range(1, len(nk)):
        h = nk[i] - nk[i - 1]
        k1 = D[i - 1]
        k2 = rhs(x + .5 * h * k1[0], y + .5 * h * k1[1], p + .5 * h * k1[2], shape, par)
        k3 = rhs(x + .5 * h * k2[0], y + .5 * h * k2[1], p + .5 * h * k2[2], shape, par)
        k4 = rhs(x + h * k3[0], y + h * k3[1], p + h * k3[2], shape, par)
        x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        y += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        p += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        S[i] = (x, y, p); D[i] = rhs(x, y, p, shape, par)
    if not dense:
        return S[-1]
    j = np.clip(np.searchsorted(nk, NOUT, side="right") - 1, 0, len(nk) - 2)
    h = nk[j + 1] - nk[j]; s = ((NOUT - nk[j]) / h)[:, None]
    h00 = 2 * s ** 3 - 3 * s ** 2 + 1; h10 = s ** 3 - 2 * s ** 2 + s; h01 = -2 * s ** 3 + 3 * s ** 2; h11 = s ** 3 - s ** 2
    return h00 * S[j] + h10 * h[:, None] * D[j] + h01 * S[j + 1] + h11 * h[:, None] * D[j + 1]

def history(Om, shape="exp", par=0.0, phi0=0.0):
    def miss(lyi):
        s = integrate(math.exp(lyi), shape, par, phi0)
        return s[0] ** 2 + s[1] ** 2 - (1 - Om)
    y_l = math.sqrt((1 - Om) / (Om * (1 + Z_START) ** 3 + 1 - Om))
    lyi, _ = root(miss, math.log(y_l) - 3, math.log(y_l) + 3, xtol=1e-12)
    st = integrate(math.exp(lyi), shape, par, phi0, dense=True)
    x, y = st[:, 0], st[:, 1]
    g = 1 + x * x - y * y
    lnH = cumtrapz(-1.5 * g, NOUT); lnH -= lnH[-1]
    E = np.exp(lnH)[::-1]; w = ((x * x - y * y) / (x * x + y * y))[::-1]
    return E, w

NCALL = [0]
def total(Om, shape, par, phi0):
    NCALL[0] += 1
    try:
        sn, bao = chi2s(history(Om, shape, par, phi0)[0])
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return 1e9
    return sn + bao

def fit_Om(shape, par, phi0, tol=1e-8):
    """Om refitted at fixed shape parameter (search.py used xatol = 1e-4 on Om in [0.2, 0.45]; every best Om it
    found lies in 0.297..0.332, so the polish searches 0.28..0.34 and checks it did not stop at an end)."""
    Om, c, _ = brent_min(lambda O: total(O, shape, par, phi0), 0.28, 0.34, tol=tol)
    if c < 1e8 and not 0.2805 < Om < 0.3395: raise RuntimeError(f"Om fit at the edge of its bracket: {Om}")
    return Om, c

# ------------------------------------------------------------------ 1. exponential slope profile
R = json.load(open(os.path.join(HERE, "results.json")))
sp = R["slope_profile"]
lams = np.array(sp["lam"]); ch = np.array(sp["chi2"])
say("")
say("1. Exponential slope from all data (profile chi2, Om refitted at each slope)")
t1 = time.time()
c0 = fit_Om("exp", 0.0, 0.0)[1]; c8 = fit_Om("exp", 0.8, 0.0)[1]
say(f"   copy agrees with stored profile: chi2(0) {c0:.5f} vs {ch[0]:.5f}, chi2(0.8) {c8:.5f} vs {ch[8]:.5f}")
prof = lambda l: fit_Om("exp", l, 0.0)[1]
lb, cb, nb = brent_min(prof, 0.7, 0.9, tol=1e-5)
lo, nlo = root(lambda l: prof(l) - (cb + 1), 0.5, 0.7, xtol=1e-5)
hi, nhi = root(lambda l: prof(l) - (cb + 1), 0.9, 1.1, xtol=1e-5)
d_l = c0 - cb
t_1 = time.time() - t1
say(f"   best slope      {sp['best']:.4f} (grid + straight-line interpolation) -> {lb:.4f} (polished)")
say(f"   1-sigma range   {sp['one_sigma'][0]:.3f} to {sp['one_sigma'][1]:.3f} -> {lo:.3f} to {hi:.3f}")
say(f"   chi2 at best    {ch.min():.4f} -> {cb:.4f};  Lambda (slope 0) worse by {ch[0] - ch.min():.3f} -> {d_l:.3f}")
geo = sorted({round(math.sqrt(p / q), 6) for p in range(1, 9) for q in range(1, 9)})
hits_old = [g for g in geo if sp["one_sigma"][0] <= g <= sp["one_sigma"][1]]
hits = [g for g in geo if lo <= g <= hi]
say(f"   simple slopes sqrt(p/q), p,q <= 8, inside 1 sigma: {len(hits_old)} of {len(geo)} -> {len(hits)} of {len(geo)}")
say(f"   sqrt(2/3) = {math.sqrt(2/3):.4f} is {abs(lb - math.sqrt(2/3)):.4f} from the polished best slope")
say(f"   cost {t_1:.1f} s ({NCALL[0]} cosmology solves)")

# ------------------------------------------------------------------ 2. hilltop and cosine free fits
# A first polish (Brent in the shape parameter with quintessence.py's history() unchanged) moved the cosine period
# 1.036 -> 0.9965 (delta chi2 3.33 -> 3.80) and the hilltop width 1.055 -> 0.835 (0.14 -> 0.42), but both stopped at
# a wall.  The wall is not physics: history() brackets the starting field between 1/20 and 20 times the
# cosmological-constant value, and when the top end of that bracket rolls the field to the bottom of the circle (or
# past the zero of the hilltop) the probe raises and the whole history is discarded, even though the solution
# itself is fine.  The grid in search.py met the same wall (every grid point below it scored 1e9).  So the polish
# below shrinks the top of the bracket until the probe runs (the solution found is the same where both work), then
# runs Brent in the shape parameter.  This is a fix to the search, not just more digits.
say("")
say("2. Hump-shaped potentials, shape parameter fitted (search.py: 12-point grid, Om refitted at each)")
stored = {}
for rd in R["rounds"]:
    for g in rd:
        if g["id"] in ("quint:cosine:free", "quint:hilltop:free"):
            stored[g["id"].split(":")[1]] = g["extra"]["full"]
PHI0 = {"hilltop": 0.1, "cosine": 1.0}
BRACKET = {"cosine": (0.82, 0.98), "hilltop": (0.46, 0.60)}   # from a coarse look with the fixed bracket

def history_fixed(Om, shape, par, phi0):
    def miss(lyi):
        s = integrate(math.exp(lyi), shape, par, phi0)
        return s[0] ** 2 + s[1] ** 2 - (1 - Om)
    yl = math.log(math.sqrt((1 - Om) / (Om * (1 + Z_START) ** 3 + 1 - Om)))
    w = 3.0
    while True:
        try:
            miss(yl + w); break
        except ValueError:
            w *= 0.7
            if w < 0.05: raise
    lyi, _ = root(miss, yl - 3, yl + w, xtol=1e-12)
    st = integrate(math.exp(lyi), shape, par, phi0, dense=True)
    x, y = st[:, 0], st[:, 1]; g = 1 + x * x - y * y
    lnH = cumtrapz(-1.5 * g, NOUT); lnH -= lnH[-1]
    return np.exp(lnH)[::-1], ((x * x - y * y) / (x * x + y * y))[::-1]

def total_fixed(Om, shape, par, phi0):
    NCALL[0] += 1
    try:
        return sum(chi2s(history_fixed(Om, shape, par, phi0)[0]))
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return 1e9

lcdm = c0
for shape in ("cosine", "hilltop"):
    t2 = time.time(); NCALL[0] = 0
    st = stored[shape]
    c_chk = total(st["Om"], shape, st["par"], PHI0[shape])
    c_chk2 = total_fixed(st["Om"], shape, st["par"], PHI0[shape])
    fitf = lambda p: brent_min(lambda O: total_fixed(O, shape, p, PHI0[shape]), 0.28, 0.34, tol=1e-6)[:2]
    pb, cbest, _ = brent_min(lambda p: fitf(p)[1], *BRACKET[shape], tol=1e-4)
    Omb = fitf(pb)[0]
    c_old = total(Omb, shape, pb, PHI0[shape])
    E, w = history_fixed(Omb, shape, pb, PHI0[shape])
    say(f"   {shape:8s} copy agrees: chi2 at stored point {c_chk:.4f} (fixed bracket {c_chk2:.4f}) vs {st['chi2']:.4f}")
    say(f"   {shape:8s} parameter {st['par']:.4f} -> {pb:.4f};  Om {st['Om']:.4f} -> {Omb:.4f};  w0 {st['w0']:.3f} -> {w[0]:.3f};  "
        f"delta chi2 vs Lambda {st['dchi2']:.2f} -> {lcdm - cbest:.2f} (polished, fixed bracket)")
    say(f"   {shape:8s} the unchanged history() at that point returns {'failure (1e9)' if c_old >= 1e8 else f'{c_old:.4f}'};"
        f"  cost {time.time() - t2:.1f} s ({NCALL[0]} solves)")

# ------------------------------------------------------------------ 3. Standard Model running
say("")
say("3. Standard Model running with measured inputs (search.py: read off a 2000-point grid)")
K = 1.0 / (16 * math.pi ** 2); MPL_FULL = 1.22089e19
def inputs_at_mt(Mh, Mt, as_mz, MW=80.3692):
    a3 = (as_mz - 0.1184) / 0.0007
    lam = 0.12604 + 0.00206 * (Mh - 125.15) - 0.00004 * (Mt - 173.34)
    yt = 0.93690 + 0.00556 * (Mt - 173.34) - 0.00042 * a3
    g3 = 1.1666 + 0.00314 * a3
    g2 = 0.64779 + 0.00004 * (Mt - 173.34) + 0.00011 * (MW - 80.384) / 0.014
    gY = 0.35830 + 0.00011 * (Mt - 173.34) - 0.00020 * (MW - 80.384) / 0.014
    return np.array([gY, g2, g3, yt, lam])
def betas(c):
    gY, g2, g3, yt, lam = c
    Y2, W2, G2, T2 = gY * gY, g2 * g2, g3 * g3, yt * yt
    bY = 41 / 6 * gY ** 3; bW = -19 / 6 * g2 ** 3; bG = -7 * g3 ** 3
    bt = yt * (4.5 * T2 - 8 * G2 - 2.25 * W2 - 17 / 12 * Y2)
    bl = (24 * lam ** 2 - 6 * T2 ** 2 + 3 / 8 * (2 * W2 ** 2 + (W2 + Y2) ** 2) + lam * (12 * T2 - 9 * W2 - 3 * Y2))
    out = K * np.array([bY, bW, bG, bt, bl])
    bY2 = gY ** 3 * (199 / 18 * Y2 + 9 / 2 * W2 + 44 / 3 * G2 - 17 / 6 * T2)
    bW2 = g2 ** 3 * (3 / 2 * Y2 + 35 / 6 * W2 + 12 * G2 - 3 / 2 * T2)
    bG2 = g3 ** 3 * (11 / 6 * Y2 + 9 / 2 * W2 - 26 * G2 - 2 * T2)
    bt2 = yt * (-12 * T2 ** 2 + 6 * lam ** 2 - 12 * lam * T2 + 36 * G2 * T2 - 108 * G2 ** 2)
    bl2 = (-312 * lam ** 3 - 144 * lam ** 2 * T2 - 3 * lam * T2 ** 2 + 30 * T2 ** 3 + 80 * lam * G2 * T2 - 32 * G2 * T2 ** 2)
    return out + K * K * np.array([bY2, bW2, bG2, bt2, bl2])
def rk4_to(c, t0, t1, h=0.01):
    n = max(1, int(math.ceil(abs(t1 - t0) / h))); h = (t1 - t0) / n
    for _ in range(n):
        k1 = betas(c); k2 = betas(c + .5 * h * k1); k3 = betas(c + .5 * h * k2); k4 = betas(c + h * k3)
        c = c + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return c
t3 = time.time()
MH, MT, AS = 125.20, 172.57, 0.1180
fd = R["flat_dial"]
c0_ = inputs_at_mt(MH, MT, AS)
t_end = math.log(MPL_FULL / MT)
lam_end = rk4_to(c0_, 0.0, t_end)[4]
say(f"   copy agrees: lambda(M_Pl) {lam_end:+.6f} vs stored {fd['lambda_MPl_measured_inputs']:+.6f}")
# lambda = 0: bracket one grid step either side of the stored crossing
tg = math.log(fd["lambda_negative_from"] / MT); dt = t_end / 1999
ca = rk4_to(c0_, 0.0, tg - 2 * dt)
f_l = lambda t: rk4_to(ca, tg - 2 * dt, t, h=0.002)[4]
tz, _ = root(f_l, tg - 2 * dt, tg + 2 * dt, xtol=1e-12)
mu_z = MT * math.exp(tz)
tb = math.log(fd["beta_min_at"] / MT)
cb_ = rk4_to(c0_, 0.0, tb - 2 * dt)
f_b = lambda t: betas(rk4_to(cb_, tb - 2 * dt, t, h=0.002))[4]
tbz, _ = root(f_b, tb - 2 * dt, tb + 2 * dt, xtol=1e-12)
mu_b = MT * math.exp(tbz); lam_b = rk4_to(cb_, tb - 2 * dt, tbz, h=0.002)[4]
say(f"   lambda turns negative at mu = {fd['lambda_negative_from']:.3e} GeV (first grid point past zero) -> {mu_z:.4e} GeV (polished)")
say(f"   beta_lambda smallest at   mu = {fd['beta_min_at']:.3e} GeV (grid argmin |beta|) -> {mu_b:.4e} GeV where beta = 0 (polished);"
    f" lambda there {fd['lambda_at_beta_min']:+.4f} -> {lam_b:+.4f}")
say(f"   cost {time.time() - t3:.1f} s")
say("")
say(f"total wall time {time.time() - T0:.1f} s, CPU {time.process_time():.1f} s")
open(os.path.join(HERE, "polish_output.txt"), "w").write("\n".join(OUT) + "\n")
