"""Polish step for the fitted numbers in expansion_fits.py and sparc_fits.py.

What those scripts searched, and how:
  expansion_fits.py  Nelder-Mead (xatol 1e-5, fatol 1e-4) for every model with free numbers, three data
                     combinations each (supernovae, BAO, both).
  sparc_fits.py      one-number gravity laws: 61-point grid, then bounded Brent (default xatol 1e-5);
                     the constant-background law: Nelder-Mead (xatol 1e-5, fatol 1e-4);
                     mass-to-light free per galaxy: bounded Brent per galaxy (xatol 1e-5), and for the
                     exp-map law a coarse grid of ten scales; NFW halos: least_squares per galaxy, 9 starts.
The polish (numpy only; scipy is not installed here, so the integrals and the bounded Brent search are
small numpy copies, the latter an exact port of scipy's routine):
  * Newton steps (central-difference gradient and Hessian) from every stored answer, until the step is
    below 1e-10; a coordinate sitting on its bound is held there;
  * per-galaxy mass-to-light: the stored-tolerance search is reproduced, then refined to xatol 1e-12
    (an optimum on the 1.5 bound is set exactly on it);
  * the exp-map free-M/L scale: a finer scan (0.25 in log10 s) over the whole range, with a local
    refinement around its best point, to see whether the coarse grid missed a better scale;
  * NFW: the per-galaxy answers were not stored, so each galaxy is re-solved from the same 9 starts with
    Levenberg-Marquardt to 1e-14 and the sums compared.
All of these are fits to measured data with real noise: digits beyond those on the page carry no
information, and the polish mostly confirms that the optimizers had converged.
Run: python3 polish.py > polish_output.txt
"""
import json, os, math, time
import numpy as np
np.seterr(all="ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.perf_counter()

def fminbound(func, x1, x2, xatol=1e-5, maxfun=500):
    """Port of scipy.optimize._minimize_scalar_bounded (bounded Brent)."""
    sqrt_eps = np.sqrt(2.2e-16); golden_mean = 0.5 * (3.0 - np.sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a); nfc, xf = fulc, fulc
    rat = e = 0.0; x = xf; fx = func(x); num = 1
    ffulc = fnfc = fx; xm = 0.5 * (a + b)
    tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
    while np.abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = 1
        if np.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc); q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r; q = 2.0 * (q - r)
            if q > 0.0: p = -p
            q = np.abs(q); r = e; e = rat
            if (np.abs(p) < np.abs(0.5 * q * r)) and (p > q * (a - xf)) and (p < q * (b - xf)):
                rat = (p + 0.0) / q; x = xf + rat
                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = np.sign(xm - xf) + ((xm - xf) == 0); rat = tol1 * si
            else:
                golden = 1
        if golden:
            e = (a - xf) if xf >= xm else (b - xf)
            rat = golden_mean * e
        si = np.sign(rat) + (rat == 0)
        x = xf + si * np.maximum(np.abs(rat), tol1)
        fu = func(x); num += 1
        if fu <= fx:
            if x >= xf: a = xf
            else: b = xf
            fulc, ffulc = nfc, fnfc; nfc, fnfc = xf, fx; xf, fx = x, fu
        else:
            if x < xf: a = x
            else: b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc; nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu
        xm = 0.5 * (a + b); tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
        if num >= maxfun: break
    return xf

def lm(fun, x0, lo=None, hi=None, xtol=1e-14, ftol=1e-15, maxit=500):
    """Small Levenberg-Marquardt (numpy only; scipy is not installed here). Central-difference Jacobian,
    box bounds by clipping. Returns (x, residuals, iterations)."""
    x = np.array(x0, float); n = len(x)
    lo = np.full(n, -np.inf) if lo is None else np.array(lo, float)
    hi = np.full(n, np.inf) if hi is None else np.array(hi, float)
    x = np.clip(x, lo, hi); r = fun(x); cost = r @ r; lam = 1e-6
    for it in range(1, maxit + 1):
        J = np.empty((len(r), n))
        for i in range(n):
            h = 1e-6 * max(1.0, abs(x[i])); e = np.zeros(n); e[i] = h
            J[:, i] = (fun(x + e) - fun(x - e)) / (2 * h)
        g = J.T @ r; A = J.T @ J
        free = ~(((x <= lo) & (g > 0)) | ((x >= hi) & (g < 0)))   # coordinates not pinned at a bound
        improved = False
        for _ in range(30):
            dx = np.zeros(n)
            Af = A[np.ix_(free, free)]
            try:
                dx[free] = np.linalg.solve(Af + lam * np.diag(np.diag(Af) + 1e-300), -g[free])
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = np.clip(x + dx, lo, hi); rn = fun(xn); cn = rn @ rn
            if np.isfinite(cn) and cn <= cost:
                improved = True; break
            lam *= 10
        if not improved:
            break
        step = np.max(np.abs(xn - x) / (np.abs(x) + 1e-12)); drop = (cost - cn) / max(cost, 1e-300)
        x, r, cost = xn, rn, cn; lam = max(lam / 10, 1e-12)
        if step < xtol or drop < ftol:
            break
    return x, r, it

def newton_polish(f, x0, bounds, h=1e-5, maxit=30, tol=1e-10):
    """Newton on a smooth objective with central-difference derivatives; bound-pinned coordinates held."""
    x = np.array(x0, float); n = len(x); lo = np.array([b[0] for b in bounds]); hi = np.array([b[1] for b in bounds])
    fx = f(x)
    for it in range(1, maxit + 1):
        hs = h*np.maximum(1, np.abs(x)); g = np.zeros(n); H = np.zeros((n, n))
        for i in range(n):
            ei = np.zeros(n); ei[i] = hs[i]
            fp, fm = f(x + ei), f(x - ei); g[i] = (fp - fm)/(2*hs[i]); H[i, i] = (fp - 2*fx + fm)/hs[i]**2
            for j in range(i):
                ej = np.zeros(n); ej[j] = hs[j]
                H[i, j] = H[j, i] = (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej))/(4*hs[i]*hs[j])
        free = ~(((x - lo) < 1e-12) & (g > 0) | ((hi - x) < 1e-12) & (g < 0))
        dx = np.zeros(n)
        try: dx[free] = -np.linalg.solve(H[np.ix_(free, free)], g[free])
        except np.linalg.LinAlgError: break
        lam = 1.0
        while lam > 1e-6:
            xn = np.clip(x + lam*dx, lo, hi); fn = f(xn)
            if fn <= fx + 1e-12*abs(fx): break
            lam /= 2
        else:
            break
        step = np.max(np.abs(xn - x)); x, fx = xn, fn
        if step < tol: break
    return x, fx, it

# =====================================================================  expansion_fits.py
C = 299792.458
P = os.path.join(HERE, "data", "pantheon")
tab = np.genfromtxt(os.path.join(P, "PantheonPlusSH0ES.dat"), names=True, dtype=None, encoding=None)
zhd, zhel, mb = tab["zHD"], tab["zHEL"], tab["m_b_corr"]
raw = np.loadtxt(os.path.join(P, "PantheonPlusSH0ES_STAT+SYS.cov"))
n = int(raw[0]); cov = raw[1:].reshape(n, n)
mask = zhd > 0.01
zhd, zhel, mb = zhd[mask], zhel[mask], mb[mask]
cov = cov[np.ix_(mask, mask)]; Cinv = np.linalg.inv(cov); ones = np.ones(len(mb)); S11 = ones @ Cinv @ ones
ZG = np.concatenate([[0.0], np.logspace(-4, np.log10(3.0), 4000)])
def cumtrapz0(y, x): return np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*np.diff(x))])
def chi_of(model, p, z):
    kind, f = model
    if kind == "E":
        chi = cumtrapz0(1.0/f(ZG, p), ZG); v = np.interp(z, ZG, chi); return v, v
    return f(z, p)
def sn_chi2(model, p):
    chi, dm = chi_of(model, p, zhd)
    if np.any(~np.isfinite(dm)) or np.any(dm <= 0): return 1e12
    d = mb - 5*np.log10((1 + zhel)*dm)
    return float(d @ Cinv @ d - (ones @ Cinv @ d)**2/S11)
B = os.path.join(HERE, "data", "desi"); bz, bv, bq = [], [], []
for l in open(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_mean.txt")):
    if l.startswith("#") or not l.strip(): continue
    a = l.split(); bz.append(float(a[0])); bv.append(float(a[1])); bq.append(a[2])
bz = np.array(bz); bv = np.array(bv); bCinv = np.linalg.inv(np.loadtxt(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_cov.txt")))
def bao_chi2(model, p):
    kind, f = model; chi, dm = chi_of(model, p, bz); E = f(bz, p)
    DM = dm; DH = 1/E; DV = (bz*DM**2*DH)**(1/3)
    t = np.array([{"DM_over_rs": DM[i], "DH_over_rs": DH[i], "DV_over_rs": DV[i]}[q] for i, q in enumerate(bq)])
    if np.any(~np.isfinite(t)): return 1e12
    hrd = (t @ bCinv @ bv)/(t @ bCinv @ t); d = bv - hrd*t
    return float(d @ bCinv @ d)
def E_lcdm(z, p):  return np.sqrt(p[0]*(1 + z)**3 + 1 - p[0])
def E_wcdm(z, p):  return np.sqrt(p[0]*(1 + z)**3 + (1 - p[0])*(1 + z)**(3*(1 + p[1])))
def E_power(z, p): return (1 + z)**(1 - p[0])
def E_cube(z, p):  return np.sqrt(np.maximum(1 + p[0]*3*z + p[1]*3*z**2 + p[2]*z**3, 1e-6))
def E_scaled_cube(z, p): return (1 + p[0]*z)**1.5
MODELS = {
    "Flat LCDM (cosmological constant)": (("E", E_lcdm), [(0.01, 1.0)]),
    "Flat wCDM (dark energy with w free)": (("E", E_wcdm), [(0.01, 1.0), (-3, 0.3)]),
    "Fused power law, chi = ((1+z)^b - 1)/b": (("E", E_power), [(-1.0, 1.5)]),
    "Cube of scaled redshift, E^2 = (1 + l z)^3": (("E", E_scaled_cube), [(0.0, 3.0)]),
    "Circle-operation cube, free weight per layer": (("E", E_cube), [(-2, 3), (-2, 3), (-2, 3)]),
}
R = json.load(open(os.path.join(HERE, "expansion_results.json")))
t_load = time.perf_counter() - T0
print(f"=== expansion_fits.py (Pantheon+ {mask.sum()} SNe, DESI {len(bv)} BAO numbers; data loaded in {t_load:.1f} s)")
t = time.perf_counter(); polished = {}
for label, (model, bounds) in MODELS.items():
    for which in ("sn", "bao", "sn+bao"):
        obj = lambda q: (sn_chi2(model, q) if "sn" in which else 0.0) + (bao_chi2(model, q) if "bao" in which else 0.0)
        p0 = R["fits"][label][which]["params"]; c0 = R["fits"][label][which]["chi2"]
        c0_here = obj(p0)
        p1, c1, it = newton_polish(obj, p0, bounds)
        polished[(label, which)] = p1
        pin = [i for i, (x, (lo, hi)) in enumerate(zip(p1, bounds)) if min(x - lo, hi - x) < 1e-9]
        print(f"  {label[:44]:44s} [{which:6s}] params {' '.join(f'{x:+.6f}' for x in p0)} -> {' '.join(f'{x:+.6f}' for x in p1)}"
              f"   chi2 {c0:.4f} -> {c1:.4f} (d = {c1-c0:+.1e}; recomputed here {c0_here:.4f}){'  on bound' if pin else ''}")
t_exp = time.perf_counter() - t
print(f"  expansion polish time {t_exp:.1f} s")
p0 = R["fits"]["Flat LCDM (cosmological constant)"]["sn"]["params"][0]; p1 = polished[("Flat LCDM (cosmological constant)", "sn")][0]
print(f"  page: Omega_m from supernovae {p0:.3f} -> {p1:.3f}" + ("  (no change)" if f"{p0:.3f}" == f"{p1:.3f}" else "  CHANGED"))
b0 = R["fits"]["Fused power law, chi = ((1+z)^b - 1)/b"]["sn"]["params"][0]; b1 = polished[("Fused power law, chi = ((1+z)^b - 1)/b", "sn")][0]
print(f"  page: power-law beta from supernovae {b0:.2f} -> {b1:.2f}" + ("  (no change)" if f"{b0:.2f}" == f"{b1:.2f}" else "  CHANGED"))
U = json.load(open(os.path.join(HERE, "uncertainties.json")))
key = "Circle-operation cube, free weight per layer"; Cc = np.array(U[f"{key} [sn+bao]"]["cov"]); J = np.array([[1, -1, 0], [0, 1, -1]])
for tag, p in (("before", np.array(R["fits"][key]["sn+bao"]["params"])), ("after", polished[(key, "sn+bao")])):
    d = J @ p; print(f"  page: layer weights {tag}: {' '.join(f'{x:.2f}' for x in p)}; equal-weights test chi2 = {d @ np.linalg.inv(J @ Cc @ J.T) @ d:.2f} (stored covariance)")

# =====================================================================  sparc_fits.py
t = time.perf_counter()
D = os.path.join(HERE, "data", "sparc"); KMS2_PER_KPC = 1e6/3.0857e19; UD, UB = 0.5, 0.7
lines = open(os.path.join(D, "SPARC_Lelli2016c.mrt")).read().splitlines()
start = max(i for i, l in enumerate(lines) if l.startswith("-----")) + 1
meta = {}
for l in lines[start:]:
    if not l.strip(): continue
    f = l.split(); meta[f[0]] = dict(D=float(f[2]), inc=float(f[5]), Q=int(f[17]))
gals = []
for name, m in meta.items():
    if m["Q"] >= 3 or m["inc"] < 30: continue
    a = np.atleast_2d(np.loadtxt(os.path.join(D, "rotmod", f"{name}_rotmod.dat"), comments="#"))
    Rr, V, eV, Vg, Vd, Vb = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
    ok = (Rr > 0) & (V > 0) & (eV > 0)
    gals.append(dict(name=name, R=Rr[ok], V=V[ok], eV=eV[ok], Vg=Vg[ok], Vd=Vd[ok], Vb=Vb[ok]))
def pieces(g, ud=UD, ub=UB):
    R_ = g["R"]
    return (g["Vg"]*np.abs(g["Vg"])/R_*KMS2_PER_KPC, ud*g["Vd"]**2/R_*KMS2_PER_KPC, ub*g["Vb"]**2/R_*KMS2_PER_KPC)
def newton_g(gg, gd, gb, p): return gg + gd + gb
def madd_exp(gg, gd, gb, p):
    s = 10**p[0]; x = (gg + gd + gb)/s; return s*np.expm1(np.clip(x, -50, 50))
def madd_pieces(gg, gd, gb, p):
    s = 10**p[0]; return s*((1 + gg/s)*(1 + gd/s)*(1 + gb/s) - 1)
def madd_background(gg, gd, gb, p):
    k, s = p[0], 10**p[1]; return (1 + k)*(gg + gd + gb) + k*s
def mond_simple(gg, gd, gb, p):
    a0 = 10**p[0]; gN = np.maximum(gg + gd + gb, 1e-16); return gN/2 + np.sqrt(gN**2/4 + gN*a0)
def rar(gg, gd, gb, p):
    a0 = 10**p[0]; gN = np.maximum(gg + gd + gb, 1e-16); return gN/(-np.expm1(-np.sqrt(gN/a0)))
GLOBAL = {
    "Newton (visible matter only)": (newton_g, []),
    "The circle operation of many small pieces (exp map)": (madd_exp, [(-14, -2)]),
    "The circle operation of gas, disk, bulge": (madd_pieces, [(-14, -2)]),
    "The circle operation with a constant background": (madd_background, [(-1, 5), (-14, -8)]),
    "MOND, simple mu (= minus circle inverse)": (mond_simple, [(-12, -8)]),
    "MOND / RAR (McGaugh 2016)": (rar, [(-12, -8)]),
}
PC = [pieces(g) for g in gals]
def chi2_gal(gi, f, p, ud=UD, ub=UB):
    g = gals[gi]
    if ud == UD and ub == UB: gg, gd, gb = PC[gi]
    else: gg, gd, gb = pieces(g, ud, ub)
    Vm = np.sqrt(np.maximum(f(gg, gd, gb, p), 0)*g["R"]/KMS2_PER_KPC)
    return float(np.sum(((g["V"] - Vm)/g["eV"])**2))
def total_chi2(f, p): return sum(chi2_gal(i, f, p) for i in range(len(gals)))
S = json.load(open(os.path.join(HERE, "sparc_results.json")))
print(f"=== sparc_fits.py ({len(gals)} galaxies; stored {S['n_galaxies']})")
for label, (f, bounds) in GLOBAL.items():
    if not bounds: continue
    st = S["fixed_ml"][label]; p0 = st["params"]
    tt = time.perf_counter()
    p1, c1, it = newton_polish(lambda q: total_chi2(f, q), p0, bounds, h=1e-4)
    sc = lambda v: f"{10**v:.3g}"
    print(f"  {label[:44]:44s} params {' '.join(f'{x:+.7f}' for x in p0)} -> {' '.join(f'{x:+.7f}' for x in p1)}"
          f"  (scale {sc(p0[-1])} -> {sc(p1[-1])})  chi2 {st['chi2']:.3f} -> {c1:.3f}  [{time.perf_counter()-tt:.1f} s]")

# free mass-to-light per galaxy
def fit_free_ml(gi, f, p, xatol):
    u = fminbound(lambda u: chi2_gal(gi, f, p, u, u*UB/UD), 0.05, 1.5, xatol=xatol)
    return chi2_gal(gi, f, p, u, u*UB/UD), u
def free_ml_polished(gi, f, p):
    c0, u0 = fit_free_ml(gi, f, p, 1e-5)
    lo, hi = max(0.05, u0 - 1e-3), min(1.5, u0 + 1e-3)
    u1 = fminbound(lambda u: chi2_gal(gi, f, p, u, u*UB/UD), lo, hi, xatol=1e-12)
    cands = [(chi2_gal(gi, f, p, u1, u1*UB/UD), u1)] + [(chi2_gal(gi, f, p, b, b*UB/UD), b) for b in (0.05, 1.5) if abs(b - u0) < 1e-3]
    return (c0, u0), min(cands)
for label in ["Newton (visible matter only)", "MOND / RAR (McGaugh 2016)"]:
    tt = time.perf_counter(); f, _ = GLOBAL[label]; st = S["free_ml"][label]; p = st["params"]
    rows = [free_ml_polished(i, f, p) for i in range(len(gals))]
    cb = sum(r[0][0] for r in rows); ca = sum(r[1][0] for r in rows)
    mb_, ma = np.median([r[0][1] for r in rows]), np.median([r[1][1] for r in rows])
    print(f"  free M/L {label[:34]:34s} chi2 stored {st['chi2']:.3f}, reproduced {cb:.3f} -> polished {ca:.3f};"
          f" median M/L {mb_:.6f} -> {ma:.6f}  [{time.perf_counter()-tt:.1f} s]")
# exp-map scale with free M/L: finer scan than the coarse ten-point grid
tt = time.perf_counter(); f, _ = GLOBAL["The circle operation of many small pieces (exp map)"]
tot_free = lambda q: sum(fit_free_ml(i, f, [q], 1e-5)[0] for i in range(len(gals)))
grid = np.arange(-14, -2 + 1e-9, 0.25); vals = np.array([tot_free(x) for x in grid])
i = int(np.argmin(vals)); st = S["free_ml"]["The circle operation of many small pieces (exp map)"]
xb = fminbound(tot_free, grid[max(i-1, 0)], grid[min(i+1, len(grid)-1)], xatol=1e-3)
interior = [(vals[j], grid[j]) for j in range(1, len(grid)-1) if vals[j] < vals[j-1] and vals[j] < vals[j+1]]
print(f"  free M/L exp-map scale: coarse grid best log10 s = {st['params'][0]:+.2f} (chi2 {st['chi2']:.3f}); fine scan best {grid[i]:+.2f}"
      f" (chi2 {vals[i]:.3f}), refined {xb:+.3f} (chi2 {tot_free(xb):.3f}); interior local minima of the scan: "
      + (", ".join(f"{x:+.2f} ({v:.1f})" for v, x in interior) if interior else "none")
      + f"; Newton with free M/L {S['free_ml']['Newton (visible matter only)']['chi2']:.3f}  [{time.perf_counter()-tt:.1f} s]")
cpol = sum(free_ml_polished(i, f, [xb])[1][0] for i in range(len(gals)))
print(f"  free M/L exp-map at log10 s = {xb:+.3f}, per-galaxy M/L polished to 1e-12: chi2 {st['chi2']:.3f} -> {cpol:.3f}")

# NFW halos: re-solve each galaxy from the same 9 starts
tt = time.perf_counter(); H0_KPC = 0.073
def v_nfw2(Rr, V200, c):
    R200 = V200/(10*H0_KPC); x = Rr/R200; m = lambda y: np.log1p(y) - y/(1 + y)
    return V200**2*m(c*x)/(x*m(c))
def fit_nfw(gi, free_ml=False):
    g = gals[gi]; gg, gd, gb = pieces(g, 1.0, 1.0); Rr = g["R"]
    bg, bd, bb = gg*Rr/KMS2_PER_KPC, gd*Rr/KMS2_PER_KPC, gb*Rr/KMS2_PER_KPC
    def res(q):
        V200, c = 10**q[0], 10**q[1]; u_d, u_b = (q[2], q[2]*UB/UD) if free_ml else (UD, UB)
        V2 = bg + u_d*bd + u_b*bb + v_nfw2(Rr, V200, c)
        return (g["V"] - np.sqrt(np.maximum(V2, 0)))/g["eV"]
    best = None
    for lv in (1.5, 2.0, 2.5):
        for lc in (0.3, 0.8, 1.3):
            q0 = [lv, lc] + ([UD] if free_ml else []); lo = [1.0, 0.0] + ([0.05] if free_ml else []); hi = [2.8, 2.0] + ([1.5] if free_ml else [])
            q, r, _ = lm(res, q0, lo, hi); c2 = float(r @ r)
            if best is None or c2 < best: best = c2
    return best
for tag, key, fm in (("fixed M/L", "fixed_ml", False), ("free M/L", "free_ml", True)):
    c2 = sum(fit_nfw(i, fm) for i in range(len(gals)))
    print(f"  NFW {tag}: chi2 stored {S[key]['NFW dark-matter halo (2 per galaxy)' if not fm else 'NFW dark-matter halo']['chi2']:.3f} -> re-solved {c2:.3f}")
print(f"  NFW time {time.perf_counter()-tt:.1f} s")
print(f"  sparc polish time {time.perf_counter()-t:.1f} s")
print(f"total wall time {time.perf_counter()-T0:.1f} s")
