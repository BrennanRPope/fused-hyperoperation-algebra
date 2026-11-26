"""The dial as a slowly rolling field today: expansion history against Pantheon+ supernovae and DESI DR2 BAO.

Data loaders follow ../dark-sector/expansion_fits.py (same files, same marginalisation: supernovae fix only
the shape of the distance curve, M and H0 are marginalised; BAO keeps one nuisance number c/(H0 r_d)).

Field equations (units reduced Planck mass = 1, N = ln a, matter + field, radiation ignored below z = 10):
  x = phi'/(sqrt6 H),  y = sqrt(V)/(sqrt3 H),  Omega_phi = x^2 + y^2,  w = (x^2 - y^2)/(x^2 + y^2)
  x' = -3x + sqrt(3/2) lam y^2 + 1.5 x (1 + x^2 - y^2)      (exponential V = V0 exp(-lam phi))
  y' = -sqrt(3/2) lam x y + 1.5 y (1 + x^2 - y^2)
  (ln H)' = -1.5 (1 + x^2 - y^2)
For other shapes lam is replaced by Gamma-dependent lam(phi) = -V'/V, carried as a third variable.
"""
import os, math
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq, minimize_scalar, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DS = os.path.join(HERE, "..", "dark-sector", "data")

# ---------- Pantheon+ ----------
P = os.path.join(DS, "pantheon")
tab = np.genfromtxt(os.path.join(P, "PantheonPlusSH0ES.dat"), names=True, dtype=None, encoding=None)
zhd, zhel, mb = tab["zHD"], tab["zHEL"], tab["m_b_corr"]
raw = np.loadtxt(os.path.join(P, "PantheonPlusSH0ES_STAT+SYS.cov"))
n = int(raw[0]); cov = raw[1:].reshape(n, n)
mask = zhd > 0.01
zhd, zhel, mb = zhd[mask], zhel[mask], mb[mask]
cov = cov[np.ix_(mask, mask)]
Cinv = np.linalg.inv(cov)
ones = np.ones(len(mb))
S11 = ones @ Cinv @ ones

# ---------- DESI DR2 BAO ----------
B = os.path.join(DS, "desi")
bz, bv, bq = [], [], []
for l in open(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_mean.txt")):
    if l.startswith("#") or not l.strip():
        continue
    a = l.split(); bz.append(float(a[0])); bv.append(float(a[1])); bq.append(a[2])
bz = np.array(bz); bv = np.array(bv)
bCinv = np.linalg.inv(np.loadtxt(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_cov.txt")))

ZG = np.concatenate([[0.0], np.logspace(-4, np.log10(3.0), 1500)])
Z_START = 10.0


def _rhs(N, s, shape, par):
    x, y, phi = s
    if shape == "exp":
        lam = par
    elif shape == "hilltop":          # V = V0 (1 - phi^2/(2 mu^2)) near the top of the hump; lam = -V'/V
        mu = par
        den = 1 - phi ** 2 / (2 * mu ** 2)
        if den < 1e-3:
            raise ValueError("rolled past the zero of V")
        lam = (phi / mu ** 2) / den
    elif shape == "cosine":           # V = V0 (1 + cos(phi/f)); lam = sin/(f(1+cos))
        f = par
        den = 1 + math.cos(phi / f)
        if den < 1e-4:
            raise ValueError("at the bottom of the circle")
        lam = math.sin(phi / f) / (f * den)
    else:
        raise ValueError(shape)
    if abs(x) > 1.2 or abs(y) > 1.2 or not math.isfinite(lam) or abs(lam) > 1e3:
        raise ValueError("field ran away")
    g = 1 + x * x - y * y
    return [-3 * x + math.sqrt(1.5) * lam * y * y + 1.5 * x * g,
            -math.sqrt(1.5) * lam * x * y + 1.5 * y * g,
            math.sqrt(6) * x]


def history(Om, shape="exp", par=0.0, phi0=0.0):
    """E(z) on the grid ZG for a thawing field (x = 0 at z = 10) with Omega_phi(0) = 1 - Om."""
    Ni, Ns = -math.log(1 + Z_START), -np.log(1 + ZG[::-1])

    def solve(yi):
        sol = solve_ivp(_rhs, (Ni, 0.0), [0.0, yi, phi0], args=(shape, par), t_eval=np.append(Ns[:-1], 0.0),
                        rtol=1e-8, atol=1e-11)
        return sol

    def miss(lyi):
        s = solve(math.exp(lyi))
        return s.y[0, -1] ** 2 + s.y[1, -1] ** 2 - (1 - Om)

    # Omega_phi at z=10 for a cosmological constant is (1-Om)/(Om*11^3 + 1-Om); bracket around it
    y_l = math.sqrt((1 - Om) / (Om * (1 + Z_START) ** 3 + 1 - Om))
    lo, w = math.log(y_l) - 3, 3.0
    # The top of the bracket (20x the Lambda value) can roll a hump-shaped field off its hump, and then the probe
    # itself fails even though the solution is fine.  Shrink the top until the probe runs (fix of Nov 26: before it,
    # whole regions of the hilltop and cosine fits scored 1e9 and the grid stopped at that wall).
    while True:
        try:
            m_hi = miss(math.log(y_l) + w)
            if not math.isfinite(m_hi): raise ValueError("probe not finite")
            break
        except Exception:
            w *= 0.7
            if w < 0.05: raise
    hi = math.log(y_l) + w
    lyi = brentq(miss, lo, hi, xtol=1e-12)
    s = solve(math.exp(lyi))
    x, y = s.y[0], s.y[1]
    g = 1 + x * x - y * y
    lnH = cumulative_trapezoid(-1.5 * g, s.t, initial=0.0)
    lnH -= lnH[-1]                      # H(0) = H0
    E = np.exp(lnH)[::-1]              # back to increasing z (ZG order)
    w = ((x * x - y * y) / (x * x + y * y))[::-1]
    phi = s.y[2][::-1]
    history.phi = phi                   # field value on the grid (reduced Planck units), z increasing
    return E, w


def chi2s(E):
    inv = 1.0 / E
    chi = cumulative_trapezoid(inv, ZG, initial=0.0)
    dm = np.interp(zhd, ZG, chi)
    mu = 5 * np.log10((1 + zhel) * dm)
    d = mb - mu
    sn = float(d @ Cinv @ d - (ones @ Cinv @ d) ** 2 / S11)
    DM = np.interp(bz, ZG, chi); Eb = np.interp(bz, ZG, E)
    DH = 1 / Eb; DV = (bz * DM ** 2 * DH) ** (1 / 3)
    t = np.array([{"DM_over_rs": DM[i], "DH_over_rs": DH[i], "DV_over_rs": DV[i]}[q] for i, q in enumerate(bq)])
    hrd = (t @ bCinv @ bv) / (t @ bCinv @ t)
    r = bv - hrd * t
    return sn, float(r @ bCinv @ r)


def fit_grid(shape, which, grid, phi0):
    """For the hump-shaped potentials: scan the shape parameter on a grid, Om fitted at each point."""
    best, ib, grid = None, 0, list(grid)
    for i, p in enumerate(grid):
        r = fit(shape, float(p), which, phi0=phi0)
        if best is None or r["chi2"] < best["chi2"]:
            best, ib = r, i
    # polish (Nov 26): Brent on the shape parameter between the best grid point's neighbours, Om refitted at each
    lo, hi = grid[max(ib - 1, 0)], grid[min(ib + 1, len(grid) - 1)]
    if hi > lo:
        rr = minimize_scalar(lambda p: fit(shape, float(p), which, phi0=phi0)["chi2"], bounds=(lo, hi),
                             method="bounded", options={"xatol": 1e-4})
        r = fit(shape, float(rr.x), which, phi0=phi0)
        if r["chi2"] < best["chi2"]:
            best = r
    best["k_fit"] = 2
    return best


def fit(shape="exp", par=None, which="both", par_bounds=(0.0, 3.0), phi0=0.0):
    """Best fit over Om (and par if par is None).  which = 'sn', 'bao' or 'both'.  Returns dict."""
    def total(Om, p):
        try:
            sn, bao = chi2s(history(Om, shape, p, phi0)[0])
        except Exception:
            return 1e9, 1e9, 1e9
        return {"sn": sn, "bao": bao, "both": sn + bao}[which], sn, bao

    if par is not None:
        r = minimize_scalar(lambda Om: total(Om, par)[0], bounds=(0.2, 0.45), method="bounded",
                            options={"xatol": 1e-4})
        c, sn, bao = total(r.x, par)
        return {"Om": float(r.x), "par": par, "chi2": c, "sn": sn, "bao": bao, "k_fit": 1}
    best = None
    for p0 in np.linspace(par_bounds[0], par_bounds[1], 7):
        r = minimize(lambda v: total(v[0], v[1])[0], [0.3, p0], method="Nelder-Mead",
                     options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 300})
        if best is None or r.fun < best.fun:
            best = r
    Om, p = best.x
    c, sn, bao = total(Om, p)
    return {"Om": float(Om), "par": float(p), "chi2": c, "sn": sn, "bao": bao, "k_fit": 2}


if __name__ == "__main__":
    E, w = history(0.3, "exp", 0.0)
    print("LCDM via field code, chi2 (sn, bao):", chi2s(E))
    E, w = history(0.3, "exp", 1.0)
    print("exp lam=1: w0 =", w[0], "chi2:", chi2s(E))
