"""Cosmic expansion: Pantheon+ supernovae and DESI DR2 BAO against LCDM and fused-algebra expansion laws.

Run:  python3 -I expansion_fits.py   (writes expansion_results.json)

Every model gives E(z) = H(z)/H0 (or a comoving distance directly).  Supernovae fix only the SHAPE of the distance-redshift
curve: the absolute magnitude M and H0 are degenerate and are marginalised analytically.  BAO fixes the shape too, with the
product H0*r_d as one nuisance number.
"""
import json, os, math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
C = 299792.458

# ---------- Pantheon+ ----------
P = os.path.join(HERE, "data", "pantheon")
cols = open(os.path.join(P, "PantheonPlusSH0ES.dat")).readline().split()
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

ZG = np.concatenate([[0.0], np.logspace(-4, np.log10(3.0), 4000)])

def chi_of(model, p, z):
    """Dimensionless comoving distance  chi = (H0/c) * D_C, and transverse distance D_M (curvature aware)."""
    kind, f = model
    if kind == "E":
        inv = 1.0 / f(ZG, p)
        chi = cumulative_trapezoid(inv, ZG, initial=0.0)
        return np.interp(z, ZG, chi), np.interp(z, ZG, chi)
    return f(z, p)

def sn_chi2(model, p):
    chi, dm = chi_of(model, p, zhd)
    if np.any(~np.isfinite(dm)) or np.any(dm <= 0):
        return 1e12
    mu = 5 * np.log10((1 + zhel) * dm)
    d = mb - mu
    return float(d @ Cinv @ d - (ones @ Cinv @ d) ** 2 / S11)

# ---------- DESI DR2 BAO ----------
B = os.path.join(HERE, "data", "desi")
bz, bv, bq = [], [], []
for l in open(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_mean.txt")):
    if l.startswith("#") or not l.strip():
        continue
    a = l.split(); bz.append(float(a[0])); bv.append(float(a[1])); bq.append(a[2])
bz = np.array(bz); bv = np.array(bv)
bcov = np.loadtxt(os.path.join(B, "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
bCinv = np.linalg.inv(bcov)

def bao_pred(model, p, hrd):
    """hrd = c / (H0 r_d), the one BAO nuisance number."""
    kind, f = model
    chi, dm = chi_of(model, p, bz)
    E = f(bz, p) if kind == "E" else EH[model](bz, p)
    DM = hrd * dm; DH = hrd / E; DV = (bz * DM ** 2 * DH) ** (1 / 3)
    return np.array([{"DM_over_rs": DM[i], "DH_over_rs": DH[i], "DV_over_rs": DV[i]}[q] for i, q in enumerate(bq)])

def bao_chi2(model, p):
    # chi2 is quadratic in hrd: minimise analytically
    t = bao_pred(model, p, 1.0)
    if np.any(~np.isfinite(t)):
        return 1e12, float("nan")
    hrd = (t @ bCinv @ bv) / (t @ bCinv @ t)
    d = bv - hrd * t
    return float(d @ bCinv @ d), float(hrd)

# ---------- models ----------
def E_lcdm(z, p):  return np.sqrt(p[0] * (1 + z) ** 3 + 1 - p[0])
def E_wcdm(z, p):  return np.sqrt(p[0] * (1 + z) ** 3 + (1 - p[0]) * (1 + z) ** (3 * (1 + p[1])))
def E_eds(z, p):   return (1 + z) ** 1.5
def E_power(z, p): return (1 + z) ** (1 - p[0])          # fused power law: chi = ((1+z)^b - 1)/b = b-th circle power / b
def E_cube(z, p):                                         # circle-operation cube of redshift with its own weight per layer
    return np.sqrt(np.maximum(1 + p[0] * 3 * z + p[1] * 3 * z ** 2 + p[2] * z ** 3, 1e-6))

def E_scaled_cube(z, p):                                  # cube of the scaled redshift: E^2 = (1 + l z)^3 = 1 + (l z)(+)(l z)(+)(l z)
    return (1 + p[0] * z) ** 1.5

def D_milne(z, p):                                        # empty open universe: D_M = sinh(ln(1+z)) = ((1+z)^2-1)/(2(1+z))
    chi = np.log1p(z); return chi, np.sinh(chi)
def D_coast(z, p):                                        # flat coasting (a ~ t): chi = ln(1+z), the linear (log) variable
    chi = np.log1p(z); return chi, chi

MODELS = {
    "Flat LCDM (cosmological constant)": (("E", E_lcdm), [0.3], [(0.01, 1.0)]),
    "Flat wCDM (dark energy with w free)": (("E", E_wcdm), [0.3, -1.0], [(0.01, 1.0), (-3, 0.3)]),
    "Matter only, no dark energy (EdS)": (("E", E_eds), [], []),
    "Empty open universe (Milne)": (("D", D_milne), [], []),
    "Flat coasting, chi = ln(1+z)": (("D", D_coast), [], []),
    "Fused power law, chi = ((1+z)^b - 1)/b": (("E", E_power), [0.3], [(-1.0, 1.5)]),
    "Cube of scaled redshift, E^2 = (1 + l z)^3": (("E", E_scaled_cube), [0.3], [(0.0, 3.0)]),
    "Circle-operation cube, free weight per layer": (("E", E_cube), [0.3, 0.3, 0.3], [(-2, 3), (-2, 3), (-2, 3)]),
}
EH = {MODELS["Empty open universe (Milne)"][0]: lambda z, p: 1 + z,
      MODELS["Flat coasting, chi = ln(1+z)"][0]: lambda z, p: 1 + z}

def fit(model, p0, bounds, which):
    def obj(q):
        if any(not (lo <= x <= hi) for x, (lo, hi) in zip(q, bounds)):
            return 1e12
        s = 0.0
        if "sn" in which: s += sn_chi2(model, q)
        if "bao" in which: s += bao_chi2(model, q)[0]
        return s
    if not p0:
        return obj([]), []
    best = None
    for jitter in (0.0, 0.1, -0.1):
        r = minimize(obj, [x + jitter for x in p0], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-4, maxiter=6000))
        if best is None or r.fun < best.fun:
            best = r
    return float(best.fun), best.x.tolist()

def q_of_z(model, p, z):
    """Deceleration parameter q = -1 + (1+z) dlnE/dz (negative = accelerating)."""
    kind, f = model
    if kind == "D":
        return np.zeros_like(z)            # both coasting models have q = 0
    h = 1e-4
    return -1 + (1 + z) * (np.log(f(z + h, p)) - np.log(f(z - h, p))) / (2 * h)

if __name__ == "__main__":
    out = dict(n_sn=int(mask.sum()), n_bao=len(bv), fits={})
    print(f"Pantheon+: {mask.sum()} SNe with z>0.01;  DESI DR2: {len(bv)} BAO numbers")
    for label, (model, p0, bounds) in MODELS.items():
        row = {}
        for which, npts in (("sn", int(mask.sum())), ("bao", len(bv)), ("sn+bao", int(mask.sum()) + len(bv))):
            c2, p = fit(model, p0, bounds, which)
            nuis = (1 if "sn" in which else 0) + (1 if "bao" in which else 0)
            k = len(p) + nuis
            row[which] = dict(chi2=c2, params=p, k=k, bic=c2 + k * math.log(npts), aic=c2 + 2 * k)
        _, p = fit(model, p0, bounds, "sn+bao")
        z = np.array([0.0, 0.5, 1.0, 2.0])
        row["q_at_z_0_0.5_1_2"] = q_of_z(model, p, z).tolist()
        row["bao_hrd"] = bao_chi2(model, p)[1]
        out["fits"][label] = row
        print(f"{label:42s} SN chi2={row['sn']['chi2']:8.1f} p={['%.3f' % x for x in row['sn']['params']]}  "
              f"BAO chi2={row['bao']['chi2']:7.1f}  joint={row['sn+bao']['chi2']:8.1f} p={['%.3f' % x for x in row['sn+bao']['params']]}")
    json.dump(out, open(os.path.join(HERE, "expansion_results.json"), "w"), indent=1)
