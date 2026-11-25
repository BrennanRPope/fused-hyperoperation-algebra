"""Galaxy rotation curves (SPARC) against Newton, fused (the circle operation) gravity, MOND and NFW dark-matter halos.

Run:  python3 -I sparc_fits.py  (writes sparc_results.json and sparc_points.npz next to this file)

Accelerations are in m/s^2.  g_N = Newtonian pull of the visible matter (gas + stars), g_obs = V_obs^2 / R.
Models map g_N -> g_model and predict V = sqrt(R * g_model).
"""
import json, os, math
import numpy as np
from scipy.optimize import minimize_scalar, minimize, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data", "sparc")
KMS2_PER_KPC = 1e6 / 3.0857e19          # (km/s)^2/kpc in m/s^2
UD, UB = 0.5, 0.7                         # stellar mass-to-light at 3.6 micron (disk, bulge), McGaugh+2016

# ---------- load ----------
lines = open(os.path.join(D, "SPARC_Lelli2016c.mrt")).read().splitlines()
start = max(i for i, l in enumerate(lines) if l.startswith("-----")) + 1
meta = {}
for l in lines[start:]:
    if not l.strip():
        continue
    f = l.split()
    meta[f[0]] = dict(D=float(f[2]), inc=float(f[5]), Q=int(f[17]))

gals = []
for name, m in meta.items():
    if m["Q"] >= 3 or m["inc"] < 30:
        continue
    a = np.loadtxt(os.path.join(D, "rotmod", f"{name}_rotmod.dat"), comments="#")
    a = np.atleast_2d(a)
    R, V, eV, Vg, Vd, Vb = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
    ok = (R > 0) & (V > 0) & (eV > 0)
    gals.append(dict(name=name, R=R[ok], V=V[ok], eV=eV[ok], Vg=Vg[ok], Vd=Vd[ok], Vb=Vb[ok]))

def pieces(g, ud=UD, ub=UB):
    """Newtonian accelerations of gas, disk, bulge (m/s^2). Gas V^2 is signed (can be negative near the centre)."""
    R = g["R"]
    gg = g["Vg"] * np.abs(g["Vg"]) / R * KMS2_PER_KPC
    gd = ud * g["Vd"] ** 2 / R * KMS2_PER_KPC
    gb = ub * g["Vb"] ** 2 / R * KMS2_PER_KPC
    return gg, gd, gb

# ---------- models: g_N pieces -> g_model ----------
def newton(gg, gd, gb, p):            return gg + gd + gb
def madd_exp(gg, gd, gb, p):           # circle operation over infinitely many small mass pieces = exp map of the formal group
    s = 10 ** p[0]; x = (gg + gd + gb) / s
    return s * np.expm1(np.clip(x, -50, 50))
def madd_pieces(gg, gd, gb, p):        # circle operation of the three named pieces: gas (+) disk (+) bulge
    s = 10 ** p[0]
    return s * ((1 + gg / s) * (1 + gd / s) * (1 + gb / s) - 1)
def madd_background(gg, gd, gb, p):    # local field (+) a constant background k on scale s: g = (1+k) gN + k s
    k, s = p[0], 10 ** p[1]
    return (1 + k) * (gg + gd + gb) + k * s
def mond_simple(gg, gd, gb, p):        # mu(x) = x/(1+x) = minus the circle inverse of x
    a0 = 10 ** p[0]; gN = np.maximum(gg + gd + gb, 1e-16)
    return gN / 2 + np.sqrt(gN ** 2 / 4 + gN * a0)
def rar(gg, gd, gb, p):                # McGaugh+2016:  g = gN / (1 - exp(-sqrt(gN/a0)))
    a0 = 10 ** p[0]; gN = np.maximum(gg + gd + gb, 1e-16)
    return gN / (-np.expm1(-np.sqrt(gN / a0)))

GLOBAL = {
    "Newton (visible matter only)": (newton, []),
    "The circle operation of many small pieces (exp map)": (madd_exp, [(-14, -2)]),
    "The circle operation of gas, disk, bulge": (madd_pieces, [(-14, -2)]),
    "The circle operation with a constant background": (madd_background, [(-1, 5), (-14, -8)]),
    "MOND, simple mu (= minus circle inverse)": (mond_simple, [(-12, -8)]),
    "MOND / RAR (McGaugh 2016)": (rar, [(-12, -8)]),
}

def chi2_gal(g, f, p, ud=UD, ub=UB):
    gg, gd, gb = pieces(g, ud, ub)
    gm = f(gg, gd, gb, p)
    Vm = np.sqrt(np.maximum(gm, 0) * g["R"] / KMS2_PER_KPC)
    return np.sum(((g["V"] - Vm) / g["eV"]) ** 2)

def total_chi2(f, p):
    return sum(chi2_gal(g, f, p) for g in gals)

# ---------- NFW halo, two free numbers per galaxy ----------
H0_KPC = 0.073                         # km/s/kpc  (73 km/s/Mpc, SPARC's distance scale)
def v_nfw2(R, V200, c):
    R200 = V200 / (10 * H0_KPC)
    x = R / R200
    m = lambda y: np.log1p(y) - y / (1 + y)
    return V200 ** 2 * m(c * x) / (x * m(c))

def fit_nfw(g, ud=UD, ub=UB, free_ml=False):
    gg, gd, gb = pieces(g, 1.0, 1.0)
    R = g["R"]; base_g = gg * R / KMS2_PER_KPC; base_d = gd * R / KMS2_PER_KPC; base_b = gb * R / KMS2_PER_KPC
    def res(q):
        V200, c = 10 ** q[0], 10 ** q[1]
        u_d, u_b = (q[2], q[2] * UB / UD) if free_ml else (ud, ub)
        V2 = base_g + u_d * base_d + u_b * base_b + v_nfw2(R, V200, c)
        return (g["V"] - np.sqrt(np.maximum(V2, 0))) / g["eV"]
    best = None
    for lv in (1.5, 2.0, 2.5):
        for lc in (0.3, 0.8, 1.3):
            q0 = [lv, lc] + ([UD] if free_ml else [])
            lo = [1.0, 0.0] + ([0.05] if free_ml else []); hi = [2.8, 2.0] + ([1.5] if free_ml else [])
            r = least_squares(res, q0, bounds=(lo, hi))
            c2 = float(np.sum(r.fun ** 2))
            if best is None or c2 < best[0]:
                best = (c2, r.x.tolist())
    return best

def fit_free_ml(g, f, p):
    """Best stellar mass-to-light per galaxy (disk; bulge scaled with it), 0.05..1.5."""
    r = minimize_scalar(lambda u: chi2_gal(g, f, p, u, u * UB / UD), bounds=(0.05, 1.5), method="bounded")
    return float(r.fun), float(r.x)

if __name__ == "__main__":
    npts = sum(len(g["R"]) for g in gals); ng = len(gals)
    print(f"{ng} galaxies (Q<3, inc>=30), {npts} points")
    out = dict(n_galaxies=ng, n_points=npts, fixed_ml={}, free_ml={})

    # global fits, fixed mass-to-light
    for label, (f, bounds) in GLOBAL.items():
        if not bounds:
            c2 = total_chi2(f, []); p = []
        elif len(bounds) == 1:
            grid = np.linspace(*bounds[0], 61)
            vals = [total_chi2(f, [x]) for x in grid]
            i = int(np.argmin(vals))
            lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
            r = minimize_scalar(lambda x: total_chi2(f, [x]), bounds=(lo, hi), method="bounded")
            c2, p = float(r.fun), [float(r.x)]
        else:
            best = None
            for k0 in (0.05, 0.3, 1.0, 3.0):
                for s0 in (-12.0, -10.5, -9.0):
                    r = minimize(lambda q: total_chi2(f, q), [k0, s0], method="Nelder-Mead",
                                 options=dict(xatol=1e-5, fatol=1e-4, maxiter=8000, initial_simplex=[[k0, s0], [k0 * 1.5 + 0.1, s0], [k0, s0 + 0.5]]))
                    if best is None or r.fun < best.fun:
                        best = r
            c2, p = float(best.fun), best.x.tolist()
        k = len(p)
        out["fixed_ml"][label] = dict(chi2=c2, params=p, k=k, bic=c2 + k * math.log(npts),
                                      chi2_per_dof=c2 / (npts - k))
        print(f"{label:48s} chi2={c2:12.1f}  k={k:3d}  params={['%.3f' % x for x in p]}")

    # NFW: 2 per galaxy
    nfw = [fit_nfw(g) for g in gals]
    c2 = sum(b[0] for b in nfw); k = 2 * ng
    out["fixed_ml"]["NFW dark-matter halo (2 per galaxy)"] = dict(chi2=c2, k=k, bic=c2 + k * math.log(npts),
                                                                   chi2_per_dof=c2 / (npts - k))
    print(f"{'NFW halo':48s} chi2={c2:12.1f}  k={k}")

    # free mass-to-light per galaxy
    for label in ["Newton (visible matter only)", "The circle operation of many small pieces (exp map)", "MOND / RAR (McGaugh 2016)"]:
        f, _ = GLOBAL[label]; p = out["fixed_ml"][label]["params"]
        if label.startswith("The circle operation"):          # refit the scale together with the per-galaxy M/L (coarse grid)
            trials = [[x] for x in np.linspace(-11, -2, 10)]
            scored = [(sum(fit_free_ml(g, f, q)[0] for g in gals), q) for q in trials]
            p = min(scored)[1]
        res = [fit_free_ml(g, f, p) for g in gals]
        c2 = sum(r[0] for r in res); k = ng + len(p)
        out["free_ml"][label] = dict(chi2=c2, k=k, params=list(map(float, p)), bic=c2 + k * math.log(npts), chi2_per_dof=c2 / (npts - k),
                                     median_ml=float(np.median([r[1] for r in res])),
                                     frac_at_upper_bound=float(np.mean([r[1] > 1.49 for r in res])))
        print(f"free M/L  {label:38s} chi2={c2:12.1f} k={k} median M/L={out['free_ml'][label]['median_ml']:.2f}")
    nfw3 = [fit_nfw(g, free_ml=True) for g in gals]
    c2 = sum(b[0] for b in nfw3); k = 3 * ng
    out["free_ml"]["NFW dark-matter halo"] = dict(chi2=c2, k=k, bic=c2 + k * math.log(npts), chi2_per_dof=c2 / (npts - k))
    print(f"free M/L  NFW  chi2={c2:.1f} k={k}")

    # Circle-operation exp map: how strong must the scale be? chi2 profile vs log10 a*
    grid = np.linspace(-12, -2, 41)
    out["madd_exp_profile"] = dict(log_s=grid.tolist(), chi2=[total_chi2(madd_exp, [x]) for x in grid])

    # where does the missing pull live? fraction of points with g_obs > 2 g_N, binned by g_N
    gN_all, gobs_all, eg_all = [], [], []
    for g in gals:
        gg, gd, gb = pieces(g)
        gN_all.append(gg + gd + gb); gobs_all.append(g["V"] ** 2 / g["R"] * KMS2_PER_KPC)
        eg_all.append(2 * g["eV"] / g["V"])
    gN_all = np.concatenate(gN_all); gobs_all = np.concatenate(gobs_all); eg_all = np.concatenate(eg_all)
    np.savez(os.path.join(HERE, "sparc_points.npz"), gN=gN_all, gobs=gobs_all, rel_err=eg_all)
    ok = gN_all > 0
    ratio = gobs_all[ok] / gN_all[ok]
    bins = np.logspace(-12.5, -8.5, 9)
    idx = np.digitize(gN_all[ok], bins)
    out["mass_discrepancy_by_gN"] = [dict(gN_lo=float(bins[i - 1]) if i > 0 else 0.0,
                                          gN_hi=float(bins[i]) if i < len(bins) else float("inf"),
                                          median_ratio=float(np.median(ratio[idx == i])), n=int(np.sum(idx == i)))
                                     for i in range(len(bins) + 1) if np.sum(idx == i) > 0]
    for b in out["mass_discrepancy_by_gN"]:
        print(f"gN {b['gN_lo']:.1e}-{b['gN_hi']:.1e}: median g_obs/g_N = {b['median_ratio']:.2f} (n={b['n']})")
    json.dump(out, open(os.path.join(HERE, "sparc_results.json"), "w"), indent=1)
