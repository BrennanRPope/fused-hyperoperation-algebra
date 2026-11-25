"""Figures for the dark-sector page.  Run after sparc_fits.py and expansion_fits.py:  python3 -I figures.py"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sparc_fits as S
import expansion_fits as X

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = ["#257de0", "#f47230", "#00b894", "#f2a600",
                                                            "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"]
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.bbox": "tight", "savefig.facecolor": "white"})
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
SR = json.load(open(os.path.join(HERE, "sparc_results.json")))
ER = json.load(open(os.path.join(HERE, "expansion_results.json")))
P = {k: v.get("params", []) for k, v in SR["fixed_ml"].items()}

# ---------- d1: the radial acceleration relation, where the missing pull lives ----------
pts = np.load(os.path.join(HERE, "sparc_points.npz"))
gN, gobs = pts["gN"], pts["gobs"]; ok = gN > 0
fig, ax = plt.subplots(figsize=(7.4, 6))
ax.scatter(gN[ok], gobs[ok], s=4, color=INK2, alpha=0.25, lw=0, label="SPARC: 153 galaxies, 3,168 points")
x = np.logspace(-12.6, -8.3, 400); z = np.zeros_like(x)
ax.plot(x, x, color=INK, lw=1.5, ls="--", label="Newton, visible matter only")
ax.plot(x, S.rar(x, z, z, P["MOND / RAR (McGaugh 2016)"]), color=BLUE, label="MOND / RAR, one constant a0")
ax.plot(x, S.madd_exp(x, z, z, P["The circle operation of many small pieces (exp map)"]), color=ORANGE,
        label="The circle operation of the mass pieces (best scale)")
ax.plot(x, S.madd_background(x, z, z, P["The circle operation with a constant background"]), color=AQUA,
        label="The circle operation with a constant background")
a0 = 10 ** P["MOND / RAR (McGaugh 2016)"][0]
ax.axvline(a0, color=BLUE, lw=0.8, ls=":")
ax.text(a0 * 1.15, 2.2e-13, "a0 = 1.2e-10 m/s²", color=INK2, fontsize=9)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(2e-13, 5e-9); ax.set_ylim(1.5e-13, 5e-9)
ax.set_xlabel("pull of the visible matter, g_N (m/s²)"); ax.set_ylabel("pull we see from the orbit, g_obs (m/s²)")
ax.set_title("The missing pull appears where gravity is weak", fontsize=12, color=INK)
ax.legend(fontsize=9, loc="upper left")
fig.savefig(os.path.join(FIG, "d1_rar.png")); plt.close(fig)

# ---------- d2: four rotation curves ----------
names = ["NGC2403", "DDO154", "NGC3198", "NGC7814"]
fig, axs = plt.subplots(2, 2, figsize=(10, 7.4))
for ax, nm in zip(axs.flat, names):
    g = next(gg for gg in S.gals if gg["name"] == nm)
    R = g["R"]; Rf = np.linspace(R.min(), R.max(), 200)
    gi = {k: np.interp(Rf, R, v) for k, v in g.items() if k not in ("name",)}
    gi["name"] = nm
    gg, gd, gb = S.pieces(gi)
    V = lambda gm: np.sqrt(np.maximum(gm, 0) * Rf / S.KMS2_PER_KPC)
    ax.errorbar(R, g["V"], yerr=g["eV"], fmt="o", ms=4, color=INK, ecolor=INK2, elinewidth=1, label="measured")
    ax.plot(Rf, V(S.madd_exp(gg, gd, gb, P["The circle operation of many small pieces (exp map)"])), color=ORANGE, label="The circle operation of pieces")
    ax.plot(Rf, V(S.newton(gg, gd, gb, [])), color=INK, ls="--", lw=1.5, label="visible matter, Newton")
    ax.plot(Rf, V(S.rar(gg, gd, gb, P["MOND / RAR (McGaugh 2016)"])), color=BLUE, label="MOND / RAR")
    c2, q = S.fit_nfw(g)
    Vn = np.sqrt(np.maximum(gg * Rf / S.KMS2_PER_KPC + (gd + gb) * Rf / S.KMS2_PER_KPC + S.v_nfw2(Rf, 10 ** q[0], 10 ** q[1]), 0))
    ax.plot(Rf, Vn, color=VIOLET, label="dark-matter halo (NFW)")
    ax.set_title(nm, fontsize=11, color=INK); ax.set_xlabel("radius (kpc)"); ax.set_ylabel("speed (km/s)")
    ax.set_ylim(0, None)
axs[0, 0].legend(fontsize=8.5, loc="lower right")
fig.suptitle("Outer stars and gas move too fast for the matter we can see", fontsize=12, color=INK)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "d2_curves.png")); plt.close(fig)

# ---------- d3: the geometry (redshift rectangle, ΛCDM cube, square = rectangle) ----------
fig, axs = plt.subplots(1, 3, figsize=(13, 4.4))
ax = axs[0]; z1, z2 = 0.6, 0.4
ax.add_patch(Rectangle((0, 0), 1, 1, fc="#dbe9fa", ec=INK2)); ax.text(0.5, 0.5, "1", ha="center", va="center")
ax.add_patch(Rectangle((1, 0), z1, 1, fc="#fde3d4", ec=INK2)); ax.text(1 + z1 / 2, 0.5, "z₁", ha="center", va="center")
ax.add_patch(Rectangle((0, 1), 1, z2, fc="#ccf1e9", ec=INK2)); ax.text(0.5, 1 + z2 / 2, "z₂", ha="center", va="center")
ax.add_patch(Rectangle((1, 1), z1, z2, fc="#fbe6b3", ec=INK2)); ax.text(1 + z1 / 2, 1 + z2 / 2, "z₁z₂", ha="center", va="center")
ax.set_xlim(-0.1, 1.8); ax.set_ylim(-0.1, 1.55); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Redshifts compose by the circle operation\n(1+z₁)(1+z₂) = 1 + (z₁ ∘₁ z₂)", fontsize=10.5, color=INK)

ax = axs[1]                                     # oblique cube: unit cube (Λ) inside the (1+z) cube (matter shell)
def box(ax, o, s, fc, ec=INK2, alpha=1.0):
    (x0, y0), (a, b, c) = o, s
    d = 0.45
    front = [(x0, y0), (x0 + a, y0), (x0 + a, y0 + b), (x0, y0 + b)]
    top = [(x0, y0 + b), (x0 + a, y0 + b), (x0 + a + d * c, y0 + b + d * c), (x0 + d * c, y0 + b + d * c)]
    side = [(x0 + a, y0), (x0 + a + d * c, y0 + d * c), (x0 + a + d * c, y0 + b + d * c), (x0 + a, y0 + b)]
    for poly, sh in ((front, 1.0), (top, 0.85), (side, 0.7)):
        ax.add_patch(Polygon(poly, closed=True, fc=fc, ec=ec, alpha=alpha * sh, lw=1))
zz = 0.55
box(ax, (0, 0), (1 + zz, 1 + zz, 1 + zz), "#fde3d4")
box(ax, (0, 0), (1, 1, 1), "#dbe9fa")
ax.text(0.5, 0.45, "unit cube\n= Λ part", ha="center", va="center", fontsize=9)
ax.text(1.75, 2.05, "shell = 3z + 3z² + z³\nweighted by Ωm", ha="center", fontsize=9)
ax.set_xlim(-0.1, 2.4); ax.set_ylim(-0.1, 2.5); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Flat ΛCDM is a circle-operation cube with dial k\nH²/H₀² = 1 + Ωm·(z ∘₁ z ∘₁ z)", fontsize=10.5, color=INK)

ax = axs[2]; gNv, a0v = 0.35, 1.6; gv = np.sqrt(gNv * a0v)
ax.add_patch(Rectangle((0, 0), a0v, gNv, fc="#dbe9fa", ec=INK2)); ax.text(a0v / 2, gNv / 2, "g_N × a0", ha="center", va="center", fontsize=9)
ax.add_patch(Rectangle((0, 0.55), gv, gv, fc="#fde3d4", ec=INK2)); ax.text(gv / 2, 0.55 + gv / 2, "g × g", ha="center", va="center", fontsize=9)
ax.text(0.05 + a0v / 2, -0.15, "same area", ha="center", fontsize=9, color=INK2)
ax.set_xlim(-0.1, 1.8); ax.set_ylim(-0.25, 1.4); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Deep MOND: square on g = rectangle g_N × a0\n(an area statement, rung 2)", fontsize=10.5, color=INK)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "d3_geometry.png")); plt.close(fig)

# ---------- d4: supernova Hubble diagram, residuals against the coasting line ----------
mods = [("Flat LCDM (cosmological constant)", BLUE, "-", "ΛCDM"),
        ("Matter only, no dark energy (EdS)", RED, "-", "matter only (no dark energy)"),
        ("Empty open universe (Milne)", YELLOW, "-", "empty universe (Milne)"),
        ("Fused power law, chi = ((1+z)^b - 1)/b", ORANGE, "-", "fused power law"),
        ("Cube of scaled redshift, E^2 = (1 + l z)^3", AQUA, "-", "cube of scaled redshift")]
zs, zh, mb, Ci = X.zhd, X.zhel, X.mb, X.Cinv
ref = 5 * np.log10((1 + zh) * np.log1p(zs))
def offset(mu):
    d = mb - mu; return (X.ones @ Ci @ d) / X.S11
r0 = mb - ref - offset(ref)
bins = np.logspace(np.log10(0.01), np.log10(2.3), 16); idx = np.digitize(zs, bins)
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[1.5, 1]))
ax = axs[0]
bz = [np.median(zs[idx == i]) for i in range(1, len(bins)) if np.sum(idx == i) > 3]
br = [np.mean(r0[idx == i]) for i in range(1, len(bins)) if np.sum(idx == i) > 3]
be = [np.std(r0[idx == i]) / np.sqrt(np.sum(idx == i)) for i in range(1, len(bins)) if np.sum(idx == i) > 3]
ax.errorbar(bz, br, yerr=be, fmt="o", color=INK, ms=5, label="Pantheon+ (1,590 supernovae, binned)")
zf = np.logspace(-2, np.log10(2.3), 300)
for key, col, ls, lab in mods:
    model, p0, b = X.MODELS[key]; p = ER["fits"][key]["sn"]["params"]
    chi, dm = X.chi_of(model, p, zs); mu = 5 * np.log10((1 + zh) * dm); off = offset(mu)
    chif, dmf = X.chi_of(model, p, zf)
    ax.plot(zf, 5 * np.log10((1 + zf) * dmf) + off - (5 * np.log10((1 + zf) * np.log1p(zf)) + offset(ref)), color=col, ls=ls, label=lab)
ax.axhline(0, color=INK2, lw=1, ls="--")
ax.set_xscale("log"); ax.set_xlabel("redshift z"); ax.set_ylabel("brightness vs. a coasting universe (mag)")
ax.set_title("Supernovae: distant ones are fainter than any decelerating model", fontsize=11, color=INK)
ax.legend(fontsize=8.5, loc="lower left")
ax = axs[1]
zq = np.linspace(0, 2.5, 200)
for key, col, ls, lab in mods:
    model, p0, b = X.MODELS[key]; p = ER["fits"][key]["sn+bao"]["params"]
    ax.plot(zq, X.q_of_z(model, p, zq), color=col, ls=ls, label=lab)
ax.axhline(0, color=INK2, lw=1, ls="--"); ax.text(1.6, -0.12, "below 0 = speeding up", fontsize=9, color=INK2)
ax.set_xlabel("redshift z (looking back in time →)"); ax.set_ylabel("deceleration q")
ax.set_title("Only ΛCDM switches from slowing to speeding up", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "d4_supernovae.png")); plt.close(fig)

# ---------- d5: DESI BAO ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
zf = np.linspace(0.05, 2.5, 200)
sig = np.sqrt(np.diag(X.bcov))
for ax, q, lab in ((axs[0], "DM_over_rs", "D_M / r_d  (sideways ruler)"), (axs[1], "DH_over_rs", "D_H / r_d  (line-of-sight ruler)")):
    sel = np.array([qq == q for qq in X.bq])
    scale = (lambda zz: zz ** (-2 / 3)) if q == "DM_over_rs" else (lambda zz: zz ** 0)
    for key, col, ls, l2 in mods:
        model, p0, b = X.MODELS[key]; p = ER["fits"][key]["sn+bao"]["params"]; hrd = ER["fits"][key]["bao_hrd"]
        chi, dm = X.chi_of(model, p, zf)
        E = model[1](zf, p) if model[0] == "E" else X.EH[model](zf, p)
        y = hrd * dm if q == "DM_over_rs" else hrd / E
        ax.plot(zf, y * scale(zf), color=col, ls=ls, label=l2)
    ax.errorbar(X.bz[sel], X.bv[sel] * scale(X.bz[sel]), yerr=sig[sel] * scale(X.bz[sel]), fmt="o", color=INK, ms=5, label="DESI DR2")
    ax.set_xlabel("redshift z"); ax.set_ylabel(lab + (" × z^-2/3" if q == "DM_over_rs" else ""))
axs[0].legend(fontsize=8.5)
fig.suptitle("Baryon acoustic rulers (DESI DR2), each model fit jointly with the supernovae", fontsize=12, color=INK)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "d5_bao.png")); plt.close(fig)
print("figures written:", sorted(os.listdir(FIG)))
