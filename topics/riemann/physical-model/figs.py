import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
O = json.load(open("phys.json")); B = np.load("burgers.npz"); C = np.load("collapse.npz")
snaps = np.load("../known-tools/rewind_snaps.npz")

# -------- Fig A: how a second zero moves near one at the centre ----------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 6.4))
Y, X = np.mgrid[-2:2:200j, -2:2:200j]; Z = X + 1j * Y
fields = [(1 / np.conj(Z), "A1. Charges pushing apart (2D Coulomb, overdamped):\nline charges in oil, superconductor vortices.\nPush outward in every direction", "C0"),
          (1j / np.conj(Z), "A2. Point vortices in a fluid:\nthe same field turned 90°.\nThey circle each other, never approach or flee", "C2"),
          (1 / Z, "A3. ζ's zeros under heat flow: the charge field MIRRORED\nin the wire. Push apart along the wire, pull together\nacross it, and swirl at 45°", "C3")]
for a, (W, title, c) in zip(ax, fields):
    a.streamplot(X, Y, W.real, W.imag, color=c, density=1.3, linewidth=0.8, arrowsize=0.9)
    a.axhline(0, color="k", lw=1.5); a.plot(0, 0, "ko", ms=8)
    a.set_aspect("equal"); a.set_xlim(-2, 2); a.set_ylim(-2, 2)
    a.set_title(title, fontsize=10); a.set_xlabel("along the wire (the critical line)")
ax[0].set_ylabel("across the wire (off the line)")
ax[2].annotate("pushed away", xy=(1.5, 0), xytext=(1.0, 0.45), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax[2].annotate("pulled in", xy=(0, 1.2), xytext=(0.35, 1.6), fontsize=9, arrowprops=dict(arrowstyle="->"))
plt.tight_layout(rect=(0, 0, 1, 0.97)); plt.savefig("figures/figA_fields.png", dpi=130); plt.close()

# -------- Fig B: the threshold in physical terms, and beads on a wire ----------------------------------
fig, ax = plt.subplots(1, 3, figsize=(18, 5.3))
a = ax[0]
S = 1.0; d = np.linspace(0.05, 3, 400)
a.plot(d, 4 / d, color="C3", label="the pair's own push (forward) / pull (rewind): 4/d")
a.plot(d, 2 * S * d, color="C0", label="the crowd's squeeze, a spring: 2·S·d")
ds = math.sqrt(2 / S); a.axvline(ds, color="k", ls=":"); a.plot(ds, 4 / ds, "ko")
a.text(ds + 0.06, 4.0, "threshold d* = √(2/S):\nthe two forces balance", fontsize=9)
a.fill_between(d[d < ds], 0, 12, color="C3", alpha=0.07); a.text(0.12, 10.5, "inside d*: on rewind the pair\nfalls together (bound)", fontsize=8.5, color="C3")
a.fill_between(d[d > ds], 0, 12, color="C0", alpha=0.07); a.text(1.75, 10.5, "outside d*: the crowd\ntears them apart", fontsize=8.5, color="C0")
a.set_ylim(0, 12); a.set_xlabel("gap d between two neighbouring zeros (S = 1)"); a.set_ylabel("force")
a.set_title("B1. The threshold is a Hill radius: where a pair's mutual pull\njust beats the crowd's tidal stretch (S = Σ 1/distance² to the rest)", fontsize=10)
a.legend(fontsize=8, loc=(0.42, 0.62))
a = ax[1]
rr = snaps["d0"] / np.sqrt(2 / snaps["S"])
a.hist(rr[(rr < 3)], bins=60, color="0.6")
a.axvline(1, color="k", ls=":"); a.text(1.55, 170, f"{O['hill']['share_inside']*100:.1f}% of our pairs\nsit inside their Hill radius", fontsize=9)
a.set_xlabel("gap ÷ Hill radius √(2/S), for each neighbouring pair of our 5,500 zeros"); a.set_ylabel("number of pairs")
a.set_title("B2. Our zeros sit right at the edge: the median gap is\n0.99 Hill radii, so about half are 'bound' and half are not", fontsize=10)
a = ax[2]
a.loglog(C["l1"], C["g1"], "o", color="C3", ms=4, label=f"1/r push (ζ's zeros, line charges, SC vortices): slope {O['collapse']['one_over_r']['exponent']:.3f}")
a.loglog(C["l2"], C["g2"], "s", color="C0", ms=4, label=f"1/r² push (charged beads on a wire): slope {O['collapse']['one_over_r2']['exponent']:.3f}")
a.set_xlabel("time left before two beads touch (rewound, in thick oil)"); a.set_ylabel("gap")
a.set_title("B3. Would beads on a wire do it? Only with the 1/r law.\nOrdinary charged beads (1/r²) collapse as (time left)^⅓, not ^½", fontsize=10)
a.legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.savefig("figures/figB_threshold_beads.png", dpi=130); plt.close()

# -------- Fig C: the exact physical home: Burgers shocks; and the dipole wire -------------------------------
fig, ax = plt.subplots(1, 3, figsize=(18, 5.3))
a = ax[0]
for y, c in [("0.8", "C0"), ("0.4", "C2"), ("0.2", "C1"), ("0.08", "C3")]:
    a.plot(B["xs"], B[f"u_{y}"], color=c, label=f"zero pair at 0.9 ± {y}i")
a.set_ylim(-14, 14); a.set_xlim(-0.5, 2.3)
a.set_xlabel("position x"); a.set_ylabel("gas velocity u(x) = −2 ∂ₓ log φ")
a.set_title("C1. Exact model: zeros of a heat solution are the POLES of a\nBurgers gas flow (Cole–Hopf). A pair nearing the line = a steepening front", fontsize=10)
a.legend(fontsize=8)
a = ax[1]
ys = np.array([0.8, 0.4, 0.2, 0.08]); sl = np.array([O["burgers_front_slope"][k] for k in ("0.8", "0.4", "0.2", "0.08")])
a.loglog(ys, sl, "o", color="C3", ms=8, label="measured steepest slope of the front")
yy = np.logspace(-1.2, 0, 50); a.loglog(yy, 4 / yy**2, "k--", lw=1, label="4 / (distance off the line)²")
a.set_xlabel("how far the zero pair sits off the line"); a.set_ylabel("steepest slope of the gas velocity")
a.set_title("C2. A zero landing on the line is a shock forming:\nthe front's slope blows up as 4/y². Heat-flow forward = Burgers run backward", fontsize=10)
a.legend(fontsize=8)
a = ax[2]
ph = np.linspace(0, 2 * math.pi, 300)
a.plot(np.degrees(ph), -np.cos(2 * ph), color="C3", lw=2, label="2D dipoles pointing along the wire: energy ∝ −cos 2φ / r²")
a.plot(np.degrees(ph[::15]), -np.cos(2 * ph[::15]), "ko", ms=4, label="the zeros' second-order (Calogero–Moser) energy, −Re 1/(z−w)²")
a.axhline(0, color="k", lw=0.5)
a.set_xticks([0, 90, 180, 270, 360]); a.set_xlabel("direction from one dipole to the other (0° = along the wire)")
a.set_ylabel("pair energy × r²")
a.set_title("C3. With inertia, the zeros are frictionless beads with dipole forces\n(attract along the wire, repel across), launched at exactly zero total energy", fontsize=10)
a.legend(fontsize=8, loc="lower right")
plt.tight_layout(); plt.savefig("figures/figC_exact_models.png", dpi=130); plt.close()
print("ok")
