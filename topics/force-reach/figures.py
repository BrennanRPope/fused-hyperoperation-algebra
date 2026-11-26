"""Pictures for the force-reach DRAFT. Run: python3 -I figures.py -> figures/"""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import special

R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})
hbarc = 197.3269804  # MeV fm
C = {"m0": "#1f77b4", "yuk": "#d62728", "tube": "#2ca02c", "pi": "#ff7f0e", "grey": "#777777"}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"figures/{name}.png", bbox_inches="tight")
    plt.close(fig)


# f1 headline: effective dimension against distance
r = np.logspace(-4, 1, 600)  # fm
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(r, 3 + 0 * r, color=C["m0"], lw=2.5, label="gravity, light (massless carrier): stays 3")
for name, mass, col in (("W", 80369.2, C["yuk"]), ("pion", 139.57, C["pi"])):
    mm = mass / hbarc
    ax.plot(r, 3 + (mm * r) ** 2 / (1 + mm * r), color=col, lw=2.5,
            label=f"{'weak force (W, Z)' if name=='W' else 'nuclear force (pion)'}: climbs, reach {hbarc/mass:.3g} fm")
a_c, sig = 4 / 3 * 0.3 * hbarc, 440.0 ** 2 / hbarc  # MeV fm, MeV/fm
F = a_c / r ** 2 + sig
Deff = 1 - np.gradient(np.log(F), np.log(r))
ax.plot(r, Deff, color=C["tube"], lw=2.5, label="quarks (flux tube): falls to 1")
ax.axvspan(1.1, 1.3, color=C["tube"], alpha=0.15)
ax.text(1.18, 5.3, "string breaks\n(real quarks)\n~1.1–1.25 fm", ha="center", fontsize=9, color=C["tube"])
ax.set_xscale("log"); ax.set_ylim(0, 7); ax.set_xlim(1e-4, 10)
ax.set_xlabel("distance r (fm)"); ax.set_ylabel("effective dimension  D_eff = 1 − d ln F / d ln r")
ax.set_title("f1. One reading for every force's reach: how many directions the flux still spreads into")
ax.axhline(1, color=C["grey"], ls=":", lw=1); ax.axhline(3, color=C["grey"], ls=":", lw=1)
ax.text(1.2e-4, 3.12, "3 = flux spreads over a sphere (1/r²)", fontsize=9, color=C["grey"])
ax.text(1.2e-4, 1.12, "1 = flux kept in a tube (constant pull)", fontsize=9, color=C["grey"])
ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(0, 0.93))
save(fig, "f1_effective_dimension")

# f2 the dial plane: massless line, screened side, closed side
k = np.linspace(-4, 4, 500); d = np.linspace(-3, 5, 500)
Kg, Dg = np.meshgrid(k, d)
Delta = Kg ** 2 - 4 * Dg
reach = np.where(Delta > 0, 2 / np.sqrt(np.abs(Delta) + 1e-300), np.nan)
mass_closed = np.where(Delta < 0, np.sqrt(np.abs(Delta)), np.nan)
fig, ax = plt.subplots(figsize=(8, 6))
cs = ax.contourf(Kg, Dg, np.log10(reach), levels=np.linspace(-0.8, 1.2, 21), cmap="Reds_r")
cb = fig.colorbar(cs, ax=ax); cb.set_label("screened side: log10 reach 2/√Δ")
ax.contourf(Kg, Dg, mass_closed, levels=20, cmap="Blues", alpha=0.6)
ax.plot(k, k ** 2 / 4, "k", lw=3)
ax.annotate("massless line Δ = 0 (d = k²/4)\ninfinite reach", xy=(2.0, 1.0), xytext=(0.6, -1.2), fontsize=10,
            bbox=dict(fc="white", ec="k"), arrowprops=dict(arrowstyle="->", lw=1.5))
ax.text(-3.7, -2.4, "Δ > 0 (line class)\ncarrier mass m = √Δ/2\nYukawa: screened", color=C["yuk"], fontsize=10,
        bbox=dict(fc="white", ec=C["yuk"]))
ax.text(-1.2, 3.6, "Δ < 0 (circle class)\nclosed direction: rungs\nmasses n√(−Δ), rung 0 massless", color=C["m0"],
        fontsize=10, bbox=dict(fc="white", ec=C["m0"]))
ax.set_xlabel("k (homothety dial)"); ax.set_ylabel("d (curvature dial)")
ax.set_title("f2. Reach is a reading of k and d: m² = Δ/4 = minus the curvature reading")
save(fig, "f2_dial_plane")

# f3 three ways field lines spread from a point (Jacobi R'' = -K R) and the pull that results
rr = np.linspace(0.001, 4, 400)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
for Kv, lab, col in ((1.0, "K > 0 (Δ<0): sin, refocuses", C["m0"]), (0.0, "K = 0 (Δ=0): r, Coulomb", "k"),
                     (-1.0, "K < 0 (Δ>0): sinh, spreads fast", C["yuk"])):
    Rr = np.sin(rr) if Kv > 0 else (rr if Kv == 0 else np.sinh(rr))
    axs[0].plot(rr, Rr, color=col, lw=2.5, label=lab)
    m = (Rr > 0.02)
    axs[1].semilogy(rr[m], 1 / (4 * np.pi * Rr[m] ** 2), color=col, lw=2.5, label=lab)
axs[0].axvline(np.pi / 2, color=C["m0"], ls=":"); axs[0].text(1.6, 0.2, "equator: area stops\ngrowing only here", fontsize=9,
                                                            color=C["m0"])
axs[0].set_ylim(-0.2, 4); axs[0].set_xlabel("distance r"); axs[0].set_ylabel("spread R(r) of field lines")
axs[0].set_title("A: spread obeys R'' = −K R, the law's own linear form")
axs[0].legend(fontsize=9)
axs[1].set_xlabel("distance r"); axs[1].set_ylabel("pull 1/(4πR²)")
axs[1].set_title("B: negative curvature screens: pull ~ e^{−2r}"); axs[1].set_ylim(1e-4, 10)
fig.suptitle("f3. Curvature of space alone already sets a reach (no carrier mass needed)")
save(fig, "f3_spreads")

# f4 Yukawa = flux x screening, and every dimension decays at the same rate
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
r = np.linspace(0.05, 6, 400); m = 1.0
axs[0].semilogy(r, 1 / (4 * np.pi * r), color="k", lw=2, label="flux part 1/(4πr)  (dimension dial, D=3)")
axs[0].semilogy(r, np.exp(-m * r), color=C["yuk"], lw=2, ls="--", label="screening part e^{−mr}  (m² = Δ/4)")
axs[0].semilogy(r, np.exp(-m * r) / (4 * np.pi * r), color=C["yuk"], lw=3, label="Yukawa = product")
axs[0].set_xlabel("r (units of 1/m)"); axs[0].set_ylabel("potential"); axs[0].legend(fontsize=9)
axs[0].set_title("A: two dials multiply (their log slopes add: −1/r − m)")
for D, col in ((2.5, "#9467bd"), (3, C["yuk"]), (4, "#8c564b"), (5, "#e377c2")):
    nu = D / 2 - 1
    G = (2 * np.pi) ** (-D / 2) * (m / r) ** nu * special.kv(nu, m * r)
    axs[1].semilogy(r, G * r ** ((D - 1) / 2), color=col, lw=2, label=f"D = {D}")
axs[1].set_xlabel("r (units of 1/m)"); axs[1].set_ylabel("G_D × r^{(D−1)/2}")
axs[1].set_title("B: after the power, every D falls at the same rate m"); axs[1].legend(fontsize=9)
fig.suptitle("f4. A massive carrier's field: the dimension dial gives the power, Δ gives the reach")
save(fig, "f4_yukawa_two_parts")

# f5 closed direction: rung 0 massless, rungs n massive
P = 2 * np.pi / np.sqrt(3); mu = np.sqrt(3)
r = np.logspace(-2, 1, 400)
G = 1 / np.tanh(np.pi * r / P) / (4 * np.pi * r * P)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [1.6, 1]})
axs[0].loglog(r, G, color="k", lw=3, label="field on open space × closed circle")
axs[0].loglog(r, 1 / (4 * np.pi ** 2 * r ** 2), color=C["grey"], ls="--", label="near: 4D, 1/r²")
axs[0].loglog(r, 1 / (4 * np.pi * r * P), color=C["m0"], ls="--", label="far: rung 0 only, 3D Coulomb 1/r")
for n in (1, 2):
    axs[0].loglog(r, 2 * np.exp(-n * mu * r) / (4 * np.pi * r * P), color=C["yuk"], alpha=0.8 - 0.3 * n,
                  label=f"rung {n}: e^(−{n}√3 r)/r, screened")
axs[0].set_ylim(1e-5, 1e3); axs[0].set_xlabel("distance r"); axs[0].set_ylabel("potential")
axs[0].legend(fontsize=9); axs[0].set_title("A: circle class k = 1, d = 1 (Δ = −3), circumference 2π/√3")
for n in range(0, 5):
    axs[1].hlines(n * mu, 0, 1, color=C["m0"] if n == 0 else C["yuk"], lw=3)
    axs[1].text(1.05, n * mu, f"rung {n}: mass {n}√3" + ("  (infinite reach)" if n == 0 else f", reach {1/(n*mu):.2f}"),
                va="center", fontsize=9)
axs[1].set_xlim(0, 3); axs[1].axis("off"); axs[1].set_title("B: the ladder seen from open space")
fig.suptitle("f5. A closed direction makes one massless carrier and a tower of short-range ones")
save(fig, "f5_closed_direction")

# f6 tube: pull of a charge in a square tube, D_eff 3 -> 1
z = np.array(R["tube_z"]); Fz = np.array(R["tube_F"])
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axs[0]
ax.add_patch(plt.Rectangle((-0.5, -3), 1, 6, fc="#e8f4e8", ec=C["tube"], lw=2))
ax.plot(0, 0, "ko", ms=8)
for th in np.linspace(-1.3, 1.3, 9):
    for sgn in (1, -1):
        t = np.linspace(0, 3, 200)
        x0 = 0.45 * np.sin(th)
        xs = x0 * (1 - np.exp(-3 * t)); ys = sgn * t
        ax.plot(xs, ys, color=C["tube"], lw=1)
ax.set_xlim(-2, 2); ax.set_ylim(-3, 3); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("A: flux kept inside a tube of side a\n(walls let none escape)")
axs[1].semilogx(z, 1 - np.gradient(np.log(Fz), np.log(z)), color=C["tube"], lw=3)
axs[1].axhline(3, color=C["grey"], ls=":"); axs[1].axhline(1, color=C["grey"], ls=":")
axs[1].set_xlabel("distance along the tube z / a"); axs[1].set_ylabel("D_eff")
axs[1].set_title("B: near: a sphere (3); far: one open direction (1)\npull → flux/(2·area), potential linear")
axs[1].text(0.012, 2.6, "sideways rungs screened\n(rate 2π/a): only rung 0 reaches far", fontsize=9)
fig.suptitle("f6. Confinement on the dimension dial: D = 1 gives a constant pull (string tension)")
save(fig, "f6_tube")

# f7 London tube: same mass sets width and tension
rho = np.linspace(0.01, 8, 400)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
for mm, col in ((0.5, "#9467bd"), (1.0, C["tube"]), (2.0, C["yuk"])):
    B = mm ** 2 / (2 * np.pi) * special.k0(mm * rho)
    axs[0].plot(rho, B, color=col, lw=2.5, label=f"m = {mm}: width (rms) 2/m = {2/mm:g}")
    axs[0].axvline(2 / mm, color=col, ls=":")
axs[0].set_ylim(0, 1.2); axs[0].set_xlabel("distance from the tube's centre ρ"); axs[0].set_ylabel("flux density B")
axs[0].legend(fontsize=9); axs[0].set_title("A: flux Φ = 1 squeezed by sideways screening, tail e^{−mρ}/√ρ")
ms = np.linspace(0.05, 3, 200); xi = 0.01
sig = ms ** 2 / (4 * np.pi) * (ms * xi) * special.k0(ms * xi) * special.k1(ms * xi)
axs[1].plot(ms, sig, color=C["tube"], lw=3, label="tension σ(m), core 0.01")
axs[1].plot(ms, ms ** 2 / (4 * np.pi) * np.log(1 / (ms * xi)), color=C["grey"], ls="--", label="(m²/4π) ln(1/(m ξ))")
axs[1].set_xlabel("screening mass m"); axs[1].set_ylabel("energy per length σ"); axs[1].legend(fontsize=9)
axs[1].set_title("B: the same m sets the tension, σ ~ m² × log")
fig.suptitle("f7. A tube made by screening: one mass gives both the width and the string tension")
save(fig, "f7_london_tube")

# f8 summary: two dials (and the switch)
fig, ax = plt.subplots(figsize=(9, 5.5))
ax.set_xlim(-0.3, 3.3); ax.set_ylim(0.3, 3.7)
ax.set_xlabel("carrier mass m  (a reading: m² = Δ/4 from k, d)")
ax.set_ylabel("open directions the flux can use  (dimension dial D)")
ax.set_xticks([0, 1.5, 3]); ax.set_xticklabels(["0", "small", "large"])
ax.set_yticks([1, 2, 3])
pts = [(0, 3, "gravity, light\nreach: unlimited", C["m0"]), (2.6, 3, "W, Z (weak)\nreach 0.002 fm", C["yuk"]),
       (1.3, 3, "pion (nuclear)\nreach 1.4 fm", C["pi"]), (0, 1, "quarks in a tube\nconstant pull", C["tube"])]
for x, y, t, c in pts:
    ax.plot(x, y, "o", ms=16, color=c); ax.text(x + 0.12, y - 0.05, t, fontsize=10, va="top")
ax.annotate("", xy=(0.15, 1.0), xytext=(1.5, 1.0), arrowprops=dict(arrowstyle="->", color=C["grey"], lw=1.5,
                                                                    ls="--"))
ax.plot(1.5, 1.0, "s", ms=12, color=C["grey"])
ax.text(1.62, 1.1, "dual-superconductor model: a mass m that\nscreens sideways, acting on flux that must be\nkept, makes the tube (width 2/m, tension ~ m²)",
        fontsize=9, color=C["grey"], va="top")
ax.set_title("f8. Two dials, not one: mass (from Δ) and open directions (D); a keep/let-end switch links them")
save(fig, "f8_two_dials")
print("ok")
