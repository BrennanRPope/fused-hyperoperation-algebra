"""DRAFT picture for round 25. Run after checks.py. Blue = the law, orange = the weak force's twist and its energy, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json")); FG = R["figure"]
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()

a = ax[0]
f = np.linspace(-0.25, 2.25, 400)
a.plot(f, 2 * f - f**2, color=BL, lw=2.4, label="bending of the 3-sphere this way of carrying gives (sectional curvature 2f − f²)")
a.plot(f, 2 * f - 2, color=GR, lw=1.5, ls="--", label="twist (torsion 2f − 2)")
a.plot(f, (f**2 - 2 * f)**2 / 2, color=OR, lw=2.4, label="energy of the weak force's field A = f σᵃTₐ at rest (barrier)")
lam = f / 2; a.plot(f, lam**2 * (3 - 2 * lam), color=OR, lw=1.5, ls=":", label="Chern–Simons number λ²(3 − 2λ), λ = f/2")
for x0, txt in ((0, "f = 0\nflat, twisted\n(left multiplication)\nweak vacuum, CS 0"), (2, "f = 2\nflat, twisted\n(right multiplication)\nweak vacuum, CS 1")):
    a.plot([x0], [0], "o", color=GR, ms=8); a.text(x0, -1.05, txt, ha="center", fontsize=8.5, color=GR)
a.plot([1], [1], "o", color=BL, ms=11); a.plot([1], [0.5], "o", color=OR, ms=9)
a.annotate("f = 1: the law's own transport\n(Levi-Civita, no twist, curvature 1)\n= the 'meron': top of the barrier, CS ½", (1, 1), (0.25, 1.3), fontsize=9, color=BL,
           arrowprops=dict(arrowstyle="->", color=BL))
a.set_xlabel("f, the setting for how a direction is carried on the 3-sphere (∇_X Y = f [X, Y]/2)"); a.set_ylim(-1.35, 2.75)
a.axhline(0, color=GR, lw=0.5); a.legend(fontsize=8.5, loc="upper left")
a.set_title("A. The law sits halfway between the two flat twists (K1, K2)")

b = ax[1]
ds = FG["ds"]; t = np.array(ds["t"]); fv = np.array(ds["f"]); et = np.array(ds["eta"]); E4 = np.array(ds["E4"])
b.plot(t, fv, color=OR, lw=2.4, label="the locked twist f(t) on the law's de Sitter (cosmic time)")
b.plot(t, 1 + 0 * t, color=BL, lw=1, ls=":", label="f = 1, the law's own transport (hilltop)")
b.plot(t, et / math.pi + 0.0, color=GR, lw=1.5, ls="--", label="conformal time η/π: only π of it in all of history")
b2 = b.twinx(); b2.plot(t, E4 / E4[len(E4) // 2], color=BL, lw=1.2, alpha=0.6); b2.set_ylim(0.9, 1.1)
b2.set_ylabel("a⁴ρ, the twist's energy × a⁴ (kept to 1e-9; blue line)", color=BL)
b.text(0.03, 0.40, f"field equation in conformal time: ψ'' = 2ψ(1 − ψ²), f = 1 + ψ\nso it makes {R['K2']['periods_in_all_time']:.2f} of a swing, then freezes as space grows.\nIts pressure is ρ/3 in every direction, like light.", transform=b.transAxes, fontsize=9, color=GR)
b.set_xlabel("cosmic time t (q = 1)"); b.set_ylabel("f"); b.legend(fontsize=8.5, loc="upper left")
b.set_title("B. The weak force's field locked to the law's expanding slices (K2, K3)")

c = ax[2]
q = FG["q"]
for nm, col, ls in (("ball", OR, "-"), ("shell", BL, "-")):
    B = FG["bubbles"][nm]; r = np.array(B["r"]); s = r > 0.03
    c.plot(r[s], np.array(B["qloc"])[s], color=col, lw=1.6, ls="--", label=f"{nm}: local dial q + 8πρ/3 (what the law would set from the twist at r)")
    c.plot(r[s], np.array(B["qbar"])[s], color=col, lw=2.4, ls=ls, label=f"{nm}: averaged dial q + 2m/r³ (the curvature across the sphere at r)")
c.axhline(q, color=GR, lw=1, ls=":", label=f"the law's de Sitter dial q = {q}")
c.set_xlim(0, 3.2); c.set_yscale('log'); c.set_ylim(0.03, 12); c.set_xlabel("radius r"); c.set_ylabel("dial (curvature)")
c.text(1.75, 0.6, "same total twist energy in both\n(M = %.2f): outside, the averaged dial\nfalls as q + 2M/r³ whatever the shape" % FG["bubbles"]["ball"]["M"], fontsize=9, color=GR)
c.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 1.0))
c.set_title("C. A bubble of twist raises the dial, unevenly (K5)")

d = ax[3]
for nm, col in (("ball", OR), ("shell", BL)):
    B = FG["bubbles"][nm]; r = np.array(B["r"]); s = r > 0.03
    tide_e = 1.5 * (np.array(B["qbar"]) - np.array(B["qloc"]))
    d.plot(r[s], tide_e[s], color=col, lw=2.4, label=f"{nm}: 3/2 × (averaged dial − local dial), from the twist's energy")
    d.plot(r[s][::3], np.array(B["tide"])[s][::3], "o", color=col, ms=3.5, mfc="none", label=f"{nm}: tide read from the metric")
M = FG["bubbles"]["ball"]["M"]; rr = np.linspace(1.2, 3.2, 200)
d.plot(rr, 3 * M / rr**3, color=GR, lw=1.6, ls="--", label="3M/r³: the textbook tide of a mass M (Newton, Schwarzschild)")
d.axhline(0, color=GR, lw=0.5)
d.text(0.97, 0.05, "no tide inside the hollow shell (Newton's shell theorem)\nno tide anywhere for an even twist (K3, K5a): Schur's lemma\ntide = −(r/2) × slope of the averaged dial (K4)", transform=d.transAxes, fontsize=9, color=GR, va="bottom", ha="right")
d.set_xlim(0, 3.2); d.set_xlabel("radius r"); d.set_ylabel("tidal stretch: curvature across − along the radius")
d.legend(fontsize=8, loc="upper right")
d.set_title("D. Twist → dial → tide: the tide is the dial's slope (K4, K5)")

fig.suptitle("DRAFT round 25: the weak force's twist and gravity's tidal pull, together", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("figures/f32_twist_dial_tide.png", bbox_inches="tight")
