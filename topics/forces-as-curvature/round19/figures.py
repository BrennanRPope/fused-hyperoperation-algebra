"""DRAFT picture for round 19. Run after checks.py. Blue = the law, orange = textbook special relativity, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 4, figsize=(23, 6.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1, 1.1]))

a = ax[0]
labels = ["3D law vs\nB·Y·B", "plain X·Y vs\nB(x⊕y)B⁻¹", "half kick² vs\nx's 4-velocity", "association\nup to the twist",
          "interval vs\nMilne form", "for contrast:\nplain X·Y vs law"]
vals = [R["sandwich"]["sandwich"], R["sandwich"]["group_vs_law"], R["boost"]["square"], R["wigner"]["gyration"],
        R["milne"]["interval"], R["sandwich"]["plain_vs_law"]]
a.bar(range(6), [max(v, 1e-17) for v in vals], color=[BL] * 5 + [GR])
for i, v in enumerate(vals):
    a.text(i, max(v, 1e-17) * 3, f"{v:.0e}" if v < 0.01 else f"{v:.2f}", ha="center", fontsize=8.5)
a.set_yscale("log"); a.set_ylim(1e-17, 30)
a.set_xticks(range(6)); a.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
a.set_ylabel("largest relative mismatch over 2,000 random cases")
a.set_title("A. At q = −1 the law is the symmetric product of\nboosts, B·Y·B: relativistic velocity addition")

b = ax[1]
W = R["wigner"]["right_angle_pairs"]
xx = np.linspace(0, 1.6, 20); b.plot(xx, xx, color=GR, ls="--", lw=1, label="equal")
b.plot([w["wigner"] for w in W], [w["law_turn"] for w in W], "o", color=BL, ms=10, mfc="none", label="the law's twist")
for w in W:
    b.annotate(f"γ = {w['gamma1']:.2f}, {w['gamma2']:.2f}", (w["wigner"], w["law_turn"]), textcoords="offset points", xytext=(8, -10), fontsize=8, color=BL)
b.set_xlabel("textbook Wigner angle, cos θ = (γ₁ + γ₂)/(1 + γ₁γ₂)  (radians)")
b.set_ylabel("angle of the law's twist U (radians)")
b.legend(fontsize=8.5, loc="upper left")
b.set_title("B. The law's twist at q = −1 is the\nThomas–Wigner rotation (boosts at right angles)")

c = ax[2]
S = R["signature"]
c.loglog(S["small_q"], S["deviation"], "o-", color=BL, ms=8, label="law minus plain vector sum")
c.loglog(S["small_q"], np.array(S["small_q"]) * S["deviation"][0] / S["small_q"][0], "--", color=GR, label="in proportion to q")
c.set_xlabel("size of the dial |q|"); c.set_ylabel("largest difference, 2,000 pairs")
c.legend(fontsize=8.5, loc="upper left")
c.text(0.40, 0.27, "q > 0: the law's moves keep w² + x² + y² + z²\n         (4D rotations, Spin(4); rounds 12–18)\n"
       "q = 0: plain addition (Galilean)\n"
       "q < 0: the law's moves keep t² − x² − y² − z²\n         (the Lorentz group, Spin(1,3))\n"
       f"each breaks the other's form by ≥ {S['cross']:.2f}", transform=c.transAxes, fontsize=8.5, color=BL, va="top")
c.set_title("C. The dial's sign picks the group:\nrotations at q > 0, Lorentz boosts at q < 0")

d = ax[3]
row = R["row"]; t0, x0 = np.array(row["t"]), np.array(row["x"])
ll = np.linspace(0, 3.6, 2); d.plot(ll, ll, color=OR, lw=1.2); d.plot(-ll, ll, color=OR, lw=1.2, label="light cone")
for tau, al in ((1, 1.0), (2, 0.75), (3, 0.5)):
    rr = np.linspace(-1.9, 1.9, 200)
    d.plot(tau * np.sinh(rr), tau * np.cosh(rr), color=GR, lw=0.9)
    d.plot(tau * x0, tau * t0, "o", color=BL, ms=5, alpha=al)
    d.text(0.08, tau - 0.25, f"τ = {tau}", color=GR, fontsize=8.5)
for xi, ti in zip(x0, t0):
    d.plot([0, 3.4 * xi], [0, 3.4 * ti], color=BL, lw=0.5, alpha=0.5)
d.set_xlim(-3.6, 3.6); d.set_ylim(0, 3.6)
d.set_xlabel("space x"); d.set_ylabel("time t")
d.legend(fontsize=8.5, loc="upper left")
d.text(0.0, -0.17, f"blue dots: the law's kicks at q = −1, one per {row['step']:.3f} of rapidity, on slices τ = 1, 2, 3;\n"
       f"lines: their worldlines. Interval = τ₁² + τ₂² − 2τ₁τ₂ cosh d, d the law's distance (mismatch {R['milne']['interval']:.0e})",
       transform=d.transAxes, fontsize=8.5, color=BL)
d.set_title("D. Flat spacetime is the cone over the law's\nq = −1 space (the Milne form)")

fig.tight_layout()
fig.savefig("figures/f26_lorentz_from_the_law.png", bbox_inches="tight")
