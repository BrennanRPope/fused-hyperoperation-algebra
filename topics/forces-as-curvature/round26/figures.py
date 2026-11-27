"""DRAFT picture for round 26. Run after checks.py. Blue = the law and its gravity, orange = the weak force's twist, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR, PU = "#e07b00", "#1f6fc2", "#777777", "#7a3fb0"
R = json.load(open("results.json")); FG = R["figure"]
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()

# A. the law's lowest gravitational wave: frozen, grows through the waist, freezes again
a = ax[0]
w = FG["w1"]; e = np.array(w["eta"]); H = np.array(w["H"])
a.plot(e, H, color=BL, lw=2.6, label="squash H = h₁ − h₂ of the law's 3-sphere slices (n = 3 wave)")
a.set_yscale("symlog", linthresh=1.0)
a.axhline(0, color=GR, lw=0.6)
a.axvspan(-math.pi / 2, -1.2, color=GR, alpha=0.08); a.axvspan(1.2, math.pi / 2, color=GR, alpha=0.08)
i = int(np.argmax(np.abs(H)))
a.annotate(f"largest squash {abs(H[i]):.0f}× the frozen one,\naround the waist (η = 0), where the slices are smallest", (e[i], H[i]), (0.12, -40), fontsize=9.5, color=BL,
           arrowprops=dict(arrowstyle="->", color=BL))
a.annotate(f"frozen again at {R['W1']['frozen_value']:+.4f}", (e[-1], H[-1]), (0.7, 0.3), fontsize=9.5, color=BL,
           arrowprops=dict(arrowstyle="->", color=BL))
a.text(-1.42, 1.6, "frozen in\nthe far past", fontsize=9, color=GR)
a.text(-1.0, 30, "H'' + 2ℋH' + 8H = 0 (conformal time η)\nlaw's de Sitter: a = sec η; all of history is η ∈ (−π/2, π/2)\nso it crosses zero only twice (η ≈ −1.49 and 0)", fontsize=9, color=GR)
a.text(1.22, -600, "outside the horizon:\nslices grow faster\nthan it can swing", fontsize=9, color=GR)
a.set_xlabel("conformal time η (q = 1): only π of it in all of history"); a.set_ylabel("squash H (linear near 0, logarithmic beyond ±1)")
a.set_ylim(-6000, 6000)
a.legend(fontsize=8.5, loc="lower right", bbox_to_anchor=(1.0, 0.0))
a.set_title("A. The law's lowest gravitational wave swings about once, then freezes (W1)")

# B. electric and magnetic tide of that wave
b = ax[1]
w = FG["w4"]; e = np.array(w["eta"]); E = np.array(w["E"]); B = np.array(w["B"])
b.semilogy(e, E, color=BL, lw=2.4, label="electric tide |E| (stretch and squeeze)")
b.semilogy(e, B, color=PU, lw=2.4, label="magnetic tide |B| = 3|ḣ|/a (needs motion)")
b.axvspan(-0.6, 0.6, color=GR, alpha=0.08)
b.text(-0.55, 2e-4, f"swinging: rms |B|/|E| = {R['W4']['ratio_swinging']:.2f}", fontsize=9.5, color=GR)
b.text(0.75, 3e-9, f"frozen: |B|/|E| falls to\n{R['W4']['ratio_end']:.0e} at the end of time;\nboth fade as the slices grow", fontsize=9.5, color=GR)
b.set_ylim(1e-12, 1e5); b.set_xlabel("conformal time η (q = 1)"); b.set_ylabel("size of each tide (physical units)")
b.legend(fontsize=8.5, loc="upper right")
b.set_title("B. A wave carries both halves of the tide; frozen, only the electric half is left (W2, W4)")

# C. a ripple in the twist makes a wave only on the law's field
c = ax[2]
w = FG["w3"]; e = np.array(w["eta"]); Hh = np.array(w["H"]); U = np.array(w["U"])
c.plot(e, U / U[0], color=OR, lw=2.2, label="twist ripple U = u₁ − u₂ (÷ its start 0.01), on the law's f = 1 field")
c2 = c.twinx()
c2.plot(e, Hh * 1e6, color=BL, lw=2.6, label="gravitational wave H it makes (× 10⁶), κ = 0.002")
c2.axhline(0, color=GR, lw=0.6)
c2.plot(e, 0 * e, color=GR, lw=1.6, ls="--", label="the same ripple on an empty weak vacuum, f = 0 or f = 2: H = 0 exactly")
c2.set_ylabel("H × 10⁶", color=BL); c2.set_ylim(-9.5, 1.6); c.set_ylim(-1.75, 1.15)
c.set_xlabel("conformal time η (closed universe, Λ = 3q, q = 0.3)"); c.set_ylabel("U / U₀", color=OR)
c.text(0.03, 0.05, "coupling of ripple to wave ∝ f(f − 2)(f + 2) and f(f − 2):\nzero at the two vacua, 3 and −6 at the law's f = 1\nwave size ∝ κ = 4πG/g² (ratios %.3f, %.3f)" % tuple(R["W3"]["ratios"]),
       transform=c.transAxes, fontsize=9, color=GR)
h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels()
c.legend(h1 + h2, l1 + l2, fontsize=8.5, loc="upper right")
c.set_title("C. A twist ripple makes a gravitational wave only on the law's field (W3)")

# D. a spinning lump drags frames; magnetic tide needs mass and spin
d = ax[3]
r = np.logspace(np.log10(0.6), np.log10(12), 300)
J, M = 0.1, 0.25
d.loglog(r, 2 * J / r**3, color=BL, lw=2.4, label="frame dragging w = 2J/r³ (J = 0.1), the same for every q")
d.loglog(r, 288 * J * M / r**7, color=PU, lw=2.4, label="magnetic tide invariant R·*R = 288 J M cos θ / r⁷ (θ = 0, M = 0.25)")
d.loglog(r, 48 * M**2 / r**6, color=GR, lw=1.6, ls="--", label="electric tide of the same lump at rest, 48 M²/r⁶ (round 23)")
d.text(0.04, 0.06, "(r⁴ w')' = 0 from Einstein's equations with Λ = 3q:\nq drops out of the drag and of R·*R\nno spin (J = 0) or no mass (M = 0): no magnetic tide",
       transform=d.transAxes, fontsize=9, color=GR)
d.set_xlabel("radius r"); d.set_ylabel("size")
d.legend(fontsize=8.5, loc="upper right")
d.set_title("D. A slowly spinning lump in the law's de Sitter (W2)")

fig.suptitle("Round 26 (DRAFT): the moving half of the tide, and gravitational waves on the law's 3-sphere", fontsize=14, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig("figures/f33_moving_tide_waves.png")
print("wrote figures/f33_moving_tide_waves.png")
