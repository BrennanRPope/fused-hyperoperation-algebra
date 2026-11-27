import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

F = np.load("frames.npz"); T = np.load("toy_tracks.npz"); M = np.load("merger.npz"); Dp = np.load("dispersion.npz")
d = json.load(open("dialwick.json")); P = json.load(open("kappa_pairs.json"))
blue, red, grey, green, purple = "#2b6cb0", "#c53030", "#718096", "#2f855a", "#6b46c1"
roots = F["roots"][:5]

# ------------------------------------------------------------------ Figure A: the quarter turn of the d dial
fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
a = ax[0]
a.axhline(0, color="k", lw=0.8); a.axvline(0, color="k", lw=0.5)
a.plot(roots, 0 * roots, "o", color=red, ms=8, label="zeros of H₀ = poles of the shock field")
a.plot(-roots, 0 * roots, "o", color=red, ms=8)
for i, r in enumerate(roots[:3]):
    a.annotate(f"{r:.2f}", (r, 0), (r - 4, 3 + 3.5 * i), fontsize=9, color=red, arrowprops=dict(arrowstyle="-", color=red, lw=0.6))
a.set_xlim(-72, 72); a.set_ylim(-20, 20)
a.set_xlabel("Re z"); a.set_ylabel("Im z")
a.set_title("A1  Before the turn (d = +1, ζ's own frame)\nshock field has poles ON the real line: unphysical, needs viscosity −1")
a.legend(loc="lower left", fontsize=9)
a.text(0, -15, "RH: every zero on this line", ha="center", fontsize=10)

a = ax[1]
im = a.pcolormesh(F["WX"], F["WY"], F["logG"], shading="auto", cmap="viridis")
a.plot(0 * roots, roots, "o", color="w", mec=red, ms=8, label="zeros of G₀ = poles of the shock field")
a.plot(0 * roots, -roots, "o", color="w", mec=red, ms=8)
a.axhline(0, color="w", lw=1.2); a.text(2.5, 2.5, "real line: no poles,\nφ > 0 here", color="w", fontsize=9)
a.set_xlabel("Re w  (w = −i z)"); a.set_ylabel("Im w")
a.set_title("A2  After the quarter turn x → i x (d → −d)\npoles stacked on ONE vertical line (the centreline)")
plt.colorbar(im, ax=a, label="log₁₀ |G₀(w)|")
a.legend(loc="lower left", fontsize=8)

a = ax[2]
xs = F["xs"]
for key, col, lab in [("um", purple, "t = −0.3 (rewound)"), ("u0", "k", "t = 0"), ("up", green, "t = +0.3")]:
    a.plot(xs, F[key], color=col, lw=1.8, label=lab)
a.axhline(0, color=grey, lw=0.6)
a.set_xlabel("x  (real line of the turned frame = the field h)"); a.set_ylabel("shock field u = −2 G′/G = −2 × mean spin")
a.set_title("A3  The shock field there is smooth for EVERY t\nφ = G_t(x) > 0 always (from Φ > 0), RH or not")
a.legend(fontsize=9, loc="upper right")
a.text(-58, -0.55, f"min G over t ∈ [−1, 1], |x| ≤ 60:  {d['physical_every_t']['min_G_over_grid']:.2f} > 0", fontsize=9)
a.text(-58, -0.68, "the three curves nearly overlap: t barely moves the field here", fontsize=9, color=grey)
plt.tight_layout(); plt.savefig("figures/figA_quarter_turn.png", dpi=130); plt.close()

# ------------------------------------------------------------------ Figure B: what RH becomes in the shock frame
fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
c = d["toy"]["c"]; Lam = d["toy"]["Lambda"]; t = T["t"]; z = T["z"]; w = T["w"]
a = ax[0]
before = t < Lam
a.plot(z.real[before], z.imag[before], color=red, lw=2, label="t < Λ: pair off the line")
a.plot(z.real[before], -z.imag[before], color=red, lw=2)
a.plot(z.real[~before], 0 * z.real[~before], color=blue, lw=2, label="t > Λ: both on the line")
a.plot(2 * math.pi - z.real[~before], 0 * z.real[~before], color=blue, lw=2)
a.plot([math.pi], [0], "k*", ms=13, label=f"landing at t = Λ = ln c = {Lam:.3f}")
a.axhline(0, color="k", lw=0.6)
a.set_xlabel("Re z"); a.set_ylabel("Im z"); a.set_xlim(1.2, 5.1)
a.set_title(f"B1  Toy H_t = c + eᵗ cos z  (c = {c})\nζ's frame: a zero pair lands on the line")
a.legend(fontsize=8, loc="upper left")

a = ax[1]
a.plot(w.real[before], w.imag[before], color=red, lw=2, label="t < Λ: two pole columns at x = ±a")
a.plot(-w.real[before], w.imag[before], color=red, lw=2)
a.plot(0 * w.real[~before], w.imag[~before], color=blue, lw=2, label="t > Λ: one column, poles spread up/down")
a.plot(0 * w.real[~before], -2 * math.pi - w.imag[~before], color=blue, lw=2)
a.plot([0], [-math.pi], "k*", ms=13, label="columns merge at t = Λ")
a.axvline(0, color=grey, ls="--", lw=0.8); a.text(0.05, -2.3, "centreline", color=grey, fontsize=9)
a.set_xlabel("Re w  (shock frame)"); a.set_ylabel("Im w")
a.set_title("B2  Same toy after the turn:\nRH = 'all pole columns have merged onto the centreline'")
a.legend(fontsize=8, loc="lower left")

a = ax[2]
u = np.linspace(0, 1.0, 400)
def Phi(u):
    s = 0
    for n in range(1, 10):
        s = s + (2 * np.pi**2 * n**4 * np.exp(9 * u) - 3 * np.pi * n**2 * np.exp(5 * u)) * np.exp(-np.pi * n**2 * np.exp(4 * u))
    return s
for tt, col, lab in [(-15.0, purple, "t = −15: antiferromagnetic tilt (rewind)"), (0.0, "k", "t = 0: ζ's spin, RH ⇔ Lee–Yang here"),
                     (15.0, green, "t = +15: ferromagnetic tilt (forward)")]:
    a.plot(u, np.exp(tt * u**2) * Phi(u), color=col, lw=1.8, label=lab)
a.set_xlabel("spin value u"); a.set_ylabel("weight e^{t u²} Φ(u)")
a.set_title("B3  The physical reading: a 'spin' with weights Φ(u)\nshock field = −2 × its magnetisation in field x")
a.legend(fontsize=8)
a.text(0.45, 0.55 * a.get_ylim()[1], "toy: spin 0 (weight c), ±1 (weight eᵗ/2)\nLee–Yang  ⇔  t ≥ ln c", fontsize=9)
plt.tight_layout(); plt.savefig("figures/figB_lee_yang_columns.png", dpi=130); plt.close()

# ------------------------------------------------------------------ Figure C: the other turns
fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
a = ax[0]
q = Dp["q"]; r = Dp["rate"]
cols = {-1.0: red, -0.5: "#dd6b20", 0.0: grey, 0.5: "#3182ce", 1.0: blue}
for kap, col in cols.items():
    lab = {-1.0: "metric d = −1 (ζ's law)", 1.0: "metric d = +1 (real charges)"}.get(kap, f"metric d = {kap:+.1f}")
    a.plot(q, kap * r, color=col, lw=1.8, label=lab)
a.axhline(0, color="k", lw=0.6)
a.set_xlabel("ripple wavenumber q"); a.set_ylabel("growth rate of a vertical ripple")
a.set_title("C1  Half turn (metric only): the line's pull\nrate = d × 2(πq − q²/2): attracts for d < 0, repels for d > 0")
a.legend(fontsize=8)

a = ax[1]
for kap, col in cols.items():
    tr = P[str(kap)]
    a.plot(tr["t"], tr["y"], color=col, lw=1.8, label=f"d = {kap:+.1f}")
a.set_xlabel("heat time t"); a.set_ylabel("height of an off-line pair in a crystal")
a.set_title("C2  Same pair, five settings: captured only for d < 0\n(flip exactly at d = 0; only d = −1 is a function's zeros)")
a.legend(fontsize=8)

a = ax[2]
th = np.linspace(0, 2 * np.pi, 400)
a.plot(np.cos(th), np.sin(th), color=grey, lw=1)
pts = [(-1, 0, red, "θ = 0: viscosity −1\nζ's frame, backward heat"),
       (0, 1, purple, "θ = π/4: viscosity i\nSchrödinger (twist page)"),
       (1, 0, blue, "θ = π/2: viscosity +1\nquarter-turn frame, physical")]
for x, y, col, lab in pts:
    a.plot(x, y, "o", color=col, ms=11)
    a.annotate(lab, (x, y), (-0.55 if x == 0 else (-1.5 if x < 0 else 0.35), 0.62 if x == 0 else -0.32), fontsize=9, color=col)
a.add_patch(FancyArrowPatch((-0.92, 0.38), (0.38, 0.92), connectionstyle="arc3,rad=-0.45", arrowstyle="->", mutation_scale=15, color="k"))
a.set_aspect("equal"); a.set_xlim(-1.6, 1.6); a.set_ylim(-1.1, 1.4)
a.set_xlabel("Re (viscosity)"); a.set_ylabel("Im (viscosity)")
a.set_title("C3  Viscosity phase = complex zoom by e^{iθ}\nviscosity(θ) = −e^{−2iθ}: one circle, the quarter turn is its far end")
plt.tight_layout(); plt.savefig("figures/figC_other_turns.png", dpi=130); plt.close()
print("ok")
