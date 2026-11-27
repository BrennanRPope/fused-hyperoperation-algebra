"""DRAFT picture for round 28. Run after checks.py. Blue = closed (the law's 3-spheres), green = flat, orange = open, grey = reference."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BL, GN, OR, GR, RD = "#1f6fc2", "#2e8b57", "#e07b00", "#777777", "#c0392b"
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
R = json.load(open("results.json"))
fig = plt.figure(figsize=(20, 13))
gs = fig.add_gridspec(2, 3)
a = fig.add_subplot(gs[0, 0], projection="3d")
H = 1.0
# hyperboloid -X0^2 + X1^2 + X2^2 = 1
T, P = np.meshgrid(np.linspace(-1.5, 1.5, 40), np.linspace(0, 2 * np.pi, 60))
a.plot_surface(np.cosh(T) * np.cos(P), np.cosh(T) * np.sin(P), np.sinh(T), color=GR, alpha=0.12, linewidth=0)
ph = np.linspace(0, 2 * np.pi, 200)
for t in (-1.0, 0.0, 1.0):   # closed: circles X0 = sinh t
    a.plot(np.cosh(t) * np.cos(ph), np.cosh(t) * np.sin(ph), np.sinh(t) + 0 * ph, color=BL, lw=2.2 if t == 0 else 1.4)
r = np.linspace(-1.6, 1.6, 200)
for t in (-0.6, 0.2, 0.8):   # flat: X0 = sinh t + r^2 e^t/2, X1 = cosh t - r^2 e^t/2, X2 = e^t r
    X0 = np.sinh(t) + r**2 * np.exp(t) / 2; X1 = np.cosh(t) - r**2 * np.exp(t) / 2; X2 = np.exp(t) * r
    m = np.abs(X0) < 2.2
    a.plot(X1[m], X2[m], X0[m], color=GN, lw=1.6)
ch = np.linspace(-1.3, 1.3, 200)
for t in (0.5, 0.9, 1.3):    # open: X0 = sinh t cosh chi, X1 = sinh t sinh chi, X2 = cosh t
    X0 = np.sinh(t) * np.cosh(ch); X1 = np.sinh(t) * np.sinh(ch); X2 = np.cosh(t) + 0 * ch
    m = np.abs(X0) < 2.2
    a.plot(X1[m], X2[m], X0[m], color=OR, lw=1.6)
a.set_xlabel("X1"); a.set_ylabel("X2"); a.set_zlabel("X0 (time)"); a.view_init(elev=16, azim=-58)
a.text2D(0.0, 0.02, "blue: closed slices (circles; the law's 3-spheres)\ngreen: flat slices\norange: open slices", transform=a.transAxes, fontsize=9.5)
a.set_title("A. One de Sitter spacetime, cut three ways (S1)")

b = fig.add_subplot(gs[0, 1])
t = np.linspace(-2.2, 2.2, 400)
b.plot(t, np.cosh(t), color=BL, lw=2.4, label="closed, k = +1: a = cosh(Ht)/H (stands still once, at t = 0)")
b.plot(t, np.exp(t), color=GN, lw=2.4, label="flat, k = 0: a = e^(Ht)")
tp = t[t > 0.02]; b.plot(tp, np.sinh(tp), color=OR, lw=2.4, label="open, k = −1: a = sinh(Ht)/H")
b.set_ylim(0, 6); b.set_xlabel("time t (units 1/H)"); b.set_ylabel("size of space a(t)"); b.legend(fontsize=9, loc="upper left")
b.text(0.6, 0.35, "all three: Λ = 3H², no matter,\nsame spacetime curvature everywhere", fontsize=9.5, color=GR)
b.set_title("B. The same Λ, three sizes of space (S1)")

c = fig.add_subplot(gs[0, 2])
c.plot(t, 1 / np.cosh(t)**2, color=BL, lw=2.4, label="closed: q(t) = H²/cosh²(Ht) > 0")
c.plot(t, 0 * t, color=GN, lw=2.4, label="flat: q = 0")
c.plot(tp, -1 / np.sinh(tp)**2, color=OR, lw=2.4, label="open: q(t) = −H²/sinh²(Ht) < 0")
c.axhline(1, color=RD, lw=1.2, ls="--", label="Λ/3 = H²")
c.plot([0], [1], "o", color=RD, ms=9)
c.annotate("the only place q = Λ/3:\nthe closed slice at its moment of rest", (0, 1), (0.5, 1.35), fontsize=9.5, color=RD, arrowprops=dict(arrowstyle="->", color=RD))
c.set_ylim(-2, 1.8); c.set_xlabel("time t (units 1/H)"); c.set_ylabel("the law's dial on the slice, q = k/a²"); c.legend(fontsize=9, loc="lower left")
c.set_title("C. The slice's dial can be any sign, and drifts (S2)")

d = fig.add_subplot(gs[1, 0:2]); d.axis("off")
boxes = [
 (0.5, 0.9, "Λ > 0 (measured)", GR),
 (0.5, 0.7, "the law: Λ = 3q and slice curvature = q\nholds only at a moment of rest (rounds 16, 20, 25)", BL),
 (0.5, 0.5, "our universe has no moment of rest (round 27, G2)", RD),
 (0.5, 0.3, "with matter: k = a₀²(8πGρ₀/3 + Λ/3 − H₀²)\nany sign, set by today's data (S3)", GR),
 (0.5, 0.1, "closed space only if one adds a premise:\n'matter's space composes by the q > 0 law'", OR),
]
for x, y, s, col in boxes:
    d.text(x, y, s, ha="center", va="center", fontsize=10.5, color=col, bbox=dict(boxstyle="round,pad=0.6", fc="white", ec=col, lw=1.8))
for y in (0.86, 0.66, 0.46, 0.26):
    d.annotate("", (0.5, y - 0.07), (0.5, y), arrowprops=dict(arrowstyle="->", color=GR, lw=1.6))
d.set_xlim(0, 1); d.set_ylim(0, 1)
d.set_title("D. Why the law doesn't force closed space")

e = fig.add_subplot(gs[1, 2])
S4 = R["S4"]
vals = [("the law (f = 1): torsion exactly 0", 1e-60, BL), ("flat twisted slices (f = 0 or 2)\non today's cosmic 3-sphere", S4["cosmic_flat_twist_GeV"], OR),
        ("spin-polarized torsion bounds\n(Kostelecký, Russell, Tasson 2008)", S4["bound_GeV"], RD)]
for i, (nm, v, col) in enumerate(vals):
    if i > 0: e.barh(i, np.log10(v) + 61, left=-61, height=0.38, color=col, alpha=0.75)
    e.text(-60.5, i + 0.24, nm, fontsize=9.5, color=col, va="bottom")
e.text(-60.5, -0.05, "exactly zero: no bar on a log scale", fontsize=9, color=BL, va="center")
e.text(np.log10(S4["cosmic_flat_twist_GeV"]) + 0.5, 1, f"{S4['cosmic_flat_twist_GeV']:.1e} GeV", va="center", fontsize=9, color=OR)
e.text(np.log10(S4["bound_GeV"]) + 0.5, 2, "~1e-31 GeV", va="center", fontsize=9, color=RD)
e.set_xlim(-61, -20); e.set_yticks([]); e.set_ylim(-0.4, 3.1)
e.set_xlabel("torsion, log10(GeV)")
e.text(-60.5, 2.9, "torsion of the family ∇_X Y = f[X,Y]/2 is 2(f − 1)", fontsize=9.5, color=GR)
e.set_title("E. Audit C4: the law has no torsion to couple to spin (S4)")

fig.suptitle("Round 28 (DRAFT): does the law force space to be closed, or fix a spin–twist coupling?", fontsize=14, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig("figures/f35_slicing.png")
print("wrote figures/f35_slicing.png")
