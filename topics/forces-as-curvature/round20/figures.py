"""DRAFT picture for round 20. Run after checks.py. Blue = the law, orange = light, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
ETA = np.array([-1.0, 1, 1, 1])
def oplus(a, b, qq=1.0):          # the same law as checks.py, Lorentzian dot
    d = lambda u, v: (u * v * ETA).sum(-1, keepdims=True)
    ab, aa, bb = d(a, b), d(a, a), d(b, b)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig = plt.figure(figsize=(23, 6.2))
ax = [fig.add_subplot(1, 4, 1, projection="3d"), fig.add_subplot(1, 4, 2), fig.add_subplot(1, 4, 3), fig.add_subplot(1, 4, 4)]

a = ax[0]
uu = np.linspace(-1.6, 1.6, 40); vv = np.linspace(0, 2 * np.pi, 60); U_, V_ = np.meshgrid(uu, vv)
a.plot_surface(np.cosh(U_) * np.cos(V_), np.cosh(U_) * np.sin(V_), np.sinh(U_), color=GR, alpha=0.12, linewidth=0)
G = R["grid"]
for i in sorted(set(g["i"] for g in G)):
    row = sorted([g for g in G if g["i"] == i], key=lambda g: g["j"])
    a.plot([g["W"] for g in row], [g["X"] for g in row], [g["T"] for g in row], color=BL, lw=1.0)
for j in sorted(set(g["j"] for g in G)):
    col = sorted([g for g in G if g["j"] == j], key=lambda g: g["i"])
    a.plot([g["W"] for g in col], [g["X"] for g in col], [g["T"] for g in col], color=BL, lw=0.6, alpha=0.6)
s = np.linspace(-1.5, 1.5, 2)
for sg in (1, -1):
    a.plot(np.ones(2), sg * s, s, color=OR, lw=2.2)
a.text(1, 1.5, 1.6, "light rays: straight lines", color=OR, fontsize=8.5)
a.set_xlabel("W"); a.set_ylabel("X (space)"); a.set_zlabel("T (time)")
a.view_init(elev=14, azim=-35); a.set_box_aspect((1, 1, 1))
a.set_title("A. 2D de Sitter grown by the law: rows are\nspace kicks (circles), columns time kicks")

b = ax[1]
G = R["grid"]
ch = {(g["i"], g["j"]): (g["T"] / (1 + g["W"]), g["X"] / (1 + g["W"])) for g in G if abs(g["j"]) <= 3}
def brk(seq):                     # break a line where it leaves through the chart's edge (a jump through infinity)
    out = [seq[0]]
    for p in seq[1:]:
        if abs(p[0] - out[-1][0]) + abs(p[1] - out[-1][1]) > 2.5 or max(abs(p[0]), abs(p[1])) > 6:
            out.append((np.nan, np.nan))
        out.append(p)
    return out
for i in sorted(set(k[0] for k in ch)):
    row = brk([ch[(i, j)] for j in sorted(k[1] for k in ch if k[0] == i)])
    b.plot([p[1] for p in row], [p[0] for p in row], "-", color=BL, lw=1.0 if i else 2.0)
for j in sorted(set(k[1] for k in ch)):
    col = brk([ch[(i, j)] for i in sorted(k[0] for k in ch if k[1] == j)])
    b.plot([p[1] for p in col], [p[0] for p in col], "-", color=BL, lw=0.6, alpha=0.6)
b.plot([ch[(0, 0)][1]], [ch[(0, 0)][0]], "o", color="k", ms=6)
xx = np.linspace(-3, 3, 400)
for sg in (1, -1):
    b.plot(xx, sg * np.sqrt(1 + xx**2), color=GR, lw=1.5, ls="--")
    b.plot(xx, sg * xx, color=OR, lw=2)
b.text(-0.95, -1.85, "dashed: chart edge t² − x² = 1,\nthe law's points at infinity", color=GR, fontsize=8.5, bbox=dict(fc="w", ec="none", alpha=0.85))
b.text(1.2, 0.55, "light:\nstraight lines", color=OR, fontsize=8.5, bbox=dict(fc="w", ec="none", alpha=0.85))
b.text(-1.95, 1.62, "thick row: space kicks from rest;\nthin lines: time kicks (q = 1)", color=BL, fontsize=8.5, bbox=dict(fc="w", ec="none", alpha=0.85))
b.set_xlim(-2, 2); b.set_ylim(-2, 2); b.set_aspect("equal")
b.set_xlabel("space x (the law's chart)"); b.set_ylabel("time t (the law's chart)")
b.set_title("B. The same grid in the law's own chart: space at +q,\ntime at −q (columns), light untouched (45°)")

c = ax[2]
rows = R["slices"]
tt = np.linspace(0, 1.7, 100)
c.plot(tt, np.cosh(tt), color=GR, ls="--", label="cosh(√q·t)/√q (rounds 14–16's slabs, q = 1)")
c.plot([r["tau"] for r in rows], [r["radius"] for r in rows], "o", color=BL, ms=10, mfc="none", label="space slice through the law's time kick")
c.set_xlabel("time from the bounce (the law's rapidity)"); c.set_ylabel("radius of the space slice")
c.legend(fontsize=8.5, loc="upper left")
E = R["einstein"]
c.text(0.04, 0.55, f"the law's metric: Ricci = 3q g (symbolic)\nEinstein with Λ = {E['Lambda'].replace('*', '')}, scalar curvature {E['scalar'].replace('*', '')}\n"
       "q > 0 de Sitter, q < 0 anti-de Sitter", transform=c.transAxes, fontsize=9, color=BL, va="top")
c.set_title("C. The law's spacetime solves Einstein with\nΛ = 3q and expands like round 14's slabs")

d = ax[3]
labels = ["on de Sitter", "adding x is linear", "keeps the 5D form", "twist keeps\n−t² + x² + y² + z²", "space = round 12", "time = law at −q",
          "light = plain sum", "anti-de Sitter (q < 0)", "slices"]
vals = [R["isometry"]["on_ds"], R["isometry"]["linear"], R["isometry"]["form"], R["twist"]["form"], R["directions"]["space"],
        R["directions"]["time"], R["directions"]["null"], R["sign"]["ads"], max(r["chart_round_trip"] for r in rows)]
d.bar(range(9), [max(v, 1e-17) for v in vals], color=BL)
for i, v in enumerate(vals):
    d.text(i, max(v, 1e-17) * 3, f"{v:.0e}" if v > 0 else "0", ha="center", fontsize=8)
d.set_yscale("log"); d.set_ylim(1e-17, 1)
d.set_xticks(range(9)); d.set_xticklabels(labels, fontsize=8, rotation=35, ha="right")
d.set_ylabel("largest mismatch (random cases)")
d.text(0.5, 0.88, f"the twist mixes time with space: largest entry {R['twist']['boost_part']:.2f}", transform=d.transAxes, ha="center", fontsize=8.5, color=BL)
d.set_title("D. The exact checks (D1–D6)\nall sit at rounding")

fig.tight_layout()
fig.savefig("figures/f27_de_sitter_from_the_law.png", bbox_inches="tight")
