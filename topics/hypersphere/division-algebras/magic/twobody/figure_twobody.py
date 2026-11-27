"""DRAFT figure D5: the pair rule for two comparable black holes. Run after checks_twobody.py."""
import os, sys, itertools; sys.path.insert(0, '..')
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stu2 import two_body
from orbit1pn import stu_params, measure
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "..", "..", "figures")
KNOWN, OURS, GREY, GREEN = "#d9822b", "#2b6cb0", "#bbbbbb", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
fig, axs = plt.subplots(1, 4, figsize=(24, 6.6), gridspec_kw={"width_ratios": [1.15, 0.9, 0.9, 1.0]})
ax = axs[0]
lines = [("four charges vs one opposite charge", [1, 1, 1, 1], [1]*4, [4, 0, 0, 0], [-1]*4, KNOWN),
         ("four charges vs four opposite", [1, 1, 1, 1], [1]*4, [1, 1, 1, 1], [-1]*4, OURS),
         ("two charges vs two opposite", [1, 1, 0, 0], [1]*4, [1, 1, 0, 0], [-1]*4, GREEN),
         ("one charge vs one charge (any directions)", [1, 0, 0, 0], [1]*4, [1, 0, 0, 0], [-1]*4, "k")]
X = np.linspace(0, 1, 11)
for lab, kA, sA, kB, sB, col in lines:
    T, TA, TB, _ = two_body(kA, sA, kB, sB)
    ax.plot([0, 1], [float(TB), float(TA)], color=col, lw=2.6, label=f"{lab}: {float(TB):g} → {float(TA):g}")
ax.axvline(0.5, color=GREY, ls="--"); ax.text(0.51, 3.25, "equal masses", color="#666", fontsize=10)
ax.annotate("at equal masses each body\nbrings half its own pairs' turn:\n3 × ½ = 3/2", xy=(0.5, 1.5), xytext=(0.58, 0.42), fontsize=10, arrowprops=dict(arrowstyle="->"))
ax.annotate("Caron-Huot & Zahraee's half-BPS\npairs: no turning at any mass", xy=(0.2, 0.0), xytext=(0.03, 0.25), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
ax.set_xlim(0, 1); ax.set_ylim(-0.15, 3.5)
ax.set_xlabel("share of the total mass carried by the first body  (0 = it is the light visitor, 1 = it is the heavy one)")
ax.set_ylabel("turn per lap ÷ (π K²/L²)"); ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(0.0, 0.83))
ax.set_title("A. Two comparable black holes: the turn slides in a\nstraight line between the two 'small visitor' answers")
def paths(kA, sA, kB, sB, laps=16):
    T, TA, TB, _ = two_body(kA, sA, kB, sB); par = stu_params(kA, sA, kB, sB, TA, TB)
    _, out = measure(par, Jscale=6.0, ecc=0.3, laps=laps)
    r, ph = out.y[0], out.y[2]; XA = float(par["XA"])
    x, y = r*np.cos(ph), r*np.sin(ph)
    return (1 - XA)*x, (1 - XA)*y, -XA*x, -XA*y, T
for ax, (kA, sA, kB, sB), ttl in [(axs[1], ([1, 0, 0, 0], [1]*4, [1, 0, 0, 0], [-1]*4), "B. Equal masses, one charge each:\nthe orbits close"),
                                  (axs[2], ([1, 1, 1, 1], [1]*4, [4, 0, 0, 0], [-1]*4), "C. Equal masses, four charges vs one:\nthe orbits turn (T = 3/2)")]:
    xa, ya, xb, yb, T = paths(kA, sA, kB, sB)
    s = np.max(np.hypot(xa, ya))
    ax.plot(xa/s, ya/s, color=OURS, lw=1.5, label="first body"); ax.plot(xb/s, yb/s, color=KNOWN, lw=1.5, label="second body")
    ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(loc="lower center", fontsize=9.5, bbox_to_anchor=(0.5, -0.17), ncol=2); ax.set_title(ttl)
for a_ in axs[1:3]: a_.text(0, -1.08, "16 laps of the first-order two-body equations,\nstrong field so the turn is visible", ha="center", fontsize=8.8, color="#555")
ax = axs[3]; ax.axis("off"); ax.set_xlim(0, 4.2); ax.set_ylim(-1.2, 2.7)
def body(ox, oy, n, col, lab):
    pts = [(0, 0), (1, 0), (0.5, 0.87), (0.5, 0.3)][:n]; P = [(ox + 0.6*x, oy + 0.6*y) for x, y in pts]
    for i, j in itertools.combinations(range(n), 2): ax.plot([P[i][0], P[j][0]], [P[i][1], P[j][1]], color=col, lw=2.5)
    for p in P: ax.plot(*p, "o", color=col, ms=11)
    ax.text(ox + 0.3, oy - 0.3, lab, ha="center", va="top", fontsize=10)
body(0.3, 1.2, 4, OURS, "first body\n4 charges, 6 pairs\nhalf the mass")
body(2.6, 1.4, 1, KNOWN, "second body\n1 charge, no pairs\nhalf the mass")
ax.annotate("", xy=(2.45, 1.45), xytext=(1.05, 1.45), arrowprops=dict(arrowstyle="<->", color="#555"))
ax.text(1.75, 1.55, "pull", ha="center", fontsize=9.5, color="#555")
ax.text(0.05, -0.1, "Turn = (each body's mass share) × (the turn its own\npairs of charges give a small visitor), summed.\nA body with one charge adds nothing, so a pair of\nsingle-charge black holes never turns, at any masses.\nAny body with two or more charges turns the orbit.", fontsize=10, va="top")
ax.set_title("D. The rule for two bodies: each body's own\npairs, weighted by its share of the mass")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD5-pair-rule-for-two-comparable-black-holes.png"), dpi=130); print("ok")
