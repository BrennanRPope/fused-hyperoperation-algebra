"""DRAFT round 36 picture: f43_charges.png, from results.json."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE, GREEN = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8", "#2e9e4f"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.18)
fig.suptitle("Round 36 (DRAFT): where hypercharge and colour could come from", fontsize=16)

# A: the law as a rotation
ax = fig.add_subplot(gs[0, 0])
t = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(t), np.sin(t), color=GREY, lw=1)
th = 1.0
ax.annotate("", xy=(np.cos(th), np.sin(th)), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.5))
ax.annotate("", xy=(1, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=GREY, lw=2.5))
ax.text(1.04, -0.08, "1", fontsize=15, color=GREY); ax.text(np.cos(th) + 0.03, np.sin(th) + 0.03, "x", fontsize=15, color=BLUE)
yv = np.array([np.cos(2.4), np.sin(2.4)]); R2 = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]); yr = R2 @ yv
ax.annotate("", xy=yv, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.5))
ax.annotate("", xy=yr, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=PURPLE, lw=2.5))
ax.text(yv[0] - 0.12, yv[1] + 0.05, "y", fontsize=15, color=ORANGE); ax.text(yr[0] - 0.33, yr[1] + 0.02, "x ⊕ y", fontsize=15, color=PURPLE)
aa = np.linspace(0, th, 30); ax.plot(0.3 * np.cos(aa), 0.3 * np.sin(aa), color=BLUE); ax.text(0.33, 0.12, "turn by the angle\nfrom 1 to x", color=BLUE, fontsize=11)
aa = np.linspace(2.4, 2.4 + th, 30); ax.plot(0.55 * np.cos(aa), 0.55 * np.sin(aa), color=PURPLE)
ax.text(-1.55, -1.05, "every direction across this plane is left alone", color=GREY, fontsize=11)
C1 = R["C1"]
txt = "\n".join(f"{k}: sandwich x^½ y x^½ = this turn to {v['rotation_vs_sandwich']:.0e}" for k, v in C1.items() if k != "sedenions")
txt += f"\nsedenions: sandwich falls off the sphere by up to {C1['sedenions']['sandwich_off_sphere']:.2f};\nthe turn still works, on any sphere"
ax.text(-1.55, 1.05, txt, fontsize=10.5, color=GREEN, va="bottom")
ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.2, 1.75); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("A. The law needs no multiplication: it is a turn in the plane of 1 and x (C1)")

# B: symmetry size against sphere size
ax = fig.add_subplot(gs[0, 1])
rows = R["C2"]["rows"]; ss = np.arange(2, 12)
ax.plot(ss, ss * (ss - 1) / 2, color=GREY, lw=1, ls="--", label="all rotations keeping 1 fixed: SO(sphere dimension)")
ax.plot([r["sphere"] for r in rows], [r["dim"] for r in rows], "o", color=BLUE, ms=11, label="symmetries of the law, found numerically")
ax.plot([7], [R["C2"]["g2"]], "s", color=ORANGE, ms=11, label="symmetries of octonion multiplication (G2)")
ax.annotate("3-sphere: SO(3) (round 1)", xy=(3, 3), xytext=(3.3, 14), color=BLUE, fontsize=11, arrowprops=dict(arrowstyle="->", color=BLUE))
ax.annotate("7-sphere: SO(7); the law\nignores the octonions", xy=(7, 21), xytext=(4.6, 30), color=ORANGE, fontsize=11, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.annotate("10-sphere: SO(10), the textbook\nunified group (45)", xy=(10, 45), xytext=(7.9, 33), color=GREEN, fontsize=11, arrowprops=dict(arrowstyle="->", color=GREEN))
ax.set_xlabel("dimension of the sphere the law lives on"); ax.set_ylabel("number of independent symmetries")
ax.set_title("B. The law's symmetries are every rotation keeping 1 fixed (C2)\nmirrors too, so it still picks no hand")
ax.legend(fontsize=10.5, loc="upper left")

# C: one generation from SO(10)
ax = fig.add_subplot(gs[1, 0])
st = R["C3"]["states"]; groups = {}
for s in st:
    key = (s["name"], s["Y"], s["Q"]); groups.setdefault(key, []).append(s)
col = {"triplet": BLUE, "antitriplet": ORANGE, "singlet": GREEN}
for (nm, Y, Q), ss_ in groups.items():
    y, q = float(Fr(Y)), float(Fr(Q)); c = col[ss_[0]["colour"]]
    ax.scatter(y, q, s=180 * len(ss_), color=c, alpha=0.8, edgecolor="k")
    lab = {"Q": "Q (u)" if q > 0 else "Q (d)", "L": "L (ν)" if q == 0 else "L (e)"}.get(nm, nm)
    ax.text(y + 0.06, q + 0.06, f"{lab} ×{len(ss_)}", fontsize=11, color=c)
ax.axhline(0, color=GREY, lw=0.5); ax.axvline(0, color=GREY, lw=0.5)
ax.set_xlabel("hypercharge Y = T3R + (B − L)/2   (fixed by SO(10), nothing to choose)"); ax.set_ylabel("electric charge Q = T3L + Y")
tr = R["C3"]["traces"]
ax.text(-1.25, 0.82, f"anomalies: Tr Y = {tr['Y']}, Tr Y³ = {tr['Y3']},\nTr T3² Y = {tr['T2Y']}, Tr T_colour² Y = {tr['Tc2Y']};\n{R['C3']['doublets']} weak doublets (even)", fontsize=11, color=PURPLE)
ax.scatter([], [], s=150, color=BLUE, label="colour triplet"); ax.scatter([], [], s=150, color=ORANGE, label="colour antitriplet"); ax.scatter([], [], s=150, color=GREEN, label="colour singlet")
ax.legend(fontsize=10.5, loc="lower right")
ax.set_xlim(-1.3, 1.4); ax.set_ylim(-1.3, 1.3)
ax.set_title("C. One hand of SO(10)'s spinor is exactly one generation (C3)\nthe 16 states with the measured charges (circle size = number of states)")

# D: running couplings
ax = fig.add_subplot(gs[1, 1])
g = np.array(R["C4"]["grid"]); run = np.array(R["C4"]["run"])
for i, (lab, c) in enumerate((("hypercharge (unified normalisation)", BLUE), ("weak SU(2)", ORANGE), ("colour SU(3)", GREEN))):
    ax.semilogx(g, run[:, i], color=c, lw=2.5, label=lab)
for k, c in (("meet12", BLUE), ("meet13", PURPLE), ("meet23", GREEN)):
    ax.axvline(R["C4"][k], color=c, ls=":", lw=1)
ax.axvline(1.22e19, color=GREY, ls="--", lw=1); ax.text(1.4e18, 55, "Planck\nmass", color=GREY, fontsize=10.5)
ax.text(3e6, 11, f"SO(10) fixes sin²θ_W = {R['C4']['sin2_unif']} where it holds;\nforcing one meeting gives {R['C4']['s2_pred']:.3f} at the Z mass,\nmeasured {R['C4']['s2_meas']}", color=RED, fontsize=11)
ax.set_xlabel("energy (GeV)"); ax.set_ylabel("1 / coupling strength")
ax.set_title("D. With only the Standard Model's particles, the three couplings never meet (C4)\nSO(10) needs a stage in between, which the law doesn't supply yet")
ax.legend(fontsize=10.5, loc="lower right", bbox_to_anchor=(1.0, 0.18))
fig.savefig(os.path.join(HERE, "figures", "f43_charges.png"), dpi=110, bbox_inches="tight")
