"""DRAFT figures for the exact split of 1 - 6M/r. Run after checks_six.py: python3 figures_six.py (numpy, matplotlib)."""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
KNOWN, OURS, GREY, RED = "#d9822b", "#2b6cb0", "#888888", "#c53030"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
res = json.load(open(os.path.join(HERE, "results_six.json")))
def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=140, bbox_inches="tight"); plt.close(fig); print("wrote", name)

# ---------------------------------------------------------------- O6: the split
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axs[0]
x = np.linspace(0, 1 / 3, 400)
v2 = x / (1 - 2 * x)
ax.plot(x, 1 - 2 * x, color=GREY, lw=2.5, label="space factor 1 − 2M/r (same for both)")
ax.plot(x, 1 - v2, color=KNOWN, lw=2, label="gyroscope's velocity factor 1 − v²")
ax.plot(x, np.where(4 * v2 <= 1.6, 1 - 4 * v2, np.nan), color=OURS, lw=2, label="orbit's velocity factor 1 − (2v)²")
ax.plot(x, 1 - 3 * x, color=KNOWN, lw=1.5, ls="--", label="gyroscope total 1 − 3M/r (product)")
ax.plot(x, 1 - 6 * x, color=OURS, lw=1.5, ls="--", label="orbit total 1 − 6M/r (product)")
geo = [g for g in res["geodesics"] if g[0] == "Schwarzschild"]
xs = [1 / g[1] for g in geo]; Ts = [1 / (1 + g[3] / (2 * math.pi)) ** 2 for g in geo]
ax.plot(xs, Ts, "o", color="black", ms=6, label="simulated orbits (geodesics)")
ax.axhline(0, color="black", lw=0.6)
ax.axvline(1 / 6, color=OURS, lw=0.8, ls=":"); ax.text(1 / 6 + 0.003, 0.82, "r = 6M: innermost stable orbit\norbit's 2v reaches c", color=OURS, fontsize=9)
ax.axvline(1 / 3, color=KNOWN, lw=0.8, ls=":"); ax.text(1 / 3 - 0.003, 0.55, "r = 3M: photon sphere\ngyroscope's v\nreaches c", color=KNOWN, fontsize=9, ha="right")
ax.set_xlim(0, 0.345); ax.set_ylim(-0.65, 1.05)
ax.set_xlabel("field strength M/r"); ax.set_ylabel("factor (radial ÷ angular frequency, squared)")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("A. Both splits are exact: the orbit and the gyroscope share\nthe space factor; the orbit's velocity factor runs at 2v")
ax = axs[1]
v = math.sqrt(1 / 10 / (1 - 2 / 10))          # r = 10M
th = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(th), np.sin(th), ":", color=GREY); ax.text(0, 1.04, "edge = light speed", ha="center", color=GREY, fontsize=9)
ax.fill(2 * v * np.cos(th), 2 * v * np.sin(th), color=OURS, alpha=0.18)
ax.plot(2 * v * np.cos(th), 2 * v * np.sin(th), color=OURS, lw=2)
ax.fill(v * np.cos(th), v * np.sin(th), color=KNOWN, alpha=0.35)
ax.plot(v * np.cos(th), v * np.sin(th), color=KNOWN, lw=2)
ax.annotate("", xy=(v, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=KNOWN))
ax.annotate("", xy=(-2 * v * 0.7071, 2 * v * 0.7071), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=OURS))
ax.text(v / 2, -0.09, "v", color=KNOWN, ha="center"); ax.text(-0.33, 0.42, "2v", color=OURS)
ax.text(-1.05, -1.22, f"at r = 10M, v = {v:.3f}c\norange: the gyroscope's disc (Thomas), radius v\nblue: the orbit's disc, radius 2v, four times the area",
        fontsize=9.5)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.35, 1.15)
ax.set_title("B. In relativity's velocity space: the orbit\nturns by the disc of radius 2v")
fig.suptitle("Figure O6. Einstein's 1 − 6M/r splits exactly: (1 − 2M/r) for curved space × (1 − (2v)²) for a velocity disc at twice the orbital speed", y=1.02)
save(fig, "figO6-the-six-splits.png")

# ---------------------------------------------------------------- O7: one rule
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5))
ax = axs[0]
rows = [("Newton (any force)", "T = N", "N = 3 + r·(force)′/(force)\n= 1 for 1/r², so orbits close", GREY),
        ("electric, special relativity", "T = N − v²", "exact, any potential:\nthe velocity disc of radius v", KNOWN),
        ("gravity, general relativity", "T = (N − 4v²)/B", "exact, any static round metric:\nspace factor 1/B × disc of radius 2v", OURS)]
for i, (name, law, note, col) in enumerate(rows):
    y = 2 - i
    ax.text(0.0, y, name, fontsize=11.5, fontweight="bold", va="center")
    ax.text(0.47, y, law, fontsize=15, color=col, va="center", fontweight="bold")
    ax.text(0.47, y - 0.36, note, fontsize=9, va="center", color="#333333")
ax.text(0, -0.85, "T = (radial frequency ÷ angular frequency)²;  advance per lap = 2π(1/√T − 1)\n"
        "v = local orbital speed (c = 1);  N is Newton's apsidal number (Bertrand)", fontsize=9, color=GREY)
ax.set_xlim(0, 1); ax.set_ylim(-1.2, 2.5); ax.axis("off")
ax.set_title("A. One rule, three forces (circular orbits, small wobbles)")
ax = axs[1]
pts = res["sr_points"]
ax.plot([p[0] for p in pts], [p[1] for p in pts], "s", color=KNOWN, label="electric type: 5 potentials × 2 radii")
cols = {"Schwarzschild": OURS, "clock warp only (B = 1)": "#6b46c1", "Reissner–Nordström Q = 0.6": "#2f855a"}
for name, col in cols.items():
    g = [q for q in res["geodesics"] if q[0] == name]
    Tp = [1 / (1 + q[4] / (2 * math.pi)) ** 2 for q in g]; Ts = [1 / (1 + q[3] / (2 * math.pi)) ** 2 for q in g]
    ax.plot(Tp, Ts, "o", color=col, label=f"gravity: {name}")
lim = [-1, 7.2]
ax.plot(lim, lim, color="black", lw=0.7)
ax.set_xlabel("T predicted by the rule"); ax.set_ylabel("T measured from the simulated orbit")
ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("B. Every simulated case lands on the line")
save(fig, "figO7-one-rule-three-forces.png")
