"""DRAFT. Labeled pictures for the wall-shapes study. Reads walls.json, depth.json, rank.json.
Writes figures/W1..W6 (png)."""
import json, math, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

S = ["#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"]
TXT, TXT2, GRID = "#191919", "#666666", "#e3e3e3"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": TXT2, "axes.labelcolor": TXT, "xtick.color": TXT2,
                     "ytick.color": TXT2, "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False})
walls = json.load(open("walls.json")); depth = json.load(open("depth.json")); rank = json.load(open("rank.json"))
SH = ["I", "II", "III", "IV"]
SHN = {"I": "I  turn by i", "II": "II  many scales", "III": "III  depth", "IV": "IV  no daylight"}
BR = ["positivity", "growth", "number theory", "topology"]

# ---------- W1 board ----------
fig, ax = plt.subplots(figsize=(10, 9.5))
order = sorted(range(len(walls)), key=lambda i: (BR.index(walls[i]["break_type"]), walls[i]["name"]))
for row, i in enumerate(order):
    w = walls[i]; y = len(order) - 1 - row
    for j, s in enumerate(SH):
        if s in w["shapes"]:
            ax.scatter(j, y, s=140, color=S[j], zorder=3, edgecolor="white", linewidth=1.5)
    ax.text(4.0, y, w["break_type"], va="center", color=S[4 + BR.index(w["break_type"])] if False else TXT2, fontsize=9)
ax.set_yticks(range(len(order))); ax.set_yticklabels([walls[i]["name"] for i in order][::-1], fontsize=8.5)
ax.set_xticks(range(4)); ax.set_xticklabels([SHN[s] for s in SH], rotation=20, ha="right")
ax.text(4.0, len(order) - 0.2, "kind of break", fontweight="bold", color=TXT, fontsize=9)
ax.set_xlim(-0.6, 5.3); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title("W1. Every wall, its shape(s) and its kind of break (25 walls)", loc="left")
fig.tight_layout(); fig.savefig("figures/W1_wall_board.png", dpi=150); plt.close(fig)

# ---------- W2 shape x break counts ----------
fig, ax = plt.subplots(figsize=(8, 4))
left = np.zeros(4)
for b_i, b in enumerate(BR):
    vals = np.array([sum(1 for w in walls if s in w["shapes"] and w["break_type"] == b) for s in SH])
    bars = ax.barh(range(4), vals, left=left, color=S[b_i], edgecolor="white", linewidth=2, label=b)
    for k, v in enumerate(vals):
        if v: ax.text(left[k] + v / 2, k, str(v), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    left += vals
ax.set_yticks(range(4)); ax.set_yticklabels([SHN[s] for s in SH]); ax.invert_yaxis()
ax.set_xlabel("number of walls (a wall can carry two shapes)")
ax.legend(frameon=False, loc="lower right", title="kind of break")
ax.set_title("W2. Shape I (turn by i) is almost all positivity breaks; III (depth) mixes kinds", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("figures/W2_shapes_by_break.png", dpi=150); plt.close(fig)

# ---------- W3 phase budget ----------
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2))
titles = ["depth 1: one factor", "depth 2: two factors", "depth 3+: three or more", "sector law"]
for n, ax in zip([1, 2, 3], axs[:3]):
    ax.set_aspect("equal"); ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.2)
    reach = min(n * 180, 360)
    ax.add_patch(Wedge((0, 0), 2.1, 0, reach if reach < 360 else 359.999, color=S[0], alpha=0.18, lw=0))
    ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
    if n == 2: ax.plot([0, 2.1], [0, 0], color=S[0], lw=2.5, ls=(0, (1, 1.5)))   # missed ray
    ax.scatter([1.3], [0], s=90, color=S[5], zorder=4, edgecolor="white")
    ax.scatter([-1.3], [0], s=90, color=S[7], zorder=4, edgecolor="white")
    ax.annotate("W > 0", (1.3, 0), (1.15, 0.45), color=TXT, fontsize=9)
    ax.annotate("W < 0", (-1.3, 0), (-1.75, 0.45), color=TXT, fontsize=9)
    ax.scatter([0], [0], s=30, color=TXT, zorder=5)
    verdict = {1: "both missed: stable", 2: "W > 0 missed: stable iff W >= 0", 3: "both hit: stable only at W = 0"}[n]
    ax.set_title(f"{titles[n-1]}\n{verdict}", fontsize=9.5)
    ax.text(0, -2.05, "shaded = reachable values of u1...un,\nevery ui in the upper half-plane", ha="center", fontsize=7.5, color=TXT2)
    ax.set_xticks([]); ax.set_yticks([])
ax = axs[3]
al = np.linspace(0.05, 1.0, 200)
for n, c in [(1, S[0]), (2, S[1]), (3, S[2]), (4, S[6])]:
    pass
nn = np.linspace(1, 8, 300)
ax.plot(nn, 2 / nn, color=S[5], lw=2, label="W > 0 safe below:  n alpha = 2 pi")
ax.plot(nn, 1 / nn, color=S[7], lw=2, label="W < 0 safe below:  n alpha = pi")
ax.scatter([2], [1], color=TXT, zorder=5); ax.annotate("half-plane, depth 2", (2, 1), (2.45, 1.02), fontsize=8.5)
ax.set_xlabel("depth n (number of factors)"); ax.set_ylabel("sector opening alpha / pi")
ax.set_ylim(0, 1.15); ax.legend(frameon=False, fontsize=8, loc="center right"); ax.set_title(titles[3] + "\ninputs with arg in (0, alpha)", fontsize=9.5)
ax.grid(color=GRID, lw=0.6)
fig.suptitle("W3. Phase budget: arguments add under composition, so the turn-by-i positivity survives only to depth 2 (u = 1 + k a, W = 1 + k c)",
             x=0.01, ha="left", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig("figures/W3_phase_budget.png", dpi=150); plt.close(fig)

# ---------- W4 Lee-Yang zeros ----------
def zeros(N, terms):
    coef = np.zeros(N + 1)
    for s in itertools.product([1, -1], repeat=N):
        E = sum(J * np.prod([s[i] for i in A]) for J, A in terms); coef[s.count(-1)] += math.exp(E)
    return np.roots(coef[::-1])
P4 = list(itertools.combinations(range(4), 2))
cases = [("pairs only\nJ2 = 0.5", [(0.5, p) for p in P4]),
         ("one 4-body term only\nJ4 = 0.5 (parity)", [(0.5, (0, 1, 2, 3))]),
         ("pairs + 4-body\nJ2 = 0.5, J4 = 0.3", [(0.5, p) for p in P4] + [(0.3, (0, 1, 2, 3))]),
         ("4-body, J4 -> infinity\nzeros at sqrt2 -+ 1", [(12.0, (0, 1, 2, 3))])]
fig, axs = plt.subplots(1, 4, figsize=(14, 4.4))
for ax, (t, terms) in zip(axs, cases):
    r = zeros(4, terms)
    ax.add_patch(Circle((0, 0), 1, fill=False, color=TXT2, lw=1, ls="--"))
    on = np.abs(np.abs(r) - 1) < 1e-7
    ax.scatter(r.real[on], r.imag[on], s=70, color=S[5], edgecolor="white", zorder=3, label="on |z| = 1")
    ax.scatter(r.real[~on], r.imag[~on], s=70, color=S[7], edgecolor="white", zorder=3, label="off the circle")
    ax.set_aspect("equal"); ax.set_xlim(-2.7, 2.7); ax.set_ylim(-2.7, 2.7)
    ax.axhline(0, color=GRID, lw=0.6); ax.axvline(0, color=GRID, lw=0.6)
    ax.set_title(t, fontsize=9.5); ax.set_xlabel("Re z"); ax.set_ylabel("Im z")
axs[0].legend(frameon=False, fontsize=8, loc="lower left")
fig.suptitle("W4. Lee-Yang zeros, 4 spins: pairs keep the circle, a lone 4-body parity term never does, pairs rescue it",
             x=0.01, ha="left", fontweight="bold")
fig.tight_layout(rect=[0, 0.03, 1, 0.92]); fig.savefig("figures/W4_lee_yang_zeros.png", dpi=150); plt.close(fig)

# ---------- W5 (J2, J4) map ----------
g = depth["D8_grid"]; G = np.array(g["max_abs_log_modulus"]); J2 = np.array(g["J2"]); J4 = np.array(g["J4"])
fig, ax = plt.subplots(figsize=(7, 5))
ok = (G < 1e-7).astype(float)
ax.contourf(J2, J4, ok, levels=[-0.5, 0.5, 1.5], colors=["#f6f6f6", "#cfe3f8"])
ax.contour(J2, J4, ok, levels=[0.5], colors=[S[0]], linewidths=2)
ax.axvline(math.log(3) / 8, color=S[1], lw=1.5, ls="--")
ax.text(math.log(3) / 8 + 0.02, 1.05, "J2 = ln3 / 8: the pair floor as J4 -> infinity", color=TXT, fontsize=8.5)
ax.text(0.7, 0.3, "Lee-Yang holds\n(every zero on |z| = 1)", color=TXT, fontsize=10)
ax.text(-0.55, 0.3, "fails", color=TXT2, fontsize=10)
ax.axhline(0, color=GRID, lw=0.6)
ax.set_xlabel("pair coupling J2 (all 6 pairs)"); ax.set_ylabel("4-body parity coupling J4")
ax.set_title("W5. 4 spins: a fixed pair amount carries any ferro 4-body term\n(exact for this symmetric model by Grace-Walsh-Szego)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("figures/W5_pair_floor_map.png", dpi=150); plt.close(fig)

# ---------- W6 rank test ----------
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for c, (name, d) in zip([S[0], S[2], S[1], S[7]], rank.items()):
    lv = [r["level"] for r in d["rows"]]; nm = [r["neighbour_maps"] for r in d["rows"]]
    ax.plot(lv, nm, marker="o", ms=7, color=c, lw=2, label=name)
    ax.text(lv[-1] + 0.15, nm[-1], f"{nm[-1]:,}", va="center", fontsize=8.5, color=TXT)
ax.set_yscale("log"); ax.set_xlabel("level (cut-off 3^-level)"); ax.set_ylabel("distinct neighbour maps (log)")
ax.legend(frameon=False, fontsize=8.5, loc="upper left"); ax.grid(color=GRID, lw=0.6)
ax.set_title("W6. Shape II: one stretch factor stays finite type; independent factors keep growing", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("figures/W6_rank_test.png", dpi=150); plt.close(fig)
print("figures written")
