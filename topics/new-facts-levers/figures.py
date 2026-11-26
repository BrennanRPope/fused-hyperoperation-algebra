"""Pictures for topics/new-facts-levers (DRAFT). Run: python3 figures.py -> figures/*.png
Colour key: blue = ours (new or not found in print), orange = already known, grey = background, red = fails."""
import json, os, itertools
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
BLUE, ORANGE, GREY, RED, GREEN = "#1f5fbf", "#e07b00", "#8a8a8a", "#c0392b", "#2e8b57"
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "figure.dpi": 130})
os.makedirs("figures", exist_ok=True)
def save(fig, name):
    fig.tight_layout(); fig.savefig("figures/" + name, dpi=130); plt.close(fig); print("figures/" + name)
TH = json.load(open("thresholds.json")); FL = json.load(open("flowers.json")); FA = json.load(open("families.json"))
RT = json.load(open("rtn_check.json"))

# ---------- fig0: the lever board ----------
levers = [
 ("Corner-sum threshold\nk ≥ √6 − 3", "numbers not found;\nread off known cones", "small: a GHZ detector\nfrom 5 regions on", "none"),
 ("Union half √2 − 1\n(3-party sharing starts)", "not found; 2-line\ncorollary (folklore?)", "no: a known phase\ntransition, re-labelled", "none"),
 ("Fuzzy monogamy\n(soft min cuts keep MMI)", "not found in\n4 searches", "yes: geometry's rule\nsurvives finite D (annealed)", "none"),
 ("Lee–Yang on the ℓp ball\nn ≤ 1 + p²", "rule not found;\npartly proved", "no", "Riemann: a test family\nfor Newman's route (weak)"),
 ("Rewind ladder k = n(n+1);\ntunnelling k = −2K", "a reading of\nknown pieces", "no: shared SU(2)\nCasimir only", "none"),
 ("e from the cube's\ncorner simplex", "known (Freudenthal,\nrandom-sum puzzle)", "no", "none"),
]
fig, ax = plt.subplots(figsize=(15, 8.6)); ax.axis("off")
cols = ["Lever", "New?", "Quantum gravity / ToE", "Millennium problems"]
xs = [0.0, 0.27, 0.49, 0.75]; ws = [0.26, 0.21, 0.25, 0.25]
for x, w, c in zip(xs, ws, cols): ax.text(x + w / 2, 0.97, c, ha="center", va="center", fontsize=16, weight="bold")
def colour(col, txt):
    t = txt.lower()
    if col == 1: return BLUE if "not found" in t or "numbers" in t else ORANGE
    if t.startswith("yes"): return GREEN
    if t.startswith("small") or "weak" in t: return "#d4a017"
    return GREY
for i, row in enumerate(levers):
    y = 0.86 - i * 0.145
    for j, (x, w, txt) in enumerate(zip(xs, ws, row)):
        c = "#ffffff" if j == 0 else colour(j, txt)
        ax.add_patch(FancyBboxPatch((x, y - 0.06), w, 0.12, boxstyle="round,pad=0.005", fc=c, ec="k", alpha=0.25 if j else 1))
        ax.text(x + w / 2, y, txt, ha="center", va="center", fontsize=12.5)
ax.text(0.5, 0.02, "Green = a lead worth following · gold = small or weak · grey = no bridge · blue = not found in print · orange = known",
        ha="center", fontsize=13, color="k")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
save(fig, "fig0-lever-board.png")

# ---------- fig1: threshold ladder, holographic vs all quantum ----------
fig, ax = plt.subplots(figsize=(12, 6.4))
ns = [3, 4, 5, 6]
q = [TH[str(n)]["vN"] for n in ns]
hol_exact = {3: TH["3"]["vN"], 4: TH["4"]["vN"], 5: (np.sqrt(5) - 3) / 2}
for n in ns:
    ax.plot([n, n], [-0.6, 0.0], color=GREY, lw=0.6)
ax.plot(ns, q, "o-", color=ORANGE, ms=11, lw=2.5, label="every quantum state (exact: an LP bound met by a real state)")
ax.plot(list(hol_exact), list(hol_exact.values()), "s-", color=BLUE, ms=11, lw=2.5, label="geometric (holographic) states, exact")
ax.plot([6], [FL["6"][0]], "s", color=BLUE, mfc="white", ms=11, mew=2.5, label="geometric, best known example (flower graph)")
ax.fill_between([4.6, 5.4], [hol_exact[5]] * 2, [q[2]] * 2, color=GREEN, alpha=0.35)
ax.annotate("geometry window:\nevery geometric state ≥ 0,\nsome quantum state < 0", xy=(5.05, -0.368), xytext=(5.35, -0.50),
            fontsize=13, arrowprops=dict(arrowstyle="->"))
labels = {3: "√6 − 3 = −0.5505", 4: "(√5 − 3)/2 = −1/φ² = −0.3820", 5: "−0.3538 (quantum)\n−0.3820 (geometric)", 6: "−0.2882 (quantum)"}
for n in ns: ax.text(n - 0.08, q[ns.index(n)] + 0.025, labels[n], fontsize=11.5, ha="right")
ax.set_xticks(ns); ax.set_xlabel("number of named regions n (plus the rest of the system)")
ax.set_ylabel("threshold k*: the corner sum is ≥ 0 for all k ≥ k*")
ax.set_title("Where the weighted corner sum turns safe: geometry only shows up from five regions on")
ax.set_ylim(-0.6, -0.2); ax.set_xlim(2.6, 6.6); ax.legend(loc="lower right")
save(fig, "fig1-threshold-ladder.png")

# ---------- fig2: the two worst states at n = 5 ----------
k = np.linspace(-0.5, -0.25, 400)
pt = (1 + k) * (4 * k**2 + 12 * k + 4)
mds = k**4 + 10 * k**3 + 20 * k**2 + 20 * k + 5
fig, ax = plt.subplots(figsize=(11, 6))
ax.axhline(0, color="k", lw=1)
ax.plot(k, pt / 4, color=BLUE, lw=3, label="perfect tensor on 4 regions (geometric worst), ÷ 4")
ax.plot(k, mds / 5, color=ORANGE, lw=3, label="Reed–Solomon GHZ state on 6 parties (quantum worst), ÷ 5")
ax.axvspan((np.sqrt(5) - 3) / 2, -0.35383, color=GREEN, alpha=0.3, label="geometry window")
ax.set_xlabel("k"); ax.set_ylabel("corner sum (per unit of single-region entropy)")
ax.set_title("Five regions: the GHZ-type state goes negative where no geometric state can")
ax.legend(loc="upper left"); ax.set_ylim(-0.25, 0.35)
save(fig, "fig2-five-region-window.png")

# ---------- fig3: the two extremal states ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 6.2))
ax = axs[0]; ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(Circle((0, 0), 0.18, color=BLUE)); 
for i, lab in enumerate("ABCD"):
    a = np.pi / 4 + i * np.pi / 2; x, y = np.cos(a), np.sin(a)
    ax.plot([0, x], [0, y], color=BLUE, lw=4); ax.add_patch(Circle((x, y), 0.14, fc="white", ec="k", lw=2))
    ax.text(x, y, lab, ha="center", va="center", fontsize=16)
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.5, 1.4)
ax.set_title("Geometric worst: perfect tensor (star graph)\nany 2 regions hold 2 units: as entangled as possible", fontsize=13)
ax = axs[1]; ax.set_aspect("equal"); ax.axis("off")
pts = [(np.cos(np.pi / 2 + i * np.pi / 3), np.sin(np.pi / 2 + i * np.pi / 3)) for i in range(6)]
for (x1, y1), (x2, y2) in itertools.combinations(pts, 2): ax.plot([x1, x2], [y1, y2], color=ORANGE, alpha=0.25, lw=1.5)
for (x, y), lab in zip(pts, ["A", "B", "C", "D", "E", "rest"]):
    ax.add_patch(Circle((x, y), 0.16, fc="white", ec="k", lw=2)); ax.text(x, y, lab, ha="center", va="center", fontsize=13)
ax.text(0, 0, "shared\nsecret\n(a, b)", ha="center", va="center", fontsize=13, color=ORANGE)
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.5, 1.4)
ax.set_title("Quantum worst: Reed–Solomon GHZ state\nany 2 parties know the whole secret; more add nothing", fontsize=13)
save(fig, "fig3-extremal-states.png")

# ---------- fig4: large n trend ----------
fig, ax = plt.subplots(figsize=(11, 6))
nn = np.arange(4, 13)
ax.plot(nn, [FL[str(n)][0] for n in nn], "s-", color=BLUE, lw=2.5, ms=9, label="geometric: best flower graph (lower bound on k*)")
ax.plot(nn, [FA["MDS"][str(n)][0] for n in nn], "o-", color=ORANGE, lw=2.5, ms=9, label="quantum: Reed–Solomon GHZ state, pure on n regions")
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("number of regions n"); ax.set_ylabel("largest k where the corner sum can go negative")
ax.set_title("Both thresholds creep to 0, so no single k separates geometry at every size")
ax.legend(loc="lower right"); ax.set_ylim(-0.42, 0.02)
save(fig, "fig4-large-n.png")

# ---------- fig5: fuzzy monogamy on real random tensor networks ----------
keys = list(RT.keys())
fig, axs = plt.subplots(1, 2, figsize=(14, 6))
x = np.arange(len(keys)); w = 0.27
axs[0].bar(x - w, [RT[k_]["annealed_I3_renyi2"] for k_ in keys], w, color=BLUE, label="annealed Rényi-2 (the theorem)")
axs[0].bar(x, [RT[k_]["quenched_mean_I3_renyi2"] for k_ in keys], w, color="#7fa7e0", label="quenched Rényi-2, average")
axs[0].bar(x + w, [RT[k_]["mean_I3_vN"] for k_ in keys], w, color=GREY, label="von Neumann, average")
axs[0].axhline(0, color="k", lw=1); axs[0].set_xticks(x); axs[0].set_xticklabels([k_.replace(",", "\n") for k_ in keys], fontsize=10)
axs[0].set_ylabel("tripartite information I₃ (≤ 0 is monogamy)"); axs[0].legend(fontsize=11, loc="lower right")
axs[0].set_title("On average, monogamy holds at bond dimension 2 and 3")
axs[1].bar(x - w / 2, [RT[k_]["frac_samples_I3_renyi2_pos"] for k_ in keys], w, color=RED, label="Rényi-2")
axs[1].bar(x + w / 2, [RT[k_]["frac_samples_I3_vN_pos"] for k_ in keys], w, color="#e8a09a", label="von Neumann")
axs[1].set_xticks(x); axs[1].set_xticklabels([k_.replace(",", "\n") for k_ in keys], fontsize=10)
axs[1].set_ylabel("share of single samples that break monogamy"); axs[1].legend()
axs[1].set_title("Single samples break it, less as D grows")
save(fig, "fig5-fuzzy-monogamy-networks.png")

# ---------- fig6: lever x Millennium grid ----------
probs = ["Quantum\ngravity", "Riemann", "Yang–Mills", "Navier–\nStokes", "P vs NP", "Hodge", "BSD"]
names = ["Corner-sum\nthreshold", "Union half\n√2 − 1", "Fuzzy\nmonogamy", "Lee–Yang\nℓp rule", "Rewind ladder /\ntunnelling", "e from the\ncorner simplex"]
grid = np.array([[1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
from matplotlib.colors import ListedColormap
fig, ax = plt.subplots(figsize=(12, 6.5))
ax.imshow(grid, cmap=ListedColormap(["#eeeeee", "#f2d16b", "#7fc79a"]), vmin=0, vmax=2)
sym = {0: "no bridge", 1: "small", 2: "lead"}
for i in range(grid.shape[0]):
    for j in range(grid.shape[1]): ax.text(j, i, sym[grid[i, j]], ha="center", va="center", fontsize=11)
ax.set_xticks(range(len(probs))); ax.set_xticklabels(probs); ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
ax.set_title("Each new fact against each goal, judged on its own merits")
save(fig, "fig6-goal-grid.png")
