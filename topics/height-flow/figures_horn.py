"""Pictures for the invariants section (DRAFT). Reads horn_invariants.json; writes figures/f9, f10."""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
H = json.load(open(os.path.join(HERE, "horn_invariants.json")))

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = "#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
                     "legend.frameon": False, "axes.titleweight": "bold", "axes.titlesize": 11})
def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)
def frac(s):
    a, b = (s.split("/") + ["1"])[:2]; return float(a)/float(b)

# ---------- f9: iterative residue against height and against the plane dial ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
s = np.linspace(0.15, 3.2, 400)
ax[0].plot(s, 1/(3*s), color=BLUE, label="résit(E^s) = 1/(3s)  (exact, from the formal series)")
ax[0].plot([1, 2], [1/3, 1/6], "o", ms=9, color=ORANGE, mec="white", mew=2, label="checked directly: E (1/3), E∘E (1/6)")
ax[0].axhline(0, color=INK2, lw=0.8)
ax[0].set_xlabel("height s"); ax[0].set_ylabel("iterative residue  résit")
ax[0].set_ylim(-0.05, 1.6); ax[0].legend(fontsize=8.5)
ax[0].set_title("A. Height × residue = 1/3 at every height")
x = np.linspace(-1.2, 3.2, 200)
ax[1].plot(x, x/3, color=BLUE, label="plane dial: résit = Δ/(3k²)")
pts = [(1.0, "e^x − 1  (d = 0)", ORANGE), (2.0, "k = 1, d = −1/4", AQUA), (0.0, "Δ = 0: Möbius,\ntrivial horns", VIOLET)]
for xx, lab, c in pts:
    ax[1].plot([xx], [xx/3], "o", ms=9, color=c, mec="white", mew=2, zorder=3)
    ax[1].annotate(lab, (xx, xx/3), xytext=(10, -22), textcoords="offset points", fontsize=8.5, color=INK)
ax[1].axhline(0, color=INK2, lw=0.8); ax[1].axvline(0, color=INK2, lw=0.8)
ax[1].set_xlabel("Δ / k²   (Δ = k² − 4d; the curvature is −Δ/4)"); ax[1].set_ylabel("iterative residue  résit")
ax[1].legend(fontsize=8.5, loc="upper left")
ax[1].set_title("B. The plane dial moves the residue; the corner dial does not")
save(fig, "f9_residue_height_and_dial.png")

# ---------- f10: half-rung gap and the same-fading-rate comparison ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
A1 = float(H["E"]["abs_A1"]); ImW = float(H["half_rung_gap_point"]["ImW"])
ss = np.linspace(0, 2, 400)
ax[0].plot(ss, 2*A1*np.abs(np.sin(np.pi*ss))*math.exp(-2*math.pi*ImW), color=BLUE,
           label="predicted from the horn map: 2|A₁| |sin πs| e^(−2π Im W)")
g = H["half_rung_gap"]
xs = sorted(g, key=float)
ax[0].plot([float(k) for k in xs], [float(g[k][0]) for k in xs], "o", ms=7, color=ORANGE, mec="white", mew=1.5,
           label="measured: attracting-side E^s minus repelling-side E^s")
ax[0].set_xlabel("height s"); ax[0].set_ylabel("gap |D| at one test point")
ax[0].legend(fontsize=8.2, loc="upper center"); ax[0].set_ylim(0, 9.5e-7)
ax[0].set_title("A. Fractional heights: the two halves of E^s disagree\nby exactly the horn coefficient (zero at whole s)")
T = H["table"]
short = {"e^x − 1 (our height-1 map)": "e^x − 1 (ours)", "Möbius x/(1 − x/2)": "Möbius",
         "x + x²/2 + x³/6 (e^x − 1 cut at x³)": "e^x − 1 cut at x³", "x·e^(x/2)": "x·e^(x/2)",
         "−log(1 − x)": "−log(1 − x)", "x + x²/2": "x + x²/2",
         "plane dial k = 1, d = −1/4": "plane k=1, d=−1/4", "plane dial k = 1, d = 1/4 (Δ = 0)": None}
floor = 1e-3
for name, r, a1, sf in T:
    lab = short.get(name, name)
    if lab is None: continue
    y = float(sf); y = floor if y < 1e-12 else y
    c = ORANGE if "ours" in lab else (VIOLET if lab == "Möbius" else BLUE)
    ax[1].plot([frac(r)], [y], "o", ms=9, color=c, mec="white", mew=2, zorder=3)
    off = {"e^x − 1 cut at x³": (-40, 10), "−log(1 − x)": (8, -4), "Möbius": (8, 2)}.get(lab, (8, 4))
    ax[1].annotate(lab + (": size 0, horn is a translation (drawn at floor)" if lab == "Möbius" else ""), (frac(r), y), xytext=off,
                   textcoords="offset points", fontsize=8.5, color=INK)
ax[1].set_yscale("log"); ax[1].set_xlim(-0.55, 1.45)
ax[1].set_xlabel("iterative residue  résit"); ax[1].set_ylabel("horn size  √|A₁·A₋₁|  (shift-free)")
ax[1].set_title("B. Same fading rate (−2/n): résit separates most,\nthe horn map separates the rest")
save(fig, "f10_gap_and_comparison.png")
