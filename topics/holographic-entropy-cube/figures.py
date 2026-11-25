"""Labeled pictures for the exploration. Reads results.json (run checks.py first)."""
import json, math, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc, Polygon
from functools import lru_cache

R = json.load(open("results.json"))
FIG = "figures/"
plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "figure.dpi": 130})
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#3b6fb6", "#e08a2c", "#3a9a5b", "#c8423b", "#7d5bb5", "#888888"
PC = {"A": BLUE, "B": ORANGE, "C": GREEN, "O": PURPLE}

# ---------------------------------------------------------------- fig 1: the entropy cube
def iso(p):  # 3D -> 2D oblique projection
    x, y, z = p
    return np.array([x + 0.45 * y, z + 0.32 * y])

fig, ax = plt.subplots(figsize=(7.4, 6.2))
labels = {(0,0,0): "∅", (1,0,0): "A", (0,1,0): "B", (0,0,1): "C",
          (1,1,0): "AB", (1,0,1): "AC", (0,1,1): "BC", (1,1,1): "ABC"}
V = list(labels)
for a, b in itertools.combinations(V, 2):
    if sum(abs(i - j) for i, j in zip(a, b)) == 1:
        pa, pb = iso(a), iso(b)
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="#bbbbbb", lw=2, zorder=1)
for v, lab in labels.items():
    p = iso(v); r = sum(v)
    if r == 0:
        col, txt, sign = "#dddddd", "S = 0", ""
    else:
        col = BLUE if r % 2 == 1 else RED
        sign = "+" if r % 2 == 1 else "−"
        txt = f"{sign} S$_{{{lab}}}$"
    ax.scatter(*p, s=1500, color=col, alpha=0.18 if r == 0 else 0.25, zorder=2, edgecolor=col, lw=2)
    ax.text(p[0], p[1], txt if r else "∅\n(S = 0)", ha="center", va="center", fontsize=13, zorder=3,
            color="black", weight="bold" if r else None)
    ax.text(p[0], p[1] - 0.2, f"weight (−1)$^{{{max(r-1,0)}}}$" if r else "", ha="center", va="center",
            fontsize=9, color=GRAY)
ax.set_title("Fig 1. Subsystem entropies sit on the corners of the 3-cube.\n"
             "The project's corner sum at k = −1 is the tripartite information I₃", fontsize=12)
ax.text(0.72, -0.48, "I₃ = S$_A$+S$_B$+S$_C$ − S$_{AB}$−S$_{AC}$−S$_{BC}$ + S$_{ABC}$\n"
        "= Σ over corners of k$^{\\,\\mathrm{rank}-1}$ S$_T$ at k = −1", ha="center", fontsize=11,
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax.text(1.75, 1.55, "blue = odd corners (+)\nred = even corners (−)", fontsize=10, color=GRAY)
ax.set_xlim(-0.45, 2.05); ax.set_ylim(-0.8, 1.75); ax.axis("off")
fig.tight_layout(); fig.savefig(FIG + "fig1_entropy_cube.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 2: RT surfaces in the hyperbolic disk
def geod(t1, t2, eps=1e-3):
    return 2 * math.log((2 / eps) * abs(math.sin((t2 - t1) / 2)))

def rt(endpoints):
    pts = sorted(endpoints)
    @lru_cache(None)
    def best(i, j):
        if i > j: return 0.0, ()
        out = (math.inf, ())
        for m in range(i + 1, j + 1, 2):
            a, pa = best(i + 1, m - 1); b, pb = best(m + 1, j)
            tot = geod(pts[i], pts[m]) + a + b
            if tot < out[0]: out = (tot, ((pts[i], pts[m]),) + pa + pb)
        return out
    return best(0, len(pts) - 1)

def geodesic_path(t1, t2, n=120):
    """Points of the hyperbolic geodesic (circle orthogonal to the unit circle) between boundary angles."""
    p1, p2 = np.exp(1j * t1), np.exp(1j * t2)
    d = (t2 - t1) % (2 * math.pi)
    if abs(d - math.pi) < 1e-6:
        s = np.linspace(0, 1, n); return p1 + (p2 - p1) * s
    mid = (t1 + d / 2)
    c = np.exp(1j * mid) / math.cos(d / 2)  # center of orthogonal circle
    r = abs(math.tan(d / 2))
    a1, a2 = np.angle(p1 - c), np.angle(p2 - c)
    da = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
    ang = a1 + da * np.linspace(0, 1, n)
    return c + r * np.exp(1j * ang)

def disk(ax, title):
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=1.5, color="black"))
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11)

def interval(ax, name, t1, t2, col):
    th = np.linspace(t1, t2, 80)
    ax.plot(np.cos(th), np.sin(th), color=col, lw=7, solid_capstyle="butt")
    tm = (t1 + t2) / 2
    ax.text(1.14 * math.cos(tm), 1.14 * math.sin(tm), name, color=col, fontsize=14, weight="bold",
            ha="center", va="center")

def chord(ax, t1, t2, col, lw=2.5, ls="-"):
    z = geodesic_path(t1, t2); ax.plot(z.real, z.imag, color=col, lw=lw, ls=ls)

fig, axs = plt.subplots(1, 3, figsize=(14, 5.2))
iv = {"A": (0.25, 1.35), "B": (1.95, 3.05), "C": (3.75, 5.35)}
disk(axs[0], "(a) one region each: S = (length of\nthe shortest bulk curve) / 4G")
for p, (a, b) in iv.items():
    interval(axs[0], p, a, b, PC[p]); chord(axs[0], a, b, PC[p])
# (b) AB far apart: disconnected phase
ivb = {"A": (0.3, 1.0), "B": (2.4, 3.1)}
disk(axs[1], "(b) A and B far apart: two separate\ncurves, so I(A:B) = 0")
for p, (a, b) in ivb.items(): interval(axs[1], p, a, b, PC[p])
for t1, t2 in rt(list(ivb["A"]) + list(ivb["B"]))[1]: chord(axs[1], t1, t2, "black")
# (c) AB close: connected phase, shade the entanglement wedge
ivc = {"A": (0.2, 1.45), "B": (1.7, 3.0)}
disk(axs[2], "(c) A and B close: the curves reconnect\n(entanglement wedge, shaded), I(A:B) > 0")
for p, (a, b) in ivc.items(): interval(axs[2], p, a, b, PC[p])
ch = rt(list(ivc["A"]) + list(ivc["B"]))[1]
for t1, t2 in ch: chord(axs[2], t1, t2, "black")
outer = geodesic_path(0.2, 3.0); inner = geodesic_path(1.45, 1.7)
th1 = np.linspace(0.2, 1.45, 60); th2 = np.linspace(1.7, 3.0, 60)
poly = np.concatenate([np.exp(1j * th1), inner, np.exp(1j * th2), outer[::-1]])
axs[2].add_patch(Polygon(np.c_[poly.real, poly.imag], color=GRAY, alpha=0.18))
fig.suptitle("Fig 2. Ryu–Takayanagi in the hyperbolic disk (a slice of anti-de Sitter space): "
             "entanglement entropy = area of the minimal surface", fontsize=13, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(FIG + "fig2_rt_disk.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 3: the contraction map proof of MMI, and the failure of its reverse
fig, axs = plt.subplots(1, 2, figsize=(14, 6.4))
m = R["B"]["MMI S_AB + S_BC + S_AC >= S_A + S_B + S_C + S_ABC"]["map"]
bnd = {"000": "O", "101": "A", "110": "B", "011": "C"}
ax = axs[0]
V = list(itertools.product((0, 1), repeat=3))
for a, b in itertools.combinations(V, 2):
    if sum(abs(i - j) for i, j in zip(a, b)) == 1:
        pa, pb = iso(a), iso(b)
        fa, fb = m["".join(map(str, a))], m["".join(map(str, b))]
        dimg = sum(x != y for x, y in zip(fa, fb))
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=GREEN if dimg == 1 else "#cfcfcf",
                lw=3 if dimg == 1 else 2, ls="-" if dimg == 1 else "--", zorder=1)
for v in V:
    key = "".join(map(str, v)); p = iso(v)
    col = PC[bnd[key]] if key in bnd else GRAY
    ax.scatter(*p, s=700, color=col, alpha=0.3, edgecolor=col, lw=2, zorder=2)
    ax.text(p[0], p[1] + 0.03, key, ha="center", va="center", fontsize=11, weight="bold", zorder=3)
    ax.text(p[0] + 0.07, p[1] - 0.13, "→ " + m[key] + (f"  ({bnd[key]})" if key in bnd else ""), ha="left",
            fontsize=9.5, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9), color=col if key in bnd else "black", zorder=3)
ax.set_title("(a) MMI: S$_{AB}$+S$_{BC}$+S$_{AC}$ ≥ S$_A$+S$_B$+S$_C$+S$_{ABC}$\n"
             "Domain: 3-cube, bits = (AB, BC, AC).  Image bits = (A, B, C, ABC)\n"
             "The one and only contraction: all four inner corners go to 0001", fontsize=11)
ax.text(0.7, -0.62, "Every one of the 12 edges lands on an edge of the image 4-cube.\n"
        "No edge is stretched, so every bulk geometry obeys MMI.", ha="center", fontsize=10,
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax.set_xlim(-0.45, 2.0); ax.set_ylim(-0.85, 1.7); ax.axis("off")

ax = axs[1]  # reversed MMI: the corner 0001 has nowhere to go
# show the 3-cube image space (bits AB, BC, AC) and the four points it must stay within distance 1 of
targets = {"000": "O", "101": "A", "110": "B", "011": "C"}
for a, b in itertools.combinations(V, 2):
    if sum(abs(i - j) for i, j in zip(a, b)) == 1:
        pa, pb = iso(a), iso(b); ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="#cfcfcf", lw=2, zorder=1)
for v in V:
    key = "".join(map(str, v)); p = iso(v)
    near = [targets[t] for t in targets if sum(x != y for x, y in zip(key, t)) <= 1]
    col = PC[targets[key]] if key in targets else GRAY
    ax.scatter(*p, s=900, color=col, alpha=0.3, edgecolor=col, lw=2, zorder=2)
    ax.text(p[0], p[1] + 0.03, key, ha="center", va="center", fontsize=11, weight="bold", zorder=3)
    ax.text(p[0], p[1] - 0.17, "near " + ",".join(near), ha="center", fontsize=9, zorder=3)
ax.set_title("(b) The reverse (false) inequality. Shown: its image 3-cube, bits (AB, BC, AC).\n"
             "The domain corner 0001 (only the ABC bit on) is one step from A, B, C and O,\n"
             "so its image must be one step from all four colored corners. None is.", fontsize=11)
ax.text(0.7, -0.62, "So an edge must stretch. Bao–Furuya–Naskar (2025): a stretched edge\n"
        "builds a bulk geometry that breaks the inequality.", ha="center", fontsize=10,
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax.set_xlim(-0.45, 2.0); ax.set_ylim(-0.85, 1.7); ax.axis("off")
fig.suptitle("Fig 3. Proof by hypercube: an entropy inequality holds for every geometry exactly when "
             "one cube maps onto another without stretching an edge", fontsize=13)
fig.tight_layout(); fig.savefig(FIG + "fig3_contraction_map.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 4: who passes the cube test
fig, ax = plt.subplots(figsize=(11, 5.6))
data = [("random graphs, min cut (1500)", np.array(R["C_I3_graphs"]), PURPLE),
        ("AdS₃ geodesics, 3 random intervals (3000)", np.array(R["E_I3_ads"]), BLUE),
        ("random 4-qubit states, Haar (4000)", np.array(R["D_I3_haar"]), GRAY),
        ("sparse 4-qubit states, 2–4 terms (4000)", np.array(R["D_I3_sparse"]), ORANGE)]
for i, (lab, x, col) in enumerate(data):
    x = np.clip(x, -4, 1.5)
    y = i + 0.32 * (np.random.default_rng(i).random(len(x)) - 0.5)
    ax.scatter(x, y, s=5, alpha=0.35, color=col)
    viol = (np.array(data[i][1]) > 1e-9).mean() * 100
    ax.text(1.62, i, f"{viol:.1f}% break MMI", va="center", fontsize=10, color=RED if viol > 0 else GREEN)
ax.axvline(0, color="black", lw=1.5)
ax.axvspan(0, 1.55, color=RED, alpha=0.06)
ax.text(0.75, 3.55, "I₃ > 0: no geometry\ncan carry this pattern", ha="center", color=RED, fontsize=10)
ax.text(-2.0, 3.55, "I₃ ≤ 0 (MMI)", ha="center", color=GREEN, fontsize=10)
for name, val, col in [("GHZ", R["D_named"]["GHZ"], RED), ("W", R["D_named"]["W"], RED)]:
    ax.plot(val, 2.5, marker="v", ms=11, color=col); ax.text(val, 2.72, name, ha="center", color=col, weight="bold")
ax.set_yticks(range(4)); ax.set_yticklabels([d[0] for d in data])
ax.set_xlim(-4.1, 2.4); ax.set_ylim(-0.6, 3.95)
ax.set_xlabel("tripartite information I₃ (the k = −1 corner sum), in bits or area units (clipped at −4)")
ax.set_title("Fig 4. Geometric entropies (min cuts, AdS geodesics) never cross zero;\nquantum states can, and then no geometry can carry them")
fig.tight_layout(); fig.savefig(FIG + "fig4_who_passes.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 5: cutoff
sw = R["E_sweep"]
L = [math.log(1 / s["eps"]) for s in sw]
fig, ax = plt.subplots(figsize=(8.6, 5))
ax.plot(L, [s["S_A"] for s in sw], color=BLUE, lw=2.5, label="S$_A$ (one corner): grows without bound")
ax.plot(L, [s["S_ABC"] for s in sw], color=GRAY, lw=2, label="S$_{ABC}$ (one corner): grows faster")
ax.plot(L, [s["I2_AB"] for s in sw], color=ORANGE, lw=2.5, label="I(A:B) (k = −1 sum on a square face): flat")
ax.plot(L, [s["I3"] for s in sw], color=GREEN, lw=2.5, label="I₃ (k = −1 sum on the whole cube): flat")
ax.axhline(0, color="black", lw=0.8)
ax.set_xlabel("log(1/cutoff)  →  finer and finer resolution at the boundary")
ax.set_ylabel("entropy (units of c/6)")
ax.set_title("Fig 5. Single corners diverge; the k = −1 corner sums do not")
ax.legend(fontsize=9.5, loc="upper left"); fig.tight_layout(); fig.savefig(FIG + "fig5_cutoff.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 6: why k = -1
k = np.linspace(-2.2, 1.2, 300)
fig, ax = plt.subplots(figsize=(7.6, 4.6))
for n, col in [(2, ORANGE), (3, GREEN), (4, BLUE)]:
    ax.plot(k, n * (1 + k) ** (n - 1), color=col, lw=2.5, label=f"n = {n}:  n(1+k)$^{{{n-1}}}$")
pts = R["A_additive_vs_k"]
ax.scatter([float(a) for a in pts], list(pts.values()), color="black", zorder=3, s=25, label="checks.py, n = 3")
ax.axvline(-1, color=RED, ls="--"); ax.axhline(0, color="black", lw=0.8)
ax.text(-1.05, 9, "k = −1:\nthe only zero", ha="right", color=RED)
ax.set_ylim(-4, 14); ax.set_xlabel("dial k"); ax.set_ylabel("cube sum of a one-body entropy")
ax.set_title("Fig 6. Only at k = −1 does the cube sum forget\nthe per-region (cutoff) pieces")
ax.legend(fontsize=9.5); fig.tight_layout(); fig.savefig(FIG + "fig6_why_k_minus_1.png"); plt.close(fig)
print("figures written")
