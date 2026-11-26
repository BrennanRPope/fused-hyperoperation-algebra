"""Pictures for the dial map (DRAFT). Run: python3 figures.py  ->  figures/*.png"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
from results import DIALS, RESULTS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
DCOL = {d[0]: d[2] for d in DIALS}
NEWCOL = {"new": "#1a7f37", "reading": "#222222", "known": "#8a8a8a"}
BINNAME = {"op": "Uses ∘k itself (or its cube of corners)",
           "plane": "Uses the plane law, or composition of maps with different centres",
           "pic": "Uses only the dial picture: the shape family or a reading sets the question"}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, x, y, w, h, text, fc, ec="#444", fs=10.5, weight="normal", color="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc=fc, ec=ec, lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, weight=weight,
            color=color, wrap=True)


# ---------------------------------------------------------------- fig 1: the order of maths
def fig1():
    fig, ax = plt.subplots(figsize=(18, 12.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.995, "Where the theory sits in the order of mathematics", ha="center", va="top",
            fontsize=18, weight="bold")
    ax.text(0.5, 0.962, "Set theory and category theory order any field in the same five steps. "
            "Each row shows the matching piece here and the result that backs it.",
            ha="center", va="top", fontsize=11.5, color="#444")
    levels = [
        ("1  Sets", "Elements, subsets,\nthe power set.",
         "The box: ∏(1 + tₓ) has one term per subset.\nThe cube's corners are the power set;\nt is each element's odds.",
         "Proved in Lean: 48 theorems\n(fractional-powerset-dials).\nOpen, ours: a 'half power set'\nneeds the absence of choice.",
         "#fff4e6"),
        ("2  Objects", "A set with structure.\nA monoid is a\none-object category.",
         "(ℝ, ∘k) is a monoid: arrows are numbers,\n0 is the identity. The shift x ↦ 1 + kx carries\n∘k to ordinary ×, so every dial k ≠ 0 is the\nsame object in new units.",
         "Proved in Lean:\nmathlib-map (19 theorems),\norder-dial (33 theorems).",
         "#e8f1fb"),
        ("3  Arrows", "Structure-preserving maps\ncarry facts from one\nobject to another.",
         "Maps out of the dial's group land in other fields:\nvelocity addition, 'A or B' in probability,\nSugeno measures, fermion and boson point\nprocesses, Hirzebruch genera, prime products.",
         "Each one is known maths.\nThe single family holding\nthem is the placement\n(README, Where it unifies).",
         "#eaf6ea"),
        ("4  Universal\n    properties", "An object pinned down\nby how every other\nobject maps to it.",
         "∘₁ is the categorical product of partial functions\n(pointed sets carry it to ×). The Jacobson radical\nis exactly the elements that never meet a wall\nat any dial.",
         "Proved in Lean: category-theory\n(36 theorems), mathlib-map.\nKnown, not built here: Lazard's\nuniversal formal group law.",
         "#f3eafa"),
        ("5  Families\n    (moduli)", "A space whose points\nare objects; invariants\nare what isomorphism keeps.",
         "The dials k, d, p, C, t, r are coordinates on\nspaces of laws and shapes. Readings such as\nΔ = k² − 4d are invariants: the sign of Δ is\nexactly what sorts the laws into three kinds.",
         "Known: Coleman–McGuinness\n1991; Lie and Klein. The rule\n'dials in, readings out' is this\nhabit (affine-names).",
         "#fdeef0"),
    ]
    xs = [(0.005, 0.11), (0.125, 0.19), (0.33, 0.38), (0.725, 0.27)]
    h, gap, ytop = 0.148, 0.024, 0.885
    for x, t in [(0.22, "How foundations order maths"),
                 (0.52, "Fused Hyperoperation Algebra +\nAffine Composition Theory"),
                 (0.86, "What backs it")]:
        ax.text(x, 0.915, t, ha="center", va="center", fontsize=12.5, weight="bold")
    for i, (lvl, found, ours, backing, col) in enumerate(levels):
        y = ytop - (i + 1) * h - i * gap
        (x0, w0), (x1, w1), (x2, w2), (x3, w3) = xs
        box(ax, x0, y, w0, h, lvl, col, fs=13, weight="bold")
        box(ax, x1, y, w1, h, found, "white", fs=11)
        box(ax, x2, y, w2, h, ours, col, fs=11)
        box(ax, x3, y, w3, h, backing, "white", ec="#999", fs=10.3, color="#333")
        if i < len(levels) - 1:
            ax.annotate("", xy=(x0 + w0 / 2, y - gap + 0.004), xytext=(x0 + w0 / 2, y - 0.004),
                        arrowprops=dict(arrowstyle="->", color="#666", lw=1.6))
    save(fig, "fig1-order-of-maths.png")


# ---------------------------------------------------------------- fig 2: the main map
def fig2():
    order = ["op", "plane", "pic"]
    rows = []
    for b in order:
        rows.append(("HEADER", b))
        rows += [r for r in RESULTS if r[2] == b]
    n = len(rows)
    fig, ax = plt.subplots(figsize=(17, 0.47 * n + 2.2))
    ax.set_xlim(-0.05, 1.62); ax.set_ylim(-1.5, n + 1.2); ax.axis("off")
    ys = {}
    yy = n
    for r in rows:
        if r[0] == "HEADER":
            ax.text(0.42, yy, BINNAME[r[1]], fontsize=12.5, weight="bold", va="center",
                    color={"op": "#0b3d91", "plane": "#7a1f1f", "pic": "#4b4b4b"}[r[1]])
            ax.plot([0.42, 1.60], [yy - 0.42, yy - 0.42], color="#bbb", lw=0.8)
        else:
            ys[r] = yy
        yy -= 1
    # dial hubs on the left, spread evenly
    nd = len(DIALS)
    hub_y = {d[0]: n - 1 - i * (n - 2) / (nd - 1) for i, d in enumerate(DIALS)}
    for letter, name, col in DIALS:
        y = hub_y[letter]
        ax.scatter([0.235], [y], s=1500, color=col, zorder=5)
        ax.text(0.235, y, letter, ha="center", va="center", color="white", fontsize=17, weight="bold", zorder=6)
        ax.text(0.195, y, name, ha="right", va="center", fontsize=10.5, color=col, weight="bold")
    # edges
    for r, y in ys.items():
        for letter in r[1]:
            ax.plot([0.235, 0.40], [hub_y[letter], y], color=DCOL[letter], lw=1.3, alpha=0.55, zorder=1)
    # result nodes
    for r, y in ys.items():
        name, dials, b, status, new, folder = r
        if status == "lean":
            ax.scatter([0.405], [y], marker="s", s=120, color="#222", zorder=4)
        elif status == "checked":
            ax.scatter([0.405], [y], marker="o", s=120, color="#222", zorder=4)
        else:
            ax.scatter([0.405], [y], marker="o", s=120, facecolor="white", edgecolor="#222", lw=1.6, zorder=4)
        ax.text(0.43, y, name, va="center", fontsize=10.6, color=NEWCOL[new],
                weight="bold" if new == "new" else "normal")
        ax.text(1.60, y, folder, va="center", ha="right", fontsize=9, color="#666", family="DejaVu Sans Mono")
    # legend (two lines)
    ly, ly2 = -0.45, -1.05
    ax.scatter([0.43], [ly], marker="s", s=90, color="#222"); ax.text(0.445, ly, "proved in Lean", va="center", fontsize=10)
    ax.scatter([0.62], [ly], marker="o", s=90, color="#222"); ax.text(0.635, ly, "checked (exact or numeric)", va="center", fontsize=10)
    ax.scatter([0.90], [ly], marker="o", s=90, facecolor="white", edgecolor="#222", lw=1.5)
    ax.text(0.915, ly, "partly proved, rest open", va="center", fontsize=10)
    ax.text(0.43, ly2, "green bold = not found in print (write-up's own search)", va="center", fontsize=10, color=NEWCOL["new"], weight="bold")
    ax.text(0.95, ly2, "black = known result, our placement", va="center", fontsize=10)
    ax.text(1.28, ly2, "grey = known", va="center", fontsize=10, color=NEWCOL["known"])
    ax.text(0.0, n + 0.9, "The dial map: each dial, the results it drives, and how much of the fused operation each one really uses",
            fontsize=16, weight="bold")
    save(fig, "fig2-dial-map.png")


# ---------------------------------------------------------------- fig 3: one centre, two centres, curvature
def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.4))
    # (a) one centre
    ax = axs[0]
    x = np.linspace(-2, 3, 50)
    c = -1.0  # centre at -1/k with k = 1
    for a in [-0.6, -0.3, 0.5, 1.0, 2.0]:
        ax.plot(x, (1 + a) * x + a, color="#1f77b4", alpha=0.8)
    ax.plot(x, x, "k--", lw=1)
    ax.scatter([c], [c], s=80, color="#1f77b4", zorder=5)
    ax.annotate("one centre, the wall −1/k", (c, c), (-1.9, 2.2), arrowprops=dict(arrowstyle="->"), fontsize=10)
    ax.set_title("(a) Fused algebra: scalings about one centre.\nThey commute; composing them is a ∘k b.", fontsize=11.5)
    ax.set_xlabel("before"); ax.set_ylabel("after"); ax.set_xlim(-2, 3); ax.set_ylim(-2.5, 4)
    # (b) two centres
    ax = axs[1]
    f = lambda v: 1.1 * v          # raise by 10% (centre 0)
    g = lambda v: v + 500          # grant of 500 (centre at infinity)
    start = 1000
    pts1 = [start, f(start), g(f(start))]
    pts2 = [start, g(start), f(g(start))]
    ax.plot([0, 1, 2], pts1, "-o", color="#1f77b4", label="raise, then grant: 1600")
    ax.plot([0, 1, 2], pts2, "-o", color="#d62728", label="grant, then raise: 1650")
    ax.annotate("", xy=(2.05, 1650), xytext=(2.05, 1600), arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text(2.1, 1625, "50 = the raise on\nthe grant (a corner)", va="center", fontsize=10)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["start", "after one", "after two"])
    ax.set_xlim(-0.2, 2.9); ax.legend(loc="upper left", fontsize=10)
    ax.set_title("(b) Different centres: order matters by an oriented\narea. This is where Rulecraft's 'history' lives.", fontsize=11.5)
    # (c) curvature: three kinds
    ax = axs[2]
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(t) - 2.4, np.sin(t), color="#2ca02c", lw=2)
    xx = np.linspace(-1, 1, 100)
    ax.plot(xx, xx ** 2 - 0.6, color="#555", lw=2)
    u = np.linspace(-1.1, 1.1, 100)
    ax.plot(2.4 + np.cosh(u) - 1, np.sinh(u), color="#d62728", lw=2)
    ax.plot(2.4 - np.cosh(u) + 0.3, np.sinh(u), color="#d62728", lw=2)
    ax.text(-2.4, -1.25, "Δ < 0\nturning\n(tan addition,\nrungs)", ha="center", va="top", fontsize=9.5)
    ax.text(0, -1.25, "Δ = 0\none double wall\n(addition in\ndisguise)", ha="center", va="top", fontsize=9.5)
    ax.text(2.55, -1.25, "Δ > 0\ntwo walls\n(relativity,\nscreening)", ha="center", va="top", fontsize=9.5)
    ax.set_xlim(-3.8, 3.9); ax.set_ylim(-3.0, 1.5); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(c) The curvature dial d adds the plane law; the sign of\nΔ = k² − 4d sorts every law into three kinds.", fontsize=11.5)
    fig.suptitle("Affine Composition Theory = the fused algebra (one centre) + curvature + maps with different centres",
                 fontsize=14, weight="bold", y=1.03)
    save(fig, "fig3-one-centre-two-centres.png")


# ---------------------------------------------------------------- fig 4: the (k, d) plane with placements
def fig4():
    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    k = np.linspace(-3, 3, 400)
    ax.fill_between(k, k ** 2 / 4, 3, color="#2ca02c", alpha=0.10)
    ax.fill_between(k, -3, k ** 2 / 4, color="#d62728", alpha=0.07)
    ax.plot(k, k ** 2 / 4, color="#333", lw=2)
    ax.text(-2.85, 2.6, "Δ < 0: no real wall, turning", color="#1d6b1d", fontsize=11)
    ax.text(0.7, -2.7, "Δ > 0: two real walls, saturating", color="#9b1c1c", fontsize=11)
    ax.text(-2.85, 2.05, "Δ = 0 on the curve k² = 4d", fontsize=10.5)
    kk = np.linspace(-2, 3, 10)
    ax.plot(kk, kk - 1, color="#9467bd", lw=1.8, ls="--")
    ax.text(-0.5, -2.35, "Hirzebruch's χ_y genera\n(one line, d = k − 1)", color="#9467bd", fontsize=10)
    pts = [
        (0, 0, "addition", (-0.75, 0.15)),
        (1, 0, "circle operation ∘₁ (Mathlib's\ncircle op); Todd genus", (0.1, 0.12)),
        (-1, 0, "k = −1: 'A or B'; same-target\nRulecraft rules; holographic\ncorner sum", (-1.95, -0.6)),
        (0, -1, "velocity addition (c = 1); signature\ngenus; barrier amplitudes (K = −1)", (0.1, -0.4)),
        (0, 1, "tan addition; beam-splitter\namplitudes (K = +1)", (0.1, 0.12)),
        (2, 1, "Euler characteristic\ngenus (on the parabola)", (-0.2, 0.25)),
    ]
    for x, y, lab, off in pts:
        ax.scatter([x], [y], s=70, color="#111", zorder=5)
        ax.text(x + off[0], y + off[1], lab, fontsize=9.6)
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.axhline(0, color="#ccc", lw=0.8); ax.axvline(0, color="#ccc", lw=0.8)
    ax.set_xlabel("homothety k"); ax.set_ylabel("curvature dial d")
    ax.set_title("The (k, d) plane: every law (a + b + k·ab)/(1 − d·ab) is one point\n"
                 "(the classification is Coleman–McGuinness 1991; the placements are from landscape/)", fontsize=12)
    save(fig, "fig4-kd-plane.png")


# ---------------------------------------------------------------- fig 5: how Rulecraft became a lead
def fig5():
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.8))
    # 1 corner
    ax = axs[0]
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, fc="#eee", ec="k"))
    ax.add_patch(plt.Rectangle((1, 0), 0.5, 1, fc="#cfe2f3", ec="k"))
    ax.add_patch(plt.Rectangle((0, 1), 1, 0.4, fc="#cfe2f3", ec="k"))
    ax.add_patch(plt.Rectangle((1, 1), 0.5, 0.4, fc="#f4b183", ec="k"))
    ax.text(1.25, 1.2, "spillover\n= corner", ha="center", va="center", fontsize=9.5)
    ax.text(0.5, 0.5, "1", ha="center", va="center")
    ax.set_xlim(-0.1, 1.7); ax.set_ylim(-0.1, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("1. Nov 25: three book ideas fit the\nbox (spillover, wisdom, the wall −1).\nTwo were closed: 'order can't matter'.",
                 fontsize=10.5)
    # 2 two orders
    ax = axs[1]
    ax.bar([0, 1], [1600, 1650], color=["#1f77b4", "#d62728"], width=0.6)
    ax.set_ylim(1500, 1680); ax.set_xticks([0, 1]); ax.set_xticklabels(["raise\nthen grant", "grant\nthen raise"])
    ax.set_title("2. Fresh look: a reform is a stretch\nplus a shift. Mixed reforms don't\ncommute (gap = oriented area).", fontsize=10.5)
    # 3 strength and target
    ax = axs[2]
    x = np.linspace(0, 1, 50)
    for s, tg, col in [(0.3, 0.8, "#1f77b4"), (0.5, 0.2, "#d62728"), (0.4, 0.8, "#2ca02c")]:
        ax.plot(x, x + s * (tg - x), color=col, lw=2)
    ax.plot(x, x, "k--", lw=1)
    ax.set_xlabel("slider before"); ax.set_ylabel("after")
    ax.set_title("3. Nov 26: each rule = strength +\ntarget. Same target: fuse at k = −1.\nBest order: sort by target.", fontsize=10.5)
    # 4 loop
    ax = axs[3]
    th = np.linspace(0, 8 * np.pi, 400)
    r = np.exp(-th / 9)
    ax.plot(0.5 + 0.4 * r * np.cos(th), 0.5 + 0.4 * r * np.sin(th), color="#9467bd", lw=1.8)
    ax.scatter([0.5], [0.5], color="k", zorder=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("slider A"); ax.set_ylabel("slider B")
    ax.set_title("4. Rules that read each other: one\nnumber, the loop product b·c, decides\nsettle, spiral, overshoot or run away.", fontsize=10.5)
    fig.suptitle("How Rulecraft became a lead: 'history' means whether the order rules are passed in changes where a "
                 "society ends up", fontsize=13.5, weight="bold", y=1.06)
    save(fig, "fig5-rulecraft-path.png")


# ---------------------------------------------------------------- fig 6: how much of the fused op is used
def fig6():
    from collections import Counter
    fig, ax = plt.subplots(figsize=(12, 3.6))
    bins = ["op", "plane", "pic"]
    cols = {"lean": "#222", "checked": "#777", "partial": "#ddd"}
    left = 0
    for b in bins:
        c = Counter(r[3] for r in RESULTS if r[2] == b)
        start = left
        for s in ["lean", "checked", "partial"]:
            if c[s]:
                ax.barh([0], [c[s]], left=left, color=cols[s], edgecolor="white")
                ax.text(left + c[s] / 2, 0, str(c[s]), ha="center", va="center",
                        color="white" if s != "partial" else "black", fontsize=12, weight="bold")
                left += c[s]
        ax.text((start + left) / 2, 0.55, {"op": "uses ∘k itself", "plane": "plane law /\nwider composition",
                                            "pic": "dial picture only"}[b] + f"\n({left - start})",
                ha="center", va="bottom", fontsize=11, weight="bold")
        ax.axvline(left, color="k", lw=2)
    ax.set_xlim(0, len(RESULTS)); ax.set_ylim(-0.6, 1.5); ax.axis("off")
    ax.text(0, -0.55, "black = proved in Lean    grey = checked    light = partly proved", fontsize=10)
    ax.set_title(f"Is the fused operation really used?  {len(RESULTS)} headline results, sorted", fontsize=13, weight="bold")
    save(fig, "fig6-how-much-is-used.png")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6):
        f()
    print("figures written to", OUT)
