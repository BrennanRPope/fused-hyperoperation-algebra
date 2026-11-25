"""Figures for cube-weights (DRAFT). Run: python3 figures.py"""
import itertools, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch

COL = {0: "#d1495b", 1: "#2e86ab", 2: "#3caa5a"}   # direction colors
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

def cube_xy(S, origin=(0, 0), scale=1.0):
    v = {0: np.array([1.0, 0.0]), 1: np.array([0.0, 1.0]), 2: np.array([0.45, 0.35])}
    p = np.array(origin, float)
    for i in S: p = p + scale * v[i]
    return p

def draw_cube(ax, origin=(0, 0), scale=2.0, corner_label=None, edge_label=None, edge_lw=2.2,
              dot=60, fs=10, highlight_face=None, path=None, lab_off=None):
    corners = [frozenset(S) for r in range(4) for S in itertools.combinations(range(3), r)]
    if highlight_face is not None:
        S, i, j = highlight_face
        pts = [cube_xy(S, origin, scale), cube_xy(S | {i}, origin, scale),
               cube_xy(S | {i, j}, origin, scale), cube_xy(S | {j}, origin, scale)]
        ax.add_patch(Polygon(pts, closed=True, color="#f6c85f", alpha=0.45, zorder=0))
    for S in corners:
        for i in range(3):
            if i in S: continue
            p, q = cube_xy(S, origin, scale), cube_xy(S | {i}, origin, scale)
            ax.plot([p[0], q[0]], [p[1], q[1]], color=COL[i], lw=edge_lw, zorder=1)
            if edge_label:
                t = edge_label(S, i)
                if t:
                    m = (p + q) / 2
                    ax.text(m[0], m[1], t, fontsize=fs - 1, color=COL[i], ha="center", va="center",
                            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9), zorder=3)
    if path:
        S = frozenset(); pts = [cube_xy(S, origin, scale)]
        for i in path:
            S = S | {i}; pts.append(cube_xy(S, origin, scale))
        pts = np.array(pts)
        ax.plot(pts[:, 0], pts[:, 1], color="black", lw=5, alpha=0.25, zorder=0.5, solid_capstyle="round")
    for S in corners:
        p = cube_xy(S, origin, scale)
        ax.scatter([p[0]], [p[1]], s=dot, color="black", zorder=4)
        if corner_label:
            t = corner_label(S)
            off = lab_off(S) if lab_off else (-0.12 if (len(S) and 2 not in S and 1 not in S) else 0.12, 0.0)
            ax.text(p[0] + off[0], p[1] + off[1] + 0.14, t, fontsize=fs, ha="center", va="bottom", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.15", fc="#fffdf5", ec="#999", lw=0.6))

def name(S):
    return "∅" if not S else "{" + ",".join(str(i + 1) for i in sorted(S)) + "}"

def save(fig, fn):
    fig.savefig(fn, dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", fn)

# ---------------------------------------------------------------- fig1 parts
fig, ax = plt.subplots(figsize=(8.5, 6.2))
draw_cube(ax, scale=3.0, corner_label=name, highlight_face=(frozenset(), 0, 1), path=[2, 0, 1], fs=11)
o = lambda S: cube_xy(frozenset(S), (0, 0), 3.0)
kw = dict(fontsize=11, arrowprops=dict(arrowstyle="->", color="#444"))
ax.annotate("corner (a subset)\nweighted already: k^(m-1),\nLucas U_m, soliton ∏A_ij", xy=o({0, 1, 2}), xytext=(5.0, 5.2), **kw)
ax.annotate("edge = one step:\nadd one input", xy=(o({0}) + o({0, 2})) / 2, xytext=(4.3, 0.2), **kw)
ax.annotate("square face:\ntwo inputs meet", xy=(1.5, 1.3), xytext=(-2.6, 2.0), **kw)
ax.annotate("path ∅ → full = one ordering\n(here 3, then 1, then 2)", xy=o({2}), xytext=(-2.8, 4.6), **kw)
for i, lbl in enumerate(["direction 1", "direction 2", "direction 3"]):
    ax.plot([], [], color=COL[i], lw=3, label=lbl)
ax.legend(loc="lower left", frameon=False)
ax.set_title("The parts of the cube that could carry a weight (n = 3)", fontsize=13)
ax.set_xlim(-3, 7.4); ax.set_ylim(-0.6, 6.4); ax.axis("off"); ax.set_aspect("equal")
save(fig, "fig1-cube-parts.png")

# ---------------------------------------------------------------- fig2 map
rows = [
    ("Corners", "Yes: k^(m-1), Lucas U_m(k,d),\nsoliton ∏A_ij, principal minors", "Boolean lattice, e_m (Newton–Girard)", "Done: only k^(m-1) composes\nas one law (check A2)", "used"),
    ("Edges", "Yes: k is the edge weight;\none dial per element", "discrete connection / link variable", "No: a flat edge field = corner\nweights (B4); curved = order dial", "dead"),
    ("Squares\n(interaction)", "Yes: soliton pair weight A_ij", "Hirota phase shift; 2-body term", "Not on numbers: fails\nassociativity (A3)", "dead"),
    ("Squares\n(reordering cost)", "Partly: order dial's one q", "quantum plane, many q_ij\n(Manin; Artin–Schelter–Tate)", "LEAD: many q_ij make soliton\ncorners from mult-add (C4–C7)", "lead"),
    ("Higher faces", "Via principal minors\n(sibling stream)", "Möbius / Harsanyi dividends,\ninteraction information", "Not here: product weights\nfactor (D1); minors = sibling", "sib"),
    ("Paths", "Yes: permutohedron", "inversions, q-factorial\n(MacMahon)", "Same q as above, no new dial", "used"),
    ("Orientations", "No", "Koszul complex; N-complexes\n(Mayer 1942, Kapranov 1996)", "Known: holonomy q gives d^N = 0\n(F2); a home, not a dial", "known"),
    ("Möbius function", "Yes at k = −1", "Rota 1964; Hall; Fubini\npolynomials", "No new dial: μ_k = permutohedron\nfaces at −k (G1–G3)", "used"),
    ("Dial plane d", "Yes (Lucas weights)", "Cayley–Hamilton", "Reading: d folds each square\nback two ranks (H1)", "used"),
    ("q-orbit", "Yes (Jacobi triple product)", "Cauchy q-binomial theorem", "Ties to squares: geometric\ndials = pair weight q (I1)", "used"),
]
color = {"used": "#d9ead3", "dead": "#f4cccc", "lead": "#ffe599", "sib": "#cfe2f3", "known": "#ead1dc"}
fig, ax = plt.subplots(figsize=(14, 9.2)); ax.axis("off")
heads = ["Part of the cube", "Already in the work?", "Existing name (literature)", "New dial?"]
xs = [0.0, 0.17, 0.43, 0.70, 1.0]
for c, h in enumerate(heads):
    ax.text((xs[c] + xs[c + 1]) / 2, 1.0, h, ha="center", va="bottom", fontsize=12, fontweight="bold")
H = 0.095
for r, row in enumerate(rows):
    y = 0.985 - (r + 1) * H
    for c in range(4):
        fc = color[row[4]] if c == 3 else ("#f3f3f3" if c == 0 else "white")
        ax.add_patch(FancyBboxPatch((xs[c] + 0.003, y + 0.004), xs[c + 1] - xs[c] - 0.006, H - 0.008,
                                    boxstyle="round,pad=0.003", fc=fc, ec="#bbb", lw=0.7))
        ax.text((xs[c] + xs[c + 1]) / 2, y + H / 2, row[c], ha="center", va="center", fontsize=10,
                fontweight="bold" if c == 0 else "normal")
leg = [("used", "in the work / a reading"), ("lead", "lead worth following"), ("dead", "dead end"),
       ("sib", "belongs to the sibling stream"), ("known", "known elsewhere, not yet used")]
for i, (kk, t) in enumerate(leg):
    ax.add_patch(FancyBboxPatch((0.0 + i * 0.205, -0.04), 0.02, 0.025, boxstyle="round,pad=0.002", fc=color[kk], ec="#999"))
    ax.text(0.025 + i * 0.205, -0.028, t, va="center", fontsize=9.5)
ax.set_xlim(0, 1); ax.set_ylim(-0.06, 1.04)
ax.set_title("Map: every part of the cube that can carry a weight (DRAFT)", fontsize=14)
save(fig, "fig2-map.png")

# ---------------------------------------------------------------- fig3 edges
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
draw_cube(ax, scale=3.0, edge_label=lambda S, i: f"k·a{i + 1}",
          corner_label=lambda S: "1" if not S else (f"k^{len(S)}·" if len(S) > 1 else "k·") + "".join(f"a{i + 1}" for i in sorted(S)), fs=10)
ax.set_title("The dial k is a weight on every edge.\nA corner multiplies the edges on any route to it.", fontsize=11.5)
ax.axis("off"); ax.set_aspect("equal"); ax.set_xlim(-1, 5); ax.set_ylim(-0.6, 5.3)
ax = axs[1]; ax.axis("off"); ax.set_aspect("equal")
def square(ax, ox, title, top_l, top_r, col):
    P = [np.array([ox, 0]), np.array([ox + 2, 0]), np.array([ox + 2, 2]), np.array([ox, 2])]
    ax.add_patch(Polygon(P, closed=True, fc=col, alpha=0.25, ec="none"))
    for a_, b_, c_ in [(0, 1, COL[0]), (1, 2, COL[1]), (0, 3, COL[1]), (3, 2, COL[0])]:
        ax.annotate("", xy=P[b_], xytext=P[a_], arrowprops=dict(arrowstyle="-|>", color=c_, lw=2.2))
    ax.text(ox + 1, -0.35, "right then up", ha="center", fontsize=9.5, color="#555")
    ax.text(ox - 0.15, 1, "up then\nright", ha="right", va="center", fontsize=9.5, color="#555")
    ax.text(ox + 1, 2.25, top_l, ha="center", fontsize=10)
    ax.text(ox + 1, 2.95, title, ha="center", fontsize=11.5, fontweight="bold")
    ax.text(ox + 1, 1.0, top_r, ha="center", va="center", fontsize=10.5)
square(ax, 0.6, "Flat square", "routes agree", "closes:\nholonomy 1", "#3caa5a")
square(ax, 4.0, "Curved square", "routes differ by a·b·[P, Q]", "misses:\ncurvature", "#d1495b")
ax.text(3.3, -1.25, "Flat edge weights carry exactly the corner weights (nothing new).\n"
        "Curvature lives on squares: with matrices it is the order dial's commutator.", ha="center", fontsize=10.5)
ax.set_xlim(-0.6, 6.8); ax.set_ylim(-1.7, 3.4)
save(fig, "fig3-edges-and-curvature.png")

# ---------------------------------------------------------------- fig4 associativity
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]
m = np.arange(1, 7); kv = 1.3
ax.semilogy(m, kv**(m - 1), "o-", color="#2e86ab", lw=2, label="chain weight k^(m−1)  (one edge factor per step)")
ax.semilogy(m, kv**(m * (m - 1) / 2), "s-", color="#d1495b", lw=2, label="pair weight A^C(m,2)  (one factor per pair)")
ax.set_xlabel("corner rank m (how many inputs meet)"); ax.set_ylabel("corner weight (k = A = 1.3)")
ax.set_title("Chain weights grow like a line in m; pair weights like m²", fontsize=11.5)
ax.legend(frameon=False, fontsize=9.5); ax.grid(alpha=0.3)
ax = axs[1]
Av = np.linspace(-0.6, 1.6, 400)
ax.plot(Av, Av**2 * (1 - Av), color="#d1495b", lw=2.4, label="pair weight: defect = A²(1 − A)·xyz")
ax.axhline(0, color="#2e86ab", lw=2.4, label="chain weight k^(m−1): defect 0 for every k")
ax.scatter([0, 1], [0, 0], s=70, color="black", zorder=5)
ax.annotate("A = 0: plain addition", xy=(0, 0), xytext=(-0.55, 0.32), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
ax.annotate("A = 1: mult-add", xy=(1, 0), xytext=(1.05, 0.32), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("constant pair weight A"); ax.set_ylabel("(x∘y)∘z  −  three-input corner sum")
ax.set_title("A pair weight on plain numbers is not one associative law", fontsize=11.5)
ax.legend(frameon=False, fontsize=9.5, loc="lower left"); ax.grid(alpha=0.3)
save(fig, "fig4-which-weights-compose.png")

# ---------------------------------------------------------------- fig5 reverse
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
def lbl_fwd(S):
    return "1" if not S else "".join(f"a{i + 1}" for i in sorted(S))
def lbl_rev(S):
    if not S: return "1"
    pr = "".join(f"q{i + 1}{j + 1}·" for i, j in itertools.combinations(sorted(S), 2))
    return pr + "".join(f"a{i + 1}" for i in sorted(S))
for ax, f, t in [(axs[0], lbl_fwd, "Order 1, 2, 3 (already sorted):\nno pair weights"),
                 (axs[1], lbl_rev, "Order 3, 2, 1 (every pair reversed):\neach corner gets all its pair weights")]:
    draw_cube(ax, scale=3.2, corner_label=f, fs=9.5)
    ax.set_title(t, fontsize=11.5); ax.axis("off"); ax.set_aspect("equal"); ax.set_xlim(-1.3, 5.6); ax.set_ylim(-0.6, 5.6)
fig.suptitle("Inputs that commute up to a factor (a_j a_i = q_ij · a_i a_j), written in sorted order, k = 1.\n"
             "Put q_ij = ((p_i − p_j)/(p_i + p_j))² and the right cube is the KdV 3-soliton (check C7).", fontsize=11.5, y=1.03)
save(fig, "fig5-reorder-makes-pair-weights.png")

# ---------------------------------------------------------------- fig6 paths and Koszul
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axs[0]; ax.axis("off"); ax.set_aspect("equal")
perms = [(1, 2, 3), (2, 1, 3), (2, 3, 1), (3, 2, 1), (3, 1, 2), (1, 3, 2)]
def inv(p): return sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j])
ang = np.pi / 2 + np.arange(6) * np.pi / 3
P = np.c_[np.cos(ang), np.sin(ang)] * 2
ax.add_patch(Polygon(P, closed=True, fc="#cfe2f3", ec="#2e86ab", lw=2))
for (x_, y_), p in zip(P, perms):
    e = inv(p); w = "1" if e == 0 else ("q" if e == 1 else f"q^{e}")
    ax.scatter([x_], [y_], s=60, color="black", zorder=3)
    ax.text(x_ * 1.28, y_ * 1.28, f"{''.join(map(str, p))}\n{w}", ha="center", va="center", fontsize=11)
ax.text(0, 0, "1 + 2q + 2q² + q³\n= (1)(1+q)(1+q+q²)\n= [3]_q!", ha="center", va="center", fontsize=11)
ax.set_title("Each route through the 3-cube is an ordering (a corner of\nthe hexagon). Weight it by q^(out-of-order pairs).", fontsize=11)
ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.1, 3.1)
ax = axs[1]
import importlib.util
def koszul(n, qq, xv):
    subs = [S for r in range(n + 1) for S in itertools.combinations(range(n), r)]
    idx = {S: r for r, S in enumerate(subs)}
    D = np.zeros((2**n, 2**n), dtype=complex)
    for S in subs:
        for i in S:
            T = tuple(s for s in S if s != i)
            D[idx[T], idx[S]] += qq**sum(1 for j in S if j < i) * xv[i]
    return D
n = 6; xv = np.linspace(0.7, 1.3, n)
grid = np.zeros((5, 6))
for a_, Nn in enumerate(range(2, 7)):
    Dq = koszul(n, np.exp(2j * np.pi / Nn), xv)
    for b_, mm in enumerate(range(1, 7)):
        grid[a_, b_] = np.abs(np.linalg.matrix_power(Dq, mm)).max()
im = ax.imshow(grid > 1e-9, cmap=matplotlib.colors.ListedColormap(["#d9ead3", "#f4cccc"]), aspect="auto")
for a_ in range(5):
    for b_ in range(6):
        ax.text(b_, a_, "0" if grid[a_, b_] < 1e-9 else "≠ 0", ha="center", va="center", fontsize=11)
ax.set_xticks(range(6)); ax.set_xticklabels([f"d^{m}" for m in range(1, 7)])
ax.set_yticks(range(5)); ax.set_yticklabels([f"q = e^(2πi/{N})" for N in range(2, 7)])
ax.set_title("Orientation: each square's edges disagree by q.\nThen d^N = 0 exactly when q is an N-th root of 1 (6-cube)", fontsize=11)
save(fig, "fig6-paths-and-orientation.png")

# ---------------------------------------------------------------- fig7 mobius
def S2(n_, m_): return sum((-1)**j * math.comb(m_, j) * (m_ - j)**n_ for j in range(m_ + 1)) // math.factorial(m_)
def mu(r, kk): return sum((-kk)**m_ * math.factorial(m_) * S2(r, m_) for m_ in range(r + 1))
fig, ax = plt.subplots(figsize=(9, 5.4))
kk = np.linspace(-1.2, 1.6, 400)
for r, c in zip(range(1, 6), ["#999", "#2e86ab", "#3caa5a", "#d1495b", "#8e44ad"]):
    ax.plot(kk, [mu(r, x) for x in kk], color=c, lw=2, label=f"n = {r}")
ax.set_yscale("symlog", linthresh=1)
ax.axvline(1, color="black", lw=0.8, ls=":"); ax.axvline(-1, color="black", lw=0.8, ls=":")
ax.text(1.03, 300, "k = 1: ±1\n(Euler: Rota's μ)", fontsize=10)
ax.text(-1.17, 300, "k = −1: all faces\n1, 3, 13, 75, 541", fontsize=10)
ax.set_xlabel("dial k"); ax.set_ylabel("μ_k(∅, full set)  (symlog scale)")
ax.set_title("The Möbius function with the dial: the permutohedron's face counts read at −k", fontsize=12)
ax.legend(frameon=False, ncol=5, loc="lower center"); ax.grid(alpha=0.3)
save(fig, "fig7-mobius-with-dial.png")

# ---------------------------------------------------------------- fig8 q-orbit
fig, ax = plt.subplots(figsize=(10.5, 5.4))
dial = {0: "1", 1: "q", 2: "q²"}
draw_cube(ax, scale=3.0, edge_label=lambda S, i: f"×{dial[i]}",
          corner_label=lambda S: {0: "1", 1: "q", 2: "q²", 3: "q³"}[sum(S)] if S else "1", fs=10)
txt = ("Direction dials in a q-orbit: 1, q, q²\n\nAdd each layer:\n"
       "  0 inputs:  1\n  1 input:   1 + q + q²\n  2 inputs:  q·(1 + q + q²)\n  3 inputs:  q³\n\n"
       "= q^C(m,2) · [3 choose m]_q\n\nA constant pair weight q per pair\nappears from dials on directions.")
ax.text(5.4, 2.6, txt, fontsize=11, va="center", family="DejaVu Sans")
ax.set_title("The q-orbit (solitons page, Jacobi triple product) seen on the cube", fontsize=12)
ax.axis("off"); ax.set_aspect("equal"); ax.set_xlim(-0.8, 10.5); ax.set_ylim(-0.6, 5.2)
save(fig, "fig8-q-orbit.png")
