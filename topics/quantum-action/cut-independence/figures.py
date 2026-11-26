"""Figures for the cut-independence page. Run after checks.py:  python3 figures.py  -> figures/*.png"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch, Circle

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
R = json.load(open(os.path.join(HERE, "results.json")))
M = R["_meta"]
BLUE, ORANGE, GREEN, PURPLE, GRAY, RED = "#2a6fdb", "#e8833a", "#2e9e5b", "#8a4fcf", "#888888", "#cc3d3d"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=140, bbox_inches="tight"); plt.close(fig)
def qn(n, r): return math.sin(n * math.pi / r) / math.sin(math.pi / r)

# ---------------------------------------------------------------- fig1: the chain
fig, ax = plt.subplots(figsize=(14, 4.6)); ax.axis("off")
steps = [("Swap rule\nyx = q·xy", "our face weight", BLUE),
         ("Its symmetry\nmust use the\nsame q", "Q2, exact", BLUE),
         ("Gluing rule\n(twisted\nproduct rule)", "Q3, exact", BLUE),
         ("Recoupling\nnumbers =\nq-6j symbols", "Q4, 342 numbers", PURPLE),
         ("Pentagon\nidentity", "Q6", PURPLE),
         ("2–3 and 1–4\nmoves", "Q7, Q8", GREEN),
         ("Same answer\nhowever space\nis cut", "Q9, 4 cuttings", GREEN)]
x0 = 0.02; w = 0.118; gap = 0.022
for i, (t, sub, c) in enumerate(steps):
    x = x0 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, 0.30), w, 0.46, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=c, alpha=0.14, ec=c, lw=2, transform=ax.transAxes))
    ax.text(x + w / 2, 0.53, t, ha="center", va="center", fontsize=11.5, transform=ax.transAxes)
    ax.text(x + w / 2, 0.20, sub, ha="center", va="center", fontsize=10, color=c, transform=ax.transAxes)
    if i < len(steps) - 1:
        ax.annotate("", (x + w + gap - 0.002, 0.53), (x + w + 0.002, 0.53), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=2, color=GRAY))
ax.text(0.5, 0.94, "Every link holds: the swap rule's q forces the whole of 3D quantum gravity's cut-independence",
        ha="center", fontsize=13.5, transform=ax.transAxes, weight="bold")
ax.text(0.5, 0.04, "Needs one extra choice: q a root of unity, q = e^(iπ/r). That is what makes the sum finite "
        "(and it is the positive cosmological constant).", ha="center", fontsize=11, color=GRAY, transform=ax.transAxes)
save(fig, "fig1-the-chain.png")

# ---------------------------------------------------------------- fig2: plane and its symmetry
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw={"width_ratios": [1.25, 1]})
a = ax[0]; a.axis("off")
rows = {0: ["1"], 1: ["1", "1"], 2: ["1", "1+q", "1"], 3: ["1", "1+q+q²", "1+q+q²", "1"],
        4: ["1", "1+q+q²+q³", "1+q+2q²+q³+q⁴", "1+q+q²+q³", "1"]}
for n, items in rows.items():
    for k, t in enumerate(items):
        x = 0.5 + (k - n / 2) * 0.235; y = 0.86 - n * 0.17
        a.text(x, y, t, ha="center", va="center", fontsize=11 if n < 4 else 9.5,
               bbox=dict(boxstyle="round,pad=0.3", fc=BLUE if 0 < k < n else "white", alpha=0.15 if 0 < k < n else 1, ec=BLUE),
               transform=a.transAxes)
a.text(0.5, 0.0, "(x + y)ⁿ when every swap y·x → x·y costs one factor of q.\nSet q = 1 and it is Pascal's triangle.",
       ha="center", fontsize=11, transform=a.transAxes)
a.set_title("The swap rule's plane counts with q-numbers", fontsize=13)
a = ax[1]
p = np.linspace(0.4, 1.9, 200); qv = 1.3
a.plot(p, p - qv, color=RED, lw=2.5, label="what E leaves behind:  (p − q)·x²")
a.plot(p, -(p - qv), color=ORANGE, lw=2.5, ls="--", label="what F leaves behind (sign flipped): (p − q)·y²")
a.axhline(0, color=GRAY, lw=1); a.axvline(qv, color=BLUE, lw=1.5, ls=":")
a.plot([qv], [0], "o", color=BLUE, ms=11)
a.text(qv - 0.06, 0.42, "only p = q\nkeeps yx = q·xy", color=BLUE, fontsize=11, ha="right")
a.set_xlabel("the symmetry's own parameter p   (swap rule fixed at q = 1.3)")
a.set_ylabel("left over after acting on yx − q·xy")
a.legend(loc="lower right", fontsize=9.5, frameon=False)
a.set_title("A symmetry of the plane has no freedom:\nit must carry the swap rule's q", fontsize=13)
save(fig, "fig2-plane-and-its-symmetry.png")

# ---------------------------------------------------------------- fig3: recoupling = q-6j
fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1, 1.05]})
a = ax[0]; a.axis("off")
def tree(a, x0, y0, left=True, lab="j₁₂"):
    leaves = [(x0 - 0.3, y0 + 0.45), (x0, y0 + 0.45), (x0 + 0.3, y0 + 0.45)]
    mid = (x0 - 0.15, y0 + 0.2) if left else (x0 + 0.15, y0 + 0.2)
    root = (x0, y0)
    pairs = [(leaves[0], mid), (leaves[1], mid), (leaves[2], root)] if left else [(leaves[0], root), (leaves[1], mid), (leaves[2], mid)]
    for (u, v) in pairs: a.plot([u[0], v[0]], [u[1], v[1]], color=BLUE, lw=2.5)
    a.plot([mid[0], root[0]], [mid[1], root[1]], color=PURPLE, lw=3)
    a.plot([root[0], root[0]], [root[1], root[1] - 0.12], color=BLUE, lw=2.5)
    for i, u in enumerate(leaves): a.text(u[0], u[1] + 0.05, f"j{'₁₂₃'[i]}", ha="center", fontsize=12)
    a.text(mid[0] + (-0.07 if left else 0.07), mid[1] - 0.06, lab, color=PURPLE, fontsize=12, ha="center")
    a.text(root[0], root[1] - 0.18, "J", ha="center", fontsize=12)
tree(a, 0.25, 0.25, True, "j₁₂"); tree(a, 1.05, 0.25, False, "j₂₃")
a.annotate("", (0.78, 0.48), (0.52, 0.48), arrowprops=dict(arrowstyle="<->", lw=2, color=GRAY))
a.text(0.65, 0.53, "recouple", ha="center", color=GRAY)
a.set_xlim(-0.15, 1.45); a.set_ylim(-0.05, 0.85)
a.text(0.65, -0.02, "Glue three copies of the plane two ways: first 1 with 2, or first 2 with 3.\n"
       "The overlap between the two ways is the recoupling number.", ha="center", fontsize=11)
a.set_title("Gluing three pieces, two ways", fontsize=13)
a = ax[1]
pairs = np.array(M["Q4pairs"])
a.plot([0, 1], [0, 1], color=GRAY, lw=1)
a.scatter(pairs[:, 1], pairs[:, 0], s=22, color=PURPLE, alpha=0.75, zorder=3)
a.set_xlabel("quantum gravity's q-6j:  √([2j₁₂+1][2j₂₃+1]) · |6j|_q")
a.set_ylabel("overlap built from the plane and its gluing rule")
a.set_title(f"All {len(pairs)} recoupling numbers land on the line\n(q = 1.3; worst difference "
            f"{float(R['Q4 recoupling = q-6j']['detail'].split('to ')[-1]):.0e})", fontsize=13)
a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.set_aspect("equal")
save(fig, "fig3-recoupling-is-the-q6j.png")

# ---------------------------------------------------------------- fig4: Pachner moves
fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))
def P(*pts): return np.array(pts)
a = ax[0]; a.axis("off")
top, bot = P(0.0, 1.0), P(0.0, -1.0); tri = [P(-0.8, 0.05), P(0.75, 0.2), P(0.1, -0.25)]
for dx, three in ((-1.1, False), (1.1, True)):
    T = top + [dx, 0]; B = bot + [dx, 0]; tr = [t + [dx, 0] for t in tri]
    for t in tr:
        a.plot(*zip(T, t), color=BLUE, lw=2); a.plot(*zip(B, t), color=BLUE, lw=2)
    for i in range(3):
        a.plot(*zip(tr[i], tr[(i + 1) % 3]), color=BLUE, lw=2, ls="-" if not three else ":")
    if not three:
        a.add_patch(Polygon(tr, closed=True, fc=ORANGE, alpha=0.3, ec="none"))
        a.text(dx, -1.35, "2 tetrahedra\nsharing a triangle", ha="center")
    else:
        a.plot(*zip(T, B), color=GREEN, lw=3)
        a.text(dx, -1.35, "3 tetrahedra\naround a new edge", ha="center")
a.text(0, 0.0, "=", fontsize=30, ha="center", va="center")
a.text(1.1, 1.18, "sum over the new edge's spin", ha="center", color=GREEN, fontsize=10.5)
a.set_xlim(-2.1, 2.1); a.set_ylim(-1.7, 1.4); a.set_aspect("equal")
w23 = max(float(v["detail"].split()[0]) if False else 0 for v in [R["Q7 Pachner 2-3 move"]])
a.set_title("2–3 move: holds to 1e-15 at every tested q = e^(iπ/r)\n(it is exactly the pentagon identity)", fontsize=12.5)
a = ax[1]; a.axis("off")
V = [P(-0.9, -0.6), P(0.9, -0.6), P(0.0, 0.95), P(0.15, -0.15)]
for dx, inner in ((-1.15, False), (1.15, True)):
    W = [v + [dx, 0] for v in V]
    for i in range(4):
        for j in range(i + 1, 4):
            a.plot(*zip(W[i], W[j]), color=BLUE, lw=2, ls=":" if (i, j) in ((0, 1),) else "-")
    if inner:
        c = sum(W) / 4 + [0, 0.02]
        for v in W: a.plot(*zip(c, v), color=GREEN, lw=2.5)
        a.plot(*c, "o", color=GREEN, ms=9)
        a.text(dx, -1.05, "4 tetrahedra around a new point\n÷ w² for the new point", ha="center")
    else:
        a.text(dx, -1.05, "1 tetrahedron", ha="center")
a.text(0, 0.1, "=", fontsize=30, ha="center", va="center")
a.set_xlim(-2.2, 2.2); a.set_ylim(-1.35, 1.25); a.set_aspect("equal")
a.set_title("1–4 move: holds when each new point is weighted 1/w²,\nw² = Σ[2j+1]² = r / (2 sin²(π/r))", fontsize=12.5)
save(fig, "fig4-pachner-moves.png")

# ---------------------------------------------------------------- fig5: four cuttings, one answer
fig, ax = plt.subplots(figsize=(12, 5))
names = ["2 tetrahedra\n(4 points)", "5 tetrahedra\n(boundary of a 4-simplex)", "8 tetrahedra\n(after a 1–4 move)", "9 tetrahedra\n(then a 2–3 move)"]
cols = [BLUE, PURPLE, ORANGE, GREEN]
rs = sorted(M["Q9"], key=int)
for i, rr in enumerate(rs):
    zs = M["Q9"][rr]
    for k, z in enumerate(zs):
        ax.bar(i * 5 + k, z, color=cols[k], label=names[k] if i == 0 else None)
    ax.plot([i * 5 - 0.6, i * 5 + 3.6], [1 / (int(rr) / (2 * math.sin(math.pi / int(rr))**2))] * 2, color="k", lw=1.5, ls="--")
    ax.text(i * 5 + 1.5, zs[0] + 0.02, f"r = {rr}\nZ = {zs[0]:.6f}", ha="center", fontsize=10.5)
ax.set_xticks([]); ax.set_ylabel("the whole state sum Z for the 3-sphere")
ax.set_ylim(0, 0.62)
ax.legend(frameon=False, fontsize=9.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.03))
ax.set_title("Four different ways to cut the 3-sphere into tetrahedra give the same number to 12 digits\n"
             "(dashed: 1/w², the known Turaev–Viro value)", fontsize=12.5)
save(fig, "fig5-four-cuttings-one-answer.png")

# ---------------------------------------------------------------- fig6: why a root of unity
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
a = ax[0]; r = 6
kk = 2 * math.cos(math.pi / (2 * r))
u = [0.0, 1.0]
for n in range(1, 4 * r + 1): u.append(kk * u[n] - u[n - 1])
n = np.arange(0, 4 * r + 1); u = np.array(u[:4 * r + 1])
a.axhline(0, color=GRAY, lw=1)
a.bar(n, u, color=[(GREEN if 0 < x < 2 * r else RED) if x % 2 == 0 else "#cfd6e6" for x in n])
a.plot(n, u, color=BLUE, lw=1)
a.text(2 * r + 0.4, 0.5, f"step 2r = {2*r}: zero", ha="left", color=RED)
a.set_xlabel("step n of the dial sequence  u(n+1) = k·u(n) − u(n−1)")
a.set_ylabel("u(n)")
a.set_title(f"Our dial's own sequence at q = e^(iπ/{r}), k = q^½ + q^-½.\nEven steps (bold) are the quantum dimensions [2j+1]; "
            f"they stop at 2j+1 = r", fontsize=12)
a = ax[1]
rr = np.arange(3, 300)
w2 = np.array([sum(qn(k + 1, x)**2 for k in range(x - 1)) for x in rr])
a.loglog(rr, w2, color=PURPLE, lw=2.5, label="w² = Σ [2j+1]²  (vertex weight)")
a.loglog(rr, rr**3 / (2 * math.pi**2), color=GRAY, lw=1.5, ls="--", label="r³ / (2π²)")
a.set_xlabel("r   (q = e^(iπ/r); q → 1 as r → ∞)"); a.set_ylabel("w²")
a.legend(frameon=False)
a.set_title("At q = 1 the vertex weight is infinite\n(Ponzano–Regge needs it removed by hand)", fontsize=12)
save(fig, "fig6-why-a-root-of-unity.png")

# ---------------------------------------------------------------- fig7: the dial picture of Lambda
fig, ax = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1.35, 1]})
a = ax[0]
k = np.linspace(0, 2.6, 400)
dp = 1 - k**2 / 4
a.plot(k, dp, color=BLUE, lw=2.5)
a.axhline(0, color=GRAY, lw=1)
a.fill_between(k, dp, 0, where=dp > 0, color=GREEN, alpha=0.15)
a.fill_between(k, dp, 0, where=dp < 0, color=RED, alpha=0.12)
for r_ in range(2, 11):
    kk = 2 * math.cos(math.pi / (2 * r_)); a.plot(kk, 1 - kk**2 / 4, "o", color=GREEN, ms=6)
    if r_ in (2, 3, 4): a.text(kk - 0.03, 1 - kk**2 / 4 + 0.05, f"r={r_}", ha="center", fontsize=9, color=GREEN)
a.plot(2, 0, "o", color="k", ms=8); a.text(2.03, 0.05, "q = 1: straight\n(k² = 4d), Λ = 0\nPonzano–Regge", fontsize=9.5)
a.text(0.2, 0.3, "q = e^(iπ/r): angle side\n(k² < 4d), Λ > 0\nd' = sin²(π/2r)", color=GREEN, fontsize=11)
a.text(1.0, -0.5, "q real: radius side\n(k² > 4d), Λ < 0", color=RED, fontsize=11)
a.set_xlabel("k = q^½ + q^-½   (with d = 1; face weight q = β/α)"); a.set_ylabel("d' = d − k²/4   (after the ruler change)")
a.set_title("The sign of d − k²/4 sorts 3D quantum gravity\nby the sign of its cosmological constant", fontsize=12.5)
a = ax[1]; a.set_aspect("equal"); a.axis("off")
r_ = 4; N = 4 * r_
th = np.linspace(0, 2 * math.pi, 300); a.plot(np.cos(th), np.sin(th), color=GRAY, lw=1)
pts = [(math.cos(i * math.pi / (2 * r_)), math.sin(i * math.pi / (2 * r_))) for i in range(N)]
for i, (x, y) in enumerate(pts):
    a.plot(x, y, "o", color=GREEN if i < 2 * r_ else ORANGE, ms=8)
for i in range(N):
    a.annotate("", pts[(i + 1) % N], pts[i], arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))
a.text(0, 0.08, f"each step turns by π/(2r) = π/{2*r_}", ha="center", fontsize=11)
a.text(0, -0.15, f"after 2r = {2*r_} steps: −1;  after {N}: back home", ha="center", fontsize=10.5, color=GRAY)
a.set_title(f"The dial's step matrix [[k, −1], [1, 0]]\ncloses up exactly at a root of unity (r = {r_})", fontsize=12.5)
save(fig, "fig7-the-dial-sorts-lambda.png")

# ---------------------------------------------------------------- fig8: scorecard
fig, ax = plt.subplots(figsize=(13, 5.2)); ax.axis("off")
rows = [("Swap rule's q weight = Manin's quantum plane", "the identification is ours (earlier pass)", BLUE),
        ("Plane's symmetry is U_q(sl2), same q", "known: Manin 1988", GRAY),
        ("Gluing → q-6j; pentagon; orthogonality", "known: Kirillov–Reshetikhin 1989", GRAY),
        ("2–3 and 1–4 moves; same Z for every cutting", "known: Turaev–Viro 1992", GRAY),
        ("Root of unity ↔ Λ > 0 ; real q ↔ Λ < 0", "known (3D gravity literature)", GRAY),
        ("All of it rebuilt from our swap rule, checked numerically", "ours: this page (11 of 11)", GREEN),
        ("Angle / straight / radius sides of our dial = Λ > 0 / = 0 / < 0", "ours: a dictionary entry, not a new result", ORANGE)]
for i, (t, s, c) in enumerate(rows):
    y = 0.9 - i * 0.125
    ax.text(0.02, y, t, fontsize=12, va="center")
    ax.text(0.62, y, s, fontsize=11.5, va="center", color=c, weight="bold" if c != GRAY else "normal")
ax.set_title("What is known and what is ours", fontsize=14)
save(fig, "fig8-scorecard.png")
print("figures written:", sorted(os.listdir(FIG)))
