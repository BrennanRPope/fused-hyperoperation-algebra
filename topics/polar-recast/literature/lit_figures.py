"""Pictures for the literature pass and the super-log question (DRAFT, Nov 26 2026).
Colour key: orange = known math, blue = ours (bridge or reading), red = gap / not rigid."""
import itertools, cmath, json, random, textwrap
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
OR, BL, GR, RD, GY = "#e07b00", "#1f6fd1", "#2a9d55", "#c0392b", "#888888"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
def save(fig, name): fig.tight_layout(); fig.savefig(f"figures/{name}.png", dpi=130); plt.close(fig)

# ---------- fig14: the plane dial is Hirzebruch's two-parameter genus; q = ratio of roots
fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
ax = axs[0]
K = np.linspace(-3, 3, 400); ax.plot(K, K**2/4, color=GR, lw=2.5, label="k² = 4d: straight (Euler characteristic, q = 1)")
ax.fill_between(K, K**2/4, 3, color=OR, alpha=0.12); ax.fill_between(K, -1.5, K**2/4, color=BL, alpha=0.10)
ax.text(0, 2.3, "angle: k² < 4d\nq on the unit circle", ha="center", color=OR, fontweight="bold")
ax.text(1.9, -0.9, "radius: k² > 4d\nq real", ha="center", color=BL, fontweight="bold")
pts = [(1, 0, "Todd genus\n(our core law)", "q = 0"), (0, 1, "signature", "q = −1"), (2, 1, "Euler char.", "q = 1"), (3, 1, "", "q ≈ 0.15"), (1, 2, "", "q on circle")]
for kk, dd, lab, ql in pts:
    ax.plot(kk, dd, "ko"); ax.annotate(f"{lab}\n{ql}".strip(), (kk, dd), xytext=(6, 6), textcoords="offset points", fontsize=8)
ax.set_xlim(-3, 3.6); ax.set_ylim(-1.5, 3); ax.set_xlabel("k (corner dial)"); ax.set_ylabel("d (plane dial)")
ax.legend(fontsize=7, loc="lower left"); ax.set_title("Plane dial = Hirzebruch's two-parameter genus\n(x+y+kxy)/(1−dxy); q = β/α, roots of t²−kt+d")
ax = axs[1]
ax.add_patch(Circle((0, 0), 1, fill=False, color=OR, lw=2)); ax.axhline(0, color=GY, lw=0.6); ax.axvline(0, color=GY, lw=0.6)
for kk, dd, c in [(0, 1, OR), (1, 2, OR), (0.6, 1, OR), (3, 1, BL), (2.5, 1, BL), (2, 1, GR), (1, 0, BL)]:
    disc = cmath.sqrt(kk*kk-4*dd); a, b = (kk+disc)/2, (kk-disc)/2; q = b/a
    ax.plot(q.real, q.imag, "o", color=c, ms=8)
    if c != BL: ax.annotate(f"k={kk}, d={dd}", (q.real, q.imag), xytext=(5, 5), textcoords="offset points", fontsize=7)
ax.annotate("radius dials (k=1,d=0 gives q=0;\nk=3,d=1 and k=2.5,d=1 just right of it)", (0.15, 0), xytext=(-0.9, 0.55), fontsize=7, color=BL, arrowprops=dict(arrowstyle="->", color=BL))
ax.set_aspect("equal"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.3, 1.3)
ax.set_title("Where q lands: angle dials on the circle\n(a turn, as in the swap rule), radius dials on the line")
ax = axs[2]
# Gr(2,4) staircase paths in a 2x2 box with q^area
paths = list(itertools.combinations(range(4), 2))
for idx, S in enumerate(paths):
    ox, oy = (idx % 3)*3.0, -(idx//3)*3.2
    x0, y0, pts_ = 0, 0, [(0, 0)]; area = 0
    for step in range(4):
        if step in S: y0 += 1
        else: x0 += 1; area += y0
        pts_.append((x0, y0))
    xs, ys = zip(*pts_)
    ax.fill_between([ox + p for p in xs], [oy + p for p in ys], oy + 2, step=None, color=BL, alpha=0.0)
    for i in range(2):
        for j in range(2): ax.add_patch(plt.Rectangle((ox+i, oy+j), 1, 1, fill=False, color=GY, lw=0.5))
    ax.plot([ox + p for p in xs], [oy + p for p in ys], color=BL, lw=2.5)
    ax.text(ox+1, oy-0.5, f"q^{area}", ha="center", fontsize=9, color=BL)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.5, 8.5); ax.set_ylim(-4.2, 2.4)
ax.set_title("Face weights = genus of the Grassmannian Gr(2,4)\n1 + q + 2q² + q³ + q⁴ (checked by fixed points)")
save(fig, "fig14-plane-dial-is-a-genus")

# ---------- fig15: every 5-region rule splits into upward corner intervals from triples
tf = json.load(open("literature/tf_check.json"))
ex = tf["examples"]["cyclic"]
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.6, 1]})
ax = axs[0]
REG = "ABCDE"; pos = {}
for g in (3, 4, 5):
    subs = ["".join(c) for c in itertools.combinations(REG, g)]
    for i, sname in enumerate(subs): pos[sname] = ((i - (len(subs)-1)/2)*1.1, g)
for a in pos:
    for b in pos:
        if len(b) == len(a)+1 and set(a) <= set(b): ax.plot(*zip(pos[a], pos[b]), color="#dddddd", lw=0.6, zorder=1)
cols = [BL, OR, GR, RD, "#8e44ad", "#16a085"]
in_rule = set(ex["corners"])
for sname, (px, py) in pos.items():
    ax.add_patch(Circle((px, py), 0.28, color="#f3f3f3" if sname not in in_rule else "white", ec=GY, lw=0.8, zorder=2))
    ax.text(px, py, sname, ha="center", va="center", fontsize=6.5, zorder=4)
for ci, (T, X) in enumerate(ex["tiling"]):
    members = ["".join(sorted(T + "".join(Y))) for kk in range(len(X)+1) for Y in itertools.combinations(X, kk)]
    for m in members:
        px, py = pos[m]; ax.add_patch(Circle((px, py), 0.28, color=cols[ci % 6], alpha=0.35, zorder=3))
    for a in members:
        for b in members:
            if len(b) == len(a)+1 and set(a) <= set(b): ax.plot(*zip(pos[a], pos[b]), color=cols[ci % 6], lw=2.5, zorder=2)
    lab = f"−I3({T[0]}:{T[1]}:{T[2]}" + (f" | {X})" if X else ")")
    ax.text(-6.3, 5.6 - 0.32*ci, lab, color=cols[ci % 6], fontsize=9, fontweight="bold")
ax.set_xlim(-6.5, 6.5); ax.set_ylim(2.4, 6.0); ax.axis("off"); ax.set_aspect("equal")
ax.set_title("A cyclic rule's corners split into five upward intervals, each starting at a 3-region corner\n(each interval is one conditional tripartite information)")
ax = axs[1]
names = list(tf["result"].keys()); ok = [tf["result"][n][0] for n in names]; tot = [tf["result"][n][1] for n in names]
ax.barh(names, tot, color=OR, alpha=0.3, label="rules in family"); ax.barh(names, ok, color=BL, alpha=0.8, height=0.5, label="split into intervals")
for i, (o, t) in enumerate(zip(ok, tot)): ax.text(t+1, i, f"{o}/{t}", va="center", fontsize=8)
ax.invert_yaxis(); ax.set_xlabel("rules"); ax.legend(fontsize=8, loc="upper right")
ax.set_title("All 357 rules split this way.\nKnown up to 6 regions; open in general\n(Hubeny & Liu 2026 conjecture)")
save(fig, "fig15-rules-split-into-tripartite-intervals")

# ---------- fig16: the super log is not any shape; it is the change between shapes
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
ax = axs[0]
S = np.linspace(-1.5, 2.5, 200)
ax.plot(S, S, color=BL, lw=2.5, label="best-fitting conic: k² − 4d = s (exact)")
ax.axhline(0, color=GR, lw=1); ax.fill_between(S, -1.8, 0, color=OR, alpha=0.10); ax.fill_between(S, 0, 2.8, color=BL, alpha=0.08)
ax.text(-1.0, -1.25, "below addition:\nlooks like a circle (angle)", color=OR, fontsize=9)
ax.text(0.9, 2.1, "above addition:\nlooks like a hyperbola (radius)", color=BL, fontsize=9)
for sv, lab in [(0, "addition\n(straight)"), (1, "a+b+ab\n(exact conic)"), (2, "rung 2")]:
    ax.plot(sv, sv, "ko"); ax.annotate(lab, (sv, sv), xytext=(6, -22), textcoords="offset points", fontsize=8)
ax.set_xlabel("rung height s"); ax.set_ylabel("discriminant k² − 4d of closest plane-dial shape"); ax.set_ylim(-1.8, 2.8)
ax.legend(fontsize=8, loc="upper left"); ax.set_title("Fit a plane-dial shape to the rung at height s\n(match bend and third-order term)")
ax = axs[1]
ax.plot(S, S*(S-1)/48, color=RD, lw=2.5); ax.axhline(0, color=GY, lw=0.8)
for sv in (0, 1): ax.plot(sv, 0, "o", color=GR, ms=8)
ax.annotate("only heights 0 and 1\nare exactly a shape", (1, 0), xytext=(1.2, 0.035), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("rung height s"); ax.set_ylabel("fourth-order gap  s(s−1)/48")
ax.set_title("The rung leaves every conic at fourth order\n(so no single shape carries the super log)")
ax = axs[2]
vals = {"plane dial, angle\n(k=1, d=2)": [-1.0]*5, "plane dial, radius\n(k=3, d=1)": [8.0]*5,
        "rung 2\n(e^(e^t−1) − 1)": [3.493483, 3.497419, 3.486762, 3.492620, 3.470510]}
for i, (lab, v) in enumerate(vals.items()):
    base = v[0] if i < 2 else 3.5
    ax.scatter([i]*5, [x - base for x in v], s=50, color=[OR, BL, RD][i], zorder=3)
ax.axhline(0, color=GY, lw=0.8); ax.set_xticks(range(3)); ax.set_xticklabels(list(vals.keys()), fontsize=8)
ax.set_ylabel("genus of CP² by fixed points, minus its true value")
ax.set_title("Same sum, five random torus weights:\nshapes give one number (rigid), rung 2 wobbles")
save(fig, "fig16-super-log-is-the-change-between-shapes")

# ---------- fig17: literature board
rows = [
 ("Composing adds radii", "Formal group law logarithm; Tsallis q-sum = multiplicative law (Lazard; Tsallis 1988; Tempesta)", "Flip = the known k ↔ −k duality (Tsallis q ↔ 2−q)", "No lead: textbook"),
 ("Rungs as stacked radii", "Regular iteration, Écalle–Voronin, Szekeres; level-index numbers (Clenshaw–Olver 1984)", "Level-index uses a rough linear in-between; ours is the regular one", "Small: a canonical level-index"),
 ("Plane dial (angle / radius)", "Hirzebruch's two-parameter Todd genus; rigid genera are exactly these (Krichever; Musin 2011)", "Dial ↔ swap-rule q = ratio of roots: ours as a bridge", "Lead: plugs our dial into genus theory"),
 ("Swap rule = d × area", "Gauss–Bonnet, Berry phase, quantum torus, q-binomial area (Pólya 1969), Inönü–Wigner", "Face weights = genera of Grassmannians (known); tied to our dial (ours)", "Folded into the dial lead"),
 ("Right shape per family", "Moment polytopes; marginal polytopes; generalized permutohedra (Postnikov); elliptope", "Permutohedron = flag variety's moment polytope", "No new lead alone"),
 ("Entropy one-sign", "Already known: Hernández-Cuenca, Hubeny & Jia, Prop. 3 (2023)", "Tripartite form: known ≤ 6 regions, open in general (Hubeny & Liu 2026)", "Small: corner-interval picture of an open conjecture"),
 ("Lee–Yang, amoebas, solitons", "Lee–Yang 1952; amoebas (Gelfand–Kapranov–Zelevinsky, Mikhalkin); Maslov; ultradiscrete KdV, box-ball", "Soliton heights close to Kodama–Williams", "No new lead"),
 ("Super log and shapes", "Abel function of e^x − 1 (Szekeres, Écalle); divergence (Baker)", "No shape carries it; it is the change between shapes (ours)", "Answers the question; no new lead from shapes"),
]
fig, ax = plt.subplots(figsize=(15, 8.6)); ax.axis("off")
hdr = ["Simplification", "Already known (by whom)", "Close, or what looks ours", "Lead?"]
xs = [0.0, 0.17, 0.53, 0.82]; ws = [0.16, 0.35, 0.28, 0.18]
for x, h in zip(xs, hdr): ax.text(x, 0.98, h, fontweight="bold", fontsize=11, transform=ax.transAxes)
for r, row in enumerate(rows):
    y = 0.92 - r*0.115
    lead = row[3]; col = BL if lead.startswith("Lead") else (GR if lead.startswith("Answers") or lead.startswith("Small") or lead.startswith("Folded") else GY)
    ax.add_patch(FancyBboxPatch((xs[3]-0.005, y-0.085), 0.18, 0.095, boxstyle="round,pad=0.003", color=col, alpha=0.15, transform=ax.transAxes))
    for c, (x, w, txt) in enumerate(zip(xs, ws, row)):
        colr = ["k", OR, "#333333", col][c]
        ax.text(x, y, textwrap.fill(txt, int(w*118)), va="top", fontsize=8.6, color=colr, transform=ax.transAxes, fontweight="bold" if c in (0, 3) else None)
ax.set_title("What the literature already has for each polar simplification (orange = known; middle column says which parts look ours)", fontsize=12)
save(fig, "fig17-literature-board")
print("done")
