"""Draws the dependency-map pictures. Run: python3 -I figures.py (needs matplotlib, networkx)."""
import os, sys, textwrap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch, Rectangle
from matplotlib.lines import Line2D
import networkx as nx
from data import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def box(ax, x, y, w, h, text, fc, ec="#333", fs=10, tc="black", lw=1.2, alpha=1):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=lw, alpha=alpha))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc)

# ------------------------------------------------------------------ d1 scorecard
def d1():
    groups = []
    for p in PREREQS:
        if p[1] not in groups: groups.append(p[1])
    rows = []
    for g in groups:
        rows.append(("H", g))
        rows += [("P", p) for p in PREREQS if p[1] == g]
    n = len(rows)
    fig, ax = plt.subplots(figsize=(17, 0.7 * n + 1.8))
    ax.set_xlim(0, 17); ax.set_ylim(-0.6, n + 0.9); ax.axis("off")
    ax.text(0.1, n + 0.45, "What any quantum-gravity answer (or G dial) has to meet, and where our results stand",
            fontsize=15, weight="bold")
    for i, (kind, r) in enumerate(rows):
        y = n - 1 - i
        if kind == "H":
            ax.text(0.1, y, r, fontsize=12.5, weight="bold", color="#1b2631", va="center")
            ax.plot([0.1, 16.9], [y - 0.32, y - 0.32], color="#bbb", lw=0.8)
            continue
        pid, g, name, corner, st, have, rests, lit = r
        lab, col = STATUS_STYLE[st]
        ax.add_patch(Rectangle((0.25, y - 0.22), 0.45, 0.44, fc=col, ec="none"))
        ax.text(0.85, y + 0.17, f"{pid}  {name}", fontsize=10.5, va="center", weight="bold")
        ax.text(0.85, y - 0.2, textwrap.shorten(have, 150, placeholder=" ..."), fontsize=9, va="center", color="#444")
        ax.text(16.9, y + 0.17, lab, fontsize=9.5, va="center", ha="right", color=col, weight="bold")
    handles = [Patch(color=c, label=l) for k, (l, c) in STATUS_STYLE.items()]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=6, frameon=False, fontsize=10)
    ax.text(0.1, n + 0.05, "bold: the prerequisite.  grey: what our results do with it.  right: verdict", fontsize=10, color="#666", style="italic")
    save(fig, "d1-prerequisite-scorecard.png")

# ------------------------------------------------------------------ d2 cube
def d2():
    # Bronstein cube in oblique projection: axes c (right), G (up), hbar (depth)
    P = lambda c, g, h: (c * 6 + h * 2.6, g * 5 + h * 2.2)
    corners = {"": (0,0,0), "c": (1,0,0), "G": (0,1,0), "h": (0,0,1), "cG": (1,1,0),
               "ch": (1,0,1), "Gh": (0,1,1), "cGh": (1,1,1)}
    names = {"": "Newtonian\nmechanics", "c": "special\nrelativity", "G": "Newtonian\ngravity",
             "h": "quantum\nmechanics", "cG": "general\nrelativity", "ch": "quantum\nfield theory",
             "Gh": "Newtonian\nquantum gravity", "cGh": "QUANTUM\nGRAVITY"}
    fig, ax = plt.subplots(figsize=(17, 10.5)); ax.axis("off"); ax.set_aspect("equal")
    for a, ca in corners.items():
        for b, cb in corners.items():
            if sum(abs(x - y) for x, y in zip(ca, cb)) == 1 and a < b:
                (x1, y1), (x2, y2) = P(*ca), P(*cb)
                # which constant this edge switches on
                d = [i for i in range(3) if ca[i] != cb[i]][0]
                col = ["#1e8449", "#c0392b", "#2471a3"][d]
                ax.plot([x1, x2], [y1, y2], color=col, lw=3.5 if d != 1 else 3.5,
                        ls="-" if d != 1 else (0, (4, 3)), alpha=0.75, zorder=1)
    by_corner = {}
    for p in PREREQS:
        by_corner.setdefault(p[3], []).append(p)
    # each corner's list sits in a side column, joined to its corner by a thin leader line
    lists = {"Gh": (-1.6, 7.9, "right"), "G": (-1.6, 5.4, "right"), "h": (-1.6, 2.6, "right"),
             "cGh": (10.4, 9.7, "left"), "cG": (10.4, 5.9, "left"), "ch": (10.4, 3.1, "left"),
             "c": (10.4, 0.6, "left")}
    for k, (lx, ly, ha) in lists.items():
        items = by_corner.get(k, [])
        if not items: continue
        x, y = P(*corners[k])
        ymid = ly - 0.21 * (len(items) - 1)
        ax.plot([x, lx + (0.1 if ha == "right" else -0.1)], [y, ymid], color="#95a5a6", lw=1, zorder=0.5)
        for j, p in enumerate(items):
            lab, col = STATUS_STYLE[p[4]]
            ax.text(lx, ly - j * 0.42, f"{p[0]} {textwrap.shorten(p[2], 40, placeholder='...')}",
                    fontsize=8.8, ha=ha, va="center", color="white", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.22", fc=col, ec="none"))
    for k, c3 in corners.items():
        x, y = P(*c3)
        box(ax, x, y, 2.0 if k != "cGh" else 2.2, 0.75, names[k], "#ffffff" if k != "cGh" else "#fadbd8",
            ec="#1b2631", fs=8.6, lw=1.6 if k != "cGh" else 2.4)
    ax.text(*P(0.5, -0.13, 0), "switch on c  (stable deformation: proved)", color="#1e8449",
            ha="center", va="center", fontsize=10.5, weight="bold")
    ax.text(-0.45, 2.55, "switch on G (no stability theorem)", color="#c0392b", ha="center", va="center",
            rotation=90, fontsize=9.5, weight="bold")
    xh, yh = P(0, 0, 0.5)
    ax.text(xh + 0.45, yh - 0.35, "switch on hbar\n(stable: proved)", color="#2471a3", ha="left", va="center",
            fontsize=10.5, weight="bold")
    handles = [Patch(color=c, label=l) for k, (l, c) in STATUS_STYLE.items()]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=10)
    ax.set_title("The prerequisites placed on the Bronstein cube (which constants each one involves)\n"
                 "green and blue edges are proven deformations; the red dashed G edges are not",
                 fontsize=13.5, weight="bold")
    ax.set_xlim(-6.2, 15.6); ax.set_ylim(-2.6, 10.2)
    save(fig, "d2-prerequisites-on-the-cube.png")

# ------------------------------------------------------------------ d3 dependency graph
def d3():
    fig, ax = plt.subplots(figsize=(25, 16)); ax.axis("off")
    X = lambda layer: layer * 4.3
    YS = 1.3
    pos = {k: (X(v[1]), v[2] * YS) for k, v in NODES.items()}
    needs = {}
    for k, (lab, targets) in INPUTS.items():
        for t in targets:
            needs.setdefault(t, []).append(lab.replace("\n", " "))
    W, H = 3.35, 1.05
    for a, b, t, note in EDGES:
        lab, col, ls, lw = EDGE_STYLE[t]
        (x1, y1), (x2, y2) = pos[a], pos[b]
        if NODES[a][1] == NODES[b][1]:
            rad = 0.45 if y2 < y1 else -0.45
            p1, p2 = (x1 - W/2, y1), (x2 - W/2, y2)
            cs = f"arc3,rad={rad}"
        else:
            p1, p2 = (x1 + W/2, y1), (x2 - W/2, y2)
            cs = "arc3,rad=0.05"
        ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=15, color=col, lw=lw, ls=ls,
                     connectionstyle=cs, zorder=1, alpha=0.85))
    layer_fc = ["#eaf2f8", "#e8f8f5", "#fef9e7", "#f4ecf7", "#fdf2e9", "#fadbd8"]
    for k, (lab, layer, y, folder) in NODES.items():
        x, yy = pos[k]
        ec = "#c0392b" if k == "GD" else "#34495e"
        box(ax, x, yy + (0.12 if k in needs else 0), W, H, lab, layer_fc[layer], ec=ec, fs=9.3,
            lw=2.4 if k in ("LED", "GD") else 1.2)
        if k in needs:
            txt = "needs: " + "; ".join(needs[k])
            ax.text(x, yy - H/2 + 0.02, textwrap.fill(txt, 46), ha="center", va="top", fontsize=7.8,
                    color="white", weight="bold", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.18", fc="#c0392b", ec="none"))
            ax.text(x, yy - H/2 - 0.42, folder, ha="center", va="top", fontsize=7, color="#7b7d7d", style="italic")
        else:
            ax.text(x, yy - H/2 - 0.08, folder, ha="center", va="top", fontsize=7, color="#7b7d7d", style="italic")
    for i, t in enumerate(["foundation", "structure", "physics results", "deeper results", "bridges and ledger", "open"]):
        ax.text(X(i), 9.6 * YS + 0.6, t, ha="center", fontsize=13.5, weight="bold", color="#1b2631")
    handles = [Line2D([0], [0], color=c, lw=lw + 0.5, ls=ls, label=l) for k, (l, c, ls, lw) in EDGE_STYLE.items()]
    handles.append(Patch(color="#c0392b", label="red tag: outside input, put in by hand"))
    ax.legend(handles=handles, loc="lower right", fontsize=11.5, frameon=True, title="how each link is backed",
              title_fontsize=12.5)
    ax.set_xlim(-4.0, X(5) + 2.2); ax.set_ylim(-1.6 * YS, 10.0 * YS + 0.6)
    ax.set_title("How the results on the path to quantum gravity depend on each other (arrow = 'is used by')",
                 fontsize=17, weight="bold")
    save(fig, "d3-dependency-graph.png")

# ------------------------------------------------------------------ d4 backing counts
def d4():
    from collections import Counter
    cnt = Counter(t for _, _, t, _ in EDGES)
    nin = sum(len(v[1]) for v in INPUTS.values())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw={"width_ratios": [1.1, 1]})
    keys = ["L", "N", "C", "P"]
    vals = [cnt[k] for k in keys] + [nin]
    labs = [EDGE_STYLE[k][0] for k in keys] + ["outside input put in by hand"]
    cols = [EDGE_STYLE[k][1] for k in keys] + ["#c0392b"]
    a1.barh(labs[::-1], vals[::-1], color=cols[::-1])
    for i, v in enumerate(vals[::-1]):
        a1.text(v + 0.3, i, str(v), va="center", fontsize=12, weight="bold")
    a1.set_title(f"All {len(EDGES)} links on the path, plus {nin} input links\n({len(INPUTS)} outside inputs, some used twice)", fontsize=13, weight="bold")
    a1.spines[["top", "right"]].set_visible(False); a1.set_xlabel("number of links")
    # where the Lean links sit: by layer of the target
    layer_names = ["foundation", "structure", "physics", "deeper", "bridges/ledger", "open"]
    import numpy as np
    M = np.zeros((6, 4))
    for a, b, t, _ in EDGES:
        M[NODES[b][1], keys.index(t)] += 1
    left = np.zeros(6)
    for j, k in enumerate(keys):
        a2.barh(layer_names[::-1], M[::-1, j], left=left[::-1], color=EDGE_STYLE[k][1], label=EDGE_STYLE[k][0])
        left += M[:, j]
    a2.set_title("Lean thins out as the claims get physical\n(links counted by the layer they lead into)", fontsize=13, weight="bold")
    a2.legend(fontsize=9, frameon=False, loc="lower right")
    a2.spines[["top", "right"]].set_visible(False); a2.set_xlabel("number of links")
    save(fig, "d4-how-links-are-backed.png")

# ------------------------------------------------------------------ d5 Lean islands
def d5():
    # rows: each module that imports Base starts a row; its own importers extend the row to the right
    rows = [["Thermodynamics", "Relativity", "TwistedBoosts", "Electromagnetism"], ["GroupsLie", "Clifford"],
            ["Cube", "E"], ["Pi"], ["OrderDial"], ["GeneralRelativity"], ["GRBridges"], ["GravityProbeB"],
            ["BinaryPulsar"]]
    W, H, DX, DY = 2.7, 0.72, 3.4, 0.98
    pos = {}
    n = len(rows)
    for r, row in enumerate(rows):
        for c, m in enumerate(row):
            pos[m] = (3.6 + c * DX, (n - 1 - r) * DY)
    pos["Base"] = (0, (n - 1) * DY / 2)
    isl = list(ISLAND_NOTES)
    for j, m in enumerate(isl):
        pos[m] = (18.2, (n - 1) * DY - 0.4 - j * 1.75)
    fig, ax = plt.subplots(figsize=(20, 9.4)); ax.axis("off")
    def arrow(p, q, rad=0.0):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13, color="#566573", lw=1.3,
                                     connectionstyle=f"arc3,rad={rad}", zorder=1))
    bx, by = pos["Base"]
    for row in rows:
        x, y = pos[row[0]]
        arrow((bx + W/2, by), (x - W/2, y))
        for a_, b_ in zip(row, row[1:]):
            (x1, y1), (x2, y2) = pos[a_], pos[b_]
            arrow((x1 + W/2, y1), (x2 - W/2, y2))
    for m, (x, y) in pos.items():
        nth = LEAN[m][0]
        island = m in ISLAND_NOTES
        fc = "#fadbd8" if island else ("#d6eaf8" if m == "Base" else "#e8f8f5")
        box(ax, x, y, W, H, f"{m}\n{nth} theorems", fc, ec="#c0392b" if island else "#34495e",
            fs=9.5, lw=2 if island else 1.2)
        if island:
            ax.text(x + W/2 + 0.15, y, ISLAND_NOTES[m], fontsize=9.5, va="center", color="#922b21")
    top = (n - 1) * DY + 0.75
    ax.add_patch(Rectangle((16.5, -0.75), 5.9, top + 0.55, fc="none", ec="#c0392b", ls="--", lw=1.6))
    ax.text(19.45, top + 0.35, "islands: import nothing,\nre-define what Base already proves", ha="center",
            fontsize=11.5, color="#c0392b", weight="bold")
    ax.text(-1.4, -1.35, "Arrows: 'imported by' (an import that a longer chain already carries is not drawn). "
            "All 20 modules are in the library's build list, and the README records each first compile.\n"
            "Every path theorem is an identity over Int, Rat or Gaussian integers; none states a limit, a convergence, "
            "a measure or a complex exponential.", fontsize=10.5, color="#333")
    ax.set_xlim(-1.6, 22.6); ax.set_ylim(-1.8, top + 1.2)
    ax.set_title("The Lean modules behind the quantum-gravity path: the newest five stand alone", fontsize=15, weight="bold")
    save(fig, "d5-lean-islands.png")

# ------------------------------------------------------------------ d6 G yardstick
def d6():
    nr, nc = len(YARD_ROWS), len(YARD_COLS)
    fig, ax = plt.subplots(figsize=(18, 1.25 * nr + 1.8)); ax.axis("off")
    ax.set_xlim(0, 18); ax.set_ylim(0, nr + 1.3)
    xs = [4.2, 7.65, 11.1, 15.0]; ws = [3.3, 3.3, 3.3, 4.2]
    for j, c in enumerate(YARD_COLS):
        ax.text(xs[j], nr + 0.55, {"hbar": "ħ  (order dial)", "c": "c  (signature dial)",
                                    "Lambda": "Λ  (curvature dial)", "G": "G  (?)"}[c],
                ha="center", fontsize=14, weight="bold", color="#c0392b" if c == "G" else "#1b2631")
    for i, (rname, cells) in enumerate(YARD_ROWS):
        y = nr - 1 - i + 0.5
        ax.text(0.15, y, textwrap.fill(rname, 22), va="center", fontsize=11.5, weight="bold")
        for j, cell in enumerate(cells):
            st = YARD_STATUS[i][j]
            col = STATUS_STYLE[st][1]
            ax.add_patch(FancyBboxPatch((xs[j] - ws[j]/2, y - 0.45), ws[j], 0.9, boxstyle="round,pad=0.02",
                         fc=col, ec="none", alpha=0.18 if st == "met" else 0.3))
            ax.add_patch(Rectangle((xs[j] - ws[j]/2, y - 0.45), 0.12, 0.9, fc=col, ec="none"))
            ax.text(xs[j], y, textwrap.fill(cell, 34 if j < 3 else 40), ha="center", va="center", fontsize=9.2)
    ax.set_title("The yardstick: what c, ħ and Λ have as dials, and what G is missing", fontsize=15, weight="bold")
    save(fig, "d6-g-dial-yardstick.png")

# ------------------------------------------------------------------ d7 weak links
def d7():
    fig, axs = plt.subplots(1, 3, figsize=(19, 6.4))
    # (a) one-number matches are ruler changes
    import numpy as np
    ax = axs[0]
    x = np.linspace(-0.6, 3, 300)
    for k, c, lab in [(1, "#2471a3", "the circle operation k=1"), (-0.5, "#1e8449", "kappa-Poincare (k=-1/kappa)"),
                      (-0.25, "#d35400", "Newton binding (k=-G/rc^2)")]:
        xx = x[1 + k * x > 0.02]
        ax.plot(xx, np.log(1 + k * xx) / k, color=c, lw=2.4, label=lab)
    ax.plot(x, x, color="k", lw=1.2, ls="--", label="plain addition")
    ax.set_title("(a) Every one-number law is plain addition\nafter a ruler change log(1+kx)/k", fontsize=12, weight="bold")
    ax.set_xlabel("quantity x"); ax.set_ylabel("ruler reading"); ax.legend(fontsize=9, frameon=False)
    ax.text(0.02, -0.25, "So a 1D match (velocity addition, collinear kappa-Poincare,\nNewton binding, Magueijo-Smolin) "
            "carries no evidence by itself.\nOnly the matrix (non-commuting) matches count.", transform=ax.transAxes,
            fontsize=9.5, color="#922b21", va="top")
    # (b) three different 'imaginary dials'
    ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(-2, 10)
    ax.text(5, 9.5, "(b) 'Turning the dial to i' names three\ndifferent moves", ha="center", fontsize=12, weight="bold")
    items = [("k = -i/ħ on the action", "phase per step; forced by\nlossless steps (chain-resolves)", "#2471a3"),
             ("k -> i k on a Cayley-Klein label", "boost <-> rotation;\nWick rotation of signature (c)", "#1e8449"),
             ("turn rate a -> i m", "checkerboard: a continuation\nin mass, not in time", "#d35400")]
    for i, (t, s, c) in enumerate(items):
        y = 7.0 - i * 2.6
        box(ax, 3.0, y, 5.2, 1.6, t, c, ec=c, fs=10.5, tc="white")
        ax.text(5.9, y, s, va="center", fontsize=9.6)
    ax.text(5, -0.9, "They are related (Nelson 1964; GJKS 1984), not identical.\nThe ledger's 'c: Wick flips the dial' and 'ħ: direction\n"
            "forced' lean on reading them as one dial.", ha="center", fontsize=9.5, color="#922b21")
    # (c) scope of 'G has no dial'
    ax = axs[2]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(5, 9.5, "(c) 'G has no dial' was tested\non a narrow class", ha="center", fontsize=12, weight="bold")
    ax.add_patch(plt.Circle((5, 4.6), 3.9, fc="#f4ecf7", ec="#8e44ad", lw=1.5))
    ax.text(5, 7.7, "all deformations a G dial could be", ha="center", fontsize=9.8, color="#6c3483")
    ax.add_patch(plt.Circle((3.8, 4.4), 1.6, fc="#d6eaf8", ec="#2471a3", lw=1.5))
    ax.text(3.8, 4.4, "finite Lie algebras,\nde Sitter family\n(dial-map's\nH² test)", ha="center", va="center", fontsize=8.8)
    ax.add_patch(plt.Circle((6.9, 4.0), 1.35, fc="#fdebd0", ec="#d35400", lw=1.5, ls="--"))
    ax.text(6.9, 4.0, "hypersurface-\ndeformation\nalgebroid\n(untested)", ha="center", va="center", fontsize=8.8)
    ax.text(5, -0.1, "Deformed GR (Bojowald-Paily 2012-13) bends this algebroid\nby a function β that depends on density in Planck units\n"
            "(G and ħ), and can flip the signature: a lead to test.", ha="center", fontsize=9.3, color="#6e2c00")
    fig.suptitle("Three links that are weaker than they read", fontsize=15, weight="bold", y=1.03)
    save(fig, "d7-weak-links.png")

# ------------------------------------------------------------------ d8 assumptions matrix
def d8():
    from assumptions import WRITEUPS, CATS
    import numpy as np
    nr, nc = len(WRITEUPS), len(CATS)
    fig, ax = plt.subplots(figsize=(18, 0.55 * nr + 3.2)); ax.axis("off")
    ax.set_xlim(-6.5, nc + 2.2); ax.set_ylim(-1.2, nr + 3.8)
    for j, (cid, cname) in enumerate(CATS):
        ax.text(j + 0.3, nr + 0.15, textwrap.fill(cname, 24), rotation=35, ha="left", va="bottom", fontsize=9, weight="bold", rotation_mode="anchor")
    for i, w in enumerate(WRITEUPS):
        y = nr - 1 - i
        ax.text(-6.3, y + 0.5, w["name"], va="center", fontsize=10, weight="bold")
        for j, (cid, cname) in enumerate(CATS):
            v = w["cats"].get(cid)
            if v is None:
                fc = "#fbfcfc"
            elif v == "known":
                fc = "#aed6f1"
            elif v == "ours":
                fc = "#f5b041"
            elif v == "input":
                fc = "#e6b0aa"
            elif v == "yes":
                fc = "#abebc6"
            elif v == "no":
                fc = "#f1948a"
            else:
                fc = "#d5d8dc"
            ax.add_patch(Rectangle((j + 0.05, y + 0.05), 0.9, 0.9, fc=fc, ec="#ccc"))
            mark = {"known": "K", "ours": "O", "input": "in", "yes": "✓", "no": "✗"}.get(v, "")
            ax.text(j + 0.5, y + 0.5, mark, ha="center", va="center", fontsize=9)
    handles = [Patch(color="#aed6f1", label="K: rests on a known (published) result"),
               Patch(color="#f5b041", label="O: rests on a step of ours"),
               Patch(color="#e6b0aa", label="in: an input put in by hand"),
               Patch(color="#abebc6", label="✓ section present"), Patch(color="#f1948a", label="✗ section missing")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.07), ncol=5, frameon=False, fontsize=9.5)
    ax.set_title("What each write-up on the path assumes (full text in assumptions.md)", fontsize=15, weight="bold")
    save(fig, "d8-assumptions-matrix.png")

if __name__ == "__main__":
    for f in (d1, d2, d3, d4, d5, d6, d7, d8):
        f(); print("drew", f.__name__)
