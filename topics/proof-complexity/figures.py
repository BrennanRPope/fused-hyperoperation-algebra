"""Figures for the proof-complexity page. Run: python3 figures.py (reads yard.json, tree.tsv)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "yard.json")))
by = {(r["result"], r["route"]): r for r in rows}
results = sorted({r["result"] for r in rows})
short = {"R1 subset expansion": "R1 subset\nexpansion", "R2 inclusion-exclusion": "R2 inclusion–\nexclusion",
         "R3 Moebius heart": "R3 Möbius on\nthe cube", "R4 binomial": "R4 binomial\n(1+x)ⁿ",
         "R5 geometric series": "R5 geometric\nseries", "R6 tangent addition": "R6 tangent\naddition"}
COL = {"standard": "#5b6c8f", "operation": "#d0453b", "control": "#2ca58d"}
NAME = {"standard": "Standard (Mathlib)", "operation": "Through the operation", "control": "Control: same plan, no operation"}

def cost(res, route, vs):
    r = by[(res, route)]
    return r["own_lines"] + r[f"unshared_lines_vs_{vs}"]

# Figure 1: route-specific cost
fig, ax = plt.subplots(figsize=(11, 5.2))
x = np.arange(len(results)); w = 0.27
for j, route in enumerate(["standard", "operation", "control"]):
    vals = []
    for res in results:
        if (res, route) not in by: vals.append(np.nan); continue
        vs = "operation" if route == "standard" else "standard"
        vals.append(cost(res, route, vs))
    b = ax.bar(x + (j - 1) * w, vals, w, color=COL[route], label=NAME[route])
    for xi, v in zip(x + (j - 1) * w, vals):
        if not np.isnan(v): ax.text(xi, v * 1.06, f"{int(v)}", ha="center", va="bottom", fontsize=8)
ax.set_yscale("log"); ax.set_ylim(30, 20000)
ax.set_xticks(x); ax.set_xticklabels([short[r] for r in results])
ax.set_ylabel("lines of Lean this route needs\nthat the standard route does not (log scale)")
ax.set_title("Route-specific cost: the proof itself plus every lemma under it that the other route doesn't need\n"
             "(standard is measured against the operation route; operation and control against standard)", fontsize=10)
ax.axvspan(-0.5, 2.5, color="#2ca58d", alpha=0.06); ax.axvspan(2.5, 5.5, color="#d0453b", alpha=0.06)
ax.text(1, 9000, "operation route cheaper than Mathlib…\n…but the control is cheaper still", ha="center", fontsize=9)
ax.text(4, 9000, "operation route costlier than Mathlib", ha="center", fontsize=9)
ax.legend(loc="upper right", bbox_to_anchor=(1, 0.86), fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig1-route-cost.png"), dpi=150); plt.close(fig)

# Figure 2: surface length and proof-term size
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, key, lab in [(axs[0], "own_lines", "lines in the proof itself"),
                     (axs[1], "term", "proof-term size (distinct sub-terms, log scale)")]:
    for j, route in enumerate(["standard", "operation", "control"]):
        vals = [by[(res, route)][key] if (res, route) in by else np.nan for res in results]
        ax.bar(x + (j - 1) * w, vals, w, color=COL[route], label=NAME[route])
        for xi, v in zip(x + (j - 1) * w, vals):
            if not np.isnan(v): ax.text(xi, v * (1.06 if key == "term" else 1) + (0 if key == "term" else 0.4),
                                        f"{int(v)}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([short[r] for r in results], fontsize=8)
    ax.set_ylabel(lab)
    if key == "term": ax.set_yscale("log")
axs[0].set_title("Surface length: Mathlib's proof is never longer", fontsize=10)
axs[1].set_title("Kernel size of the proof: Mathlib's is always smallest", fontsize=10)
axs[0].legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig2-surface-length.png"), dpi=150); plt.close(fig)

# Figure 3: head to head, operation vs control
fig, ax = plt.subplots(figsize=(9.5, 4.4))
res5 = [r for r in results if (r, "control") in by]
x5 = np.arange(len(res5)); w2 = 0.36
o = [cost(r, "operation", "control") for r in res5]; c = [cost(r, "control", "operation") for r in res5]
fw = [by[(r, "operation")]["framework_lines"] for r in res5]
ax.bar(x5 - w2 / 2, o, w2, color=COL["operation"], label="operation route: its own lines + what the control doesn't need")
ax.bar(x5 - w2 / 2, fw, w2, color="#8a1f17", label="…of which the framework's setup lemmas (chain, lift, corners)")
ax.bar(x5 + w2 / 2, c, w2, color=COL["control"], label="control: its own lines + what the operation route doesn't need")
for xi, v in zip(x5 - w2 / 2, o): ax.text(xi, v + 4, str(v), ha="center", fontsize=8)
for xi, v in zip(x5 + w2 / 2, c): ax.text(xi, v + 4, str(v), ha="center", fontsize=8)
ax.set_xticks(x5); ax.set_xticklabels([short[r] for r in res5], fontsize=8)
ax.set_ylabel("lines of Lean")
ax.set_title("Same proof plan, with and without the operation: the operation only adds cost", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.set_ylim(0, max(o) * 1.45)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig3-operation-vs-control.png"), dpi=150); plt.close(fig)

# Figure 4: dependency trees for R3 (two levels of named lemmas)
tree = {}
for line in open(os.path.join(HERE, "tree.tsv")):
    root, d, n, lines, kids = line.rstrip("\n").split("\t")
    kids = [k.strip() for k in kids.strip("[]").split(",") if k.strip()]
    tree.setdefault(root, {})[n] = (int(d), int(lines), kids)
roots = [("Finset.sum_powerset_neg_one_pow_card_of_nonempty", "Standard (Mathlib)"),
         ("HCT.Results.op_sum_powerset_neg_one_pow_card_of_nonempty", "Through the operation"),
         ("HCT.Control.ctl_sum_powerset_neg_one_pow_card_of_nonempty", "Control: same plan, no operation")]
def colour(n):
    if n.startswith("HCT.Results") or n.startswith("HCT.Control.ctl_sum") or n == "Finset.sum_powerset_neg_one_pow_card_of_nonempty": return "#333333"
    if n.startswith("HCT.Control"): return "#2ca58d"
    if n.startswith("HCT."): return "#d0453b"
    if "choose" in n or "apply_card" in n or "powerset_neg_one_pow_card" in n or "nsmul" in n: return "#5b6c8f"
    if "powerset_insert" in n or "induction" in n or "prod_insert" in n: return "#2ca58d"
    return "#999999"
fig, axs = plt.subplots(1, 3, figsize=(19, 7.5))
for ax, (root, title) in zip(axs, roots):
    t = tree[root]
    lvl1 = t[root][2]
    ax.set_xlim(-0.2, 3.6); ax.set_ylim(-0.5, 0.5 + max(14, sum(max(1, len([k for k in t.get(n, (0, 0, []))[2] ][:6])) for n in lvl1)))
    ax.axis("off"); ax.set_title(title, fontsize=11)
    ytop = ax.get_ylim()[1] - 1
    ax.text(0, ytop / 2, root.split(".")[-1].replace("_of_nonempty", "\n_of_nonempty"), fontsize=7.5, ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white", ec=colour(root)), color=colour(root))
    y = ytop
    for n in lvl1:
        kids = t.get(n, (0, 0, []))[2][:6]
        span = max(1, len(kids))
        yc = y - (span - 1) / 2
        ax.plot([1.05, 1.25], [ytop / 2, yc], color="#cccccc", lw=0.8)
        ax.text(1.27, yc, f"{n.split('.')[-1]}  ({t[n][1]})", fontsize=7, va="center", color=colour(n), fontweight="bold" if n.startswith("HCT") else None)
        for i, k in enumerate(kids):
            yk = y - i
            ax.plot([2.3, 2.45], [yc, yk], color="#dddddd", lw=0.6)
            kl = t[k][1] if k in t else ""
            ax.text(2.47, yk, f"{k.split('.')[-1]}" + (f"  ({kl})" if kl != "" else ""), fontsize=6.3, va="center", color=colour(k))
        if len(t.get(n, (0, 0, []))[2]) > 6:
            ax.text(2.47, y - 6, f"… +{len(t[n][2]) - 6} more", fontsize=6, color="#999999", va="center")
            span += 1
        y -= span + 0.4
axs[0].text(0, -0.4, "Blue: counting subsets by size and alternating binomial sums.\n(Numbers in brackets: lines of that lemma.)", fontsize=8, color="#5b6c8f", transform=axs[0].transAxes)
axs[1].text(0, -0.4, "Red: the framework's own lemmas (chain, lift, corners).\nGreen: peeling one element off a set.", fontsize=8, color="#d0453b", transform=axs[1].transAxes)
axs[2].text(0, -0.4, "Green: the same peeling, done directly on the product,\nwith no chain or lift in between.", fontsize=8, color="#2ca58d", transform=axs[2].transAxes)
fig.suptitle("R3, the signed count of a cube's corners is 0: two levels of the lemmas each proof cites", fontsize=12)
fig.tight_layout(rect=(0, 0.06, 1, 0.96)); fig.savefig(os.path.join(HERE, "fig4-dependency-trees-R3.png"), dpi=150); plt.close(fig)

# Figure 5: scoreboard
yard = [("Own lines", lambda r, rt: by[(r, rt)]["own_lines"]),
        ("Proof-term size", lambda r, rt: by[(r, rt)]["term"]),
        ("Lemmas cited\n(kernel)", lambda r, rt: by[(r, rt)]["direct"]),
        ("Dependency\ndepth", lambda r, rt: by[(r, rt)]["depth"]),
        ("Route-specific\ncost", lambda r, rt: cost(r, rt, "operation" if rt == "standard" else "standard"))]
fig, ax = plt.subplots(figsize=(10.5, 4.6))
ax.set_xlim(0, len(yard)); ax.set_ylim(0, len(results)); ax.invert_yaxis(); ax.axis("off")
for i, res in enumerate(results):
    ax.text(-0.05, i + 0.5, short[res].replace("\n", " "), ha="right", va="center", fontsize=9)
    for j, (lab, f) in enumerate(yard):
        routes = [rt for rt in ("standard", "operation", "control") if (res, rt) in by]
        vals = {rt: f(res, rt) for rt in routes}
        best = min(vals.values()); winners = [rt for rt in routes if vals[rt] == best]
        colr = COL[winners[0]] if len(winners) == 1 else "#cccccc"
        ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, color=colr, alpha=0.85))
        txt = " / ".join(str(vals[rt]) for rt in routes)
        ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=7.5, color="white")
for j, (lab, _) in enumerate(yard):
    ax.text(j + 0.5, -0.15, lab, ha="center", va="bottom", fontsize=9)
ax.set_title("Who wins each yardstick (lower is better). Cell: standard / operation / control. Colour = winner.\n"
             "Red (the operation) never wins on its own; grey = tie.", fontsize=10, pad=34)
for k, rt in enumerate(["standard", "operation", "control"]):
    ax.add_patch(plt.Rectangle((k * 1.7, len(results) + 0.25), 0.25, 0.4, color=COL[rt], clip_on=False))
    ax.text(k * 1.7 + 0.32, len(results) + 0.45, NAME[rt], va="center", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig5-scoreboard.png"), dpi=150); plt.close(fig)
print("figures written")
