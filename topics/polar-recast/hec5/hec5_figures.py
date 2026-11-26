"""Pictures for the 5-region sign check (DRAFT). Orange = known, blue = ours, red = carved (weight below zero)."""
import json, itertools, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OR, BL, RD, GR, GY = "#e07b00", "#1f6fd1", "#c0392b", "#2a9d55", "#888888"
rep = json.load(open("hec5_signs.json"))
CORN = {g: ["".join(J) for J in itertools.combinations("ABCDE", g)] for g in range(1, 6)}
def corner_panel(ax, coeffs, title):
    for g in range(1, 6):
        items = CORN[g]
        for i, J in enumerate(items):
            x = i - (len(items) - 1)/2; c = coeffs.get(J, 0)
            col = "#eeeeee" if c == 0 else (RD if c == -1 else ("#7b1010" if c < -1 else GR))
            ax.scatter(x, g, s=95, color=col, edgecolor="k", lw=0.4, zorder=3)
            if c: ax.text(x, g - 0.42, J, ha="center", fontsize=5)
    ax.set_ylim(0.4, 5.5); ax.set_xlim(-5.3, 5.3); ax.set_yticks(range(1, 6)); ax.set_yticklabels(["1", "2", "3", "4", "5"], fontsize=7)
    ax.set_xticks([]); ax.set_title(title, fontsize=9)
fig, axs = plt.subplots(2, 4, figsize=(15, 6.6))
for ax, k in zip(axs.flat, map(str, range(1, 9))):
    r = rep[k]; ex = r["example"]
    key = [kk for kk in ex if "grade 3" in kk] or list(ex)
    corner_panel(ax, ex[key[0]], f"{k}. {r['name']}: {r['instances']} rules\n" + ("all one sign" if all("one sign" in s for s in r["stats"] if "sign" in s) else "MIXED"))
fig.supylabel("number of regions in the corner", fontsize=9)
fig.suptitle("All 372 five-region holographic rules on the cube's corners (one example per family). "
             "Red = weight −1, dark red = −2, grey = untouched.\nOutside subadditivity (panel 1), no rule touches a single region or a pair, and none has a plus sign.", fontsize=11)
fig.tight_layout(); fig.savefig("../figures/fig12-every-5-region-rule-on-corners.png", dpi=130); plt.close(fig)

# fig13: net carving: GHZ vs perfect tensor on one cyclic rule
group = ["ABE", "BCE", "CDE", "ABCE", "BCDE"]
y_pt = {3: 0, 4: -2, 5: 6}; y_ghz = {3: 1, 4: -1, 5: 1}
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, (nm, y, col) in zip(axs, [("6-qubit perfect tensor (holographic)", y_pt, BL), ("6-qubit GHZ (not holographic)", y_ghz, RD)]):
    vals = [y[len(J)] for J in group]; tot = sum(vals)
    ax.bar(range(len(group)), vals, color=[GR if v > 0 else (RD if v < 0 else GY) for v in vals])
    ax.bar(len(group) + 0.6, tot, color=BL if tot <= 0 else RD, hatch="//")
    ax.axhline(0, color="k", lw=0.8); ax.set_xticks(list(range(len(group))) + [len(group) + 0.6]); ax.set_xticklabels(group + ["SUM"])
    ax.set_title(f"{nm}\nsum over the rule's corners = {tot:+d} → " + ("passes (net carved)" if tot <= 0 else "fails (net added)"), fontsize=10)
axs[0].set_ylabel("corner weight")
fig.suptitle("One cyclic rule: the corners it names must add to zero or less.\nSingle corners may be positive (the perfect tensor's top corner is +6), but each named group must net out carved.", fontsize=10)
fig.tight_layout(); fig.savefig("../figures/fig13-net-carved-not-every-corner.png", dpi=130); plt.close(fig)
print("ok")
