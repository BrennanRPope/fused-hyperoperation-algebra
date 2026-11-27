"""DRAFT picture for round 8. Run after checks.py. Orange = known (textbook formula), blue = from loops."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [2, 1]})
rows = R["einstein"]["rows"]
xx = np.arange(len(rows))
ax[0].bar(xx - 0.2, [b for _, _, b in rows], 0.4, color=OR, label="Einstein's equations: −√g G^{μν} (textbook formula)")
ax[0].bar(xx + 0.2, [a for _, a, _ in rows], 0.4, color=BL, label="change in the action read from loops")
ax[0].axhline(0, color=GR, lw=0.8)
ax[0].set_xticks(xx); ax[0].set_xticklabels([f"g{c[0]}{c[1]}" for c, _, _ in rows])
ax[0].set_xlabel("metric component nudged"); ax[0].set_ylabel("δS / δ(nudge)")
ax[0].set_title(f"A. Gravity: all 10 of Einstein's equations from varying the loop action\n(generic 4D metric; worst difference {R['einstein']['worst_rel']:.0e} of the largest)")
ax[0].legend(fontsize=8)
rs = R["maxwell"]["rows"]
x2 = np.arange(len(rs))
ax[1].bar(x2 - 0.2, [b for _, _, b in rs], 0.4, color=OR, label="Maxwell: κ² ∂_μF_{μν}")
ax[1].bar(x2 + 0.2, [a for _, a, _ in rs], 0.4, color=BL, label="5D loop action, varied")
ax[1].axhline(0, color=GR, lw=0.8)
ax[1].set_xticks(x2); ax[1].set_xticklabels([f"A{c}" for c, _, _ in rs])
ax[1].set_xlabel("potential component nudged")
ax[1].set_title("B. Light: the same 5D loop action,\nvaried in A, gives Maxwell's equations")
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/f15_einstein_and_maxwell_from_loops.png", bbox_inches="tight")
print("f15 written")
