"""DRAFT picture for round 10. Run after checks.py. Blue = the law's bookkeeping, orange = ordinary spinor bookkeeping."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
L = np.array(R["loops"])
aa = np.linspace(0, 4 * math.pi, 100)
ax[0].plot(aa, aa / 2, "-", color=BL, lw=1, alpha=0.6)
ax[0].plot(L[:, 1], L[:, 5] / 2, "o", color=BL, ms=8, label="law: divide by the complex corner → phase = area / 2")
ax[0].plot(aa, 0 * aa, "-", color=OR, lw=2, label="ordinary spinor: divide by its size → no common phase")
ax[0].axvline(4 * math.pi, color=GR, ls=":", lw=1); ax[0].text(4 * math.pi - 0.2, 0.5, "whole sphere:\n½ × 4π = 1 turn", ha="right", fontsize=8, color=GR)
ax[0].set_xlabel("area enclosed by the loop on the velocity sphere")
ax[0].set_ylabel("common phase picked up (beyond spin)")
ax[0].set_title("A. The law's corner carries half the twist,\nso its objects have charge ½ under the twist")
ax[0].legend(fontsize=8, loc="upper left")
lab = ["law's bookkeeping\n(charge ½, flux 1)", "ordinary spinor\n(charge 0)"]
Lc = [R["law_counts"][0], R["spin_counts"][0]]; Rc = [R["law_counts"][1], R["spin_counts"][1]]
xx = np.arange(2)
ax[1].bar(xx - 0.2, Lc, 0.4, color=BL, label="left-handed massless")
ax[1].bar(xx + 0.2, Rc, 0.4, color=OR, label="right-handed massless")
for i in range(2):
    ax[1].text(xx[i] - 0.2, Lc[i] + 0.03, str(Lc[i]), ha="center", color=BL)
    ax[1].text(xx[i] + 0.2, Rc[i] + 0.03, str(Rc[i]), ha="center", color=OR)
ax[1].set_ylim(0, 1.25)
ax[1].set_xticks(xx); ax[1].set_xticklabels(lab); ax[1].set_yticks([0, 1])
ax[1].set_ylabel("massless modes on the hidden velocity sphere")
ax[1].set_title("B. Charge ½ gives exactly one\nmassless one-handed mode")
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/f17_charge_from_the_corner.png", bbox_inches="tight")
print("f17 written")
