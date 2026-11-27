"""DRAFT picture for round 7. Run after checks.py. Orange = known, blue = ours/computed here, grey = reference."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw={"width_ratios": [1.3, 1]})
tt = np.array(R["flow_t"]); fl = R["flow"]
for t_, f_ in zip(tt, fl):
    ax[0].plot([t_] * len(f_), f_, ".", color=BL, ms=1.2)
ax[0].axhline(0, color=GR, lw=0.8)
for c in R["zero_mode_twists"]:
    ax[0].axvline(c, color=OR, lw=0.8, ls=":")
ax[0].axvline(0, color=GR, lw=1.2)
ax[0].text(0.05, -5.6, "no twist:\ngap 3/2", fontsize=8, color=GR)
ax[0].text(1.05, -5.6, "flat twist:\n2 massless", fontsize=8, color=OR)
ax[0].set_ylim(-5.9, 5.9); ax[0].set_xlabel("twist strength t (0 = round sphere, ±1 = twisted flat)")
ax[0].set_ylabel("Dirac eigenvalue × r")
ax[0].set_title("A. Twisting the hidden 3-sphere moves matter's masses;\nat special twists some become massless")
lab = ["twist t=1\n(S³×S¹)", "twist t=5/3\n(S³×S¹)", "flux 1\n(S²)", "flux 2\n(S²)", "flux 3\n(S²)"]
o = {round(t, 3): (a, b) for t, a, b in R["S3xS1_counts"]}
p = {nf: (a, b) for nf, a, b in R["S2_flux_counts"]}
L_ = [o[1.0][0], o[round(1 + 2 / 3, 3)][0], p[1][0], p[2][0], p[3][0]]
R_ = [o[1.0][1], o[round(1 + 2 / 3, 3)][1], p[1][1], p[2][1], p[3][1]]
xx = np.arange(5)
ax[1].bar(xx - 0.2, L_, 0.4, color=BL, label="left-handed massless")
ax[1].bar(xx + 0.2, R_, 0.4, color=OR, label="right-handed massless")
ax[1].set_xticks(xx); ax[1].set_xticklabels(lab, fontsize=9)
ax[1].set_ylabel("number of massless modes")
ax[1].set_title("B. Twist makes them in left–right pairs;\nflux through a hidden 2-sphere makes one hand only")
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/f14_twist_and_handedness.png", bbox_inches="tight")
print("f14 written")
