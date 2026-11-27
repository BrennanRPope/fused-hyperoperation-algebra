"""DRAFT pictures for round 5. Run after checks.py. Orange = known, blue = ours, grey = control."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
rho = R["rho"]

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
# A: the field term
sw = np.array(R["sweep"])
ff = np.linspace(0, sw[:, 0].max() * 1.05, 100)
ax[0].plot(ff, -(rho**2 / 8) * ff**2, color=OR, label="Kaluza–Klein formula −(ρ²/16) FᵃFᵃ")
ax[0].plot(sw[:, 0], sw[:, 1], "o", color=BL, ms=8, label="hidden SU(3) space: from loops (12×12 chain)")
ax[0].plot([R["torus_c1c2"]], [R["torus_R"]], "s", color=GR, ms=8, label="hidden flat 8-torus: from loops (0)")
ax[0].axhline(0, color=GR, lw=0.8, ls=":")
ax[0].set_xlabel("c₁·c₂  (field λ₄ dx, λ₅ dy: only the k·ab part, no slopes)")
ax[0].set_ylabel("12D curvature minus the hidden space's own")
ax[0].set_title("A. Gluons pulling on gluons = the twist of the hidden space\n(¾ of this field lands in colour 8, which SU(2) lacks)")
ax[0].legend(fontsize=8, loc="lower left")
# B: colour factor read off the hidden space's curvature
groups = ["U(1)\ncircle", "SU(2)\n3-sphere", "SU(3)\n8D group space"]
ca_geo = [0.0, 2.0, R["C_A_from_curvature"]]       # U(1): flat circle; SU(2): 6/rho^2 x rho^2/3 (round 4); SU(3): this round
xx = np.arange(3)
ax[1].bar(xx, ca_geo, 0.5, color=[GR, BL, BL])
for i, v in enumerate(ca_geo):
    ax[1].text(i, v + 0.08, f"{v:.3f}" if i == 2 else f"{v:.0f}", ha="center", color=BL if i else GR)
for j, (k, (v, e)) in enumerate(R["lep"].items()):
    ax[1].errorbar(2.33 + 0.1 * j, v, yerr=e, fmt="o", color=OR, capsize=4, label=f"{k} (LEP four-jet)")
ax[1].set_xticks(xx); ax[1].set_xticklabels(groups)
ax[1].set_ylabel("C_A = (hidden curvature) × ρ² / (number of hidden directions)")
ax[1].set_ylim(0, 4)
ax[1].set_title("B. The self-coupling strength is the hidden space's curvature\n(light: flat circle, no self-coupling)")
ax[1].legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig("figures/f12_su3_hidden_space.png", bbox_inches="tight")
print("f12 written")
