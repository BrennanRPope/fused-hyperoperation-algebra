"""fig5: the heat-type bounds as jets. Run: python3 jet_figure.py"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
etaC = np.linspace(0, 0.999, 400)
gs = [0, 0.25, 0.5, 0.75, 1.0]
cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(gs)))
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
for g, c in zip(gs, cols):
    e = etaC / (2 - g * etaC)
    axs[0].plot(etaC, e, color=c, lw=2, label=f"γ = {g}")
    axs[1].plot(etaC, e / (1 + (g - 1) * e), color=c, lw=6 - 4 * g, alpha=0.8)
    axs[2].plot(etaC, 1 + (g - 1) * e, color=c, lw=2)
axs[0].set_title("Scalar view: five different bounds\nη* = η_C / (2 − γη_C)", fontsize=11)
axs[0].legend(fontsize=9); axs[0].text(0.03, 0.72, "γ = 1: upper bound\n(velocity point)\nγ = 0: η_C/2\n(parity point)", fontsize=9)
axs[1].set_title("Jet view, velocity layer: all five\nland on one curve tanh(β/2)  (J1)", fontsize=11)
axs[1].text(0.05, 0.75, "u = η* / (1 + (γ − 1)η*)\n= η_C / (2 − η_C)\n= tanh(β/2), β = −ln(1 − η_C)", fontsize=9)
axs[2].set_title("Jet view, scale layer: γ lives only here\ns = 1 + (γ − 1)η*", fontsize=11)
axs[2].text(0.05, 0.62, "stages with different γ multiply\nas commuting jets: no turn (J3)", fontsize=9)
for a in axs:
    a.set_xlabel("Carnot efficiency η_C"); a.set_xlim(0, 1); a.grid(alpha=0.3)
axs[0].set_ylabel("efficiency at max power η*"); axs[0].set_ylim(0, 1); axs[1].set_ylim(0, 1); axs[2].set_ylim(0, 1.05)
fig.suptitle("Fig 5. The heat-type family as one velocity carried with different scales", fontsize=12)
fig.tight_layout(); fig.savefig("fig5-jet-view.png", dpi=150)
print("wrote fig5-jet-view.png")
