"""Figure 27: Lee-Yang at p = 4 above a = 4 reduces to (1 - x^2)^c (1 + x^2)^theta with theta in [0, 1)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = json.load(open("even_p_check.json"))["E4"]
BLUE, ORANGE, GREY = "#2b6cb0", "#dd6b20", "#888888"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 6.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.32)
ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
txt = ("At p = 4, with c = m + θ (m whole, 0 ≤ θ < 1):\n\n"
       "   (1 − x⁴)^c = (1 − x²)^c (1 + x²)^θ · (1 + x²)^m\n\n"
       "Each factor 1 + x² turns into the operator 1 − D²,\n"
       "which keeps real zeros (Laguerre). So:\n\n"
       "   real zeros for (1 − x²)^c (1 + x²)^θ\n"
       "        ⇒ real zeros at p = 4, step c.\n\n"
       "θ = 0 is the round ball (Bessel): the whole-a case.\n"
       "θ = 1 is one operator step from it: real again.\n"
       "Open: every θ strictly between, for every c > 4.\n\n"
       "No such reduction at p = 6 or 8: their extra factors\n"
       "(1 + x² + x⁴, 1 + x⁴) give operators with non-real roots.")
ax.text(0, 0.98, txt, va="top", fontsize=10.3, family="DejaVu Sans")
ax.set_title("The reduction (proved)", fontsize=12.5, loc="left")
th = np.linspace(0, 1, 21)
for k, alpha in enumerate(["8", "20"]):
    ax = fig.add_subplot(gs[0, 1 + k])
    Z = np.array(D[alpha]["zeros"])
    cmap = plt.get_cmap("viridis")
    for j in range(Z.shape[1]):
        ax.plot(th, Z[:, j] - Z[0, j], color=cmap(j/Z.shape[1]), lw=1.8, label=f"zero {j + 1}" if j in (0, Z.shape[1] - 1) else None)
    ax.axhline(0, color=GREY, lw=0.8)
    ax.set_xlabel("θ (fractional part of c)"); ax.set_ylabel("how far each zero has moved since θ = 0")
    ax.legend(fontsize=9.5, loc="lower left")
    ax.set_title(f"c = {alpha}: each of the first {Z.shape[1]} zeros slides down,\nby at most {abs(Z[-1] - Z[0]).max():.1f}; "
                 f"neighbours stay ≥ {D[alpha]['min_gap']:.2f} apart", fontsize=11.5)
fig.suptitle("Figure 27. Lee–Yang at even roundness above a = p: p = 4 reduces to one cleaner family (DRAFT)", fontsize=14, y=1.0)
fig.savefig("figures/fig27-lee-yang-even-p.png", dpi=110, bbox_inches="tight"); print("saved")
