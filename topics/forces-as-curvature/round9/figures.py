"""DRAFT picture for round 9. Run after checks.py. Orange = known (geometry/textbook), blue = from the law, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
pf = np.array(R["per_face_sample"])
ax[0].plot(pf[:, 2] * 1e3, pf[:, 1] * 1e3, "o", color=BL, ms=3, label="sampled faces (both charts)")
lim = [pf[:, 2].min() * 1e3 * 0.98, pf[:, 2].max() * 1e3 * 1.02]
ax[0].plot(lim, lim, "-", color=OR, lw=1, label="twist = area")
ax[0].set_xlabel("area of the velocity triangle on the sphere (×10⁻³)")
ax[0].set_ylabel("twist from the law (×10⁻³ rad)")
ax[0].set_title(f"A. q = +1: each of {R['faces']} triangles' twist is its area;\nall together: {R['total_twist']/(2*math.pi):.6f} whole turns")
ax[0].legend(fontsize=8)
d = np.array(R["disc"])
ax[1].plot(d[:, 0], -d[:, 3] / (2 * math.pi), "o-", color=BL, label="q = −1 (relativity): open disc |u| < R")
ax[1].axhline(R["flux_number"], color=OR, lw=2, label="q = +1 (rotations): closed sphere, exactly 2")
ax[1].set_yscale("log"); ax[1].set_xlabel("how far out in velocity space (R; 1 = light speed)")
ax[1].set_ylabel("total twist (whole turns)")
ax[1].set_title("B. Only the closed end of the dial gives\na whole-number flux")
ax[1].legend(fontsize=8)
ch = R["charges"]
xx = np.arange(len(ch))
ax[2].bar(xx - 0.2, [c[2] for c in ch], 0.4, color=BL, label="left-handed massless")
ax[2].bar(xx + 0.2, [c[3] for c in ch], 0.4, color=OR, label="right-handed massless")
ax[2].set_xticks(xx); ax[2].set_xticklabels([f"charge {c[0]:+.1f}\n(flux {c[1]:+d})" for c in ch], fontsize=9)
ax[2].set_ylabel("number of massless modes"); ax[2].set_yticks([0, 1, 2])
ax[2].set_title("C. Matter charged under the twist on the\nhidden velocity sphere comes out one-handed")
ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/f16_flux_from_the_twist.png", bbox_inches="tight")
print("f16 written")
