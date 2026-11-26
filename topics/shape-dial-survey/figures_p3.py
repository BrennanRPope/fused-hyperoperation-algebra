"""Figure 21: Lee-Yang at p = 3, proved (run after p3_certificate.py)."""
import json
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = json.load(open("p3_certificate.json"))
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
mp.mp.dps = 40
def coeffs(c, K=140):
    c = mp.mpf(c); return [(-1)**k/mp.factorial(2*k)*mp.beta(mp.mpf(2*k + 1)/3, c + 1)/3 for k in range(K)]
def ev(cf, s):
    s2 = mp.mpf(s)**2; acc = mp.mpf(0)
    for ck in reversed(cf): acc = acc*s2 + ck
    return acc
fig, axs = plt.subplots(2, 2, figsize=(15.5, 11), gridspec_kw=dict(wspace=0.24, hspace=0.38))
# A: the climb
ax = axs[0, 0]
ax.add_patch(plt.Rectangle((0, 0.2), 1, 0.6, color=GREEN, alpha=0.35))
ax.add_patch(plt.Rectangle((1, 0.2), 1, 0.6, color=BLUE, alpha=0.30))
ax.add_patch(plt.Rectangle((2, 0.2), 1, 0.6, color=BLUE, alpha=0.55))
ax.add_patch(plt.Rectangle((3, 0.2), 1.2, 0.6, color=RED, alpha=0.25))
ax.text(0.5, 0.5, "Pólya\n(1918)", ha="center", va="center", fontsize=11)
ax.text(1.5, 0.5, "step 1\nb in (0, 1]", ha="center", va="center", fontsize=11)
ax.text(2.5, 0.5, "step 2\nb in (1, 2]", ha="center", va="center", fontsize=11, color="white")
ax.text(3.6, 0.5, "breaks\n(proved)", ha="center", va="center", fontsize=11)
for x0 in (0.5, 1.5):
    ax.annotate("", xy=(x0 + 1, 0.92), xytext=(x0, 0.92), arrowprops=dict(arrowstyle="->", lw=1.6, connectionstyle="arc3,rad=-0.4"))
ax.text(1.5, 1.22, "each step: G at b + 1 alternates at the zeros of G at b,\nthen Krein's theorem leaves no room for complex zeros",
        ha="center", fontsize=9.5)
ax.axvline(3, color="k", lw=1.5); ax.text(3.02, 0.08, "a = 3, n = 10 = 1 + p²", fontsize=9.5)
ax.set_xlim(-0.1, 4.3); ax.set_ylim(0, 1.6); ax.set_yticks([]); ax.set_xlabel("a  (ball dimension n = 3a + 1)")
ax.set_title("p = 3: real zeros for every a ≤ 3, in two whole steps\nup from Pólya's case")
# B: zeros of G_b against b, coloured by the sign of G at b + 1
ax = axs[0, 1]; bs = np.round(np.arange(0.0, 2.0001, 0.1), 3)
for b in bs:
    cb, ca = coeffs(b), coeffs(b + 1); grid = np.arange(0.5, 40, 0.05); v = [ev(cb, t) for t in grid]
    zs = [float(mp.findroot(lambda x: ev(cb, x), (mp.mpf(grid[i]), mp.mpf(grid[i + 1])), solver="illinois", verify=False))
          for i in range(len(v) - 1) if v[i]*v[i + 1] < 0]
    for z in zs:
        sg = ev(ca, z) > 0; ax.scatter(b, z, color=BLUE if sg else ORANGE, s=14)
ax.axhspan(21, 40, color=PURPLE, alpha=0.08); ax.axhline(25, color=GREY, lw=1, ls="--"); ax.text(2.04, 25, "boxes\nreach 25", fontsize=8.5, color=GREY, va="center")
ax.text(1.0, 31, "far zone (s ≥ 21): explicit inequalities", ha="center", color=PURPLE, fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.9))
ax.text(1.0, 13.5, "near zone: every zero boxed,\n200 boxes in b", ha="center", fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.9))
ax.scatter([], [], color=BLUE, s=14, label="G at b + 1 > 0 there"); ax.scatter([], [], color=ORANGE, s=14, label="G at b + 1 < 0 there")
ax.set_xlabel("b"); ax.set_ylabel("zeros of G at b"); ax.set_ylim(0, 40); ax.legend(fontsize=9, loc="lower right")
ax.set_title("Every zero of G at b, with the sign of G at b + 1:\ncolours alternate up every column")
# C: the mechanism at b = 2 (a = 3)
ax = axs[1, 0]; cb, ca = coeffs(2.0), coeffs(3.0); ts = np.linspace(2, 34, 1200)
gb = np.array([float(ev(cb, t)) for t in ts]); ga = np.array([float(ev(ca, t)) for t in ts])
ax.plot(ts, gb*ts**3/18, color=GREY, lw=1.6, label="G at b = 2 (scaled by s³/18)")
ax.plot(ts, ga*ts**4/162, color=BLUE, lw=2, label="G at a = 3 (scaled by s⁴/162)")
zi = [i for i in range(len(gb) - 1) if gb[i]*gb[i + 1] < 0]
ax.scatter(ts[zi], ga[zi]*ts[zi]**4/162, color=RED, s=40, zorder=3, label="G at a = 3, at the zeros of G at 2")
ax.axhline(0, color="k", lw=0.8); ax.axhline(1/9, color=RED, lw=0.8, ls=":"); ax.axhline(-1/9, color=RED, lw=0.8, ls=":")
ax.text(2.2, 1.55, "dotted: centre's kink, ±1/9", fontsize=9, color=RED)
ax.set_xlabel("s"); ax.set_ylim(-1.6, 1.9); ax.legend(fontsize=9, loc="upper right", ncol=1)
ax.set_title("Why it alternates: at each zero of G at b, G at b + 1 sits near a crest,\nclear of the kink's ±1/9 band")
# D: certified far-zone margins over b at s = 21
ax = axs[1, 1]
bm = [np.mean(f["b"]) for f in D["far"]["boxes"]]
for k, (lab, col) in enumerate([("|ε| < 1 (rim wave beats the kink)", GREEN), ("one simple zero per half-turn", BLUE),
                                ("G at b + 1 keeps its predicted sign", ORANGE)]):
    ax.plot(bm, [f["margins"][k] for f in D["far"]["boxes"]], "o-", color=col, ms=3, lw=1.6, label=lab)
ax.axhline(0, color="k", lw=1.2); ax.set_ylim(-0.1, 1.1); ax.set_xlabel("b (40 interval boxes)")
ax.set_ylabel("certified margin at s = 21 (holds for all s ≥ 21)")
ax.set_title("Far zone: all three inequalities hold with room to spare,\nproved with interval arithmetic")
ax.legend(fontsize=9, loc="lower left")
fig.suptitle("Figure 21. Lee–Yang at roundness 3 is proved for every a ≤ 3 (DRAFT, computer-assisted)", fontsize=14, y=0.995)
fig.savefig("figures/fig21-lee-yang-p3-proof.png", dpi=110, bbox_inches="tight")
print("saved")
