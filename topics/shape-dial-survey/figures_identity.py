"""Figure 18: the step identity (run after identity_check.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
exec(open("interlace_check.py").read().split("if __name__")[0])
K = json.load(open("identity_check.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(18.5, 5.4), gridspec_kw=dict(wspace=0.28, width_ratios=[1.05, 1.1, 1]))
# A: Rolle picture at p = 3, a = 1.5 -> 2.5
ax = axs[0]; p, a = 3.0, 1.5; kap = (a + 2)*p + 1 - p   # kappa for the pair (a, a+1) = (a+1)p+1
kap = (a + 1)*p + 1
s = np.linspace(0.05, 26, 1200)
cb = coeffs(p, a + 1); ca = coeffs(p, a)
H = np.array([float(mp.mpf(x)**kap*ev(cb, x)) for x in s]); g = np.array([float(ev(ca, x))*x**(a + 1) for x in s])
H = H/np.max(np.abs(H)); w = np.sign(H)*np.abs(H)**(1/kap*6)                     # shrink the growth so all humps show
ax.plot(s, w, color=ORANGE, lw=2, label=f"s^κ G at a+1 = {a+1:g} (growth squashed)")
ax.plot(s, g/np.max(np.abs(g[s > 5])), color=BLUE, lw=2, label=f"G at a = {a:g} (scaled)")
za, _ = zeros(p, a, polish=True); zb, _ = zeros(p, a + 1, polish=True)
for z in za: ax.axvline(z, color=BLUE, lw=.7, ls=":")
ax.scatter(zb, [0]*len(zb), color=ORANGE, s=36, zorder=4, edgecolor="white"); ax.scatter(za, [0]*len(za), color=BLUE, s=36, zorder=4, edgecolor="white")
ax.axhline(0, color="k", lw=.6); ax.set_ylim(-1.5, 1.7); ax.set_xlabel("s")
ax.set_title("p = 3: each hump of the orange curve\nhas its top at a blue zero (Rolle)")
ax.legend(fontsize=9, loc="upper left")
# B: the reduction map
ax = axs[1]; P = np.linspace(1, 7, 601)
ax.fill_between(P, 0, 1, color=BLUE, alpha=.35, lw=0)
ax.fill_between(P, np.maximum(1, P - 1), P, color="#555555", alpha=.45, lw=0)
ax.fill_between(P, 1, np.maximum(1, P - 1), color=GREY, alpha=.15, lw=0, hatch="//", edgecolor="#999999")
ax.fill_between(P, P, 8, color=RED, alpha=.15, lw=0)
for pe in [2, 4, 6]:
    ax.plot([pe, pe], [0, 8], color=GREEN if pe == 2 else "#9ae6b4", lw=4 if pe == 2 else 3, solid_capstyle="butt")
ax.scatter([4]*8, np.arange(1, 9), s=40, color=GREEN, zorder=4, edgecolor="white")
for k in range(3):
    ax.annotate("", xy=(5.6, 4.6 - k), xytext=(5.6, 5.4 - k), arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
ax.text(5.72, 3.2, "step down\n(proved)", fontsize=9.5)
ax.plot(P, P, color="k", lw=1); ax.set_xlim(1, 7); ax.set_ylim(0, 8)
ax.set_xlabel("roundness p"); ax.set_ylabel("a = (n − 1)/p"); ax.set_title("What's left: one strip of width 1 under a = p")
ax.legend(handles=[Patch(color=BLUE, alpha=.35, label="holds, proved (a ≤ 1)"),
                   Patch(color="#555555", alpha=.45, label="open: the top strip p − 1 < a ≤ p"),
                   Patch(facecolor=GREY, alpha=.25, hatch="//", edgecolor="#999999", label="follows once the strip is proved"),
                   Patch(color=RED, alpha=.15, label="breaks, proved (p not even)"),
                   Line2D([], [], color=GREEN, lw=4, label="proved: p = 2, p = 4 whole a")],
          loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=8.6, frameon=False)
# C: why the step can't climb: F at the zeros of G_a
ax = axs[2]
for key, col, off, lab in [("p=3, a=1 -> 2", BLUE, -0.18, "p = 3, a = 1 → 2 (a+1 ≤ p): signs alternate"),
                           ("p=1.5, a=1 -> 2", RED, 0.18, "p = 1.5, a = 1 → 2 (a+1 > p): one sign")]:
    F = np.array(K["lobes"][key]["F_at_zeros"]); k = np.arange(1, len(F) + 1)
    v = np.sign(F)*np.log10(1 + np.abs(F)/np.min(np.abs(F)))
    ax.bar(k + off, v, width=0.34, color=col, label=lab)
ax.axhline(0, color="k", lw=.6); ax.set_xlabel("k-th zero of G at a"); ax.set_ylabel("sign · log size of F there")
ax.set_title("Climbing a → a + 1 needs F to flip sign\nat every zero of G at a")
ax.legend(fontsize=9, loc="lower left")
fig.savefig("figures/fig18-lee-yang-step-identity.png", dpi=130, bbox_inches="tight"); print("wrote fig18")
