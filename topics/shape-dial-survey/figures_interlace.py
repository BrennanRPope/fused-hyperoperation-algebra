"""Figure 17: the sliding route for the open band (run after interlace_check.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
exec(open("interlace_check.py").read().split("if __name__")[0])
R = json.load(open("interlace_check.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw=dict(wspace=0.27))
for ax, p, title in [(axs[0], "3.0", "p = 3: zeros slide right, never meet"), (axs[1], "1.5", "p = 1.5: zeros meet just past a = p")]:
    T = R["tracks"][p]; P = float(p)
    for a, z in zip(T["a"], T["zeros"]):
        col = BLUE if a <= 1 else (GREY if a <= P else RED)
        ax.scatter([a]*len(z), z, s=7, color=col)
    ax.axvspan(1, P, color=GREY, alpha=.12, lw=0)
    ax.axvline(P, color="k", lw=1); ax.text(P + 0.03, 1.2, "a = p", fontsize=10)
    ax.set_xlabel("a = (n − 1)/p"); ax.set_ylabel("real zeros of G (0 < s ≤ 30)"); ax.set_ylim(0, 31)
    ax.set_title(title)
axs[0].text(1.05, 1.2, "open band (grey): no zeros meet", fontsize=9.5, color="#444")
axs[1].text(1.55, 15, "pairs meet and\nleave the line\n(far ones first)", color=RED, fontsize=10)
# right: alternation at p = 4, a = 2.5 and 3.5
ax = axs[2]; s = np.linspace(0.3, 30, 900)
for a, col in [(2.5, BLUE), (3.5, ORANGE)]:
    c = coeffs(4.0, a); v = np.array([float(ev(c, x)) for x in s])*s**(a + 1)
    ax.plot(s, v/np.max(np.abs(v[s > 8])), color=col, lw=2, label=f"a = {a:g}")
    z, _ = zeros(4.0, a, polish=True); ax.scatter(z, [0]*len(z), color=col, s=40, zorder=4, edgecolor="white")
ax.axhline(0, color="k", lw=.6); ax.set_ylim(-1.6, 1.6); ax.set_xlabel("s")
ax.set_ylabel("G(s)·s^(a+1), scaled"); ax.set_title("p = 4: zeros of a and a + 1 alternate")
ax.legend(loc="upper right", fontsize=10)
fig.savefig("figures/fig17-lee-yang-sliding-route.png", dpi=130, bbox_inches="tight"); print("wrote fig17")
