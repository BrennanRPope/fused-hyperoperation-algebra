"""Figure 19: the top strip through the shape dial (run after strip_check.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
exec(open("strip_check.py").read().split("t0 = time.perf_counter()")[0])
SC = json.load(open("strip_check.json"))
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(18.5, 5.3), gridspec_kw=dict(wspace=0.27))
ax = axs[0]
for p, col, targets in [(2.5, GREEN, [8, 15, 25, 40, 70]), (3.0, BLUE, [8, 15, 25, 40, 70]),
                        (5.0, ORANGE, [15, 25, 40, 70, 110]), (7.0, RED, [25, 50, 80, 110, 160, 220])]:
    a = p; rs = [crest(a, t) for t in targets]; ms = [margin(p, a, r, n=180) for r in rs]
    ax.plot(rs, ms, "o-", color=col, lw=2, label=f"p = {p:g}, a = p")
    ax.axhline(abs(np.sin(np.pi*p/2))/p**(p - 1), color=col, lw=1, ls=":")
ax.axhline(1, color="k", lw=1.2); ax.text(5.2, 1.04, "below this line the ball comparison works", fontsize=9.5)
ax.set_xscale("log"); ax.set_ylim(0, 2.2); ax.set_xlabel("radius of the circle (log scale)")
ax.set_ylabel("size of G_p − ball, relative to the ball")
ax.set_title("Far out every roundness looks like the ball\n(dotted: where each curve settles, the centre's kink)")
ax.legend(fontsize=9, loc="upper right")
ax = axs[1]; p, a, R = 5.0, 5.0, SC["far_radius"]["p=5, a=5"]
zp, _ = real_zeros(p, a, R); z2, _ = real_zeros(2.0, a, R)
ax.scatter(zp, [1]*len(zp), color=ORANGE, s=30, label=f"roundness 5: {len(zp)} zeros")
ax.scatter(z2, [0]*len(z2), color=BLUE, s=30, label=f"ball, same a: {len(z2)} zeros")
for u, v in zip(zp, z2): ax.plot([u, v], [1, 0], color=GREY, lw=.6)
ax.set_yticks([0, 1]); ax.set_yticklabels(["ball (p = 2)", "p = 5"]); ax.set_ylim(-0.6, 1.8)
ax.set_xlabel("s"); ax.axvline(R, color="k", lw=1); ax.text(R - 1, 1.55, "far zone →", ha="right", fontsize=10)
ax.set_title(f"a = 5: the zeros differ near the centre but\nmatch one for one; past r = {R:.0f} the ball takes over")
ax.legend(fontsize=9, loc="lower left")
ax = axs[2]
for key, col in [("p=3, a=3", BLUE), ("p=5, a=5", ORANGE), ("p=7, a=7", RED)]:
    L = SC["lobes"][key]["log10_abs_L"]; ax.plot(range(1, len(L) + 1), L, "o-", ms=3.5, color=col, lw=1.6, label=key.replace("p=", "p = ").replace("a=", "a = "))
ax.set_xlabel("lobe number k"); ax.set_ylabel("log10 |weighted lobe of G at a − 1|")
ax.set_title("Peeling one shell: the weighted lobes grow\nfrom the very first, so each gap gets a zero")
ax.legend(fontsize=9, loc="lower right")
fig.savefig("figures/fig19-lee-yang-top-strip.png", dpi=130, bbox_inches="tight"); print("wrote fig19")
