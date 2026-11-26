"""Figure 20: the lobe law far out (run after lobe_law_check.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
exec(open("lobe_law_check.py").read().split("t0 = time.perf_counter()")[0])
D = json.load(open("lobe_law_check.json"))
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(2, 2, figsize=(15.5, 11), gridspec_kw=dict(wspace=0.22, hspace=0.36))
# A: the zigzag closes onto 1 +/- kink share at a = p
ax = axs[0, 0]
for (key, col) in [("p=2.5", GREEN), ("p=3", BLUE), ("p=3.5", ORANGE)]:
    r = np.array(D["edge_rows"][key]); p = float(key[2:]); m = abs(np.sin(np.pi*p/2))/p**(p - 1)
    ax.plot(r[:, 0], r[:, 3], "-", color=col, lw=0.8, alpha=0.6)
    ax.scatter(r[:, 0], r[:, 3], color=col, s=12, label=f"{key}, a = p: lines at 1 ± {m:.3f}")
    for y in (1 - m, 1 + m): ax.axhline(y, color=col, lw=1, ls=":")
ax.axhline(0, color="k", lw=1.2)
ax.text(60, 0.12, "below 0 a lobe would shrink", fontsize=9.5)
ax.set_ylim(-0.2, 3.2); ax.set_xlim(0, 220)
ax.set_xlabel("zero of G at a − 1 (where the lobe ends)")
ax.set_ylabel("lobe excess ÷ the weight's share\n(L_(k+1)/L_k − 1) · z / (π a (p−1))")
ax.set_title("Lobes grow by the weight's share, give or take the centre's kink\n(zigzag closes onto the dotted lines, never below)")
ax.legend(fontsize=9, loc="upper right")
# B: why it zigzags: the integrand is a rim wave riding on the centre's constant-sign offset
ax = axs[0, 1]; p, a = 2.5, 2.5; b = a - 1
c = coeffs(p, b, 60); ts = np.linspace(4, 40, 900)
A = p**b*float(mp.gamma(a)); cc = b*float(mp.gamma(p + 1))*np.sin(np.pi*p/2)
w = np.array([float(ev(c, t)) for t in ts])*ts**a/A
ax.plot(ts, w, color=BLUE, lw=1.8, label="G at a − 1, scaled by its rim size")
ax.plot(ts, cc*ts**(a - p - 1)/A, color=RED, lw=2, label="centre kink offset (one sign)")
ax.plot(ts, np.cos(ts - np.pi*a/2), color=GREY, lw=1, ls="--", label="pure rim wave")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("t"); ax.set_ylim(-1.4, 1.4)
ax.set_title("Why: the centre's kink nudges the wave one way,\nso one sign of lobe is a bit bigger (p = 2.5, a = 2.5)")
ax.legend(fontsize=9, loc="lower right")
# C: past the edge the zigzag opens up and a lobe shrinks
ax = axs[1, 0]
for rows, lab, col in [(D["inside_rows"]["p=3, a=2.5"], "a = 2.5 (inside)", GREEN), (D["edge_rows"]["p=3"], "a = 3 (edge)", BLUE),
                       (D["past_edge"]["rows"], "a = 4.5 (past the edge)", RED)]:
    r = np.array(rows); r = r[r[:, 0] < 120]
    ax.plot(r[:, 0], r[:, 3], "o-", color=col, ms=3, lw=0.8, label=lab)
ax.axhline(0, color="k", lw=1.2); ax.set_ylim(-2.5, 4); ax.set_xlim(0, 120)
ax.set_xlabel("zero of G at a − 1"); ax.set_ylabel("lobe excess ÷ the weight's share")
ax.set_title("p = 3: inside the strip the zigzag closes; past a = p it opens\nand lobes start to shrink, just where the rule breaks")
ax.legend(fontsize=9, loc="upper right")
# D: the law's floor at a = p is the same kink share as the ball comparison
ax = axs[1, 1]; ps = np.linspace(1.02, 7, 600); m = np.abs(np.sin(np.pi*ps/2))/ps**(ps - 1)
ax.plot(ps, 1 - m, color=PURPLE, lw=2.2, label="floor of the lobe law at a = p:  1 − |sin(πp/2)| / p^(p−1)")
for g in D["grid"]:
    if abs(g["a"] - g["p"]) < 1e-9: ax.scatter(g["p"], g["min_D"], color=ORANGE, s=40, zorder=3)
ax.scatter([], [], color=ORANGE, s=40, label="smallest excess seen, all lobes (a = p)")
for pe in [2, 4, 6]: ax.scatter(pe, 1, color="k", marker="s", s=30, zorder=3)
ax.scatter([], [], color="k", marker="s", s=30, label="even p: no kink, floor exactly 1")
ax.axhline(0, color="k", lw=1.2); ax.set_ylim(-0.1, 1.7); ax.set_xlabel("roundness p")
ax.set_ylabel("lobe excess ÷ the weight's share")
ax.set_title("The floor stays above 0 for every p > 1, and every lobe\nchecked sits above it")
ax.legend(fontsize=9, loc="lower right")
fig.suptitle("Figure 20. The lobe law: far out, weighted lobes grow at every roundness up to a = p (DRAFT)", fontsize=14, y=0.995)
fig.savefig("figures/fig20-lee-yang-lobe-law.png", dpi=110, bbox_inches="tight")
print("saved")
