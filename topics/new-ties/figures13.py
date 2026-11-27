"""Round 13 pictures (DRAFT): fig41 distance to R against the pair's size, fig42 bad pairs against random stretches,
fig43 where the near misses sit."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

dd = json.load(open("dist13.json")); sf = json.load(open("safe13.json")); cc = json.load(open("control13.json"))
P = dd["per_pair"]
A = np.array([p[0] for p in P]); D = np.array([-p[4] for p in P])          # D = -log2 dist (upper-bound point)
order = np.argsort(A, kind="stable"); A, D = A[order], D[order]
rec = np.maximum.accumulate(D); n_so_far = np.arange(1, len(A) + 1)

# tail exponent s: P(dist < 2^-k) ~ 2^(-s k), fitted on both samples over k in [8, 30]
def tail(vals, ks):
    v = np.array(vals); return np.array([(v >= k).mean() for k in ks])
ks = np.arange(4, 61)
tb = tail(D, ks); tc = tail([-p[4] for p in cc["per_point"]], ks)
def fit(t, lo=8, hi=30, kk=ks):
    m = (kk >= lo) & (kk <= hi) & (t > 0); return np.polyfit(kk[m], np.log2(t[m]), 1)
fb, fc = fit(tb), fit(tc); s_b, s_c = -fb[0], -fc[0]
kf = np.arange(4, 231); tcf = tail([-p[4] for p in cc["per_point"]], kf); tbf = tail(D, kf)
ff = fit(tcf, 40, 150, kf); s_far = -ff[0]; fbf = fit(tbf, 40, 150, kf)
json.dump(dict(s_bad_k8_30=s_b, s_control_k8_30=s_c, s_control_k40_150=s_far, s_bad_k40_150=-fbf[0], intercept_far=ff[1],
               record=float(rec[-1]), pairs=len(A),
               predicted_record=float((math.log2(len(A)) + ff[1]) / s_far)), open("figures13_fit.json", "w"), indent=1)

# fig41
fig, ax = plt.subplots(figsize=(9, 5.2))
ax.scatter(A, D, s=2, color="#1f77b4", alpha=0.25, label="each bad pair: -log2 of its distance to R")
ax.step(A, rec, where="post", color="#d62728", lw=1.8, label="closest approach so far")
pred = (np.log2(n_so_far) + ff[1]) / s_far
ax.plot(A, pred, "--", color="black", lw=1.2, label="random-point prediction from the control's far tail (s = %.3f)" % s_far)
SA = [p[0] for p in sf["per_pair"] if p[3]]; SK = [p[3] for p in sf["per_pair"] if p[3]]
ax.scatter(SA, SK, s=6, marker="x", color="#2ca02c", alpha=0.5, label="safe radius from the range search, a <= 10,000")
ax.set_ylim(bottom=0); ax.set_xscale("log"); ax.set_xlabel("a (the stretch is lam^a / mu^b)"); ax.set_ylabel("-log2 distance to the line-up set R")
ax.set_title("Fig 41 (DRAFT). Bad stretches come ever closer to R, at the pace random points would")
ax.legend(fontsize=8, loc="upper left"); fig.tight_layout(); fig.savefig("figures/fig41-distance-to-R.png", dpi=130); plt.close(fig)

# fig42
fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogy(ks, tb, "o-", ms=3, color="#1f77b4", label="bad pairs, a <= %d (%d pairs), slope s = %.3f" % (dd["amax"], len(A), s_b))
ax.semilogy(ks, tc, "s-", ms=3, color="#ff7f0e", label="random stretches in the same windows (%d), slope s = %.3f" % (cc["n"], s_c))
ax.set_xlabel("k"); ax.set_ylabel("share with distance to R below 2^-k")
ax.set_title("Fig 42 (DRAFT). Bad pairs and random stretches approach R the same way")
ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig("figures/fig42-bad-vs-random.png", dpi=130); plt.close(fig)

# fig43
R = np.array([p[2] for p in P])[order]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, (lo, hi, nm) in zip(axs, [(9 / 19, 100 / 171, "K window"), (0.9, 10 / 9, "N window")]):
    m = (R >= lo) & (R <= hi)
    ax.scatter(R[m], D[m], s=2, alpha=0.3, color="#1f77b4")
    ax.set_xlabel("stretch rho (%s)" % nm)
    if nm == "K window": ax.axvline(10 / 19, color="#d62728", lw=1, ls=":"); ax.text(10 / 19, ax.get_ylim()[1] * 0.95, " 10/19 (N = (10/19) K)", color="#d62728", fontsize=8, va="top")
    else: ax.axvline(18 / 19, color="#d62728", lw=1, ls=":"); ax.text(18 / 19, ax.get_ylim()[1] * 0.95, " 18/19", color="#d62728", fontsize=8, va="top")
axs[0].set_ylabel("-log2 distance to R")
fig.suptitle("Fig 43 (DRAFT). Where the near misses sit in the two windows")
fig.tight_layout(); fig.savefig("figures/fig43-near-misses.png", dpi=130); plt.close(fig)
print(open("figures13_fit.json").read())
