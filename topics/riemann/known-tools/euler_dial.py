"""The k dial on the Euler product, F_k(s) = prod_p (1 + k p^-s) = zeta(s)^k G_k(s), G_k holomorphic, nonzero
for Re s > 1/2 when |k| < sqrt 2. Three uses on RH:
 (B1) the singularity each zero of zeta gives F_k: zero (k = 1), pole (k = -1), branch point (k not an integer);
 (B2) k = -1 is the Mertens function M(x) = sum mu(n); RH <=> M(x) = O(x^(1/2+e)); computed to 1e8;
      k = 1 is the squarefree count, whose error is blind to the zeros of zeta(s) (it sees zeta(2s));
 (C)  the second layer of the Euler product (prime squares) gives -1/2 li(sqrt x) in pi(x): the one finite-step
      term at the RH scale. (pi(x) - li(x)) log x / sqrt x has mean -1.
"""
import json, math
import numpy as np, mpmath as mp
from scipy.special import expi
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 10**8
out = {}
# --- primes and mu to 1e8 ---------------------------------------------------------------
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
primes = np.nonzero(s)[0]
del s
mu = np.ones(N + 1, np.int8); mu[0] = 0
for p in primes:
    mu[p::p] *= -1
    if p * p <= N:
        mu[p*p::p*p] = 0
xs = np.unique(np.logspace(1, 8, 4000).astype(np.int64))
M = np.empty(len(xs)); Q = np.empty(len(xs))
run_m = 0; run_q = 0; last = 0
for i, x in enumerate(xs):
    seg = mu[last + 1:x + 1]
    run_m += int(seg.sum(dtype=np.int64)); run_q += int(np.count_nonzero(seg)); last = x
    M[i] = run_m; Q[i] = run_q
# maximum of |M|/sqrt x over all x (exact, chunked)
mx = 0.0; arg = 0; run = 0
for a in range(1, N + 1, 10**7):
    b = min(N, a + 10**7 - 1)
    c = np.cumsum(mu[a:b + 1], dtype=np.int64) + run
    r = np.abs(c) / np.sqrt(np.arange(a, b + 1))
    k = int(np.argmax(r[(np.arange(a, b + 1) >= 1000) - 0 == 1] if False else r))
    if a + k >= 1000 and r[k] > mx: mx, arg = float(r[k]), a + k
    run = int(c[-1])
del mu
out["mertens"] = {"M_at_1e8": float(M[-1]), "max_abs_M_over_sqrt_x_from_1e3": mx, "at_x": arg}
q_err = Q - 6 * xs / math.pi**2
out["squarefree"] = {"err_at_1e8": float(q_err[-1]), "max_abs_err_over_x_quarter_above_1e3":
                     float(np.max(np.abs(q_err[xs >= 1000]) / xs[xs >= 1000]**0.25))}
# --- (C) second layer ---------------------------------------------------------------------
pc = np.searchsorted(primes, xs, side="right")
li = lambda x: expi(np.log(x))
lx = np.log(xs.astype(float)); sq = np.sqrt(xs.astype(float))
A = (pc - li(xs.astype(float))) * lx / sq
Bc = (pc - li(xs.astype(float)) + 0.5 * li(sq)) * lx / sq
sel = xs >= 10**4
# mean in log x (each decade weighted equally)
w = np.gradient(lx)[sel]
out["second_layer"] = {"mean_A_logweighted_1e4_1e8": float(np.sum(A[sel] * w) / w.sum()),
                       "mean_B_logweighted_1e4_1e8": float(np.sum(Bc[sel] * w) / w.sum()),
                       "A_ever_positive_1e4_1e8": bool((A[sel] > 0).any())}
print(json.dumps(out, indent=1))

# --- (B1) the phase of zeta^k around the first zero ----------------------------------------
mp.mp.dps = 15
re = np.linspace(0.2, 0.8, 90); im = np.linspace(13.4, 14.9, 90)
Zt = np.array([[complex(mp.zeta(mp.mpc(a, b))) for a in re] for b in im])
json.dump(out, open("euler_dial.json", "w"), indent=1)

fig = plt.figure(figsize=(17, 10))
for i, (k, lab) in enumerate([(1, "k = 1: ζ(s)/ζ(2s)\nthe zero of ζ stays a zero"),
                               (-1, "k = −1: 1/ζ(s)\nthe zero becomes a pole"),
                               (0.5, "k = ½: ζ(s)^½ G(s)\nthe zero becomes a branch point")]):
    a = fig.add_subplot(2, 3, i + 1)
    ph = np.angle(np.exp(k * np.log(Zt)))   # principal branch: the cut shows
    a.imshow(ph, origin="lower", extent=[re[0], re[-1], im[0], im[-1]], aspect="auto", cmap="hsv")
    a.axvline(0.5, color="k", lw=0.6, ls=":")
    a.plot(0.5, 14.134725, "wo", ms=6, mec="k")
    a.set_title(lab, fontsize=10); a.set_xlabel("Re s"); a.set_ylabel("Im s")
fig.text(0.5, 0.965, "B1. The k dial on the Euler product, ∏(1 + k p⁻ˢ) = ζ(s)^k G_k(s): colour is the phase near the first zero "
         "(white dot). For every k in (−√2, √2) except 0 and 1, a zero off the line would be a singularity,\n"
         "so RH ⇔ this fused Euler product continues to Re s > ½ cut only along [½, 1]. k = 1 is the one blind setting.",
         ha="center", va="top", fontsize=10.5)
a = fig.add_subplot(2, 3, 4)
a.semilogx(xs, M / np.sqrt(xs), lw=0.6, label="M(x)/√x  (k = −1)")
a.axhline(0, color="k", lw=0.5); a.axhline(1, color="C3", ls=":", lw=0.8); a.axhline(-1, color="C3", ls=":", lw=0.8)
a.set_title(f"B2. k = −1 is the Mertens function. RH ⇔ M(x) = O(x^(½+ε)) (Littlewood 1912).\n"
            f"Up to 10⁸, |M(x)|/√x ≤ {mx:.3f} (from x = 10³)", fontsize=10)
a.set_xlabel("x"); a.legend(fontsize=8)
a = fig.add_subplot(2, 3, 5)
a.semilogx(xs[xs >= 10], (q_err / xs**0.25)[xs >= 10], lw=0.6, color="C2", label="(Q(x) − 6x/π²)/x^¼  (k = 1)")
a.axhline(0, color="k", lw=0.5)
a.set_title("B2. k = 1 is the squarefree count Q(x). Its error is set by the zeros of ζ(2s),\n"
            "a quarter-power scale; the zeros of ζ(s) itself do not show here", fontsize=10)
a.set_xlabel("x"); a.legend(fontsize=8)
a = fig.add_subplot(2, 3, 6)
a.semilogx(xs[xs >= 100], A[xs >= 100], lw=0.6, label="(π(x) − li(x))·log x/√x")
a.semilogx(xs[xs >= 100], Bc[xs >= 100], lw=0.6, label="same with ½ li(√x) added back")
a.axhline(-1, color="C0", ls="--", lw=0.8); a.axhline(0, color="C1", ls="--", lw=0.8)
a.set_title("C. The second layer of the Euler product (prime squares) is −½ li(√x):\n"
            f"the one finite-step term at the RH scale. Log-mean 10⁴–10⁸: {out['second_layer']['mean_A_logweighted_1e4_1e8']:.2f} → "
            f"{out['second_layer']['mean_B_logweighted_1e4_1e8']:.2f}", fontsize=10)
a.set_xlabel("x"); a.legend(fontsize=8)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("figures/fig_euler_dial.png", dpi=105)
