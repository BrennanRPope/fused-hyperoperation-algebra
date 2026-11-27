"""Round 10: hyper integrals against the jet layers of topics/gravity-from-theory/jet-river.

jet-river: one step of a 3D flow, det(I + hS) = 1 + hP + h^2 Q + h^3 R (P, Q, R the elementary
symmetric functions of S's eigenvalues). Our rung-2 fused sum is the same fold, prod(1 + h*lambda).
Three comparisons:
 A. Matrix product integral: prod_i (I + h S(x_i)) -> ordered exponential; its det -> exp(int tr S).
    Q and R are O(h^2) per step and vanish in the limit; only P survives.
 B. On primes there is no h -> 0, so every layer survives: prod_p (1 + k p^-s) = sum_m e_m k^m,
    e_m = sum over squarefree n with m prime factors of n^-s.
 C. As a function of the dial k the product has only real zeros, at k = -p^s, so the layers e_m
    satisfy Newton's inequalities e_m^2 >= e_{m-1} e_{m+1} (a Polya frequency sequence).
"""
import json, numpy as np, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm

out = {}
rng = np.random.default_rng(7)

# ---- A. matrix rung-2 product integral -------------------------------------------------
A0, A1, A2 = (rng.normal(size=(3, 3)) * 0.6 for _ in range(3))
S = lambda x: A0 + A1 * math.sin(2 * x) + A2 * x
L = 1.0
def prod_integral(n, k=1.0):
    h = L / n; M = np.eye(3); dets = []
    for i in range(n):
        Si = S((i + 0.5) * h)
        M = (np.eye(3) + k * h * Si) @ M
    return M
def ordered_exp(n=20000):
    h = L / n; M = np.eye(3)
    for i in range(n):
        M = expm(h * S((i + 0.5) * h)) @ M
    return M
ref = ordered_exp(4000)
trint = sum(np.trace(S((i + 0.5) / 20000)) for i in range(20000)) / 20000
ns = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
errM, errdet, layerQ = [], [], []
for n in ns:
    M = prod_integral(n)
    errM.append(float(np.linalg.norm(M - ref)))
    errdet.append(float(abs(np.linalg.det(M) - math.exp(trint))))
    h = L / n
    # size of the h^2 Q + h^3 R part of one step, summed over the steps
    q = 0.0
    for i in range(n):
        Si = S((i + 0.5) * h); ev = np.linalg.eigvals(Si)
        P = ev.sum(); Q = ev[0]*ev[1] + ev[0]*ev[2] + ev[1]*ev[2]; R = ev.prod()
        q += abs(h*h*Q + h**3 * R)
    layerQ.append(float(q))
out["A"] = {"ns": ns, "err_vs_ordered_exp": errM, "err_det_vs_exp_int_trace": errdet,
            "sum_of_Q_R_parts": layerQ, "exp_int_trace": math.exp(trint),
            "det_ref": float(np.linalg.det(ref))}

# ---- B. layers on primes ---------------------------------------------------------------
N = 10**6
sieve = np.ones(N + 1, bool); sieve[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0].astype(float)
MMAX = 10
def layers(s):
    e = np.zeros(MMAX + 1); e[0] = 1.0
    for x in primes ** (-s):
        e[1:] = e[1:] + x * e[:-1]
    return e
z2, z4 = math.pi**2 / 6, math.pi**4 / 90
B = {}
for s in (2.0, 3.0):
    e = layers(s)
    B[str(s)] = {"e": e.tolist(),
                 "sum_k1": float(e.sum()), "sum_km1": float(sum(e[m] * (-1)**m for m in range(MMAX + 1)))}
B["2.0"]["target_k1_zeta2_over_zeta4"] = z2 / z4
B["2.0"]["target_km1_one_over_zeta2"] = 1 / z2
# direct check of layer 2 at s = 2: squarefree n = pq, p < q
pr = primes[primes < 2000]
e2_direct = sum(1/(p*q)**2 for i, p in enumerate(pr) for q in pr[i+1:])
B["e2_direct_primes_below_2000"] = e2_direct
out["B"] = B

# ---- C. Newton / Turan on the layers -----------------------------------------------------
C = {}
for s in ("2.0", "3.0"):
    e = np.array(B[s]["e"])
    ratios = [float(e[m]**2 / (e[m-1] * e[m+1])) for m in range(1, MMAX)]
    C[s] = {"turan_ratios": ratios, "all_ge_1": all(r >= 1 for r in ratios)}
out["C"] = C

json.dump(out, open("jets.json", "w"), indent=1)

# ---- figure ------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
a.loglog(ns, errM, "o-", label="rung-2 matrix product  vs  ordered exponential")
a.loglog(ns, layerQ, "s-", label="total of the h²Q + h³R jet layers over all steps")
a.loglog(ns, errdet, "^-", label="det of the product  vs  exp(∫ tr S)")
a.loglog(ns, [1.0/n for n in ns], "k:", label="1/n")
a.set_xlabel("number of steps n (h = 1/n)")
a.set_title("A. A 3×3 flow S(x): the fused product integral is the ordered exponential.\n"
            "The Q and R jet layers shrink like 1/n; in the limit only P = tr S is left", fontsize=10)
a.legend(fontsize=8)

a = ax[1]
e = np.array(B["2.0"]["e"]); m = np.arange(MMAX + 1)
a.semilogy(m, e, "o-", color="C3", label="layer e_m on primes, s = 2 (all survive)")
a.semilogy([0, 1, 2, 3], [1, 0.9, 0.25, 0.02], "s--", color="C7",
           label="a 3D jet: 1, P, Q, R, then nothing (illustrative sizes)")
a.set_xlabel("layer m")
a.set_title("B. On primes there is no step size to send to 0, so every layer stays:\n"
            "∏(1 + k p⁻²) = Σ e_m k^m, e_m = Σ over squarefree n with m prime factors of n⁻²", fontsize=10)
a.legend(fontsize=8)

a = ax[2]
ks = np.linspace(-30, 1.5, 2000)
P_small = primes[primes < 2000]
vals = np.array([np.prod(1 + k * P_small**-2.0) for k in ks])
a.plot(ks, np.sign(vals) * np.log1p(np.abs(vals)), color="C0")
for p in (2, 3, 5):
    a.axvline(-p * p, color="C1", ls=":", lw=1)
    a.text(-p * p - 0.4, -1.85, f"−{p}²", ha="center", fontsize=9, color="C1")
a.plot([-1, 1], [math.log1p(1/z2), math.log1p(z2/z4)], "ko")
a.annotate("k = −1: 1/ζ(2)", (-1, math.log1p(1/z2)), (-14, -1.1), fontsize=9, arrowprops=dict(arrowstyle="->"))
a.annotate("k = 1: ζ(2)/ζ(4)", (1, math.log1p(z2/z4)), (-14, 1.45), fontsize=9, arrowprops=dict(arrowstyle="->"))
a.axhline(0, color="k", lw=0.6)
a.set_ylim(-2, 2); a.set_xlabel("dial k")
a.set_ylabel("sign · log(1 + |∏|)")
a.set_title("C. In the dial k the prime product has only real zeros, at k = −p^s,\n"
            "so its layers obey Newton's inequalities (Turán holds at every m)", fontsize=10)
plt.tight_layout()
plt.savefig("figures/fig_hyper_vs_jets.png", dpi=110)
print(json.dumps({"A_last": [errM[-1], layerQ[-1], errdet[-1]], "B2": B["2.0"]["sum_k1"],
                  "target": z2/z4, "Bm1": B["2.0"]["sum_km1"], "t_m1": 1/z2,
                  "e2": [B["2.0"]["e"][2], e2_direct], "C": C}, indent=1))
