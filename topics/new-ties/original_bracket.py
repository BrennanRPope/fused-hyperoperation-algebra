"""Round 11 (DRAFT): the certified bracket for the original break (README section 15, part O5).

Lower bound: the reduced pieces whose class-1 runs have at most k copies of lam (k = AMAX of original.py) have no
exact overlaps (Round 10 Steps 1-4 for pairs outside the bad windows, original.py's descent for every bad pair with
a <= k). Hochman (2014), applied to their finite truncations, gives dim K >= s_k, the root of
    W_k(s) = (Lam_k - 1)(Q - 1) + Lam_k (M - 1) nu^s Q,   Lam_k = (1 - L^(k+1))/(1 - L),  L = lam^s, M = 1/(1 - mu^s),
    Q = 1/(1 - nu^s - kappa^s).
Upper bound: s_F, the root of the same weight with Lam = 1/(1 - L) (Round 9 Theorem P10 / Round 10 Step 0).
Both weights decrease in s, so one root each; each is certified by interval evaluation on both sides.
Writes original_bracket.json."""
import json, sys
from mpmath import mp, iv, mpf, log10

K = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
DPS = int(K * 0.49) + 80                     # the gap is about lam^(k s) ~ 10^(-0.481 k)
mp.dps = DPS; iv.dps = DPS

def W(s, k, ar):
    lam, mu, nu, kap = ar.mpf(19) / 200, ar.mpf(1) / 20, ar.mpf(1) / 19, ar.mpf(1) / 10
    L, Mu, N_, K_ = lam ** s, mu ** s, nu ** s, kap ** s
    Lam = 1 / (1 - L) if k is None else (1 - L ** (k + 1)) / (1 - L)
    M = 1 / (1 - Mu); Q = 1 / (1 - N_ - K_)
    return (Lam - 1) * (Q - 1) + Lam * (M - 1) * N_ * Q

def root(k):
    r = mp.findroot(lambda s: W(s, k, mp) - 1, mpf("0.4704834"))
    eps = mpf(10) ** (-(DPS - 30))
    lo, hi = r - eps, r + eps
    assert (W(iv.mpf(lo), k, iv) - 1).a > 0 and (W(iv.mpf(hi), k, iv) - 1).b < 0
    return lo, hi

if __name__ == "__main__":
    kl, kh = root(K)
    fl, fh = root(None)
    gap = fh - kl
    out = dict(k=K, dps=DPS, lower=mp.nstr(kl, 20), upper=mp.nstr(fh, 20), s_F=mp.nstr(fl, 40),
               gap=mp.nstr(gap, 5), log10_gap=float(log10(gap)),
               round10=dict(k=12, lower=0.4704833085, upper=0.4704834402, gap=1.317e-7))
    json.dump(out, open("original_bracket.json", "w"), indent=1)
    print(out)
