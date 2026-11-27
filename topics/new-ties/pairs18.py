"""Round 18 (DRAFT): Theorem 18.1, no two different reduced three-letter words share a fingerprint. Computer-assisted.

Words and notation as in pairs17.py. A reduced three-letter word is W = N u1 x2 u2 x3, u_i = l^L_i m^M_i, H_i = L_i + M_i,
x_i = K only if M_i = 0, prefix N u1 x2 != NN. T3 = 19^(L1+L2) + r3 P2 T_pre, T_pre = 19^L1 + r2 P1, K(W) = T3 / S with
S = 2^e2 5^e5 19^g the {2,5,19}-part of T3.
  Type A (not (x3 = N, u2 empty)): e2 = e5 = 0;  D1 := r2 r3 P1 P2;  T3 = D1 (1 + eps), 0 < eps <= 0.11 * 10^-H1.
  Type B (x3 = N, u2 empty, so H2 = 0): T3 = 20 19^L1 + 19 r2 P1, delta = 2^e2 5^e5 <= 40;  D1 := 19 r2 P1;
          T3 = D1 (1 + eps), 0 < eps <= 0.106 * 10^-H1.
  F2 (proved by hand, README 22.2): g = v19(T3) <= 2.
Suppose W != W' share K; by symmetry H1 <= H1'. X := the largest of the four gaps.
  (1,1)  Lambda11 = log(S' D1 / (S D1')) = log((1 + eps') / (1 + eps)), 0 < |Lambda11| < 0.11 * 10^-H1, a form in log 19, 2, 5
         with coefficients <= 6X + 6.  So H1 <= m11.
  (2,1)  If W is Type A: T3 = r3 P2 T_pre (1 + eps2), 0 < eps2 <= 0.1 * 10^-H2, and
         Lambda21 = log(S' r3 P2 T_pre / (S D1')) has 0 < |Lambda21| < 0.11 * 10^-min(H2, H1'), a form in log 19, 2, 5, T_pre,
         one form per prefix with H1 <= m11.  So H2 <= m21 or H1' <= m21.
  (2,2)  If both are Type A: 19^D = 1 modulo 2^min(v2(r3 P2), v2(r3' P2')) and modulo 5^min(v5(r3 P2), v5(r3' P2')),
         D = (L1+L2-g) - (L1'+L2'-g'), D != 0, |D| <= 2X + 2.  So min(H2, H2') <= m2.
  Cases: one of the two words is "short": its gaps add up to at most s = max(m11 + m21, m11 + m2, m21 + m2).
  Size:  100 * 20^(H1'+H2') <= T3' = S' K <= 14440 * T3 <= 14440 * 401 * 200^s, which caps the other word.
The first bound comes from Matveev's theorem (Bugeaud-Mignotte-Siksek 2006, Thm 9.4) with n = 3 and 4; each later one from
exact LLL (de Weger's lemma). Below the final bound every pair is compared exactly.
Usage: python3 pairs18.py   Writes pairs18.json (about two minutes; peak memory well under 1 GB).
"""
import json, math, time
from decimal import Decimal
from baker16 import matveev_constant, largest_ok
from pairs17 import (strip, v, runs, P_of, T2, T3, name2, name3, LN, lattice_bound)

S_MAX = 2 ** 3 * 5 * 19 ** 2          # largest {2,5,19}-part of T3 (delta <= 40, g <= 2)
T3_UP = 401                           # T3 <= 401 * 200^(H1+H2)   (1.11 * 361 = 400.71)
T3_LO = 100                           # T3 >= 100 * 20^(H1+H2)

def prefixes(Hmax):
    return [p for p in runs(Hmax) if p != (19, 0, 0)]

def is_typeB(s):
    return s == (19, 0, 0)

def flog(n, b):
    """floor(log_b n) for integers n >= 1, b >= 2, exactly"""
    k = 0
    while b ** (k + 1) <= n: k += 1
    return k

def m2_of(X):
    """the D-step: min(H2, H2') <= min((2 + floor(log2(2X+2))) / 2, 1 + floor(log5(2X+2)))"""
    n = 2 * X + 2
    return min((2 + flog(n, 2)) // 2, 1 + flog(n, 5))

def size_of(s):
    """largest Y with 100 * 20^Y <= 14440 * 401 * 200^s (exact integers)"""
    rhs = S_MAX * T3_UP * 200 ** s
    return largest_ok(lambda Y: T3_LO * 20 ** Y <= rhs, 0, 4 * s + 20)

def a_priori():
    """Matveev (BMS 2006, Thm 9.4) for prod alpha^b - 1, plus log 2 for the logarithmic form (see 21.5). ell = 1 + log(6X + 6).
    m11 <= a1 ell;  log T_pre <= 3.0 + 5.3 H1 <= aT ell;  m21 <= a21 ell^2;  m2 <= ell;  new X <= aS ell^2."""
    C3 = matveev_constant(n=3, alphas=(19, 2, 5)); C4 = matveev_constant(n=4, alphas=(19, 2, 5))
    l10 = math.log(10); e0 = 1 + math.log(12)
    a1 = C3 / l10 + 1
    aT = 3.0 / e0 + 5.3 * a1
    a21 = C4 * aT / l10 + 1
    aS = (math.log(S_MAX * T3_UP / T3_LO) / e0 ** 2 + math.log(200) * (a1 / e0 + a21)) / math.log(20)
    g = lambda X: X - aS * (1 + math.log(6 * X + 6)) ** 2
    lo, hi = 1.0, 1e300
    for _ in range(4000):
        mid = math.sqrt(lo * hi)
        if g(mid) > 0: hi = mid
        else: lo = mid
    Xs = math.ceil(hi); ell = 1 + math.log(6 * Xs + 6)
    assert g(Xs) > 0 and 2 * aS * ell / (Xs + 1) < 1          # g is increasing from Xs on
    return dict(C3=C3, C4=C4, a1=a1, aT=aT, a21=a21, aS=aS, X_star=Xs,
                check="g(X) = X - aS (1 + log(6X+6))^2 > 0 and increasing for X >= X_star")

def one_round(X, out_rounds):
    t0 = time.time(); X0 = 6 * X + 6
    m11, C3d, _ = lattice_bound(X0, [LN[19], LN[2], LN[5]], "0.11")
    if m11 is None: m11 = X
    m11 = min(m11, X)
    m21 = 0; worst = None; pre = prefixes(m11)
    for p in pre:
        h, _, _ = lattice_bound(X0, [LN[19], LN[2], LN[5], Decimal(T2(p)).ln()], "0.11")
        if h is None: h = X
        if h > m21: m21 = h; worst = p
    m21 = min(m21, X); m2 = min(m2_of(X), X)
    s = max(m11 + m21, m11 + m2, m21 + m2); hS = size_of(s)
    Xn = max(m11, m21, m2, hS)
    rec = dict(X=X, X0=X0, C_digits_3d=len(str(C3d)) - 1, m11=m11, prefixes_checked=len(pre), m21=m21,
               worst_prefix=name2(worst) if worst else None, m2=m2, short_word_sum_s=s, other_word_sum=hS, new_X=Xn,
               seconds=round(time.time() - t0, 1))
    out_rounds.append(rec); print(rec, flush=True)
    return Xn

def all_words(B):
    """canonical three-letter words with H1 + H2 <= B, as (prefix, suffix)"""
    for H1 in range(B + 1):
        for p in runs(H1):
            if p[1] + p[2] != H1 or p == (19, 0, 0): continue
            for s in runs(B - H1): yield p, s

def brute(B):
    seen = {}; clashes = []; n = 0; nB = 0
    for p, s in all_words(B):
        n += 1; nB += is_typeB(s); K = strip(T3(p, s))
        if K in seen: clashes.append((name3(p, s), name3(*seen[K])))
        else: seen[K] = (p, s)
    # planted: non-reduced words with m directly before K are reduced words in disguise (mK = lN) and must be found
    planted = 0; found = 0; wrong = []
    for p, s in all_words(B):
        for which in (0, 1):
            r, L, M = (p, s)[which]
            if r != 19 or L == 0: continue
            nr = (10, L - 1, M + 1)                              # l^(L-1) m^(M+1) K  is  l^L m^M N
            q = (nr, s) if which == 0 else (p, nr)
            planted += 1; K = strip(T3(*q))
            if K in seen:
                found += 1
                if seen[K] != (p, s): wrong.append((name3(*q), name3(*seen[K])))
    return dict(H1_plus_H2_up_to=B, three_letter_words=n, type_B_words=nB, distinct_fingerprints=len(seen), clashes=clashes,
                planted_mK_words=planted, planted_found=found, planted_matched_wrong_word=wrong)

if __name__ == "__main__":
    t0 = time.time(); out = {}
    out["a_priori"] = ap = a_priori(); print("a priori:", ap, flush=True)
    rounds = []; X = ap["X_star"]
    while True:
        Xn = one_round(X, rounds)
        if Xn >= X: break
        X = Xn
    out["rounds"] = rounds; out["final_bound"] = X
    out["exact_check"] = brute(X); print({k: v_ for k, v_ in out["exact_check"].items()}, flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open("pairs18.json", "w"), indent=1, default=str)
    print("done", out["seconds"], "s")
