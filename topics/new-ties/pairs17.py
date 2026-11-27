"""Round 17 (DRAFT): Theorem 17.1, no word with three N/K letters shares a fingerprint with a word with two. Computer-assisted.

Words. A reduced two-letter word is V = N l^L' m^M' x' (x' = K only if M' = 0), T' = 19^L' + r' P', K(V) = T' / 19^g'.
A reduced three-letter word is W = N u1 x2 u2 x3, u_i = l^L_i m^M_i, H_i = L_i + M_i, x_i = K only if M_i = 0, and not
(u1 empty and x2 = N), since N N ... = N m ... has only two N/K letters. Its psi is T3 / 19^(L1+L2) with
    T3 = 19^(L1+L2) + r3 P2 T_pre,   T_pre = 19^L1 + r2 P1   (the T of the prefix N u1 x2),
and K(W) = T3 with every factor 2, 5 and 19 removed. Suppose K(W) = K(V) = K. X := max(H1, H2, H').
Type A (not (x3 = N and u2 empty)): T3, T_pre are prime to 10, and
  (L3)  Lambda3 = log(19^(g'-g3) r2 r3 P1 P2 / (r' P')), 0 < |Lambda3| < 0.11 * 10^-min(H1, H'), coefficients <= 6X + 6
        (nonzero: otherwise 19^L2 + r3 P2 would be a product of 2s, 5s and 19s, so u2 x3 = N by Theorem 15.5).
  (D)   19^D = 1 mod 2^min(v2(r3 P2), v2(r' P')) with D = (L1+L2-g3) - (L'-g'), D != 0 (else T_pre would be a product of
        2s, 5s and 19s, so u1 x2 = N, excluded); hence min(H2, H') <= (2 + floor(log2(5X + 3))) / 2.
  (L1)  Lambda1 = log(19^g3 K / (r2 r3 P1 P2)) in (0, 0.11 * 10^-H1), a form in log 19, 2, 5, K, coefficients <= 6X + 6.
  (2ad) K 19^m = 1 mod 2^v2(r3 P2) with m = g3 - L1 - L2, |m| <= 4X + 3, and v2(r3 P2) >= 2 H2.
  (L2)  Lambda2 = log(19^g3 K / (r3 P2 T_pre)) in (0, 0.1 * 10^-H2), a form in log 19, 2, 5, K, T_pre (theory only).
Type B (x3 = N, u2 empty): T3 = 20 19^L1 + 19 r2 P1, 19^g3 * delta * K = T3 with delta = 2^e2 5^e5, and
  (LB)  Lambda_B = log(19^(g'+1-g3) r2 P1 / (delta r' P')), 0 < |Lambda_B| < 0.11 * 10^-min(H1, H') (nonzero by the 2-, 5-
        and 19-adic sizes, see README 21.3), and (LB1) log(delta 19^g3 K / (19 r2 P1)) in (0, 0.106 * 10^-H1).
Size: T' <= 19 K <= 19 T3 <= 19 * 381 * 200^(H1+H2) and T' >= 10 * 20^H'.
Case split each round: either H' <= max(m1, m2) ("V short": then K is one of finitely many and (L1)/(LB1), (2ad) bound
H1, H2), or H1 <= m1 and H2 <= m2 ("W short": then the size bound caps H'). The first bound comes from Matveev's theorem
(Bugeaud-Mignotte-Siksek 2006, Thm 9.4) with n = 3, 4, 5; each later one from exact LLL (de Weger's lemma) and an exact 2-adic
discrete logarithm. Below the final bound every pair is compared exactly.
Usage: python3 pairs17.py   Writes pairs17.json (about a few minutes).
"""
import json, math, time
from decimal import Decimal, getcontext
from baker16 import matveev_constant, largest_ok
from lattice17 import reduce_fast
getcontext().prec = 700
LN = {p: Decimal(p).ln() for p in (2, 5, 19, 10)}
LN2F = math.log(2)

def strip(a):
    for p in (2, 5, 19):
        while a % p == 0: a //= p
    return a
def v(a, p):
    k = 0
    while a % p == 0: a //= p; k += 1
    return k
def runs(Hmax):
    """(r, L, M) for a gap l^L m^M followed by x (r = 19 for N, 10 for K; K only when M = 0), H = L + M <= Hmax"""
    for H in range(Hmax + 1):
        for L in range(H + 1):
            M = H - L
            yield (19, L, M)
            if M == 0: yield (10, L, M)
P_of = lambda L, M: 2 ** (3 * L + 2 * M) * 5 ** (2 * L + M)
def T2(w):
    r, L, M = w
    return 19 ** L + r * P_of(L, M)
def T3(p, s):
    (r2, L1, M1), (r3, L2, M2) = p, s
    return 19 ** (L1 + L2) + r3 * P_of(L2, M2) * (19 ** L1 + r2 * P_of(L1, M1))
def name2(w):
    r, L, M = w
    return "N" + "l" * L + "m" * M + ("N" if r == 19 else "K")
def name3(p, s):
    return name2(p) + name2(s)[1:]
def two_letter_words(Hmax):
    return [w for w in runs(Hmax) if w != (19, 0, 0)]                        # NN excluded (fingerprint 1)
def three_letter_words(H1max, H2max):
    for p in runs(H1max):
        if p == (19, 0, 0): continue                                         # N N ... = N m ...: two letters
        for s in runs(H2max): yield p, s

def jmax(K, M):
    """max over integers |m| <= M of v2(19^m K - 1), exactly (K odd, K != 1, prime to 19). Uses the 2-adic discrete logarithm:
    19 generates the closed subgroup {1, 3 mod 8} of Z_2^*; v2(19^x - 1) = 1 for odd x and 2 + v2(x) for even x in Z_2."""
    if K % 8 in (5, 7): return 2                                             # 19^m K = 5 or 7 mod 8
    J = max(M.bit_length() + 8, 16)
    while True:
        mod = 1 << (J + 2); Kinv = pow(K, -1, mod); lam = 0
        for n in range(3, J + 3):                                            # now 19^lam = K^-1 mod 2^(n-1)
            if pow(19, lam, 1 << n) != Kinv % (1 << n): lam += 1 << (n - 3)
        t = 0                                                                # lam = dlog_19(K^-1) mod 2^J
        while t < J:
            c = lam % (1 << (t + 1)); m = -M + ((c + M) % (1 << (t + 1)))
            if m > M: break
            t += 1
        if t < J: return 1 if t == 0 else 2 + t
        J *= 2

def form_bound(lam, thr):
    """largest h >= 0 with thr * 10^-h > lam (a lower bound lam for |Lambda| and |Lambda| < thr 10^-h)"""
    if lam is None: return None
    h = math.floor(float((Decimal(thr) / lam).log10()))
    while h >= 0 and not Decimal(thr) * Decimal(10) ** (-h) > lam: h -= 1
    while Decimal(thr) * Decimal(10) ** (-(h + 1)) > lam: h += 1
    return max(h, 0)

def lattice_bound(X0, logs, thr, tries=6):
    n = len(logs); C = 10 ** (n * len(str(X0)) + 2 * n)
    for _ in range(tries):
        lam, b1sq = reduce_fast(X0, C, logs)
        if lam is not None and lam > 0: return form_bound(lam, thr), C, lam
        C *= 10 ** 6
    return None, C, None

def hW_of(m1, m2):
    """largest H' with 10 * 20^H' <= 19 * 381 * 200^(m1 + m2) (exact integers)"""
    rhs = 19 * 381 * 200 ** (m1 + m2)
    return largest_ok(lambda H: 10 * 20 ** H <= rhs, 0, 4 * (m1 + m2) + 10)

def m2_of(X):
    return (2 + (5 * X + 3).bit_length() - 1) // 2

def a_priori():
    """Matveev (BMS 2006, Thm 9.4) for Lambda = prod alpha^b - 1; for the logarithmic form |log(1+y)| <= 2|y| when small,
    so log|Lambda_log| > -C(1 + log B) - log 2. Each bound below is <= a * ell^k with ell = 1 + log(6X + 6) >= 3.48."""
    C3 = matveev_constant(n=3, alphas=(19, 2, 5)); C4 = matveev_constant(n=4, alphas=(19, 2, 5))
    C5 = matveev_constant(n=5, alphas=(19, 2, 5)); l10 = math.log(10); e0 = 1 + math.log(12)
    a1 = C3 / l10 + 1                                     # m1 <= a1 ell  (and m2 <= ell <= a1 ell)
    aK = 3.0 / e0 + 5.3 * a1                              # log K <= 3.0 + 5.3 hV <= aK ell
    a_1 = C4 * aK / l10 + 1                               # h1 <= a_1 ell^2  (A4 = log K, B = 6X + 6)
    aP = 3.0 / e0 ** 2 + 5.3 * a_1                        # log T_pre <= 3.0 + 5.3 h1 <= aP ell^2
    a2 = C5 * aK * aP / l10 + 1                           # h2 <= a2 ell^4  (A4 = log K, A5 = log T_pre)
    aW = (math.log(19 * 381 / 10) / e0 + math.log(200) * (a1 + 1)) / math.log(20)   # hW <= aW ell
    assert a2 * e0 ** 3 > max(a1, aW) and a2 * e0 ** 2 > a_1      # so F(X) <= a2 ell^4 for every X >= 1
    g = lambda X: X - a2 * (1 + math.log(6 * X + 6)) ** 4
    lo, hi = 1.0, 1e300
    for _ in range(4000):
        mid = math.sqrt(lo * hi)
        if g(mid) > 0: hi = mid
        else: lo = mid
    Xs = math.ceil(hi)
    ell = 1 + math.log(6 * Xs + 6)
    assert g(Xs) > 0 and 4 * a2 * ell ** 3 / (Xs + 1) < 1   # g increasing from Xs on (ell^3/(X+1) decreases for ell > 3)
    return dict(C3=C3, C4=C4, C5=C5, a1=a1, aK=aK, a_1=a_1, aP=aP, a2=a2, aW=aW, X_star=Xs,
                check="g(X) = X - a2 (1 + log(6X+6))^4 > 0 and increasing for X >= X_star")

def one_round(X, out_rounds, matveev_per_K_over=40000):
    t0 = time.time(); X0 = 6 * X + 6
    m1, C3d, lam3 = lattice_bound(X0, [LN[19], LN[2], LN[5]], "0.11")
    m1 = min(m1, X); m2 = min(m2_of(X), X); hV = max(m1, m2); hW = min(hW_of(m1, m2), X)
    Vs = two_letter_words(hV); Ks = sorted({strip(T2(w)) for w in Vs})
    C4 = matveev_constant(n=4, alphas=(19, 2, 5)); h1max = 0; h2max = 0; worst = None; method = "LLL per K"
    M = 4 * X + 3
    for K in Ks:
        lk = Decimal(K).ln()
        if len(Ks) > matveev_per_K_over:
            method = "Matveev per K"
            h1 = math.floor((C4 * max(float(lk), 0.16) * (1 + math.log(X0)) + LN2F + math.log(0.11)) / math.log(10))
        else:
            h1, _, _ = lattice_bound(X0, [LN[19], LN[2], LN[5], lk], "0.11")
            if h1 is None: h1 = X
        h2 = jmax(K, M) // 2
        if h1 > h1max: h1max = h1; worst = K
        h2max = max(h2max, h2)
    h1max = min(h1max, X); h2max = min(h2max, X)
    Xn = max(hV, h1max, h2max, hW)
    rec = dict(X=X, X0=X0, C_digits_3d=len(str(C3d)) - 1, m1=m1, m2=m2, V_short_hV=hV, V_words=len(Vs), distinct_K=len(Ks),
               per_K_method=method, H1_bound_V_short=h1max, worst_K_digits=len(str(worst)) if worst else 0, H2_bound_V_short=h2max,
               W_short_Hprime_bound=hW, new_X=Xn, seconds=round(time.time() - t0, 1))
    out_rounds.append(rec); print(rec, flush=True)
    return Xn

def brute(Xf):
    Vs = two_letter_words(Xf); setV = {}
    for w in Vs: setV.setdefault(strip(T2(w)), w)
    nW = 0; clashes = []
    for p, s in three_letter_words(Xf, Xf):
        nW += 1; K = strip(T3(p, s))
        if K in setV: clashes.append((name3(p, s), name2(setV[K])))
    planted = 0; planted_found = 0                 # N N u2 x3 is N m u2 x3 in disguise: must be found
    for s in runs(Xf - 1):
        planted += 1
        if strip(T3((19, 0, 0), s)) in setV: planted_found += 1
    return dict(H_up_to=Xf, three_letter_words=nW, two_letter_words=len(Vs), clashes=clashes,
                planted_NN_words=planted, planted_found=planted_found)

if __name__ == "__main__":
    t0 = time.time(); out = {}
    out["a_priori"] = ap = a_priori(); print("a priori:", ap, flush=True)
    rounds = []; X = ap["X_star"]
    while True:
        Xn = one_round(X, rounds)
        if Xn >= X: break
        X = Xn
    out["rounds"] = rounds; out["final_bound"] = X
    out["exact_check"] = brute(X); print(out["exact_check"], flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open("pairs17.json", "w"), indent=1, default=str)
    print("done", out["seconds"], "s")
