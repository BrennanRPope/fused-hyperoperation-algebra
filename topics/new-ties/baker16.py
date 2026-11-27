"""Round 16 (DRAFT): Theorem 16.1, no word with three N/K letters has psi equal to a product of powers of 2, 5 and 19
(fingerprint 1). Computer-assisted: an explicit lower bound for linear forms in logarithms (Matveev 2000, in the form of
Bugeaud, Mignotte and Siksek 2006, Thm 9.4) bounds the exponents, an exact LLL reduction (de Weger's method) cuts the bound
to a few dozen, and exact integer arithmetic checks everything below it.

The word is x1 u1 x2 u2 x3 (x1 = N without loss; leading l, m do not change psi, trailing ones multiply it by a unit).
u_i holds L_i letters l and M_i letters m (ml = lm as maps). With A_i = (200/19)^L_i 20^M_i and r_N = 19, r_K = 10:
    psi = 1 + r3 A2 (1 + r2 A1).
Case I  (x3 = N and u2 empty):  psi = 20 + 19 r2 A1.  The 2- and 5-adic sizes force psi = 2^al 5^be 19^ga with al <= 2,
        be <= 1, so  Lambda = log(psi / (19 r2 A1)) = log(1 + 20/c)  is a positive linear form in log 2, log 5, log 19
        smaller than 10^-H / 9, H = L1 + M1, with coefficients at most 3H + 5.
Case II (otherwise): psi - 1 has positive 2- and 5-adic size, so psi = 19^E. Then a2 <= v2(19^E - 1) <= 2 + log2 E,
        b2 <= v5(19^E - 1) <= 1 + log5 E, L2 = 2 b2 - a2, and the 19-adic size gives L1 <= 2. So two of the four terms of
        19^(L1+L2) psi = 19^(L1+L2) + r3 19^L1 P2 + r2 r3 P1 P2  (P_i = 2^a_i 5^b_i)  are at most polynomial in E, and
        Lambda = log(19^(L1+L2+E) / (r2 r3 P1 P2)) lies in (0, 2(1 + 380 E^2)/19^E), coefficients at most B(E).
Usage: python3 baker16.py   Writes baker16.json.
"""
import json, math, time
from fractions import Fraction as F
from decimal import Decimal, getcontext
getcontext().prec = 160
LOG = {p: Decimal(p).ln() for p in (2, 5, 19)}
ALPHA = (19, 2, 5)

def matveev_constant(n=3, D=1, alphas=ALPHA):
    # log|Lambda| > -1.4 * 30^(n+3) * n^4.5 * D^2 * (1 + log D) * (1 + log B) * A1...An,  A_i >= max(D h, |log a|, 0.16)
    # (for the alphas given; a missing A_i is left to the caller as a factor)
    A = 1.0
    for a in alphas: A *= max(D * math.log(a), abs(math.log(a)), 0.16)
    return 1.4 * 30 ** (n + 3) * n ** 4.5 * D ** 2 * (1 + math.log(D)) * A

def largest_ok(pred, lo, hi):
    """largest integer x in [lo, hi] with pred(x) true, pred monotone (true then false)"""
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if pred(mid): lo = mid
        else: hi = mid - 1
    return lo

def lll(B, delta=F(3, 4)):
    B = [list(map(int, r)) for r in B]; n = len(B); m = len(B[0])
    def gs(B):
        Bs = []; mu = [[F(0)] * n for _ in range(n)]
        for i in range(n):
            v = [F(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(F(B[i][k]) * Bs[j][k] for k in range(m)) / sum(x * x for x in Bs[j])
                v = [v[k] - mu[i][j] * Bs[j][k] for k in range(m)]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gs(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [B[k][i] - q * B[j][i] for i in range(m)]; Bs, mu = gs(B)
        nk = sum(x * x for x in Bs[k]); nk1 = sum(x * x for x in Bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1: k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(B); k = max(k - 1, 1)
    return B

def reduce_once(X0, C, alphas=ALPHA, logs=None):
    """Lower bound for |Lambda| = |x1 log a1 + ... + xn log an| over integer x != 0 with |x_i| <= X0 (de Weger's lemma).
    Lattice rows (e_i, phi_i) for i < n and (0, ..., 0, phi_n), phi_i = nearest integer to C log a_i. If the LLL-reduced
    first vector b1 has |b1|^2 / 2^(n-1) > T^2 + S, with S = (n-1) X0^2 and T = (n X0 + 2)/2, then
    |Lambda| >= (sqrt(|b1|^2 / 2^(n-1) - S) - T) / C."""
    n = len(alphas); logs = logs or [Decimal(a).ln() for a in alphas]
    phi = [int((Decimal(C) * lg).to_integral_value()) for lg in logs]
    rows = [[1 if j == i else 0 for j in range(n - 1)] + [phi[i]] for i in range(n - 1)] + [[0] * (n - 1) + [phi[-1]]]
    R = lll(rows)
    assert is_lll_reduced(R), "LLL output not reduced"
    b1sq = sum(x * x for x in R[0])
    l_sq = F(b1sq, 2 ** (n - 1))
    S = (n - 1) * X0 * X0; T = F(n * X0 + 2, 2)
    if l_sq <= T * T + S: return None, b1sq
    lam = (Decimal(int(l_sq.numerator)) / Decimal(int(l_sq.denominator)) - Decimal(S)).sqrt() - Decimal(T.numerator) / Decimal(T.denominator)
    return lam / Decimal(C), b1sq           # |Lambda| >= this

def is_lll_reduced(B, delta=F(3, 4)):
    n = len(B); Bs = []; mu = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        vv = [F(x) for x in B[i]]
        for j in range(i):
            mu[i][j] = sum(F(B[i][k]) * Bs[j][k] for k in range(len(B[i]))) / sum(x * x for x in Bs[j])
            vv = [vv[k] - mu[i][j] * Bs[j][k] for k in range(len(vv))]
        Bs.append(vv)
    size = all(abs(mu[i][j]) <= F(1, 2) for i in range(n) for j in range(i))
    lov = all(sum(x * x for x in Bs[k]) >= (delta - mu[k][k - 1] ** 2) * sum(x * x for x in Bs[k - 1]) for k in range(1, n))
    return size and lov

def v(n, p):
    k = 0
    while n % p == 0: n //= p; k += 1
    return k
def kappa_frac(q):
    def strip(a):
        for p in (2, 5, 19):
            while a % p == 0: a //= p
        return a
    return strip(q.numerator), strip(q.denominator)
def psi3(r2, L1, M1, r3, L2, M2):
    A1 = F(200, 19) ** L1 * 20 ** M1; A2 = F(200, 19) ** L2 * 20 ** M2
    return 1 + r3 * A2 * (1 + r2 * A1)

if __name__ == "__main__":
    t0 = time.time(); out = {}
    CM = matveev_constant(); out["matveev_constant"] = CM
    # ---- Case II: E log 19 - log(2(1 + 380 E^2)) < CM (1 + log B(E))
    BE = lambda E: math.ceil(4.25 * (E + 4 + 2 * math.log(max(E, 1), 5)))
    E0 = largest_ok(lambda E: E * math.log(19) - math.log(2 * (1 + 380 * E * E)) < CM * (1 + math.log(BE(E))), 4, 10 ** 16)
    # ---- Case I: H log 10 + log 9 < CM (1 + log(3H + 5))
    H0 = largest_ok(lambda H: H * math.log(10) + math.log(9) < CM * (1 + math.log(3 * H + 5)), 1, 10 ** 16)
    out["matveev_bounds"] = dict(case_II_E=E0, case_II_coefficients=BE(E0), case_I_H=H0, case_I_coefficients=3 * H0 + 5)
    print("Matveev:", out["matveev_bounds"], flush=True)
    steps = []
    for case, bound, X0of, ok in (
            ("II", E0, BE, lambda E, lam: Decimal(E) * LOG[19] - Decimal(2 * (1 + 380 * E * E)).ln() < -lam.ln()),
            ("I", H0, lambda H: 3 * H + 5, lambda H, lam: Decimal(H) * Decimal(10).ln() + Decimal(9).ln() < -lam.ln())):
        cur = bound
        for rnd in range(6):
            X0 = X0of(cur); C = 10 ** (3 * len(str(X0)) + 6)
            lam, b1sq = reduce_once(X0, C)
            if lam is None or lam <= 0:
                C *= 10 ** 6; lam, b1sq = reduce_once(X0, C)
            new = largest_ok(lambda x: ok(x, lam), 1, cur)
            steps.append(dict(case=case, round=rnd + 1, old_bound=cur, X0=X0, C_digits=len(str(C)) - 1, b1_length=float(Decimal(b1sq).sqrt()),
                              lower_bound_for_Lambda='%.4e' % lam, new_bound=new))
            print(steps[-1], flush=True)
            if new >= cur: break
            cur = new
        out["case_" + case + "_final_bound"] = cur
    out["lll_rounds"] = steps
    # ---- exact check below the bounds
    E1 = out["case_II_final_bound"]; H1 = out["case_I_final_bound"]
    found = []; checked = 0
    for E in range(1, E1 + 1):
        D = 19 ** E - 1; v2, v5 = v(D, 2), v(D, 5)
        for r3 in (19, 10):
            for L2 in range(0, 2 * v5 + 1):
                for M2 in range(0, v2 + 1):
                    a2, b2 = 3 * L2 + 2 * M2, 2 * L2 + M2
                    if a2 > v2 or b2 > v5: continue
                    if r3 == 19 and L2 == 0 and M2 == 0: continue        # that is Case I
                    for r2 in (19, 10):
                        for L1 in range(0, 3):
                            checked += 1
                            num = 19 ** (L1 + L2) * D
                            den = r3 * 2 ** a2 * 5 ** b2
                            if num % den: continue
                            R = num // den - 19 ** L1                     # must be r2 * 2^(3L1) 5^(2L1) * 20^M1
                            base = r2 * 2 ** (3 * L1) * 5 ** (2 * L1)
                            if R <= 0 or R % base: continue
                            q = R // base; M1 = 0
                            while q % 20 == 0: q //= 20; M1 += 1
                            if q == 1: found.append(("II", E, r2, L1, M1, r3, L2, M2))
    for H in range(0, H1 + 1):
        for L1 in range(0, H + 1):
            M1 = H - L1
            for r2 in (19, 10):
                checked += 1
                if kappa_frac(psi3(r2, L1, M1, 19, 0, 0)) == (1, 1): found.append(("I", r2, L1, M1))
    out["exact_check"] = dict(case_II_E_up_to=E1, case_I_H_up_to=H1, candidates_checked=checked, solutions=found)
    # ---- direct brute force on the formula as an independent check (runs up to 12)
    bf = 0; bf_units = []
    for r2 in (19, 10):
        for r3 in (19, 10):
            for L1 in range(13):
                for M1 in range(13):
                    for L2 in range(13):
                        for M2 in range(13):
                            bf += 1
                            if kappa_frac(psi3(r2, L1, M1, r3, L2, M2)) == (1, 1): bf_units.append((r2, L1, M1, r3, L2, M2))
    out["brute_force_runs_up_to_12"] = dict(words=bf, unit_psi=bf_units)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open("baker16.json", "w"), indent=1, default=str)
    print("solutions:", found, "brute force units:", bf_units, "time", out["seconds"])
