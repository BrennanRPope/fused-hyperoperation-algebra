"""Round 17 (DRAFT): a fast exact LLL and de Weger's reduction lemma, shared by pairs17.py and checks17.py.

lll_int is the integral LLL algorithm (H. Cohen, A Course in Computational Algebraic Number Theory, Algorithm 2.6.7) with
delta = 3/4: integers only, so it is exact. Its output is re-checked for LLL-reducedness by baker16.is_lll_reduced (exact
rational Gram-Schmidt) before any bound is drawn from it, so a bug here could only cost time, never correctness.
"""
from decimal import Decimal
from fractions import Fraction as F
from baker16 import is_lll_reduced

def lll_int(B):
    b = [list(map(int, r)) for r in B]; n = len(b)
    dot = lambda x, y: sum(p * q for p, q in zip(x, y))
    d = [0] * (n + 1); lam = [[0] * n for _ in range(n)]
    d[0] = 1; d[1] = dot(b[0], b[0]); k = 2; kmax = 1          # 1-based k as in Cohen
    def red(k, l):                                              # 1-based
        if 2 * abs(lam[k - 1][l - 1]) > d[l]:
            q = (2 * lam[k - 1][l - 1] + d[l]) // (2 * d[l])   # nearest integer
            b[k - 1] = [x - q * y for x, y in zip(b[k - 1], b[l - 1])]
            lam[k - 1][l - 1] -= q * d[l]
            for i in range(1, l): lam[k - 1][i - 1] -= q * lam[l - 1][i - 1]
    def swap(k):
        nonlocal kmax
        b[k - 1], b[k - 2] = b[k - 2], b[k - 1]
        for j in range(1, k - 1): lam[k - 1][j - 1], lam[k - 2][j - 1] = lam[k - 2][j - 1], lam[k - 1][j - 1]
        l = lam[k - 1][k - 2]; Bv = (d[k - 2] * d[k] + l * l) // d[k - 1]
        for i in range(k + 1, kmax + 1):
            t = lam[i - 1][k - 1]
            lam[i - 1][k - 1] = (d[k] * lam[i - 1][k - 2] - l * t) // d[k - 1]
            lam[i - 1][k - 2] = (Bv * t + l * lam[i - 1][k - 1]) // d[k]
        d[k - 1] = Bv
    while k <= n:
        if k > kmax:
            kmax = k
            for j in range(1, k + 1):
                u = dot(b[k - 1], b[j - 1])
                for i in range(1, j): u = (d[i] * u - lam[k - 1][i - 1] * lam[j - 1][i - 1]) // d[i - 1]
                if j < k: lam[k - 1][j - 1] = u
                else:
                    d[k] = u
                    if u == 0: raise ValueError("dependent rows")
        while True:
            red(k, k - 1)
            if 4 * d[k] * d[k - 2] < 3 * d[k - 1] ** 2 - 4 * lam[k - 1][k - 2] ** 2:
                swap(k); k = max(2, k - 1)
            else:
                for l in range(k - 2, 0, -1): red(k, l)
                k += 1; break
    return b

def reduce_fast(X0, C, logs):
    """De Weger's lemma, exactly as baker16.reduce_once but with lll_int. logs must be accurate to well below 1/C.
    Returns (lower bound for |Lambda| over integer x != 0 with all |x_i| <= X0, or None; |b1|^2)."""
    n = len(logs)
    phi = [int((Decimal(C) * lg).to_integral_value()) for lg in logs]
    rows = [[1 if j == i else 0 for j in range(n - 1)] + [phi[i]] for i in range(n - 1)] + [[0] * (n - 1) + [phi[-1]]]
    R = lll_int(rows)
    assert is_lll_reduced(R), "LLL output not reduced"
    b1sq = sum(x * x for x in R[0])
    l_sq = F(b1sq, 2 ** (n - 1)); S = (n - 1) * X0 * X0; T = F(n * X0 + 2, 2)
    if l_sq <= T * T + S: return None, b1sq
    lam = (Decimal(l_sq.numerator) / Decimal(l_sq.denominator) - Decimal(S)).sqrt() - Decimal(T.numerator) / Decimal(T.denominator)
    return lam / Decimal(C), b1sq
