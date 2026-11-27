"""Referee-20 engine (written from the definitions in BRIEF.md only).

Words: N x_1 ... x_n, item x = (L, M, y), y in 'NK', K only with M = 0, x_1 != (0,0,'N').
T_0 = 1, T_k = 19^(A_k) + r_k P_k T_(k-1).  kappa = T_n with 2,5,19 removed.

Tools:
  * Matveev (BMS 2006 Thm 9.4) constants,
  * de Weger reduction with python-flint LLL, certified by: output rows lie in the lattice, |det| equal
    (so it is a basis), and lambda_1^2 >= min_k d_k/d_(k-1) (Gram determinants, exact integers), which
    holds for every basis (no reducedness assumption needed),
  * exact 2-adic bound for v2(19^m Q - 1), |m| <= Mb, via Z_2^* = {+-1} x 5^(Z_2),
  * the 2-adic vanishing resolver (see REPORT.md, Lemma N').
"""
import math, functools
from fractions import Fraction as F
import mpmath as mp
import flint

LOG10 = math.log(10)

# ---------------------------------------------------------------- words
def rP(it):
    L, M, y = it
    return (19 if y == 'N' else 10) * 200 ** L * 20 ** M

def gap(it):
    return it[0] + it[1]

def g2(it):
    L, M, y = it
    return 3 * L + 2 * M + (1 if y == 'K' else 0)

def Tval(items):
    T, A = 1, 0
    for it in items:
        A += it[0]
        T = 19 ** A + rP(it) * T
    return T

def strip(n):
    for p in (2, 5, 19):
        while n % p == 0:
            n //= p
    return n

def vp(n, p):
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def items_gap(H):
    out = [(L, H - L, 'N') for L in range(H + 1)]
    out.append((H, 0, 'K'))
    return out

def items_upto(Hmax, first=False):
    out = []
    for H in range(Hmax + 1):
        for it in items_gap(H):
            if first and it == (0, 0, 'N'):
                continue
            out.append(it)
    return out

@functools.lru_cache(maxsize=None)
def items_g2(b):
    """all items with g2 <= b (g2 >= 2*gap, so gap <= b/2)"""
    if b < 0:
        return ()
    return tuple(it for it in items_upto(b // 2) if g2(it) <= b)

# S-part data (from item4.py / v19dp.py, proved there): allowed e5 per number of items
E5 = {0: (0,), 1: (0,), 2: (0, 1), 3: (0, 1)}
SMAX = {0: 1, 1: 19, 2: 40 * 19 ** 2, 3: 80 * 19 ** 3}   # max S = 2^e2 5^e5 19^g for n items

# ---------------------------------------------------------------- Matveev
def matveev_C(n, A):
    c = 1.4 * 30 ** (n + 3) * n ** 4.5
    for a in A:
        c *= a
    return c

LOGS3 = [math.log(2), math.log(5), math.log(19)]

def matveev_h(nlogs, Aalpha, B):
    """bound on h from |Lambda_log| < (1/9) 10^-h, Matveev with logs 2,5,19 (+ alpha)."""
    A = list(LOGS3) + ([max(Aalpha, 0.16)] if nlogs == 4 else [])
    c = matveev_C(nlogs, A) * (1 + math.log(B))
    # |Lambda_mult| <= 1.12 |Lambda_log| < 0.1245 * 10^-h ;  log|Lambda_mult| > -c
    return math.floor((c + math.log(0.1245)) / LOG10)

def logTmax(f, Y):
    """log of an upper bound for T_f (f items, gap sum <= Y): T <= (10/9) 19^f 200^Y"""
    return math.log(10 / 9) + f * math.log(19) + Y * math.log(200)

def derived_bound(p, Hp, q):
    """if word P (p items) has gap sum <= Hp and kappa(P) = kappa(Q), Q with q items:
    T(Q) = T(P) S_Q/S_P <= T(P) SMAX[q], T(Q) >= 10^q 20^{H_Q}."""
    x = (logTmax(p, Hp) + math.log(SMAX[q]) - q * math.log(10)) / math.log(20)
    return max(0, math.floor(x + 1e-9))

# ---------------------------------------------------------------- de Weger with flint
def _certify_and_c1sq(red, a, n):
    """red: list of int rows. check basis of lattice {(x_1..x_{n-1}, sum x_i a_i)}; return lower bound
    (Fraction) for lambda_1^2."""
    an = a[-1]
    for v in red:
        s = v[n - 1] - sum(v[j] * a[j] for j in range(n - 1))
        assert s % an == 0, 'row not in lattice'
    Mred = flint.fmpz_mat(red)
    assert abs(int(Mred.det())) == abs(an), 'not a basis'
    G = Mred * Mred.transpose()
    d = [1]
    for k in range(1, n + 1):
        sub = flint.fmpz_mat([[int(G[i, j]) for j in range(k)] for i in range(k)])
        d.append(int(sub.det()))
        assert d[-1] > 0
    return min(F(d[k], d[k - 1]) for k in range(1, n + 1))

def deweger_once(alphas, X0, C):
    """alphas: mp values of the logs; the LAST one may have coefficient 1 (it is the one attached to a_n).
    All coefficients bounded by X0 in absolute value; returns lower bound on log|Lambda| or None."""
    n = len(alphas)
    a = [int(mp.nint(C * al)) for al in alphas]
    assert a[-1] != 0
    basis = []
    for j in range(n - 1):
        v = [0] * n
        v[j] = 1
        v[n - 1] = a[j]
        basis.append(v)
    v = [0] * n
    v[n - 1] = a[-1]
    basis.append(v)
    red = [[int(x) for x in row] for row in flint.fmpz_mat(basis).lll().tolist()]
    c1sq = _certify_and_c1sq(red, a, n)
    S = (n - 1) * X0 * X0
    if c1sq <= S:
        return None, c1sq
    T = mp.mpf(n) * X0 * (mp.mpf(1) / 2 + mp.mpf(10) ** (-20))
    sq = mp.sqrt(mp.mpf(c1sq.numerator) / c1sq.denominator - S)
    if sq <= T:
        return None, c1sq
    return mp.log((sq - T) / C), c1sq

def lll_h(alpha_fn, n, X0, logk=math.log(1 / 9), extra=6):
    """bound h with |Lambda| < exp(logk) 10^-h, Lambda = sum x_i alpha_i, |x_i| <= X0, x != 0."""
    e = int(math.ceil(n * math.log10(max(X0, 2)))) + extra
    for _ in range(60):
        C = mp.mpf(10) ** e
        with mp.workdps(e + 40):
            low, c1 = deweger_once(alpha_fn(), X0, C)
            if low is not None:
                h = (mp.mpf(logk) - low) / mp.log(10)
                return int(mp.floor(h + mp.mpf(10) ** -15)), e
        e += 2
    raise RuntimeError('reduction failed')

def lll3(X0):
    return lll_h(lambda: [mp.log(2), mp.log(5), mp.log(19)], 3, X0)

def lll4(K, X0):
    """form log K + x2 log2 + x5 log5 + x19 log19, K > 1 not an S-unit."""
    return lll_h(lambda: [mp.log(2), mp.log(5), mp.log(19), mp.log(K)], 4, X0)

# ---------------------------------------------------------------- 2-adic
J2 = 512
MOD2 = 1 << (J2 + 2)
P5 = [pow(5, 1 << j, MOD2) for j in range(J2)]

@functools.lru_cache(maxsize=None)
def st2(x):
    """odd integer x (may be negative) -> (s, t mod 2^J2) with x = (-1)^s 5^t in Z_2"""
    assert x % 2 == 1
    s = 0 if x % 4 == 1 else 1
    y = (x if s == 0 else -x) % MOD2
    t, cur = 0, 1
    for j in range(J2):
        m = 1 << (j + 3)
        if cur % m != y % m:
            t += 1 << j
            cur = cur * P5[j] % MOD2
        assert cur % m == y % m
    return s, t

S19, T19 = st2(19)
assert S19 == 1 and T19 % 2 == 1
INV_T19 = pow(T19, -1, 1 << J2)

def odd_part(x):
    while x % 2 == 0:
        x //= 2
    return x

def vmax2(sQ, tQ, Mb, exclude=None):
    """max over integers m, |m| <= Mb, m != exclude, of v2(19^m Q - 1), Q = (-1)^sQ 5^tQ.
    v2(19^m Q - 1) = 1 if m + sQ odd, else 2 + v2(m - z) with z = -tQ / t19."""
    MODJ = 1 << J2
    z = (-tQ * INV_T19) % MODJ
    if (z - sQ) % 2 != 0:
        return 2            # every m of the right parity has v2(m - z) = 0
    j0 = 1
    while (1 << j0) <= 2 * Mb:
        j0 += 1
    assert j0 + 40 < J2, 'precision'
    mod = 1 << j0
    ms = z % mod
    if ms > mod // 2:
        ms -= mod
    best = j0 - 1
    if abs(ms) <= Mb:
        # other m in range: v2(m - z) = v2(m - ms) <= floor(log2(max |m - ms|))
        best = max(1, int(math.floor(math.log2(Mb + abs(ms)))))
        if ms != exclude:
            d = (ms - z) % MODJ
            if d == 0:
                return None      # vanishing to full precision: caller must handle
            best = max(best, vp(d, 2))
    return 2 + best

def is_19power(q):
    """q Fraction > 0: return k with q = 19^k, else None"""
    num, den = q.numerator, q.denominator
    k = 0
    while num % 19 == 0:
        num //= 19; k += 1
    while den % 19 == 0:
        den //= 19; k -= 1
    if num == 1 and den == 1:
        return k
    return None

class Unresolved(Exception):
    pass

def resolve_vanishing(b, bq, tP, tQ, tauP, tauQ, RP, RQ):
    """vanishing of the main terms at the 2-adic split; returns 'impossible' or 'equal' (W = W'),
    raises Unresolved otherwise.  See REPORT.md Lemma N'."""
    if b == 0 or bq == 0:
        if b == 0 and bq == 0:
            raise Unresolved('both complete')
        return 'impossible'      # remainder of one side is 0, of the other > 0
    if b != bq:
        return 'impossible'      # K(prefix of b letters) = K(prefix of b' letters), b != b' <= 3
    rho = F(RQ * tauP, RP * tauQ)
    v = {}
    for p in (2, 5, 19):
        v[p] = vp(rho.numerator, p) - vp(rho.denominator, p)
    if strip(rho.numerator) != 1 or strip(rho.denominator) != 1:
        return 'impossible'
    if v[2] - 2 * v[5] != v[19]:
        return 'impossible'
    if rho == 1 and tuple(tP) == tuple(tQ):
        return 'equal'
    raise Unresolved(('rho', rho, tP, tQ))

def adic_bounds(p, q, tP, tQ, Mb, log=None):
    """2-adic step.  P has p items with known tail tP (last len(tP) items), Q likewise.
    Returns (bP, bQ): either g2 of P's next unknown item from the back <= bP, or the same for Q <= bQ.
    (None = that side has no unknown item.)"""
    b, bq = p - len(tP), q - len(tQ)
    tauP, tauQ = Tval(tP), Tval(tQ)
    RP = 1
    for it in tP: RP *= rP(it)
    RQ = 1
    for it in tQ: RQ *= rP(it)
    t2P, t2Q = vp(tauP, 2), vp(tauQ, 2)
    oP, oQ = odd_part(tauP), odd_part(tauQ)
    sP, tPv = st2(oP)
    sQ_, tQv = st2(oQ)
    V = 0
    for e5P in E5[p]:
        for e5Q in E5[q]:
            # Q = oP 5^e5Q / (oQ 5^e5P)
            sQ = (sP + sQ_) % 2
            tQ_ = (tPv - tQv + (e5Q - e5P)) % (1 << J2)
            q19 = F(oP * 5 ** e5Q, oQ * 5 ** e5P)
            k = is_19power(q19)
            excl = None
            if k is not None:
                m0 = -k
                res = resolve_vanishing(b, bq, tP, tQ, tauP, tauQ, RP, RQ)
                excl = m0
                if log is not None:
                    log.append(('vanish', tuple(tP), tuple(tQ), e5P, e5Q, m0, res))
            vm = vmax2(sQ, tQ_, Mb, exclude=excl)
            if vm is None:
                raise Unresolved(('precision/vanishing', tP, tQ))
            V = max(V, vm)
    vR_P, vR_Q = vp(RP, 2), vp(RQ, 2)
    bP = (V + t2P - vR_P) if b > 0 else None
    bQ = (V + t2Q - vR_Q) if bq > 0 else None
    return bP, bQ
