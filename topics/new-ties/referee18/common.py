"""Shared tools for the independent referee check.

Everything here is written from the definitions in the brief:
  letters l=(19/200,0), m=(1/20,0), N=(1/19,1/19), K=(1/10,1/10)
  psi(w x) = psi(w)/s_x + d_x,  kappa = psi with 2,5,19 stripped.
"""
from fractions import Fraction as F
import math
import mpmath as mp

S_PRIMES = (2, 5, 19)


def vp(n, p):
    """p-adic valuation of a nonzero integer."""
    n = abs(n)
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def strip(n):
    """remove all factors 2,5,19 from a positive integer"""
    for p in S_PRIMES:
        while n % p == 0:
            n //= p
    return n


def kappa_of_fraction(q):
    q = F(q)
    return F(strip(q.numerator), strip(q.denominator))


def P(L, M):
    return 200 ** L * 20 ** M


def rval(x):
    return 19 if x == 'N' else 10


# ---------- two-letter words  V = N l^L m^M x ----------
def two_letter_words(Hmax):
    """all reduced two-letter words with H=L+M<=Hmax, excluding NN.
    x=K only if M=0."""
    out = []
    for H in range(Hmax + 1):
        for L in range(H + 1):
            M = H - L
            for x in 'NK':
                if x == 'K' and M > 0:
                    continue
                if x == 'N' and L == 0 and M == 0:
                    continue  # NN
                out.append((L, M, x))
    return out


def T2(L, M, x):
    return 19 ** L + rval(x) * P(L, M)


def kappa2(L, M, x):
    """fingerprint (integer) and g = v19(T)"""
    T = T2(L, M, x)
    k = strip(T)
    return k, vp(T, 19), T


# ---------- three-letter words  W = N u1 x2 u2 x3 ----------
def three_letter_words(H1max, H2max, H1min=0, H2min=0):
    """reduced three-letter words not starting NN."""
    for H1 in range(H1min, H1max + 1):
        for L1 in range(H1 + 1):
            M1 = H1 - L1
            for x2 in 'NK':
                if x2 == 'K' and M1 > 0:
                    continue
                if x2 == 'N' and H1 == 0:
                    continue  # starts NN
                for H2 in range(H2min, H2max + 1):
                    for L2 in range(H2 + 1):
                        M2 = H2 - L2
                        for x3 in 'NK':
                            if x3 == 'K' and M2 > 0:
                                continue
                            yield (L1, M1, x2, L2, M2, x3)


def T3(L1, M1, x2, L2, M2, x3):
    r2, r3 = rval(x2), rval(x3)
    Tpre = 19 ** L1 + r2 * P(L1, M1)
    return 19 ** (L1 + L2) + r3 * P(L2, M2) * Tpre


def is_typeB(w):
    L1, M1, x2, L2, M2, x3 = w
    return L2 == 0 and M2 == 0 and x3 == 'N'


# ---------- direct affine-map evaluation (independent of the T formulas) ----------
MAPS = {'l': (F(19, 200), F(0)), 'm': (F(1, 20), F(0)),
        'N': (F(1, 19), F(1, 19)), 'K': (F(1, 10), F(1, 10))}


def compose_word(word):
    """word is a string; w x = w o x (apply x first).  returns (s,t)."""
    s, t = F(1), F(0)
    for ch in word:
        sx, tx = MAPS[ch]
        # (s,t) o (sx,tx): y -> s(sx y + tx) + t
        s, t = s * sx, s * tx + t
    return s, t


def psi_word(word):
    s, t = compose_word(word)
    return t / s


def word_string2(L, M, x):
    return 'N' + 'l' * L + 'm' * M + x


def word_string3(w):
    L1, M1, x2, L2, M2, x3 = w
    return 'N' + 'l' * L1 + 'm' * M1 + x2 + 'l' * L2 + 'm' * M2 + x3


# ---------- Matveev (BMS 2006 Thm 9.4) ----------
def matveev_C(n, A, D=1):
    """constant c with log|Lambda| > -c (1+log B)."""
    c = 1.4 * 30 ** (n + 3) * n ** 4.5 * D ** 2 * (1 + math.log(D))
    for a in A:
        c *= a
    return c


LOG2, LOG5, LOG19, LOG10 = (math.log(2), math.log(5), math.log(19), math.log(10))


# ---------- 2-adic tools ----------
def v2_19pow_minus1(D):
    """v2(19^D - 1) for an integer D != 0"""
    assert D != 0
    if D % 2:
        return 1
    return 2 + vp(D, 2)


def dlog19_2adic(K, J):
    """z mod 2^J with 19^z = K in Z_2, or None if K not in closure of <19>
    (i.e. K mod 8 not in {1,3})."""
    if K % 8 not in (1, 3):
        return None
    z = 0
    for j in range(J):
        mod = 1 << (j + 3)
        if pow(19, z, mod) != K % mod:
            z += 1 << j
        assert pow(19, z, mod) == K % mod
    return z


def max_v2_K19m(K, Mb, J=None, exclude=None):
    """an upper bound for  max_{|m|<=Mb} v2(K*19^(-m) - 1)  (exact max up to the
    j0-1 slack described in the notes).  Returns integer."""
    if K % 8 not in (1, 3):
        return 2
    Mb = int(Mb)
    j0 = 1
    while (1 << (j0 - 1)) <= Mb:
        j0 += 1
    # residues of [-Mb,Mb] mod 2^j0 are distinct since 2^j0 > 2Mb
    if J is None:
        J = j0 + 64
    while True:
        z = dlog19_2adic(K, J)
        mod = 1 << j0
        s = z % mod
        if s >= mod // 2:
            s -= mod
        if abs(s) > Mb or s == exclude:
            # (exclude: the one m allowed to have v2(z-m) >= j0 is ruled out)
            vmax = j0 - 1
            break
        # v2(z - s) exactly, need it < J
        diff = (z - s) % (1 << J)
        if diff == 0:
            J *= 2
            continue
        vmax = max(j0 - 1, vp(diff, 2))
        break
    # v2(19^y - 1) = 1 if y odd, 2+v2(y) if y even
    return 2 + vmax


# ---------- integral LLL (Cohen, Algorithm 2.6.7), delta = p/q ----------
def lll_int(basis, dnum=99, dden=100):
    b = [list(map(int, v)) for v in basis]
    n = len(b)

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    d = [0] * (n + 1)  # d[0]=1, d[i] for i=1..n
    lam = [[0] * n for _ in range(n)]
    d[0] = 1
    d[1] = dot(b[0], b[0])
    k = 1  # 0-based index of current vector (Cohen's k=2)
    kmax = 0

    def red(k, l):
        dl = d[l + 1]
        if 2 * abs(lam[k][l]) > dl:
            q = (2 * lam[k][l] + dl) // (2 * dl)
            bk, bl = b[k], b[l]
            for i in range(len(bk)):
                bk[i] -= q * bl[i]
            lam[k][l] -= q * dl
            for i in range(l):
                lam[k][i] -= q * lam[l][i]

    def swap(k):
        b[k], b[k - 1] = b[k - 1], b[k]
        for j in range(k - 1):
            lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
        lm = lam[k][k - 1]
        B = (d[k - 1] * d[k + 1] + lm * lm) // d[k]
        for i in range(k + 1, kmax + 1):
            t = lam[i][k]
            lam[i][k] = (d[k + 1] * lam[i][k - 1] - lm * t) // d[k]
            lam[i][k - 1] = (B * t + lm * lam[i][k]) // d[k + 1]
        d[k] = B

    while k < n:
        if k > kmax:
            kmax = k
            for j in range(k + 1):
                u = dot(b[k], b[j])
                for i in range(j):
                    u = (d[i + 1] * u - lam[k][i] * lam[j][i]) // d[i]
                if j < k:
                    lam[k][j] = u
                else:
                    d[k + 1] = u
                    if u == 0:
                        raise ValueError('dependent vectors')
        red(k, k - 1)
        if dden * d[k + 1] * d[k - 1] < dnum * d[k] ** 2 - dden * lam[k][k - 1] ** 2:
            swap(k)
            k = max(1, k - 1)
        else:
            for l in range(k - 2, -1, -1):
                red(k, l)
            k += 1
    return b


def gram_schmidt_exact(b):
    n = len(b)
    bs = []
    mu = [[F(0)] * n for _ in range(n)]
    nrm = []
    for i in range(n):
        v = [F(x) for x in b[i]]
        for j in range(i):
            mu[i][j] = sum(F(x) * y for x, y in zip(b[i], bs[j])) / nrm[j]
            v = [vi - mu[i][j] * wj for vi, wj in zip(v, bs[j])]
        bs.append(v)
        nrm.append(sum(x * x for x in v))
    return mu, nrm


def check_lll_reduced(b, delta=F(99, 100)):
    """exact check: size reduced and Lovasz with delta. returns (ok, nrm)"""
    mu, nrm = gram_schmidt_exact(b)
    n = len(b)
    for i in range(n):
        for j in range(i):
            if abs(mu[i][j]) > F(1, 2):
                return False, nrm
    for k in range(1, n):
        if nrm[k] < (delta - mu[k][k - 1] ** 2) * nrm[k - 1]:
            return False, nrm
    return True, nrm


def deweger(alphas_mp, X0, C, delta=F(99, 100), check=True):
    """de Weger lemma.  alphas_mp: list of mpf (precision must exceed log10 C).
    Lattice: rows b_j = e_j + a_j e_n (j<n), b_n = a_n e_n, a_i = nint(C alpha_i).
    Returns (lower bound on log|Lambda| as mpf or None, info dict)."""
    n = len(alphas_mp)
    a = [int(mp.nint(C * al)) for al in alphas_mp]
    assert a[-1] != 0
    basis = []
    for j in range(n - 1):
        v = [0] * n
        v[j] = 1
        v[n - 1] = a[j]
        basis.append(v)
    v = [0] * n
    v[n - 1] = a[n - 1]
    basis.append(v)
    red = lll_int(basis)
    ok, nrm = check_lll_reduced(red, delta)
    if not ok:
        raise RuntimeError('LLL output not reduced')
    # also check same lattice: |det| preserved
    c2_a = min(nrm)
    b1sq = sum(x * x for x in red[0])
    c2_b = F(b1sq) * (delta - F(1, 4)) ** (n - 1)
    c2 = max(c2_a, c2_b)
    S = (n - 1) * X0 ** 2
    T = mp.mpf(n) * X0 * (mp.mpf(1) / 2 + mp.mpf(10) ** (-30))
    info = dict(n=n, C=C, X0=X0, c2=c2, log10c=float(mp.log10(mp.sqrt(mp.mpf(c2.numerator) / c2.denominator))))
    if c2 <= S:
        return None, info
    sq = mp.sqrt(mp.mpf(c2.numerator) / c2.denominator - S)
    if sq <= T:
        return None, info
    low = (sq - T) / C
    return mp.log(low), info


def reduce_bound(alpha_fn, n, X0, logk, logbase, C_start_exp=None, extra=8, maxtries=40):
    """Given |Lambda| <= exp(logk) * base^{-H}, Lambda = sum x_i alpha_i,
    |x_i|<=X0 (x != 0), return (newH, info) with H <= newH.
    alpha_fn() must return the alphas computed at the current mp.dps;
    precision is set to (digits of C)+50 before each call."""
    if C_start_exp is None:
        C_start_exp = int(math.ceil(n * math.log10(max(X0, 2)))) + extra
    e = C_start_exp
    for _ in range(maxtries):
        C = 10 ** e
        with mp.workdps(e + 60):
            alphas = alpha_fn()
            lowlog, info = deweger(alphas, X0, C)
            if lowlog is not None:
                Hb = (mp.mpf(logk) - lowlog) / mp.mpf(logbase)
                info['Cexp'] = e
                return int(mp.floor(Hb + mp.mpf(10) ** -20)), info
        e += 3
    raise RuntimeError('reduction failed')
