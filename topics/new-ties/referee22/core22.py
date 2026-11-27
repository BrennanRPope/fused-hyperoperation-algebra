"""Referee 22 core: letters, words, psi, canonical form, fingerprints, word ratios.
Written from the statements in BRIEF.md only.
A word is a string over 'l','m','N','K'; as a map xw = x o w (rightmost letter applied first).
"""
from fractions import Fraction as F

S = {'l': F(19, 200), 'm': F(1, 20), 'N': F(1, 19), 'K': F(1, 10)}
T = {'l': F(0), 'm': F(0), 'N': F(1, 19), 'K': F(1, 10)}
NK = set('NK')


def st(w):
    """(s_w, t_w) by composing maps; w(y) = s y + t."""
    s, t = F(1), F(0)
    for x in reversed(w):          # apply rightmost first
        s, t = S[x] * s, S[x] * t + T[x]
    return s, t


def psi(w):
    s, t = st(w)
    return t / s


def psi_terms(w):
    """Theorem 14.1 terms, in position order (leftmost N/K first)."""
    out = []
    for i, x in enumerate(w):
        if x in NK:
            p = F(1)
            for y in w[i + 1:]:
                p /= S[y]
            out.append(p)
    return out


def strip(n):
    for p in (2, 5, 19):
        while n % p == 0:
            n //= p
    return n


def fingerprint(q):
    q = F(q)
    if q == 0:
        return F(0)
    return F(strip(q.numerator), strip(q.denominator))


def canon(w):
    w = w.lstrip('lm')
    if not w:
        return ''
    w = 'N' + w[1:]
    if len(w) > 1 and w[1] == 'N':
        w = 'N' + 'm' + w[2:]
    while True:
        v = w.replace('mK', 'lN').replace('ml', 'lm')
        if v == w:
            break
        w = v
    return w.rstrip('lm')


def is_canon(w):
    return canon(w) == w


def val(n, p):
    """p-adic valuation of a nonzero rational."""
    n = F(n)
    a, b, v = n.numerator, n.denominator, 0
    while a % p == 0:
        a //= p; v += 1
    while b % p == 0:
        b //= p; v -= 1
    return v


def is_ratio_product(q):
    """Is q (a rational made of 2,5,19 only) a product of one or more word ratios
    19^L/(r 200^L 20^M), r in {19,10}?  Derivation: a product of k word ratios is
    19^(L-p) / (10^q 200^L 20^M) with p+q = k >= 1 (p = #r=19, q = #r=10)."""
    q = F(q)
    if q <= 0:
        return False
    rest = q
    e2, e5, e19 = val(q, 2), val(q, 5), val(q, 19)
    if rest != F(2) ** e2 * F(5) ** e5 * F(19) ** e19:
        return False
    A, B = -e2, -e5
    if 2 * B - A - e19 < 1:
        return False
    lo, hi = max(0, e19), min(A - B, 2 * B - A)
    return lo <= hi


def gen_canon(maxlen):
    """All canonical words of length 1..maxlen, by generating every word over the
    alphabet that starts with N and testing canon(w) == w (slow but honest) -- used
    only for small lengths; for larger lengths use gen_canon_local and cross-check."""
    from itertools import product
    out = []
    for n in range(1, maxlen + 1):
        for tail in product('lmNK', repeat=n - 1):
            w = 'N' + ''.join(tail)
            if canon(w) == w:
                out.append(w)
    return out


def gen_canon_local(maxlen):
    """Canonical words via the local description: starts with N, second letter not N,
    no 'mK', no 'ml', last letter N or K.  (Cross-checked against gen_canon.)"""
    out = []

    def rec(w):
        if w[-1] in NK:
            out.append(w)
        if len(w) == maxlen:
            return
        for x in 'lmNK':
            if len(w) == 1 and x == 'N':
                continue
            if w[-1] == 'm' and x in 'Kl':
                continue
            rec(w + x)
    rec('N')
    return out
