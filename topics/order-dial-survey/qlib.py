"""Small exact tools for the order-dial survey (Nov 26, 2026).

- The quantum plane: elements are dicts {(a, b): coeff} meaning sum coeff * x^a y^b (normal order, x left of y),
  with the swap rule y x = q x y, so (x^a y^b)(x^c y^d) = q^(b c) x^(a+c) y^(b+d).  This is the repo's swapExp
  (CutIndependence: swapExp b c = b * c).
- Multiparameter words: products of generators x_1..x_n with x_j x_i = q x_i x_j for i < j.
- Gaussian binomials, q-integers, inversions.
Coefficients are sympy expressions in the symbol q (and k, t where needed)."""
import itertools, functools
import sympy as sp

q, k, t, hb = sp.symbols('q k t hbar')


# ---------------- quantum plane ----------------
def qp_mul(A, B, qq=q):
    out = {}
    for (a, b), c1 in A.items():
        for (c, d), c2 in B.items():
            key = (a + c, b + d)
            out[key] = out.get(key, 0) + c1 * c2 * qq ** (b * c)
    return qp_clean(out)


def qp_add(*As):
    out = {}
    for A in As:
        for key, c in A.items():
            out[key] = out.get(key, 0) + c
    return qp_clean(out)


def qp_scale(A, s):
    return qp_clean({key: s * c for key, c in A.items()})


def qp_clean(A):
    out = {}
    for key, c in A.items():
        c = sp.expand(c)
        if c != 0:
            out[key] = c
    return out


def qp_pow(A, n, qq=q):
    R = {(0, 0): sp.Integer(1)}
    for _ in range(n):
        R = qp_mul(R, A, qq)
    return R


def circ(A, B, kk=k, qq=q):
    """a o_k b = a + b + k a b, on quantum-plane elements (the circle operation of the ring)."""
    return qp_add(A, B, qp_scale(qp_mul(A, B, qq), kk))


X = {(1, 0): sp.Integer(1)}
Y = {(0, 1): sp.Integer(1)}


# ---------------- q-numbers ----------------
def qint(n, qq=q):
    return sum((qq ** i for i in range(n)), sp.Integer(0))


def qfact(n, qq=q):
    r = sp.Integer(1)
    for i in range(1, n + 1):
        r *= qint(i, qq)
    return sp.expand(r)


@functools.lru_cache(maxsize=None)
def gauss(n, m, qq=q):
    """Gaussian binomial [n, m]_q as a polynomial (Pascal rule [n,m] = [n-1,m-1] + q^m [n-1,m])."""
    if m < 0 or m > n:
        return sp.Integer(0)
    if m == 0 or m == n:
        return sp.Integer(1)
    return sp.expand(gauss(n - 1, m - 1, qq) + qq ** m * gauss(n - 1, m, qq))


def inv(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def words(n, m):
    """0/1 words with m ones (y's) and n - m zeros (x's)."""
    for pos in itertools.combinations(range(n), m):
        yield tuple(1 if i in pos else 0 for i in range(n))


# ---------------- multiparameter words ----------------
def sort_word(word, qs):
    """Normal-order a word in generators 0..n-1 with x_j x_i = qs[i][j] x_i x_j (i < j).
    Returns (factor, sorted tuple)."""
    w = list(word)
    f = sp.Integer(1)
    changed = True
    while changed:
        changed = False
        for i in range(len(w) - 1):
            if w[i] > w[i + 1]:
                f *= qs[w[i + 1]][w[i]]
                w[i], w[i + 1] = w[i + 1], w[i]
                changed = True
    return f, tuple(w)
