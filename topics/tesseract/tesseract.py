"""Tesseract coordinates for Fused Hyperoperation Algebra.

A point in n dimensions gets 2**n coordinates, one per corner of the
n-cube (n = 4 is the tesseract). A corner is addressed by n switches,
one per axis: switch i is on when axis i takes part. The corner's
grade (its dimension) is how many switches are on.

Corners are stored as bitmasks: bit 0 is axis a, bit 1 is axis b, and so
on. A coordinate vector is a plain list of length 2**n, indexed by mask.

The three ideas, in code:
  from_point(a)  -> corner coordinates  c_S = product of a_i for i in S
  zeta(c)        -> box coordinates     B_T = sum of c_S over S inside T
                                             = product of (1 + a_i), i in T
  compose(v, w)  -> circle operation on corners (switches combine by OR)

∘₁ a (+) b = a + b + ab = (1+a)(1+b) - 1. In box coordinates it is
plain multiplication, corner by corner; mobius() goes back.

Everything works with ints, floats or fractions.Fraction.
"""

from fractions import Fraction
from itertools import combinations
from math import comb

AXES = "abcdefgh"


# ---------------------------------------------------------------- numbers

def multadd(x, y):
    """The circle operation of two numbers: x + y + x*y."""
    return x + y + x * y


def multadd_all(values):
    """The circle operation of a whole list, left to right (0 is the identity)."""
    total = 0
    for v in values:
        total = multadd(total, v)
    return total


def multadd_inverse(x):
    """The y with x (+) y = 0, that is -x / (1 + x). Undefined at x = -1."""
    if x == -1:
        raise ZeroDivisionError("-1 has no circle inverse (its box side is 0)")
    return -x / (1 + x) if isinstance(x, float) else Fraction(-x) / (1 + x)


# ---------------------------------------------------------------- corners

def size(n):
    return 1 << n


def grade(mask):
    """Number of switches on = the corner's dimension."""
    return bin(mask).count("1")


def label(mask, n):
    """Switch address in axis order a b c d..., e.g. '1010' for a and c."""
    return "".join("1" if mask >> i & 1 else "0" for i in range(n))


def name(mask, n):
    """Letters of the switches that are on, e.g. 'ac'; the empty corner is '1'."""
    s = "".join(AXES[i] for i in range(n) if mask >> i & 1)
    return s or "1"


def corners_by_grade(n):
    """List of lists: corners_by_grade(n)[k] holds every corner of grade k."""
    layers = [[] for _ in range(n + 1)]
    for mask in range(size(n)):
        layers[grade(mask)].append(mask)
    return layers


def subsets(mask):
    """Every sub-corner of mask (every mask whose switches are inside it)."""
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & mask


# ------------------------------------------------- ordinary <-> tesseract

def from_point(a):
    """Corner coordinates of an ordinary point a = (a_1..a_n).

    Corner S holds the product of a_i over the switches in S: the measure
    of that piece when a box with sides 1 + a_i is cut at 1 on every axis.
    """
    n = len(a)
    c = [1] * size(n)
    for mask in range(1, size(n)):
        low = mask & -mask
        i = low.bit_length() - 1
        c[mask] = c[mask ^ low] * a[i]
    return c


def to_point(c):
    """Back to ordinary coordinates: the grade-1 corners, scaled by corner 1."""
    n = size_to_n(len(c))
    if c[0] == 0:
        raise ValueError("corner 1 is 0, so there is no scale to read from")
    return [c[1 << i] / c[0] if isinstance(c[0], float) else Fraction(c[1 << i]) / c[0]
            for i in range(n)]


def is_pure(c, tol=0):
    """True when c comes from an ordinary point (every corner factors)."""
    n = size_to_n(len(c))
    if c[0] == 0:
        return all(abs(x) <= tol for x in c)
    a = to_point(c)
    rebuilt = from_point(a)
    return all(abs(c[m] - c[0] * rebuilt[m]) <= tol for m in range(len(c)))


def size_to_n(length):
    n = length.bit_length() - 1
    if length != 1 << n:
        raise ValueError("a coordinate list must have length 2**n")
    return n


def evaluate(c, x):
    """Value of the coordinates at x: sum over corners of c_S times prod x_i."""
    n = size_to_n(len(c))
    if len(x) != n:
        raise ValueError("x needs one value per axis")
    monomials = from_point(x)
    return sum(c[m] * monomials[m] for m in range(len(c)))


def total(c):
    """Sum of every corner. For a pure point this is 1 + (a (+) b (+) ...)."""
    return sum(c)


def grade_totals(c):
    """Sum of the corners in each grade: [point, lines, areas, volumes, ...]."""
    n = size_to_n(len(c))
    out = [0] * (n + 1)
    for m, v in enumerate(c):
        out[grade(m)] += v
    return out


# ------------------------------------------------- transforms and the circle operation

def zeta(c):
    """Box coordinates: B_T = sum of c_S over every S inside T."""
    b = list(c)
    n = size_to_n(len(b))
    for i in range(n):
        bit = 1 << i
        for m in range(len(b)):
            if m & bit:
                b[m] += b[m ^ bit]
    return b


def mobius(b):
    """Inverse of zeta (Mobius inversion): c_S = sum over T inside S of (-1)^|S-T| B_T."""
    c = list(b)
    n = size_to_n(len(c))
    for i in range(n):
        bit = 1 << i
        for m in range(len(c)):
            if m & bit:
                c[m] -= c[m ^ bit]
    return c


def compose(v, w):
    """The circle operation of two tesseract points: corners combine by OR of switches.

    (v * w)_T = sum of v_S * w_U over every pair with S | U == T.
    Same as mobius(zeta(v) * zeta(w)) corner by corner.
    """
    if len(v) != len(w):
        raise ValueError("both points need the same number of axes")
    out = [0] * len(v)
    for s, vs in enumerate(v):
        if vs == 0:
            continue
        for u, wu in enumerate(w):
            out[s | u] += vs * wu
    return out


def compose_fast(v, w):
    """compose() through box coordinates; O(n 2^n) instead of O(4^n)."""
    return mobius([x * y for x, y in zip(zeta(v), zeta(w))])


# ------------------------------------------------- the soliton example

def soliton_tau_weights(A):
    """Hirota N-soliton corner weights: w_S = product of A[i][j], i < j in S.

    A[i][j] is the interaction between solitons i and j. With every
    A[i][j] = 1 the point is pure: the tau function factors as
    (1 + e^eta_1)(1 + e^eta_2)...
    """
    n = len(A)
    w = [1] * size(n)
    for mask in range(size(n)):
        on = [i for i in range(n) if mask >> i & 1]
        for i, j in combinations(on, 2):
            w[mask] *= A[i][j]
    return w


def describe(c):
    """Readable listing, one line per corner, grouped by grade."""
    n = size_to_n(len(c))
    lines = []
    for k, layer in enumerate(corners_by_grade(n)):
        parts = [f"{label(m, n)} {name(m, n)}={c[m]}" for m in layer]
        lines.append(f"grade {k} ({comb(n, k)}): " + ", ".join(parts))
    return "\n".join(lines)


if __name__ == "__main__":
    a = [Fraction(1, 2), 2, 3, Fraction(-1, 3)]
    c = from_point(a)
    print(describe(c))
    print("total =", total(c), "= 1 + circle operation =", 1 + multadd_all(a))
