"""Where tesseract coordinates make things shorter, and where they don't.

Run: python3 examples.py   (needs sympy for the soliton check)
Every example prints a before/after and asserts they agree.
"""

import random
from fractions import Fraction
from itertools import product

import tesseract as t


def header(s):
    print("\n" + "=" * 70 + "\n" + s + "\n" + "=" * 70)


# ---------------------------------------------------------------- 1
def example_or_of_events():
    header("1. Chance that at least one of 4 independent things happens")
    p = [Fraction(1, 10), Fraction(1, 5), Fraction(1, 4), Fraction(1, 3)]
    # Before: inclusion-exclusion, 15 alternating terms (every non-empty corner)
    terms = []
    for mask in range(1, 16):
        prod_ = 1
        for i in range(4):
            if mask >> i & 1:
                prod_ *= p[i]
        sign = 1 if t.grade(mask) % 2 else -1
        terms.append(sign * prod_)
    before = sum(terms)
    # After: one circle-operation chain on -p, i.e. P = -((-p1) (+) (-p2) (+) ...)
    after = -t.multadd_all([-x for x in p])
    print(f"before: {len(terms)} signed terms, one per non-empty corner -> {before}")
    print(f"after : -( -p1 (+) -p2 (+) -p3 (+) -p4 ) -> {after}")
    assert before == after
    # Dependent events: put P(all of S happen) on corner S, one Mobius step.
    rng = random.Random(0)
    worlds = [tuple(rng.random() < 0.4 for _ in range(4)) for _ in range(2000)]
    w = [Fraction(sum(all(x[i] for i in range(4) if m >> i & 1) for x in worlds),
                  len(worlds)) for m in range(16)]
    none = sum((-1) ** t.grade(m) * w[m] for m in range(16))
    direct = Fraction(sum(any(x) for x in worlds), len(worlds))
    print(f"dependent events: 1 - signed corner sum = {1 - none} = direct count {direct}")
    assert 1 - none == direct


# ---------------------------------------------------------------- 2
def anf(table, n):
    """Truth table (index = input mask) -> corner weights (Mobius transform)."""
    return t.mobius(list(table))


def walsh(table, n):
    """Fourier-Walsh weights with +-1 inputs: f = sum_S f^(S) prod x_i."""
    N = 1 << n
    out = []
    for s in range(N):
        tot = 0
        for x in range(N):
            chi = (-1) ** t.grade(x & s)
            tot += table[x] * chi
        out.append(Fraction(tot, N))
    return out


def example_boolean():
    header("2. A 16-row truth table becomes a few corners")
    n = 4
    parity = [t.grade(x) % 2 for x in range(16)]
    pm = [1 - 2 * v for v in parity]          # parity as +1/-1
    W = walsh(pm, n)
    nz = {t.name(m, n): W[m] for m in range(16) if W[m]}
    print("parity of 4 bits: 16 rows -> Walsh corners", nz)
    assert nz == {"abcd": 1}
    maj3 = [1 if t.grade(x) >= 2 else 0 for x in range(8)]
    c = anf(maj3, 3)
    nz = {t.name(m, 3): c[m] for m in range(8) if c[m]}
    print("majority of 3: 8 rows -> corners", nz, "(ab+ac+bc-2abc)")
    assert nz == {"ab": 1, "ac": 1, "bc": 1, "abc": -2}
    # where the 0/1 basis is worse: parity in 0/1 corners fills everything
    c = anf(parity, 4)
    print("same parity in 0/1 corners:", {t.name(m, 4): c[m] for m in range(16) if c[m]})


# ---------------------------------------------------------------- 3
def example_kdv():
    header("3. KdV solitons: 16-term formula = one corner rule")
    import sympy as sp
    x, tt = sp.symbols("x t", real=True)
    ks = [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2)]
    n = len(ks)
    A = [[((ks[i] - ks[j]) / (ks[i] + ks[j])) ** 2 if i != j else 1
          for j in range(n)] for i in range(n)]
    w = t.soliton_tau_weights(A)
    eta = [k * x - k ** 3 * tt for k in ks]
    tau = sum(w[m] * sp.exp(sum(eta[i] for i in range(n) if m >> i & 1))
              for m in range(16))
    u = 2 * sp.diff(sp.log(tau), x, 2)
    resid = sp.diff(u, tt) + 6 * u * sp.diff(u, x) + sp.diff(u, x, 3)
    f = sp.lambdify((x, tt), resid, "mpmath")
    import mpmath
    mpmath.mp.dps = 40
    pts = [(mpmath.mpf(xv), mpmath.mpf(tv)) for xv in ("-6", "-2", "0", "1.5", "4")
           for tv in ("-1", "0", "0.7")]
    worst = max(abs(f(xv, tv)) for xv, tv in pts)
    print("tau has 16 terms; corner S weight = product of A_ij over pairs in S")
    print("KdV residual u_t + 6 u u_x + u_xxx at 15 sample points: max |r| =",
          mpmath.nstr(worst, 3))
    assert worst < 1e-25
    # control: drop the interaction (every A = 1) and it is no longer a solution
    w1 = [1] * 16
    tau1 = sum(w1[m] * sp.exp(sum(eta[i] for i in range(n) if m >> i & 1))
               for m in range(16))
    u1 = 2 * sp.diff(sp.log(tau1), x, 2)
    r1 = sp.diff(u1, tt) + 6 * u1 * sp.diff(u1, x) + sp.diff(u1, x, 3)
    f1 = sp.lambdify((x, tt), r1, "mpmath")
    ctrl = max(abs(f1(xv, tv)) for xv, tv in pts)
    print("control with every A = 1 (pure box, no interaction): max |r| =",
          mpmath.nstr(ctrl, 3))
    assert ctrl > 1e-3
    # the interaction lives only in grade 2: Mobius of log-weights
    logw = [sp.log(v) for v in w]
    m = t.mobius(logw)
    by_grade = {}
    for mask in range(16):
        g = t.grade(mask)
        val = sp.nsimplify(sp.simplify(m[mask]))
        by_grade.setdefault(g, set()).add(val == 0)
    print("Mobius of log-weights is zero at grades:",
          [g for g, s in by_grade.items() if s == {True}])
    assert by_grade[3] == {True} and by_grade[4] == {True}


# ---------------------------------------------------------------- 4
def example_rotation_is_worse():
    header("4. Where it gets worse: rotating a point")
    a, b = Fraction(3), Fraction(4)
    c = t.from_point([a, b])
    # rotate by the angle with cos = 3/5, sin = 4/5 (exact)
    co, si = Fraction(3, 5), Fraction(4, 5)
    a2, b2 = co * a - si * b, si * a + co * b
    c2 = t.from_point([a2, b2])
    print("corners before:", c, " after rotation:", c2)
    # The new 'ab' corner needs a^2 and b^2, which are not corners:
    rebuilt = co * si * (a * a - b * b) + (co * co - si * si) * c[3]
    assert c2[3] == rebuilt
    print("new ab = cos*sin*(a^2 - b^2) + (cos^2 - sin^2)*ab  -> needs a^2, b^2,")
    print("which no corner holds. Polar does this in one step (add the angle).")
    print("And -1 swallows everything: (-1) (+) b =", t.multadd(-1, Fraction(5)))


if __name__ == "__main__":
    example_or_of_events()
    example_boolean()
    example_kdv()
    example_rotation_is_worse()
    print("\nall checks passed")
