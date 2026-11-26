"""Numeric checks for the hyperoperation ladder (topics/hyperoperation-ladder).

The proofs live in lean-library/FusedAlgebra/Topics/HyperLadder.lean. These checks are evidence:
they recompute the proved facts by brute force, and they test the claims Lean does not reach yet
(right distributivity on the giant rungs, the height ladder at fractional heights).

Part A  whole rungs: lad n a b + 1 = H_n(a+1, b+1), Goodstein's sequence moved by x -> x+1.
Part B  the height ladder B_s(a,b) = E^s(E^-s a + E^-s b), E(x) = e^x - 1, regular iteration.
Part C  the dial k.
"""
import math, os, sys
from fractions import Fraction as Fr
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))

results = []
def check(name, ok, detail=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""))

# ---------------------------------------------------------------- Part A: whole rungs
@lru_cache(maxsize=None)
def hyp(n, x, y):
    """Goodstein: H_0(x,y)=y+1, H_1(x,0)=x, H_2(x,0)=0, H_n(x,0)=1 for n>=3, H_{n+1}(x,y+1)=H_n(x,H_{n+1}(x,y))."""
    if n == 0: return y + 1
    if y == 0: return x if n == 1 else (0 if n == 2 else 1)
    if n == 1: return x + y
    if n == 2: return x * y
    if n == 3: return x ** y
    return hyp(n - 1, x, hyp(n, x, y - 1))

@lru_cache(maxsize=None)
def lad(n, a, b):
    """The fused ladder by its own recursion: lad 0 a b = b+1, lad 1 a 0 = a+1, lad n a 0 = a (n>=2)."""
    if n == 0: return b + 1
    if b == 0: return a + 1 if n == 1 else a
    if n == 1: return a + b + 1
    if n == 2: return a + b + a * b
    if n == 3: return (a + 1) ** (b + 1) - 1
    return lad(n - 1, a, lad(n, a, b - 1))

def lad_raw(n, a, b):
    """No shortcuts, only the recursion (for small arguments)."""
    if n == 0: return b + 1
    if b == 0: return a + 1 if n == 1 else a
    return lad_raw(n - 1, a, lad_raw(n, a, b - 1))

SAFE = {0: 6, 1: 6, 2: 6, 3: 5, 4: 3, 5: 2}   # grid size per rung that stays computable
def grid(n):
    m = SAFE[n]
    pts = [(a, b) for a in range(m) for b in range(m)]
    if n == 4: pts = [(a, b) for (a, b) in pts if not (a >= 2 and b >= 2)] + [(1, 2), (2, 1)]
    if n == 5: pts = [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (2, 0), (1, 2)]
    return pts

ok = all(lad(n, a, b) + 1 == hyp(n, a + 1, b + 1) for n in range(6) for (a, b) in grid(n))
check("A1 every rung is Goodstein's H_n moved by x -> x+1 (n = 0..5)", ok)
ok = all(lad_raw(n, a, b) == lad(n, a, b) for n in range(4) for a in range(4) for b in range(4))
check("A1 the closed forms a+b+1, a+b+ab, (1+a)^(1+b)-1 agree with the bare recursion", ok)

def arrow(m, a, b):
    if m == 0: return a ** b
    if b == 0: return 1
    return arrow(m - 1, a, arrow(m, a, b - 1))
ok = all(hyp(m + 3, a, b) == arrow(m, a, b) for m in range(3) for a in range(1, 4) for b in range(0, 4)
         if not (m >= 1 and a >= 3 and b >= 3) and not (m == 2 and (a, b) not in [(1, 3), (2, 2), (2, 3), (3, 1), (1, 1)]))
check("A2 H_(m+3) = Knuth's m+1 up-arrows", ok, "2^3=8, 2^^3=16, 2^^^3=65536 included")
check("A2 values 2^^^3 = 65536, H_n(2,2) = 4 on every rung n >= 1",
      arrow(2, 2, 3) == 65536 and all(hyp(n, 2, 2) == 4 for n in range(1, 8)))

def comm_fails(n): return [(a, b) for (a, b) in grid(n) if (b, a) in grid(n) and lad(n, a, b) != lad(n, b, a)]
check("A3 commutative exactly at rungs 1 and 2 on the grid",
      [n for n in range(6) if not comm_fails(n)] == [1, 2])
def assoc_fail(n):
    if n >= 4: return lad(n, lad(n, 1, 0), 1) != lad(n, 1, lad(n, 0, 1))   # the Lean witness (1,0,1)
    return any(lad(n, lad(n, a, b), c) != lad(n, a, lad(n, b, c)) for a in range(4) for b in range(4) for c in range(3))
check("A3 associative exactly at rungs 1 and 2", [n for n in range(6) if not assoc_fail(n)] == [1, 2])
check("A3 the corner 1*1 = 3 on every rung n >= 1", all(lad(n, 1, 1) == 3 for n in range(1, 8)))
check("A3 0 is a right identity from rung 2 up, a left absorber from rung 3 up",
      all(lad(n, a, 0) == a for n in range(2, 6) for a in range(5)) and all(lad(n, 0, b) == 0 for n in range(3, 8) for b in range(3))
      and lad(1, 3, 0) != 3 and lad(2, 0, 3) != 0)

def left_dist(n, a, b, c): return lad(n + 1, a, lad(n, b, c)) == lad(n, lad(n + 1, a, b), lad(n + 1, a, c))
def right_dist(n, a, b, c): return lad(n + 1, lad(n, a, b), c) == lad(n, lad(n + 1, a, c), lad(n + 1, b, c))
ok = all(left_dist(n, a, b, c) for n in (0, 1) for a in range(5) for b in range(5) for c in range(5))
ok = ok and all(not left_dist(n, 1, 0, 0) for n in range(2, 7))
check("A4 rung n+1 distributes over rung n on the left only for n = 0, 1 (witness 1,0,0 kills the rest)", ok)
ok = all(right_dist(n, a, b, c) for n in (0, 1, 2) for a in range(4) for b in range(4) for c in range(4))
check("A4 rung n+1 distributes over rung n on the right for n = 0, 1, 2", ok)
check("A4 ...and fails at rung 4 over rung 3 (Lean witness 1,2,1)", not right_dist(3, 1, 2, 1))

# A5: right distributivity on the giant rungs, refuted by residues.
# Witness (a,b,c) = (1,2,1), i.e. A,B,C = 2,3,2 in Goodstein's coordinates. With C = 2 the two sides are
#   LHS = H_n(P, P),  P = H_n(2,3)        RHS = H_n(H_n(2,2), H_n(3,3)) = H_n(4, H_n(3,3)).
# For n >= 4 both are power towers: LHS = tower of P's, RHS = tower of 4's. Towers mod m are computable by
# the generalized Euler theorem: a^E = a^((E mod phi m) + phi m) (mod m) once E >= log2 m. Different residues
# mean different numbers.
from sympy import totient
BIG = 1 << 64

def tower(base, h, m):
    """(base^^h mod m, exact value or None if >= 2^64). base is a function q -> (residue mod q, exact or None)."""
    if m == 1:
        r, ex = tower(base, h, 2) if h else (0, 1)
        return 0, ex
    if h == 0: return 1 % m, 1
    br, bex = base(m)
    if h == 1: return br, bex
    ph = int(totient(m))
    er, eex = tower(base, h - 1, ph)
    if eex is not None and bex is not None:
        # exact if small enough
        if eex * max(bex.bit_length(), 1) < 64:
            v = bex ** eex
            return v % m, (v if v < BIG else None)
        return pow(bex, eex, m), None
    if eex is not None and eex < 64:
        return pow(br, eex, m), None if (bex is None or eex > 1) else bex
    return pow(br, er + ph, m), None   # exponent is huge: generalized Euler

def const(x): return lambda q: (x % q, x if x < BIG else None)
def stable_tower(base, m):
    """residue of base^^H for every H large enough: check that the value has settled."""
    vals = [tower(base, h, m)[0] for h in (40, 50, 60)]
    assert len(set(vals)) == 1
    return vals[0]

# rung 5 over rung 4 (n = 4): LHS = 16^^16, RHS = 4^^(3^27)  (3^27 is past every height we need)
moduli = [m for m in range(3, 2000)]
m4 = next(m for m in moduli if tower(const(16), 16, m)[0] != stable_tower(const(4), m))
check("A5 rung 5 over rung 4: the two sides of (1,2,1) differ mod %d" % m4, True,
      "16^^16 = %d, 4^^(3^27) = %d (mod %d)" % (tower(const(16), 16, m4)[0], stable_tower(const(4), m4), m4))
# rung 6 over rung 5 (n = 5): P = H_5(2,3) = 65536; LHS = 65536^^(huge), RHS = 4^^(huge)
m5 = next(m for m in moduli if stable_tower(const(65536), m) != stable_tower(const(4), m))
check("A5 rung 6 over rung 5: the two sides differ mod %d" % m5, True)
# rung n+1 over n for every n >= 6: P = H_n(2,3) is itself a tower of 2's of height >= 65536
two_inf = lambda q: (stable_tower(const(2), q), None)
m6 = next(m for m in moduli if stable_tower(two_inf, m) != stable_tower(const(4), m))
check("A5 rung n+1 over rung n, every n >= 6: one modulus (%d) separates the two sides" % m6, True,
      "computer check by residues, not in Lean")

# ---------------------------------------------------------------- Part B: the height ladder
import regular_iteration as RI
E = math.expm1
def Bs(s, a, b): return RI.B(s, a, b)
tol = 1e-9
check("B1 height 0 is +, height 1 is circ", abs(Bs(0, 0.7, 1.3) - 2.0) < tol and abs(Bs(1, 0.7, 1.3) - (0.7 + 1.3 + 0.91)) < 1e-8)
B2 = lambda a, b: (1 + a) * (1 + b) * math.exp(math.log1p(a) * math.log1p(b)) - 1
check("B1 height 2 is (1+a)(1+b)e^(L a L b) - 1, height -1 is log(e^a+e^b-1)",
      abs(Bs(2, 0.5, 0.8) - B2(0.5, 0.8)) < 1e-7 and abs(Bs(-1, 0.5, 0.8) - math.log(math.exp(0.5) + math.exp(0.8) - 1)) < 1e-8)
pts = [(0.3, 0.9, 1.4), (1.2, 0.4, 0.7)]
ok = all(abs(Bs(s, a, b) - Bs(s, b, a)) < 1e-8 for s in (-0.5, 0.5, 1.5) for (a, b, _) in pts)
check("B2 commutative at fractional heights -1/2, 1/2, 3/2", ok)
ok = all(abs(Bs(s, Bs(s, a, b), c) - Bs(s, a, Bs(s, b, c))) < 1e-7 for s in (-0.5, 0.5, 1.5) for (a, b, c) in pts)
check("B2 associative at fractional heights", ok)
check("B2 0 is the identity at every height", all(abs(Bs(s, 0.83, 0.0) - 0.83) < 1e-9 for s in (-0.5, 0.5, 2.5)))
fexp = lambda a, b: (1 + a) ** (1 + b) - 1
check("B3 the two ladders part after circ: height 2 is not rung 3, height 0 is not rung 1",
      abs(Bs(2, 1.0, 1.0) - fexp(1, 1)) > 0.5 and abs(Bs(0, 1.0, 1.0) - 3.0) > 0.5,
      "B_2(1,1) = %.6f vs 1*1 = 3 on rung 3" % Bs(2, 1.0, 1.0))
# Bennett 1915 (shifted by x -> 1+x): F_2(x,y) = exp(log x log y). Agrees with the E ladder only up to height 1.
bennett2 = lambda a, b: math.exp(math.log1p(a) * math.log1p(b)) - 1
check("B3 the E-ladder is not Bennett's commutative ladder at height 2", abs(B2(1, 1) - bennett2(1, 1)) > 1)
# partner distributivity: B_(s+1) distributes over (E^s)_*(sadd) = E^s(E^-s a + E^-s b + 1)
s = 0.5
part = lambda a, b: RI.Es(s, RI.Es(-s, a) + RI.Es(-s, b) + 1)
a, b, c = 0.4, 0.3, 0.6
lhs = Bs(s + 1, a, part(b, c)); rhs = part(Bs(s + 1, a, b), Bs(s + 1, a, c))
check("B4 height 3/2 distributes over its partner (E^1/2)_*(sadd)", abs(lhs - rhs) < 1e-7, "%.10f vs %.10f" % (lhs, rhs))
lhs = Bs(s + 1, a, Bs(s, b, c)); rhs = Bs(s, Bs(s + 1, a, b), Bs(s + 1, a, c))
check("B4 ...but not over height 1/2 itself", abs(lhs - rhs) > 1e-3)
lnh = lambda s, X: RI.Es(-s, X - 1)
X, Y = 1.7, 2.3
check("B5 hypernatural log at height 1/2: ln_s(XY) = B_(1-s)(ln_s X, ln_s Y)",
      abs(lnh(0.5, X * Y) - Bs(0.5, lnh(0.5, X), lnh(0.5, Y))) < 1e-8)

# ---------------------------------------------------------------- Part C: the dial k
def saddK(k, a, b): return a + b + Fr(1) / k
def circK(k, a, b): return a + b + k * a * b
def fexpK(k, a, m): return ((1 + k * a) ** m - 1) / k
ok = True
for k in (Fr(2), Fr(-3), Fr(1, 5)):
    for (a, b, c) in [(Fr(1), Fr(2), Fr(3)), (Fr(1, 2), Fr(-1, 3), Fr(2))]:
        ok &= k * saddK(k, a, b) == saddK(1, k * a, k * b)          # rescaling of rung 1
        ok &= k * circK(k, a, b) == circK(1, k * a, k * b)          # rescaling of rung 2
        ok &= k * fexpK(k, a, 3) == fexpK(1, k * a, 3)              # rescaling of rung 3
        ok &= circK(k, a, saddK(k, b, c)) == saddK(k, circK(k, a, b), circK(k, a, c))
check("C1 every dial k != 0 is the k = 1 ladder rescaled by x -> kx (rungs 1-3), distributivity kept", ok)
ks = [Fr(1, 10 ** j) for j in (2, 4, 6)]
r1 = [float(saddK(k, 1, 1)) for k in ks]; r2 = [float(circK(k, Fr(1), Fr(1))) for k in ks]; r3 = [float(fexpK(k, Fr(1), 2)) for k in ks]
check("C2 at k -> 0 the ladder collapses: rung 1 blows up, rung 2 -> +, rung 3 -> m*a",
      r1[-1] > 1e5 and abs(r2[-1] - 2) < 1e-5 and abs(r3[-1] - 2) < 1e-5, "k=1e-6: %.3g, %.8f, %.8f" % (r1[-1], r2[-1], r3[-1]))

print()
n_ok = sum(results)
print(f"{n_ok}/{len(results)} checks passed")
