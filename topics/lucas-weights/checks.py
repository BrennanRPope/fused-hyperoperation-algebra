"""Quick symbolic checks for the Lucas corner formula (sympy). Run: python3 checks.py"""
import sympy as sp
from itertools import combinations
k, d = sp.symbols('k d')
def U(m):
    a, b = sp.Integer(0), sp.Integer(1)
    if m == 0: return a
    for _ in range(m - 1): a, b = b, sp.expand(k*b - d*a)
    return b
def e(m, xs):
    return sum((sp.Mul(*c) for c in combinations(xs, m)), sp.Integer(0))
ok = 0
# 1. Corner formula: chain of n plane-law steps = sum U_m e_m / (1 - d sum U_{m-1} e_m), n = 1..5
for n in range(1, 6):
    xs = sp.symbols(f'x1:{n+1}')
    A, B = sp.Integer(0), sp.Integer(1)      # chain value A/B, start at 0
    for a in xs:
        A, B = sp.expand(A + a*(B + k*A)), sp.expand(B - d*a*A)
    num = sp.expand(sum(U(m)*e(m, xs) for m in range(n+1)))
    den = sp.expand(1 - d*sum(U(m-1)*e(m, xs) for m in range(1, n+1)))
    assert sp.expand(A - num) == 0 and sp.expand(B - den) == 0, n
    ok += 1
# 2. Pair recursion really is the law: x = A/B, y = F(x, a)
x, a = sp.symbols('x a')
F = (x + a + k*x*a) / (1 - d*x*a)
Ab, Bb = sp.symbols('A B')
assert sp.simplify((Ab + a*(Bb + k*Ab))/(Bb - d*a*Ab) - F.subs(x, Ab/Bb)) == 0; ok += 1
# 3. Cayley-Hamilton: N^m = U_m N - d U_{m-1} I, N = [[k,1],[-d,0]], m = 1..10
N = sp.Matrix([[k, 1], [-d, 0]])
for m in range(1, 11):
    assert sp.expand(N**m - (U(m)*N - d*U(m-1)*sp.eye(2))) == sp.zeros(2); ok += 1
# 4. Named points (first 8 weights U_1..U_8)
pts = {"(0,0) plain addition": (0,0), "(1,0) circle op": (1,0), "(-1,0) incl-excl": (-1,0),
       "(0,1) tangent": (0,1), "(0,-1) tanh/velocity": (0,-1), "(1,-1) Fibonacci": (1,-1),
       "(2,-1) Pell": (2,-1), "(1,-2) Jacobsthal": (1,-2), "(3,2) Mersenne": (3,2),
       "(2,1) parabola, U_m=m": (2,1), "(1,1) period 6": (1,1)}
for name, (kk, dd) in pts.items():
    print(f"{name:26s}", [U(m).subs({k: kk, d: dd}) for m in range(1, 9)])
print(f"\n{ok} symbolic checks passed")
