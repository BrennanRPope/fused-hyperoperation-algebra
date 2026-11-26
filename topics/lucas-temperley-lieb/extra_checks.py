"""A single crossing is itself a plane-dial element: sigma = A + A^-1 e in TL_2 obeys sigma^2 = k sigma - d."""
import sys, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from tl import TL
A = sp.symbols("A", nonzero=True)
delta = -A**2 - A**-2
T = TL(2, delta)
sig = {T.one: A, T.e(0): 1/A}
k, d = A - A**-3, -A**-2
diff = T.add(T.mul(sig, sig), T.add(T.scal(k, sig), {T.one: -d}), -1)
ok = all(sp.simplify(v) == 0 for v in diff.values())
print(("PASS" if ok else "FAIL"), "Kauffman crossing sigma = A + A^-1 e (loop value -A^2 - A^-2) obeys sigma^2 = k sigma - d with k = A - A^-3, d = -A^-2 (the Hecke relation)")
ratio = sp.simplify(k**2/d)
print("   invariant k^2/d =", sp.factor(ratio), "; discriminant k^2 - 4d =", sp.factor(sp.simplify(k**2 - 4*d)), "(zero only at A^4 = -1)")
