import sys, time; sys.path.insert(0, '.')
import sympy as sp
from pn2 import run
R = sp.Rational
for case, kw in [('GR', {}), ('EMD', dict(a=sp.sqrt(3))), ('EMD', dict(a=1)), ('EMD', dict(a=1/sp.sqrt(3))), ('EMD', dict(a=0))]:
    for M2 in [R(1, 1000), R(1, 3), 1]:
        t = time.time(); T, T2, T2e, rest = run(case, 1, M2, **kw); nu = M2/(1 + M2)**2
        print(case, kw, 'nu =', nu, '| T, T2, T2e =', T, T2, T2e, '| rest', rest, f'{time.time()-t:.1f}s', flush=True)
