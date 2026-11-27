"""The border of the two-base product: the jet q-MZV as value plus q-derivative.  (DRAFT checks)
Run: python3 border_checks.py   (sympy exact where possible, numpy/mpmath for limits q -> 1)

Asked by Brennan, Nov 27 05:31Z: "Yes" to "Should I test the border case next, where the jet holds a q-MZV and its
q-derivative?"
"""
import json, itertools, numpy as np, mpmath as mp, sympy as sp
from fractions import Fraction as Fr
mp.mp.dps = 30
res = []
def ok(name, c, note=""):
    res.append({"id": name, "pass": bool(c), "note": note}); print(("PASS " if c else "FAIL ") + name)

z = {s: float(mp.zeta(s)) for s in range(2, 7)}
Z22, Z31, Z21, Z32 = 0.75*z[4], z[4]/4, z[3], 3*z[2]*z[3] - 11*z[5]/2   # Euler's evaluations

# ---- numerics for Bradley's q-MZVs at q < 1 (floats) ----
def letters(s, q, Nmax):
    n = np.arange(1, Nmax + 1, dtype=np.float64); qn = np.exp(n*np.log(q))
    return (1 - q)**s * qn**(s - 1) / (1 - qn)**s
def brad(word, q):
    Nmax = int(90/(1 - q)) + 10
    acc = np.ones(Nmax)                      # running nested sum, innermost first
    for s in reversed(word):
        L = letters(s, q, Nmax) * acc
        acc = np.concatenate([[0.0], np.cumsum(L)[:-1]])
    return float(np.sum(L))
def slope_at_1(word, value):
    """d zeta_q / dq at q = 1-, by Richardson on (value - zeta_q(1-h))/h"""
    hs = [4e-3, 2e-3, 1e-3, 5e-4]
    sl = [(value - brad(word, 1 - h))/h for h in hs]
    return 2*sl[-1] - sl[-2]

# B1 away from q = 1 on the border, the jet's second slot is the ordinary derivative, so a relation that holds for all q
#    differentiates to a true relation: Bradley's q-analogue of Euler, zeta_q(2,1) = zeta_q(3), and its q-derivative
B1 = True; vals = []
for q0 in (0.3, 0.5, 0.8):
    h = 1e-6
    e = brad((2, 1), q0) - brad((3,), q0)
    de = (brad((2, 1), q0 + h) - brad((2, 1), q0 - h) - brad((3,), q0 + h) + brad((3,), q0 - h))/(2*h)
    B1 &= abs(e) < 1e-13 and abs(de) < 1e-7; vals.append((q0, e, de))
ok("B1 off q = 1 the jet adds nothing new: zeta_q(2,1) = zeta_q(3) (Bradley's q-Euler) and so do their q-derivatives (q = 0.3, 0.5, 0.8)", B1)

# B2 at the origin of the dial (both bases 1, k = d = 0, N^2 = 0) the jet of a letter is its q-derivative at q = 1:
#    d/dq [q^((s-1)n) / [n]_q^s] at q = 1 equals (s/2 - 1)/n^(s-1) + (s/2)/n^s
qs = sp.symbols("q")
B2 = all(sp.simplify(sp.diff(qs**((s-1)*n) / sum(qs**j for j in range(n))**s, qs).subs(qs, 1)
                     - (sp.Rational(s, 2) - 1)/sp.Integer(n)**(s-1) - sp.Rational(s, 2)/sp.Integer(n)**s) == 0
         for n in range(1, 9) for s in range(1, 6))
ok("B2 the jet of a letter at q = 1: d/dq z_s(n) = (s/2 - 1)/n^(s-1) + (s/2)/n^s  (n = 1..8, s = 1..5)", B2,
   "so the jet's second slot is delta = (weight/2) * identity + a lowering operator z_s -> (s/2 - 1) z_(s-1)")

# B3 the second slot obeys the second slot of the jet quasi-shuffle (N = eps, eps^2 = 0); exact on truncated sums
def Hd(word, M):     # (value, formal derivative at q = 1) of the truncated sum, exact
    v = Fr(0); dv = Fr(0)
    for ns in itertools.combinations(range(M, 0, -1), len(word)):
        a = Fr(1); da = Fr(0)
        for s, n in zip(word, ns):
            l = Fr(1, n**s); dl = Fr(s - 2, 2)*Fr(1, n**(s-1)) + Fr(s, 2)*Fr(1, n**s)
            a, da = a*l, a*dl + da*l
        v += a; dv += da
    return v, dv
B3 = True
for M in range(1, 9):
    (a, da), (b, db), (c, dc), (e, de) = Hd((2,), M), Hd((2, 2), M), Hd((4,), M), Hd((3,), M)
    B3 &= a*a == 2*b + c and 2*a*da == 2*db + dc - e
ok("B3 second slot of H(2)H(2) = 2H(2,2) + H(4) - N H(3) at the origin: 2 H(2) dH(2) = 2 dH(2,2) + dH(4) - H(3)  (M = 1..8, exact)", B3)

# B4 the jet's slope versus the true slope of the infinite q-series at q = 1-: they agree for s >= 3, and at s = 2
#    the true slope is larger by exactly 1/2
B4 = True; tab = []
for s in (2, 3, 4, 5):
    true = slope_at_1((s,), z[s]); formal = (s/2 - 1)*(z[s-1] if s > 2 else 0) + (s/2)*z[s]
    tab.append((s, true, formal, true - formal))
    B4 &= abs(true - formal - (0.5 if s == 2 else 0)) < 5e-6
ok("B4 true slope of zeta_q(s) at q = 1 equals the jet's (s/2 - 1)zeta(s-1) + (s/2)zeta(s) for s = 3, 4, 5; at s = 2 it is zeta(2) + 1/2, "
   "and the jet alone gives zeta(2)", B4, "; ".join(f"s={s}: true {t:.7f}, jet {f:.7f}, gap {g:+.7f}" for s, t, f, g in tab))

# B5 where the 1/2 comes from: the lowering coefficient (s/2 - 1) is zero at s = 2 but zeta(s - 1) has a pole there;
#    continued in s, (s/2 - 1) zeta(s - 1) -> 1/2.  In depth 2 the same residue gives (1/2) zeta(rest) for a front 2
lim = mp.limit(lambda s: (s/2 - 1)*mp.zeta(s - 1), 2)
B5a = abs(float(lim) - 0.5) < 1e-20
ok("B5a (s/2 - 1) zeta(s - 1) -> 1/2 as s -> 2: the 1/2 is a residue, zero times a pole", B5a)
t22 = slope_at_1((2, 2), Z22); p22 = 2*Z22 + z[2]/2; f22 = 2*Z22
t31 = slope_at_1((3, 1), Z31); p31 = 2*Z31 + 0.5*Z21 - 0.5*(z[2] - z[3])
t32 = slope_at_1((3, 2), Z32); p32 = 2.5*Z32 + 0.5*Z22   # w/2 zeta(3,2) + (1/2) zeta(2,2); the back 2 lowers with coefficient 0
B5b = abs(t22 - p22) < 5e-6 and abs(t22 - f22) > 0.5
B5c = abs(t31 - p31) < 5e-6 and abs(t32 - p32) < 5e-6
ok("B5b depth 2, front 2: true slope of zeta_q(2,2) = 2 zeta(2,2) + zeta(2)/2 (jet alone: 2 zeta(2,2)): %.7f vs %.7f" % (t22, p22), B5b)
ok("B5c depth 2, no front 2: true slope = jet slope for (3,1) and (3,2): %.7f vs %.7f, %.7f vs %.7f" % (t31, p31, t32, p32), B5c)

# B6 consistency with q-Euler: the rule (jet + "front 2 gives half of the rest") gives the same slope for zeta(2,1) and
#    zeta(3), whatever value T is used for the divergent zeta(1)
T, z2s, z3s = sp.symbols("T zeta2 zeta3")
rule21 = sp.Rational(3, 2)*z3s + sp.Rational(1, 2)*T + sp.Rational(-1, 2)*(T - z2s)    # w/2 zeta(2,1) + front-2 residue (1/2)zeta(1) + lower the 1: -(1/2) zeta(2,0), zeta(2,0) = zeta(1) - zeta(2)
rule3 = sp.Rational(1, 2)*z2s + sp.Rational(3, 2)*z3s
B6 = sp.simplify(rule21 - rule3) == 0 and abs(slope_at_1((2, 1), Z21) - (z[2]/2 + 1.5*z[3])) < 5e-6
ok("B6 the rule passes Bradley's q-Euler: slope of zeta_q(2,1) = slope of zeta_q(3) = zeta(2)/2 + 3 zeta(3)/2, and the divergent zeta(1) cancels", B6)

# B7 negative control: without the residue term the depth-2 front-2 slope misses by zeta(2)/2
B7 = abs(t22 - f22 - z[2]/2) < 5e-6
ok("B7 negative control: the plain jet slope of zeta_q(2,2) misses the true one by zeta(2)/2 = %.7f" % (z[2]/2), B7)

json.dump(res, open("border_checks.json", "w"), indent=1)
print(f"{sum(r['pass'] for r in res)}/{len(res)} pass")
