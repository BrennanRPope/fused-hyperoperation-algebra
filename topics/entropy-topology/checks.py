"""Checks for topics/entropy-topology (DRAFT, Nov 27 2026).
Run: python3 checks.py  (needs numpy, scipy, sympy, mpmath; about 2 minutes).
Recomputes windows_index.json and jensen.json, then tests everything on the page."""
import json, math, subprocess, sys
from math import gcd
import numpy as np, sympy as sp, mpmath as mp

subprocess.run([sys.executable, 'windows_index.py'], check=True, capture_output=True)
subprocess.run([sys.executable, 'jensen.py'], check=True, capture_output=True)

from dialtop import (dial_power_h, rp1_degree, rp1_degree_alg, newton_parts, distinct_real_roots,
                     newton_preimage_counts, newton_turning_counts, real_preimages)
from orbits import (logistic_iter, interval_index, orbits_of_exact_period, pattern, birkhoff_pair,
                    linear_type, find_henon_orbits, henon_n, winding_index)
from morse_graph import morse_graph, set_centre

results = []
def check(name, ok, info=''):
    results.append(bool(ok)); print(('PASS' if ok else 'FAIL'), name, info)

EXPECT = {'elliptic': (1, 1), 'hyperbolic': (-1, -1), 'inverse hyperbolic': (1, -1)}
rng = np.random.default_rng(20261127)

# ---------------- A. The dial square is Newton's method ----------------
x, y, k, d, e1, e2 = sp.symbols('x y k d e1 e2')
F = (2 * x + k * x ** 2) / (1 - d * x ** 2)
N = (y ** 2 - d) / (2 * y + k)
check('A1 x o x = 1/N(1/x): the dial square is Newton\'s method on q(y) = y^2 + k y + d (exact)',
      sp.simplify(F - 1 / N.subs(y, 1 / x)) == 0)
kk, dd = -(e1 + e2), e1 * e2
q = y ** 2 + kk * y + dd
def comp(a, b): return sp.cancel((a + b + kk * a * b) / (1 - dd * a * b))
X = x; okA2 = okA3 = True
for m in range(2, 7):
    X = comp(X, x)
    Gy = sp.cancel(1 / X.subs(x, 1 / y))
    g = 1 / q
    K = y + (m - 1) * sp.diff(g, y, m - 2) / sp.diff(g, y, m - 1)
    okA2 &= sp.simplify(Gy - K) == 0
    M = lambda t: (t - e1) / (t - e2)
    okA3 &= sp.simplify(M(Gy) - M(y) ** m) == 0
check('A2 m-fold self-composition = Koenig\'s root-finding method of order m on q, m = 2..6 (m = 3 is Halley)', okA2)
check('A3 under z = (y - e1)/(y - e2) every m-fold power is z -> z^m (e1, e2 the roots of q), m = 2..6', okA3)
num = sp.factor(sp.numer(sp.together(F - x)))
check('A4 fixed points of x o x: x (1 + k x + d x^2) = 0, so besides 0 they are the zeros of the dial polynomial',
      sp.expand(num + x * (1 + k * x + d * x ** 2)) == 0 or sp.expand(num - x * (1 + k * x + d * x ** 2)) == 0, str(num))
dF = sp.diff(F, x)
ok5 = sp.simplify(dF.subs(x, 0) - 2) == 0
for root in sp.solve(1 + k * x + d * x ** 2, x):
    ok5 &= sp.simplify(dF.subs(x, root)) == 0
check('A5 multiplier 2 at 0 (a repeller), multiplier 0 at each zero of 1 + kx + dx^2 (superattracting)', ok5)
Fw = sp.cancel(F.subs(d, k ** 2 / 4))
check('A6 on the wall Delta = k^2 - 4d = 0 a zero meets a pole and x o x drops to the Mobius map 2x/(1 - kx/2)',
      sp.simplify(Fw - 2 * x / (1 - k * x / 2)) == 0, str(Fw))

# ---------------- B. Lefschetz, Nielsen, Conley-Morse on the dial plane ----------------
samples = [(1, 1), (0, 1), (2, 3), (-1, 0.5), (1, 0), (0, -1), (3, 1), (-2, -1), (2, 1), (0, 0), (4, 4), (1, 0.25)]
okB1 = True; tab = []
for (kk_, dd_) in samples:
    D = kk_ * kk_ - 4 * dd_
    for m in range(1, 7):
        deg = rp1_degree(lambda u, v: dial_power_h(u, v, kk_, dd_, m))
        want = m if D < 0 else (1 if abs(D) < 1e-12 else m % 2)
        okB1 &= deg == want
    tab.append(((kk_, dd_), D))
check('B1 circle degree of the m-fold power on RP^1: m when Delta < 0, 1 on the wall, m mod 2 when Delta > 0 (12 dial points, m = 1..6)', okB1)
okB2 = True
for (kk_, dd_) in samples:
    D = kk_ * kk_ - 4 * dd_
    if abs(D) < 1e-12: continue
    deg = rp1_degree(lambda u, v: dial_power_h(u, v, kk_, dd_, 2))
    # indices: 0 has multiplier 2 -> -1; each real zero of 1 + kx + dx^2 has multiplier 0 -> +1;
    # infinity is a fixed point only when d = 0 (there x o x is a polynomial), multiplier 0 -> +1
    nzeros = (2 if D > 0 else 0) if dd_ != 0 else (1 if kk_ != 0 else 0)
    s = -1 + nzeros + (1 if dd_ == 0 and kk_ != 0 else 0)
    okB2 &= s == 1 - deg
check('B2 Lefschetz: sum of fixed-point indices of x o x = 1 - degree at every off-wall dial point (L = 1 - 2 = -1 elliptic, 1 - 0 = 1 hyperbolic)', okB2)

def count_fixed(kk_, dd_, j, grid=2_000_003):
    th = np.linspace(0, np.pi, grid, endpoint=False) + 1e-4
    u, v = np.sin(th), np.cos(th)
    U, V = u, v
    for _ in range(j):
        U, V = dial_power_h(U, V, kk_, dd_, 2)
        nrm = np.hypot(U, V); U, V = U / nrm, V / nrm
    s = U * v - V * u
    return int(np.sum(np.sign(s[:-1]) != np.sign(s[1:])))   # s is anti-periodic, so no wrap-around pair
ell = [count_fixed(1, 1, j) for j in range(1, 9)]
hyp = [count_fixed(1, -1, j) for j in range(1, 9)]
check('B3 Nielsen: x o x has exactly 2^j - 1 points of period dividing j when Delta < 0 and exactly 3 when Delta > 0 (j = 1..8)',
      ell == [2 ** j - 1 for j in range(1, 9)] and hyp == [3] * 8, 'elliptic %s, hyperbolic %s' % (ell, hyp))
ce, ch = newton_preimage_counts([1, 1, 1], 16), newton_preimage_counts([1, 1, -1], 16)
re_ = math.log((ce[-1] - ce[-2]) / (ce[-2] - ce[-3]))
check('B4 entropy of x o x: log 2 when Delta < 0, 0 when Delta > 0 (backward orbit of the repeller)',
      abs(re_ - math.log(2)) < 1e-9 and ch[-1] == ch[3], 'rate %.6f; hyperbolic counts %s' % (re_, ch[:5]))
okB5 = True; switch = []
for dd_ in np.linspace(-1, 1, 41):
    S, R, Kd, edges = morse_graph(lambda u, v: dial_power_h(u, v, 1, dd_, 2), ncell=1500, samples=14)
    D = 1 - 4 * dd_
    if D > 1e-9:
        att = sorted(round(set_centre(edges, s), 2) for s, kd in zip(S, Kd) if kd == 'attractor')
        rep = [round(set_centre(edges, s), 2) for s, kd in zip(S, Kd) if kd == 'repeller']
        zs = np.roots([dd_, 1, 1]).real if abs(dd_) > 1e-12 else np.array([-1.0, np.inf])
        want = sorted(round(z, 2) for z in zs if np.isfinite(z))
        good = len(S) == 3 and rep == [0.0] and all(any(abs(a - w) < 0.05 * (1 + abs(w)) for a in att) for w in want)
    elif D < -1e-9:
        good = len(S) == 1 and len(S[0]) == 1500
    else:
        good = True
    okB5 &= good
    switch.append((round(dd_, 3), len(S)))
check('B5 Conley-Morse graph along k = 1: two attractors at the real zeros and a repeller at 0 for d < 1/4; one Morse set, the whole circle, for d > 1/4',
      okB5, 'Morse-set counts %s' % [s for s in switch if s[0] in (-1.0, 0.2, 0.25, 0.3, 1.0)])
c = 3.0
Fv = lambda t: (2 * t) / (1 + t * t / c ** 2)          # k = 0, d = -1/c^2: velocity addition
t = 0.1; errs = []
for _ in range(12):
    t = Fv(t); errs.append(abs(c - t))
ratios = [errs[i + 1] / errs[i] ** 2 for i in range(len(errs) - 1) if errs[i] > 1e-6 and errs[i + 1] > 0]
d1 = (Fv(c + 1e-6) - Fv(c - 1e-6)) / 2e-6
check('B6 velocity addition (k = 0, d = -1/c^2): composing a speed with itself is Newton finding the light speeds +-c; the jet at c has slope 0 (quadratic convergence)',
      abs(d1) < 1e-5 and len(ratios) >= 2 and np.std(ratios[-2:]) < 0.05 * np.mean(ratios[-2:]), 'slope at c %.1e, e_{n+1}/e_n^2 -> %.4f (= 1/2c = %.4f)' % (d1, ratios[-1], 1 / (2 * c)))
tl = rng.uniform(-3, 3); acc = 0.0; nlyap = 200000
for _ in range(nlyap):
    acc += math.log(abs(2 * (1 + tl * tl) / (1 - tl * tl) ** 2))
    tl = 2 * tl / (1 - tl * tl)
check('B7 tangent law (k = 0, d = 1): the jet along a typical orbit grows like 2^n, Lyapunov exponent log 2',
      abs(acc / nlyap - math.log(2)) < 0.01, '%.5f' % (acc / nlyap))

# ---------------- C. Real Newton maps in general ----------------
bad = tot = 0
for n in range(2, 13):
    for _ in range(300):
        r = int(rng.integers(0, n + 1))
        if (n - r) % 2: r -= 1
        r = max(r, 0)
        roots = list(rng.normal(size=r) * 2)
        for _ in range((n - r) // 2):
            a, b = rng.normal(), abs(rng.normal()) + 0.2
            roots += [a + 1j * b, a - 1j * b]
        cf = np.real(np.poly(roots)) * rng.uniform(0.5, 2) * rng.choice([-1, 1])
        nm, dn, _ = newton_parts(cf)
        tot += 1; bad += rp1_degree_alg(nm, dn) != 2 - len(distinct_real_roots(cf))
check('C1 the circle degree of real Newton\'s method is 2 - r, r = number of distinct real roots (3300 random polynomials, degree 2..12)',
      bad == 0, '%d of %d' % (tot - bad, tot))
bad = tot = 0
for _ in range(200):
    rr = list(rng.normal(size=int(rng.integers(1, 4))) * 2)
    mult = [int(rng.integers(1, 4)) for _ in rr]
    cpx = int(rng.integers(0, 2))
    roots = sum(([z] * m for z, m in zip(rr, mult)), [])
    for _ in range(cpx):
        a, b = rng.normal(), abs(rng.normal()) + 0.3
        roots += [a + 1j * b, a - 1j * b]
    # reduced Newton map, exact: q'/q = sum m_i/(y - z_i) = P/Q over the distinct roots z_i, N = y - Q/P
    zs = list(rr) + [z for z in roots if abs(np.imag(z)) > 0]
    ms = list(mult) + [1] * (2 * cpx)
    if sum(ms) < 2: continue                     # q linear: Newton is a constant map
    Q = np.poly1d(np.real(np.poly(zs)))
    P = sum((m * np.poly1d(np.real(np.poly([w for j, w in enumerate(zs) if j != i]))) for i, m in enumerate(ms)), np.poly1d([0.0]))
    nm, dn = np.poly1d([1, 0]) * P - Q, P
    tot += 1; bad += rp1_degree_alg(nm, dn) != 2 - len(rr)
check('C2 repeated roots count once: degree = 2 - (number of distinct real roots) (random polynomials of degree 2 or more, many with multiple roots)', bad == 0, '%d of %d' % (tot - bad, tot))
cases = {'y^2+1': [1, 0, 1], '(y^2+1)(y^2+4)': np.poly([1j, -1j, 2j, -2j]).real,
         'quartic, 4 real': np.poly([-2, -.5, .7, 1.9]), 'quintic, 5 real': np.poly([-2, -1, .2, 1.1, 2.3]),
         'sextic, 6 real': np.poly([-2.5, -1.4, -.3, .6, 1.5, 2.7]), 'sextic, 4 real': np.poly([-2, -.5, .7, 1.9, .1 + 1j, .1 - 1j]).real,
         'sextic, 0 real': np.poly([1j, -1j, .5 + 2j, .5 - 2j, -1 + .5j, -1 - .5j]).real}
okC3 = True; info = []
for name, cf in cases.items():
    cnt = newton_preimage_counts(cf, 14, cap=60000)
    rt = math.log((cnt[-1] - cnt[-2]) / (cnt[-2] - cnt[-3]))
    r = len(distinct_real_roots(cf))
    okC3 &= abs(rt - math.log(abs(2 - r))) < 1e-6
    info.append('%s %.4f' % (name, rt))
check('C3 the backward orbit of infinity grows at exactly log|2 - r| (the Lefschetz floor) in 7 examples', okC3, '; '.join(info))
okC4 = True; info = []
for i in range(12):
    n = int(rng.integers(4, 8))
    cf = np.poly(np.sort(rng.normal(size=n) * 2))
    tc = newton_turning_counts(cf, 12, cap=40000)
    dd2 = np.diff(tc); rate = math.log(dd2[-1] / dd2[-2])
    okC4 &= rate < math.log(n - 2) + 0.03
    info.append('n=%d %.3f/%.3f' % (n, rate, math.log(n - 2)))
check('C4 all roots real: turning points of N^j grow no faster than (n-2)^j, so h = log(n - 2), the floor (Hurley-Martin), 12 random cases', okC4, ', '.join(info[:4]) + ' ...')
tc = newton_turning_counts([1, 0, -2, 2], 18, cap=60000)
dd2 = np.diff(tc); rate = math.log(dd2[-1] / dd2[-2])
check('C5 one real root (y^3 - 2y + 2): circle degree 1 gives floor 0, yet turning points grow at about 0.6 per step: the excess over the floor comes from the folds',
      0.4 < rate < 0.8 and rp1_degree_alg(*newton_parts([1, 0, -2, 2])[:2]) == 1, 'rate %.3f' % rate)

# C6: no real roots, n = 6.  Without folds (q'' has no real root) N is a monotone degree-2 covering, so h = log 2 exactly.
# With folds the turning points grow faster than the floor (lap growth = entropy, Misiurewicz-Szlenk).
cf_a = np.real(np.poly([1 + 1j, 1 - 1j, -.5 + .8j, -.5 - .8j, .3 + 2j, .3 - 2j]))
cf_b = np.real(np.poly([2 + .3j, 2 - .3j, -1 + .5j, -1 - .5j, .1 + 1j, .1 - 1j]))
tc_a = newton_turning_counts(cf_a, 8, cap=60000)
tc_b = newton_turning_counts(cf_b, 12, cap=60000); dd6 = np.diff(tc_b); rate6 = math.log(dd6[-1] / dd6[-2])
check('C6 sextic, no real roots (degree 2, floor log 2): with no folds there are no turning points (h = log 2 exactly); with folds the turning points grow at about 0.90 > log 2',
      max(tc_a) == 0 and rp1_degree_alg(*newton_parts(cf_a)[:2]) == 2 and rp1_degree_alg(*newton_parts(cf_b)[:2]) == 2 and 0.85 < rate6 < 0.95,
      'rate %.3f vs log 2 = %.3f' % (rate6, math.log(2)))

# ---------------- D. Conley index: what knot type cannot see ----------------
r = 3.835
orbs = orbits_of_exact_period(r, 3)
pats = [pattern(o, r) for o, m in orbs]; idx = sorted(int(np.sign(1 - m)) for o, m in orbs)
check('D1 period-3 twins at r = 3.835: one pattern (knot type), fixed-point indices of f^3 +1 (stable) and -1 (unstable)',
      len(orbs) == 2 and pats[0] == pats[1] and idx == [-1, 1], '%s, multipliers %s' % (pats[0], [round(m, 3) for o, m in orbs]))
okD2 = True; born = []
for r in np.linspace(3.80, 3.86, 61):
    okD2 &= interval_index(lambda z: logistic_iter(z, r, 3), 0.44, 0.58) == 0
    born.append(len(orbits_of_exact_period(r, 3)) > 0)
check('D2 continuation: the index of the isolating interval [0.44, 0.58] stays 0 across the birth at 1 + sqrt 8 (read from the two endpoints only)',
      okD2 and not born[0] and born[-1], 'orbits present from r = %.3f' % np.linspace(3.80, 3.86, 61)[born.index(True)])
W = json.load(open('windows_index.json'))
lone = [(w['n'], round(w['r'], 4)) for w in W if not w['twin_indices']]
check('D3 all 65 windows (periods 2..9): the stable orbit has index +1; 60 have a same-pattern twin, always of index -1',
      all(w['stable_index'] == 1 for w in W) and sum(1 for w in W if w['twin_indices'] == [-1]) == 60
      and all(w['twin_indices'] in ([-1], []) for w in W), 'no twin (born by period doubling): %s' % lone)
Kst = 0.9; pairs = []; skipped = 0
for n in range(1, 14):
    for m in range(0, n):
        if gcd(m, n) != 1: continue
        pr = birkhoff_pair(Kst, m, n)
        if all(o['ok'] for o in pr): pairs.append((m, n, pr))
        else: skipped += 1
orbs = [o for _, _, pr in pairs for o in pr]
okD4 = all((o['i1'], o['i2']) == EXPECT[linear_type(o['trace'])] for o in orbs)
okD4 &= all(pr[0]['perm'] == pr[1]['perm'] for _, _, pr in pairs)
okD4 &= all(pr[0]['kind'] == 'minimizing' and pr[0]['i1'] == -1 for _, _, pr in pairs)
from collections import Counter
check('D4 standard map K = 0.9: in every Birkhoff pair the two orbits share their cyclic order (braid type), yet (index of F^n, index of F^2n) is (+1,+1) elliptic, (-1,-1) hyperbolic, (+1,-1) inverse hyperbolic',
      okD4 and len(pairs) >= 45, '%d pairs (periods 1..13, %d not converged), types %s' % (len(pairs), skipped, dict(Counter(linear_type(o['trace']) for o in orbs))))
check('D5 the index is the parity of the action\'s Morse index: i(F^n) = (-1)^(Morse + 1) for all those orbits; the minimizer (Morse 0) is always hyperbolic',
      all(o['i1'] == (-1) ** (o['morse'] + 1) for o in orbs))
ca = 0.24; sa = math.sqrt(1 - ca * ca)
import warnings; warnings.filterwarnings('ignore')
H = find_henon_orbits(ca, sa, 8, seeds=60)
okD6 = True
for o in H:
    n = o['n']; O = o['orbit']
    Dm = np.linalg.norm(O[:, None] - O[None], axis=2) + 9 * np.eye(n)
    rho = min(1e-5, Dm.min() / 10)
    i1 = winding_index(lambda z: henon_n(z, ca, sa, n), O[0], rho)
    i2 = winding_index(lambda z: henon_n(z, ca, sa, 2 * n), O[0], rho / 10)
    okD6 &= (i1, i2) == EXPECT[linear_type(o['trace'])]
check('D6 area-preserving Henon map (cos a = 0.24): the same index pair classifies every orbit found', okD6 and len(H) >= 5,
      '%d orbits, %s' % (len(H), dict(Counter(linear_type(o['trace']) for o in H))))

# ---------------- E. Shub, Yomdin, Ivanov ----------------
okE1 = True
for _ in range(60):
    kk_, dd_, a = rng.normal(size=3) * 2
    deg = rp1_degree(lambda u, v: (u + a * v + kk_ * a * u, v - dd_ * a * u))
    okE1 &= deg == int(np.sign(1 + kk_ * a + dd_ * a * a))
check('E1 one composition x -> x o a is a Mobius map, a circle homeomorphism of degree sign(1 + k a + d a^2): entropy 0 for every single element (60 random)', okE1)
check('E2 Shub bound sharp on the elliptic side: real entropy of x o x = log 2 = log|degree| = log(algebraic degree), the Gromov-Lyubich ceiling; hyperbolic side: real 0 under complex log 2',
      abs(re_ - math.log(2)) < 1e-9 and ch[-1] == ch[3])

# ---------------- F. A door to the Riemann hypothesis ----------------
J = json.load(open('jensen.json'))
check('F1 Jensen data: gamma(0) = xi(1/2) = 0.4971207781883141...', abs(mp.mpf(J['gamma'][0]) - mp.mpf(J['xi_half'])) < mp.mpf('1e-25'), J['xi_half'][:20])
rows = J['rows']
check('F2 all 248 Jensen polynomials of Xi with d <= 8, n <= 30 are hyperbolic and their real Newton maps have circle degree 2 - d',
      len(rows) == 248 and all(r['real_roots'] == r['d'] for r in rows) and all(r['newton_degree'] == 2 - r['d'] for r in rows if r['d'] >= 2))
# a polynomial that loses two real roots crosses a double-root wall, where the degree jumps by 2
jumps = []
svals = np.linspace(-1.5, 4, 220) + 1e-3          # offset so no sample sits exactly on a wall
for s in svals:
    cf = np.poly([-2, -0.3, 0.5, 1.2]) + np.array([0, 0, 0, 0, s])
    jumps.append(rp1_degree_alg(*newton_parts(cf)[:2]))
check('F3 losing real roots means crossing a double-root wall: along a quartic family the Newton degree only jumps (by 2) where two roots merge',
      set(np.abs(np.diff(jumps))) <= {0, 2} and len(set(jumps)) >= 2,
      'degrees met %s; jumps at s = %s' % (sorted(set(jumps)), [round(float(svals[i]), 3) for i in np.nonzero(np.diff(jumps))[0]]))

print('%d of %d pass' % (sum(results), len(results)))
