"""DRAFT checks for the three-body Conley-index test. Run: python3 checks.py  (about a minute; needs gcc, numpy).
Uses the saved floquet.json, curated.json, twins.json, analysis.json and recomputes the key numbers live."""
import json, numpy as np, sympy as sp
import cb3, section as S, degree as D

R = []
def check(name, ok):
    R.append((name, bool(ok))); print(('PASS ' if ok else 'FAIL ') + name)

FL = [o for o in json.load(open('floquet.json')) if o['ok']]
CU = [o for o in json.load(open('curated.json')) if o['ok']]
TW = {(o['name'], o['which']): o for o in json.load(open('twins.json'))}
AN = json.load(open('analysis.json'))

# ---- A. the return map ----
v1, v2 = 0.347116888120318, 0.532724945387253
X0 = cb3.suvakov(v1, v2); E8 = cb3.energy(X0); s8 = np.array([2, 0, v1, v1])
s8, res8 = S.refine(s8, E8, 1)
check('A1 the figure-eight is a fixed point of the section map (one hit per period), residual < 1e-11', res8 < 1e-11)
J8 = S.jac(s8, E8, 1)
X, Phi = cb3.flow_var(X0, 6.325913982936183)
M8 = cb3.reduced_block(X0, Phi)
ts_sec = sorted(np.real(cb3.pair_traces(J8)[:2])); ts_var = sorted(np.real(cb3.pair_traces(M8)[:2]))
check('A2 section-map Jacobian and the reduced variational block give the same Floquet traces (to 1e-6)',
      np.allclose(ts_sec, ts_var, atol=1e-6))
ev = np.linalg.eigvals(J8)
check('A3 figure-eight is E-E: all four multipliers on the unit circle, arguments ±0.0529 and ±1.8730',
      np.allclose(np.abs(ev), 1, atol=1e-6) and np.allclose(sorted(np.abs(np.angle(ev))), [0.0529, 0.0529, 1.8730, 1.8730], atol=1e-3))
check('A4 det of the section Jacobian is 1 (volume preserving) to 1e-6', abs(np.linalg.det(J8) - 1) < 1e-6)
check('A5 692 of 695 catalog orbits integrated; 108 pass the precision margin', AN['n_all'] == 692 and AN['n_robust'] == 108)

old = json.load(open('../three-body-order/stability.json')); newd = {o['cls']: o for o in FL}
dd = np.array([o['logmaxlam'] - newd[o['cls']]['logmaxlam'] for o in old if o['cls'] in newd])
check('A6 correction: the earlier finite-difference instability numbers (three-body-order/stability.json) are too high by more than 1 (in log) for 78 of 239 orbits',
      len(dd) == 239 and np.sum(dd > 1) == 78)

# ---- B. the index ----
check('B1 parity law: i(F^n) = (-1)^#H * ((-1)^(n+1))^#IH on all 108 robust orbits, n = 1..4', AN['parity_law_holds'])
t, n = sp.symbols('t n')
lam = sp.symbols('lam', positive=True)
okB2 = True
rng = np.random.default_rng(0)
def sympl_block(kind):
    if kind == 'E': a = rng.uniform(0.2, 3); return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    if kind == 'H': l = rng.uniform(1.2, 4); return np.diag([l, 1 / l])
    if kind == 'IH': l = rng.uniform(1.2, 4); return np.diag([-l, -1 / l])
for kinds in [('E', 'E'), ('H', 'H'), ('IH', 'IH')]:
    for _ in range(20):
        M = np.zeros((4, 4)); M[:2, :2] = sympl_block(kinds[0]); M[2:, 2:] = sympl_block(kinds[1])
        P = np.linalg.qr(rng.normal(size=(4, 4)))[0]; M = P @ M @ P.T
        for k in range(1, 9):
            okB2 &= np.sign(np.linalg.det(np.eye(4) - np.linalg.matrix_power(M, k))) == 1 or \
                    abs(np.linalg.det(np.eye(4) - np.linalg.matrix_power(M, k))) < 1e-8
for _ in range(20):   # complex unstable quartet
    a, m = rng.uniform(0.2, 3), rng.uniform(0.1, 1)
    l = np.exp(m + 1j * a); vals = [l, 1 / l, np.conj(l), 1 / np.conj(l)]
    for k in range(1, 9):
        okB2 &= np.real(np.prod([1 - v ** k for v in vals])) > 0
check('B2 E-E, H-H, IH-IH and complex-unstable points have index +1 for every iterate 1..8 (random examples)', okB2)
d_ok = all(D.degree(lambda x, A=A: A @ x, np.zeros(4), 1.0, n=2)[0] == [int(np.sign(np.linalg.det(A)))] * 5
           for A in [rng.normal(size=(4, 4)) for _ in range(6)])
def sq(x):
    z1 = complex(x[0], x[1]); w = z1 * z1; return np.array([w.real, w.imag, x[2], x[3]])
d2 = D.degree(sq, np.zeros(4), 1.0, n=4)[0]
check('B3 value-only degree counter: right sign on 6 random linear maps, and degree 2 for (z1^2, z2)', d_ok and d2 == [2] * 5)
good = [o for o in CU if o['residual'] < 1e-8]
settled = [o for o in good if o['deg_F'] not in (0, None) or o['degrees']['F_cube'] not in (0, None)]
agree = all(o['idx_jac'][0] == v for o in good for v in (o['deg_F'], o['degrees']['F_cube']) if v not in (0, None))
check('B4 value-only index of F^k settles on 14 of the 15 well-refined orbits (residual < 1e-8) and always matches sign det(I - DF^k)',
      len(good) == 15 and len(settled) == 14 and agree)
B8 = np.linalg.inv(J8 - np.eye(4)); B8 /= np.abs(B8).max()
deg8, _ = D.robust_degree(lambda x: S.F(x, E8, 1)[0] - x, s8, 1e-6, B=B8)
o = TW[('I.A1', 'twin')]; sH, _ = S.refine(np.array([2, 0, o['v'][0], o['v'][0]]), o['E'], 1)
JH = S.jac(sH, o['E'], 1); BH = np.linalg.inv(JH - np.eye(4)); BH /= np.abs(BH).max()
degH, _ = D.robust_degree(lambda x: S.F(x, o['E'], 1)[0] - x, sH, 1e-6, B=BH)
check('B5 live, from map values only: figure-eight index +1, its twin (Simo H orbit) index -1', deg8 == 1 and degH == -1)
pp = AN['plusplus']
check('B6 of the 52 robust orbits with index pair (+,+) only 4 are E-E (the rest H-H 26, CU 15, IH-IH 7)',
      sum(pp.values()) == 52 and pp['E-E'] == 4 and pp['H-H'] == 26 and pp['CU'] == 15 and pp['IH-IH'] == 7)

# ---- C. twins ----
a, b = TW[('I.A1', 'cat')], TW[('I.A1', 'twin')]
check('C1 figure-eight E-E (index +1), its twin E-H (index -1)', a['klass'] == 'E-E' and b['klass'] == 'E-H' and a['idx_jac'] == 1 and b['idx_jac'] == -1)
ta = min(x[0] for x in a['tau'] if abs(x[0] - 2) < 0.1); tb = min(x[0] for x in b['tau'] if abs(x[0] - 2) < 0.1)
check('C2 the fold pair: trace 1.9972 (elliptic) for the eight, 2.0084 (hyperbolic) for its twin', abs(ta - 1.9972) < 1e-3 and abs(tb - 2.0084) < 1e-3)
names = sorted({k[0] for k in TW})
opp = [nm for nm in names if TW[(nm, 'cat')]['idx_jac'] == -TW[(nm, 'twin')]['idx_jac']]
check('C3 5 of the 7 same-word pairs have opposite indices (II.C66 and II.C75 do not)', len(opp) == 5 and set(names) - set(opp) == {'II.C66', 'II.C75'})
plus_more_stable = [nm for nm in opp if (TW[(nm, 'cat')]['idx_jac'] == 1) == (TW[(nm, 'cat')]['logmaxlam'] < TW[(nm, 'twin')]['logmaxlam'])]
check('C4 in opposite-index pairs the +1 member is the less unstable one in 3 of 5 (so +1 does not mark the stabler twin)', len(plus_more_stable) == 3)

# ---- D. isolation ----
I = AN['isolation']
check('D1 every start within 1e-3 of the figure-eight stays within 50x its distance for 300 returns',
      I['figure-eight I.A1 (E-E)']['frac_stay'][:2] == [1.0, 1.0])
check('D2 every start near the H-H and complex-unstable orbits leaves, at all four radii',
      all(v == 0.0 for v in I['I.A4 (H-H)']['frac_stay'] + I['II.C24 (CU)']['frac_stay']))

# ---- E. the dial reading of one Floquet pair ----
okE1 = True; okE2 = True; okE3 = True
for o in CU:
    for tr, ti in o['tau']:
        if abs(ti) > 1e-9: continue
        Rr = (2 - tr) / 4; R2 = (2 - (tr * tr - 2)) / 4
        okE1 &= abs(R2 - 4 * Rr * (1 - Rr)) < 1e-9 * max(1, abs(R2))
def tangent_law(x, y, d): return (x + y) / (1 - d * x * y)
for th in rng.uniform(0.1, 3.0, 10):
    x = np.tan(th / 2); y = x
    for k in range(2, 6):
        y = tangent_law(y, x, 1.0)
        okE2 &= abs(np.arctan(y) * 2 % (2 * np.pi) - (k * th) % (2 * np.pi)) < 1e-9 or abs(np.cos(2 * np.arctan(y)) - np.cos(k * th)) < 1e-9
for mu in rng.uniform(0.1, 1.5, 10):
    x = np.tanh(mu / 2); y = x
    for k in range(2, 6):
        y = tangent_law(y, x, -1.0)
        okE2 &= abs(np.arctanh(y) * 2 - k * mu) < 1e-9
    # index factor 2 - tau = 4R; hyperbolic: R = -x^2/(1-x^2) < 0 ; elliptic: R = x^2/(1+x^2) > 0
    okE3 &= ((2 - 2 * np.cosh(mu)) / 4 + x * x / (1 - x * x)) < 1e-12
for th in rng.uniform(0.1, 3.0, 10):
    x = np.tan(th / 2); okE3 &= abs((2 - 2 * np.cos(th)) / 4 - x * x / (1 + x * x)) < 1e-12
check("E1 Greene's residue of each pair squares by the logistic rule R -> 4R(1 - R) (curated orbits)", okE1)
check('E2 iterating an elliptic pair is the tangent law (d = +1) on x = tan(theta/2); a hyperbolic pair is velocity addition (d = -1) on x = tanh(mu/2)', okE2)
check('E3 a pair\'s index factor sign(2 - trace) = sign of its curvature dial d: + for elliptic, - for hyperbolic', okE3)

print('\n%d of %d pass' % (sum(ok for _, ok in R), len(R)))
