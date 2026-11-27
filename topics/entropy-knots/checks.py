"""Checks for topics/entropy-knots (DRAFT, Nov 27 2026).
Run: python3 checks.py  (needs numpy, scipy, sympy, mpmath; about 15 minutes).
Recomputes windows.json, zeros.json and posbraid.json, then tests them."""
import json, math, subprocess, sys
import numpy as np, mpmath as mp, sympy as sp

subprocess.run([sys.executable, 'windows.py'], check=True, capture_output=True)
subprocess.run([sys.executable, 'window_braids.py'], check=True, capture_output=True)
subprocess.run([sys.executable, 'zeros.py'], check=True, capture_output=True)
subprocess.run([sys.executable, 'posbraid.py'], check=True, capture_output=True)

from windows import F, pattern, pattern_matrix, pattern_entropy, kneading_exact
from braids import braid_entropy
R = json.load(open('windows.json')); Z = json.load(open('zeros.json')); PB = json.load(open('posbraid.json'))
phi = (1 + 5 ** .5) / 2
results = []
def check(name, ok, info=''):
    results.append(bool(ok)); print(('PASS' if ok else 'FAIL'), name, info)

# ---- A. One dimension: the orbit's pattern pins its stable window ----
from collections import Counter
cnt = Counter(r['n'] for r in R)
check('A1 superstable windows of period 2..9 number 1,1,2,3,5,9,16,28 (65)',
      [cnt[n] for n in range(2, 10)] == [1, 1, 2, 3, 5, 9, 16, 28], str(dict(cnt)))
dev = max(max(abs(r['h_pattern'] - r['h_map_below']), abs(r['h_pattern'] - r['h_map_above'])) for r in R)
check('A2 entropy forced by the orbit pattern = map entropy at its window (all 65)', dev < 1e-9, 'max dev %.1e' % dev)
# e^{-h} is a zero of the kneading numerator and of the Markov characteristic polynomial
okz = True
for r in R:
    if r['h_pattern'] < 1e-9: continue
    z = math.exp(-r['h_pattern'])
    kn = sum(c * z ** k for k, c in enumerate(r['kneading']))
    A = pattern_matrix(r['perm']).astype(float)
    ch = np.linalg.det(np.eye(len(A)) - z * A)
    okz &= abs(kn) < 1e-8 and abs(ch) < 1e-8
check('A3 e^-h is a zero of the kneading determinant and of det(I - tA) for every window', okz)
p3 = [r for r in R if abs(r['h_pattern'] - math.log(phi)) < 1e-9]
lo, hi = 1 + math.sqrt(8), 3.8568
check('A4 the log(phi) plateau holds the period-3 window and its period-6 and period-9 sub-windows, all inside [1+sqrt8, 3.8568]',
      sorted(r['n'] for r in p3) == [3, 6, 9] and all(lo < r['r'] < hi for r in p3), str([(r['n'], round(r['r'], 5)) for r in p3]))
# saddle-node partner: both period-3 orbits at r = 3.835 have the same pattern
mp.mp.dps = 50
rr = mp.mpf('3.835')
def f3(x): return F(rr, F(rr, F(rr, x)))
roots = []
for x0 in np.linspace(0.01, 0.99, 400):
    try:
        x = mp.findroot(lambda y: f3(y) - y, mp.mpf(x0))
    except Exception: continue
    if 0 < x < 1 and abs(F(rr, x) - x) > 1e-6 and all(abs(x - y) > 1e-12 for y in roots): roots.append(x)
orbs = []
for x in roots:
    o = sorted([x, F(rr, x), F(rr, F(rr, x))])
    if not any(abs(o[0] - q[0]) < 1e-12 for q in orbs): orbs.append(o)
mults, pats = [], []
for o in orbs:
    m_ = 1
    for y in o: m_ *= rr * (1 - 2 * y)
    mults.append(float(m_))
    x = o[0]; pts = [x, F(rr, x), F(rr, F(rr, x))]
    order = sorted(range(3), key=lambda i: pts[i]); rank = {i: k for k, i in enumerate(order)}
    perm = [0] * 3
    for i in range(3): perm[rank[i]] = rank[(i + 1) % 3]
    pats.append(perm)
check('A5 at r = 3.835 the stable and unstable period-3 orbits have the same pattern (same 1D braid type)',
      len(orbs) == 2 and pats[0] == pats[1] and min(abs(m) for m in mults) < 1 < max(abs(m) for m in mults),
      'patterns %s, multipliers %s' % (pats, [round(m, 3) for m in mults]))

# ---- B. The plane: the orbit's braid ----
check('B1 braid entropy test: sigma1 sigma2^-1 gives log(phi^2)', abs(braid_entropy([1, -2], 3) - 2 * math.log(phi)) < 1e-3)
check('B2 braid entropy test: sigma1^2 sigma2^-1 gives log(2+sqrt3)', abs(braid_entropy([1, 1, -2], 3) - math.log(2 + 3 ** .5)) < 1e-3)
check('B3 braid entropy test: sigma1 sigma2 (a rotation) gives 0', braid_entropy([1, 2], 3) < 1e-3)
check('B4 braid entropy <= pattern entropy for all 65 windows (Thurston-Nielsen floor below the interval floor)',
      all(r['h_braid'] <= r['h_pattern'] + 5e-3 for r in R))
eq = [r for r in R if r['h_pattern'] > 1e-9 and abs(r['h_braid'] - r['h_pattern']) < 2e-3]
zero = [r for r in R if r['h_pattern'] > 1e-9 and r['h_braid'] < 1e-2]
mid = [r for r in R if r['h_pattern'] > 1e-9 and r not in eq and r not in zero]
check('B5 of the 62 positive-entropy windows: braid keeps all the entropy in 8, part in 31, none in 23',
      (len(eq), len(mid), len(zero)) == (8, 31, 23), '%d / %d / %d' % (len(eq), len(mid), len(zero)))
r3 = [r for r in R if r['n'] == 3][0]
check('B6 the period-3 window orbit braids as a rotation: braid entropy 0 (period 3 forces chaos on a line, not in the plane)',
      r3['h_braid'] < 1e-3, 'braid %s' % r3['braid'])
rot = [r for r in R if r['perm'] == list(range(1, r['n'])) + [0]]
check('B7 every rotation pattern (1,2,...,n-1,0) braids with entropy 0', all(r['h_braid'] < 1e-2 for r in rot), '%d windows' % len(rot))
r5 = [r for r in R if r['perm'] == [1, 3, 4, 2, 0]][0]
lam5 = max(abs(np.roots([1, -1, -1, -1, 1])))
check('B8 the period-5 window [1,3,4,2,0] braid has the smallest 5-strand dilatation 1.72208 (Ham-Song 2007)',
      abs(r5['h_braid'] - math.log(lam5)) < 2e-3, 'h=%.5f log(lam)=%.5f' % (r5['h_braid'], math.log(lam5)))
r7 = [r for r in R if r['perm'] == [1, 3, 4, 5, 6, 2, 0]][0]
lam7 = max(abs(np.roots([1, -1, 0, -1])))
check('B9 the period-7 window [1,3,4,5,6,2,0] braid reaches the smallest 7-strand dilatation 1.46557 (Lanneau-Thiffeault 2011)',
      abs(r7['h_braid'] - math.log(lam7)) < 2e-3, 'h=%.5f log(lam)=%.5f' % (r7['h_braid'], math.log(lam7)))
# Not a locator: braid entropy is not monotone along the parameter
hb = [r['h_braid'] for r in R]
drops = sum(1 for a, b in zip(hb, hb[1:]) if b < a - 0.05)
check('B10 along r, braid entropy jumps up and down (%d drops > 0.05) while map entropy only rises' % drops,
      drops > 10 and all(b >= a - 1e-9 for a, b in zip([r['h_pattern'] for r in R], [r['h_pattern'] for r in R][1:])))

# ---- C. Zeros: trace monoids, hard-core gases and the Jones index ----
check('C1 path P_m (Temperley-Lieb generators): d^(m+1) mu(1/d^2) = U_(m+1)(d/2), m = 1..8',
      all(Z['P%d' % m]['identity'] in (True, 'True') for m in range(1, 9)))
check('C2 path P_m: growth = Beraha number 4cos^2(pi/(m+2)), all zeros real, m = 1..8',
      all(abs(Z['P%d' % m]['growth'] - Z['P%d' % m]['beraha']) < 1e-9 for m in range(1, 9)))
check('C3 cycle C_m: growth = 4cos^2(pi/(2m)), m = 3..8',
      all(abs(Z['C%d' % m]['growth'] - Z['C%d' % m]['beraha']) < 1e-9 for m in range(3, 9)))
mp_ = {tuple(p['sizes']): p for p in Z['multipartite']}
check('C4 affine centre patterns with growth below 4: [1,1],[2,1],[1,1,1],[2,2],[2,1,1] are Beraha B4,B5,B6,B8,B12',
      all(mp_[s]['is_beraha'] for s in [(1, 1), (2, 1), (1, 1, 1), (2, 2), (2, 1, 1)]) and
      [mp_[s]['beraha_n'] for s in [(1, 1), (2, 1), (1, 1, 1), (2, 2), (2, 1, 1)]] == [4, 5, 6, 8, 12])
check('C5 but [3,1] (the claw) and [4,1] are not Beraha and have complex zeros (Chudnovsky-Seymour: claw-free => real)',
      not mp_[(3, 1)]['is_beraha'] and not mp_[(4, 1)]['is_beraha'] and not mp_[(3, 1)]['real_zeros'] and not mp_[(4, 1)]['real_zeros'],
      'growth %.5f, %.5f' % (mp_[(3, 1)]['growth'], mp_[(4, 1)]['growth']))
sm = Z['strand_maps']
check('C6 strand maps (affine maps on neighbouring coordinate pairs) give exactly the trace counts of the path, m = 2..6',
      all(sm[k]['match'] in (True, 'True') for k in sm), str({k: sm[k]['composites'][:5] for k in sm}))
check('C7 positive braid words: Burau products at t = 3/7 match Deligne counts for B3+, B4+, B5+',
      all(PB[k]['match'] for k in PB))
check('C8 positive 3-braids grow like phi (entropy log phi), positive 4-braids like 2.08679',
      abs(PB['B3']['growth'] - phi) < 1e-9 and abs(PB['B4']['growth'] - 2.086791) < 1e-5)
# the braid relation lowers the growth below the trace monoid on the same generators
check('C9 braid relation lowers growth: B3+ 1.618 < 2 (free pair), B4+ 2.087 < 2.618 (path P3), B5+ 2.395 < 3 (path P4)',
      PB['B3']['growth'] < 2 and PB['B4']['growth'] < Z['P3']['growth'] and PB['B5']['growth'] < Z['P4']['growth'])
# Graeffe: mu_P3(t^2) = D(t) D(-t), D = period-3 kneading numerator 1 - t - t^2
tt = sp.symbols('t')
check('C10 mu_P3(t^2) = (1 - t - t^2)(1 + t - t^2): the trace monoid of three maps (two sharing a centre) is the period-3 window squared',
      sp.expand((1 - tt - tt**2) * (1 + tt - tt**2) - (1 - 3 * tt**2 + tt**4)) == 0 and r3['kneading'] == [1, -1, -1])

print('%d of %d pass' % (sum(results), len(results)))
