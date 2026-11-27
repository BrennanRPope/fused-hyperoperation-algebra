"""DRAFT. Checks for the Conley-Zehnder follow-up (section 8 of README.md). Run: python3 checks_cz.py
Reads cz.json and krein.json (from czcatalog.py and krein_test.py) and recomputes the key numbers live."""
import json, os, numpy as np
import cb3, cz
from czcatalog import reduced_with_basis
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'cz.json'))); KR = json.load(open(os.path.join(HERE, 'krein.json')))
res = []

def check(name, ok, detail=''):
    res.append(bool(ok)); print(('PASS' if ok else 'FAIL'), name, detail)

def margin(r):
    m = []
    for tr, ti in r['tau']:
        if abs(ti) > 1e-9: m.append(abs(ti))
        else: m += [abs(2 - tr), abs(2 + tr)]
    return min(m)

def L(r):
    mu = r['muhat']; return [mu[m - 1] - m * mu[0] for m in range(1, len(mu) + 1)]

# F1 convention: the harmonic oscillator H = (q^2 + p^2)/2 over less than a turn is a positive pair with index +1
phi = 1.1
M = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
p = cz.krein_pairs(M)[0]
mu1, _ = cz.crossing_index(np.eye(2), np.zeros((2, 2)), np.zeros((2, 2)), phi, 400)
check('F1 oscillator: Krein sign +, theta = phi, crossing index +1', p['s'] == 1 and abs(p['theta'] - phi) < 1e-12 and mu1 == 1,
      's=%d theta=%.3f mu=%d' % (p['s'], p['theta'], mu1))

# F2 Long's iteration formula with Krein-normalised angles against direct crossing counts (periodic linear systems)
rng = np.random.default_rng(3); T = 2 * np.pi; n_ok = n = 0
for trial in range(6):
    w = rng.uniform(0.15, 1.4, 2) * np.array([1, (-1) ** trial])
    S0 = np.diag([w[0], w[1], w[0], w[1]])
    B = rng.normal(size=(4, 4)); S1 = (B + B.T) / 2 * 0.08
    C = rng.normal(size=(4, 4)); S2 = (C + C.T) / 2 * 0.08
    mu1, M = cz.crossing_index(S0, S1, S2, T, 1200)
    pairs = cz.krein_pairs(M); r = sum(q['kind'] == 'E' for q in pairs)
    for m in (2, 3):
        mum, _ = cz.crossing_index(S0, S1, S2, m * T, 1200 * m)
        pred = m * (mu1 - r) + sum(2 * cz.E(m * q['theta'] / (2 * np.pi)) for q in pairs if q['kind'] == 'E') - r
        n += 1; n_ok += int(abs(mum - pred) < 1e-9)
check('F2 iteration formula = crossing count (stable test systems, both Krein types)', n_ok == n, '%d of %d' % (n_ok, n))

# F3 the round-one fixed-point index is the parity of the iterate index, all 692 orbits, laps 1-4
check('F3 (-1)^mu(m-th lap) = fixed-point index of F^m, m = 1..4', all(r['parity_ok'] for r in D) and len(D) == 692,
      '%d of %d orbits' % (sum(r['parity_ok'] for r in D), len(D)))

# F4 number of elliptic pairs = number of E's in the Floquet class
check('F4 elliptic pairs match the Floquet class', all(r['r'] == r['klass'].split('-').count('E') for r in D))

# F5 orbits with no elliptic pair: mu(m) = m mu(1) for every m <= 60; with one: the line bends for some m <= 60
noE = [r for r in D if r['r'] == 0]; hasE = [r for r in D if r['r'] > 0]
check('F5 straight line exactly when no elliptic pair', all(not any(L(r)) for r in noE) and all(any(L(r)) for r in hasE),
      '%d straight, %d bend' % (len(noE), len(hasE)))

# F6 the integers mu(1..60) alone give back the number of elliptic pairs
check('F6 decode from integers', all(r['r'] == r['r_decoded'] for r in D))

# F7 the 52 safe orbits with fixed-point index +1 at every lap: the 4 stable ones bend, the 48 others never do
pp = [r for r in D if margin(r) > 0.01 and r['rel_pal'] < 1e-6 and r['klass'] in ('E-E', 'H-H', 'IH-IH', 'CU')]
st = [r for r in pp if r['klass'] == 'E-E']
check('F7 separation inside the (+1,+1) group', len(pp) == 52 and len(st) == 4 and all(any(L(r)) for r in st)
      and not any(any(L(r)) for r in pp if r['klass'] != 'E-E'), '%d orbits, %d stable' % (len(pp), len(st)))

# F8 second lap: mu(2) - 2 mu(1) is 0 with no elliptic pair, odd with one, even with two
d2 = lambda r: r['muhat'][1] - 2 * r['muhat'][0]
check('F8 mu(2) - 2 mu(1) parity counts elliptic pairs', all(d2(r) % 2 == r['r'] % 2 for r in D)
      and all(d2(r) == 0 for r in noE))

# F9 the figure-eight live: mixed Krein signs, and the line first bends at lap 4
v1, v2 = 0.347116888120318, 0.532724945387253
X0 = cb3.suvakov(v1, v2); X, Phi = cb3.flow_var(X0, 6.325913982936183)
M, V = reduced_with_basis(X0, Phi); t1, t2, _ = cb3.pair_traces(M)
pairs = cz.pairs_from_traces(M, [t1, t2], V, cz.Jmat(6))
mu = [int(cz.muhat(pairs, m)) for m in range(1, 13)]
Lf = [mu[m - 1] - m * mu[0] for m in range(1, 13)]
check('F9 the eight: E-E, Krein signs (-, +), line bends first at lap 4', sorted(q['s'] for q in pairs) == [-1, 1]
      and Lf[:3] == [0, 0, 0] and Lf[3] != 0, 'angles %s, L = %s' % ([round(q['phi'], 4) for q in pairs], Lf[:6]))

# F10 confirmed stable orbits after refinement and their Krein types
ty = [x['krein'] for x in KR['orbits']]
check('F10 19 refined stable orbits: 10 mixed, 9 definite', len(ty) == 19 and ty.count('mixed') == 10,
      str(KR['fate']))

print('%d of %d pass' % (sum(res), len(res)))
if not all(res): raise SystemExit(1)
