"""Round 3 computations other than the 5,000-zero rewind. -> roads.json
T1 fused rewind on an exact toy polynomial (coefficients flow linearly; roots fold)
T2 the Kerr turn: imaginary ring average J0(aD) (forward) vs real ring average I0(aD) (rewind)
T3 barrier chain with random phases = Euler product: mean |zeta| and log|zeta| spread vs the chain
T4 far zone: S(T) at the first 5,000 zeros"""
import json, numpy as np, mpmath as mp
from numpy.polynomial import polynomial as P
from math import factorial
from scipy.special import ellipk, spence
from multiprocessing import Pool
R = {}

def evenpoly(roots):  # coefficients (low->high) of prod (z - r) for r in roots (+ mirrors), real
    c = np.array([1.0 + 0j])
    for r in roots: c = P.polymul(c, [-r, 1])
    return c.real

def series_op(c, coef):  # sum_k coef(k) D^(2k) c
    out = np.zeros_like(c); d = c.copy(); k = 0
    while d.size and np.any(d):
        out[:d.size] += coef(k) * d
        d = P.polyder(d, 2); k += 1
    return out

def heat(c, t): return series_op(c, lambda k: (-t) ** k / factorial(k))     # e^{-t D^2}: forward for t > 0
def J0op(c, a): return series_op(c, lambda k: (-1) ** k * (a / 2) ** (2 * k) / factorial(k) ** 2)
def I0op(c, a): return series_op(c, lambda k: (a / 2) ** (2 * k) / factorial(k) ** 2)

def pair_near(c, x0):
    r = np.roots(c[::-1]); r = r[np.argsort(abs(r - x0))][:2]
    return r, float((r[0] + r[1]).real), float((r[0] * r[1]).real)

# ---- T1
base = [1, 2, 2.3, 3.5, 5]; c0 = evenpoly(base + [-x for x in base])
ts = np.linspace(0.04, -0.06, 401); tr = []
for t in ts:
    ct = heat(c0, t); r, s, p = pair_near(ct, 2.15)
    tr.append({'t': float(t), 're': sorted(r.real.tolist()), 'im': float(abs(r.imag).max()), 's': s, 'p': p, 'disc': s * s - 4 * p,
               'coef': ct.tolist(), 'all_re': sorted(np.roots(ct[::-1]).real.tolist())})
# first t (going back) with a complex pair
tstar = next(x['t'] for x in tr if x['im'] > 1e-6)
cb = heat(heat(c0, -0.06), 0.06)
R['T1'] = {'roots': base, 'track': tr, 't_first_complex_grid': tstar, 'isolated_estimate': -(0.3 ** 2) / 8,
           'reversibility_max_coef_err': float(np.max(abs(cb - c0)) / np.max(abs(c0)))}
# exact collision time from the discriminant of the pair, by bisection on the polynomial's double root
def has_complex(t): return np.max(abs(np.roots(heat(c0, t)[::-1]).imag)) > 1e-9
lo, hi = -0.06, 0.0
for _ in range(60):
    m = (lo + hi) / 2
    (lo, hi) = (m, hi) if has_complex(m) else (lo, m)
R['T1']['t_collide'] = (lo + hi) / 2
print('T1', R['T1']['t_collide'], R['T1']['isolated_estimate'], R['T1']['reversibility_max_coef_err'])

# ---- T2
offA = [1, 3.5, 5]; cA = P.polymul(evenpoly(offA + [-x for x in offA]), P.polymul(evenpoly([2 + 0.3j, 2 - 0.3j]), evenpoly([-2 + 0.3j, -2 - 0.3j])))
cB = c0
avals = np.linspace(0, 0.8, 161)
def tracks(c, op, x0):
    out = []
    for a in avals:
        r, s, p = pair_near(op(c, a), x0); out.append({'a': float(a), 're': sorted(r.real.tolist()), 'im': float(abs(r.imag).max()), 'disc': s * s - 4 * p})
    return out
R['T2'] = {'A_imag_ring': tracks(cA, J0op, 2.0), 'A_real_ring': tracks(cA, I0op, 2.0),
           'B_imag_ring': tracks(cB, J0op, 2.15), 'B_real_ring': tracks(cB, I0op, 2.15)}
# check the ring average really is J0(aD): average f(z + i a cos th) over th by the trapezoid rule
z0, a0 = 1.7, 0.45; th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
avg_i = np.mean(P.polyval(z0 + 1j * a0 * np.cos(th), cA)).real; avg_r = np.mean(P.polyval(z0 + a0 * np.cos(th), cA)).real
R['T2']['ring_check'] = {'imag': [avg_i, float(P.polyval(z0, J0op(cA, a0)))], 'real': [avg_r, float(P.polyval(z0, I0op(cA, a0)))]}
# Gaussian smear along i = forward heat flow
sig = 0.37; x, w = np.polynomial.hermite_e.hermegauss(40); w = w / w.sum()
R['T2']['gauss_check'] = [float(np.sum(w * P.polyval(z0 + 1j * sig * x, cA)).real), float(P.polyval(z0, heat(cA, sig ** 2 / 2)))]
def land(tr): return next((x['a'] for x in tr if x['im'] < 1e-7), None)
def leave(tr): return next((x['a'] for x in tr if x['im'] > 1e-7), None)
R['T2']['A_lands_at'] = land(R['T2']['A_imag_ring']); R['T2']['B_leaves_at'] = leave(R['T2']['B_real_ring'])
R['T2']['A_real_lands'] = land(R['T2']['A_real_ring']); R['T2']['B_imag_leaves'] = leave(R['T2']['B_imag_ring'])
print('T2', R['T2']['ring_check'], R['T2']['gauss_check'], R['T2']['A_lands_at'], R['T2']['B_leaves_at'], R['T2']['A_real_lands'], R['T2']['B_imag_leaves'])

# ---- T3
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(float)
PR = primes(10 ** 7)
def chain_mean(sig):  # prod (2/pi) K(p^-sig) = sum d_{1/2}(n)^2 n^{-2 sig}
    return float(np.exp(np.sum(np.log(2 / np.pi * ellipk(PR ** (-2 * sig))))))
def chain_var(sig):   # variance of log|zeta| in the random-phase chain: (1/2) sum Li2(p^{-2 sig})
    x = PR ** (-2 * sig); return float(0.5 * np.sum(spence(1 - x)))
def absz(args):
    sig, t = args; return float(abs(mp.zeta(mp.mpc(sig, t))))
sigs_mean = [3.0, 2.0, 1.5, 1.25, 1.1]
rng = np.random.default_rng(1)
tsamp = rng.uniform(1000, 20000, 20000)
with Pool(4) as pool:
    meanz = {}
    for s in sigs_mean:
        v = np.array(pool.map(absz, [(s, t) for t in tsamp], chunksize=200)); meanz[s] = (float(v.mean()), float(v.std() / np.sqrt(v.size)))
        print('T3 mean', s, meanz[s], chain_mean(s), flush=True)
    sigs_dist = [1.5, 1.0, 0.75, 0.6]; dist = {}
    for s in sigs_dist:
        v = np.log(np.array(pool.map(absz, [(s, t) for t in tsamp[:12000]], chunksize=200)))
        dist[s] = v.tolist(); print('T3 var', s, v.var(), chain_var(s), flush=True)
# Monte Carlo of the chain itself (primes < 1e5)
PS = PR[PR < 1e5]
mc = {}
for s in [1.5, 1.0, 0.75, 0.6]:
    acc = np.zeros(6000)
    r = PS ** (-s)
    for k in range(0, PS.size, 2000):
        ph = rng.uniform(0, 2 * np.pi, (6000, min(2000, PS.size - k)))
        acc += -np.sum(np.log(np.abs(1 - r[k:k + 2000][None, :] * np.exp(1j * ph))), axis=1)
    mc[s] = acc.tolist()
R['T3'] = {'mean': {str(s): {'data': meanz[s], 'chain': chain_mean(s)} for s in sigs_mean},
           'var': {str(s): {'data': float(np.var(dist[s])), 'chain': chain_var(s)} for s in [1.5, 1.0, 0.75, 0.6]},
           'var_curve': {'sigma': np.linspace(0.52, 2, 75).tolist(), 'chain': [chain_var(s) for s in np.linspace(0.52, 2, 75)]},
           'dist': {str(s): dist[s] for s in dist}, 'mc': {str(s): mc[s] for s in mc}}

# ---- T4
Z = json.load(open('../elliptic-object/zeros.json'))['gamma']
mp.mp.dps = 20
S = [float(n - 1.5 - mp.siegeltheta(mp.mpf(g)) / mp.pi) for n, g in enumerate(Z, start=1)]
R['T4'] = {'gamma': [float(g) for g in Z], 'S': S}
print('T4 S range', min(S), max(S))
json.dump(R, open('roads.json', 'w'))
print('done')
