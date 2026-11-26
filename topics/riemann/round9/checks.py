"""Round 9 checks."""
import json, re
import numpy as np
import mpmath as mp
from math import pi, log
from stringlib import theta, gram, String, true_S, primes_upto, mean_spacing

R = []
def check(name, ok, info=''):
    R.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {info}")

gam = np.array(json.load(open('zeros.json'))['gamma'])
P = json.load(open('probe.json')); L = json.load(open('learn.json')); G = json.load(open('gram.json'))
D = json.load(open('drift.json')); PR = json.load(open('predict.json')); RW = json.load(open('rwidth.json'))

z1000 = np.load('../../primes/zeta_zeros_1000.npy')
check('zeros 1-1000 match the earlier table', np.max(np.abs(gam[:1000] - z1000[:1000])) < 1e-9)
check('5500 zeros, strictly increasing', len(gam) == 5500 and np.all(np.diff(gam) > 0))
for n in [1, 2500, 5500]:
    check(f'zero {n} agrees with mpmath.zetazero', abs(float(mp.zetazero(n).imag) - gam[n - 1]) < 1e-8)
for t in [20.0, 300.0, 5000.0]:
    check(f'theta asymptotic = siegeltheta at t={t}', abs(theta(t) - float(mp.siegeltheta(t))) < 1e-9)
g = gram(np.array([-1, 0, 125, 126]))
check('theta(g_m) = m pi', np.max(np.abs(theta(g) - np.array([-1, 0, 125, 126]) * pi)) < 1e-9)
check('g_-1 = 9.6669, g_0 = 17.8456', abs(g[0] - 9.66690808) < 1e-6 and abs(g[1] - 17.84559954) < 1e-6)
S_g = np.array(G['S_gram']); counts = np.array(G['counts']); ms = np.array(G['ms'])
check('S at Gram points is an integer', np.all(S_g == np.round(S_g)))
check('zeros in [g_m, g_m+1) = 1 + jump of S', np.all(counts == 1 + np.diff(S_g)))
i125 = int(np.nonzero(ms == 125)[0][0])
check("Hutchinson's first failure: [g125,g126) empty, [g126,g127) has two", counts[i125] == 0 and counts[i125 + 1] == 2)
check('no Rosser failure below zero 5500 (Brent: first at Gram point 13,999,525)', G['rosser_failures'] == 0)
# good iff S even, against the sign of Z at Gram points
bad = ms[S_g % 2 != 0][:20]; good = ms[S_g % 2 == 0][::400][:20]
ok = True
for m in np.concatenate([bad, good]):
    z = float(mp.siegelz(float(gram(m))))
    is_good = (-1) ** int(m) * z > 0
    ok &= (is_good == (S_g[ms == m][0] % 2 == 0))
check('good Gram point iff S(g_m) even (40 points, against siegelz)', ok)
hist = {int(k): v for k, v in G['hist'].items()}
check('interval counts 393 / 4715 / 387 / 3', hist == {0: 393, 1: 4715, 2: 387, 3: 3})
st = String(X=3000.0, taper='sharp'); E = st.levels(10)
_e = np.abs(E - gam[:10]) / mean_spacing(gam[:10])
check('prime string (X=3000) places zeros 1-10 within 0.05 mean spacings', np.max(_e) < 0.05, f"max err {np.max(_e):.4f} spacings")
for taper in ['sharp', 'fejer', 'gauss']:
    ks = sorted([k for k in P['cutoff_grid'] if k.startswith(taper)], key=lambda k: float(k.split('|')[1]))
    v = [P['cutoff_grid'][k]['test']['rms'] for k in ks]
    check(f'held-out error falls with the cutoff ({taper})', all(a > b for a, b in zip(v, v[1:])))
r = np.array(L['free_weights']['fitted'][:10]) / np.array(L['free_weights']['truth'][:10])
check('weights learned from 500 zeros within 5% of p^(-m/2)/m (first 10)', np.max(np.abs(r - 1)) < 0.05, str(np.round(r, 3)))
check('Landau dips at log 2 (all zeros) and log 3, log 5 (first 500)',
      min(abs(np.array(L['landau_minima_all']) - log(2))) < 0.01 and
      min(abs(np.array(L['landau_minima_train']) - log(3))) < 0.01 and min(abs(np.array(L['landau_minima_train']) - log(5))) < 0.01)
tr = [t['train']['rms'] for t in L['trap']]; te = [t['test']['rms'] for t in L['trap']]
check('trap: fit error falls and held-out error rises past 320 settings', tr[-1] < tr[4] / 5 and te[-1] > 2 * te[4])
check('orbit lengths must be log p to about 3e-5 at these heights', RW['0.0001']['1501-5500'] > 3 * RW['0']['1501-5500'] and RW['3e-05']['1501-5500'] < 1.3 * RW['0']['1501-5500'])
sh = [r_ for r_ in D if r_['taper'] == 'sharp']
dev = [abs(r_['s']['best_train'] - 0.5) for r_ in sh]
check('best line drifts toward 1/2 as orbits are added (sharp)', all(a >= b for a, b in zip(dev, dev[1:])) and dev[-1] < dev[0])
gs = [r_ for r_ in D if r_['taper'] == 'gauss']
check('gauss taper overshoots below 1/2 and comes back', all(r_['s']['best_train'] < 0.5 for r_ in gs) and gs[-1]['s']['best_train'] >= gs[0]['s']['best_train'])
def _rms(lab):
    pr_ = np.array(PR[lab], float); idx = np.arange(500, 5500); idx = idx[idx < len(pr_)]
    d = (pr_[idx] - gam[idx]) / mean_spacing(gam[idx]); d = d[np.isfinite(d)]
    return float(np.sqrt(np.mean(d * d)))
check('prime string beats the Gram midpoint on held-out zeros', _rms('P') < _rms('G') / 4, f"P {_rms('P'):.4f} G {_rms('G'):.4f}")
check('Rosser block hit is 1.0 on held-out zeros', PR['rosser_block_hit'] == 1.0)
check('p = 2 roughly doubles the level density', 1.6 < P['dial_p2']['density_ratio'] / (np.log(3000 / (2 * pi)) / pi / 2) < 2.4)
pr = primes_upto(100000).astype(float)
fermi = np.prod(1 + pr ** -2.0)
check('Fermi primon gas: prod(1+p^-2) = zeta(2)/zeta(4)', abs(fermi - float(mp.zeta(2) / mp.zeta(4))) < 1e-5)
st_ = P['string_theory']
check('each p-adic place matters (dropping p=2 hurts more than p=97)', st_['drop_place']['2']['test']['rms'] > st_['drop_place']['97']['test']['rms'] > P['cutoff_grid']['gauss|200']['test']['rms'])
check('Bose beats Fermi and Boltzmann', st_['statistics']['Bose (k=1)'] < st_['statistics']['Boltzmann (k=0)'] < st_['statistics']['Fermi (k=-1)'])
check('sign flip and time reversal both hurt', st_['sign_flip'] > 5 * P['cutoff_grid']['gauss|200']['test']['rms'] and st_['time_reversal'] > 2 * P['cutoff_grid']['gauss|200']['test']['rms'])
lean = open('lean/RiemannRound9.lean').read()
check('Lean: 12 theorems, no sorry', len(re.findall(r'^theorem ', lean, re.M)) == 12 and 'sorry' not in lean)
print(f"\n{sum(o for _, o in R)} of {len(R)} pass")
