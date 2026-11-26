"""DRAFT. Fits on hyper.json: is the period, or its hyperspherical versions, simpler in word length or in runs?"""
import json, os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from hyper import rf, runs
Z = json.load(open(os.path.join(H, 'hyper.json'))); C = [c for c in Z['catalog'] if c['ok']]; O = Z['ours']
def q(c):
    E = abs(c['E']); n = len(c['word'])
    return dict(Ts=c['T'] * E**1.5, taus=c['tau'] * E**-0.5, arc=c['arc'], nR=c['nRmax'], Is=c['Imean'] * E**2, n=n, nruns=len(runs(c['word'])), rf=rf(c['word']))
QC = [q(c) for c in C]; QO = [q(c) for c in O]
def loo(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float); A = X.T @ X; b = X.T @ y; r = []
    for i in range(len(y)):
        Ai = A - np.outer(X[i], X[i]); bi = b - X[i] * y[i]; r.append((X[i] @ np.linalg.solve(Ai, bi) - y[i]) / y[i])
    r = np.array(r); return dict(sd=float(np.sqrt(np.mean(r**2))), max=float(np.max(np.abs(r))), coef=np.linalg.lstsq(X, y, rcond=None)[0].tolist())
out = {}
for name in ('Ts', 'taus', 'arc', 'Is'):
    y = [x[name] for x in QC]
    out[name] = {'letters': loo([[x['n']] for x in QC], y), 'run prices': loo([x['rf'] for x in QC], y), 'number of runs': loo([[x['nruns']] for x in QC], y)}
nR = np.array([x['nR'] for x in QC]); n = np.array([x['n'] for x in QC]); nr = np.array([x['nruns'] for x in QC])
out['R breaths'] = dict(equal_letters=float(np.mean(nR == n)), equal_runs=float(np.mean(nR == nr)), equal_half_letters=float(np.mean(nR == n // 2)),
    ratio_letters_mean=float(np.mean(nR / n)), ratio_letters_sd=float(np.std(nR / n)), ratio_runs_mean=float(np.mean(nR / nr)), ratio_runs_sd=float(np.std(nR / nr)),
    corr_T_nR=float(np.corrcoef([x['Ts'] for x in QC], nR)[0, 1]), Ts_per_breath=loo([[v] for v in nR], [x['Ts'] for x in QC]))
# stability (66 polished): log of the largest Floquet multiplier per unit of each clock
lam = np.array([max(o['maxlam'], 1.0) for o in O]); L = np.log(lam)
per = {}
for name, v in [('per unit T*', [x['Ts'] for x in QO]), ('per unit tau*', [x['taus'] for x in QO]), ('per letter', [x['n'] for x in QO]),
                ('per run', [x['nruns'] for x in QO]), ('per breath of R', [x['nR'] for x in QO]), ('per unit shape-sphere length', [x['arc'] for x in QO])]:
    r = L / np.array(v, float); m = r[L > 0.01]; per[name] = dict(mean=float(np.mean(m)), cv=float(np.std(m) / np.mean(m)), n=int(len(m)))
out['stability'] = dict(per=per, n_stable=int(np.sum(L <= 0.01)), stable=[(O[i]['word'], O[i].get('match'), float(lam[i])) for i in np.nonzero(L <= 0.01)[0]])
yo = [x['Ts'] for x in QO]
coef = np.array(out['Ts']['run prices']['coef']); ro = (np.array([x['rf'] for x in QO]) @ coef - yo) / np.array(yo)
out['ours run law'] = dict(sd=float(np.sqrt(np.mean(ro**2))), within1=float(np.mean(abs(ro) < .01)), n=len(yo))
json.dump(out, open(os.path.join(H, 'hyper_fit.json'), 'w'), indent=1)
for k in ('Ts', 'taus', 'arc', 'Is'): print(k, {m: '%.4f/%.3f' % (v['sd'], v['max']) for m, v in out[k].items()})
print(json.dumps(out['R breaths'], indent=0)); print(json.dumps(out['stability'], indent=0)[:1500]); print(out['ours run law'])
