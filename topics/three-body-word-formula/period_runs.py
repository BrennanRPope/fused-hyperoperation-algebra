"""DRAFT. Run law for the period. Split the cyclic word into runs: maximal stretches whose letters all have the same
case (the shape keeps going round the collisions the same way); a case change is a reversal.
Models, each fitted on 4/5 of Li & Liao's 695 orbits and tested on the held-out 1/5 (5 folds):
  letters only:         T* = c L
  letters + reversals:  T* = b L + r R        (R = number of runs = number of reversals)
  cost per run length:  T* = sum_k c_k n_k    (n_k = runs of length k)
Outside test: our own one-lap orbits from ../three-body (read only)."""
import json, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
def runs(w):
    n = len(w); i0 = next((i for i in range(n) if w[i].isupper() != w[i - 1].isupper()), 0); u = w[i0:] + w[:i0]
    out, cur = [], 1
    for x, y in zip(u, u[1:]):
        if x.isupper() == y.isupper(): cur += 1
        else: out.append(cur); cur = 1
    out.append(cur); return out
KS = [1, 2, 3, 4, 5]
def feats(w):
    r = runs(w); return dict(L=len(w), R=len(r), n=[sum(1 for x in r if x == k) for k in KS], syz=2 * len(w) - len(r))
F = [feats(o['word']) for o in d]; Ts = np.array([o['Tstar'] for o in d])
models = {'letters only': lambda f: [f['L']], 'letters + reversals': lambda f: [f['L'], f['R']],
          'syzygies (incl. hidden 3-middle)': lambda f: [f['syz']], 'cost per run length': lambda f: f['n']}
rng = np.random.default_rng(5); idx = rng.permutation(len(d)); folds = np.array_split(idx, 5)
out = {'run_length_counts': {k: int(sum(f['n'][i] for f in F)) for i, k in enumerate(KS)}}
for nm, fn in models.items():
    X = np.array([fn(f) for f in F], float); pred = np.zeros_like(Ts)
    for k in range(5):
        te = folds[k]; tr = np.setdiff1d(idx, te); w = np.linalg.lstsq(X[tr], Ts[tr], rcond=None)[0]; pred[te] = X[te] @ w
    rel = (Ts - pred) / Ts; w = np.linalg.lstsq(X, Ts, rcond=None)[0]
    out[nm] = dict(coeffs=w.tolist(), heldout_rel_sd=float(np.std(rel)), heldout_max_rel=float(np.max(np.abs(rel))),
                   within_1pct=float(np.mean(np.abs(rel) < 0.01)), within_03pct=float(np.mean(np.abs(rel) < 0.003)), rel=rel.tolist())
# outside test on our own orbits (one lap only)
E = json.load(open(os.path.join(H, '..', 'three-body', 'results.json')))['E']
wr = np.array(out['cost per run length']['coeffs']); wlr = np.array(out['letters + reversals']['coeffs'])
ours = [(e['word'], e['Tstar'], e['match']) for e in E if abs(e['Tstar'] / e['len'] - 2.43) < 0.3]
rr = np.array([(T - np.array(feats(w)['n']) @ wr) / T for w, T, m in ours])
out['our_orbits'] = dict(n=len(ours), rel=rr.tolist(), sd=float(np.std(rr)), within_1pct=float(np.mean(np.abs(rr) < 0.01)),
                         listed=[m for w, T, m in ours], words=[w for w, T, m in ours])
json.dump(out, open(os.path.join(H, 'period_runs.json'), 'w'))
for nm in models: print(nm, {k: (np.round(v, 4) if k != 'rel' else '') for k, v in out[nm].items()})
print('ours', out['our_orbits']['sd'], out['our_orbits']['within_1pct'])
