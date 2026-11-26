"""DRAFT. Pair-type period law: T* = t_aa * n(aa) + t_ab * n(ab) + t_aB * n(aB), counting the cyclic word's adjacent
letter pairs by type (same letter / same direction / opposite direction). Fit on 4/5, test on the held-out 1/5."""
import json, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
def pairs(w):
    c = w + w[0]; n_same_letter = n_same_dir = n_opp = 0
    for x, y in zip(c, c[1:]):
        if x == y: n_same_letter += 1
        elif x.isupper() == y.isupper(): n_same_dir += 1
        else: n_opp += 1
    return [n_same_letter, n_same_dir, n_opp]
X = np.array([pairs(o['word']) for o in d], float); Ts = np.array([o['Tstar'] for o in d]); Lf = X.sum(1)
rng = np.random.default_rng(11); idx = rng.permutation(len(d)); folds = np.array_split(idx, 5)
out = {}
for nm, cols in [('length only', None), ('two types (same direction / opposite)', 'two'), ('three types', 'three')]:
    pred = np.zeros_like(Ts)
    def M(sel):
        if cols is None: return Lf[sel][:, None]
        if cols == 'two': return np.column_stack([X[sel, 0] + X[sel, 1], X[sel, 2]])
        return X[sel]
    for k in range(5):
        te = folds[k]; tr = np.setdiff1d(idx, te)
        w = np.linalg.lstsq(M(tr), Ts[tr], rcond=None)[0]; pred[te] = M(te) @ w
    rel = (Ts - pred) / Ts
    wall = np.linalg.lstsq(M(np.arange(len(d))), Ts, rcond=None)[0]
    out[nm] = dict(coeffs=wall.tolist(), heldout_rel_sd=float(np.std(rel)), heldout_max_rel=float(np.max(np.abs(rel))),
                   within_1pct=float(np.mean(np.abs(rel) < 0.01)), within_half_pct=float(np.mean(np.abs(rel) < 0.005)))
# our own 66 as an outside test (one-lap words only), using the three-type law fitted on Li-Liao
E = json.load(open(os.path.join(H, '..', 'three-body', 'results.json')))['E']
w3 = np.array(out['three types']['coeffs'])
ours = [(e['word'], e['Tstar']) for e in E if abs(e['Tstar'] / e['len'] - 2.43) < 0.3]
rel = [(T - np.array(pairs(w)) @ w3) / T for w, T in ours]
out['our_orbits_one_lap'] = dict(n=len(ours), rel_sd=float(np.std(rel)), max_rel=float(np.max(np.abs(rel))), within_1pct=float(np.mean(np.abs(np.array(rel)) < 0.01)))
json.dump(out, open(os.path.join(H, 'period_pairs.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
