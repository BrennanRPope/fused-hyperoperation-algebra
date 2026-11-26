"""DRAFT. Check the topological Kepler third law (Dmitrasinovic & Suvakov 2015; Li & Liao 2017) on 695 + our 66 orbits,
and test simple word features as predictors of the start (v1, v2) with honest 5-fold cross-validation."""
import json, os, re, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
E = json.load(open(os.path.join(H, '..', 'three-body', 'results.json')))['E']
out = {}
Lf = np.array([o['Lf'] for o in d], float); Ts = np.array([o['Tstar'] for o in d])
cls = np.array([o['cls'].split('.')[0] + '.' + o['cls'].split('.')[1][0] for o in d])
v1 = np.array([o['v1'] for o in d]); v2 = np.array([o['v2'] for o in d])
En = -2.5 + 3 * (v1**2 + v2**2)
# recompute T* from T and energy as a check
Ts2 = np.array([o['T'] for o in d]) * np.abs(En)**1.5
out['Tstar_recompute_maxrel'] = float(np.max(np.abs(Ts2 - Ts) / Ts))
A = np.column_stack([Lf, np.ones_like(Lf)]); (a, b), *_ = np.linalg.lstsq(A, Ts, rcond=None)
r = Ts - (a * Lf + b)
out['fit'] = dict(slope=a, intercept=b, rms=float(np.sqrt(np.mean(r**2))), max_abs_rel=float(np.max(np.abs(r) / Ts)),
                  ratio_mean=float(np.mean(Ts / Lf)), ratio_sd=float(np.std(Ts / Lf)), ratio_min=float(np.min(Ts / Lf)), ratio_max=float(np.max(Ts / Lf)))
# through origin
a0 = float((Lf @ Ts) / (Lf @ Lf)); r0 = Ts - a0 * Lf
out['fit_origin'] = dict(slope=a0, rms=float(np.sqrt(np.mean(r0**2))), frac_within_10pct=float(np.mean(np.abs(r0) / Ts < 0.10)), frac_within_5pct=float(np.mean(np.abs(r0) / Ts < 0.05)))
out['by_class'] = {c: dict(n=int((cls == c).sum()), ratio_mean=float(np.mean(Ts[cls == c] / Lf[cls == c])), ratio_sd=float(np.std(Ts[cls == c] / Lf[cls == c]))) for c in sorted(set(cls))}
# does the leftover depend on where the orbit starts? correlation of T*/Lf with |v| and angle
ratio = Ts / Lf; vm = np.hypot(v1, v2); ang = np.arctan2(v2, v1)
out['corr_ratio_vs_speed'] = float(np.corrcoef(ratio, vm)[0, 1]); out['corr_ratio_vs_angle'] = float(np.corrcoef(ratio, ang)[0, 1])
# our 66 (word length from our reader; power = laps)
ours = [(e['len'], e['Tstar'], e['power']) for e in E]
oL = np.array([x[0] for x in ours], float); oT = np.array([x[1] for x in ours])
out['ours_ratio_mean'] = float(np.mean(oT / oL)); out['ours_ratio_sd'] = float(np.std(oT / oL)); out['ours_n'] = len(ours)

# ---------- word features -> start, 5-fold CV
def feats(w):
    n = len(w); f = [n, w.count('a') + w.count('A'), w.count('b') + w.count('B'), w.count('a') + w.count('b'),
         sum(w[i] == w[i + 1] for i in range(n - 1)),             # doubled letters (loops round one collision twice)
         len(re.findall(r'(?=abAB|baBA|ABab|BAba)', w)),         # figure-eight syllables
         w.count('a') - w.count('A'), w.count('b') - w.count('B')]  # exponent sums (abelianized word)
    return f
X = np.array([feats(o['word']) for o in d], float)
X = np.column_stack([X, X[:, 1:] / X[:, :1], (cls == 'I.A'), (cls == 'I.B'), (cls == 'II.C')])
Y = np.column_stack([v1, v2])
rng = np.random.default_rng(1); idx = rng.permutation(len(d)); folds = np.array_split(idx, 5)
def cv(pred):
    err = np.zeros_like(Y)
    for k in range(5):
        te = folds[k]; tr = np.setdiff1d(idx, te)
        err[te] = Y[te] - pred(X[tr], Y[tr], X[te])
    return err
def mean_pred(Xtr, Ytr, Xte): return np.tile(Ytr.mean(0), (len(Xte), 1))
def lin_pred(Xtr, Ytr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    A = np.column_stack([(Xtr - mu) / sd, np.ones(len(Xtr))]); B = np.column_stack([(Xte - mu) / sd, np.ones(len(Xte))])
    lam = 1e-2; W = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ Ytr); return B @ W
def knn_pred(Xtr, Ytr, Xte, k=5):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    D = (((Xte - mu) / sd)[:, None, :] - ((Xtr - mu) / sd)[None, :, :])**2
    nn = np.argsort(D.sum(2), 1)[:, :k]; return Ytr[nn].mean(1)
res = {}
for nm, p in [('mean', mean_pred), ('linear', lin_pred), ('knn5', knn_pred)]:
    e = cv(p); dist = np.hypot(e[:, 0], e[:, 1])
    res[nm] = dict(rms_dist=float(np.sqrt(np.mean(dist**2))), median_dist=float(np.median(dist)),
                   R2_v1=float(1 - np.mean(e[:, 0]**2) / np.var(v1)), R2_v2=float(1 - np.mean(e[:, 1]**2) / np.var(v2)),
                   frac_within_0p01=float(np.mean(dist < 0.01)))
out['feature_cv'] = res
# typical distance between neighbouring catalog orbits (how sharp a prediction has to be)
P = Y; D = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1]); np.fill_diagonal(D, 9)
out['nearest_neighbour_dist_median'] = float(np.median(D.min(1)))
json.dump(out, open(os.path.join(H, 'period_law.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
