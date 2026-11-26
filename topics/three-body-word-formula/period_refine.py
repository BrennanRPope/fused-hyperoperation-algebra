"""DRAFT. Can the word itself explain the ~3% leftover in period-per-letter? 5-fold CV regressions of T*/Lf on word features."""
import json, os, re, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
def inverse(w): return w[::-1].swapcase()
def syll(w):
    """split the cyclic word into runs of alternating case (figure-eight-like) vs repeated case"""
    return w
def feats(w):
    n = len(w); c = w + w[0]
    same_case = sum(c[i].isupper() == c[i+1].isupper() for i in range(n)) / n     # e.g. 'ab' or 'BA': both letters same direction
    doubled = sum(c[i] == c[i+1] for i in range(n)) / n                           # 'aa','bb': circling one collision twice
    alt = sum(c[i].lower() != c[i+1].lower() for i in range(n)) / n               # a then b (switching collision)
    f8 = len(re.findall(r'(?=abAB|baBA|ABab|BAba|aBAb|bABa|AbaB|BabA)', w + w[:3])) / n
    return [same_case, doubled, alt, f8]
X = np.array([feats(o['word']) for o in d]); r = np.array([o['Tstar'] / o['Lf'] for o in d])
names = ['same-direction pairs', 'doubled letters', 'switches a<->b', 'figure-eight syllables']
out = {'corr': {nm: float(np.corrcoef(X[:, k], r)[0, 1]) for k, nm in enumerate(names)}}
rng = np.random.default_rng(3); idx = rng.permutation(len(d)); folds = np.array_split(idx, 5)
def cvfit(cols):
    pred = np.zeros_like(r)
    for k in range(5):
        te = folds[k]; tr = np.setdiff1d(idx, te)
        A = np.column_stack([X[tr][:, cols], np.ones(len(tr))]); w = np.linalg.lstsq(A, r[tr], rcond=None)[0]
        pred[te] = np.column_stack([X[te][:, cols], np.ones(len(te))]) @ w
    e = r - pred; return float(np.std(e)), float(np.mean(np.abs(e) / r < 0.01)), pred
base = float(np.std(r))
out['sd_ratio_plain'] = base
for cols, nm in [([0], 'same-direction'), ([1], 'doubled'), ([2], 'switches'), ([3], 'f8'), ([0, 1, 2, 3], 'all four')]:
    sd, f1, pred = cvfit(cols); out['cv_' + nm] = dict(sd=sd, frac_within_1pct=f1)
A = np.column_stack([X, np.ones(len(r))]); w = np.linalg.lstsq(A, r, rcond=None)[0]
out['fit_all_coeffs'] = dict(zip(names + ['const'], map(float, w)))
out['frac_within_1pct_plain'] = float(np.mean(np.abs(r - r.mean()) / r < 0.01))
json.dump(out, open(os.path.join(H, 'period_refine.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
