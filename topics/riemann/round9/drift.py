"""Round 9: are the small dial shifts structural or regulator artefacts?

For each regulator (taper, cutoff X), find the best height dial s (the line the orbit
sum is read on), odds t and amplitude A on the training zeros, and see whether the
optimum drifts back to the truth (s = 1/2, t = 1, A = 1) as more orbits are kept."""
import json
import numpy as np
from stringlib import String, score

gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
TR = np.arange(0, 500); TE = np.arange(500, NZ)
OUT = []
for taper, X in [('sharp', 30), ('sharp', 100), ('sharp', 300), ('sharp', 1000), ('sharp', 3000),
                 ('gauss', 50), ('gauss', 200), ('gauss', 500)]:
    row = dict(taper=taper, X=X)
    for name, vals in [('s', np.arange(0.30, 0.701, 0.025)), ('A', np.arange(0.8, 1.401, 0.025))]:
        tr, te = [], []
        for v in vals:
            st = String(taper=taper, X=float(X), **{name: float(v)})
            E = st.levels(NZ)
            tr.append(score(E, gam, TR)['rms']); te.append(score(E, gam, TE)['rms'])
        tr = np.array(tr); te = np.array(te)
        row[name] = dict(vals=vals.tolist(), train=tr.tolist(), test=te.tolist(),
                         best_train=float(vals[np.argmin(tr)]), best_test=float(vals[np.argmin(te)]))
    print(f"{taper:5s} X={X:5d}: best s train {row['s']['best_train']:.3f} test {row['s']['best_test']:.3f} | "
          f"best A train {row['A']['best_train']:.3f} test {row['A']['best_test']:.3f}")
    OUT.append(row)
json.dump(OUT, open('drift.json', 'w'))
