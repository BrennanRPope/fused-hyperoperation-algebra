"""DRAFT. Floquet-pair types of the 695 Li-Liao orbits (equal masses, zero angular momentum, collinear start)
and of the figure-eight, from exact variational Jacobians. Writes floquet.json."""
import json, sys, os, numpy as np, multiprocessing as mp
import cb3
SRC = os.path.join(cb3.HERE, '..', 'three-body-word-formula', 'data', 'li-liao-695.json')

def one(o):
    X0 = cb3.suvakov(o['v1'], o['v2']); T = o['T']
    try:
        X, Phi = cb3.flow_var(X0, T)
    except RuntimeError:
        return dict(cls=o['cls'], word=o['word'], ok=False)
    M = cb3.reduced_block(X0, Phi)
    t1, t2, info = cb3.pair_traces(M)
    lam = np.abs(np.linalg.eigvals(M))
    rec = dict(cls=o['cls'], word=o['word'], T=T, Tstar=o['Tstar'], v1=o['v1'], v2=o['v2'], E=cb3.energy(X0),
               closure=float(np.abs(X - X0).max()), tau=[[t1.real, t1.imag], [t2.real, t2.imag]],
               klass=cb3.orbit_class(t1, t2), logmaxlam=float(np.log(lam.max())),
               idx=[cb3.index_from_traces(t1, t2, n) for n in (1, 2, 3, 4)], ok=True, **info)
    rec['rel_pal'] = info['palindrome_err'] / max(1.0, float(np.abs(np.poly(M)).max()))
    return rec

if __name__ == '__main__':
    d = json.load(open(SRC))
    with mp.Pool(4) as p: out = p.map(one, d, chunksize=4)
    json.dump(out, open(os.path.join(cb3.HERE, 'floquet.json'), 'w'))
    ok = [o for o in out if o['ok']]
    from collections import Counter
    print(len(out), 'orbits,', len(ok), 'integrated')
    print(Counter(o['klass'] for o in ok))
