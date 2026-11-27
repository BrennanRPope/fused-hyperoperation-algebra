"""DRAFT. Krein's theorem on the catalog's stable orbits: a stable orbit whose two elliptic pairs have the same Krein
sign (definite) cannot turn complex-unstable under a small Hamiltonian perturbation; one with opposite signs (mixed)
can, once the perturbation is about as big as the gap between the two pairs' angles (Krein 1950; Gelfand-Lidskii
1955; Moser 1958).  Here: refine every catalog E-E orbit as a fixed point of the section map, keep the ones whose residual is below
1e-9 and that are still E-E, and record each pair's angle and Krein sign. Writes krein.json."""
import json, os, numpy as np
import cb3, cz
from czcatalog import reduced_with_basis

def refine_one(o):
    import section as S
    E = o['E']; s0 = np.array([2, 0, o['v1'], o['v1']])
    try:
        k = S.returns_per_period(s0, E, o['T'])
        s, res = S.refine(s0, E, k)
        _, Tper = S.F(s, E, k)
        return S.to_state(s, E), float(Tper), res
    except Exception:
        return None

if __name__ == '__main__':
    import multiprocessing as mp
    from collections import Counter
    d = json.load(open(os.path.join(cb3.HERE, 'cz.json')))
    fl = {o['cls']: o for o in json.load(open(os.path.join(cb3.HERE, 'floquet.json'))) if o.get('ok')}
    ee = [r for r in d if r['r'] == 2]
    with mp.Pool(4) as pool: refined = pool.map(refine_one, [fl[r['cls']] for r in ee])
    out, fate = [], Counter()
    for r, rf in zip(ee, refined):
        if rf is None: fate['no refinement'] += 1; continue
        X0, Tper, res = rf
        X, Phi = cb3.flow_var(X0, Tper)
        M, V = reduced_with_basis(X0, Phi)
        t1, t2, info = cb3.pair_traces(M); kl = cb3.orbit_class(t1, t2)
        if res > 1e-9: fate['residual > 1e-9'] += 1; continue
        fate[kl] += 1
        if kl != 'E-E': continue
        pairs = cz.pairs_from_traces(M, [t1, t2], V, cz.Jmat(6))
        ph = [p['phi'] for p in pairs]
        out.append(dict(cls=r['cls'], word=r['word'], T=Tper, residual=res,
                        krein='definite' if pairs[0]['s'] == pairs[1]['s'] else 'mixed',
                        phi=ph, s=[p['s'] for p in pairs], theta=[p['theta'] for p in pairs],
                        collision_gap=abs(ph[1] - ph[0]), edge_gap=min(min(ph), np.pi - max(ph))))
    json.dump(dict(fate=fate, orbits=out), open(os.path.join(cb3.HERE, 'krein.json'), 'w'), indent=1)
    print('catalog E-E after refinement:', dict(fate))
    print('Krein types of confirmed stable orbits:', Counter(x['krein'] for x in out))
    print('sign patterns:', Counter(tuple(sorted(x['s'])) for x in out))
