"""DRAFT. The 7 same-word twin pairs (from ../three-body-word-formula/twins.json): each orbit refined as a fixed
point of the section map, its Floquet type, and its fixed-point index (Jacobian formula and value-only count).
Twins born together at a fold must have opposite indices. Writes twins.json."""
import json, os, numpy as np, multiprocessing as mp
import cb3, section as S, degree as D

def one(args):
    name, which, (v1, v2, T) = args
    X0 = cb3.suvakov(v1, v2); E = cb3.energy(X0); s0 = np.array([2, 0, v1, v1])
    try:
        k = S.returns_per_period(s0, E, T)
        s, res = S.refine(s0, E, k)
        J = S.jac(s, E, k)
    except Exception as ex:
        return dict(name=name, which=which, ok=False, err=str(ex))
    t1, t2, info = cb3.pair_traces(J)
    G = J - np.eye(4); B = np.linalg.inv(G); B /= np.abs(B).max()
    deg, _ = D.robust_degree(lambda x: S.F(x, E, k)[0] - x, s, 1e-6, B=B)
    lam = np.abs(np.linalg.eigvals(J))
    return dict(name=name, which=which, ok=True, k=k, residual=res, E=E, T=T, v=[v1, v2],
                tau=[[t1.real, t1.imag], [t2.real, t2.imag]], klass=cb3.orbit_class(t1, t2),
                idx_jac=cb3.index_from_traces(t1, t2, 1), deg=deg, logmaxlam=float(np.log(lam.max())),
                residues=[cb3.residue(t1).real, cb3.residue(t2).real])

if __name__ == '__main__':
    tw = json.load(open(os.path.join(cb3.HERE, '..', 'three-body-word-formula', 'twins.json')))
    jobs = [(n, w, tw[n][w]) for n in tw for w in ('cat', 'twin')]
    with mp.Pool(4) as p: out = p.map(one, jobs)
    json.dump(out, open('twins.json', 'w'), indent=1)
    for o in out:
        print(o['name'], o['which'], o.get('klass'), 'res %.0e' % o.get('residual', -1), 'idx', o.get('idx_jac'), 'deg', o.get('deg'), 'logλ %.3f' % o.get('logmaxlam', 0))
