"""DRAFT. A small set of orbits, the shortest robust one of each Floquet class (plus a second E-E), refined as
fixed points of the section map F^k, with Floquet types from the section Jacobian and the fixed-point index of
F^k and F^2k computed two ways: from the Jacobian and derivative-free (degree on a small 4-cube).
Writes curated.json."""
import json, numpy as np, multiprocessing as mp
import cb3, section as S, degree as D

def margin(o):
    m = []
    for tr, ti in o['tau']:
        if abs(ti) > 1e-9: m.append(abs(ti))
        else: m += [abs(2 - tr), abs(2 + tr)]
    m.append(abs(o['disc']))
    return min(m)

def work(o):
    E = o['E']; s0 = np.array([2, 0, o['v1'], o['v1']])
    try:
        k = S.returns_per_period(s0, E, o['T'])
        s, res = S.refine(s0, E, k)
        J = S.jac(s, E, k)
    except Exception as ex:
        return dict(cls=o['cls'], ok=False, err=str(ex))
    t1, t2, info = cb3.pair_traces(J)
    klass = cb3.orbit_class(t1, t2)
    # cube half-width: well inside the linear zone (scaled so |F^k - id| stays far above noise)
    # Two boxes per iterate: a plain cube, and a box shaped by (J^n - I)^-1 so the count is well conditioned.
    # Either way the count itself uses only values of F; the shape only decides where we sample.
    out = {}
    for n_it, key in ((1, 'F'), (2, 'F2')):
        Jn = np.linalg.matrix_power(J, n_it); G = Jn - np.eye(4)
        B = np.linalg.inv(G); B /= np.abs(B).max()
        g = (lambda x, n_it=n_it: S.F(x, E, n_it * k)[0] - x)
        lam = np.abs(np.linalg.eigvals(Jn)).max()
        out[key + '_cube'] = D.robust_degree(g, s, 1e-5 / max(1.0, lam))[0]
        out[key + '_shaped'] = D.robust_degree(g, s, 1e-6, B=B)[0]
    d1, d2 = out['F_shaped'], out['F2_shaped']
    m1, m2 = out, None
    return dict(cls=o['cls'], word=o['word'], T=o['T'], E=E, k=k, s=s.tolist(), residual=res, ok=True,
                tau=[[t1.real, t1.imag], [t2.real, t2.imag]], klass=klass, klass_catalog=o['klass'],
                idx_jac=[cb3.index_from_traces(t1, t2, 1), cb3.index_from_traces(t1, t2, 2)],
                deg_F=d1, deg_F2=d2, degrees=m1,
                eig_abs=np.abs(np.linalg.eigvals(J)).tolist(), eig_arg=np.angle(np.linalg.eigvals(J)).tolist())

if __name__ == '__main__':
    d = [o for o in json.load(open('floquet.json')) if o['ok']]
    picks = []
    for kl in ['E-E', 'E-H', 'E-IH', 'H-H', 'H-IH', 'IH-IH', 'CU']:
        cand = sorted([o for o in d if o['klass'] == kl and margin(o) > 0.05 and o['v2'] > 0], key=lambda o: o['T'])
        picks += cand[:3]
    with mp.Pool(4) as p: out = p.map(work, picks)
    json.dump(out, open('curated.json', 'w'), indent=1)
    for o in out:
        if not o['ok']: print(o); continue
        print(o['cls'], o['k'], o['klass'], '(catalog %s)' % o['klass_catalog'], 'res %.0e' % o['residual'],
              'jac idx', o['idx_jac'], o['degrees'])
