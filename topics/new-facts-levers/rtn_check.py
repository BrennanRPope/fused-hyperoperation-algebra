"""Does monogamy (I3 <= 0) survive for real random tensor networks at small bond dimension?
The repo's theorem (holographic-entropy-cube, third pass) covers the ANNEALED Renyi-2 entropy (Ising free energy).
Here: exact states from random Gaussian tensors on four small networks; Renyi-2 annealed vs quenched, and von Neumann."""
import numpy as np, itertools, json
rng = np.random.default_rng(7)
def rnd(*shape): return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
def state(net, D):
    if net == 'single':   # one tensor, legs A B C O
        return rnd(D, D, D, D)
    if net == 'chain2':   # X(A,B,e) Y(C,O,e)
        return np.einsum('abe,coe->abco', rnd(D, D, D), rnd(D, D, D))
    if net == 'ring4':    # four nodes in a ring, one leg each
        return np.einsum('aps,bpq,cqr,ors->abco', rnd(D, D, D), rnd(D, D, D), rnd(D, D, D), rnd(D, D, D))
    if net == 'tree5':    # centre joined to four nodes, each with a leg
        return np.einsum('pqrs,ap,bq,cr,os->abco', rnd(D, D, D, D), rnd(D, D), rnd(D, D), rnd(D, D), rnd(D, D))
SUBS = [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
def moments(psi):
    """return Tr rho_I^2 (unnormalised) and Tr rho (norm^2), plus normalised S2, SvN for each subset"""
    D = psi.shape[0]; z0 = np.vdot(psi, psi).real; out = {}
    for I in SUBS:
        J = [i for i in range(4) if i not in I]
        M = np.transpose(psi, list(I) + J).reshape(D ** len(I), -1)
        sv = np.linalg.svd(M, compute_uv=False) ** 2
        p = sv / sv.sum(); p = p[p > 1e-15]
        out[I] = ((sv ** 2).sum(), -np.log((p ** 2).sum()), -(p * np.log(p)).sum())
    return z0, out
def I3(S):  # S: dict subset -> entropy
    return S[(0,)] + S[(1,)] + S[(2,)] - S[(0, 1)] - S[(0, 2)] - S[(1, 2)] + S[(0, 1, 2)]
res = {}
for net in ['single', 'chain2', 'ring4', 'tree5']:
    for D in [2, 3]:
        Z1 = {I: 0.0 for I in SUBS}; Z0 = 0.0; q2 = []; qv = []; T = 3000
        for _ in range(T):
            z0, m = moments(state(net, D)); Z0 += z0 ** 2
            for I in SUBS: Z1[I] += m[I][0]
            q2.append(I3({I: m[I][1] for I in SUBS})); qv.append(I3({I: m[I][2] for I in SUBS}))
        ann = I3({I: -np.log(Z1[I] / Z0) for I in SUBS})
        q2 = np.array(q2); qv = np.array(qv)
        res[f'{net},D={D}'] = dict(annealed_I3_renyi2=float(ann), quenched_mean_I3_renyi2=float(q2.mean()),
                                   frac_samples_I3_renyi2_pos=float((q2 > 1e-9).mean()),
                                   mean_I3_vN=float(qv.mean()), frac_samples_I3_vN_pos=float((qv > 1e-9).mean()))
        print(net, D, {k: round(v, 4) for k, v in res[f'{net},D={D}'].items()})
json.dump(res, open('rtn_check.json', 'w'), indent=1)
