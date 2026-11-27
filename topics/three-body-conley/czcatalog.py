"""DRAFT. Conley-Zehnder iteration data for the Li-Liao catalog: per Floquet pair, the Krein-normalised angle,
the iterate indices muhat(1..K), a decode of the number of elliptic pairs from those integers alone, a parity check
against the fixed-point index, and a strong-stability (Krein) perturbation test of the stable orbits.
Writes cz.json."""
import json, os, numpy as np, multiprocessing as mp
import cb3, cz

K = 60
SRC = os.path.join(cb3.HERE, 'floquet.json')

def reduced_with_basis(X0, Phi):
    """same reduction as cb3.reduced_block, also returning the 12x4 basis V that represents the reduced space."""
    G = cb3.grad_conserved(X0)
    _, s, Vt = np.linalg.svd(G); Kb = Vt[6:].T
    f = cb3.field(X0); xi = cb3.rot_gen(X0)
    A = np.column_stack([f, xi]); A = Kb @ np.linalg.lstsq(Kb, A, rcond=None)[0]
    Qa, _ = np.linalg.qr(A)
    P = Kb @ Kb.T - Qa @ Qa.T
    w, U = np.linalg.eigh(P); V = U[:, w > 0.5]
    B = np.column_stack([A, V])
    C = np.linalg.lstsq(B, Phi @ V, rcond=None)[0]
    return C[2:, :], V

def strong_test(M, Om, rng, eps_grid, nsamp=60):
    """fraction of random Hamiltonian perturbations exp(eps Om^-1 S) M that are still E-E, per eps."""
    Oi = np.linalg.inv(Om); out = []
    Ss = []
    for _ in range(nsamp):
        A = rng.normal(size=(4, 4)); S = (A + A.T) / 2; Ss.append(S / np.linalg.norm(S, 2))
    from scipy.linalg import expm
    for e in eps_grid:
        ok = 0
        for S in Ss:
            Mp = expm(e * Oi @ S) @ M
            lam = np.linalg.eigvals(Mp)
            ok += int(np.all(np.abs(np.abs(lam) - 1) < 1e-6))
        out.append(ok / nsamp)
    return out

def one(o):
    if not o.get('ok'): return None
    X0 = cb3.suvakov(o['v1'], o['v2'])
    X, Phi = cb3.flow_var(X0, o['T'])
    M, V = reduced_with_basis(X0, Phi)
    J = cz.Jmat(6); Om = V.T @ J @ V
    sym_err = float(np.abs(M.T @ Om @ M - Om).max() / (np.abs(Om).max() * np.abs(M).max() ** 2))
    taus = [complex(a, b) for a, b in o['tau']]
    pairs = cz.pairs_from_traces(M, taus, V, J)
    m = np.arange(1, K + 1)
    mu = [int(x) for x in cz.muhat(pairs, m)]
    rec = dict(cls=o['cls'], word=o['word'], T=o['T'], klass=o['klass'], rel_pal=o['rel_pal'], tau=o['tau'],
               idx=o['idx'], sym_err=sym_err,
               pairs=[{k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in p.items()} for p in pairs],
               muhat=mu, mean_hat=float(cz.mean_hat(pairs)), r=sum(p['kind'] == 'E' for p in pairs),
               r_decoded=cz.decode_r(mu))
    rec['parity_ok'] = all(((-1) ** mu[n - 1]) == o['idx'][n - 1] for n in (1, 2, 3, 4))
    if rec['r'] == 2:
        rec['krein'] = 'definite' if pairs[0]['s'] == pairs[1]['s'] else 'mixed'
        ph = sorted(p['phi'] for p in pairs)
        rec['collision_gap'] = float(ph[1] - ph[0])                       # distance to a Krein collision
        rec['edge_gap'] = float(min(ph[0], np.pi - ph[1]))               # distance to +1 or -1
        rec['strong'] = strong_test(M, Om, np.random.default_rng(sum(map(ord, o['cls']))),
                                    [0.003, 0.01, 0.03, 0.1, 0.3])
    return rec

if __name__ == '__main__':
    d = json.load(open(SRC))
    with mp.Pool(4) as p: out = [r for r in p.map(one, d, chunksize=4) if r]
    json.dump(out, open(os.path.join(cb3.HERE, 'cz.json'), 'w'))
    from collections import Counter
    print(len(out), 'orbits; median relative symplectic error', np.median([r['sym_err'] for r in out]))
    print('parity ok', sum(r['parity_ok'] for r in out))
    print('class vs r', Counter((r['klass'], r['r']) for r in out))
    print('r vs decoded', Counter((r['r'], r['r_decoded']) for r in out))
    print('E-E krein', Counter(r.get('krein') for r in out if r['r'] == 2))
