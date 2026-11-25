import numpy as np, json
np.set_printoptions(precision=6, suppress=True)
out = {}

def h_of_k(v, w, k):
    return v + w * np.exp(1j * k)

# 1. Mult-add winding: steps delta_m = h(k_{m+1})/h(k_m) - 1
def multadd_winding(v, w, M):
    ks = 2 * np.pi * np.arange(M + 1) / M
    h = h_of_k(v, w, ks)
    d = h[1:] / h[:-1] - 1          # the mult-add steps
    total = np.prod(1 + d) - 1      # the mult-add total of all steps
    turns = np.sum(np.angle(1 + d)) / (2 * np.pi)
    return total, turns

out['winding'] = []
for v, w in [(0.5, 1.0), (1.0, 0.5), (0.9, 1.0), (1.0, 0.9)]:
    for M in [4, 8, 64, 1000]:
        tot, turns = multadd_winding(v, w, M)
        out['winding'].append(dict(v=v, w=w, M=M, total=abs(tot), turns=round(float(turns), 12)))
        print('winding', v, w, M, abs(tot), turns)

# Zak (Berry) phase of the lower band by a discrete Wilson loop
def zak(v, w, M=2000):
    ks = 2 * np.pi * np.arange(M) / M
    vecs = []
    for k in ks:
        H = np.array([[0, h_of_k(v, w, k)], [np.conj(h_of_k(v, w, k)), 0]])
        e, U = np.linalg.eigh(H)
        vecs.append(U[:, 0])
    prod = 1
    for m in range(M):
        prod *= np.vdot(vecs[m], vecs[(m + 1) % M])
    return -np.angle(prod)
out['zak'] = {f'{v},{w}': float(zak(v, w)) for v, w in [(0.5, 1), (1, 0.5)]}
print('zak', out['zak'])

# 2. Finite open chain, A/B ordering. h[n,n]=v_n (A_n-B_n), h[n+1,n]=w_n (B_n-A_{n+1})
def blocks(vs, ws):
    N = len(vs)
    h = np.zeros((N, N))
    for n in range(N):
        h[n, n] = vs[n]
        if n + 1 < N:
            h[n + 1, n] = ws[n]
    return h

def full_H(vs, ws, onsite=None):
    h = blocks(vs, ws)
    N = len(vs)
    H = np.block([[np.zeros((N, N)), h], [h.T, np.zeros((N, N))]])
    if onsite is not None:
        H = H + np.diag(onsite)
    return H

def realspace_winding(vs, ws):
    h = blocks(vs, ws)
    N = len(vs)
    U, S, Vt = np.linalg.svd(h)
    QAB = U @ Vt
    QBA = QAB.T
    X = np.diag(np.arange(N, dtype=float))
    C = QBA @ (QAB @ X - X @ QAB)
    mid = range(N // 4, 3 * N // 4)
    return float(np.real(np.mean([C[n, n] for n in mid])))

N = 40
for v, w in [(0.5, 1.0), (1.0, 0.5)]:
    E = np.linalg.eigvalsh(full_H([v] * N, [w] * (N)))
    E = np.sort(np.abs(E))
    print('clean', v, w, 'smallest |E|', E[:4], 'gap', 2 * abs(w - v), 'rs winding', realspace_winding([v] * N, [w] * N))
    out[f'clean_{v}_{w}'] = dict(smallest=[float(x) for x in E[:4]], rsw=realspace_winding([v] * N, [w] * N))

# edge state profile, v=0.5, w=1, N=40
v, w = 0.5, 1.0
H = full_H([v] * N, [w] * N)
e, U = np.linalg.eigh(H)
idx = np.argsort(np.abs(e))[:2]
# combine the two near-zero modes into left (A) and right (B) edge states
psiA = U[:, idx[0]] + U[:, idx[1]]
psiB = U[:, idx[0]] - U[:, idx[1]]
psiA /= np.linalg.norm(psiA); psiB /= np.linalg.norm(psiB)
if np.sum(psiA[:N] ** 2) < 0.5:
    psiA, psiB = psiB, psiA
A = psiA[:N]; B = psiA[N:]
sgn = np.sign(A[0])
A *= sgn; B *= sgn
pred = np.array([(-v / w) ** n for n in range(N)])
pred /= np.linalg.norm(pred)
print('edge A amps first 8', A[:8])
print('pred      first 8', pred[:8])
print('max |B| on left state', np.max(np.abs(B[:N // 2])), 'max dev from pred', np.max(np.abs(A - pred)))
print('energies of near-zero modes', e[idx])
out['edge'] = dict(A=[float(x) for x in A[:12]], pred=[float(x) for x in pred[:12]], B=[float(x) for x in B[:12]],
                   E=[float(x) for x in e[idx]], maxdev=float(np.max(np.abs(A - pred))))

# 3. Disorder (chiral-preserving: random hoppings only)
rng = np.random.default_rng(7)
dis = []
for (v, w) in [(0.5, 1.0), (1.0, 0.5)]:
    for eps in [0.0, 0.3, 0.6, 0.9]:
        rows = []
        for trial in range(200):
            vs = v * (1 + eps * rng.uniform(-1, 1, N))
            ws = w * (1 + eps * rng.uniform(-1, 1, N))
            E = np.sort(np.abs(np.linalg.eigvalsh(full_H(vs, ws))))
            nz = int(np.sum(E < 1e-6))
            rsw = realspace_winding(vs, ws)
            pred_nu = 1 if np.mean(np.log(ws)) > np.mean(np.log(vs)) else 0
            rows.append((nz, rsw, pred_nu, E[0], E[2]))
        nzs = [r[0] for r in rows]; rsws = np.array([r[1] for r in rows])
        d = dict(v=v, w=w, eps=eps, zero_modes_counts={str(k): nzs.count(k) for k in set(nzs)},
                 rsw_mean=float(rsws.mean()), rsw_min=float(rsws.min()), rsw_max=float(rsws.max()),
                 frac_round_match=float(np.mean([round(r[1]) == r[2] for r in rows])),
                 median_E0=float(np.median([r[3] for r in rows])), median_E2=float(np.median([r[4] for r in rows])))
        dis.append(d); print(d)
out['disorder'] = dis

# 4. Break chiral symmetry: on-site potential
rng = np.random.default_rng(3)
v, w = 0.5, 1.0
for W in [0.0, 0.1, 0.3]:
    E = np.sort(np.linalg.eigvalsh(full_H([v] * N, [w] * N, onsite=W * rng.uniform(-1, 1, 2 * N))))
    near = E[np.argsort(np.abs(E))[:2]]
    print('onsite', W, 'two nearest-zero energies', near)
    out[f'onsite_{W}'] = [float(x) for x in near]

# loop data for picture
ks = np.linspace(0, 2 * np.pi, 49)
loops = []
for name, v, w in [('trivial', 1.0, 0.5), ('topological', 0.5, 1.0)]:
    for k in ks:
        z = h_of_k(v, w, k)
        loops.append(dict(case=name, k=round(float(k), 4), re=round(float(z.real), 4), im=round(float(z.imag), 4)))
out['loops'] = loops
json.dump(out, open('results.json', 'w'), indent=1)
