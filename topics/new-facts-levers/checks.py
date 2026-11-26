"""Checks for topics/new-facts-levers (DRAFT). Each line: PASS/FAIL, label, numbers.
Run: python3 checks.py   (about 30 s)"""
import itertools, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thresholds import cone, lp_min, threshold, setup, row, fvec
import mpmath as mp
results = []
def check(label, ok, info=''):
    results.append(ok); print(('PASS' if ok else 'FAIL'), label, info)

# T1-T4: corner-sum thresholds over the von Neumann outer cone (SA + SSA + purification) and SA+MMI
expect = {3: -0.5505102572, 4: -0.3819660113, 5: -0.3538307955, 6: -0.2881734520}
for n in [3, 4, 5, 6]:
    A_vn, reps, var = cone(n); A_h, _, _ = cone(n, holo=True)
    kv = threshold(lambda k: lp_min(n, k, A_vn, reps)[0]); kh = threshold(lambda k: lp_min(n, k, A_h, reps)[0])
    check(f'T{n-2} n={n} vN-LP threshold', abs(kv - expect[n]) < 1e-8, f'{kv:.10f}')
    check(f'T{n-2}b n={n} SA+MMI-LP threshold equals vN-LP', abs(kh - kv) < 1e-8, f'{kh:.10f}')

# H1: n = 5 with Hernandez-Cuenca facet orbits 3,4,5 added -> threshold drops to (sqrt5-3)/2
n = 5; N, full, reps, idx, var = setup(n); L = 'ABCDEO'
FACETS = [("ABC ADE BCDE", "A BC DE ABCDE"), ("ABC ABD ACE BDE CDE", "AB AC BD CE DE ABCDE"),
          ("ABC ABD ABE ACD ACE ADE BCE BDE CDE", "AB AC AD BE CE DE BCD ABCE ABDE ACDE")]
def mask(w, perm): return sum(1 << perm[L.index(c)] for c in w)
extra = set()
for l, r in FACETS:
    for perm in itertools.permutations(range(6)):
        rr = row(var, len(reps), [(1, mask(w, perm)) for w in l.split()] + [(-1, mask(w, perm)) for w in r.split()])
        extra.add(tuple(np.round(rr, 9)))
A_mmi, _, _ = cone(n, holo=True); A_h = np.vstack([A_mmi, -np.array(sorted(extra))])
kh5 = threshold(lambda k: lp_min(n, k, A_h, reps)[0])
check('H1 n=5 holographic (SA+MMI+facets 3-5) threshold = (sqrt5-3)/2', abs(kh5 - (np.sqrt(5) - 3) / 2) < 1e-8, f'{kh5:.10f}')
pt4 = np.array([min(bin(m & 15).count('1'), 4 - bin(m & 15).count('1')) for m in reps], float)
check('H2 perfect tensor on 4 regions is in that cone', bool(np.all(A_h @ pt4 <= 1e-9)))
mds = np.array([min(bin(m).count('1'), 2, 6 - bin(m).count('1')) for m in reps], float)
check('H3 MDS-GHZ vector is outside the holographic cone (breaks facets 4 and 5)', bool(np.max(A_h @ mds) > 0.5), f'max violation {np.max(A_h @ mds):.1f}')
check('H4 MDS-GHZ vector is inside the von Neumann cone', bool(np.all(cone(5)[0] @ mds <= 1e-9)))
k0 = -0.36
fv = float(fvec(5, k0, reps) @ mds); fp = float(fvec(5, k0, reps) @ pt4)
check('H5 at k=-0.36: MDS-GHZ corner sum < 0, perfect tensor > 0', fv < 0 < fp, f'{fv:.4f}, {fp:.4f}')

# Q1: build the MDS-GHZ state ([6,2] Reed-Solomon over F_5) and compute von Neumann entropies
q = 5; psi = np.zeros([q] * 6)
for a in range(q):
    for b in range(q):
        psi[tuple([(a + b * x) % q for x in range(5)] + [b])] += 1
psi /= np.linalg.norm(psi)
def S(I):
    I = list(I); J = [i for i in range(6) if i not in I]
    sv = np.linalg.svd(np.transpose(psi, I + J).reshape(q ** len(I), -1), compute_uv=False) ** 2
    sv = sv[sv > 1e-12]; return float(-(sv * np.log(sv)).sum() / np.log(q))
ok = all(abs(S(I) - min(len(I), 2, 6 - len(I))) < 1e-9 for r in range(1, 6) for I in itertools.combinations(range(6), r))
check('Q1 Reed-Solomon state has S_I = min(|I|, 2, 6-|I|) (units of log 5)', ok)

# P1: polish the three closed-form thresholds to 30 digits
mp.mp.dps = 30
r3 = mp.findroot(lambda k: k**2 + 6*k + 3, -0.55); r4 = mp.findroot(lambda k: k**2 + 3*k + 1, -0.38)
r5 = mp.findroot(lambda k: k**4 + 10*k**3 + 20*k**2 + 20*k + 5, -0.35)
check('P1 n=3 root = sqrt6 - 3', abs(r3 - (mp.sqrt(6) - 3)) < mp.mpf(10)**-28, mp.nstr(r3, 25))
check('P2 n=4 root = (sqrt5 - 3)/2 = -1/phi^2', abs(r4 + 1 / mp.phi**2) < mp.mpf(10)**-28, mp.nstr(r4, 25))
check('P3 n=5 quantum root (MDS-GHZ)', abs(r5 - mp.mpf('-0.3538307955')) < 1e-9, mp.nstr(r5, 25))

# R1: annealed Renyi-2 monogamy on small random tensor networks (fixed seed)
rng = np.random.default_rng(7)
def rnd(*s): return rng.standard_normal(s) + 1j * rng.standard_normal(s)
SUBS = [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
def I3(Sd): return Sd[(0,)] + Sd[(1,)] + Sd[(2,)] - Sd[(0, 1)] - Sd[(0, 2)] - Sd[(1, 2)] + Sd[(0, 1, 2)]
D = 2; Z1 = {I: 0.0 for I in SUBS}; Z0 = 0.0; pos = 0; T = 1500; q2 = []
for _ in range(T):
    p = np.einsum('abe,coe->abco', rnd(D, D, D), rnd(D, D, D)); z0 = np.vdot(p, p).real; Z0 += z0 ** 2; s2 = {}
    for I in SUBS:
        J = [i for i in range(4) if i not in I]
        sv = np.linalg.svd(np.transpose(p, list(I) + J).reshape(D ** len(I), -1), compute_uv=False) ** 2
        Z1[I] += (sv ** 2).sum(); pp = sv / sv.sum(); s2[I] = -np.log((pp ** 2).sum())
    q2.append(I3(s2))
ann = I3({I: -np.log(Z1[I] / Z0) for I in SUBS}); q2 = np.array(q2)
check('R1 two-node network D=2: annealed Renyi-2 I3 <= 0 (the theorem)', ann < 0, f'{ann:.4f}')
check('R2 same network: quenched mean I3 <= 0', q2.mean() < 0, f'{q2.mean():.4f}')
check('R3 same network: some single samples break monogamy', (q2 > 1e-9).mean() > 0.1, f'fraction {(q2 > 1e-9).mean():.3f}')

print(f'{sum(results)}/{len(results)} passed')
