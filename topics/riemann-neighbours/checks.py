"""Checks for riemann-neighbours/README.md. Runs compute.py (about 3 s), then tests each claim.
Polish checks (marked POLISH) redo a result a second, independent way or at higher precision."""
import os, runpy, itertools, numpy as np, mpmath as mp, sympy as sp, networkx as nx
HERE = os.path.dirname(os.path.abspath(__file__))
G = runpy.run_path(os.path.join(HERE, 'compute.py'), run_name='compute')
R = G['R']
results = []
def check(name, ok, note=''):
    results.append((name, bool(ok))); print(('PASS ' if ok else 'FAIL ') + name + (' : ' + note if note else ''))

# G: graphs
check('G1 Bass determinant formula (Petersen)', R['G1_bass_petersen_relerr'] < 1e-12, f"rel err {R['G1_bass_petersen_relerr']:.1e}")
check('G2 the k = -1 fold over prime cycles equals det(I - uB), a polynomial of degree 30', R['G2_fold_equals_det'] and R['G2_det_degree'] == 30 and R['G2_fold_coeffs_31_34'] == [0, 0, 0, 0])
check('G3 Petersen: nontrivial non-backtracking eigenvalues have size sqrt 2, all khat in [-2, 2]',
      R['G3_petersen_nontrivial_mu_abs'] == [round(np.sqrt(2), 8)] and all(abs(x) <= 2 for x in R['G3_petersen_khat']))
g4 = R['G4_graphs']
check('G4 two-cluster graph: lambda2 > 2 sqrt 2 and poles off the circle; random graph: none',
      g4['twocluster60']['lambda2'] > 2*np.sqrt(2) and g4['twocluster60']['off_circle'] > 0 and g4['random60']['off_circle'] == 0 and g4['random60']['ramanujan'],
      f"lambda2 = {g4['twocluster60']['lambda2']:.4f}, off-circle poles {g4['twocluster60']['off_circle']}")
check('G5 Bass formula on all three graphs', max(v['bass_relerr'] for v in g4.values()) < 1e-10)
# POLISH: prime-cycle counts by direct enumeration of closed non-backtracking walks of length 5 and 6 on Petersen (no matrices)
P = nx.petersen_graph()
def nb_closed(L):
    cnt = 0
    for v0 in P:
        stack = [(v0, None, 0)]
        while stack:
            v, prev, l = stack.pop()
            if l == L:
                cnt += (v == v0); continue
            for w in P.neighbors(v):
                if w != prev: stack.append((w, v, l+1))
    return cnt
# tailless condition: for cycles of length below girth*2 every closed NB walk of length 5, 6 is tailless (girth 5)
check('G6 POLISH: closed non-backtracking walks of length 5 and 6 counted by search match tr B^5, tr B^6',
      [nb_closed(5), nb_closed(6)] == R['G2_trace_counts'][4:6], str([nb_closed(5), nb_closed(6)]))

# B: Bartholdi's bump dial is d in disguise
check('B1 Bartholdi determinant formula holds (t = 0, 0.3, 0.7, 1.4)', R['B1_bartholdi_relerr'] < 1e-10)
u, lam, t, q = sp.symbols('u lambda t q')
fac = 1 - lam*u + (1-t)*(q + t - 1 + 1 - 1 + 0)*0  # placeholder removed below
fac = 1 - lam*u + (1-t)*((q+1) - (1-t))*u**2       # D = q+1 on a (q+1)-regular graph
check('B2 per eigenvalue the bump dial only changes d: d_eff = (1 - t)(q + t)', sp.simplify(sp.expand(fac) - (1 - lam*u + (1-t)*(q+t)*u**2)) == 0)

# W: curves
check('W1 Hasse: a_p^2 < 4p for all 300 good primes below 2000 (y^2 = x^3 + x + 1)', R['W1_all_hasse'] and R['W1_n_primes'] == 300)
check('W2 points over F_{p^2} counted directly equal p^2 + 1 - (a_p^2 - 2p) from the dial point', all(v['direct'] == v['from_dial_point'] for v in R['W2_fp2'].values()))
check('W3 degree form m^2 - a_p m n + p n^2 is positive definite (minimum 1 over a box)', all(v == 1 for v in R['W3_degree_form_min'].values()))
m, n, k, d = sp.symbols('m n k d')
check('W4 4(m^2 + k m n + d n^2) = (2m + k n)^2 + (4d - k^2) n^2 (positive definite exactly on the elliptic side)',
      sp.expand(4*(m**2 + k*m*n + d*n**2) - ((2*m + k*n)**2 + (4*d - k**2)*n**2)) == 0)

# L: Lee-Yang on the cube
check('L1 ferromagnetic 9-spin cube: every zero on |z| = 1', R['L1_ferro_root_abs_dev'] < 1e-9, f"{R['L1_ferro_root_abs_dev']:.1e}")
check('L1b antiferromagnetic: zeros leave the circle', R['L1_anti_root_abs_dev'] > 0.1)
check('L2 two directions: 1 + a(z1 + z2) + z1 z2 = 1 + a (z1 o_k z2) with k = 1/a', R['L2_two_dir_identity'])
z, a = sp.symbols('z a', real=True)
disc2 = sp.discriminant(z**2 + 2*a*z + 1, z)
check('L2b uniform field: z^2 + 2az + 1 has discriminant 4a^2 - 4, elliptic side exactly for |a| <= 1', sp.simplify(disc2 - (4*a**2 - 4)) == 0)
check('L3 three directions: weights 1, a^2, a^2, 1 are a product form only at a^2 = 1', all((v > 0) == (a_ != '1.0') for a_, v in R['L3_three_dir_product_form_residual'].items()))
# Hankel test: corner weights c_j = binom(n,j) a^{j(n-j)}/binom(n,j) = a^{j(n-j)} follow an order-2 recurrence (all a (k,d) fold can make) iff 3x3 Hankel dets vanish
ok = True; hv = []
for nn in range(4, 9):
    c = [sp.Rational(1, 2)**(j*(nn-j)) for j in range(nn+1)]
    H = sp.Matrix(3, 3, lambda i, j: c[i+j])
    hv.append(float(H.det())); ok &= H.det() != 0
check('L3b uniform coupling a = 1/2, n = 4..8: Hankel determinants nonzero, so no (k, d) fold gives these weights', ok, ', '.join(f'{x:.2e}' for x in hv))
check('L4 coupled prime cube (primes to 19), ferromagnetic, z_p = p^(1/2 - s): every zero found has Re s = 1/2',
      R['L4_ferro_max_dev_from_half'] < 1e-10 and len(R['L4_ferro_zeros']) >= 50, f"{len(R['L4_ferro_zeros'])} zeros, max |Re s - 1/2| = {R['L4_ferro_max_dev_from_half']:.1e}")
check('L4b antiferromagnetic couplings: zeros leave the line', R['L4_anti_max_dev_from_half'] > 0.1)
check('L4c without the sqrt p weights the zeros sit on Re s = 0', R['L4_unscaled_zeros_max_dev_from_0'] < 1e-10)
# POLISH: refine one ferromagnetic zero at 40 digits with mpmath, from the same weights
W = G['cube_weights'](G['aP']); logn = G['logn']
mp.mp.dps = 40
aPm = [[mp.mpf(float(x)) for x in row] for row in G['aP']]; PR = G['PR']
Wm = []
for S in G['subsets']:
    w = mp.mpf(1)
    for i in range(len(PR)):
        for j in range(len(PR)):
            if S[i] and not S[j]: w *= aPm[i][j]
    Wm.append(w)
Lm = [mp.fsum(mp.log(PR[i]) for i in range(len(PR)) if S[i]) for S in G['subsets']]
fz = lambda s: mp.fsum(w*mp.exp((mp.mpf(1)/2 - s)*l) for w, l in zip(Wm, Lm))
z0 = complex(*R['L4_ferro_zeros'][len(R['L4_ferro_zeros'])//2])
zr = mp.findroot(fz, mp.mpc(z0.real, z0.imag))
check('L5 POLISH: a ferromagnetic cube zero refined at 40 digits stays on Re s = 1/2', abs(zr.real - mp.mpf(1)/2) < mp.mpf(10)**-30, f"zero {mp.nstr(zr, 12)}")

# H: Davenport-Heilbronn
check('H1 Davenport-Heilbronn function: zero off the line at 0.808517 + 85.699348 i', abs(R['H1_DH_zero'][0] - 0.808517) < 1e-5 and abs(R['H1_DH_zero'][1] - 85.699348) < 1e-5, str(R['H1_DH_zero']))
check('H2 its completed function is real on the critical line (it has a functional equation)', all(abs(x - 1) < 1e-12 for x in R['H2_DH_phase_on_line']))
mp.mp.dps = 50
chi = {1: 1, 2: 1j, 3: -1j, 4: -1}; chib = {a_: complex(v).conjugate() for a_, v in chi.items()}
L5f = lambda s, c: mp.power(5, -s)*sum(c[a_]*mp.zeta(s, mp.mpf(a_)/5) for a_ in range(1, 5))
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
f = lambda s: (1 - 1j*kap)/2*L5f(s, chi) + (1 + 1j*kap)/2*L5f(s, chib)
zz = mp.findroot(f, mp.mpc(*R['H1_DH_zero']))
check('H3 POLISH: the off-line zero holds at 50 digits', abs(f(zz)) < mp.mpf(10)**-40 and zz.real > 0.8, mp.nstr(zz, 15))
check('D1 L(s, chi_-4): first zero on the line at t = 6.0209 (an open neighbour, checked only numerically)', abs(R['D1_L_chi4_first_zero'][0] - 0.5) < 1e-12 and abs(R['D1_L_chi4_first_zero'][1] - 6.020949) < 1e-5)

# S: Selberg's local factor
s_, l_ = sp.symbols('s lambda', real=True)
roots = sp.solve(s_**2 - s_ + l_, s_)
check('S1 s^2 - s + lambda: roots are 1/2 +- sqrt(1/4 - lambda), so on Re s = 1/2 exactly when lambda >= 1/4',
      sp.simplify(sum(roots) - 1) == 0 and sp.simplify(roots[0]*roots[1] - l_) == 0 and sp.discriminant(s_**2 - s_ + l_, s_) == 1 - 4*l_)
check('S2 Selberg factor as a dial point (k, d) = (-1, lambda): khat = -1/sqrt(lambda), |khat| <= 2 iff lambda >= 1/4',
      all((abs(-1/np.sqrt(L)) <= 2) == (L >= 0.25) for L in [0.1, 0.2, 0.25, 0.3, 3/16, 975/4096, 91.1]))

# L6: the Lee-Yang coupling is our k dial placed on edges (the repo's "dial gas", topics/lee-yang)
import itertools, fractions
Fr = fractions.Fraction
rng6 = np.random.default_rng(3)
n6 = 5
A6 = {(i, j): Fr(int(rng6.integers(2, 9)), int(rng6.integers(10, 20))) for i in range(n6) for j in range(i+1, n6)}
a6 = lambda i, j: A6[(min(i, j), max(i, j))]
zz6 = [Fr(int(rng6.integers(1, 7)), 7) for _ in range(n6)]
lhs6 = rhs6 = Fr(0)
for S in itertools.product([0, 1], repeat=n6):
    w_cut = Fr(1)
    for i in range(n6):
        for j in range(n6):
            if S[i] and not S[j]: w_cut *= a6(i, j)
    zprod = Fr(1)
    for i in range(n6):
        if S[i]: zprod *= zz6[i]
    lhs6 += zprod * w_cut
    resc = Fr(1); edge = Fr(1)
    for i in range(n6):
        if S[i]:
            for j in range(n6):
                if j != i: resc *= a6(i, j)
    for i in range(n6):
        for j in range(i+1, n6):
            if S[i] and S[j]: edge *= a6(i, j)**-2
    rhs6 += zprod * resc * edge
check('L6 Lee-Yang coupling = the k dial on each edge: cut weights equal dial-gas weights k_ij = a_ij^-2 after z_i -> z_i * prod_j a_ij (exact, 5 spins, random couplings)', lhs6 == rhs6)

npass = sum(ok for _, ok in results)
print(f'\n{npass} of {len(results)} pass')
