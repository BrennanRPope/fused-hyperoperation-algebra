"""Riemann's neighbours on the dial map: the computations behind README.md.
Writes results.json. Groups: G = graphs (Ihara), W = curves (Hasse/Weil), L = Lee-Yang cube, H = Davenport-Heilbronn,
D = Dirichlet L (open neighbour), B = Bartholdi (dial-in-disguise test)."""
import json, itertools, time, sys, numpy as np, networkx as nx, mpmath as mp
T0 = time.time()
def stamp(x): print(f'[{time.time()-T0:6.1f}s] {x}', flush=True)
rng = np.random.default_rng(20261126)
R = {}

# ---------- G. Ihara zeta of graphs ----------
def hashimoto(G):
    E = [(u, v) for u, v in G.edges()] + [(v, u) for u, v in G.edges()]
    idx = {e: i for i, e in enumerate(E)}
    B = np.zeros((len(E), len(E)))
    J = np.zeros((len(E), len(E)))
    for (a, b), i in idx.items():
        J[i, idx[(b, a)]] = 1
        for c in G.neighbors(b):
            if c != a:
                B[i, idx[(b, c)]] = 1
    return B, J, E
def ihara_check(G):
    n, m = G.number_of_nodes(), G.number_of_edges()
    A = nx.to_numpy_array(G); D = np.diag(A.sum(1)); B, J, E = hashimoto(G)
    worst = 0
    for u in [0.13, 0.31 + 0.2j, -0.27 + 0.05j]:
        lhs = np.linalg.det(np.eye(2*m) - u*B)
        rhs = (1 - u*u)**(m - n) * np.linalg.det(np.eye(n) - u*A + u*u*(D - np.eye(n)))
        worst = max(worst, abs(lhs - rhs)/abs(rhs))
    return worst, A, B, J
P = nx.petersen_graph()
w, A, B, J = ihara_check(P)
lam = np.sort(np.linalg.eigvalsh(A))[::-1]
mu = np.linalg.eigvals(B)
R['G1_bass_petersen_relerr'] = w
R['G1_petersen_eigs'] = [round(x, 10) for x in lam]
# prime cycles: N_m = tr B^m, pi(l) by Mobius inversion of N_m = sum_{d|m} d pi(d)
L = 34
Nm = []; Bk = np.eye(B.shape[0])
for mm in range(1, L+1):
    Bk = Bk @ B; Nm.append(int(round(np.trace(Bk))))
pi = {}
for l in range(1, L+1):
    s = Nm[l-1] - sum(d*pi[d] for d in pi if l % d == 0 and d < l)
    pi[l] = s // l
R['G2_trace_counts'] = Nm[:12]
R['G2_prime_cycles'] = {str(l): pi[l] for l in range(1, 21)}
# the k=-1 fold over prime cycles, prod (1 - u^l)^pi(l), truncated at degree L: it is a polynomial of degree 2m=30
# exp of the log series: log prod(1-u^l)^pi(l) = -sum_m N_m u^m / m, so n a_n = -sum_{k<=n} N_k a_{n-k} (exact integers)
coef = [1]
for nn in range(1, L+1):
    acc = -sum(Nm[k-1]*coef[nn-k] for k in range(1, nn+1))
    assert acc % nn == 0
    coef.append(acc // nn)
detpoly = np.poly(B)  # det(I - uB) = sum_k c_k u^k with c_k the char-poly coefficients (highest first)
detpoly = [int(round(x.real)) for x in detpoly]
R['G2_fold_equals_det'] = all(int(coef[i]) == (detpoly[i] if i < len(detpoly) else 0) for i in range(L+1))
R['G2_fold_coeffs_31_34'] = [int(coef[i]) for i in range(31, 35)]
R['G2_det_degree'] = len(detpoly) - 1
# Ramanujan / dial points: nontrivial eigenvalue lambda -> factor 1 - lambda u + q u^2 -> khat = -lambda/sqrt q
q = 2
R['G3_petersen_khat'] = [round(-x/np.sqrt(q), 6) for x in lam[1:]]
R['G3_petersen_nontrivial_mu_abs'] = sorted(set(round(abs(z), 8) for z in mu if abs(abs(z) - 1) > 1e-6 and abs(abs(z) - 2) > 1e-6))
# random 3-regular graphs: fraction Ramanujan, and a two-cluster graph that is not
def two_cluster(nh, seed):
    G1 = nx.random_regular_graph(3, nh, seed=seed); G2 = nx.random_regular_graph(3, nh, seed=seed+1)
    G = nx.disjoint_union(G1, G2)
    (a, b) = list(G1.edges())[0]; e2 = list(G2.edges())[0]; c, d = e2[0] + nh, e2[1] + nh
    G.remove_edge(a, b); G.remove_edge(c, d); G.add_edge(a, c); G.add_edge(b, d)
    return G
rams = []
for s in range(40):
    G = nx.random_regular_graph(3, 200, seed=1000+s)
    ev = np.sort(np.linalg.eigvalsh(nx.to_numpy_array(G)))[::-1]
    rams.append(max(abs(ev[1]), abs(ev[-1])) <= 2*np.sqrt(2) + 1e-9)
R['G3_random_200_ramanujan_fraction'] = sum(rams)/len(rams)
Grand = nx.random_regular_graph(3, 60, seed=7); Gbad = two_cluster(30, 11)
out = {}
for name, G in [('petersen', P), ('random60', Grand), ('twocluster60', Gbad)]:
    w2, A2, B2, J2 = ihara_check(G)
    ev = np.sort(np.linalg.eigvalsh(A2))[::-1]
    mu2 = np.linalg.eigvals(B2)
    poles = [1/z for z in mu2 if abs(z) > 1e-9]
    out[name] = dict(bass_relerr=w2, eig=list(map(float, ev)), poles=[[p.real, p.imag] for p in poles],
                     lambda2=float(ev[1]), ramanujan=bool(max(abs(ev[1]), abs(ev[-1])) <= 2*np.sqrt(2) + 1e-9),
                     off_circle=int(sum(1 for p in poles if abs(abs(p) - 1/np.sqrt(2)) > 1e-6 and abs(abs(p) - 1) > 1e-6 and abs(abs(p) - 0.5) > 1e-6)))
R['G4_graphs'] = out

stamp('# ---------- B')
# ---------- B. Bartholdi bump-count parameter: a new dial or d in disguise? ----------
def bartholdi(G, t, u):
    n, m = G.number_of_nodes(), G.number_of_edges()
    A = nx.to_numpy_array(G); D = np.diag(A.sum(1)); B, J, E = hashimoto(G)
    lhs = np.linalg.det(np.eye(2*m) - u*(B + t*J))
    rhs = (1 - (1-t)**2*u*u)**(m - n) * np.linalg.det(np.eye(n) - u*A + (1-t)*(D - (1-t)*np.eye(n))*u*u)
    return abs(lhs - rhs)/abs(rhs)
R['B1_bartholdi_relerr'] = max(bartholdi(G, t, u) for G in [P, Grand] for t in [0.0, 0.3, 0.7, 1.4] for u in [0.11, 0.2+0.1j])
# per-eigenvalue factor 1 - lambda u + (1-t)(q+t) u^2 : t only moves d (d_eff = (1-t)(q+t)); with u -> u/(something)? record d_eff
R['B2_deff'] = {str(t): (1-t)*(2+t) for t in [0, 0.25, 0.5, 1.0]}

stamp('# ---------- W')
# ---------- W. Hasse/Weil: Frobenius as a dial point, degree form = norm form ----------
def ap(p, a=1, b=1):
    sq = [0]*p
    for x in range(p): sq[x*x % p] += 1
    N = 1 + sum(sq[(x**3 + a*x + b) % p] for x in range(p))
    return p + 1 - N
def primes(n):
    s = [True]*(n+1); s[0] = s[1] = False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = [False]*len(s[i*i::i])
    return [i for i in range(n+1) if s[i]]
good = [p for p in primes(2000) if p > 3 and (4 + 27) % p != 0]
aps = {p: ap(p) for p in good}
R['W1_all_hasse'] = all(aps[p]**2 < 4*p for p in good)
R['W1_n_primes'] = len(good)
# over F_{p^2} (p = 3 mod 4, i^2 = -1): count points directly, compare with deg(1 - F^2) = p^2 + 1 - (a_p^2 - 2p)
def count_fp2(p, a=1, b=1):
    els = [(x, y) for x in range(p) for y in range(p)]
    def mul(u, v): return ((u[0]*v[0] - u[1]*v[1]) % p, (u[0]*v[1] + u[1]*v[0]) % p)
    def add(u, v): return ((u[0]+v[0]) % p, (u[1]+v[1]) % p)
    sqc = {}
    for y in els:
        s = mul(y, y); sqc[s] = sqc.get(s, 0) + 1
    N = 1
    for x in els:
        rhs = add(add(mul(mul(x, x), x), mul((a % p, 0), x)), (b % p, 0))
        N += sqc.get(rhs, 0)
    return N
fp2 = {}
for p in [7, 11, 19, 23, 43, 47]:
    a = aps[p]
    fp2[str(p)] = dict(direct=count_fp2(p), from_dial_point=p*p + 1 - (a*a - 2*p))
R['W2_fp2'] = fp2
# degree form deg(m - nF) = m^2 - a m n + p n^2 is positive definite iff a^2 < 4p; min over small (m,n) != 0
R['W3_degree_form_min'] = {str(p): min(m*m - aps[p]*m*n + p*n*n for m in range(-40, 41) for n in range(-40, 41) if (m, n) != (0, 0)) for p in [5, 7, 101, 499, 1999]}
R['W3_ap'] = {str(p): aps[p] for p in [5, 7, 101, 499, 1999]}
R['W4_khat'] = [-aps[p]/np.sqrt(p) for p in good]
R['W4_primes'] = good

stamp('# ---------- L')
# ---------- L. Lee-Yang on the cube ----------
def ly_poly_uniform(a):
    n = a.shape[0]; c = np.zeros(n+1)
    for S in itertools.product([0, 1], repeat=n):
        w = 1.0
        for i in range(n):
            for j in range(n):
                if S[i] and not S[j]: w *= a[i, j]
        c[sum(S)] += w
    return c
def sym(n, lo, hi):
    a = np.ones((n, n))
    for i in range(n):
        for j in range(i+1, n): a[i, j] = a[j, i] = rng.uniform(lo, hi)
    return a
n = 9
ferro = sym(n, 0.05, 0.95); anti = sym(n, 1.2, 2.5)
rf = np.roots(ly_poly_uniform(ferro)[::-1]); ra = np.roots(ly_poly_uniform(anti)[::-1])
R['L1_ferro_root_abs_dev'] = float(max(abs(abs(rf) - 1)))
R['L1_anti_root_abs_dev'] = float(max(abs(abs(ra) - 1)))
R['L1_roots'] = dict(ferro=[[z.real, z.imag] for z in rf], anti=[[z.real, z.imag] for z in ra])
# two directions: 1 + a(z1+z2) + z1 z2 = 1 + a (z1 o_{1/a} z2); uniform field: z^2 + 2a z + 1, elliptic iff |a| <= 1
R['L2_two_dir_identity'] = all(abs((1 + a*(x+y) + x*y) - (1 + a*(x + y + x*y/a))) < 1e-12 for a in [0.3, 0.9, 1.7] for x, y in [(0.2, -1.3), (2.0, 0.7)])
# three directions: corner weights a^{|S|(3-|S|)} = 1, a^2, a^2, 1 are of product form (after z -> lam z) only if a^2 = 1
R['L3_three_dir_product_form_residual'] = {str(a): abs((a*a)**2 - (a*a)*1) for a in [0.3, 0.5, 0.9, 1.0]}
# coupled prime cube: z_p = p^(1/2 - s) on primes up to 19, random ferromagnetic couplings: zeros on Re s = 1/2
PR = [2, 3, 5, 7, 11, 13, 17, 19]
aP = sym(len(PR), 0.1, 0.9); aA = sym(len(PR), 1.3, 2.2)
subsets = list(itertools.product([0, 1], repeat=len(PR)))
def cube_weights(a):
    W = []
    for S in subsets:
        w = 1.0
        for i in range(len(PR)):
            for j in range(len(PR)):
                if S[i] and not S[j]: w *= a[i, j]
        W.append(w)
    return np.array(W)
logn = np.array([sum(np.log(PR[i]) for i in range(len(PR)) if S[i]) for S in subsets])
def cubeZ(W, shift):
    f = lambda s: np.sum(W * np.exp((shift - s) * logn))
    f.d = lambda s: np.sum(-logn * W * np.exp((shift - s) * logn))
    return f
def newton(f, z):
    for _ in range(60):
        dz = f(z)/f.d(z); z -= dz
        if abs(dz) < 1e-14: break
    return z
def find_zeros(f, sig0, sig1, t1, ns=50, nt=700):
    zs = []
    S = np.linspace(sig0, sig1, ns); T = np.linspace(0.5, t1, nt)
    V = np.array([[abs(f(complex(s, t))) for s in S] for t in T])
    for i in range(1, nt-1):
        for j in range(1, ns-1):
            if V[i, j] <= V[i-1:i+2, j-1:j+2].min():
                z = newton(f, complex(S[j], T[i]))
                if sig0 - .2 < z.real < sig1 + .2 and 0 < z.imag < t1 and abs(f(z)) < 1e-9 and all(abs(z - y) > 1e-6 for y in zs):
                    zs.append(z)
    return zs
stamp('zf = find_zero')
zf = find_zeros(cubeZ(cube_weights(aP), 0.5), -1.5, 2.5, 30)
za = find_zeros(cubeZ(cube_weights(aA), 0.5), -1.5, 2.5, 30)
z0 = find_zeros(cubeZ(cube_weights(aP), 0.0), -2.0, 2.0, 30)
R['L4_ferro_zeros'] = [[z.real, z.imag] for z in zf]
R['L4_ferro_max_dev_from_half'] = float(max(abs(z.real - 0.5) for z in zf))
R['L4_anti_zeros'] = [[z.real, z.imag] for z in za]
R['L4_anti_max_dev_from_half'] = float(max(abs(z.real - 0.5) for z in za))
R['L4_unscaled_zeros_max_dev_from_0'] = float(max(abs(z.real) for z in z0))

stamp('# ---------- H')
# ---------- H. Davenport-Heilbronn: functional equation, no Euler product, zero off the line ----------
mp.mp.dps = 30
def L5(s, chi):
    return mp.power(5, -s) * sum(chi[a] * mp.zeta(s, mp.mpf(a)/5) for a in range(1, 5))
chi = {1: 1, 2: 1j, 3: -1j, 4: -1}; chib = {a: complex(v).conjugate() for a, v in chi.items()}
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
fDH = lambda s: (1 - 1j*kap)/2 * L5(s, chi) + (1 + 1j*kap)/2 * L5(s, chib)
zDH = mp.findroot(fDH, mp.mpc(0.808, 85.699))
R['H1_DH_zero'] = [float(zDH.real), float(zDH.imag)]
# functional equation check: Lambda(s) = (pi/5)^(-(s+1)/2) Gamma((s+1)/2) f(s) is real on the critical line up to a constant phase
Lam = lambda s: mp.power(mp.pi/5, -(s+1)/2) * mp.gamma((s+1)/2) * fDH(s)
ph = [mp.arg(Lam(mp.mpc(0.5, t))) for t in [10, 20, 33.3]]
R['H2_DH_phase_on_line'] = [float(mp.cos(2*x)) for x in ph]
# a zero on the line it does have, near t ~ 5 (scan sign of Re(Lam e^{-i phase}))
stamp('# ---------- D')
# ---------- D. Dirichlet L(s, chi_-4): open neighbour, first zero on the line ----------
chi4 = lambda s: mp.power(4, -s) * (mp.zeta(s, mp.mpf(1)/4) - mp.zeta(s, mp.mpf(3)/4))
z4 = mp.findroot(chi4, mp.mpc(0.5, 6.02))
R['D1_L_chi4_first_zero'] = [float(z4.real), float(z4.imag)]
import os
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results.json'), 'w'), indent=1, default=float)
for k, v in R.items():
    if k not in ('L1_roots', 'G4_graphs', 'W4_khat', 'L4_ferro_zeros', 'L4_anti_zeros'): print(k, v)
print('G4', {k: (v['lambda2'], v['ramanujan'], v['off_circle'], v['bass_relerr']) for k, v in R['G4_graphs'].items()})
print('ferro zeros', len(zf), 'anti zeros', len(za))
