"""Round 6 robustness: (1) exact theta inversion (independent of the engine) at the p = 2 zeros of the curvature runs;
(2) square p = 2.5 and hexagon p = 3 recomputed with a changed small-x fit (extra power (p-1)/2, smaller x0);
(3) curvature runs at p = 3, 4, inf recomputed with the fit powers -1/2 and (p-1)/2 added."""
import json, numpy as np, mpmath as mp
from lattice_engine import *
mp.mp.dps = 30
res = {}

def make_exact_mp(M, dps=45, nodes=500, R=30):
    """Z_M(s) = sum' (m^T M m)^-s by Epstein's exact theta inversion, in mpmath at high precision (no code shared
    with the engine). Theta is tabulated once at Gauss-Legendre nodes in u = log x on [0, log X]."""
    with mp.workdps(dps):
        M = mp.matrix(M); Mi = M**-1; dt = mp.det(M)
        def qvals(A):
            out = []
            for a in range(-R, R + 1):
                for c in range(-R, R + 1):
                    if a or c: out.append(A[0, 0]*a*a + (A[0, 1] + A[1, 0])*a*c + A[1, 1]*c*c)
            return out
        gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
        tab = []
        for A in (M, Mi):
            q = qvals(A); X = mp.mpf(110)/(mp.pi*min(q)); L = mp.log(X)
            # degree 9 Gauss-Legendre rule (3*2^8 nodes) on [0, L]
            pts = gl.calc_nodes(9, mp.mp.prec)
            us = [(L/2)*(x + 1) for x, w in pts]; ws = [(L/2)*w for x, w in pts]
            th = [sum(mp.exp(-mp.pi*mp.exp(u)*qq) for qq in q) for u in us]
            tab.append((us, [w*t for w, t in zip(ws, th)]))
    def Z(s):
        with mp.workdps(dps):
            s = mp.mpc(s)
            (u1, f1), (u2, f2) = tab
            xi = -1/s - 1/(mp.sqrt(dt)*(1 - s)) + mp.fsum(f*mp.exp(s*u) for u, f in zip(u1, f1)) \
                 + mp.fsum(f*mp.exp((1 - s)*u) for u, f in zip(u2, f2))/mp.sqrt(dt)
            return xi*mp.power(mp.pi, s)/mp.gamma(s)
    return Z
d = json.load(open('lattice_zeros.json'))
import sympy
char = lambda D: [int(sympy.kronecker_symbol(D, n)) if n > 0 else 0 for n in range(abs(D))]
c20, c4, c5 = char(-20), char(-4), char(5)
mp.mp.dps = 20
genus20 = lambda s: mp.zeta(s)*mp.dirichlet(s, c20) + mp.dirichlet(s, c4)*mp.dirichlet(s, c5)
chk = {}
# disc -20: exact genus formula Z_{x^2+5y^2} = zeta L(chi_-20) + L(chi_-4) L(chi_5); refine each engine zero on it
rows = []
for z in d['disc-20_square']['2.0']['zeros']:
    r = mp.findroot(genus20, mp.mpc(*z))
    rows.append(dict(engine=z, genus_formula=[float(r.real), float(r.imag)], diff=float(abs(r - mp.mpc(*z)))))
chk['disc-20_square'] = rows
print('disc-20 max diff', max(r['diff'] for r in rows), 'off-line', [r['genus_formula'] for r in rows if abs(r['genus_formula'][0] - 0.5) > 1e-6], flush=True)
# disc -23: exact theta inversion for the p = 2 gauge's quadratic form
b1 = np.array([1.0, 0.0]); b2 = np.array([0.5, np.sqrt(23)/2]); V = [b1, b2, b2 - b1]
Nv = np.array([2*v/np.dot(v, v) for v in V] + [-2*v/np.dot(v, v) for v in V]); B = np.array([b1, b2])
G = B @ Nv.T; M23 = ((G @ G.T)/Nv.shape[0]).tolist()
Zx = make_exact_mp(M23)
rows = []
for z in d['disc-23_voronoi']['2.0']['zeros']:
    r = mp.findroot(Zx, mp.mpc(*z), tol=1e-24)
    rows.append(dict(engine=z, exact=[float(r.real), float(r.imag)], diff=float(abs(r - mp.mpc(*z)))))
chk['disc-23_voronoi'] = dict(matrix=M23, rows=rows)
print('disc-23 max diff', max(r['diff'] for r in rows), 'off-line', [r['exact'] for r in rows if abs(r['exact'][0] - 0.5) > 1e-6], flush=True)
res['p2_curvature_exact_check'] = chk

def zeros_of(B, N, p, **kw):
    zs, _ = find_zeros(LatticeZeta(B, N, p, **kw)); return np.array([[z.real, z.imag] for z in zs])

def match(z1, z2):
    """largest distance from each zero of z1 to the nearest of z2 (and counts)"""
    dd = [np.min(np.hypot(*(z2 - a).T)) for a in z1]
    return dict(n1=len(z1), n2=len(z2), max_dist=float(np.max(dd)))

rob = {}
for name, (B, N), p in [('square', square_setup(), 2.5), ('hexagon', hex_setup(), 3.0), ('square', square_setup(), 2.1)]:
    base = np.array(d[name][str(p)]['zeros'])
    alt = zeros_of(B, N, p, x0=1e-3, fit_pows=[0, 0.5, (p - 1)/2, 1, 1.5, 2, 2.5])
    rob[f'{name}_{p}'] = match(base, alt); print(name, p, rob[f'{name}_{p}'], flush=True)
res['fit_change'] = rob

curv = {}
B5 = np.array([[1.0, 0.0], [0.0, np.sqrt(5.0)]])
b1 = np.array([1.0, 0.0]); b2 = np.array([0.5, np.sqrt(23)/2]); V = [b1, b2, b2 - b1]
Nv = np.array([2*v/np.dot(v, v) for v in V] + [-2*v/np.dot(v, v) for v in V])
for name, B, N in [('disc-20_square', B5, square_setup()[1]), ('disc-23_voronoi', np.array([b1, b2]), Nv)]:
    out = {}
    for p in [3.0, 4.0, np.inf]:
        pw = [-0.5, 0, 0.5, 1, 1.5, 2] + ([] if np.isinf(p) else [(p - 1)/2])
        Z = LatticeZeta(B, N, p, x0=1e-3, fit_pows=pw, fit_span=20)
        zs, _ = find_zeros(Z); z = np.array([[q.real, q.imag] for q in zs])
        key = 'inf' if np.isinf(p) else str(p)
        old = np.array(d[name][key]['zeros'])
        out[key] = dict(zeros=z.tolist(), fit_resid=Z.fit_resid, vs_first_run=match(z, old) if len(z) else None)
        print(name, key, 'resid %.1e' % Z.fit_resid, out[key]['vs_first_run'], 're range', z[:, 0].min().round(3), z[:, 0].max().round(3), flush=True)
    curv[name] = out
res['curvature_refit'] = curv
json.dump(res, open('robustness.json', 'w'), indent=1)
print('done')
