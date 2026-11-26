"""Round 8, item 3b-d: holes and twists.

(b) Twist on the transverse ladder.  For an orientation-preserving closed geodesic of length l the
    Selberg factor is prod_k (1 - e^{-(s+k) l}); its trace-formula weight is 1/(2 sinh(l/2)) =
    e^{-l/2} sum_k e^{-k l}.  A glide reflection flips the transverse direction: linearised return map
    eigenvalues -e^{+-l}, so |det(1-P)|^{1/2} = 2 cosh(l/2) = e^{l/2} sum_k (-1)^k e^{-kl}... i.e. the
    ladder gets the sign (-1)^k:  Z^-_l(s) = prod_k (1 - (-1)^k e^{-(s+k) l}).
    The orientation double cover sees the glide as one geodesic of length 2l, and
    Z^-_l(s) * Z^+_l(s) = Z_{2l}(s)   with Z^+_l = prod_k (1 + (-1)^k e^{-(s+k) l})   (rung by rung (1-x)(1+x)=1-x^2).
(c) Flat check of the same principle: Klein bottle heat trace minus half the torus trace equals
    (1/2) * (1/sqrt(4 pi tau)) * sum_j exp(-(j+1/2)^2/(4 tau)): each glide of length |j+1/2| enters with
    weight 1/|det(1 - (-1))| = 1/2.
(d) Quantum gravity: GUE (orientable JT sum) vs GOE (crosscaps included, Stanford-Witten) for zeta's spacings.
"""
import json, numpy as np, mpmath as mp
mp.mp.dps = 30
out = {}

# (b) ladder identity
rng = np.random.default_rng(1)
worst = 0
for _ in range(20):
    l = rng.uniform(0.5, 4); s = complex(rng.uniform(0.1, 2), rng.uniform(-20, 20))
    K = 200
    k = np.arange(K)
    zm = np.prod(1 - (-1.0)**k*np.exp(-(s+k)*l))
    zp = np.prod(1 + (-1.0)**k*np.exp(-(s+k)*l))
    z2 = np.prod(1 - np.exp(-(s+k)*2*l))
    worst = max(worst, abs(zm*zp - z2)/abs(z2))
print(f"(b) Z^-_l Z^+_l = Z_2l: worst relative error over 20 random (l, s): {worst:.1e}")
out['ladder_identity_err'] = worst
# weights: log-derivative expansion check  d/ds log Z^-_l = sum_m l e^{-s m l} w_m, w_m = 1/(1-e^{-ml}) (m even), 1/(1+e^{-ml}) (m odd)
l = 1.3; s = mp.mpc(0.7, 3.1)
lhs = mp.diff(lambda x: mp.log(mp.nprod(lambda k: 1 - (-1)**int(k)*mp.e**(-(x+k)*l), [0, 300])), s)
rhs = mp.nsum(lambda m: l*mp.e**(-s*m*l)*(1/(1-mp.e**(-m*l)) if int(m) % 2 == 0 else 1/(1+mp.e**(-m*l))), [1, mp.inf])
print(f"    log-derivative: {mp.nstr(lhs, 12)} vs orbit sum with sinh/cosh weights {mp.nstr(rhs, 12)}")
out['weight_check'] = float(abs(lhs - rhs))

# (c) flat Klein bottle heat trace
def theta_T(tau, R=60):
    m = np.arange(-R, R+1); e = np.exp(-4*np.pi**2*m**2*tau)
    return e.sum()**2
def theta_K(tau, R=60):
    m = np.arange(-R, R+1); e1 = np.exp(-4*np.pi**2*m**2*tau)
    n = np.arange(1, R+1); e2 = np.exp(-4*np.pi**2*n**2*tau)
    return e1.sum()*e2.sum() + np.exp(-4*np.pi**2*(m[m % 2 == 0])**2*tau).sum()
errs = []
for tau in [0.005, 0.02, 0.1, 0.5]:
    lhs = theta_K(tau) - 0.5*theta_T(tau)
    j = np.arange(-60, 60)
    rhs = 0.5/np.sqrt(4*np.pi*tau)*np.exp(-(j+0.5)**2/(4*tau)).sum()
    errs.append([tau, float(lhs), float(rhs)])
    print(f"(c) tau={tau}: Klein - torus/2 = {lhs:.10f}; glide sum with weight 1/2 = {rhs:.10f}")
out['flat_glide'] = errs

# (d) GOE vs GUE
S = lambda d: np.sinc(d)
def E_pm(s, sign, m=40):
    x, w = np.polynomial.legendre.leggauss(m)
    res = []
    for L in np.atleast_1d(s):
        if L == 0: res.append(1.0); continue
        a = L/2; y = (x+1)*a/2; ww = w*a/2
        K = S(y[:, None]-y[None, :]) + sign*S(y[:, None]+y[None, :])
        res.append(np.linalg.det(np.eye(m) - np.sqrt(ww)[:, None]*K*np.sqrt(ww)[None, :]))
    return np.array(res)
sg = np.linspace(0, 4, 401); h = sg[1]-sg[0]
Ep, Em = E_pm(sg, +1), E_pm(sg, -1)
import sys; sys.path.insert(0, '../round7')
from rmt import gue_E
print(f"(d) check E_+ E_- = E_GUE: max diff {np.abs(Ep*Em - gue_E(sg)).max():.1e}")
p_goe = np.gradient(np.gradient(Ep, h), h); p_gue = np.gradient(np.gradient(Ep*Em, h), h)
for nm, p in (('GOE', p_goe), ('GUE', p_gue)):
    print(f"    {nm}: integral {np.trapezoid(p, sg):.4f}, mean {np.trapezoid(sg*p, sg):.4f}")
d7 = json.load(open('../round7/rmt.json'))
def ks(sample, pdf):
    cdf = np.cumsum(pdf)*h; cdf /= cdf[-1]
    xs = np.sort(sample); emp = np.arange(1, len(xs)+1)/len(xs)
    return float(np.max(np.abs(emp - np.interp(xs, sg, cdf))))
sp = np.array(d7['zeta']['spacings'])
kg, ko = ks(sp, p_gue), ks(sp, p_goe)
small = float(np.mean(sp < 0.25))
pg, po = [float(np.trapezoid(p[sg <= 0.25], sg[sg <= 0.25])) for p in (p_gue, p_goe)]
print(f"    zeta spacings (n={len(sp)}, heights 10000-11600): KS to GUE {kg:.3f}, KS to GOE {ko:.3f}")
print(f"    fraction of spacings < 0.25: zeta {small:.4f}, GUE {pg:.4f}, GOE {po:.4f}")
out.update({'sgrid': sg.tolist(), 'p_goe': p_goe.tolist(), 'p_gue': p_gue.tolist(), 'zeta_spacings': sp.tolist(),
            'ks_gue': kg, 'ks_goe': ko, 'small': [small, pg, po]})

# topology table: integer picks -> chi, area (Gauss-Bonnet, curvature -1), max small eigenvalues (BMM 2014: < -chi)
tab = []
for g in range(2, 6): tab.append({'surface': f'orientable genus {g}', 'chi': 2-2*g, 'area_over_pi': -2*(2-2*g), 'max_small': -(2-2*g)})
for k in range(3, 9): tab.append({'surface': f'{k} crosscaps', 'chi': 2-k, 'area_over_pi': -2*(2-k), 'max_small': -(2-k)})
out['table'] = tab
for r in tab: print("   ", r)
json.dump(out, open('topology.json', 'w'))
print('saved')
