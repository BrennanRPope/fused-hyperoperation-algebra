"""Round 4: the hypernatural logarithm ln_s(X) = E^(-s)(X - 1), E(x) = e^x - 1 (fractional-rungs/e-every-height).
Facts used: ln_s(1 + u) = u - (s/2) u^2 + O(u^3) for every height s.
(1) Laguerre at height s: normalise f at x0 (g = f/f(x0)); then -(ln_s g)''(x0) = (s f'^2 - f f'')/f^2.
    s = 1: Laguerre's expression; s = 1 - 1/n: Laguerre's inequality for real-rooted degree-n polynomials.
    At an extremum (f' = 0) every height gives -f''/f, so D = -8 f/f'' does not depend on s.
(2) The critical height h*(x) = f f''/f'^2 = 1 - 1/N_eff(x), N_eff = (sum a_k)^2 / sum a_k^2, a_k = 1/(x - z_k):
    with real zeros, h* < 1 everywhere (Cauchy-Schwarz); a complex pair close to the line pushes h* above 1.
    We measure h* along the line for Xi (zeros 1..4900, from the zeros, spot-checked from Xi itself).
(3) Products: ln_s(XY) = B_(1-s)(ln_s X, ln_s Y); only s = 1 turns the Euler product into a sum. -> hyperlog.json"""
import json, numpy as np, mpmath as mp, sympy as sp
res = {}
# (1) series of ln_s(1+u): E^-s(u) = u - (s/2)u^2 + ... ; check by the flow equation at small s numerically
# E^s(x) = x + s*(x^2/2) + ...: the generator of the flow at x is V(x) with V(x) = x^2/2 + x^3/12 - ... ; first terms suffice:
x, s, u = sp.symbols('x s u')
V = x**2/2 + x**3/12  # vector field of the flow E^s near 0 (first two terms of the iterative logarithm of e^x - 1)
# E^(-s)(u) to second order: u - s V(u) + O(s^2 u^3) -> coefficient of u^2 is -s/2 for every s (higher orders change u^3 only)
res['ln_s_second_order'] = '-s/2'
f0, f1, f2 = sp.symbols('f0 f1 f2', positive=True)
g1, g2 = f1/f0, f2/f0
lag_s = sp.simplify(-(1*g2 + (-s)*g1**2))           # -(ln_s g)'' with ln_s'(1)=1, ln_s''(1)=-s
res['laguerre_height_s'] = str(sp.simplify(lag_s*f0**2))
res['extremum_height_free'] = bool(sp.simplify(lag_s.subs(f1, 0) + f2/f0) == 0)
# degree-n polynomial: (n-1) f'^2 - n f f'' >= 0 is height s = 1 - 1/n; check on a random real-rooted polynomial
rng = np.random.default_rng(5); ok = True
for n in [3, 5, 8]:
    r = np.sort(rng.normal(size=n)); P = np.poly(r); P1 = np.polyder(P); P2 = np.polyder(P1)
    xs = np.linspace(r[0]-2, r[-1]+2, 4001)
    val = (n-1)*np.polyval(P1, xs)**2 - n*np.polyval(P, xs)*np.polyval(P2, xs)
    ok &= val.min() > -1e-9*np.abs(val).max()
res['laguerre_degree_n_holds'] = bool(ok)
# (2) critical height along the line for Xi, from the zeros
G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
w = np.linspace(G[-1], G[-1]*400, 200001); rho = np.log(w/(2*np.pi))/(2*np.pi)
# the tail beyond the last zero is smooth in x: tabulate it once and interpolate
xt = np.linspace(0, G[-1]*0.995, 3001)
T1t = np.array([np.trapezoid((1/(x0-w)+1/(x0+w))*rho, w) for x0 in xt])
T2t = np.array([np.trapezoid((1/(x0-w)**2+1/(x0+w)**2)*rho, w) for x0 in xt])
def sums(xv):
    d1 = xv[:, None]-G[None, :]; d2 = xv[:, None]+G[None, :]
    S1 = np.sum(1/d1 + 1/d2, axis=1); S2 = np.sum(1/d1**2 + 1/d2**2, axis=1)
    return S1+np.interp(xv, xt, T1t), S2+np.interp(xv, xt, T2t)
per = []
for i in range(0, 4900):
    a, b = G[i], G[i+1]
    xv = a + (b-a)*np.linspace(0.01, 0.99, 49)
    S1, S2 = sums(xv); h = 1 - S2/S1**2
    j = int(np.argmax(h))
    per.append((float(a), float(b), float(h[j]), float((xv[j]-a)/(b-a))))
per = np.array(per)
res['hstar_max'] = float(per[:, 2].max()); res['hstar_argmax_gamma'] = float(per[np.argmax(per[:, 2]), 0])
res['Neff_max'] = float(1/(1-per[:, 2].max()))
res['hstar_median_of_interval_max'] = float(np.median(per[:, 2]))
res['per_interval'] = per[::7].tolist()
blocks = [(lo, lo+700) for lo in range(0, 4900, 700)]
res['blocks'] = [dict(gamma_lo=float(G[lo]), gamma_hi=float(G[hi]), hmax=float(per[lo:hi, 2].max()), hmed=float(np.median(per[lo:hi, 2])),
                      Neff_med=float(np.median(1/(1-per[lo:hi, 2])))) for lo, hi in blocks]
# spot-check h* from Xi itself (Z derivatives and the gamma factor) at a few points
mp.mp.dps = 20
def lF(g, k):
    wv = mp.mpc(0.25, g/2)
    return 2*mp.re((-1)**(k-1)*mp.factorial(k-1)/mp.mpc(g, 0.5)**k) + mp.re((1j/2)**k*mp.psi(k-1, wv))
def hstar_xi(g):
    g = mp.mpf(g); Z = [mp.siegelz(g, derivative=k) for k in range(3)]; l1, l2 = lF(g, 1), lF(g, 2)
    X0 = Z[0]; X1 = Z[1]+l1*Z[0]; X2 = Z[2]+2*l1*Z[1]+(l2+l1**2)*Z[0]
    return float(X0*X2/X1**2)
spots = []
for i in [0, 10, 100, 1000, 3000]:
    a, b = G[i], G[i+1]; xv = np.array([a+0.3*(b-a)]); S1, S2 = sums(xv)
    spots.append((float(xv[0]), float(1-S2[0]/S1[0]**2), hstar_xi(xv[0])))
res['spot_checks'] = spots
# (3) the product law: ln_s(XY) vs B_(1-s); only s = 1 additive -- recorded as a statement (proved in e-every-height)
res['product_law'] = 'ln_s(XY) = B_(1-s)(ln_s X, ln_s Y); additive only at s = 1'
json.dump(res, open('hyperlog.json', 'w'))
print({k: v for k, v in res.items() if k not in ('per_interval',)})

# refine: the 25 intervals with the largest sampled h*, re-sampled at 2,001 points; plus the stretch 0 < x < gamma_1
top = np.argsort(per[:, 2])[-25:]
ref = []
for i in top:
    a, b = G[i], G[i+1]; xv = a + (b-a)*np.linspace(1e-4, 1-1e-4, 2001); S1, S2 = sums(xv); ref.append(float((1-S2/S1**2).max()))
xv = np.linspace(0.05, G[0]*(1-1e-4), 2001); S1, S2 = sums(xv); h0 = float((1-S2/S1**2).max())
res['hstar_max_refined'] = max(max(ref), h0); res['hstar_below_first_zero'] = h0
res['Neff_max_refined'] = 1/(1-res['hstar_max_refined'])
json.dump(res, open('hyperlog.json', 'w'))
print('refined', res['hstar_max_refined'], h0, res['Neff_max_refined'])
