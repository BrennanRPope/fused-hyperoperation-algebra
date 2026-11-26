"""Round 5, pair 3: the height dial s across the roundness (lp) kernels.

Kernel e^(-|u|^p) on the roundness dial; its cosine transform
  F_p(x) = int e^(-|u|^p) cos(x u) du = (2/p) sum_n (-1)^n Gamma((2n+1)/p) x^(2n)/(2n)!
(entire for p > 1). Riemann's Phi decays like exp(-pi e^(4u)), the far end of this dial.

Height s on Laguerre's expression (round 4): -(ln_s g)'' ~ s f'^2 - f f''. The critical height
h*(x) = f f''/f'^2 = 1 - 1/N_eff. With only real zeros h* < 1 everywhere; a non-real pair pushes it
above 1. So the height dial at s = 1 is a detector for zeros leaving the line, and the largest h* in
each gap says which Laguerre degree the kernel passes there (s = 1 - 1/n).
"""
import json, numpy as np, mpmath as mp

res = {}
def series_fd(p, x, nmax=None):
    """f, f', f'' of F_p at real x with adaptive precision"""
    x = mp.mpf(x)
    # size of the largest term (log10) to set the precision
    lg = lambda n: (mp.loggamma((2*n + 1)/mp.mpf(p)) + 2*n*mp.log(abs(x) + 1e-30) - mp.loggamma(2*n + 1))/mp.log(10)
    n = 0; best = -1e9
    while True:
        val = lg(n)
        best = max(best, val)
        if n > 5 and val < best - 40 and val < -40: break
        n += 1
    N = n + 5
    with mp.workdps(int(best) + 40):
        f = f1 = f2 = mp.mpf(0)
        for k in range(N):
            c = (-1)**k*mp.gamma((2*k + 1)/mp.mpf(p))/mp.factorial(2*k)
            f += c*x**(2*k)
            if k >= 1: f1 += c*2*k*x**(2*k - 1)
            if k >= 1: f2 += c*2*k*(2*k - 1)*x**(2*k - 2)
        return [float(2/mp.mpf(p)*v) for v in (f, f1, f2)]

pvals = [2.5, 3.0, 4.0, 6.0, 8.0]
X = 22.0
xs = np.arange(0.02, X, 0.01)
out = {}
for p in pvals:
    F = np.array([series_fd(p, x) for x in xs])
    f, f1, f2 = F[:, 0], F[:, 1], F[:, 2]
    zi = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    zeros = [float(xs[i] - f[i]*(xs[i+1]-xs[i])/(f[i+1]-f[i])) for i in zi]
    with np.errstate(divide='ignore', invalid='ignore'):
        hstar = f*f2/(f1*f1)
    # largest h* in each gap (and below the first zero)
    edges = [0.0] + zeros + [X]
    gaps = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (xs > a) & (xs < b) & np.isfinite(hstar)
        if m.sum() > 3:
            j = np.argmax(np.where(m, hstar, -np.inf))
            gaps.append(dict(a=a, b=b, hmax=float(hstar[j]), at=float(xs[j])))
    lag_fail = [g for g in gaps if g['hmax'] > 1]
    out[str(p)] = dict(zeros=zeros, gaps=gaps, n_zeros=len(zeros),
                       hmax_overall=float(max(g['hmax'] for g in gaps)),
                       laguerre_fails_in=len(lag_fail),
                       xs=xs[::5].tolist(), h=np.clip(hstar, -5, 5)[::5].tolist())
    print(p, len(zeros), out[str(p)]['hmax_overall'], len(lag_fail), flush=True)

# non-real zeros where the detector fires (p = 2.5, 3): Newton from the h* peak
mp.mp.dps = 40
def Fc(p, z, N=400):
    z = mp.mpc(z)
    return 2/mp.mpf(p)*mp.fsum((-1)**k*mp.gamma((2*k+1)/mp.mpf(p))/mp.factorial(2*k)*z**(2*k) for k in range(N))
cz = {}
for p in [2.5, 3.0]:
    found = []
    for g in out[str(p)]['gaps']:
        if g['hmax'] > 1:
            for y0 in [0.3, 0.8, 1.5]:
                try:
                    with mp.workdps(120):
                        r = mp.findroot(lambda z: Fc(p, z, 600), mp.mpc(g['at'], y0))
                    rr = (round(float(r.real), 6), round(abs(float(r.imag)), 6))
                    if rr[1] > 1e-6 and all(abs(rr[0]-a)+abs(rr[1]-b) > 1e-4 for a, b in found) and rr[0] < X + 2:
                        found.append(rr)
                except Exception:
                    pass
    cz[str(p)] = sorted(found)
res['nonreal_zeros'] = cz
# Gaussian (p = 2): F = sqrt(pi) e^(-x^2/4), h* = 1 - 2/x^2 -> 1, no zeros
res['gaussian_hstar'] = '1 - 2/x^2'
res['kernels'] = out
res['xi_round4_hmax'] = 0.7048
json.dump(res, open('pair3.json', 'w'), indent=1)
print(cz)
