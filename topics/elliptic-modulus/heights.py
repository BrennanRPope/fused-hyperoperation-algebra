"""Fractional heights (regular iteration) of the slope-chart elliptic generator G_m = sn dn / cn, on the real line.
G^s(x) = G^{+-n}( P_s( G^{-+n}(x) ) ), where the inner iterates run toward 0 along the real axis and P_s is the
truncated formal series of the height-s map (numseries.flow of the Julia field). Double precision; checked against
mpmath in polish.py."""
import numpy as np
from scipy.special import ellipj
import mpmath as mp
from numseries import slope, julia, flow
def G(x, m):
    sn, cn, dn, _ = ellipj(x, m); return sn*dn/cn
def Gp(x, m):
    sn, cn, dn, _ = ellipj(x, m); return (1 - 2*m*sn**2 + m*sn**4)/cn**2
def Ginv(t, m, x0=None):
    """Inverse on the principal sheet (-K, K): the law's own logarithm (height 1 of the hypernatural log)."""
    t = np.asarray(t, float)
    # L(t) = int_0^t dtau / sqrt(1 + 2 kappa tau^2 + tau^4); good start: Newton from arctan-ish guess
    x = np.arctan(t) if x0 is None else x0
    from scipy.special import ellipk
    K = ellipk(m)
    x = np.clip(x*(2*K/np.pi), -0.999*K, 0.999*K)
    for _ in range(60):
        dx = (G(x, m) - t)/Gp(x, m)
        x = np.clip(x - dx, -0.9999999*K, 0.9999999*K)
        if np.all(np.abs(dx) < 1e-15): break
    return x
_cache = {}
def series_coeffs(m, s, N=25):
    key = (round(m, 12), round(s, 12), N)
    if key not in _cache:
        mp.mp.dps = 40
        V = julia(slope(mp.mpf(m), N), N)
        _cache[key] = [float(c) for c in flow(V, mp.mpf(s), N)]
    return _cache[key]
def P(x, m, s):
    c = series_coeffs(m, s); return sum(cj*x**j for j, cj in enumerate(c))
def height(x, m, s, small=None):
    """G_m^s(x) for real x on the principal branch."""
    x = np.atleast_1d(np.asarray(x, float)).copy()
    if small is None: small = 0.04 if abs(1 - 2*m) > 0.15 else 0.3
    attract = (1 - 2*m) < 0          # cubic coefficient (1-2m)/3 < 0: forward iterates go to 0 on the real line
    out = np.empty_like(x)
    for i, xi in enumerate(x):
        if xi == 0: out[i] = 0; continue
        y, n = xi, 0
        while abs(y) > small:
            y = G(y, m) if attract else float(Ginv(y, m)); n += 1
            if n > 200000: raise RuntimeError('no convergence')
        z = P(y, m, s)
        for _ in range(n):
            z = float(Ginv(z, m)) if attract else G(z, m)
        out[i] = z
    return out
if __name__ == "__main__":
    for m in [0.2, 0.5, 0.8]:
        x = np.array([0.3, 0.6, 0.9])
        h = height(x, m, 0.5); hh = height(h, m, 0.5)
        print(m, 'half∘half - G =', np.max(np.abs(hh - G(x, m))), ' G^1 - G =', np.max(np.abs(height(x, m, 1.0) - G(x, m))))
