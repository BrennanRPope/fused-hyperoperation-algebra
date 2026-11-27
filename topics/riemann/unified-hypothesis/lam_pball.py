"""De Bruijn-Newman constant of the p-ball Lee-Yang kernel w(u) = (1-|u|^p)^a on [-1,1].
H_lam(z) = int_0^1 e^{lam u^2} (1-u^p)^a cos(z u) du.  Lambda = inf{lam : H_lam has only real zeros}.
Complex zeros are counted with the argument principle on the rectangle [0,S] x [-Y,Y]."""
import numpy as np, json, sys, time
from scipy.special import roots_jacobi
from math import pi, sin, log

def nodes(p, a, n=1200):
    x, w = roots_jacobi(n, a, 0.0)          # weight (1-x)^a on [-1,1]
    u = (x + 1) / 2
    wt = w * 0.5 ** (a + 1)                  # (1-u)^a du
    g = np.where(u < 1, ((1 - u ** p) / np.where(u < 1, 1 - u, 1)) ** a, p ** a)
    return u, wt * g

class H:
    def __init__(self, p, a, n=1200):
        self.u, self.w = nodes(p, a, n)
    def __call__(self, lam, z, chunk=4000):
        z = np.asarray(z, complex)
        ww = self.w * np.exp(lam * self.u ** 2)
        out = np.empty(z.shape, complex)
        zf = z.ravel(); of = out.ravel()
        for i in range(0, len(zf), chunk):
            of[i:i + chunk] = np.cos(np.outer(zf[i:i + chunk], self.u)) @ ww
        return out

def real_zeros(h, lam, S, ds=0.01):
    s = np.arange(0.005, S, ds)
    v = h(lam, s).real
    idx = np.where(np.sign(v[1:]) != np.sign(v[:-1]))[0]
    return s[idx] + ds / 2, s, v

def count(h, lam, S=150.0, Y=4.0, step=0.02):
    rz, s, v = real_zeros(h, lam, S)
    # right edge: midway between the last two real zeros below S (or at S)
    R = 0.5 * (rz[-1] + rz[-2]) if len(rz) >= 2 else S
    nreal = int(np.sum(rz < R))
    nb = int(R / step) + 1; nv = int(2 * Y / step) + 1
    bot = np.linspace(0, R, nb) - 1j * Y
    right = R + 1j * np.linspace(-Y, Y, nv)
    top = np.linspace(R, 0, nb) + 1j * Y
    left = 1j * np.linspace(Y, -Y, nv)
    path = np.concatenate([bot, right[1:], top[1:], left[1:]])
    f = h(lam, path)
    ph = np.unwrap(np.angle(f))
    ntot = int(round((ph[-1] - ph[0]) / (2 * pi)))
    return ntot, nreal, R

def has_complex(h, lam, n0=None, S=150.0, **kw):
    ntot, nreal, R = count(h, lam, S=S, **kw)
    bad = ntot > nreal or nreal < 2 or R < 0.8 * S or (n0 is not None and nreal < n0 - 1)
    return bad, ntot, nreal

def Lambda(p, a, lo=-12.0, hi=0.0, tol=0.02, **kw):
    h = H(p, a)
    kw['n0'] = count(h, 0.0, S=kw.get('S', 150.0))[1]
    c_hi = has_complex(h, hi, **kw)
    if c_hi[0]:
        return dict(status='complex at lam=0', ntot=c_hi[1], nreal=c_hi[2])
    c_lo = has_complex(h, lo, **kw)
    if not c_lo[0]:
        return dict(status=f'real down to {lo}', Lambda_upper=lo)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if has_complex(h, mid, **kw)[0]:
            lo = mid
        else:
            hi = mid
    return dict(status='found', Lambda=0.5 * (lo + hi), bracket=[lo, hi])

if __name__ == '__test__':
    import mpmath as mp
    p, a = 2.5, 1.5
    h = H(p, a)
    for z in [3.0, 40.0, 120.0, 30 + 3j]:
        ex = mp.quad(lambda u: mp.e ** (-2 * u * u) * (1 - u ** p) ** a * mp.cos(z * u), [0, 0.5, 1])
        print(z, h(-2.0, [z])[0], complex(ex))

if __name__ == '__main__':
    grid = []
    for p in [1.5, 2, 2.5, 3, 3.5, 4, 5, 6]:
        for a in [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]:
            if a <= p:
                grid.append((p, a))
    grid += [(2, 3), (2, 4), (2, 6), (4, 5), (4, 6), (6, 8), (1, 0.5), (1, 1), (1.2, 0.6), (1.2, 1.2)]
    import os
    out = json.load(open('lam_pball.json')) if os.path.exists('lam_pball.json') else {}
    t0 = time.time()
    for p, a in grid:
        if f"{p}|{a}" in out:
            continue
        r = Lambda(p, a, lo=-14.0, tol=0.05)
        B = abs(sin(pi * p / 2)) / p ** (p - 1)
        r['far_threshold_at_a_eq_p'] = (log(B) if (a == p and B > 1e-12) else None)
        out[f"{p}|{a}"] = r
        print(p, a, r, round(time.time() - t0), flush=True)
        json.dump(out, open('lam_pball.json', 'w'))
