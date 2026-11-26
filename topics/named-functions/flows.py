"""Height flows for every shape's composition law (generic regular iteration).

A composition law a * b with a generator G (Aczel): a * b = G(G^-1(a) + G^-1(b)), G(0) = 0, G'(0) = 1.
The rung ladder of that law climbs by G: rung s operation = (G^s)_+ (+).  Rung 0 is +, rung 1 is the law itself.
The hypernatural logarithm of the law at height s is G^-s (read on the lifted variable for the cube).

G^s for real s needs a smooth choice: regular iteration at the parabolic fixed point 0
(Szekeres 1958, Ecalle).  Classical method: walk toward 0 with n steps, apply the formal
series of G^s there (its coefficients are polynomials in s), walk back out.
Polish (cheap): run at n and 2n, keep the 2n value, report |difference| as the error estimate.
"""
import mpmath as mp
mp.mp.dps = 40
ORDER = 16

# ---------- truncated power series (lists, index = power, c[0] = 0) ----------
def s_mul(a, b, N=ORDER):
    out = [mp.mpf(0)]*(N+1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(0, N+1-i):
            if b[j] != 0: out[i+j] += x*b[j]
    return out

def s_compose(f, g, N=ORDER):
    """f(g(x)), g[0] = 0"""
    out = [mp.mpf(0)]*(N+1)
    gp = [mp.mpf(0)]*(N+1); gp[0] = mp.mpf(1)
    for i in range(1, N+1):
        gp = s_mul(gp, g, N)
        if f[i] != 0:
            out = [o + f[i]*c for o, c in zip(out, gp)]
    return out

def s_revert(f, N=ORDER):
    """compositional inverse of f = x + ..., by Newton-free fixed point on coefficients"""
    g = [mp.mpf(0)]*(N+1); g[1] = mp.mpf(1)
    for n in range(2, N+1):
        c = s_compose(f, g, N)
        g[n] -= c[n]
    return g

class Flow:
    """Regular iteration of a map f tangent to the identity at 0."""
    def __init__(self, f, finv, series, name=""):
        self.f, self.finv, self.name = f, finv, name
        self.ser = series
        # integer iterates 0..ORDER (series), coefficient i is a polynomial of degree <= i-1 in s
        its = [[mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)]*(ORDER-1)]
        for _ in range(ORDER):
            its.append(s_compose(series, its[-1]))
        self.its = its
        self._cache = {}
        # direction: sign of first nonzero correction term a_m x^m
        m = next(i for i in range(2, ORDER+1) if abs(series[i]) > mp.mpf(10)**(-30))
        self.m, self.am = m, series[m]

    def coeffs(self, s):
        key = mp.nstr(s, 30)
        if key in self._cache: return self._cache[key]
        cs = [mp.mpf(0)]*(ORDER+1)
        for i in range(1, ORDER+1):
            pts = range(0, i)
            val = mp.mpf(0)
            for j in pts:
                w = mp.mpf(1)
                for q in pts:
                    if q != j: w *= (s - q)/(j - q)
                val += self.its[j][i]*w
            cs[i] = val
        self._cache[key] = cs
        return cs

    def P(self, s, y):
        cs = self.coeffs(s)
        out = mp.mpf(0)
        for i in range(ORDER, 0, -1):
            out = (out + cs[i])*y
        return out

    def repelling(self, x):
        return self.am * x**(self.m - 1) > 0   # f moves x away from 0

    def raw(self, s, x, n):
        s = mp.mpf(s); x = mp.mpf(x)
        if x == 0: return mp.mpf(0)
        if self.repelling(x):
            y = x
            for _ in range(n): y = self.finv(y)
            y = self.P(s, y)
            for _ in range(n): y = self.f(y)
        else:
            y = x
            for _ in range(n): y = self.f(y)
            y = self.P(s, y)
            for _ in range(n): y = self.finv(y)
        return y

    def __call__(self, s, x, n=400, polish=True):
        if not polish: return self.raw(s, x, n)
        a = self.raw(s, x, n); b = self.raw(s, x, 2*n)
        self.last_err = abs(a - b)
        return b

# ---------- the maps ----------
def E_k(k):
    """exponential of the circle operation: (e^{kx}-1)/k"""
    k = mp.mpf(k)
    f = lambda x: mp.expm1(k*x)/k
    finv = lambda y: mp.log1p(k*y)/k
    ser = [mp.mpf(0)] + [k**(i-1)/mp.factorial(i) for i in range(1, ORDER+1)]
    return Flow(f, finv, ser, f"E_{k}")

def T_kappa(kappa, t):
    kappa = mp.mpf(kappa)
    if kappa > 0: r = mp.sqrt(kappa); return mp.tan(r*t)/r
    if kappa < 0: r = mp.sqrt(-kappa); return mp.tanh(r*t)/r
    return mp.mpf(t)

def T_kappa_inv(kappa, x):
    kappa = mp.mpf(kappa)
    if kappa > 0: r = mp.sqrt(kappa); return mp.atan(r*x)/r
    if kappa < 0: r = mp.sqrt(-kappa); return mp.atanh(r*x)/r
    return mp.mpf(x)

def E_kd(k, d):
    """exponential of the two-dial law (a+b+kab)/(1-dab):  E = sigma(T_kappa(t)),
    sigma(x) = x/(1-(k/2)x),  kappa = d - k^2/4 = -Delta/4 (the curvature reading)."""
    k = mp.mpf(k); d = mp.mpf(d); kap = d - k**2/4
    f = lambda t: (lambda T: T/(1 - (k/2)*T))(T_kappa(kap, t))
    finv = lambda x: T_kappa_inv(kap, x/(1 + (k/2)*x))
    # series of L' = 1/(1+kt+dt^2)  -> L -> revert
    inv = [mp.mpf(0)]*(ORDER+1); inv[0] = mp.mpf(1)
    for n in range(1, ORDER+1):
        inv[n] = -(k*inv[n-1] + (d*inv[n-2] if n >= 2 else 0))
    L = [mp.mpf(0)] + [inv[n-1]/n for n in range(1, ORDER+1)]
    ser = s_revert(L)
    fl = Flow(f, finv, ser, f"E_({k},{d})")
    fl.k, fl.d, fl.kappa = k, d, kap
    return fl

def law_kd(k, d, a, b):
    return (a + b + k*a*b)/(1 - d*a*b)

def L_kd_numeric(k, d, x):
    return mp.quad(lambda t: 1/(1 + k*t + d*t*t), [0, x])
