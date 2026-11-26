"""Round 9: the Hilbert-Polya string probe. Shared machinery.

A "string" here is a semiclassical spectral model: a smooth counting law (from the
classical dynamics, with an extension angle) plus a sum over closed orbits of length
m*log p with weights.  Its levels E_n solve the quantisation condition

    F(E) = Nbar(E) + S_model(E) = n - 1/2,      Nbar(E) = theta(E)/pi + 1,

which at the truth (all orbits, weights p^{-m/2}/m) is the Franca-LeClair equation
theta(t_n) + arg zeta(1/2 + 0 + i t_n) = (n - 3/2) pi.

Dials (only the seven model dials): k homothety, t odds, s height (read as the line
sigma on which the orbit sum is taken), d curvature, q order, r ratio, p roundness.
Settings: everything else the string can adjust (cutoff X, taper, extension angle,
amplitude, Maslov phase, growth law, free weights/frequencies).
"""
import numpy as np
from math import pi, log

TWO_PI = 2 * pi


# ---------- smooth part ----------
def theta(t):
    t = np.asarray(t, float)
    return t / 2 * np.log(t / TWO_PI) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def dtheta(t):
    t = np.asarray(t, float)
    return 0.5 * np.log(t / TWO_PI) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)


def gram(m):
    """Gram points g_m: theta(g_m) = m*pi (m may be fractional, m >= -1)."""
    m = np.asarray(m, float)
    t = TWO_PI * np.exp(1 + np.real(_lambertw((m + 0.125) / np.e)))
    for _ in range(30):
        t = t - (theta(t) - m * pi) / dtheta(t)
    return t


def _lambertw(x):
    from scipy.special import lambertw
    return lambertw(x)


def mean_spacing(t):
    return TWO_PI / np.log(np.asarray(t, float) / TWO_PI)


# ---------- orbits ----------
def primes_upto(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def factor_exponents(n, P):
    ex = []
    for p in P:
        if p * p > n:
            break
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            ex.append(e)
    if n > 1:
        ex.append(1)
    return ex


def orbit_table(X, q=0.0):
    """Orbits n <= X: returns (n, log n, coefficient c(n), prime-power info).
    c(n) at q = 0: 1/m for n = p^m (log zeta).  At q = 1: (Omega-1)!/prod e_i! (log 1/(1-P)).
    In between: log-linear path (1-q) a(n) + q b(n)  [not the exact q-multinomial dial]."""
    from math import factorial
    P = primes_upto(max(int(X), 2))
    if q == 0.0:
        ns, cs, ps, ms = [], [], [], []
        for p in P:
            m, pm = 1, int(p)
            while pm <= X:
                ns.append(pm); cs.append(1.0 / m); ps.append(p); ms.append(m)
                m += 1; pm *= int(p)
        n = np.array(ns, float)
        o = np.argsort(n)
        return dict(n=n[o], logn=np.log(n[o]), c=np.array(cs)[o], pp=np.ones(len(n), bool),
                    p=np.array(ps, float)[o], m=np.array(ms, float)[o])
    rows = []
    for n in range(2, int(X) + 1):
        ex = factor_exponents(n, P)
        Om = sum(ex)
        a = (1.0 / ex[0]) if len(ex) == 1 else 0.0
        if q == 0.0 and a == 0.0:
            continue
        b = factorial(Om - 1) / np.prod([factorial(e) for e in ex])
        c = (1 - q) * a + q * b
        if c == 0:
            continue
        pp = (len(ex) == 1)
        p = 0
        if pp:
            # recover p
            for pr in P:
                if n % pr == 0:
                    p = pr
                    break
        rows.append((n, log(n), c, pp, p, ex[0] if pp else Om))
    n = np.array([r[0] for r in rows], float)
    return dict(n=n, logn=np.log(n), c=np.array([r[2] for r in rows]),
                pp=np.array([r[3] for r in rows]), p=np.array([r[4] for r in rows], float),
                m=np.array([r[5] for r in rows], float))


DEFAULT = dict(k=1.0, t=1.0, s=0.5, d=0.0, q=0.0, r=1.0, p=None,
               X=200.0, taper='fejer', alpha=0.0, A=1.0, phi=0.0, growth='xp',
               growth_c=None)


def weights(orb, P):
    """Orbit weights w(n) including all dials and settings except the phase."""
    k, t, s, d = P['k'], P['t'], P['s'], P['d']
    m = orb['m']
    c = orb['c'].copy()
    pp = orb['pp']
    # curvature d: local factor 1 - u + d u^2, roots alpha+beta = 1, alpha*beta = d
    if d != 0.0:
        disc = complex(1 - 4 * d) ** 0.5
        al, be = (1 + disc) / 2, (1 - disc) / 2
        c = np.where(pp, np.real(al ** m + be ** m) / m, c)
    c = np.where(pp, c * k ** (m - 1) * t ** m, c)
    w = c * orb['n'] ** (-s)
    # taper
    X = P['X']
    if P['taper'] == 'fejer':
        w = w * np.clip(1 - orb['logn'] / log(X), 0, None)
    elif P['taper'] == 'gauss':
        w = w * np.exp(-(orb['logn'] / log(X)) ** 2 * 2)
    elif P['taper'] == 'sharp':
        w = w * (orb['n'] <= X)
    return P['A'] * w


class String:
    def __init__(self, **kw):
        self.P = dict(DEFAULT)
        self.P.update(kw)
        Xb = self.P['X'] if self.P['taper'] != 'gauss' else self.P['X'] ** 1.6
        self.orb = orbit_table(Xb, self.P['q'])
        self.w = weights(self.orb, self.P)
        keep = np.abs(self.w) > 1e-12
        self.w = self.w[keep]
        self.L = self.orb['logn'][keep] * self.P['r']
        self.free = None  # optional (amps, freqs, phases) replacing the prime sum

    # smooth count
    def Nbar(self, E):
        g = self.P['growth']
        if g == 'xp':
            base = theta(E) / pi + 1
        elif g == 'xp_lead':      # Berry-Keating leading term only
            base = E / TWO_PI * (np.log(E / TWO_PI) - 1) + 7 / 8
        elif g == 'oscillator':   # N ~ c E  (harmonic oscillator, 1D)
            base = self.P['growth_c'] * E
        elif g == 'box':          # N ~ c E^{1/2}  (particle in a box)
            base = self.P['growth_c'] * np.sqrt(E)
        elif g == 'billiard':     # N ~ c E  (2D billiard, Weyl), same as oscillator in form
            base = self.P['growth_c'] * E
        elif g == 'gammaC':       # archimedean factor of a complex place: (2pi)^-s Gamma(s)
            from scipy.special import loggamma
            base = (np.imag(loggamma(0.5 + 1j * np.asarray(E))) - np.asarray(E) * np.log(TWO_PI)) / pi + 1
        elif g == 'union':
            raise NotImplementedError
        elif g == 'nopi':         # Gamma(s/2) without the pi^{-s/2}
            base = (theta(E) + np.asarray(E) / 2 * np.log(pi)) / pi + 1
        return base + self.P['alpha'] / TWO_PI

    def dNbar(self, E):
        g = self.P['growth']
        if g in ('xp',):
            return dtheta(E) / pi
        if g == 'xp_lead':
            return np.log(np.asarray(E) / TWO_PI) / TWO_PI
        if g in ('oscillator', 'billiard'):
            return self.P['growth_c'] + 0 * np.asarray(E)
        if g == 'box':
            return 0.5 * self.P['growth_c'] / np.sqrt(E)

    def S(self, E, chunk=400):
        E = np.atleast_1d(np.asarray(E, float))
        out = np.empty_like(E)
        if self.free is not None:
            a, f, ph = self.free
            for i in range(0, len(E), chunk):
                e = E[i:i + chunk, None]
                out[i:i + chunk] = (a * np.sin(e * f + ph)).sum(1)
            return out
        for i in range(0, len(E), chunk):
            e = E[i:i + chunk, None]
            out[i:i + chunk] = -(self.w * np.sin(e * self.L + self.P['phi'])).sum(1) / pi
        return out

    def F(self, E):
        return self.Nbar(E) + self.S(E)

    def levels(self, nmax, nmin=1, iters=48):
        """E_n: a crossing of F = n - 1/2 found by bisection in [g_{n-3.5}, g_{n+0.5}]
        (for the xp growth law); other growth laws bracket from Nbar."""
        n = np.arange(nmin, nmax + 1, dtype=float)
        if self.P['growth'] in ('xp', 'xp_lead', 'nopi_off'):
            lo = gram(np.maximum(n - 3.5, -0.9))
            hi = gram(n + 0.5)
        else:
            # invert Nbar roughly with a wide bracket
            lo = np.full_like(n, 1.0)
            hi = np.full_like(n, 1.0)
            Eg = np.linspace(1, 1e5, 200001)
            Nb = self.Nbar(Eg)
            lo = np.interp(n - 3.0, Nb, Eg)
            hi = np.interp(n + 2.0, Nb, Eg)
            lo = np.maximum(lo, 1.0)
        target = n - 0.5
        flo = self.F(lo) - target
        fhi = self.F(hi) - target
        ok = (flo < 0) & (fhi > 0)
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            fm = self.F(mid) - target
            neg = fm < 0
            lo = np.where(neg, mid, lo)
            hi = np.where(neg, hi, mid)
        E = 0.5 * (lo + hi)
        E[~ok] = np.nan
        return E

    def crossings(self, Emin, Emax, step=0.02):
        """Count up/down crossings of each half-integer level (well-posedness check)."""
        Eg = np.arange(Emin, Emax, step)
        Fv = self.F(Eg)
        lev = np.floor(Fv + 0.5)  # index n with n-1/2 <= F < n+1/2
        d = np.diff(lev)
        return Eg, Fv, d


# ---------- truth from the zeros ----------
def true_S(E, gam):
    E = np.asarray(E, float)
    N = np.searchsorted(gam, E, side='right')
    return N - theta(E) / pi - 1


def unfold(x):
    return theta(np.asarray(x)) / pi + 1


# ---------- signatures ----------
GUE_R2 = lambda x: 1 - (np.sinc(x)) ** 2


def spacings(levels):
    u = unfold(np.sort(levels))
    return np.diff(u)


def number_variance(levels, Ls, nwin=4000, rng=None):
    rng = np.random.default_rng(0) if rng is None else rng
    u = np.sort(unfold(levels))
    out = []
    for L in Ls:
        starts = rng.uniform(u[0], u[-1] - L, nwin)
        cnt = np.searchsorted(u, starts + L) - np.searchsorted(u, starts)
        out.append(cnt.var())
    return np.array(out)


def pair_corr(levels, xmax=3.0, nb=30):
    u = np.sort(unfold(levels))
    diffs = []
    for j in range(1, 8):
        dd = u[j:] - u[:-j]
        diffs.append(dd[dd < xmax])
    dd = np.concatenate(diffs)
    h, edges = np.histogram(dd, bins=nb, range=(0, xmax))
    dens = h / (len(u) * (edges[1] - edges[0]))
    return 0.5 * (edges[1:] + edges[:-1]), dens


def ks2(a, b):
    a = np.sort(a); b = np.sort(b)
    allv = np.concatenate([a, b])
    ca = np.searchsorted(a, allv, side='right') / len(a)
    cb = np.searchsorted(b, allv, side='right') / len(b)
    return float(np.max(np.abs(ca - cb)))


def score(E, gam, idx):
    """Closeness score and signatures on the index set idx (0-based)."""
    En = E[idx]; g = gam[idx]
    good = np.isfinite(En)
    dl = (En[good] - g[good]) / mean_spacing(g[good])
    rms = float(np.sqrt(np.mean(dl ** 2))) if good.any() else np.nan
    # Gram interval agreement: same Gram interval index as the true zero
    gi_true = np.floor(theta(g[good]) / pi)
    gi_mod = np.floor(theta(En[good]) / pi)
    gram_hit = float(np.mean(gi_true == gi_mod))
    # nearest-zero hit: predicted level's nearest true zero is the right one
    near = float(np.mean(np.abs(dl) < 0.5))
    sp_m = spacings(En[good]); sp_z = spacings(g)
    ks = ks2(sp_m, sp_z)
    Ls = np.array([0.5, 1.0, 2.0])
    nv_m = number_variance(En[good], Ls); nv_z = number_variance(g, Ls)
    x, pc_m = pair_corr(En[good]); _, pc_z = pair_corr(g)
    return dict(rms=rms, bias=float(np.mean(dl)) if good.any() else np.nan,
                osc=float(np.std(dl)) if good.any() else np.nan,
                gram_hit=gram_hit, near=near, ks=ks,
                nv_err=float(np.max(np.abs(nv_m - nv_z))),
                pc_err=float(np.sqrt(np.mean((pc_m - pc_z) ** 2))),
                missing=float(np.mean(~good)))
