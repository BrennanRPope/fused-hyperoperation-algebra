"""Round 7, item 3: random-matrix statistics of zeros on the critical line.

Compare (unfolded by the exact smooth phase theta(t)/pi):
  zeta                   one Euler product (degree 1)         heights 10000-11600
  L(s, chi_-4)           one Euler product (degree 1)         heights 1000-3000
  x^2+y^2  = 4 zeta L(chi_-4)   product of two Euler products   heights 500-1400
  x^2+5y^2 = zeta L(chi_-20) + L(chi_-4) L(chi_5)   no Euler product  heights 500-1400
  Davenport-Heilbronn    no Euler product, degree 1           heights 500-2300
against GUE (Gaudin, exact Fredholm determinant), two independent GUE superposed, Poisson.
Also the count deficit [theta(T2)-theta(T1)]/pi - N_on, which counts zeros off the line.
"""
import json, math
import numpy as np
from scipy.special import loggamma
from fastl import hurwitz, dirichlet_series_periodic, char_coeffs, zeta

SQ5 = math.sqrt(5)
XI = (math.sqrt(10 - 2*SQ5) - 2)/(SQ5 - 1)
DH = [1, XI, -XI, -1, 0]
C4, C5, C20 = char_coeffs(-4), char_coeffs(5), char_coeffs(-20)

def fam(name):
    if name == 'zeta':
        th = lambda t: (loggamma((0.5+1j*t)/2)).imag - t/2*math.log(math.pi)
        f = lambda s: zeta(s)
    elif name == 'L-4':
        th = lambda t: (loggamma((1.5+1j*t)/2)).imag + t/2*math.log(4/math.pi)
        f = lambda s: dirichlet_series_periodic(s, C4)
    elif name == 'x2+y2':
        th = lambda t: (loggamma(0.5+1j*t)).imag + t*math.log(2/(2*math.pi))
        f = lambda s: 4*zeta(s)*dirichlet_series_periodic(s, C4)
    elif name == 'x2+5y2':
        th = lambda t: (loggamma(0.5+1j*t)).imag + t*math.log(math.sqrt(20)/(2*math.pi))
        f = lambda s: zeta(s)*dirichlet_series_periodic(s, C20) + dirichlet_series_periodic(s, C4)*dirichlet_series_periodic(s, C5)
    elif name == 'DH':
        th = lambda t: (loggamma((1.5+1j*t)/2)).imag + t/2*math.log(5/math.pi)
        f = lambda s: dirichlet_series_periodic(s, DH)
    return th, f

def hardy(name, t):
    th, f = fam(name)
    v = np.exp(1j*th(t))*f(0.5+1j*t)
    return v

def zeros_on_line(name, T1, T2, step=0.02):
    t = np.arange(T1, T2, step)
    v = hardy(name, t)
    imrel = float(np.max(np.abs(v.imag))/np.max(np.abs(v.real)))
    r = v.real
    idx = np.where(np.sign(r[:-1]) * np.sign(r[1:]) < 0)[0]
    a, b = t[idx].copy(), t[idx+1].copy()
    fa = r[idx].copy()
    for _ in range(30):
        m = (a+b)/2
        fm = hardy(name, m).real
        left = np.sign(fm) == np.sign(fa)
        a = np.where(left, m, a); fa = np.where(left, fm, fa)
        b = np.where(left, b, m)
    return (a+b)/2, imrel

# ---- GUE spacing via Bornemann's Fredholm determinant ----
def gue_E(s, m=40):
    x, w = np.polynomial.legendre.leggauss(m)
    out = []
    for L in np.atleast_1d(s):
        if L == 0: out.append(1.0); continue
        y = (x+1)*L/2; ww = w*L/2
        d = y[:, None]-y[None, :]
        with np.errstate(invalid='ignore', divide='ignore'):
            K = np.where(d == 0, 1.0, np.sin(np.pi*d)/(np.pi*d))
        M = np.eye(m) - np.sqrt(ww)[:, None]*K*np.sqrt(ww)[None, :]
        out.append(np.linalg.det(M))
    return np.array(out)

def spacing_densities(sgrid):
    h = sgrid[1]-sgrid[0]
    E = gue_E(sgrid)
    p_gue = np.gradient(np.gradient(E, h), h)
    E2 = gue_E(sgrid/2)**2
    p_sup = np.gradient(np.gradient(E2, h), h)
    return p_gue, p_sup, np.exp(-sgrid)

def ks(sample, sgrid, pdf):
    cdf = np.cumsum(pdf)*(sgrid[1]-sgrid[0]); cdf /= cdf[-1]
    xs = np.sort(sample)
    emp = np.arange(1, len(xs)+1)/len(xs)
    model = np.interp(xs, sgrid, cdf)
    return float(np.max(np.abs(emp-model)))

def pair_corr(x, maxd=3.0, nb=60):
    d = []
    for i in range(len(x)):
        j = i+1
        while j < len(x) and x[j]-x[i] < maxd:
            d.append(x[j]-x[i]); j += 1
    h, e = np.histogram(d, bins=nb, range=(0, maxd))
    return (h/(len(x)*(e[1]-e[0]))).tolist(), ((e[1:]+e[:-1])/2).tolist()

if __name__ == '__main__':
    runs = {'zeta': (10000, 11600), 'L-4': (1000, 3000), 'x2+y2': (500, 1400),
            'x2+5y2': (500, 1400), 'DH': (500, 2300)}
    sgrid = np.linspace(0, 4, 401)
    p_gue, p_sup, p_poi = spacing_densities(sgrid)
    out = {'sgrid': sgrid.tolist(), 'p_gue': p_gue.tolist(), 'p_sup': p_sup.tolist(), 'p_poi': p_poi.tolist()}
    for name, (T1, T2) in runs.items():
        if name == 'x2+y2':
            # product of two Euler products: merge the zeros of each factor (a scan misses near-coincident pairs)
            g1, i1 = zeros_on_line('zeta', T1, T2)
            g2, i2 = zeros_on_line('L-4', T1, T2)
            g, imrel = np.sort(np.concatenate([g1, g2])), max(i1, i2)
        else:
            g, imrel = zeros_on_line(name, T1, T2)
        th = fam(name)[0]
        x = th(g)/math.pi
        sp = np.diff(x)
        expected = (th(np.array([T2]))[0]-th(np.array([T1]))[0])/math.pi
        deficit = expected - len(g)
        k_gue, k_sup, k_poi = ks(sp, sgrid, p_gue), ks(sp, sgrid, p_sup), ks(sp, sgrid, p_poi)
        small = float(np.mean(sp < 0.25))
        spr = sp/sp.mean()
        ksr = [ks(spr, sgrid, p_gue), ks(spr, sgrid, p_sup), ks(spr, sgrid, p_poi)]
        smallr = float(np.mean(spr < 0.25))
        pc, pcx = pair_corr(x)
        print(f"{name:7s} T=[{T1},{T2}] zeros on line {len(g):5d}, phase reality {imrel:.1e}, "
              f"expected {expected:8.1f}, deficit {deficit:+.1f} | mean spacing {sp.mean():.4f} | "
              f"KS: GUE {k_gue:.4f} superposed {k_sup:.4f} Poisson {k_poi:.4f} | P(s<0.25) {small:.4f}"
              f" | rescaled to mean 1: KS GUE {ksr[0]:.4f} sup {ksr[1]:.4f} Poi {ksr[2]:.4f}, P(s<0.25) {smallr:.4f}")
        out[name] = {'T': [T1, T2], 'n': len(g), 'deficit': deficit, 'spacings': sp.tolist(),
                     'ks': [k_gue, k_sup, k_poi], 'small': small, 'ks_rescaled': ksr, 'small_rescaled': smallr, 'pc': pc, 'pcx': pcx,
                     'first': g[:5].tolist()}
    # reference small-spacing probabilities
    h = sgrid[1]-sgrid[0]; i = int(0.25/h)
    out['ref_small'] = {'gue': float(np.sum(p_gue[:i])*h), 'sup': float(np.sum(p_sup[:i])*h), 'poi': float(1-math.exp(-0.25))}
    print("reference P(s<0.25): GUE %.4f superposed %.4f Poisson %.4f" % tuple(out['ref_small'].values()))
    print("KS 5%% critical value ~ 1.36/sqrt(n)")
    json.dump(out, open('rmt.json', 'w'))
