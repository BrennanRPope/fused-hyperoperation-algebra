"""Every zero of the Klein and odd-half zetas with -3 < Re < 2.5, 1 < Im < 600, plus the near-line
zeros (-0.9 < Re < 2.4) in the window 1900 < Im < 2000, plus real zeros in (-3, 1).
Spurious grid cells are dropped unless Newton lands on |F| < 1e-8 |F(z + 0.05)|."""
import json, numpy as np
from scipy.optimize import brentq
from zetas import Z
from zerofind import find_zeros

def clean(z, kind):
    if len(z) == 0: return z
    r = np.abs(Z(z, kind)); ref = np.abs(Z(z + 0.05, kind))
    return z[r < 1e-8*ref]

out = {}
for kind in ('torus', 'klein', 'odd'):
    lo = clean(find_zeros(kind, -3, 2.5, 1, 600, ds=0.02, dt=0.02, band=25), kind)
    hi = clean(find_zeros(kind, -0.9, 2.4, 1900, 2000, ds=0.02, dt=0.01, band=25), kind)
    f = lambda x: Z(np.array([x + 0j]), kind)[0].real
    xs = np.linspace(-2.999, 0.999, 4001); v = np.array([f(x) for x in xs])
    real = [brentq(f, xs[i], xs[i+1]) for i in range(len(xs)-1) if np.sign(v[i]) != np.sign(v[i+1]) and abs(xs[i] - round(xs[i])) > 1e-9 or (np.sign(v[i]) != np.sign(v[i+1]))]
    real = sorted(set(round(r, 10) for r in real))
    out[kind] = {'low': [[c.real, c.imag] for c in lo], 'high': [[c.real, c.imag] for c in hi], 'real': real}
    print(kind, 'zeros 1<Im<600:', len(lo), ' 1900<Im<2000 near line:', len(hi), ' real zeros:', np.round(real, 6).tolist(), flush=True)
json.dump(out, open('zeros.json', 'w'))
