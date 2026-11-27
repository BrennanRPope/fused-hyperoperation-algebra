"""Lee-Yang rule, a starting point that works at every p (Nov 27).

s*G_a(s) = Im F(s), F(w) = E exp(i w X), X = U^(1/p), U ~ Beta(1, a)  (X is Kumaraswamy(p, a), density
a p x^(p-1) (1 - x^p)^(a-1)). Hermite-Biehler: if F has no zeros in Im w <= 0, then Re F and Im F have only real,
interlacing zeros, so G_a has only real zeros. This file checks that criterion numerically.
Run: python3 hb_check.py   (writes hb_check.json)
"""
import json, os, math, time
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
HERE = os.path.dirname(os.path.abspath(__file__)); t0 = time.time()
res = []
def ok(name, cond, detail):
    res.append(dict(name=name, passed=bool(cond), detail=detail)); print(('PASS ' if cond else 'FAIL ') + name + ' | ' + detail, flush=True)

NODES, WTS = leggauss(5000); XN = (NODES+1)/2; WN = WTS/2
def make_F(p, a):
    rho = WN*a*p*XN**(p-1)*(1-XN**p)**(a-1)
    def F(w):
        w = np.atleast_1d(np.asarray(w, dtype=complex)); out = np.empty(w.shape, complex)
        for i in range(0, len(w), 1500):
            out[i:i+1500] = np.exp(1j*np.outer(w[i:i+1500], XN)) @ rho
        return out
    dF = lambda w: np.exp(1j*w*XN) @ (rho*1j*XN)
    return F, dF
def height(p, a, w):  # two-end balance: B e^{iw}(iw)^-a = -A(-iw)^-p, A = a Gamma(p+1), B = p^a Gamma(a+1)
    return math.log(p**a*math.gamma(a+1)/(a*math.gamma(p+1))) + (p-a)*math.log(abs(w))
def lower_count(F, R, H, step=0.02):
    """zeros of F in the box -R < Re w < R, -H < Im w < 0, by the argument principle (path runs clockwise, hence the sign)"""
    seg = lambda z0, z1: z0 + (z1-z0)*np.linspace(0, 1, max(50, int(abs(z1-z0)/step)))
    path = np.concatenate([seg(-R, R), seg(R, R-1j*H), seg(R-1j*H, -R-1j*H), seg(-R-1j*H, -R)])
    v = F(path); ph = np.unwrap(np.angle(v)); return -int(round((ph[-1]-ph[0])/(2*np.pi))), float(np.min(np.abs(v[:len(path)//4])))
def zeros_near_line(p, a, F, dF, R):
    out = []
    for x in np.arange(6, R, 1.0):
        w = x + 1j*height(p, a, x)
        h0 = w.imag; conv = False
        for _ in range(60):
            d = dF(w)
            if not np.isfinite(d) or d == 0: break
            step = F(w)[0]/d; w = w - step
            if abs(w.imag - h0) > 40: break
            if abs(step) < 1e-4: conv = True; break
        if conv and 0 < w.real < R and abs(F(w)[0]) < 1e-3*np.min(np.abs(F(np.array([w+0.7, w-0.7, w+0.7j, w-0.7j])))) \
           and all(abs(w-z) > 1e-5 for z in out): out.append(w)
    return sorted(out, key=lambda z: z.real)

# B1 identity sG = Im F (series at 40 digits against the quadrature F)
mp.mp.dps = 40; gaps = []
for p, a, s in [(3, 2.5, 7.3), (9.5, 9.5, 31.0), (12.5, 12, 55.0)]:
    G = mp.quad(lambda x: (1-x**p)**a*mp.cos(s*x), [0, 0.5, 0.9, 1])
    F, _ = make_F(p, a); gaps.append(abs(float(s*G) - F(s)[0].imag))
ok('B1 s G_a(s) = Im F(s), F the characteristic function of X = U^(1/p), U ~ Beta(1, a)', max(gaps) < 1e-9, f'max gap {max(gaps):.1e} at (p,a,s) = (3,2.5,7.3), (9.5,9.5,31), (12.5,12,55)')

# B2 the Hermite-Biehler consequence: Re F and Im F zeros on the real line interlace (p = 9.5, a = 9.5)
F, dF = make_F(9.5, 9.5); s = np.linspace(0.01, 200, 400001); v = F(s)
zr = s[:-1][np.sign(v.real[:-1]) != np.sign(v.real[1:])]; zi = s[:-1][np.sign(v.imag[:-1]) != np.sign(v.imag[1:])]
allz = sorted([(z, 'r') for z in zr] + [(z, 'i') for z in zi]); alt = all(allz[i][1] != allz[i+1][1] for i in range(len(allz)-1))
ok('B2 consequence at p = 9.5, a = 9.5: zeros of Re F and Im F on (0, 200) alternate, as Hermite-Biehler says', alt and len(zi) > 50,
   f'{len(zr)} zeros of Re F, {len(zi)} of Im F (= real zeros of G), alternating: {alt}')

# B3 no zeros of F in the lower half-plane box, across the top strip, p from 1.2 to 20.5
rows = []
for p in [1.2, 1.5, 2.5, 3, 5, 7, 9.5, 10.5, 12.5, 15.3, 20.5]:
    for a in ([p/2, p] if p < 2 else [p-0.75, p]):
        F, dF = make_F(p, a); n, mn = lower_count(F, 250, 15); rows.append((p, a, n, mn))
ok('B3 F has no zeros in the box |Re w| < 250, -15 < Im w < 0, at both ends of the top strip, p = 1.2 to 20.5 (even and odd alike)',
   all(r[2] == 0 for r in rows), '; '.join(f'p={r[0]},a={r[1]:g}: {r[2]}' for r in rows))

# B4 the zeros of F run along the predicted height; at a = p it settles at (p - 1) ln p
hrows = []
for p in [2.5, 3, 5, 7, 9.5, 10.5, 12.5, 15.3]:
    F, dF = make_F(p, p); z = zeros_near_line(p, p, F, dF, 250)
    hrows.append((p, len(z), float(min(q.imag for q in z)), float(z[-1].imag), (p-1)*math.log(p)))
ok('B4 at a = p every zero of F lies above the real line, and the far ones settle at height (p - 1) ln p',
   all(r[2] > 0 and abs(r[3]-r[4]) < 0.03*r[4] + 0.05 for r in hrows),
   '; '.join(f'p={r[0]}: {r[1]} zeros, lowest {r[2]:.2f}, farthest {r[3]:.2f}, (p-1)ln p {r[4]:.2f}' for r in hrows))

# B5 above a = p the zeros sink like -(a - p) ln|w| and cross the real line; the criterion fails where the rule breaks
F, dF = make_F(3, 3.6); z = zeros_near_line(3, 3.6, F, dF, 250)
slope = np.polyfit(np.log([q.real for q in z[5:]]), [q.imag for q in z[5:]], 1)[0]
F2, _ = make_F(2.5, 5); n2, _ = lower_count(F2, 250, 15)
ok('B5 above a = p: at p = 3, a = 3.6 the heights fall with slope -(a - p) = -0.6 per ln|w|; at p = 2.5, a = 5 zeros already sit below the line',
   abs(slope+0.6) < 0.05 and n2 > 0, f'fitted slope {slope:.3f}; p=2.5, a=5: {n2} zeros in the lower box')

# B6 even p above a = p: G's zeros are real (proved, section N3 even) but F's cross below: the criterion is sufficient only
F, dF = make_F(4, 6); n4, _ = lower_count(F, 250, 15)
ok('B6 even p = 4, a = 6: F has zeros below the line (crossing near |w| = 143, as predicted) though G is real-rooted there',
   n4 > 0 and abs(math.exp(height(4, 6, 1)/2) - 143) < 2, f'{n4} zeros in the lower box; predicted crossing |w| = {math.exp(height(4, 6, 1)/2):.0f}')

print(f'\n{sum(r["passed"] for r in res)} of {len(res)} pass, {time.time()-t0:.0f} s')
json.dump(dict(results=res, strip=rows, heights=hrows), open(os.path.join(HERE, 'hb_check.json'), 'w'), indent=1)
