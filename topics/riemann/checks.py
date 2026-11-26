"""Checks for the Riemann DRAFT. Run: python3 checks.py > checks_output.txt
Reads the JSON files written by cube.py, weil.py, ladder.py, up_step_*.py and dbn_scan.py, recomputes the cheap
parts, and adds a polishing pass (marked POLISH) that re-derives a few numbers a second way."""
import json, time, math, numpy as np, mpmath as mp
from up_step_grid import F, Xi
T0 = time.time(); results = []
def check(name, ok, detail):
    results.append(ok); print(('PASS' if ok else 'FAIL'), name, '|', detail, flush=True)

# ---------- A. prime cube ----------
cube = json.load(open('cube.json'))
gap = max(abs(r['fold_lift'] - r['product']) for r in cube['fold_rows'])
check('A1 fold of ∘ₖ over p^-s lifts to the Euler product (10 rows)', gap < 1e-12, 'max gap %.1e' % gap)
gap2 = max(abs(r['product'] - r['zeta_k_G'])/abs(r['product']) for r in cube['fold_rows'])
check('A2 product = zeta^k * G_k (primes < 2e5)', gap2 < 2e-6, 'max rel gap %.1e (prime truncation)' % gap2)
mp.mp.dps = 30
r1 = [r for r in cube['fold_rows'] if r['s'] == 2 and r['k'] == 1][0]['product']
rm = [r for r in cube['fold_rows'] if r['s'] == 2 and r['k'] == -1][0]['product']
check('A3 k=1 is zeta(s)/zeta(2s), k=-1 is 1/zeta(s) at s=2', abs(r1 - float(mp.zeta(2)/mp.zeta(4))) < 1e-5 and abs(rm - float(1/mp.zeta(2))) < 1e-5,
      'k=1 %.7f vs %.7f; k=-1 %.7f vs %.7f' % (r1, mp.zeta(2)/mp.zeta(4), rm, 1/mp.zeta(2)))
# local zeros 1 + k p^-s = 0 at Re s = log|k|/log p
worst = 0; inside = []
for k in (0.5, 1, 2, 3):
    for p in (2, 3, 5, 7, 11):
        s = mp.log(k)/mp.log(p) + 1j*mp.pi/mp.log(p)   # p^s = -k
        worst = max(worst, abs(1 + k*mp.mpf(p)**(-s)))
        sig = math.log(k)/math.log(p)
        if 0 < sig < 1: inside.append((k, p, round(sig, 3)))
check('A4 local zeros of the cube sit on Re s = log|k|/log p; none in the strip for |k| <= 1', worst < 1e-25 and all(k > 1 for k, p, s in inside),
      'residual %.1e; inside the strip: %s' % (worst, inside))

# ---------- B. corner-sign walk ----------
N = 10**7
mu = np.ones(N+1, dtype=np.int8); mu[0] = 0
isp = np.ones(N+1, dtype=bool); isp[:2] = False
for p in range(2, int(N**0.5)+1):
    if isp[p]: isp[p*p::p] = False
for p in np.nonzero(isp)[0]:
    mu[p::p] *= -1
    if p*p <= N: mu[p*p::p*p] = 0
M = np.cumsum(mu.astype(np.int64))
mx = (np.abs(M[100:]) / np.sqrt(np.arange(100, N+1))).max()
check('B1 Mertens walk M(10^7) = 1037 and |M(x)| < 0.6 sqrt(x) for 100 <= x <= 10^7', M[N] == 1037 and mx < 0.6, 'M(1e7)=%d, max |M|/sqrt x = %.4f' % (M[N], mx))

# ---------- C. local factors as dial points ----------
weil = json.load(open('weil.json'))
g1 = weil['genus1']; g2 = weil['genus2']
check('C1 genus 1 (y^2=x^3+x+1): every dial point (-a_p, p) is elliptic, Delta <= 0', all(r['Delta'] <= 0 for r in g1),
      '%d primes < 3000, max Delta/4p = %.4f' % (len(g1), max((r['Delta'])/(4*r['p']) for r in g1)))
check('C2 genus 2 (y^2=x^5-x+1): L-polynomial splits into two real quadratics, both elliptic', all(r['AB_real'] and r['DeltaA'] <= 1e-9 and r['DeltaB'] <= 1e-9 for r in g2),
      '%d primes < 400' % len(g2))
check('C3 point counts obey |N1 - p - 1| <= 4 sqrt p (genus 2)', all(abs(r['N1'] - r['p'] - 1) <= 4*math.sqrt(r['p']) for r in g2), 'consistency of the counts')
th = np.array([r['theta'] for r in g1]); hist, edges = np.histogram(th, bins=8, range=(0, np.pi))
st = np.array([(edges[i+1]-edges[i] - (np.sin(2*edges[i+1]) - np.sin(2*edges[i]))/2)/np.pi for i in range(8)])*len(th)
chi2 = float(((hist - st)**2/st).sum())
check('C4 Frobenius angles follow the Sato–Tate curve (2/pi) sin^2', chi2 < 20, 'chi^2 = %.1f on 7 dof' % chi2)

# ---------- D. one step up the ladder ----------
up = json.load(open('up_step.json'))
def PhiT(u):  # Titchmarsh's Phi (factor 2), Xi(t) = 2 int_0^inf Phi cos
    u = mp.mpf(u); return 2*mp.nsum(lambda n: (2*mp.pi**2*n**4*mp.e**(9*u/2) - 3*mp.pi*n**2*mp.e**(5*u/2))*mp.e**(-mp.pi*n**2*mp.e**(2*u)), [1, mp.inf])
M1 = mp.quad(lambda t: mp.re(Xi(t)), mp.linspace(0, 120, 61))
check('D1 the push: M_1 = int_0^inf Xi = pi Phi(0) > 0', abs(M1 - mp.pi*PhiT(0)) < 1e-15, 'M_1 = %s, pi Phi(0) = %s' % (mp.nstr(M1, 15), mp.nstr(mp.pi*PhiT(0), 15)))
I = mp.mpf(0); mn = mp.inf
for a in range(0, 100):
    I += mp.quad(lambda t: mp.re(Xi(t)), [a, a+1]); mn = min(mn, I/(a+1))
tail = mp.quad(lambda t: abs(Xi(t)), [100, 140])
check('D2 F_1(x) = (1/x) int_0^x Xi > 0 on [1, 100]; beyond, the wave is far below the push', mn > 0 and tail < 1e-10, 'min F_1 = %.4f, int_100^140 |Xi| = %.1e' % (mn, tail))
z1 = up['1']['zeros']; z2 = up['2']['zeros']
check('D3 every zero found of F_1 and F_2 (Re < 60) is off the line, |Im| > 8', len(z1) >= 9 and len(z2) >= 9 and min(abs(z[1]) for z in z1 + z2) > 8,
      'F_1: %d zeros, Im from %.2f to %.2f; F_2: %d zeros' % (len(z1), min(z[1] for z in z1), max(z[1] for z in z1), len(z2)))
# POLISH: argument principle on the box [0.5, 60] x [-8, 8] (F_1 is real on the real line, so the box sees every near-line zero)
mp.mp.dps = 20
def arg_count(f, x0, x1, y0, y1, n=240):
    pts = [mp.mpc(x0 + (x1-x0)*i/n, y0) for i in range(n)] + [mp.mpc(x1, y0 + (y1-y0)*i/n) for i in range(n)] + \
          [mp.mpc(x1 - (x1-x0)*i/n, y1) for i in range(n)] + [mp.mpc(x0, y1 - (y1-y0)*i/n) for i in range(n)]
    vals = [f(z) for z in pts]; vals.append(vals[0]); tot = 0
    for a, b in zip(vals[:-1], vals[1:]): tot += float(mp.arg(b/a))
    return tot/(2*math.pi)
t1 = time.time(); cnt = arg_count(lambda z: F(z, 1), 0.5, 60, -8, 8, n=150)
check('D4 POLISH argument principle: F_1 has no zeros in [0.5, 60] x [-8, 8]', abs(cnt) < 0.1, 'winding %.3f (%.0f s)' % (cnt, time.time()-t1))
moves = []
for z in z1[:3]:
    z0 = mp.mpc(z[0] + 0.7, z[1] - 0.5); zz = mp.findroot(lambda w: F(w, 1), z0); moves.append(abs(complex(zz) - complex(*z)))
check('D5 POLISH zeros of F_1 re-found from shifted starts', max(moves) < 1e-8, 'max move %.1e' % max(moves))

# ---------- E. heat flow (de Bruijn–Newman), Rodgers–Tao convention ----------
dbn = json.load(open('dbn.json'))
g = [float(mp.zetazero(n).imag) for n in range(1, 19)]
ch0 = dbn['0.00']['changes']
check('E1 at t = 0 the 18 sign changes of H_0 in [0,150] are the first 18 zeros, x = 2 gamma', len(ch0) == 18 and max(abs(c + 0.125 - 2*gm) for c, gm in zip(ch0, g)) < 0.25,
      'max offset %.3f (grid 0.25)' % max(abs(c + 0.125 - 2*gm) for c, gm in zip(ch0, g)))
ts = sorted(dbn, key=float); counts = [len(dbn[t]['changes']) for t in ts]
drops = [(ts[i], counts[i-1], counts[i]) for i in range(1, len(ts)) if counts[i] != counts[i-1]]
check('E2 for t >= -1.2 all 18 stay real; going down, they leave in pairs', all(len(dbn[t]['changes']) == 18 for t in ts if float(t) >= -1.2) and all(abs(a - b) == 2 for _, a, b in drops),
      'count changes at t = %s' % [d[0] for d in drops])
# POLISH: H_t with 384 nodes instead of 768 at t = 0.2, x = 2 gamma_1 region; refine a zero with findroot from a complex start
mp.mp.dps = 45
def PhiRT(u):
    s = mp.mpf(0)
    for n in range(1, 8): s += (2*mp.pi**2*n**4*mp.e**(9*u) - 3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u))
    return s
U = mp.mpf('1.3'); pts = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(8, mp.mp.prec)
us = [(x+1)*U/2 for x, w in pts]; ws = [w*U/2 for x, w in pts]; Ph = [PhiRT(u) for u in us]
H = lambda t, z: mp.fsum(w*p*mp.e**(t*u*u)*mp.cos(z*u) for u, w, p in zip(us, ws, Ph))
zr = mp.findroot(lambda z: H(mp.mpf('0.2'), z), mp.mpc(28, 0.5))
check('E3 POLISH H_0.2 zero refined from a start off the line lands back on it', abs(mp.im(zr)) < 1e-20, 'zero at %s' % mp.nstr(zr, 12))
zr0 = mp.findroot(lambda z: H(0, z), mp.mpf(28.2))
check('E4 POLISH half the nodes still give x_1 = 2 gamma_1', abs(zr0 - 2*mp.zetazero(1).imag) < 1e-25, '%s vs %s' % (mp.nstr(zr0, 20), mp.nstr(2*mp.zetazero(1).imag, 20)))

# ---------- F. kernel ladder ----------
lad = json.load(open('ladder.json'))
even_ok = all(len(lad[p]['sign_changes']) >= 15 for p in ('4', '6'))
odd = {p: (len(lad[p]['sign_changes']), lad[p]['ratio_at_40']) for p in ('2.5', '3', '3.5', '5')}
check('F1 even p keeps its zeros coming; non-even p stops, with the tail equal to the kink push', even_ok and all(abs(r - 1) < 0.01 for n, r in odd.values()),
      'even: %s zeros in [0,40]; non-even (zeros, tail/push at 40): %s' % ([len(lad[p]['sign_changes']) for p in ('4', '6')], odd))
mp.mp.dps = 40
def Gp(p, x):
    p = mp.mpf(p); Uu = mp.mpf(110)**(1/p); n = max(4, int(abs(x)*Uu/mp.pi) + 2)
    return mp.quad(lambda u: mp.exp(-u**p)*mp.cos(x*u), mp.linspace(0, Uu, n+1))
a = mp.findroot(lambda x: Gp(3, x), 7.75); b = mp.findroot(lambda x: Gp(3, x), 8.05)
check('F2 POLISH p = 3: the last two real zeros are a close pair (about to leave the line)', abs(a - b) > 0.05 and abs(a - b) < 0.5, 'zeros %s and %s' % (mp.nstr(a, 10), mp.nstr(b, 10)))

print('\n%d of %d pass, %.0f s' % (sum(results), len(results), time.time() - T0))
