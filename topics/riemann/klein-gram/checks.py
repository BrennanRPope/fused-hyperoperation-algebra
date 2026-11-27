"""Checks for the Gram-test and prime-drift DRAFT (torus, Klein bottle, odd half). Numeric evidence, not proofs."""
import json, math, numpy as np, mpmath as mp
from zetas import Z, parts, theta_T, gram_points
from riemann_siegel import Z_rs
import gram as GR
from drift import coeffs, log_derivative_coeffs

mp.mp.dps = 25
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('PASS' if ok else 'FAIL'), name, info)

def L4m(s): return (mp.zeta(s, 0.25) - mp.zeta(s, 0.75))*mp.power(4, -s)
def ZK_mp(s, sign=-1): return 2*mp.zeta(s)*L4m(s) + sign*(1 - mp.power(2, 1 - 2*s))*mp.zeta(2*s)

# C1 sector formulas against direct lattice sums at s = 3
R = 300; m = np.arange(-R, R + 1)
def direct(sector, s=3.0):
    tot = 0.0
    for n in range(0, R + 1):
        lam = (m*m + n*n).astype(float)
        if n > 0: tot += np.sum(lam**(-s))
        else:
            mm = m[(m != 0) & ((m % 2 == 0) if sector == 'klein' else (m % 2 == 1))]
            tot += np.sum((mm*mm).astype(float)**(-s))
    return tot
for sec in ('klein', 'odd'):
    d = direct(sec); f = Z(np.array([3.0 + 0j]), sec)[0].real
    check(f'C1 {sec} sector formula = lattice sum at s=3', abs(d - f) < 1e-5, f'{d:.8f} vs {f:.8f}')

# C2 Klein + odd half = torus
s = np.array([0.3 + 17j, 1.7 + 240j, -1.2 + 55j])
check('C2 Z_K + Z_O = Z_T', np.max(np.abs(Z(s, 'klein') + Z(s, 'odd') - Z(s, 'torus'))/np.abs(Z(s, 'torus'))) < 1e-12)

# C3 fast evaluator against mpmath
err = max(abs(complex(Z(np.array([c]), 'klein')[0]) - complex(ZK_mp(mp.mpc(c.real, c.imag))))/abs(complex(ZK_mp(mp.mpc(c.real, c.imag))))
          for c in (0.5 + 700j, 1.3 + 1500j, -2.5 + 900j))
check('C3 vectorised Klein zeta matches mpmath (t up to 1500)', err < 1e-6, f'max rel err {err:.1e}')

# C4-C6 Gram counts to height 600 by the argument principle
ns, gs = gram_points(theta_T, 10.0, 600.0)
H = [GR.horizontal(g) for g in gs]; V = {x: GR.vertical(x, gs) for x in GR.XS}
xs = GR.XS[::-1]; cnt = {}; worst = 0.0
for j, name in enumerate(['right', 'middle', 'far_left']):
    a, b = xs[j+1], xs[j]
    for k in ('torus', 'klein', 'odd'):
        bottom = -np.array([H[i][j][k] for i in range(len(gs)-1)]); top = np.array([H[i+1][j][k] for i in range(len(gs)-1)])
        w = (bottom + V[b][k] + top - V[a][k])/(2*np.pi); worst = max(worst, np.max(np.abs(w - np.rint(w))))
        cnt[(name, k)] = np.rint(w).astype(int)
check('C4 every winding number is an integer to 0.05', worst < 0.05, f'worst {worst:.3f}')
tor_sc = GR.sign_change_counts(theta_T, lambda s: parts(s)[0], gs, h=0.002)
tot_strip = int(cnt[('middle', 'torus')].sum() + cnt[('right', 'torus')].sum() + cnt[('far_left', 'torus')].sum())
check('C5 torus: zeros in the strip = sign changes on the line (all on the line to height 600)',
      tot_strip == int(tor_sc.sum()) and cnt[('right', 'torus')].sum() == 0 and cnt[('far_left', 'torus')].sum() == 0,
      f'{tot_strip} vs {int(tor_sc.sum())}')
one = {k: np.mean((cnt[('middle', k)] + cnt[('right', k)]) == 1) for k in ('torus', 'klein', 'odd')}
check('C6 Gram share with one zero: Klein < torus < odd half (to height 600)', one['klein'] < one['torus'] < one['odd'],
      ', '.join(f'{k} {v:.3f}' for k, v in one.items()))
mean_cnt = {k: np.mean(cnt[('middle', k)] + cnt[('right', k)]) for k in ('klein', 'odd')}
check('C7 Klein and odd half keep one zero per Gram interval on average right of Re = -1',
      all(abs(v - 1) < 0.02 for v in mean_cnt.values()), str({k: round(v, 4) for k, v in mean_cnt.items()}))
fl = cnt[('far_left', 'klein')].sum()/(gs[-1] - gs[0])
check('C8 far-left band density = log 2 / pi per unit height (Klein)', abs(fl/(math.log(2)/math.pi) - 1) < 0.08,
      f'{fl:.4f} vs {math.log(2)/math.pi:.4f}')

# C9 mean Gram value = 2 a_T(1)/2 -/+ 1: torus 8, Klein 3, odd half 5
ns2, g2 = gram_points(theta_T, 10.0, 2000.0)
FT = parts(0.5 + 1j*g2); ph = np.exp(1j*theta_T(g2)); sg = (-1.0)**ns2
mv = {'torus': np.mean(sg*np.real(ph*FT[0])), 'klein': np.mean(sg*np.real(ph*(FT[0]/2 - FT[1]))),
      'odd': np.mean(sg*np.real(ph*(FT[0]/2 + FT[1])))}
check('C9 mean Gram values 8, 3, 5 (torus, Klein, odd half) over 3473 Gram points',
      abs(mv['torus'] - 8) < 0.1 and abs(mv['klein'] - 3) < 0.1 and abs(mv['odd'] - 5) < 0.1,
      ', '.join(f'{k} {v:.3f}' for k, v in mv.items()))

# C10 lattice identity N_K = N_T/2 -/+ [floor sqrt X odd]
X = 200000; A = coeffs(X)
NT, NK, NO = (np.cumsum(A[k]) for k in ('torus', 'klein', 'odd'))
r = np.array([math.isqrt(n) % 2 for n in range(X + 1)], float)
check('C10 lattice count: N_K = N_T/2 - [floor sqrt X odd], N_O = N_T/2 + [...] for X <= 2e5',
      np.max(np.abs(NK - (NT/2 - r))) == 0 and np.max(np.abs(NO - (NT/2 + r))) == 0)

# C11 torus log-derivative coefficients = von Mangoldt of Z[i]
cT = log_derivative_coeffs(A['torus'][:2001])
def lam_gauss(n):
    from sympy import factorint
    f = factorint(n)
    if len(f) != 1: return 0.0
    (p, e), = f.items()
    if p == 2: return math.log(2)
    if p % 4 == 1: return 2*math.log(p)
    return 2*math.log(p) if e % 2 == 0 else 0.0
err = max(abs(cT[n] - lam_gauss(n)) for n in range(2, 2001))
check('C11 torus prime count = Gaussian-prime von Mangoldt (n <= 2000)', err < 1e-8, f'max err {err:.1e}')

# C12 Klein has no Euler product: composite "primes"
cK = log_derivative_coeffs(A['klein'][:101])
check('C12 Klein log-derivative is nonzero at 10 = 2*5 (no Euler product)', abs(cK[10]) > 1, f'c(10) = {cK[10]:.4f}')

# C13 rightmost Klein zero (Re > 1) confirmed by mpmath; odd half zero at s = 0
Zs = json.load(open('zeros.json'))
zk = np.array([complex(a, b) for a, b in Zs['klein']['low']]); top = zk[np.argmax(zk.real)]
r0 = mp.findroot(ZK_mp, mp.mpc(top.real, top.imag))
check('C13 a Klein zero with Re > 1 (mpmath)', mp.re(r0) > 1.2 and abs(ZK_mp(r0)) < 1e-15, f'rho = {complex(r0):.6f}')
check('C14 odd half vanishes at s = 0', abs(ZK_mp(mp.mpf(0), +1)) < 1e-20)

# C15 explicit formula: Klein's integrated prime drift = minus the sum over its zeros
Xc = 1_000_000; aK = coeffs(Xc)['klein']; psi = np.cumsum(log_derivative_coeffs(aK))
p1 = np.cumsum(psi - np.arange(Xc + 1))
grid = np.array([2e5, 4e5, 6e5, 8e5, 1e6]); zz = zk[zk.real > -1]
pred = -2*np.real(np.sum(grid[:, None]**(zz[None, :] + 1)/(zz*(zz + 1))[None, :], axis=1))
ratio = p1[grid.astype(int)]/pred
check('C15 explicit formula rebuilds Klein drift from its zeros to 2% (X = 2e5..1e6)', np.all(np.abs(ratio - 1) < 0.02), np.round(ratio, 4).tolist())
check('C16 Klein prime drift outgrows X (|psi - X| > X at X = 1e6)', abs(psi[Xc] - Xc) > Xc, f'{psi[Xc] - Xc:.3e}')

# C17 Riemann-Siegel evaluator
e = max(abs(Z_rs(t)[0] - float(mp.siegelz(t))) for t in (5000.5, 3e5 + 0.25))
check('C17 Riemann-Siegel Z(t) matches mpmath', e < 1e-4, f'max err {e:.1e}')

print(f'{sum(res)}/{len(res)} pass')
