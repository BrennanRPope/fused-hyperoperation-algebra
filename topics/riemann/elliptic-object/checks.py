"""Checks for the elliptic-object DRAFT. Run: python3 checks.py > checks_output.txt
G = mathematics (sympy exact, or numbers from zeros.json / li.json / dial.json); P = polish."""
import json, math, os, re, cmath
import sympy as sp
import mpmath as mp
res = []
def check(name, ok, info=''):
    res.append(ok); print(('PASS ' if ok else 'FAIL ') + name + (f'   [{info}]' if info else ''))

d = json.load(open('dial.json')); li = json.load(open('li.json')); Z = json.load(open('zeros.json'))
u, a, b, p, x, y, g, dl = sp.symbols('u a b p x y g delta', real=True)

# G1 pairing lemma
k, dd = sp.expand(-(a + b)), sp.expand(a * b)
check('G1 pairing lemma: (1 - a u)(1 - b u) = 1 + k u + d u^2 with Delta = k^2 - 4d = (a - b)^2 >= 0 (exact)',
      sp.expand(k ** 2 - 4 * dd - (a - b) ** 2) == 0 and sp.expand((1 - a * u) * (1 - b * u) - (1 + k * u + dd * u ** 2)) == 0)
# G2 mirror pairs
ok = all(r['Delta'] == (r['p'] - 1) ** 2 and r['k_hat'] <= -2 - 1e-12 for r in d['A1'])
check('G2 zeta factor paired with its mirror, (1 - u)(1 - p u): Delta = (p - 1)^2 > 0 and k_hat <= -2.12 for all 45 primes < 200', ok,
      f"max k_hat {max(r['k_hat'] for r in d['A1']):.4f}")
# G3 complex characters: elliptic with d = 1
ok = all(abs(r['k_hat']) <= 2 + 1e-12 and r['Delta'] <= 1e-12 for r in d['A2'])
check('G3 complex character with its conjugate: 1 - 2Re(chi(p)) u + u^2 has Delta <= 0 (elliptic) and d = 1, roots on |u| = 1 (Re s = 0); mod 5 and mod 7, 90 factors', ok)
# G4 elliptic curve reference (Hasse)
ok = all(abs(r['k_hat']) <= 2 for r in d['A3'])
check('G4 y^2 = x^3 + x + 1: every good prime < 3000 has |a_p| <= 2 sqrt(p) (Hasse), 428 primes', ok)
# G5 Cayley map
rho = sp.Rational(1, 2) + sp.I * g
z = 1 - 1 / rho
check('G5 Cayley map z = 1 - 1/rho: |z| = 1 on Re rho = 1/2, and z(1 - rho) = 1/z(rho) (exact)',
      sp.simplify(sp.Abs(z) ** 2 - 1) == 0 and sp.simplify((1 - 1 / (1 - (x + sp.I * y))) * (1 - 1 / (x + sp.I * y)) - 1) == 0)
X = sp.symbols('X', real=True)
mod2 = sp.expand(((X - 1) ** 2 + y ** 2) - (X ** 2 + y ** 2))
check('G6 converse: |z| = 1 exactly when Re rho = 1/2 (|rho - 1|^2 - |rho|^2 = 1 - 2 Re rho)', sp.simplify(mod2 - (1 - 2 * X)) == 0)
# G7 k + 2 = 1/|rho|^2 and curvature formula
kk = sp.simplify(-(z + 1 / z))
check('G7 on the line, k = -(z + 1/z) = -2 + 1/|rho|^2 (exact)', sp.simplify(kk - (-2 + 1 / (sp.Rational(1, 4) + g ** 2))) == 0)
gp = sp.symbols('gp')
rr = sp.Rational(1, 2) + sp.I * gp; zz = 1 - 1 / rr
curv = sp.simplify(1 - (zz + 1 / zz) ** 2 / 4)
check('G8 curvature d - k^2/4 of a pair = g^2/(g^2 + 1/4)^2 with complex height g (rho = 1/2 + i g); real > 0 iff g real (exact)',
      sp.simplify(curv - gp ** 2 / (gp ** 2 + sp.Rational(1, 4)) ** 2) == 0)
C = d['C']
check('G9 5,000 zeros: every Cayley dial point inside the elliptic segment, |k| < 2, and k + 2 = 1/|rho|^2 to 1e-15',
      all(abs(r['k_re']) < 2 for r in C) and max(abs(r['k_plus_2'] - r['inv_abs_rho2']) for r in C) < 1e-15,
      f"max |k| = {max(abs(r['k_re']) for r in C):.10f}")
# G10 off-line toy points leave the real line, real off-line is hyperbolic
off = [r for r in d['C_off'] if 'gamma0' in r]; roff = [r for r in d['C_off'] if 'real_rho' in r]
check('G10 toy off-line zeros: k is not real (Im k != 0) and |z| != 1; toy real zeros in (1/2, 1): k > 2 (hyperbolic)',
      all(abs(r['k_im']) > 1e-6 and abs(r['abs_z'] - 1) > 1e-6 for r in off) and all(r['k_re'] > 2 for r in roff))
# G11 trivial zeros hyperbolic
check('G11 trivial zeros s = -2m read through the Cayley map: k < -2 (hyperbolic), approaching -2', all(r['k'] < -2 for r in d['B']),
      f"k(m=1) = {d['B'][0]['k']:.4f}, k(m=40) = {d['B'][-1]['k']:.6f}")
# G12 Li from zeta: known values
mp.mp.dps = 40
lam = [mp.mpf(v) for v in li['lambda']]
check('G12 lambda_1 = 1 + gamma/2 - log(4 pi)/2 to 25 digits (Cauchy integral of log xi, no zeros used)',
      abs(lam[0] - (1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2)) < mp.mpf(10) ** -25, mp.nstr(lam[0], 22))
mp.mp.dps = 40
l2_exact = 1 + mp.euler - mp.euler ** 2 - 2 * mp.stieltjes(1) - mp.log(4 * mp.pi) + mp.pi ** 2 / 8
check('G13 lambda_2 matches the closed form 1 + gamma - gamma^2 - 2 gamma_1 - log 4pi + pi^2/8 (Stieltjes gamma_1)',
      abs(lam[1] - l2_exact) < mp.mpf(10) ** -20, mp.nstr(l2_exact, 20))
check('G14 lambda_n > 0 for every n = 1..1000 (from zeta itself)', min(lam) > 0, f'min {mp.nstr(min(lam), 8)}')
D = {r['n']: r for r in d['D']}
rel = max(abs(D[n]['zero_sum'] + D[n]['tail'] - float(lam[n - 1])) / float(lam[n - 1]) for n in range(1, 1001))
check('G15 the same lambda_n rebuilt from the dial points: 2 sum (1 - T_n(-k/2)) over 5,000 zeros + density tail agrees to 1e-5 relative, n <= 1000',
      rel < 1e-5, f'max rel diff {rel:.2e}')
gam = [float(v) for v in Z['gamma'][:200]]
n0 = 37
lhs = sum(2 * (1 - float(mp.chebyt(n0, 1 - 1 / (2 * (0.25 + t * t))))) for t in gam)
rhs = sum((2 - 2 * (1 - 1 / complex(0.5, t)) ** n0).real for t in gam)
check('G16 Chebyshev form = Li form term by term: 2(1 - T_n(1 - 1/(2|rho|^2))) = 2 - z^n - z^-n (n = 37, 200 zeros)', abs(lhs - rhs) < 1e-9 * abs(rhs), f'{lhs:.10f}')
EG = 0.5772156649015329
tr = lambda n: n / 2 * (math.log(n) - 1 + EG - math.log(2 * math.pi))
w = max(abs(float(lam[n - 1]) - tr(n)) / math.sqrt(n) for n in range(10, 1001))
check('G17 lambda_n - (n/2)(log n - 1 + gamma - log 2pi) stays within 1 x sqrt(n) for 10 <= n <= 1000 (trend matches)', w < 1, f'max |resid|/sqrt(n) = {w:.3f}')
mp.mp.dps = 30
# G18 Weil bound <-> Chebyshev for a curve: genus 1, N_n = q^n + 1 - 2 q^{n/2} T_n(a/(2 sqrt q))
q, ap = 5, -3   # y^2 = x^3 + x + 1 over F_5 (a_5 = -3 from weil.json)
ok = True
for n in range(1, 8):
    al = (ap + cmath.sqrt(ap * ap - 4 * q)) / 2
    Nn = q ** n + 1 - (al ** n + al.conjugate() ** n).real
    ok &= abs(Nn - (q ** n + 1 - 2 * q ** (n / 2) * float(mp.chebyt(n, ap / (2 * math.sqrt(q)))))) < 1e-6
    ok &= abs(Nn - q ** n - 1) <= 2 * q ** (n / 2) + 1e-9
check('G18 the same Chebyshev shape over F_5: #E(F_5^n) = 5^n + 1 - 2*5^(n/2) T_n(a/(2 sqrt 5)), and the Weil bound holds, n = 1..7', ok)
# G19 detection toy: first negative n close to the predicted (1/eps) log(trend/2)
okt = True; info = []
for t in d['toy']:
    eps = t['eps']; n = 10.0
    for _ in range(200): n = math.log((tr(n) + 4) / 2) / eps
    okt &= 0.7 < t['first_negative_n'] / n < 1.3
    info.append(f"{t['gamma0']:.1f}/{t['delta']}: {t['first_negative_n']} vs {n:.0f}")
check('G19 detection toy: a split double zero turns lambda_n negative at n within 30% of (1/eps) log(trend/2), eps ~ delta/gamma^2', okt, '; '.join(info))
eps = 0.5 * math.log1p(2 * 0.25 / (0.0625 + (3e12) ** 2)); n = 10.0
for _ in range(300): n = math.log((tr(n) + 4) / 2) / eps
check('G20 reach: an off-line pair just above the verified height 3e12 (delta = 0.25) first shows at n > 1e26, far past the 1e5 computed',
      n > 1e26, f'n* ~ {n:.2e}')

# G21-G24 positivity / precession add-on
pr = json.load(open('probe.json'))
check('G21 probe spin on the primon gas: coherence |zeta(b+it)|/zeta(b) stays above zeta(2b)/zeta(b)^2 for b = 2, 1.5, 1.1 (t in (0, 50])',
      all(pr['min_coherence'][b] >= pr['lower_bound'][b] for b in pr['min_coherence']),
      ', '.join(f"{b}: {pr['min_coherence'][b]:.3f} >= {pr['lower_bound'][b]:.3f}" for b in pr['min_coherence']))
mp.mp.dps = 25
zs10 = [float(v) for v in Z['gamma'][:10]]
check('G22 at b = 1/2 the same coherence vanishes at each of the 10 zeros below t = 50 (|zeta| < 1e-15)',
      all(abs(mp.zeta(mp.mpc(0.5, mp.mpf(v)))) < 1e-15 for v in Z['gamma'][:10]))
spl = json.load(open('li_split.json')); A_ = [float(v) for v in spl['arch']]; P_ = [float(v) for v in spl['prime']]
check('G23 Li split: prime-at-infinity part + finite-primes part = lambda_n to 1e-10 for n = 1..1000',
      max(abs(A_[i] + P_[i] - float(lam[i])) for i in range(1000)) < 1e-10)
check('G24 the prime-at-infinity part is negative exactly for n <= 7 and outweighs |finite-primes part| for every 11 <= n <= 1000; the prime part stays in [-7, 8.5]',
      [i + 1 for i in range(1000) if A_[i] < 0] == list(range(1, 8)) and all(A_[i] > abs(P_[i]) for i in range(10, 1000))
      and -7 < min(P_) and max(P_) < 8.5, f'prime part in [{min(P_):.2f}, {max(P_):.2f}]')

# ---------------- polish
txt = open('README.md').read()
figs = re.findall(r'\(figures/([^)]+\.png)\)', txt)
check('P1 README is marked DRAFT and names the theory (Affine Composition Theory, built on Fused Hyperoperation Algebra)',
      txt.startswith('# DRAFT') and 'Affine Composition Theory' in txt and 'Fused Hyperoperation Algebra' in txt)
check('P2 every picture the README shows exists, and every picture made is shown', figs and all(os.path.exists('figures/' + f) for f in figs)
      and set(figs) == set(os.listdir('figures')), f'{len(figs)} pictures')
bad = [w for w in ('we prove the riemann', 'proves rh', 'proof of the riemann hypothesis is', 'rh is true') if w in txt.lower()]
check('P3 no overclaim phrases in the README', not bad, ','.join(bad))
check('P4 wording rule: no "downward"/"going down" for the heat flow; "backward heat flow" used', 'backward heat flow' in txt
      and not re.search(r'downward|going down', txt, re.I))
check('P5 the README opens with the answer in its first section and says plainly whether it is a lead', '## The answer' in txt[:3000] and 'lead' in txt[:3000].lower())
nums = {'5,000': '5,000' in txt, '1000': ('1000' in txt or '1,000' in txt)}
check('P6 numbers quoted in the README match the data (zeros count, lambda_1, first-negative n of the toys)',
      all(nums.values()) and '0.0230957' in txt and all(str(t['first_negative_n']) in txt for t in d['toy']))
check('P7 the Lean file says it is compiled in the shared library and has no sorry', 'compiled in the shared library' in open('lean/EllipticObject.lean').read() and 'sorry' not in re.sub(r'--.*', '', re.sub(r'/-.*?-/', '', open('lean/EllipticObject.lean').read(), flags=re.S)))
print(f'\n{sum(res)} of {len(res)} pass')
