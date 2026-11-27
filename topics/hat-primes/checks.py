"""Checks for the hat-and-primes DRAFT. Needs numpy, sympy, matplotlib (for point-in-polygon) only.
Prints PASS/FAIL per check and a final count."""
import json, time, numpy as np, sympy as sp
from matplotlib.path import Path
import hatgen as h

t0 = time.time(); results = []
def check(name, ok, detail=''):
    results.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name + (f'  [{detail}]' if detail else ''))

# ---------- the generator builds real hat tilings ----------
L = h.levels(4)
counts = [len(h.hats_of(q[0])) for q in L]
check('1. hats per H supertile, levels 0-4, are 4, 25, 169, 1156, 7921', counts == [4, 25, 169, 1156, 7921], str(counts))
patch = h.constructPatch(*L[4])
KEEP = set(h.NEW_H_CH + h.NEW_P_CH + h.NEW_F_CH + h.NEW_T_CH)
Ts = [T for i, (T, g) in enumerate(patch.children) if i in KEEP for T, k in g.hats(T)]
o = h.transPt(Ts[0], h.HAT[0])
corners = np.array([h.transPt(T, p) - o for T in Ts for p in h.HAT])
xy = np.array([h.to_hex(p) for p in corners])
check('2. every corner of all 129,314 hats lies on one triangular lattice', len(Ts) == 129314 and np.abs(xy - np.round(xy)).max() < 1e-6,
      f'{len(Ts)} hats, max offset {np.abs(xy - np.round(xy)).max():.1e}')
# no overlaps: sample the reference hat at every point of the fine lattice (1/6)L strictly inside it. That lattice is
# carried to itself by every motion that places a hat (60-degree turns, reflections, lattice shifts), so all hats
# are sampled on one common grid, and two hats overlap exactly when they share a sample point.
def strictly_inside(poly, pts, eps=1e-9):
    ok = Path(poly).contains_points(pts)
    a, b = poly, np.roll(poly, -1, axis=0)
    for (x1, y1), (x2, y2) in zip(a, b):
        dx, dy = x2 - x1, y2 - y1
        t = np.clip(((pts[:, 0] - x1)*dx + (pts[:, 1] - y1)*dy)/(dx*dx + dy*dy), 0, 1)
        d = np.hypot(pts[:, 0] - (x1 + t*dx), pts[:, 1] - (y1 + t*dy))
        ok &= d > eps
    return ok
g = np.mgrid[-36:60, -36:60].reshape(2, -1).T / 6.0
pts = np.array([h.hexPt(a, b) for a, b in g])
P = pts[strictly_inside(np.array(h.HAT), pts)]
def place(T): return np.stack([T[0]*P[:, 0] + T[1]*P[:, 1] + T[2], T[3]*P[:, 0] + T[4]*P[:, 1] + T[5]], 1)
def shared(Tlist):
    k = np.round(np.vstack([place(T) for T in Tlist]) * 1200).astype(np.int64)
    return len(k) - np.unique(k, axis=0).shape[0]
sh = shared(Ts)
sh_scaffold = shared([T for T, k in h.hats_of(patch)])   # all 29 patch children, scaffolding included
check('3. no two hats overlap (no shared sample point); the same test does catch the one known overlap in the scaffolding patch',
      sh == 0 and sh_scaffold > 0, f'{len(P)} samples per hat; shared: {sh} in the 129,314 hats, {sh_scaffold} in the scaffolding patch')

# ---------- the substitution and its numbers ----------
KINDS = ['H', 'T', 'P', 'F']
def row(chs): return [sum(1 for c in chs if h.PATCH_KINDS[c] == k) for k in KINDS]
A = sp.Matrix([row(h.NEW_H_CH), row(h.NEW_T_CH), row(h.NEW_P_CH), row(h.NEW_F_CH)])
x = sp.symbols('x')
check('4. substitution matrix [[3,1,3,3],[1,0,0,0],[2,0,1,2],[2,0,1,3]] (matches Schirmann et al. 2023)',
      A.tolist() == [[3, 1, 3, 3], [1, 0, 0, 0], [2, 0, 1, 2], [2, 0, 1, 3]])
cp = sp.expand(A.charpoly(x).as_expr())
check('5. characteristic polynomial = (x-1)(x+1)(x^2-7x+1); largest root phi^4', sp.expand(cp - (x - 1)*(x + 1)*(x**2 - 7*x + 1)) == 0
      and sp.simplify(((1 + sp.sqrt(5))/2)**4 - (7 + 3*sp.sqrt(5))/2) == 0)
hv = sp.Matrix([4, 1, 2, 2]); eH = sp.Matrix([[1, 0, 0, 0]])
w = (A**3 - 8*A**2 + 8*A - sp.eye(4)) * hv
fibsq = [sp.fibonacci(2*n + 3)**2 for n in range(3)]
first = [(eH * A**n * hv)[0] for n in range(3)]
check('6. hats in a level-n H supertile = F(2n+3)^2 for every n (proof: both obey s(n+3)=8s(n+2)-8s(n+1)+s(n), same start)',
      (A*w == -w) and (eH*w)[0] == 0 and first == fibsq,
      f'(eH w = {(eH*w)[0]}; first terms {first})')
lam = (7 + 3*sp.sqrt(5))/2
u = [v for (val, m, vecs) in A.T.eigenvects() if sp.simplify(val - lam) == 0 for v in vecs][0]; u = u / sum(u)
refl = u[0]; tot = (sp.Matrix([[4, 1, 2, 2]]) * u)[0]
ratio = sp.radsimp((tot - refl) / refl)
r6 = [(eH * A**n * hv)[0] for n in (8,)][0]; f6 = [(eH * A**n * sp.Matrix([1, 0, 0, 0]))[0] for n in (8,)][0]
check('7. unreflected : reflected hats -> phi^4 exactly (6.8541)', sp.simplify(ratio - lam) == 0 and abs(float((r6 - f6)/f6) - float(lam)) < 1e-5,
      f'level 8: {float((r6 - f6)/f6):.7f}')
check('8. metatile shares H, T, P, F = 1/3, 7/6-sqrt5/2, sqrt5-2, 3/2-sqrt5/2 (all in Q(sqrt5))',
      all(sp.simplify(a - b) == 0 for a, b in zip(u, [sp.Rational(1, 3), sp.Rational(7, 6) - sp.sqrt(5)/2, sp.sqrt(5) - 2, sp.Rational(3, 2) - sp.sqrt(5)/2])))
# periods mod p
def period(p):
    Ap = np.array(A.tolist(), dtype=np.int64) % p; X = Ap.copy(); I = np.eye(4, dtype=np.int64)
    for n in range(1, 4*p*p):
        if (X == I).all(): return n
        X = (X @ Ap) % p
bad = []
for p in sp.primerange(7, 400):
    T = period(p)
    ok = ((p - 1) % T == 0) if p % 5 in (1, 4) else ((p + 1) % T == 0)
    if not ok: bad.append((p, T))
check('9. counts mod p repeat with period dividing p-1 (p = +-1 mod 5) or p+1 (p = +-2 mod 5), all primes 7..397; p=5 period 10',
      not bad and period(5) == 10 and period(2) == 6 and period(3) == 6, f'exceptions {bad}')

# ---------- the hat on the Eisenstein lattice ----------
anc = np.array([h.transPt(T, h.HAT[0]) for T in Ts])
cen = np.array([h.hat_polygon(T).mean(0) for T in Ts]); c0 = cen.mean(0)
o2 = anc[np.argmin(((anc - c0)**2).sum(1))]
A0 = np.round(np.array([h.to_hex(p - o2) for p in anc])).astype(np.int64)
X2, Y2 = 2*A0[:, 0] + A0[:, 1], A0[:, 1] - A0[:, 0]
check('10. every hat anchor is a hexagon centre of the kite grid (a multiple of 2(1+w6))', (X2 % 6 == 0).all() and (Y2 % 6 == 0).all())
Ah = np.stack([X2 // 6, Y2 // 6], 1)
nrm = lambda P: P[:, 0]**2 + P[:, 0]*P[:, 1] + P[:, 1]**2
def classes(P, p, r):
    if r is None: return (P[:, 0] % p)*p + (P[:, 1] % p), p*p
    return (P[:, 0] + P[:, 1]*r) % p, p
tests = [(2, None), (3, 2), (5, None)]
for p in sp.primerange(7, 44):
    if p % 3 == 1: tests += [(p, r) for r in range(p) if (r*r - r + 1) % p == 0]
worst = 0; worst_chi = 0
for p, r in tests:
    c, n = classes(Ah, p, r); cnt = np.bincount(c, minlength=n); e = len(Ah)/n
    worst = max(worst, np.abs(cnt/e - 1).max()); worst_chi = max(worst_chi, ((cnt - e)**2/e).sum()/(n - 1))
check('11. hat anchors are even over the residue classes of every Eisenstein prime of norm <= 43 (within 2%, chi2/dof < 1)',
      worst < 0.02 and worst_chi < 1, f'worst class off by {worst:.4f}, worst chi2/dof {worst_chi:.3f}')
# prime and square-free shares, anchors vs empty hexagon centres in the same region
Aset = set(map(tuple, Ah.tolist()))
nb = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
Rg = np.array(sorted({(a + dx, b + dy) for a, b in Aset for dx, dy in nb})); Rg = Rg[(Rg != 0).any(1)]
NM = int(nrm(Rg).max()) + 10
isp = np.ones(NM + 1, bool); isp[:2] = False
for i in range(2, int(NM**0.5) + 1):
    if isp[i]: isp[i*i::i] = False
def prime_mask(P):
    n = nrm(P); q = np.round(np.sqrt(n)).astype(np.int64); qq = np.maximum(q, 1)
    return isp[n] | ((q*q == n) & isp[q] & (q % 3 == 2) & (P[:, 0] % qq == 0) & (P[:, 1] % qq == 0))
def sqfree_mask(P):
    bad = np.zeros(len(P), bool)
    for p in sp.primerange(2, int(NM**0.5) + 2):
        if p == 3: bad |= (P[:, 0] % 3 == 0) & (P[:, 1] % 3 == 0)
        elif p % 3 == 2: bad |= (P[:, 0] % (p*p) == 0) & (P[:, 1] % (p*p) == 0)
        else:
            for r in [t for t in range(p) if (t*t - t + 1) % p == 0]:
                q = p*p; r2 = (r - (r*r - r + 1)*pow(2*r - 1, -1, q)) % q
                bad |= (P[:, 0] + P[:, 1]*r2) % q == 0
    return ~bad
hole = np.array([tuple(p) not in Aset for p in Rg.tolist()])
pm, sm = prime_mask(Rg), sqfree_mask(Rg)
L2 = float(sum((1 if n % 3 == 1 else -1 if n % 3 == 2 else 0)/n**2 for n in range(1, 300001)))
inv_zK2 = 1/((np.pi**2/6)*L2)
check('12. hats sit on primes exactly as often as empty hexagon centres do (prime share within 0.3 points)',
      abs(pm[~hole].mean() - pm[hole].mean()) < 0.003, f'anchors {pm[~hole].mean():.4f}, empty {pm[hole].mean():.4f}')
check('13. square-free share of hat anchors = 1/zeta_K(2) = 0.77809 for K = Q(sqrt-3) (within 0.002)',
      abs(sm[~hole].mean() - inv_zK2) < 0.002, f'{sm[~hole].mean():.5f} vs {inv_zK2:.5f}')
# the union (k = -1) composite of the local removal rates gives the same number
Az = Ah[(Ah != 0).any(1)]; rates = []
for p in sp.primerange(2, int(nrm(Az).max()**0.5) + 2):
    if p == 3: rates.append(((Az[:, 0] % 3 == 0) & (Az[:, 1] % 3 == 0)).mean())
    elif p % 3 == 2: rates.append(((Az[:, 0] % (p*p) == 0) & (Az[:, 1] % (p*p) == 0)).mean())
    else:
        for r in [t for t in range(p) if (t*t - t + 1) % p == 0]:
            q = p*p; r2 = (r - (r*r - r + 1)*pow(2*r - 1, -1, q)) % q
            rates.append(((Az[:, 0] + Az[:, 1]*r2) % q == 0).mean())
uc = 0.0
for a in rates: uc = uc + a - uc*a
meas = sqfree_mask(Az).mean()
check('14. 1 - (union composite a o_{-1} b = a + b - ab of the local removal rates) = measured square-free share (within 0.002)',
      abs((1 - uc) - meas) < 0.002, f'{1 - uc:.5f} vs {meas:.5f}')

# ---------- spectra ----------
N = 10**6
isq = np.ones(N + 1, bool)
for p in sp.primerange(2, int(N**0.5) + 1): isq[p*p::p*p] = False
sqf = np.nonzero(isq[1:])[0] + 1
def formula(q):
    f = sp.factorint(q)
    if any(e >= 3 for e in f.values()): return 0.0
    v = 6/np.pi**2
    for p in f: v *= -1.0/(p*p - 1)
    return v
errs = [abs(np.exp(2j*np.pi*sqf*a/q).sum().real/N - formula(q)) for a, q in [(1, 2), (1, 3), (1, 4), (1, 6), (2, 9), (5, 12), (1, 25), (7, 36), (1, 8), (1, 27), (1, 7)]]
check('15. square-free integers: peak amplitude at a/q = (6/pi^2) * prod over p | q of -1/(p^2-1), 0 if q not cube-free (Baake-Moody-Pleasants)',
      max(errs) < 1e-4, f'max error {max(errs):.1e} at N = 10^6')
z = np.load('../primes/zeta_zeros_1000.npy')
xs = np.linspace(0.3, 3.6, 8000); Ssum = np.cos(np.outer(xs, z)).sum(1)
mins = [xs[i] for i in range(1, len(xs) - 1) if Ssum[i] < Ssum[i-1] and Ssum[i] < Ssum[i+1] and Ssum[i] < -60]
ppw = [m for m in range(2, 40) if len(sp.factorint(m)) == 1]
check('16. the first 1,000 zeta zeros: every deep dip of sum cos(gamma x) on [0.3, 3.6] sits at log of a prime power (within 2e-4)',
      len(mins) >= 12 and all(min(abs(np.log(m) - x0) for m in ppw) < 2e-4 for x0 in mins), f'{len(mins)} dips')
# Bragg peaks: same height at 19k and 129k hats, positions in Q(sqrt5)/8
def hexc(Tl):
    a = np.array([h.transPt(T, h.HAT[0]) for T in Tl]); oo = a[0]
    q = np.round(np.array([h.to_hex(p - oo) for p in a])).astype(np.int64)
    return np.stack([(2*q[:, 0] + q[:, 1])//6, (q[:, 1] - q[:, 0])//6], 1)
p3 = h.constructPatch(*L[3])
T3 = [T for i, (T, g) in enumerate(p3.children) if i in KEEP for T, k in g.hats(T)]
U3, U4 = hexc(T3), Ah
def Iat(U, q): return float(np.abs(np.exp(-2j*np.pi*(U[:, 0]*q[0] + U[:, 1]*q[1])).sum())**2/len(U)**2)
s5 = 5**0.5
q_pred = [((4 - s5)/4, (9 - s5)/8), ((1 + s5)/8, s5/4)]
hs3 = [Iat(U3, q) for q in q_pred]; hs4 = [Iat(U4, q) for q in q_pred]
check('17. Bragg spots at the exact Q(sqrt5) points ((4-sqrt5)/4, (9-sqrt5)/8) and ((1+sqrt5)/8, sqrt5/4): same height (within 2%) at 19k and 129k hats',
      all(b > 0.015 and abs(a/b - 1) < 0.02 for a, b in zip(hs3, hs4)), f'{[round(v, 5) for v in hs3]} vs {[round(v, 5) for v in hs4]}')
check('18. log 2, log 3, log 5, log 7 have no integer relation with coefficients up to 6 (they are independent over Q)',
      not any(abs(a*np.log(2) + b*np.log(3) + c*np.log(5) + d*np.log(7)) < 1e-9
              for a in range(-6, 7) for b in range(-6, 7) for c in range(-6, 7) for d in range(-6, 7) if (a, b, c, d) != (0, 0, 0, 0)))
print(f'{sum(results)} of {len(results)} pass ({time.time() - t0:.0f} s)')
