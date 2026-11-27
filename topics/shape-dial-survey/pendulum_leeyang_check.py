"""Pendulum picture -> Lee-Yang work: tests of the three ideas (Nov 27).

Dial gas: Z = sum_S z^|S| k^e(S) on the square-lattice torus = Ising with k = e^{4K};
zero field is z = k^-2 (ring: z = k^-1). Onsager's modulus k_O = 1/sinh^2(2K) = 4k/(k-1)^2.
Pendulum: m = E/E_top, reading G(u) = sn dn / cn, mirror G(u) G(K-u) = 1.
Run: python3 pendulum_leeyang_check.py  (writes pendulum_leeyang_check.json)
"""
import json, pickle, os
import numpy as np, mpmath as mp
mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
LY = os.path.join(HERE, '..', 'lee-yang')
res = []
def ok(name, cond, detail):
    res.append(dict(name=name, passed=bool(cond), detail=detail))
    print(('PASS ' if cond else 'FAIL ') + name + ' | ' + detail)

kO = lambda k: 4*k/(k-1)**2          # Onsager modulus from the dial k = e^{4K}
# ---------------- Idea 1: the dial gas transition is the separatrix ----------------
kc = 3 + 2*mp.sqrt(2); ka = 3 - 2*mp.sqrt(2)
Kc = mp.log(kc)/4
ok('L1 k_c = 3+2sqrt2 <=> sinh2K = 1 <=> k_O = 1 (m = 1); 3-2sqrt2 is K = -K_c',
   abs(mp.sinh(2*Kc)-1) < 1e-28 and abs(kO(kc)-1) < 1e-28 and abs(kO(ka)-1) < 1e-28
   and abs(mp.log(ka)/4 + Kc) < 1e-28,
   f'sinh2K_c-1 = {mp.nstr(mp.sinh(2*Kc)-1,3)}, k_O(k_c)-1 = {mp.nstr(kO(kc)-1,3)}, k_O(k_a)-1 = {mp.nstr(kO(ka)-1,3)}')

def onsager_ss(k):
    """Nearest-neighbour <ss> on the infinite lattice, with Onsager's elliptic integral
    written as (1+kappa) x the pendulum quarter period at m = kappa^2, kappa = min(k_O, 1/k_O)."""
    k = mp.mpf(k); K = mp.log(k)/4; ko = kO(k); kap = min(ko, 1/ko)
    Kpend = mp.ellipk(kap**2)                       # pendulum: time to the turning point
    Konsager = (1+kap)*Kpend                        # Landen
    return mp.coth(2*K)/2*(1 + 2/mp.pi*(2*mp.tanh(2*K)**2-1)*Konsager)

# Landen check: K(k1) = (1+k_O) K(k_O), k1 = 2 sinh2K / cosh^2 2K
gaps = []
for k in [1.5, 2, 4, 5.8, 6, 10, 50]:
    k = mp.mpf(k); K = mp.log(k)/4; k1 = 2*mp.sinh(2*K)/mp.cosh(2*K)**2
    ko = kO(k); kap = min(ko, 1/ko)
    gaps.append(abs(mp.ellipk(k1**2) - (1+kap)*mp.ellipk(kap**2)))
ok('L2a Onsager integral = (1 + k_O) x pendulum quarter period at m = k_O^2 (Landen), both sides of k_c',
   max(gaps) < 1e-25, f'max gap {mp.nstr(max(gaps),3)} over k = 1.5..50')

C = pickle.load(open(os.path.join(LY, 'torus_counts.pkl'), 'rb'))
def torus_ss(L, k):
    c = C[L].astype(object); N = L*L
    z = mp.mpf(k)**-2; Z = mp.mpf(0); Ze = mp.mpf(0)
    for mm in range(c.shape[0]):
        for e in range(c.shape[1]):
            if c[mm, e]:
                w = c[mm, e]*z**mm*mp.mpf(k)**e; Z += w; Ze += e*w
    return 2*(Ze/Z)/N - 1                           # <e>/N = (1 + <ss>)/2
rows = []
for k in [2, 3, 10, 20]:
    on = onsager_ss(k); tor = [torus_ss(L, k) for L in (4, 5, 6)]
    rows.append((k, float(on), [float(t) for t in tor]))
far = max(abs(r[2][2]-r[1]) for r in rows if r[0] != 3)
near = [abs(t-rows[1][1]) for t in rows[1][2]]
ok('L2b dial gas torus (z = k^-2) matches Onsager-via-pendulum <ss>: within 1e-3 at L=6 for k = 2, 10, 20; closing in at k = 3 (nearer k_c)',
   far < 1e-3 and near[0] > near[1] > near[2], '; '.join(f'k={r[0]}: Onsager {r[1]:.5f}, torus L=4,5,6 ' + ','.join(f'{t:.5f}' for t in r[2]) for r in rows))

# Dictionary: Yang's magnetisation, Kramers-Wannier, over-the-top period
errs = []
for k in [6, 8, 12, 40]:
    k = mp.mpf(k); K = mp.log(k)/4; m = kO(k)**2
    M = (1 - mp.sinh(2*K)**-4)**(mp.mpf(1)/8)
    Kd = mp.asinh(1/mp.sinh(2*K))/2                # Kramers-Wannier dual coupling
    errs += [abs(M**8 - (1-m)), abs(kO(mp.e**(4*Kd)) - 1/kO(k))]
ok('L3a ordered side: Yang M^8 = 1 - m; Kramers-Wannier dual sends m -> 1/m (swinging <-> over the top)',
   max(errs) < 1e-25, f'max err {mp.nstr(max(errs),3)}')
# Onsager's elliptic part is the same at m and 1/m; the pendulum's over-the-top period is 2K(1/m)/sqrt(m)
e2 = []
for m in [mp.mpf('0.3'), mp.mpf('0.8')]:
    k1 = lambda ko: 2*mp.sqrt(ko)/(1+ko)
    e2.append(abs(k1(mp.sqrt(m)) - k1(1/mp.sqrt(m))))
ok('L3b Onsager k1 = 2 sqrt(k_O)/(1+k_O) is the same at k_O and 1/k_O (one integral serves both phases)',
   max(e2) < 1e-28, f'max err {mp.nstr(max(e2),3)}')

# Anisotropic couplings run along the pendulum's imaginary time; isotropic = halfway
e3 = []
for ko in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('0.999')]:
    Kp = mp.ellipk(1-ko**2)                         # imaginary quarter period K'
    s = lambda u: (-1j*mp.ellipfun('sn', 1j*u, m=ko**2))
    for u in [Kp*0.13, Kp*0.37, Kp*0.5, Kp*0.81]:
        e3.append(abs(s(u)*s(Kp-u) - 1/ko))
    e3.append(abs(s(Kp/2)**2 - 1/ko))               # halfway: both couplings equal, sinh^2 2K = 1/k_O
u0 = mp.mpf('0.6'); tanerr = abs(-1j*mp.ellipfun('sn', 1j*u0, m=mp.mpf(1)) - mp.tan(u0))
ok('L4 couplings sinh2K1 = -i sn(iu, k_O), sinh2K2 = same at K\'-u: product 1/k_O; equal (isotropic) at u = K\'/2; tan at k_c',
   max(e3) < 1e-25 and tanerr < 1e-25, f'max err {mp.nstr(max(e3),3)}, tan err {mp.nstr(tanerr,3)}')

# Fisher zeros (complex k, zero field) against |k_O| = 1, the separatrix in complex energy
F = pickle.load(open(os.path.join(LY, 'fisher_zeros.pkl'), 'rb'))
fz = {}
for L in (4, 5, 6):
    z = np.array(F[L]); d = np.abs(np.abs(4*z/(z-1)**2) - 1)
    # compare: the same number of points scattered over the same region
    rng = np.random.default_rng(L); r = np.abs(z).max()
    rnd = rng.uniform(-r, r, (4000,)) + 1j*rng.uniform(-r, r, (4000,))
    dr = np.abs(np.abs(4*rnd/(rnd-1)**2) - 1)
    fz[L] = (float(np.median(d)), float(np.median(dr)), len(z))
med = [fz[L][0] for L in (4, 5, 6)]
ok('L5 dial gas Fisher zeros close in on |k_O| = 1 (|m| = 1, the separatrix in complex energy) as L grows, well inside random scatter',
   med[0] > med[1] > med[2] and all(v[0] < 0.5*v[1] for v in fz.values()),
   '; '.join(f'L={L}: median ||k_O|-1| {v[0]:.3f} vs random {v[1]:.3f} ({v[2]} zeros)' for L, v in fz.items()))

# ---------------- Idea 2: unit-size lines ----------------
# Pendulum reading in complex time: where is |G| = 1, and where are the zeros of G^{2N} + 1?
# G = sn dn / cn has periods 2K, 2iK' and order 2 (poles at K and iK'), so G^{2N} = -1 has 4N roots per cell.
def halfway_roots(m, N=3):
    K = mp.ellipk(m); Kp = mp.ellipk(1-m)
    G = lambda u: mp.ellipfun('sn', u, m=m)*mp.ellipfun('dn', u, m=m)/mp.ellipfun('cn', u, m=m)
    roots = set()
    for j in range(2*N):
        t = mp.exp(1j*mp.pi*(2*j+1)/(2*N))
        for x in np.linspace(-0.95, 0.95, 10):
            for y in np.linspace(0.03, 1.97, 10):
                try:
                    r = mp.findroot(lambda u: G(u)-t, K*x+1j*Kp*y)
                except Exception:
                    continue
                if abs(G(r)-t) > 1e-15: continue
                xr = (float(mp.re(r)/K)+1) % 2 - 1; yr = float(mp.im(r)/Kp) % 2
                roots.add((round(xr, 6), round(yr, 6)))
    on = [abs(abs(x)-0.5) < 1e-6 or abs(y % 1 - 0.5) < 1e-6 for x, y in roots]
    # size check on the two kinds of halfway line
    sz = max(max(abs(abs(G(K/2+1j*Kp*y))-1), abs(abs(G(K*x+1j*Kp/2))-1)) for x, y in [(0.1, 0.3), (0.7, 1.4), (-0.4, 0.9)])
    return len(roots), sum(on), float(sz)
hw = {float(m): halfway_roots(m) for m in [mp.mpf('0.2'), mp.mpf('0.7'), mp.mpf('0.95')]}
ok('L6a pendulum: |G| = 1 on the halfway lines Re u = +-K/2 (real time) and Im u = K\'/2 (imaginary time); all 4N zeros of G^2N + 1 lie on them (N = 3)',
   all(v[0] == 12 and v[1] == 12 and v[2] < 1e-15 for v in hw.values()),
   '; '.join(f'm={m}: {v[0]} roots found, {v[1]} on halfway lines, max ||G|-1| {v[2]:.1e}' for m, v in hw.items()))
# 1D Ising ring: Z = det^{N/2} (g^N + g^-N), g = lambda+/sqrt(det): the same two-term mirror form
def ring_poly(N, k):
    """Z(w) for the lattice-gas ring, w = z k (so the Lee-Yang circle is |w| = 1). Coeffs high->low."""
    # transfer matrix in w: T = [[1, sqrt(w/k)], [sqrt(w/k), w]] up to the factor; use exact polynomial algebra
    # weights z^n k^{n_i n_j} = w^n k^{n_i n_j - n}; on the ring sum(n_i n_j) - sum(n) = -(number of occupied runs)
    from math import comb
    P = np.zeros(N+1)
    P[0] += 1; P[N] += 1                            # empty ring, full ring (zero runs)
    for n in range(1, N):                           # n occupied, r runs: count = N/r C(n-1,r-1) C(N-n-1,r-1)
        for r in range(1, min(n, N-n)+1):
            P[n] += N/r*comb(n-1, r-1)*comb(N-n-1, r-1)*float(k)**(-r)
    return P[::-1]
e5 = []; edge = []
for N in [8, 16, 32]:
    for k in [2, 5, 20]:
        w = np.array([complex(r) for r in mp.polyroots([mp.mpf(float(c)) for c in ring_poly(N, k)], maxsteps=200, extraprec=200)]); th = np.angle(w)
        x = np.sort(np.cos(th/2)**2)                # cos^2(theta/2), sign-free
        c = np.sqrt(1-1/k); psi = np.pi*(2*np.arange(N)+1)/(2*N)
        xp = np.sort((c*np.cos(psi))**2)
        e5.append(max(np.max(np.abs(np.abs(w)-1)), np.max(np.abs(x-xp))))
ok('L6b 1D Ising ring: every Lee-Yang zero has |w| = 1 and cos(theta/2) = sqrt(1-1/k) cos(pi(2j+1)/2N) (equal-time samples of a swing)',
   max(e5) < 1e-10, f'max err {max(e5):.1e} over N = 8,16,32 and k = 2,5,20')

# ---------------- Idea 3: crowding at the arc's end ----------------
# Edge exponent: with x = cos(theta/2) = c cos(psi), the zeros' density in x is 1/sqrt(c^2 - x^2): exponent -1/2
N = 2000; c = np.sqrt(1-1/5); psi = np.pi*(2*np.arange(N)+1)/(2*N); x = c*np.cos(psi)
dd = np.array([1e-4, 1e-3, 1e-2]); cnt = np.array([np.sum(x > c-d) for d in dd])
slope = np.polyfit(np.log(dd), np.log(cnt), 1)[0]
ok('L7a 1D crowding: count within d of the arc end grows like d^{1/2} (density exponent sigma = -1/2, the turning-point linger)',
   abs(slope-0.5) < 0.03, f'fitted count exponent {slope:.3f}')
ok('L7b the swing size is set by temperature, not system size: c = sqrt(1-1/k) -> 1 only as k -> infinity; N only refines the sampling',
   abs(np.sqrt(1-1/20)-np.sqrt(1-1/20)) == 0 and np.sqrt(1-1/2) < np.sqrt(1-1/20) < 1, 'c(k=2) = 0.707, c(k=20) = 0.975, c independent of N (L6b)')
# antiferro k < 1: m = 1 - 1/k < 0 (the Wick band) and the zeros leave the circle
w = np.roots(ring_poly(16, 0.5))
ok('L7c antiferro k = 1/2 (m = 1 - 1/k < 0): ring zeros leave |w| = 1 and lie on the negative real axis',
   np.max(np.abs(w.imag)) < 1e-6 and np.all(w.real < 0) and np.max(np.abs(np.abs(w)-1)) > 0.1,
   f'max |Im w| {np.max(np.abs(w.imag)):.1e}, |w| range {np.abs(w).min():.3f}..{np.abs(w).max():.3f}')
# p-ball single spin: zeros evenly spaced (spacing -> pi), no crowding anywhere
mp.mp.dps = 90
p, a = 4, 3
cf = [mp.gamma(a+1)*mp.gamma(mp.mpf(2*j+1)/p)/(p*mp.factorial(2*j)*mp.gamma(mp.mpf(2*j+1)/p+a+1)) for j in range(220)]
def Gs(s):
    w = -mp.mpf(s)**2; acc = mp.mpf(0)
    for cc in reversed(cf): acc = acc*w + cc
    return acc
grid = [mp.mpf(i)/10 for i in range(1, 601)]; vals = [Gs(s) for s in grid]
zs = [mp.findroot(Gs, (grid[i], grid[i+1]), solver='anderson', tol=mp.mpf(10)**-40) for i in range(len(grid)-1) if vals[i]*vals[i+1] < 0]
sp = [float(zs[i+1]-zs[i]) for i in range(len(zs)-1)]
ok('L7d p-ball single spin (p = 4, a = 3): zeros evenly spaced, spacing closing on pi from above; no arc, no crowding',
   len(zs) > 10 and abs(sp[-1]-np.pi) < 0.03 and all(sp[i+1] < sp[i] for i in range(len(sp)-4, len(sp)-1)), f'{len(zs)} zeros in (0, 60); last spacings ' + ', '.join(f'{v:.4f}' for v in sp[-3:]))
mp.mp.dps = 30

print(f'\n{sum(r["passed"] for r in res)} of {len(res)} pass')
json.dump(dict(results=res, onsager_rows=rows, fisher=fz), open(os.path.join(HERE, 'pendulum_leeyang_check.json'), 'w'), indent=1, default=str)
