"""Checks for the similarity-ratio survey (Nov 26, 2026). Run from this folder: python3 checks.py > checks_output.txt
Each line prints PASS/FAIL with a label; the last line is 'X of Y pass'. numpy, scipy, sympy, mpmath only."""
import math, time, random
from collections import Counter
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import mpmath as mp

results = []
def check(label, ok, note=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {label}" + (f"  [{note}]" if note else ""))

def moran(rs, s0=0.5):
    """Similarity dimension: the real root of sum r_i^s = 1 (Moran 1946). An output of the ratios."""
    from scipy.optimize import brentq
    return brentq(lambda s: sum(r**s for r in rs) - 1, 1e-9, 60)

lg = math.log
# ---------------------------------------------------------------- C1: two pieces, uniform ratio r
print("== C1  two pieces at ratio r: segment, Cantor, half box, overlap")
grid = np.linspace(0.02, 0.98, 49)
check("A1 Moran root = log 2 / log(1/r) along the whole dial (49 settings)",
      max(abs(moran([r, r]) - lg(2)/lg(1/r)) for r in grid) < 1e-10)
check("A2 stops: r = 1/4 -> D = 1/2 (half box), r = 1/3 -> log2/log3 (Cantor), r = 1/2 -> 1 (segment)",
      abs(moran([.25, .25]) - .5) < 1e-12 and abs(moran([1/3, 1/3]) - lg(2)/lg(3)) < 1e-12 and abs(moran([.5, .5]) - 1) < 1e-12,
      f"{moran([1/3,1/3]):.6f}")
gaps = [(r, 1 - 2*r) for r in (0.25, 1/3, 0.45, 0.5, 0.55, 0.618)]
check("A3 the two images [0, r] and [1-r, 1] have gap 1 - 2r: open below 1/2, touching at 1/2, overlapping above",
      all((g > 0) == (r < 0.5) and (abs(g) < 1e-15) == (r == 0.5) for r, g in gaps))

def boxdim(r, n=14):
    """Box-count slope of the level-n pieces of {x -> r x, x -> r x + 1 - r}, at the scale of a level-n piece."""
    lefts = np.array([0.0])
    for _ in range(n):
        lefts = np.concatenate([r*lefts, r*lefts + 1 - r])
    def count(eps):
        # number of eps-grid boxes met by the union of intervals [x, x + r^n]
        a = np.floor(lefts/eps).astype(np.int64); b = np.floor((lefts + r**n)/eps).astype(np.int64)
        boxes = set()
        for lo, hi in zip(a, b):
            boxes.update(range(lo, hi + 1))
        return len(boxes)
    e1, e2 = r**(n - 4), r**n
    return (lg(count(e2)) - lg(count(e1)))/lg(e1/e2)
bd = {r: boxdim(r) for r in (0.25, 0.4)}
check("A4 below 1/2 the box-count slope matches D (open set condition holds; Hutchinson 1981)",
      all(abs(bd[r] - lg(2)/lg(1/r)) < 0.02 for r in bd), ", ".join(f"r={r}: {bd[r]:.4f} vs {lg(2)/lg(1/r):.4f}" for r in bd))
b6 = boxdim(0.6, 16)
check("A5 above 1/2 the set is the whole segment: box slope 1, while log2/log(1/r) = 1.357 > 1 (OSC fails)",
      abs(b6 - 1) < 0.02 and lg(2)/lg(1/0.6) > 1.3, f"slope {b6:.4f}")

print("== C1b  the classic stops with N pieces")
check("B1 Sierpinski (3 at 1/2) log3/log2, Koch (4 at 1/3) log4/log3, n-cube (2^n at 1/2) n for n = 1..4",
      abs(moran([.5]*3) - lg(3)/lg(2)) < 1e-12 and abs(moran([1/3]*4) - lg(4)/lg(3)) < 1e-12
      and all(abs(moran([.5]*2**n) - n) < 1e-10 for n in range(1, 5)))
check("B2 turning r with N = 2^n fixed slides the cube's D = n log2/log(1/r) (n = 3, 20 settings)",
      max(abs(moran([r]*8) - 3*lg(2)/lg(1/r)) for r in np.linspace(.05, .9, 20)) < 1e-9)

# ---------------------------------------------------------------- C2: box layers along D
print("== C2  box layers C(D, k) slide as D moves")
def layer(D, k): return sp.binomial(sp.Rational(D), k) if not isinstance(D, sp.Basic) else sp.binomial(D, k)
half = [sp.binomial(sp.Rational(1, 2), k) for k in range(6)]
check("L1 D = 1/2 layers 1, 1/2, -1/8, 1/16, -5/128, 7/256 (FractalBoxes.half_layers)",
      half == [1, sp.Rational(1, 2), sp.Rational(-1, 8), sp.Rational(1, 16), sp.Rational(-5, 128), sp.Rational(7, 256)])
check("L2 D = 3 layers 1, 3, 3, 1, 0 (cube_layers)", [sp.binomial(3, k) for k in range(5)] == [1, 3, 3, 1, 0])
mp.mp.dps = 30
t0 = time.time()
err = 0
for D in (0.25, 0.5, 0.631, 1.585, 2.5):
    for a in (mp.mpf('0.5'), mp.mpf('0.9')):
        s = mp.nsum(lambda k: mp.binomial(D, k)*a**k, [0, mp.inf])
        err = max(err, abs(s - (1 + a)**D))
tpol1 = time.time() - t0
check("L3 POLISH (30 digits): sum_k C(D,k) a^k = (1+a)^D at D = 1/4, 1/2, 0.631, 1.585, 5/2 and a = 0.5, 0.9",
      err < mp.mpf(10)**-25, f"max error {mp.nstr(err, 3)}, {tpol1:.2f} s")
# the Moran box: side 1/r, so (1 + a)^D = (1/r)^D = N exactly
check("L4 Moran box: a box of dimension D with side 1/r holds (1/r)^D = N pieces, for N = 2, 3, 8 and many r",
      max(abs((1/r)**moran([r]*N) - N) for N in (2, 3, 8) for r in np.linspace(.05, .9, 15)) < 1e-9)
def terms(D, a, ks): return [abs(float(mp.binomial(D, k)*a**k)) for k in ks]
tc = terms(lg(2)/lg(3), 2, (20, 40, 80))      # Cantor: a = 1/r - 1 = 2
t6 = terms(lg(2)/lg(1/0.6), 1/0.6 - 1, (20, 40, 80))   # r = 0.6: a = 2/3
check("L5 the Moran box's layer series converges only for r >= 1/2 (a = 1/r - 1 <= 1): Cantor terms grow, r = 0.6 terms shrink",
      tc[0] < tc[1] < tc[2] and t6[0] > t6[1] > t6[2], f"Cantor |term| at k=20,40,80: {tc[0]:.2e}, {tc[1]:.2e}, {tc[2]:.2e}")
check("L6 joining at a common ratio adds D: two half boxes (2 at 1/4) make 4 at 1/4, D = 1 (half_join_stick); three 1/3-boxes (2 at 1/8) make 8 at 1/8, D = 1 (third_join_stick)",
      abs(moran([.25]*4) - 1) < 1e-12 and abs(moran([1/8]*2) - 1/3) < 1e-12 and abs(moran([1/8]*8) - 1) < 1e-12)
check("L7 D = 1/2 joined with D = 3/2 at a common ratio gives D = 2 (half_and_three_halves): 2 at 1/4 times 8 at 1/4 = 16 at 1/4",
      abs(moran([.25]*2) + moran([.25]*8) - moran([.25]*16)) < 1e-12 and abs(moran([.25]*16) - 2) < 1e-12)

# ---------------------------------------------------------------- C3: the weight r^D and the k dial
print("== C3  the weight r^s: k at a fixed s, frozen at the true D, odds when one piece turns")
t, u, w = sp.symbols('t u w')
circ = lambda x, y, k: x + y + k*x*y
check("W1 w (t o_w u) = (w t) o_1 (w u): a weight is the homothety k = w (circK_scale, ZoomDial.zoom_is_twist)",
      sp.expand(w*circ(t, u, w) - circ(w*t, w*u, 1)) == 0)
check("W2 uniform turning: at the Moran dimension every piece weighs r^D = 1/N, whatever r (N = 2, 3, 5; 30 settings)",
      max(abs(r**moran([r]*N) - 1/N) for N in (2, 3, 5) for r in np.linspace(.03, .95, 30)) < 1e-10)
check("W3 at a fixed trial exponent s the weight r^s does move (s = 1/2: 0.316 at r = 0.1, 0.707 at r = 0.5): that is k moving",
      abs(0.1**0.5 - 0.316) < 1e-3 and abs(0.5**0.5 - 0.7071) < 1e-4)
oks = True
for r1 in (0.1, 0.2, 0.3, 0.45):
    r2 = 0.4; D = moran([r1, r2]); p1, p2 = r1**D, r2**D; odds = p1/p2
    oks &= abs(p1 + p2 - 1) < 1e-12 and abs(odds - (r1/r2)**D) < 1e-12 and abs(odds/(1 + odds) - p1) < 1e-12
check("W4 one piece turned alone: the weights r_i^D stay on the simplex (sum 1); their odds are (r1/r2)^D and the chance t/(1+t) is the weight",
      oks)
oks = all(abs(mp.log(sum(mp.e**(moran(rs)*mp.log(r)) for r in rs))/moran(rs)) < 1e-12
          for rs in ([.25, .25], [.3, .5], [.1, .2, .4]))
check("W5 Bowen's equation: the soft maximum of log r_i at temperature 1/D is exactly 0 (free energy zero)", oks)
x = lambda s, rs: sum(r**s for r in rs)
check("W6 depth-n words weigh x(s)^n; the k = -1 word series 1/(1 - x) converges for s > D and hits its wall at s = D (zoom-dial)",
      abs(x(0.7, [1/3, 1/3])**3 - sum((a*b*c)**0.7 for a in [1/3]*2 for b in [1/3]*2 for c in [1/3]*2)) < 1e-14
      and x(0.7, [1/3]*2) < 1 and abs(x(lg(2)/lg(3), [1/3]*2) - 1) < 1e-14)

# ---------------------------------------------------------------- C4: Kraft and Renyi
print("== C4  Kraft's sum is Moran's wall at D = 1; Renyi entropy becomes Renyi dimension")
random.seed(3)
def rand_tree(depth):
    if depth == 0 or random.random() < 0.3: return None
    return (rand_tree(depth - 1), rand_tree(depth - 1))
def tree_sum(tr, r1, r2, s):
    if tr is None: return 1.0
    return r1**s*tree_sum(tr[0], r1, r2, s) + r2**s*tree_sum(tr[1], r1, r2, s)
def kraft(tr):
    if tr is None: return Fr(1)
    return (kraft(tr[0]) + kraft(tr[1]))/2
trees = [rand_tree(9) for _ in range(40)]
check("K1 every full binary tree has Kraft sum exactly 1 (InformationFull.kraft_one), 40 random trees, exact fractions",
      all(kraft(tr) == 1 for tr in trees))
ok, off = True, True
for r1, r2 in ((.5, .5), (.2, .6), (.3, .3), (.05, .9), (.7, .45)):
    D = moran([r1, r2])
    ok &= all(abs(tree_sum(tr, r1, r2, D) - 1) < 1e-10 for tr in trees)
    off &= any(abs(tree_sum(tr, r1, r2, D + 0.1) - 1) > 1e-3 for tr in trees)
check("K2 with branch ratios (r1, r2) at every node, every full tree weighs 1 at the Moran D (and not at D + 0.1): Kraft slides to Moran",
      ok and off, "5 ratio pairs x 40 trees; the induction is kraft_one's with 1/2 replaced by r_i^D")
check("K3 c+1 branches at ratio 1/(c+1) give D = 1 for every c (dary_full's tree); at ratio r they give log(c+1)/log(1/r)",
      all(abs(moran([1/(c+1)]*(c+1)) - 1) < 1e-10 for c in range(1, 7)) and abs(moran([0.2]*3) - lg(3)/lg(5)) < 1e-10)
p = np.array([0.5, 0.3, 0.2]); r = 0.25
def renyi(p, q): return float(np.sum(p*np.log(1/p))) if q == 1 else float(np.log(np.sum(p**q))/(1 - q))
def tau(q, p, rs):
    from scipy.optimize import brentq
    return brentq(lambda T: float(np.sum(p**q*np.array(rs)**(-T))) - 1, -50, 50)
ok = all(abs(tau(q, p, [r]*3)/(q - 1) - renyi(p, q)/lg(1/r)) < 1e-9 for q in (0, 0.5, 2, 3))
check("K4 uniform ratio: the Renyi dimension D_q = tau(q)/(q-1) equals H_q(p)/log(1/r) (q = 0, 1/2, 2, 3)", ok)
rs = [0.2, 0.3, 0.45]; D = moran(rs); pn = np.array(rs)**D
check("K5 natural weights p_i = r_i^D make every D_q equal D (a monofractal), uneven ratios",
      all(abs(tau(q, pn, rs)/(q - 1) - D) < 1e-9 for q in (0, 0.5, 2, 4)))
pu = np.array([0.6, 0.3, 0.1]); h = 1e-5
D1 = float(np.sum(pu*np.log(pu))/np.sum(pu*np.log(rs)))
check("K6 unequal ratios (multifractal formalism, Cawley-Mauldin 1992): tau(0) = -D_Moran and tau'(1) = information dimension sum p log p / sum p log r",
      abs(tau(0, pu, rs) + D) < 1e-9 and abs((tau(1 + h, pu, rs) - tau(1 - h, pu, rs))/(2*h) - D1) < 1e-6, f"D_1 = {D1:.5f}, D_0 = {D:.5f}")

# ---------------------------------------------------------------- C5: Gentile ladder r = 1/(m+1)
print("== C5  Stefan-Boltzmann 7/8 and the cap ladder r = 1/(m+1)")
X = sp.symbols('x')
check("G1 at most m copies: (1-x)(1 + x + ... + x^m) = 1 - x^(m+1), m = 1..6 (SizeDialResidual.at_most_m)",
      all(sp.expand((1 - X)*sum(X**j for j in range(m + 1)) - (1 - X**(m + 1))) == 0 for m in range(1, 7)))
t0 = time.time()
def share(beta, m, energy=False):
    wgt = (lambda y: y**beta) if energy else (lambda y: y**(beta - 1))
    if energy:   # occupation n = 1/(e^y-1) - (m+1)/(e^{(m+1)y}-1); weight y^beta (energy density in beta dims)
        cap = mp.quad(lambda y: wgt(y)*(1/(mp.e**y - 1) - (m + 1)/(mp.e**((m + 1)*y) - 1)), [0, mp.inf])
        bose = mp.quad(lambda y: wgt(y)/(mp.e**y - 1), [0, mp.inf])
    else:        # log of the partition function: -log(1-e^-y) + log(1-e^-(m+1)y)
        cap = mp.quad(lambda y: wgt(y)*(-mp.log(1 - mp.e**-y) + mp.log(1 - mp.e**(-(m + 1)*y))), [0, mp.inf])
        bose = mp.quad(lambda y: wgt(y)*(-mp.log(1 - mp.e**-y)), [0, mp.inf])
    return cap/bose
err = max(abs(share(b, m) - (1 - mp.mpf(m + 1)**-b)) for b in (1, 2, 3) for m in (1, 2, 3))
tpol2 = time.time() - t0
check("G2 POLISH (30 digits): capped-at-m log Z / Bose log Z = 1 - (m+1)^(-beta) = 1 - r^beta, beta = 1, 2, 3, m = 1, 2, 3",
      err < mp.mpf(10)**-25, f"max error {mp.nstr(err, 3)}, {tpol2:.2f} s")
e1, e2 = share(3, 1, True), share(3, 2, True)
check("G3 energy density in 3-D: 7/8 at m = 1 (Stefan-Boltzmann's fermion/boson ratio), 26/27 at m = 2, 63/64 at m = 3",
      abs(e1 - mp.mpf(7)/8) < 1e-25 and abs(e2 - mp.mpf(26)/27) < 1e-25 and abs(share(3, 3, True) - mp.mpf(63)/64) < 1e-25)
ok = True
for m in (2, 3):
    for s in (2, 3, 4):
        em = sum((1 - m*(n % m == 0))*mp.mpf(n)**-s for n in range(1, 6001))  # partial sum over whole periods
        # tail correction: sum over a full period of (1 - m[m|n]) is 0, so the tail is small; compare to closed form
        ok &= abs(em - (1 - mp.mpf(m)**(1 - s))*mp.zeta(s)) < 1e-6
check("G4 Dirichlet form: sum (1 - m[m|n]) n^-s = (1 - m^(1-s)) zeta(s), so the leftover share is m^(-beta) at s = beta + 1 (m = 2 is eta)",
      ok)
check("G5 the ladder's stops: r = 1/2 (Fermi), 1/3, 1/4, ... -> 0 (Bose, share 0); share r^beta",
      all(abs((1 - mp.mpf(1)/(m + 1)**3) - (1 - (mp.mpf(1)/(m + 1))**3)) < 1e-30 for m in range(1, 20)))

# ---------------------------------------------------------------- C6: the clock
print("== C6  the clock r of the height flow is the d dial in disguise for rho")
rho = lambda k, d: 1/3 - 4*d/(3*k**2)
ok = all(abs(rho(k, k**2*(1 - lg(rr))/4) - lg(rr)/3) < 1e-14 for k in (0.5, 1, 2) for rr in (2, math.e, 6, math.exp(9/5)))
check("T1 rho = (ln r)/3 (newton-g Q1) equals the plane-law rho = 1/3 - 4d/(3k^2) at d = k^2 (1 - ln r)/4, for k = 1/2, 1, 2", ok)
check("T2 the parabola d = k^2/4 has rho = 0 for every k: the only clock-free value (newton-g Q4)",
      all(abs(rho(k, k**2/4)) < 1e-15 for k in (0.3, 1, 7)))

# ---------------------------------------------------------------- fractal strings
print("== Z  fractal strings and complex dimensions (Lapidus-van Frankenhuijsen)")
mp.mp.dps = 20
f = lambda s: 2*mp.mpf(1)/3**s - 1
D = mp.log(2)/mp.log(3)
check("Z1 Cantor: complex dimensions D + 2 pi i n / log 3 are zeros of 1 - 2 (1/3)^s (n = 1..5)",
      all(abs(f(D + 2j*mp.pi*n/mp.log(3))) < 1e-15 for n in range(1, 6)))
zc = lambda s: sum(2**n*mp.mpf(3)**(-(n + 1)*s) for n in range(400))
check("Z2 Cantor string lengths 3^-(n+1) (2^n of them): geometric zeta = 1/(3^s - 2), total length 1 at s = 1",
      abs(zc(1) - 1) < 1e-15 and abs(zc(mp.mpf('0.9')) - 1/(3**mp.mpf('0.9') - 2)) < 1e-12)
lens = [mp.mpf(3)**-(n + 1) for n in range(8) for _ in range(2**n)]
s = 2
direct = sum(sum((kk/l)**-s for kk in range(1, 3000)) for l in lens[:31])
fact = sum(1/mp.mpf(kk)**s for kk in range(1, 3000))*sum(l**s for l in lens[:31])
check("Z3 spectral zeta of a string = zeta(s) x geometric zeta (Lapidus): checked on the first 31 Cantor-string intervals at s = 2",
      abs(direct - fact) < 1e-12)
fu = lambda s: mp.mpf(0.3)**s + mp.mpf(0.5)**s - 1
Du = mp.findroot(fu, 0.8)
mins = min(abs(fu(Du + 1j*tt)) for tt in np.linspace(0.5, 120, 2400))
fl = lambda s: mp.mpf(0.25)**s + mp.mpf(0.5)**s - 1
Dl = mp.findroot(fl, 0.7)
check("Z4 one piece turned alone: lattice (1/4, 1/2) has complex dimensions on the line Re s = D, every 2 pi/log 2; non-lattice (0.3, 0.5) has none on that line up to height 120",
      abs(fl(Dl + 2j*mp.pi/mp.log(2))) < 1e-15 and mins > 1e-3, f"non-lattice min |x(s)-1| on the line = {mp.nstr(mins, 3)}")

# ---------------------------------------------------------------- N: in-between laws
print("== N  laws at in-between settings")
check("N1 two pieces at r = 2^(-1/D) have similarity dimension D: D = 1/2 (r = 1/4), 1 (1/2), 2 (0.7071), 3 (0.7937)",
      all(abs(moran([2**(-1/Dd)]*2) - Dd) < 1e-10 for Dd in (0.5, 1, 2, 3)))
lam = 2**-0.5
def nu_hat(xi, terms=80): return np.prod([np.cos(math.pi*xi*lam**j) for j in range(terms)])
sinc = lambda z: 1.0 if z == 0 else math.sin(math.pi*z)/(math.pi*z)
ok = all(abs(nu_hat(xi) - sinc(2*xi)*sinc(2*xi*lam)) < 1e-9 for xi in (0.1, 0.37, 1.3, 2.9))
check("N2 at r = 1/sqrt 2 the natural measure's Fourier transform is sinc(2 xi) sinc(sqrt2 xi): digits split into even and odd, so it is uniform * uniform (a trapezoid density), the square seen edge-on (Wintner 1935)",
      ok)
n = 6
seg = {sum(Fr(e, 2**(i + 1)) for i, e in enumerate(eps)) for eps in __import__('itertools').product((0, 1), repeat=2*n)}
cq = [sum(Fr(3*e, 4**(j + 1)) for j, e in enumerate(eps)) for eps in __import__('itertools').product((0, 1), repeat=n)]
mix = {Fr(2, 3)*a + Fr(1, 3)*b for a in cq for b in cq}
check("N3 the joined half boxes, geometrically: (2/3) C_1/4 + (1/3) C_1/4 = the dyadic segment, level by level (exact, level 6: 4096 points)",
      seg == mix and len(seg) == 4**n)
diag = {(a + b)/2 for a in cq for b in cq}
check("N4 but (1/2)(C_1/4 + C_1/4) has only 3^n points at level n (exact overlap): a dust of D = log3/log4 = 0.79, not 1",
      len(diag) == 3**n)
def sumset_maxgap(r, n=9):
    lefts = np.array([0.0])
    for _ in range(n):
        lefts = np.unique(np.round(np.concatenate([r*lefts, r*lefts + (1 - r), r*lefts + 2*(1 - r)]), 14))
    L = 2*r**n
    ends = np.sort(lefts); gap = np.max(np.maximum(ends[1:] - (ends[:-1] + L), 0)) if len(ends) > 1 else 0
    return gap
gp = {r: sumset_maxgap(r) for r in (0.3, 1/3, 0.4)}
check("N5 C_r + C_r is the 3-map set with shifts 0, 1-r, 2(1-r): gaps at r = 0.3, an interval from r = 1/3 on (cover condition 2r >= 1 - r)",
      gp[0.3] > 1e-3 and gp[1/3] < 1e-12 and gp[0.4] < 1e-12, ", ".join(f"r={r:.3f}: max gap {g:.2e}" for r, g in gp.items()))

def garsia(pw, n_max):
    dist = Counter({(0,)*len(pw[0]): 1}); H = []
    for i in range(n_max):
        li = pw[i]; new = Counter()
        for v, c in dist.items():
            new[v] += c; new[tuple(a + b for a, b in zip(v, li))] += c
        dist = new; tot = 2**(i + 1)
        H.append(-sum(c/tot*math.log(c/tot) for c in dist.values()))
        if i < 12: counts.append(len(dist))
    return H
counts = []
pw = [(1, 0), (0, 1)]
for i in range(2, 40): pw.append(tuple(a - b for a, b in zip(pw[i - 2], pw[i - 1])))
t0 = time.time(); Hg = garsia(pw, 24); tg = time.time() - t0
fib = [1, 1]
for _ in range(20): fib.append(fib[-1] + fib[-2])
check("N6 golden r = 1/phi: exact overlaps (1 = r + r^2); n digits give F(n+3) - 1 distinct points, not 2^n",
      counts[:10] == [fib[nn + 2] - 1 for nn in range(1, 11)], f"{counts[:8]}")
lamg = (math.sqrt(5) - 1)/2
dimg = (Hg[-1] - Hg[-2])/lg(1/lamg)
check("N7 golden: dimension of the natural measure = Garsia entropy / log phi = 0.99571 < 1 (Alexander-Zagier 1991)",
      abs(dimg - 0.99571) < 5e-5, f"{dimg:.6f}, {tg:.1f} s")
counts = []
pw3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
for i in range(3, 40): pw3.append(tuple(a - b - c for a, b, c in zip(pw3[i - 3], pw3[i - 2], pw3[i - 1])))
lam3 = [z.real for z in np.roots([1, 1, 1, -1]) if abs(z.imag) < 1e-12][0]
H3 = garsia(pw3, 21)
dim3 = (H3[-1] - H3[-2])/lg(1/lam3)
check("N8 tribonacci r = 0.5437 (also Pisot): dimension 0.98041 < 1 (Grabner-Kirschenhofer-Tichy's 0.9804)",
      abs(dim3 - 0.9804) < 2e-4, f"{dim3:.6f}")

print(f"\npolish cost: {tpol1 + tpol2:.2f} s for checks L3 and G2 (mpmath, 30 digits)")
print(f"{sum(results)} of {len(results)} pass")
