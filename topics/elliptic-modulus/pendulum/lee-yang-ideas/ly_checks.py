"""Checking the three Lee-Yang ideas from the pendulum (Nov 27). Reads ../../../lee-yang/ (dialgas.py, fisher_zeros.pkl)
without changing it. Run from this folder: python3 ly_checks.py -> ly_checks_output.txt, ly_checks.json."""
import sys, math, cmath, pickle, json, itertools
import numpy as np, mpmath as mp
sys.path.insert(0, '../../../lee-yang')
from dialgas import cycle_counts, brute_counts, torus_edges
LY = '../../../lee-yang/'
res = []
def rec(name, ok, val):
    res.append((name, bool(ok), str(val))); print(('PASS' if ok else 'FAIL'), name, val)

def ring_poly_zeros(n, k, dps=60):
    """exact ring partition function sum_m,e C[m][e] z^m k^e, roots in z, returned as zeta = k z."""
    C = cycle_counts(n); mp.mp.dps = dps
    coeffs = [sum(int(C[mm][e]) * mp.mpf(k)**e for e in range(C.shape[1])) for mm in range(n + 1)]
    roots = mp.polyroots(coeffs[::-1], maxsteps=400, extraprec=400)
    return [complex(r) * k for r in roots]

# ---- idea 3 first (it carries idea 2): the ring's zeros as a pendulum
# L1 closed form <=> pendulum form: theta = pi - arg(zeta) satisfies sin(theta/2) = sqrt(m) cos(psi_j/2), m = 1 - 1/k,
#     psi_j = pi(2j+1)/n equally spaced.
worst = 0
for k, n in [(3.0, 21), (1.5, 21), (8.0, 25)]:
    m = 1 - 1/k
    zs = ring_poly_zeros(n, k)
    pred = sorted(2*math.asin(math.sqrt(m)*abs(math.cos(math.pi*(2*j+1)/(2*n)))) for j in range(n))
    got = sorted(math.pi - abs(cmath.phase(z)) for z in zs)
    # each predicted angle appears for zeta and its conjugate
    pred2 = sorted(pred + pred)[::1]
    pred2 = sorted(set(round(p, 12) for p in pred))
    gotu = sorted(set(round(g, 9) for g in got))
    worst = max(worst, max(min(abs(g - p) for p in pred2) for g in got))
rec('L1 ring zeros: sin(theta/2) = sqrt(m) cos(psi_j/2), theta = pi - arg zeta, m = 1 - 1/k', worst < 1e-8, worst)

# L2 the Yang-Lee edge is a pendulum amplitude: sin^2(theta0/2) = m = 1 - 1/k, so cos(theta0) = 2/k - 1 (= the README's arccos(1-2/k) edge in phi)
worst = 0
for k in [1.5, 3.0, 6.0, 20.0]:
    th0 = 2*math.asin(math.sqrt(1 - 1/k)); phi_e = math.acos(1 - 2/k)
    worst = max(worst, abs((math.pi - th0) - phi_e))
rec('L2 edge: pi - theta0 = arccos(1 - 2/k) (the Lee-Yang page\'s gap edge)', worst < 1e-12, worst)

# L3 density: zeros per unit angle = cos(theta/2) x pendulum dwell density 1/(2 sqrt(m - sin^2(theta/2))), normalised
def zero_angles(k, n):
    m = 1 - 1/k
    return np.array([2*math.asin(math.sqrt(m)*abs(math.cos(math.pi*(2*j+1)/(2*n)))) for j in range(n)])
ks_worst = 0
for k in [1.5, 3.0, 10.0, 100.0]:
    m = 1 - 1/k; th0 = 2*math.asin(math.sqrt(m)); a = np.sort(zero_angles(k, 20001))
    # cdf from density cos(t/2)/(2 sqrt(m - sin^2(t/2))): substitute y = sin(t/2)/sqrt(m): cdf = (2/pi) asin(y)
    xs = np.linspace(0, th0*0.999999, 400)
    cdf_pred = (2/math.pi)*np.arcsin(np.sin(xs/2)/math.sqrt(m))
    cdf_emp = np.searchsorted(a, xs)/len(a)
    ks_worst = max(ks_worst, np.max(np.abs(cdf_pred - cdf_emp)))
rec('L3 zero density = cos(theta/2) x pendulum dwell density (KS distance, n = 20001)', ks_worst < 2e-4, ks_worst)

# L3b near the edge both pile up like (theta0 - theta)^(-1/2)
k = 3.0; m = 1 - 1/k; th0 = 2*math.asin(math.sqrt(m)); a = np.sort(zero_angles(k, 200001))
eps = np.array([1e-2, 3e-3, 1e-3]); cnt = np.array([np.sum(a > th0 - e) for e in eps])
slope = np.polyfit(np.log(eps), np.log(cnt), 1)[0]   # count within eps ~ eps^(1/2)
rec('L3b edge pile-up: zeros within eps of the edge ~ eps^(1/2)', abs(slope - 0.5) < 0.02, slope)

# L4 the separatrix: as k -> infinity (m -> 1) the arc closes and the density goes flat (cos x dwell = 1/2)
k = 1e8; m = 1 - 1/k; a = np.sort(zero_angles(k, 20001)); xs = np.linspace(0.05, 3.0, 200)
flat = np.max(np.abs(np.searchsorted(a, xs)/len(a) - xs/math.pi))
rec('L4 k -> infinity: zeros uniform on the whole circle (separatrix)', flat < 2e-3, flat)

# L5 the three energy bands of the ring: k > 1 swinging (on the circle, short of zeta = +1); k < 0 over the top (m > 1:
#     the circle all the way to zeta = +1, the rest on the positive real axis); 0 < k < 1 Wick image (m < 0: negative real axis)
ok = True; notes = []
for k in [3.0]:
    zs = ring_poly_zeros(41, k); ok &= all(abs(abs(z) - 1) < 1e-9 and z.real < 1 - 1e-6 for z in zs)
for k in [-1.0, -3.0]:
    zs = ring_poly_zeros(41, k); m = 1 - 1/k
    on = [z for z in zs if abs(abs(z) - 1) < 1e-9]; real = [z for z in zs if abs(z.imag) < 1e-9 and z.real > 0 and abs(abs(z)-1) > 1e-9]
    # predicted: off the circle exactly where sqrt(m) cos(psi/2) > 1
    npred = sum(1 for j in range(41) if math.sqrt(m)*abs(math.cos(math.pi*(2*j+1)/82)) > 1)
    ok &= len(on) + len(real) == len(zs) and len(real) == npred; notes.append((k, len(real), npred))
for k in [0.5]:
    zs = ring_poly_zeros(41, k); ok &= all(abs(z.imag) < 1e-9 and z.real < 0 for z in zs)
rec('L5 bands: k>1 swing, k<0 over the top (spill to +real exactly where sqrt(m)cos > 1), 0<k<1 negative real', ok, notes)

# ---- idea 2: the unit circle. The Lee-Yang circle |zeta| = 1 is the pendulum's circle of real angles (theta = pi - arg zeta):
#      the wall zeta = -1 is the pendulum at rest at the bottom, the zero-field point zeta = +1 is the top. The halfway line
#      |t| = 1 has no role: test it - the zeros' slope readings t = sin(theta)/theta' are not of size 1.
k = 3.0; m = 1 - 1/k; th = zero_angles(k, 41)
tval = np.sin(th)/(2*np.sqrt(np.maximum(m - np.sin(th/2)**2, 1e-300)))
rec('I2 halfway-line idea: zeros do not sit at |t| = 1 (spread of log10|t|)', np.std(np.log10(tval)) > 0.3, float(np.std(np.log10(tval))))

# ---- idea 1: the 2D dial gas is the Ising model, and its transition is the separatrix of Onsager's modulus
# D1 lattice-gas identity on a 3x3 torus: z^n k^e = (zk)^(N/2) (z^(1/2) k)^M k^(sum ss / 4)  with k = e^(4J), h = ln(sqrt(z) k)
L = 3; N = 9; ed = torus_edges(L, L)
worst = 0
for z, k in [(0.3, 2.0), (1.7, 0.6), (0.11, 5.0)]:
    Zd = sum(z**bin(S).count('1') * k**sum(1 for a, b in ed if (S >> a) & 1 and (S >> b) & 1) for S in range(1 << N))
    J = math.log(k)/4; h = math.log(math.sqrt(z)*k)
    Zi = 0
    for S in range(1 << N):
        s = [1 if (S >> i) & 1 else -1 for i in range(N)]
        Zi += math.exp(J*sum(s[a]*s[b] for a, b in ed) + h*sum(s))
    worst = max(worst, abs(Zd/((z*k)**(N/2)*Zi) - 1))
rec('D1 dial gas = Ising: k = e^(4J), zero field at z = 1/k^2 (3x3 torus, exact sum)', worst < 1e-12, worst)

# D2 Onsager's modulus in dial terms: k_O = 1/sinh^2(2J) = 4k/(k-1)^2; k_O = 1 exactly at k = 3 +- 2 sqrt 2
f = lambda k: 4*k/(k-1)**2
rec('D2 k_O = 4k/(k-1)^2 equals 1 at k = 3 +- 2sqrt2', max(abs(f(3+2*math.sqrt(2))-1), abs(f(3-2*math.sqrt(2))-1)) < 1e-12, (f(3+2*math.sqrt(2)), f(3-2*math.sqrt(2))))

# D3 Kramers-Wannier in dial terms: k_O -> 1/k_O; ordered side (k > k_c) has k_O < 1 (swinging), disordered (1 < k < k_c) k_O > 1 (over the top)
kc = 3 + 2*math.sqrt(2)
def dual(k):   # solve 4k*/(k*-1)^2 = (k-1)^2/(4k) on k* > 1
    q = (k-1)**2/(4*k); a = 1; b = -(2 + 4/q); c = 1   # q (k*-1)^2 = 4k*  -> k*^2 - (2+4/q)k* + 1 = 0
    r = (-b + math.sqrt(b*b - 4*a*c))/2; return r
pairs = [(k, dual(k)) for k in [6.0, 10.0, 30.0]]
ok = all((k > kc) == (d < kc) and abs(f(k)*f(d) - 1) < 1e-12 for k, d in pairs) and abs(dual(kc) - kc) < 1e-9
rec('D3 duality k <-> k* swaps k_O and 1/k_O (swinging <-> over the top), fixed at k_c', ok, [(round(k, 3), round(d, 4)) for k, d in pairs])

# D4 Fisher zeros of the finite tori against the self-dual curve |k_O| = 1, i.e. |k - 1|^2 = 4|k|
Fz = pickle.load(open(LY + 'fisher_zeros.pkl', 'rb'))
dist = lambda z: abs(abs(f(z)) - 1)
rng = np.random.default_rng(1)
out = {}
for Lf in (4, 5, 6):
    zs = np.array(Fz[Lf]); d = np.median([dist(z) for z in zs])
    # null: same moduli, random phases
    null = np.median([dist(abs(z)*cmath.exp(1j*rng.uniform(0, 2*math.pi))) for z in zs for _ in range(50)])
    out[Lf] = (round(float(d), 3), round(float(null), 3))
trend = out[4][0] > out[6][0]
rec('D4 Fisher zeros approach |k_O| = 1 (median | |k_O|-1 |, data vs random-phase null)', trend and all(v[0] < v[1] for v in out.values()), out)

npass = sum(r[1] for r in res); print(f'{npass} of {len(res)} pass')
json.dump(res, open('ly_checks.json', 'w'), indent=1)
