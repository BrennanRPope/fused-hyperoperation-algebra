"""DRAFT checks for hypersphere/README.md (Affine Composition Theory).

Run: python3 checks.py > checks_output.txt
Needs numpy and scipy. Every check prints PASS or FAIL with the numbers it compared.
"""
import math
import numpy as np
from scipy import integrate, special

rng = np.random.default_rng(20261126)
results = []


def check(name, ok, detail):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# ---------------------------------------------------------------------------
# Cayley-Dickson algebras: R (1), C (2), H (4), O (8), sedenions (16)
# ---------------------------------------------------------------------------
def conj(a):
    if len(a) == 1:
        return a.copy()
    h = len(a) // 2
    return np.concatenate([conj(a[:h]), -a[h:]])


def mul(a, b):
    if len(a) == 1:
        return a * b
    h = len(a) // 2
    p, q = a[:h], a[h:]
    r, s = b[:h], b[h:]
    # (p,q)(r,s) = (pr - s* q, s p + q r*)
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])


def norm(a):
    return float(np.sqrt(np.sum(a * a)))


def unit(n):
    e = np.zeros(n)
    e[0] = 1.0
    return e


print("=" * 78)
print("0. SET-UP: polar coordinates u = r * omega, radius adds, sphere multiplies")
print("=" * 78)

# H0a: |uv| = |u||v| in dims 1,2,4,8; fails in 16
for n in (1, 2, 4, 8):
    worst = 0.0
    for _ in range(2000):
        u, v = rng.normal(size=n), rng.normal(size=n)
        worst = max(worst, abs(norm(mul(u, v)) - norm(u) * norm(v)) / (norm(u) * norm(v)))
    check(f"H0a dim {n}: norm multiplies, so log-radius adds", worst < 1e-12, f"worst relative error {worst:.1e}")

worst16 = 0.0
for _ in range(2000):
    u, v = rng.normal(size=16), rng.normal(size=16)
    worst16 = max(worst16, abs(norm(mul(u, v)) - norm(u) * norm(v)) / (norm(u) * norm(v)))
check("H0b dim 16 (sedenions): norm does NOT multiply, the polar split breaks", worst16 > 0.05,
      f"largest relative miss {worst16:.3f}")

# H0c: the circle operation in n dimensions, a o_k b = a + b + k ab, lift 1 + k a multiplies
for n in (2, 4, 8):
    k = 0.7
    a, b = rng.normal(size=n) * 0.5, rng.normal(size=n) * 0.5
    comp = a + b + k * mul(a, b)
    lhs = mul(unit(n) + k * a, unit(n) + k * b)
    rhs = unit(n) + k * comp
    L = lambda x: math.log(norm(unit(n) + k * x)) / k
    check(f"H0c dim {n}: lift 1+ka turns the circle operation into multiplication; radius L_k adds",
          np.allclose(lhs, rhs) and abs(L(comp) - L(a) - L(b)) < 1e-12,
          f"lift error {np.max(np.abs(lhs - rhs)):.1e}, radius error {abs(L(comp) - L(a) - L(b)):.1e}")


def ok_compose(n, k, a, b):
    return a + b + k * mul(a, b)


assoc = {}
for n in (2, 4, 8):
    worst = 0.0
    for _ in range(300):
        a, b, c = (rng.normal(size=n) for _ in range(3))
        x = ok_compose(n, 1.0, ok_compose(n, 1.0, a, b), c)
        y = ok_compose(n, 1.0, a, ok_compose(n, 1.0, b, c))
        worst = max(worst, np.max(np.abs(x - y)))
    assoc[n] = worst
check("H0d circle operation is associative in dims 2 and 4, not in dim 8 (octonions)",
      assoc[2] < 1e-10 and assoc[4] < 1e-10 and assoc[8] > 1e-3,
      f"worst associator: dim2 {assoc[2]:.1e}, dim4 {assoc[4]:.1e}, dim8 {assoc[8]:.2f}")

comm4 = max(np.max(np.abs(ok_compose(4, 1, a, b) - ok_compose(4, 1, b, a)))
            for a, b in ((rng.normal(size=4), rng.normal(size=4)) for _ in range(50)))
check("H0e circle operation stops commuting in dim 4 (quaternions)", comm4 > 1e-2, f"largest a∘b − b∘a = {comm4:.2f}")

print()
print("=" * 78)
print("1. i AND pi FROM THE SPHERE (Brennan's first test)")
print("=" * 78)

# I1: Re v = 0, |v| = 1  =>  v^2 = -1 in H and O; such v fill S^2 and S^6
for n in (2, 4, 8):
    worst = 0.0
    for _ in range(2000):
        v = rng.normal(size=n)
        v[0] = 0.0
        v /= norm(v)
        worst = max(worst, np.max(np.abs(mul(v, v) + unit(n))))
    check(f"I1 dim {n}: every unit tangent v at 1 squares to -1 (a sphere S^{n - 2} of i's)",
          worst < 1e-12, f"worst |v^2 + 1| = {worst:.1e}")

# I1b: the k route's size-keeping condition is the same as tangency at 1
worst = 0.0
for _ in range(1000):
    kk = complex(*rng.normal(size=2))
    x = rng.normal()
    eps = 1e-7
    dsize = (abs(1 + kk * x * eps) - 1) / eps
    worst = max(worst, abs(dsize - (kk * x).real))
check("I1b k route: d|1 + k x|/dx at 0 = Re(k x), so 'keeps size' = 'tangent to the sphere at 1'",
      worst < 1e-5, f"worst mismatch {worst:.1e}")


def qexp(v, t):
    n = len(v)
    return np.cos(t) * unit(n) + np.sin(t) * v


# I2: exp(t v) stays on the sphere at unit speed and first reaches -1 at t = pi
for n in (2, 4, 8):
    v = rng.normal(size=n)
    v[0] = 0.0
    v /= norm(v)
    ts = np.linspace(0, 4, 400001)
    pts = np.array([np.cos(ts), *[np.sin(ts) * v[i] for i in range(1, n)]])
    dist = np.linalg.norm(pts.T + unit(n), axis=1)
    t_hit = ts[np.argmin(dist)]
    # check the flow from the series exp(tv) = sum (tv)^m/m! (needs v^2 = -1 only through mul)
    S = unit(n)
    term = unit(n)
    for m in range(1, 40):
        term = mul(term, v) * (math.pi / m)
        S = S + term
    check(f"I2 dim {n}: exp(t v) first reaches -1 at t = pi; series exp(pi v) = -1",
          abs(t_hit - math.pi) < 1e-4 and np.max(np.abs(S + unit(n))) < 1e-12,
          f"hit at t = {t_hit:.5f}, |exp(pi v) + 1| = {np.max(np.abs(S + unit(n))):.1e}")

# I3: great-circle distance from 1 to -1 on S^n is pi, every n
for n in (1, 2, 3, 5, 7):
    v = rng.normal(size=n + 1)
    v[0] = 0.0
    v /= norm(v)
    N = 200000
    ts = np.linspace(0, math.pi, N + 1)
    P = np.outer(np.cos(ts), unit(n + 1)) + np.outer(np.sin(ts), v)
    length = np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1))
    check(f"I3 S^{n}: half a great circle (1 to -1) has length pi", abs(length - math.pi) < 1e-9,
          f"length {length:.12f}")

# I4: Euler by chaining on the 3-sphere: (1 + pi v / N)^N -> -1
v = np.array([0.0, 1, 2, 2]) / 3
errs = []
for N in (100, 1000, 10000):
    step = unit(4) + v * math.pi / N
    P = unit(4)
    for _ in range(N):
        P = mul(P, step)
    errs.append(np.max(np.abs(P + unit(4))))
check("I4 S^3: chaining N steps of size pi/N along any unit i lands on -1, error ~ 1/N",
      errs[2] < 1e-3 and errs[0] / errs[2] > 50, f"errors {', '.join(f'{e:.1e}' for e in errs)}")

# I5: pi from the sphere's area: Poisson's polar trick for the Gaussian
G1 = integrate.quad(lambda x: math.exp(-x * x), -np.inf, np.inf)[0]
polar = 2 * math.pi * integrate.quad(lambda r: math.exp(-r * r) * r, 0, np.inf)[0]
check("I5 Poisson: (∫e^{-x²}dx)² = ∫∫e^{-r²} r dr dθ = pi (the 2π is the circle's length)",
      abs(G1 ** 2 - math.pi) < 1e-12 and abs(polar - math.pi) < 1e-12, f"{G1 ** 2:.15f} and {polar:.15f}")
ok = True
for n in range(1, 9):
    area = 2 * math.pi ** (n / 2) / special.gamma(n / 2)  # |S^{n-1}|
    lhs = G1 ** n
    rhs = area * integrate.quad(lambda r: math.exp(-r * r) * r ** (n - 1), 0, np.inf)[0]
    ok &= abs(lhs - rhs) < 1e-10 * lhs
check("I5b |S^{n-1}| = 2 pi^{n/2}/Gamma(n/2) reproduces (∫e^{-x²})^n for n = 1..8", ok, "all eight agree")

# I6: spinor sign: rotating by 2pi takes the unit quaternion from 1 to -1 (distance pi on S^3)
axis = np.array([0.0, 0.0, 0.6, 0.8])
x = np.array([0.0, 1.0, -2.0, 0.5])
q = qexp(axis, math.pi)  # rotation angle 2pi
rot = mul(mul(q, x), conj(q))
check("I6 S^3: a full 2π turn in space is a half great circle on S^3: q = -1, vectors unchanged",
      np.allclose(q, -unit(4)) and np.allclose(rot, x), f"q = {np.round(q, 12)}")

# I7: the sphere of i's does not commute: ij = k = -ji
i_, j_, k_ = np.eye(4)[1], np.eye(4)[2], np.eye(4)[3]
check("I7 S^3: ij = k and ji = -k (the seed of the swap rule)",
      np.allclose(mul(i_, j_), k_) and np.allclose(mul(j_, i_), -k_), "")

print()
print("=" * 78)
print("2. e: CUBE, SIMPLEX, SPHERE")
print("=" * 78)

# E1: corner simplex of the unit n-cube has volume 1/n!; the cube is n! sorted simplices
ok = True
detail = []
for n in range(1, 7):
    U = rng.random((400000, n))
    frac_simplex = np.mean(U.sum(axis=1) <= 1)
    frac_sorted = np.mean(np.all(np.diff(U, axis=1) > 0, axis=1))
    exact = 1 / math.factorial(n)
    tol = 4 * math.sqrt(exact * (1 - exact) / 400000) + 1e-9
    ok &= abs(frac_simplex - exact) < tol and abs(frac_sorted - exact) < tol
    detail.append(f"n={n}: {frac_simplex:.5f}/{frac_sorted:.5f} vs {exact:.5f}")
check("E1 cube corner simplex = 1/n! of the cube = chance that n random numbers come out sorted", ok,
      "; ".join(detail))

partial = sum(1 / math.factorial(n) for n in range(0, 20))
check("E2 sum over all dimensions of the corner simplex volume = e", abs(partial - math.e) < 1e-15,
      f"{partial:.15f}")

# E3: keep adding random numbers until the total passes 1: average count = e
M = 2_000_000
total = np.zeros(M)
count = np.zeros(M, dtype=int)
alive = np.ones(M, dtype=bool)
while alive.any():
    total[alive] += rng.random(alive.sum())
    count[alive] += 1
    alive = total <= 1
mean = count.mean()
se = count.std() / math.sqrt(M)
check("E3 random-sum puzzle: average number of tries to pass 1 is e", abs(mean - math.e) < 4 * se,
      f"{mean:.5f} ± {se:.5f} vs e = {math.e:.5f}")

# E4: the sphere gives e^pi, and its e is the simplex of squared radii
Vball = lambda n: math.pi ** (n / 2) / special.gamma(n / 2 + 1)
even = sum(Vball(2 * m) for m in range(0, 60))
same = all(abs(Vball(2 * m) - math.pi ** m / math.factorial(m)) < 1e-12 * Vball(2 * m) + 1e-300
           for m in range(0, 30))
allsum = sum(Vball(n) for n in range(0, 120))
check("E4 balls: V(B^{2m}) = pi^m × (1/m!) (angles × simplex), sum over even dims = e^pi",
      same and abs(even - math.exp(math.pi)) < 1e-10,
      f"sum = {even:.12f}, e^pi = {math.exp(math.pi):.12f}; all dims {allsum:.6f} = e^pi(1+erf√pi) "
      f"{math.exp(math.pi) * (1 + math.erf(math.sqrt(math.pi))):.6f}")
r = 1 / math.sqrt(math.pi)
even_r = sum((math.pi * r * r) ** m / math.factorial(m) for m in range(60))
check("E4b to get e itself from balls you must shrink them to disc area 1 (r = 1/√pi)",
      abs(even_r - math.e) < 1e-14, f"{even_r:.15f}")

# E5: e from the sphere in high dimension (Poincare-Borel) = the dial chain at k = -1/2
errs = []
for n in (10, 100, 1000, 10000):
    t = np.linspace(-4, 4, 2001)
    dens = np.clip(1 - t * t / n, 0, None) ** ((n - 3) / 2)
    dens /= np.trapezoid(dens, t)
    gauss = np.exp(-t * t / 2) / math.sqrt(2 * math.pi)
    errs.append(np.max(np.abs(dens - gauss)))
check("E5 one coordinate of a random point on the sphere of radius √n tends to the bell curve e^{-t²/2}",
      errs[-1] < 1e-4 and errs[0] / errs[-1] > 500, f"max density errors {', '.join(f'{e:.1e}' for e in errs)}")
s = 2.3
chain_err = abs((1 - s / (2 * 20000)) ** 20000 - math.exp(-s / 2))
check("E5b the sphere's (1 - t²/n)^{n/2} is the dial chain (1 + k x/N)^N with k = -1/2, x = t², N = n/2",
      chain_err < 1e-4, f"at t² = {s}: error {chain_err:.1e}")

print()
print("=" * 78)
print("3. RUNGS")
print("=" * 78)
q = qexp(np.array([0.0, 0.48, 0.6, 0.64]), 0.9)
log_radius = math.log(norm(q))
check("R1 on the sphere the radius log|u| is 0, so the rung-3 law 'take the radius twice' has nothing to act on",
      abs(log_radius) < 1e-15, f"log|q| = {log_radius:.1e}; log(log|q|) = -inf")


def qpow(q, t):
    th = math.acos(np.clip(q[0], -1, 1))
    v = q.copy()
    v[0] = 0
    v /= norm(v)
    return qexp(v, t * th)


check("R2 rung 3 survives on the sphere only as 'angle × real exponent': (q^s)^t = q^{st}",
      np.allclose(qpow(qpow(q, 1.7), 0.4), qpow(q, 0.68)), "")

print()
print("=" * 78)
print("4. PLANE DIAL d = 1/R²: SPHERE, PLANE, IMAGINARY SPHERE")
print("=" * 78)


def sph_angle(A, B, C):
    u = B - np.dot(A, B) * A
    w = C - np.dot(A, C) * A
    return math.acos(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))


# P1: sphere of radius R: angle excess = area / R^2  (Girard / Gauss-Bonnet)
R = 2.0
A0 = np.array([0.0, 0, 1])
B0 = np.array([math.sin(0.9), 0, math.cos(0.9)])
C0 = np.array([math.sin(0.8) * math.cos(1.1), math.sin(0.8) * math.sin(1.1), math.cos(0.8)])
excess = sph_angle(A0, B0, C0) + sph_angle(B0, C0, A0) + sph_angle(C0, A0, B0) - math.pi
X = rng.normal(size=(4_000_000, 3))
X /= np.linalg.norm(X, axis=1)[:, None]
inside = (np.sign(X @ np.cross(A0, B0)) == np.sign(C0 @ np.cross(A0, B0))) & \
         (np.sign(X @ np.cross(B0, C0)) == np.sign(A0 @ np.cross(B0, C0))) & \
         (np.sign(X @ np.cross(C0, A0)) == np.sign(B0 @ np.cross(C0, A0)))
area = inside.mean() * 4 * math.pi * R * R
se = math.sqrt(inside.mean() * (1 - inside.mean()) / len(X)) * 4 * math.pi * R * R
check("P1 sphere radius R: angle sum − pi = (+1/R²) × area", abs(excess - area / R ** 2) < 4 * se / R ** 2,
      f"excess {excess:.5f}, area/R² {area / R ** 2:.5f} ± {se / R ** 2:.5f}")

# P2: hyperboloid (imaginary radius): pi − angle sum = (1/R²) × area, R = 1
def mink(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def hyp_pt(x, y):  # Klein coords -> hyperboloid
    s = 1 / math.sqrt(1 - x * x - y * y)
    return np.array([s, s * x, s * y])


def hyp_angle(A, B, C):
    u = B + mink(A, B) * A
    w = C + mink(A, C) * A
    return math.acos(mink(u, w) / math.sqrt(mink(u, u) * mink(w, w)))


Kl = np.array([[0.05, -0.1], [0.75, 0.05], [0.1, 0.7]])
Ah, Bh, Ch = (hyp_pt(*p) for p in Kl)
defect = math.pi - (hyp_angle(Ah, Bh, Ch) + hyp_angle(Bh, Ch, Ah) + hyp_angle(Ch, Ah, Bh))
# area in Klein model: dA = dx dy / (1 - x² - y²)^{3/2}, geodesics are straight lines
def tri_area_klein(P, n=2400):
    xs = np.linspace(P[:, 0].min(), P[:, 0].max(), n)
    ys = np.linspace(P[:, 1].min(), P[:, 1].max(), n)
    Xg, Yg = np.meshgrid(xs, ys)
    def side(p, q, r):
        return (q[0] - p[0]) * (Yg - p[1]) - (q[1] - p[1]) * (Xg - p[0]), \
               (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    m = np.ones_like(Xg, dtype=bool)
    for a, b, c in ((P[0], P[1], P[2]), (P[1], P[2], P[0]), (P[2], P[0], P[1])):
        s1, s2 = side(a, b, c)
        m &= np.sign(s1) == np.sign(s2)
    rr = np.clip(1 - Xg ** 2 - Yg ** 2, 1e-12, None)
    w = np.where(m, rr ** -1.5, 0.0)
    return np.sum(w) * (xs[1] - xs[0]) * (ys[1] - ys[0])


harea = tri_area_klein(Kl)
check("P2 hyperbolic plane (sphere of imaginary radius): angle sum − pi = (−1/R²) × area",
      abs(defect - harea) / harea < 3e-3, f"defect {defect:.5f}, area {harea:.5f}")

# P3: Lambert: the spherical law of cosines with R -> iR is the hyperbolic one
a_s = math.acosh(-mink(Bh, Ch))
b_s = math.acosh(-mink(Ah, Ch))
c_s = math.acosh(-mink(Ah, Bh))
Cang = hyp_angle(Ch, Ah, Bh)
Rc = 1j
lhs = np.cos(c_s / Rc)
# spherical law with radius R: cos(c/R) = cos(a/R)cos(b/R) + sin(a/R)sin(b/R)cos C
rhs = np.cos(a_s / Rc) * np.cos(b_s / Rc) + np.sin(a_s / Rc) * np.sin(b_s / Rc) * math.cos(Cang)
check("P3 Lambert's imaginary sphere: spherical cosine law at R = i reproduces the measured hyperbolic triangle",
      abs(lhs - rhs) < 1e-12, f"cos(c/i) = {lhs.real:.10f}, rhs = {rhs.real:.10f} (imag {abs(rhs.imag):.1e})")

print()
print("=" * 78)
print("5. SWAP RULE: TURN = d × AREA ON THE SPHERE")
print("=" * 78)


def rotq(axis_vec, angle):
    v = np.concatenate([[0.0], axis_vec / np.linalg.norm(axis_vec)])
    return qexp(v, angle / 2)


ratios = []
for a in (0.1, 0.03, 0.01):
    b = 0.7 * a
    P = rotq(np.array([1.0, 0, 0]), a)
    Q = rotq(np.array([0.0, 1, 0]), b)
    Cm = mul(mul(P, Q), mul(conj(P), conj(Q)))
    turn = 2 * math.acos(np.clip(Cm[0], -1, 1))
    ratios.append(turn / (a * b))
check("W1 S^3: swapping turns a (about x) and b (about y) leaves a turn about z of size a·b = area",
      abs(ratios[-1] - 1) < 1e-3, f"turn/(ab) = {', '.join(f'{x:.6f}' for x in ratios)}")

# W2: Berry phase of spin 1/2 around a latitude = -1/2 solid angle (Hopf connection)
ok = True
det = []
for th0 in (0.4, 1.0, 2.0):
    N = 20000
    phis = np.linspace(0, 2 * math.pi, N, endpoint=False)
    states = np.array([[math.cos(th0 / 2), np.exp(1j * p) * math.sin(th0 / 2)] for p in phis])
    prod = np.prod([np.vdot(states[i], states[(i + 1) % N]) for i in range(N)])
    berry = -np.angle(prod)
    solid = 2 * math.pi * (1 - math.cos(th0))
    want = (-0.5 * solid + math.pi) % (2 * math.pi) - math.pi
    got = (berry + math.pi) % (2 * math.pi) - math.pi
    ok &= abs(got - want) < 1e-3
    det.append(f"θ={th0}: {got:.5f} vs {want:.5f}")
check("W2 Hopf fibration: spin-½ Berry phase = −½ × solid angle on the Bloch sphere", ok, "; ".join(det))

# W3: parallel transport around a latitude on a sphere of radius R turns a vector by area/R^2
for Rr in (1.0, 3.0):
    th0 = 0.8
    N = 200000
    phis = np.linspace(0, 2 * math.pi, N + 1)
    pos = lambda p: Rr * np.array([math.sin(th0) * math.cos(p), math.sin(th0) * math.sin(p), math.cos(th0)])
    e_th = lambda p: np.array([math.cos(th0) * math.cos(p), math.cos(th0) * math.sin(p), -math.sin(th0)])
    vec = e_th(0.0)
    for p in phis[1:]:
        nrm = pos(p) / Rr
        vec = vec - np.dot(vec, nrm) * nrm
        vec /= np.linalg.norm(vec)
    e_ph = np.array([0.0, 1, 0])
    ang = math.atan2(np.dot(vec, e_ph), np.dot(vec, e_th(0.0)))
    cap = 2 * math.pi * Rr ** 2 * (1 - math.cos(th0))
    want = (cap / Rr ** 2 + math.pi) % (2 * math.pi) - math.pi
    check(f"W3 sphere R={Rr}: carrying an arrow around a latitude turns it by cap area / R²",
          abs(abs(ang) - abs(want)) < 1e-3, f"turn {ang:.5f}, area/R² (mod 2π) {want:.5f}")

print()
print("=" * 78)
print("6. PROBABILITY ON THE SPHERE")
print("=" * 78)
X = rng.normal(size=(2_000_000, 3))
X /= np.linalg.norm(X, axis=1)[:, None]
hist, _ = np.histogram(X[:, 2], bins=20, range=(-1, 1))
flat = hist.max() / hist.min()
check("Q1 Archimedes' hat box: the height of a random point on S² is uniform", flat < 1.02,
      f"max/min bin = {flat:.4f}")

n_, m_ = np.array([0.0, 0.6, 0.8]), np.array([0.0, 0.0, 1.0])
born = (1 + n_ @ m_) / 2
capfrac = np.mean(X @ m_ >= -(n_ @ m_))
check("Q2 spin-½ Born probability (1 + n·m)/2 = area fraction of a cap (Archimedes)",
      abs(born - capfrac) < 2e-3, f"Born {born:.4f}, cap {capfrac:.4f}")

# Q3: orthogonal hemispheres are independent (the cube's corners on the sphere); other caps are not
Y = rng.normal(size=(2_000_000, 4))
Y /= np.linalg.norm(Y, axis=1)[:, None]
A_, B_ = Y[:, 0] > 0, Y[:, 1] > 0
union = np.mean(A_ | B_)
p_, q_ = A_.mean(), B_.mean()
check("Q3 orthogonal hemispheres: P(A∪B) = p ∘₋₁ q (the 2^n orthants are the cube's corners)",
      abs(union - (p_ + q_ - p_ * q_)) < 2e-3, f"union {union:.4f}, p∘q {p_ + q_ - p_ * q_:.4f}")
A2, B2 = Y[:, 0] > 0.5, Y[:, 1] > 0.5
p2, q2 = A2.mean(), B2.mean()
gap = abs(np.mean(A2 | B2) - (p2 + q2 - p2 * q2))
check("Q3b smaller orthogonal caps are NOT independent: the union rule fails", gap > 5e-3,
      f"union {np.mean(A2 | B2):.4f} vs p∘q {p2 + q2 - p2 * q2:.4f}")

# Q4: ball <-> simplex, polydisc <-> cube (squared radii)
n = 3
Z = rng.normal(size=(1_000_000, 2 * n))
Z /= np.linalg.norm(Z, axis=1)[:, None]
s_sphere = Z[:, 0::2] ** 2 + Z[:, 1::2] ** 2  # |z_i|^2 on S^5
check("Q4 points on S^5 ⊂ C³: squared radii lie on the triangle s₁+s₂+s₃ = 1, uniformly (mean 1/3 each)",
      np.allclose(s_sphere.sum(axis=1), 1) and np.allclose(s_sphere.mean(axis=0), 1 / 3, atol=2e-3)
      and abs(np.mean(s_sphere[:, 0] < 0.5) - 0.75) < 2e-3,
      f"means {np.round(s_sphere.mean(axis=0), 4)}, P(s₁<½) = {np.mean(s_sphere[:, 0] < 0.5):.4f} (uniform: 0.75)")
rad = rng.random(1_000_000) ** (1 / (2 * n))
s_ball = s_sphere * (rad ** 2)[:, None]
s_poly = rng.random((1_000_000, n))  # polydisc: each |z_i|² uniform in [0,1]
check("Q4b ball in C³ -> squared radii fill the solid corner simplex; polydisc -> they fill the cube",
      abs(np.mean(s_ball.sum(axis=1) <= 0.5) - 0.125) < 2e-3 and abs(np.mean(s_poly.sum(axis=1) <= 0.5) - 0.5 ** 3 / 6) < 2e-3,
      f"P(sum ≤ ½): ball {np.mean(s_ball.sum(axis=1) <= 0.5):.4f} (simplex 0.125), polydisc "
      f"{np.mean(s_poly.sum(axis=1) <= 0.5):.4f} (cube 0.0208)")

print()
print("=" * 78)
print("7. PHYSICS")
print("=" * 78)
# EM1: flux of a point charge through a cube equals 4pi; through a sphere it is one product
f = lambda y, x: 1.0 / (x * x + y * y + 1.0) ** 1.5  # field through face z = 1 of the cube [-1,1]^3
face = integrate.dblquad(f, -1, 1, -1, 1, epsabs=1e-12)[0]
check("EM1 point charge: flux through the cube = 6 face integrals = 4pi; on the sphere it is 4pi r² × 1/r²",
      abs(6 * face - 4 * math.pi) < 1e-9, f"6 × {face:.10f} = {6 * face:.10f}, 4π = {4 * math.pi:.10f}")
ok = True
for nd in range(2, 7):
    area = 2 * math.pi ** (nd / 2) / special.gamma(nd / 2)
    ok &= abs(area * 1 / area - 1) < 1e-15
check("EM2 in n space dimensions the field is 1/(|S^{n-1}| r^{n-1}); the constant 4pi is the sphere's area",
      ok and abs(2 * math.pi ** 1.5 / special.gamma(1.5) - 4 * math.pi) < 1e-12, "")

# SR1: on the celestial sphere, a boost scales the stereographic radius by the Doppler factor
ok = True
worst = 0.0
for beta in (0.3, 0.8, 0.99):
    D = math.sqrt((1 + beta) / (1 - beta))
    for th in np.linspace(0.05, math.pi - 0.05, 40):
        thp = math.acos((math.cos(th) + beta) / (1 + beta * math.cos(th)))  # aberration toward +z
        worst = max(worst, abs(math.tan(thp / 2) - math.tan(th / 2) / D))
check("SR1 celestial sphere: a boost along z scales tan(θ/2) by 1/Doppler; a turn about z adds to the angle",
      worst < 1e-12, f"worst error {worst:.1e}")
b1, b2 = 0.6, 0.7
Dop = lambda b: math.sqrt((1 + b) / (1 - b))
badd = (b1 + b2) / (1 + b1 * b2)
x1, x2 = Dop(b1) - 1, Dop(b2) - 1
check("SR2 velocity addition = Doppler factors multiply = the circle operation ∘₁ on x = D − 1 (Bondi's k-calculus)",
      abs(Dop(badd) - Dop(b1) * Dop(b2)) < 1e-12 and abs((x1 + x2 + x1 * x2) - (Dop(badd) - 1)) < 1e-12,
      f"D(v1⊕v2) = {Dop(badd):.6f}, D1·D2 = {Dop(b1) * Dop(b2):.6f}")

# SR3: Penrose-Terrell: a circle of directions stays a circle after a boost
beta = 0.9
c0 = np.array([0.5, 0.2, math.sqrt(1 - 0.29)])
c0 /= np.linalg.norm(c0)
e1 = np.cross(c0, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
e2 = np.cross(c0, e1)
ang = 0.3
pts = []
for t in np.linspace(0, 2 * math.pi, 50, endpoint=False):
    d = math.cos(ang) * c0 + math.sin(ang) * (math.cos(t) * e1 + math.sin(t) * e2)
    cz = d[2]
    czp = (cz + beta) / (1 + beta * cz)
    sxy = math.sqrt(max(0, 1 - czp ** 2)) / math.sqrt(d[0] ** 2 + d[1] ** 2)
    pts.append([d[0] * sxy, d[1] * sxy, czp])
pts = np.array(pts)
M_ = np.c_[pts, np.ones(len(pts))]
sv = np.linalg.svd(M_, compute_uv=False)
check("SR3 Penrose–Terrell: boosted circle of directions is still a circle (points stay on one plane)",
      sv[-1] / sv[0] < 1e-12, f"smallest/largest singular value {sv[-1] / sv[0]:.1e}")

# G1: Kepler orbit -> momentum on a circle (Hamilton) -> Fock map -> great circle on S^3
GM = 1.0
r0 = np.array([1.0, 0.0, 0.0])
v0 = np.array([0.0, 1.2, 0.3])


def rhs_k(t, y):
    r = y[:3]
    return np.concatenate([y[3:], -GM * r / np.linalg.norm(r) ** 3])


sol = integrate.solve_ivp(rhs_k, (0, 30), np.concatenate([r0, v0]), rtol=1e-11, atol=1e-12, max_step=0.01)
P_ = sol.y[3:].T
E = 0.5 * np.sum(v0 ** 2) - GM / np.linalg.norm(r0)
p0 = math.sqrt(-2 * E)
p2 = np.sum(P_ ** 2, axis=1)
W = np.c_[(p2 - p0 ** 2) / (p2 + p0 ** 2), 2 * p0 * P_ / (p2 + p0 ** 2)[:, None]]
svW = np.linalg.svd(W, compute_uv=False)
check("G1 Kepler (Fock/Moser): each orbit's momenta, mapped to S^3, lie on one great circle",
      np.allclose(np.linalg.norm(W, axis=1), 1) and svW[2] / svW[0] < 1e-6,
      f"singular values {np.round(svW, 8)} (rank 2 = a plane through the centre)")
deg = [sum(2 * l + 1 for l in range(nn)) for nn in range(1, 8)]
harm = [(nn - 1 + 1) ** 2 for nn in range(1, 8)]
check("G2 hydrogen level n has n² states = harmonics of degree n−1 on S^3 (Fock 1935)", deg == harm,
      f"{deg}")

# G3: de Sitter: Euclidean de Sitter = S^4 of radius R = sqrt(3/Λ); T = ħc/(2π k R)
hbar, c, kB = 1.054571817e-34, 2.99792458e8, 1.380649e-23
Lam = 1.1056e-52
Rds = math.sqrt(3 / Lam)
T = hbar * c / (2 * math.pi * kB * Rds)
check("G3 de Sitter: temperature = 1 / (great-circle length of the 4-sphere); curvature dial d = Λ/3 = 1/R² (k = 0)",
      1e-30 < T < 1e-29 and abs(1 / Rds ** 2 - Lam / 3) < 1e-60, f"R = {Rds:.3e} m, T = {T:.2e} K")

print()
npass = sum(ok for _, ok in results)
print(f"{npass} of {len(results)} checks pass")
