"""New-law pass: sweep the laws at the in-between settings of the roundness dial (smooth <-> angular, the p dial),
the shape dial (round <-> triangle), and their combination, and test every candidate law that appears there.
Each numerical search ends with a cheap polish; its cost is recorded.  Writes new_laws.json."""
import json, time, math
from fractions import Fraction
from contextlib import contextmanager
import numpy as np
from scipy import special, optimize

R, POL, RESULTS = {}, {}, []
def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok))); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
@contextmanager
def timed(name):
    t0 = time.perf_counter(); yield; POL[name] = time.perf_counter() - t0
rng = np.random.default_rng(7)

# ---------------------------------------------------------------------------
# N1 (roundness dial, and its crossing with the shape dial): Bell's bound as a floor
#   With Alice's observables in a symmetric shape K and Bob's in its dual, the best CHSH score is
#   sup over unit vectors x, y of ||x + y|| + ||x - y|| in the norm whose ball is the dual of K.
#   Candidate law: this is never below 2 sqrt 2, for EVERY symmetric shape, and equals 2^(1 + 1/min(p, q)) on the p dial
#   in every dimension.
# ---------------------------------------------------------------------------
def sum_constant_2d(gauge, grid=1500):
    t = np.linspace(0, 2*np.pi, grid, endpoint=False); D = np.column_stack([np.cos(t), np.sin(t)])
    U = D/gauge(D)[:, None]
    S = gauge((U[:, None, :] + U[None, :, :]).reshape(-1, 2)).reshape(grid, grid)
    Dd = gauge((U[:, None, :] - U[None, :, :]).reshape(-1, 2)).reshape(grid, grid)
    i, j = np.unravel_index(np.argmax(S + Dd), S.shape)
    return float((S + Dd)[i, j]), t[i], t[j]
def polish_2d(gauge, a0, b0):
    def f(v):
        x = np.array([np.cos(v[0]), np.sin(v[0])]); y = np.array([np.cos(v[1]), np.sin(v[1])])
        x = x/gauge(x[None])[0]; y = y/gauge(y[None])[0]
        return -(gauge((x + y)[None])[0] + gauge((x - y)[None])[0])
    r = optimize.minimize(f, [a0, b0], method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
    return -r.fun
def polygon_gauge(V):                       # gauge of the symmetric polygon with vertices V (and -V), via its facets
    from scipy.spatial import ConvexHull
    P = np.vstack([V, -V]); h = ConvexHull(P)
    A = h.equations[:, :2]; b = -h.equations[:, 2]
    g = lambda Z: np.max(Z @ (A/b[:, None]).T, axis=-1)
    g.verts = P[h.vertices]
    return g
def exact_polygon(g):
    """For a polygon the sum ||x+y|| + ||x-y|| is convex in x and in y, so its max over the ball sits at vertex pairs."""
    Vt = g.verts; S = Vt[:, None, :] + Vt[None, :, :]; Dd = Vt[:, None, :] - Vt[None, :, :]
    return float(np.max(g(S.reshape(-1, 2)) + g(Dd.reshape(-1, 2))))
floor_vals = []
t0 = time.perf_counter()
for trial in range(300):
    m = rng.integers(2, 9)                  # 4 to 16 vertices
    ang = np.sort(rng.uniform(0, np.pi, m)); rad = rng.uniform(0.3, 1.0, m)
    V = np.column_stack([rad*np.cos(ang), rad*np.sin(ang)])
    try: g = polygon_gauge(V)
    except Exception: continue
    floor_vals.append(exact_polygon(g))
grid_time = time.perf_counter() - t0
# smooth random shapes too: random symmetric "superellipse-like" norms (sum of |<z, d_k>|^r)
for trial in range(100):
    k = rng.integers(2, 6); Dk = rng.normal(size=(k, 2)); r_ = rng.uniform(1.05, 6)
    g = lambda Z, Dk=Dk, r_=r_: np.sum(np.abs(Z @ Dk.T)**r_, axis=-1)**(1/r_)
    v, a0, b0 = sum_constant_2d(g, grid=500)
    t1 = time.perf_counter(); floor_vals.append(max(v, polish_2d(g, a0, b0))); POL["N1 polish (smooth shapes)"] = POL.get("N1 polish (smooth shapes)", 0) + time.perf_counter() - t1
floor_vals = np.array(floor_vals)
lowest = int(np.argmin(floor_vals))
R["bell_floor"] = dict(n_shapes=len(floor_vals), min=float(floor_vals.min()), max=float(floor_vals.max()),
                       hist=np.histogram(floor_vals, bins=24, range=(2.7, 4.05))[0].tolist())
check("N1a no symmetric shape (400 random polygons and smooth norms) scores below 2 sqrt 2",
      floor_vals.min() >= 2*np.sqrt(2) - 2e-3, f"lowest {floor_vals.min():.4f}, highest {floor_vals.max():.4f}; 2 sqrt 2 = {2*np.sqrt(2):.4f}")
# the floor is reached by more than the disk: regular 8-gon, and the p dial's formula in 3 dimensions
oct_g = polygon_gauge(np.column_stack([np.cos(np.pi*np.arange(4)/4), np.sin(np.pi*np.arange(4)/4)]))
with timed("N1 polish (octagon)"):
    v, a0, b0 = sum_constant_2d(oct_g); v_oct = polish_2d(oct_g, a0, b0)
check("N1b the regular octagon also sits exactly on the floor (so the floor does not single out the round shape)",
      abs(v_oct - 2*np.sqrt(2)) < 1e-9, f"{v_oct:.12f}")
def lp3(p):
    def nrm(z): return np.sum(np.abs(z)**p)**(1/p)
    def f(w):
        x, y = w[:3], w[3:]; x = x/nrm(x); y = y/nrm(y); return -(nrm(x + y) + nrm(x - y))
    best = 0
    for _ in range(60):
        r = optimize.minimize(f, rng.normal(size=6), method="Nelder-Mead", options=dict(maxiter=3000, xatol=1e-10, fatol=1e-12))
        best = max(best, -r.fun)
    return best
with timed("N1 polish (3D p-balls)"):
    d3 = {p: lp3(p) for p in [1.25, 1.5, 3.0, 5.0]}
form = lambda p: 2**(1 + 1/min(p, p/(p - 1)))
R["bell_3d"] = {str(p): v for p, v in d3.items()}
check("N1c in 3 dimensions the p-ball still scores 2^(1 + 1/min(p, q)): one formula in every dimension",
      all(abs(d3[p] - form(p)) < 1e-6 for p in d3), str({p: round(v, 8) for p, v in d3.items()}))
# along the shape dial the score goes the other way: Janotta et al. (2011) odd polygons (triangle = classical trit = 2;
# pentagon 6/sqrt5 = 2.683), all below 2 sqrt 2.  We quote these; the saddle picture is our reading.
R["bell_saddle_shape_dial"] = {"3": 2.0, "5": 6/np.sqrt(5), "inf": 2*np.sqrt(2)}

# ---------------------------------------------------------------------------
# N2 (roundness dial at rational settings): Fermat for every rational p
#   Candidate law: |x|^p + |y|^p = 1 with p = m/n in lowest terms has infinitely many rational points only for
#   m = 1 or 2 (on p >= 1: p = 1 and p = 2), plus p = infinity; for every other rational p only the 4 trivial points.
#   Proof: x rational and x^(m/n) rational with gcd(m, n) = 1 force x = u^n, so the curve becomes u^m + v^m = 1 (Fermat m).
# ---------------------------------------------------------------------------
def is_perfect_power(fr, n):
    if fr < 0: return None
    def iroot(k):
        r = round(k**(1/n))
        for c in (r - 1, r, r + 1):
            if c >= 0 and c**n == k: return c
        return None
    a, b = iroot(fr.numerator), iroot(fr.denominator)
    return None if a is None or b is None else Fraction(a, b)
def rational_points_rational_p(m, n, D=60):
    """All rational (x, y), x, y >= 0, denominators <= D**n after the forced substitution, on x^(m/n) + y^(m/n) = 1.
    Exhaustive over u = a/c with c <= D: x = u^n, then y must be v^n with u^m + v^m = 1 exactly."""
    pts = set()
    for c in range(1, D + 1):
        for a in range(0, c + 1):
            u = Fraction(a, c); rest = 1 - u**m
            v = is_perfect_power(rest, m)
            if v is not None: pts.add((u**n, v**n))
    return pts
# brute-force confirmation that no non-n-th-power x slips through: scan x = a/c directly and test y numerically
def brute(p, D=80):
    hits = 0
    for c in range(1, D + 1):
        for a in range(1, c):
            y = (1 - (a/c)**p)**(1/p); fr = Fraction(y).limit_denominator(D)
            if fr.denominator <= D and abs(float(fr) - y) < 1e-13: hits += 1
    return hits
cases = [(3, 2), (4, 3), (5, 2), (5, 3), (5, 4), (7, 4), (6, 5), (7, 2), (9, 4)]
with timed("N2 exact search"):
    counts = {f"{m}/{n}": len(rational_points_rational_p(m, n)) for m, n in cases}
    counts["1/1"] = len(rational_points_rational_p(1, 1)); counts["2/1"] = len(rational_points_rational_p(2, 1))
    br = {f"{m}/{n}": brute(m/n) for m, n in cases}
R["fermat_rational_p"] = counts; R["fermat_brute"] = br
check("N2a every rational p = m/n with m >= 3 tested has only the trivial points (here: x or y = 0)",
      all(counts[f"{m}/{n}"] == 2 for m, n in cases) and all(v == 0 for v in br.values()),
      f"exact {counts}; brute-force interior hits {br}")
check("N2b p = 1 and p = 2 keep filling up (here " + str(counts['1/1']) + " and " + str(counts['2/1']) + " points in one quadrant)",
      counts["1/1"] > 1000 and counts["2/1"] > 10)
# Lean coverage: our Lean proves exponents 3 and 4, hence every m divisible by 3 or 4
covered = [f"{m}/{n}" for m, n in cases if m % 3 == 0 or m % 4 == 0]
R["fermat_lean_covered"] = covered

# ---------------------------------------------------------------------------
# N3 (roundness dial): where Lee-Yang breaks
#   A spin uniform in the n-dimensional p-ball, read along one axis, has weight (1 - |x|^p)^a with a = (n - 1)/p.
#   Its field function is G(s) = int_0^1 (1 - x^p)^a cos(s x) dx; Lee-Yang = every zero of G is real.
#   Candidate law: for p not an even whole number, Lee-Yang holds only if a <= p, i.e. n <= 1 + p^2.
#   Why: far out, G is a wave from the rim, ~ s^-(a+1), plus a steady push from the centre, where |x|^p has a kink,
#   ~ -a Gamma(p+1) cos(pi(p+1)/2) s^-(p+1).  If a > p the push wins, G stops changing sign, and the zeros leave the line.
#   At even whole p the kink is absent and nothing pushes.
#   The earlier quadrature scan (checks.py J3) looked only at |Re s| < 28 and missed this; here G is summed exactly.
# ---------------------------------------------------------------------------
import mpmath as mp
def ly_coeffs(p, a, K, dps):
    mp.mp.dps = dps; p = mp.mpf(p); a = mp.mpf(a)
    return [(-1)**k/mp.factorial(2*k)*mp.beta((2*k + 1)/p, a + 1)/p for k in range(K)]
def ly_eval(c, s):
    s2 = mp.mpmathify(s)**2; acc = mp.mpf(0)
    for ck in reversed(c): acc = acc*s2 + ck
    return acc
def far_values(p, a, S0, S1, step=0.25):
    c = ly_coeffs(p, a, int(1.5*S1) + 60, int(S1/2) + 40)
    return [ly_eval(c, S0 + i*step) for i in range(int((S1 - S0)/step) + 1)]
def changes(v): return sum(1 for i in range(len(v) - 1) if v[i]*v[i + 1] < 0)
def crossover(p, a):
    """s where the centre's push matches the rim's wave: beyond it, G keeps one sign (only when a > p)."""
    push = a*math.gamma(p + 1)*abs(math.sin(math.pi*p/2)); wave = math.gamma(a + 1)*p**a
    return (wave/push)**(1/(a - p))
expected = round(60/np.pi)
t0 = time.perf_counter()
far = {}
for p in [1.0, 1.25, 1.5, 1.8, 3.0, 5.0]:
    for da in [-0.5, -0.1]:
        far[f"p={p:g}, a={p + da:g}"] = dict(window=[150, 210], changes=changes(far_values(p, p + da, 150, 210)))
    for da in [0.5, 1.0, 2.0]:
        a = p + da; sx = crossover(p, a)
        if 1.3*sx + 60 > 700:
            far[f"p={p:g}, a={a:g}"] = dict(window=None, crossover=sx); continue
        S0 = max(20.0, 1.3*sx); v = far_values(p, a, S0, S0 + 60)
        sign_pred = math.copysign(1, math.sin(math.pi*p/2))
        far[f"p={p:g}, a={a:g}"] = dict(window=[S0, S0 + 60], crossover=sx, changes=changes(v),
                                        sign_matches_push=all(math.copysign(1, float(x)) == sign_pred for x in v))
for a in [2.0, 6.0, 10.0]:
    for p in [2.0, 4.0]:
        far[f"p={p:g}, a={a:g}"] = dict(window=[150, 210], changes=changes(far_values(p, a, 150, 210)))
far_time = time.perf_counter() - t0
R["lee_yang_far"] = far; R["lee_yang_far_expected"] = expected
for k_, v in far.items(): print("     ", k_, v)
pv = lambda k_: float(k_.split(",")[0][2:]); av = lambda k_: float(k_.split("a=")[1])
odd = [k_ for k_ in far if pv(k_) not in (2.0, 4.0)]
ok_below = all(far[k_]["changes"] >= expected - 1 for k_ in odd if av(k_) < pv(k_))
tested_above = [k_ for k_ in odd if av(k_) > pv(k_) and far[k_]["window"] is not None]
ok_above = all(far[k_]["changes"] == 0 and far[k_]["sign_matches_push"] for k_ in tested_above)
check("N3a below a = p the zeros keep coming on the line far out; above a = p, past the predicted crossover, G keeps the "
      "sign the centre's kink predicts and has no real zeros", ok_below and ok_above and len(tested_above) >= 6,
      f"{len(tested_above)} cases above a = p tested past their crossover")
check("N3b at even whole p (2 and 4) the zeros keep coming for every a tried (2, 6, 10): no kink, no push",
      all(far[k_]["changes"] >= expected - 1 for k_ in far if pv(k_) in (2.0, 4.0)), "")
# full count inside a disc (argument principle on |s| = X) against real zeros on [-X, X], for a <= p: no zero off the line
def zeros_in_disc(p, a, X=60.0, npts=2400):
    c = ly_coeffs(p, a, int(1.5*X) + 80, int(X/2) + 50)
    zs = [X*mp.expj(2*mp.pi*k/npts) for k in range(npts + 1)]
    vals = [ly_eval(c, z) for z in zs]
    ang = [float(mp.arg(v)) for v in vals]
    tot = float(np.sum(np.diff(np.unwrap(ang))))/(2*np.pi)
    rv = [ly_eval(c, -X + i*0.05) for i in range(int(2*X/0.05) + 1)]
    real = sum(1 for i in range(len(rv) - 1) if rv[i]*rv[i + 1] < 0)
    return round(tot), real
with timed("N3 polish (exact zero count in a disc)"):
    disc = {f"p={p}, a={a:g}": zeros_in_disc(p, a) for p, a in [(1.0, 0.9), (1.5, 1.4), (1.8, 1.7), (3.0, 2.9), (1.5, 2.5)]}
R["lee_yang_disc"] = disc
print("    zeros in |s| < 60 (all, real):", disc)
check("N3c polish: for a <= p every zero inside |s| < 60 is real; for a > p (p = 1.5, a = 2.5) some are not",
      all(v[0] == v[1] for k_, v in disc.items() if k_ != "p=1.5, a=2.5") and disc["p=1.5, a=2.5"][0] > disc["p=1.5, a=2.5"][1], "")
# in spin terms: uniform n-dim p-ball, a = (n - 1)/p <= p  <=>  n <= 1 + p^2
R["lee_yang_max_n"] = {str(p): int(np.floor(1 + p*p)) for p in [1.0, 1.25, 1.5, 1.8, 3.0]}
# ---------------------------------------------------------------------------
# N4 (shape dial, round <-> triangle): the in-between laws found are the turn (simplex page) and the ordering fraction;
# neither gives a new closed form.  Recorded from checks.py (no recomputation here).
# N5 (roundness dial): in-between laws that are immediate from known formulas:
#   volume sums E_{1/p}(2 Gamma(1 + 1/p)) (from the known volume formula); replicator at any p is the same law
#   (u = x^(1/p) is only a change of variable).  Verify the replicator statement.
# ---------------------------------------------------------------------------
f_ = rng.normal(size=5); x0 = rng.dirichlet(np.ones(5))
def rep(t, x): return x*(f_ - x @ f_)
from scipy.integrate import solve_ivp
solx = solve_ivp(rep, [0, 2], x0, rtol=1e-11, atol=1e-13).y[:, -1]
gaps = []
for p in [1.3, 2.0, 3.7]:
    def up(t, u, p=p):
        x = u**p; return u*(f_ - x @ f_)/p
    solu = solve_ivp(up, [0, 2], x0**(1/p), rtol=1e-11, atol=1e-13).y[:, -1]
    gaps.append(float(np.max(np.abs(solu**p - solx))))
check("N5 the replicator read through u = x^(1/p) is the same law at every p (no new in-between law)",
      max(gaps) < 1e-8, f"max gap {max(gaps):.1e}")

R["polish_seconds"] = POL; R["lee_yang_far_seconds"] = far_time; R["bell_grid_seconds"] = grid_time
print("    polish cost (s):", {k_: round(v, 3) for k_, v in POL.items()}, "total", round(sum(POL.values()), 2))
json.dump(R, open("new_laws.json", "w"), indent=1, default=float)
print(f"\n{sum(ok for _, ok in RESULTS)} of {len(RESULTS)} checks pass")
