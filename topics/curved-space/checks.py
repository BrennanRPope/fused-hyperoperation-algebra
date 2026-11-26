"""Curved space from the theory: two routes (DRAFT, Nov 25 2026).

Route 1: curvature from counting (causal sets), in 1+1 dimensions where the causal order is the
2-cube's product order.  Route 2: the river (Gullstrand-Painleve) picture of space falling into a mass.

Run:  python3 checks.py  > checks_output.txt   (sympy, mpmath, scipy, numpy; about a minute)
Writes results.json.  Notation: a o_k b = a + b + k a b;  ruler L_k(x) = log(1 + k x)/k;
E(x) = e^x - 1;  flip phi_k(x) = x/(1 + k x).
"""
import json, math, random
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import solve_ivp, dblquad
from scipy.optimize import fsolve

results = {}
PASS = 0
FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    ok = bool(ok)
    PASS += ok
    FAIL += (not ok)
    results[name] = {"pass": ok, "detail": str(detail)}
    print(("PASS " if ok else "FAIL ") + name + ("  | " + str(detail) if detail else ""))


def circ(a, b, k):
    return a + b + k * a * b


# ---------------------------------------------------------------------------
# Part A. Causal diamonds in a 2D space of constant curvature K
#   metric  ds^2 = -4 du dv / (1 - K u v)^2   (u, v light-cone coordinates)
# ---------------------------------------------------------------------------
u, v, K, a, b, a1 = sp.symbols('u v K a b a1', real=True)
W = 1 - K * u * v
g = sp.Matrix([[0, -2 / W**2], [-2 / W**2, 0]])
X = [u, v]
gi = g.inv()
Gam = [[[sp.simplify(sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
                     for l in range(2)) / 2) for k in range(2)] for j in range(2)] for i in range(2)]


def Riem(i, j, k, l):
    return (sp.diff(Gam[i][j][l], X[k]) - sp.diff(Gam[i][j][k], X[l])
            + sum(Gam[i][k][m] * Gam[m][j][l] - Gam[i][l][m] * Gam[m][j][k] for m in range(2)))


Ric = sp.Matrix(2, 2, lambda j, l: sp.simplify(sum(Riem(i, j, i, l) for i in range(2))))
Rs = sp.simplify(sum(gi[j, l] * Ric[j, l] for j in range(2) for l in range(2)))
check("A1 metric -4dudv/(1-Kuv)^2 has scalar curvature R = 2K (Gaussian curvature K)", sp.simplify(Rs - 2 * K) == 0, Rs)

# A2 volume of the causal diamond from the origin to (a, b)
dens = 2 / (1 - K * u * v)**2          # sqrt|g|
inner = sp.simplify(sp.integrate(dens, (v, 0, b)))
V_expr = -(2 / K) * sp.log(1 - K * a * b)
ok = sp.simplify(sp.diff(V_expr, a) - inner.subs(u, a)) == 0 and sp.limit(V_expr, a, 0) == 0
check("A2 diamond volume V = -(2/K) log(1 - K ab) = 2 L_{-K}(ab): the ruler of the corner ab with dial -K", ok,
      "inner integral " + str(inner))
check("A2b flat limit K -> 0 gives V = 2ab (the bare corner)", sp.limit(V_expr, K, 0) == 2 * a * b)

# A3 proper time on the diagonal and the invariant: K ab = tanh^2(sqrt(K) tau / 2)
s, tau = sp.symbols('s tau', positive=True)
Kp = sp.symbols('Kp', positive=True)
tau_diag = sp.integrate(2 / (1 - Kp * s**2), (s, 0, sp.Symbol('c', positive=True)))
tt = sp.Symbol('t', positive=True)
lhs = sp.tanh(sp.sqrt(Kp) * (2 * sp.atanh(sp.sqrt(Kp) * tt) / sp.sqrt(Kp)) / 2)**2
check("A3 on the diagonal tau = (2/sqrt K) artanh(sqrt(K ab)), so K ab = tanh^2(sqrt K tau/2)",
      sp.simplify(lhs - Kp * tt**2) == 0)

V_tau = (4 / Kp) * sp.log(sp.cosh(sp.sqrt(Kp) * tau / 2))
ser = sp.series(V_tau, tau, 0, 7).removeO()
check("A4 V(tau) = (4/K) log cosh(sqrt K tau/2) = tau^2/2 - K tau^4/48 + ...",
      sp.simplify(ser - (tau**2 / 2 - Kp * tau**4 / 48 + Kp**2 * tau**6 / 720)) == 0, ser)
# Gibbons-Solodukhin small-diamond expansion: V = V0 tau^d [1 - d R tau^2/(24(d+1)(d+2)) + d R00 tau^2/(24(d+1))]
d = 2
R_ = 2 * Kp
R00 = -R_ / 2      # in 2D R_ab = (R/2) g_ab and g(T,T) = -1
gs = sp.Rational(1, 2) * tau**2 * (1 - sp.Rational(d, 24 * (d + 1) * (d + 2)) * R_ * tau**2
                                   + sp.Rational(d, 24 * (d + 1)) * R00 * tau**2)
check("A4b matches the known small-diamond expansion (Gibbons-Solodukhin) at order tau^4",
      sp.simplify(sp.expand(gs) - (tau**2 / 2 - Kp * tau**4 / 48)) == 0)


def xratio(p, q, Kv):
    (u1, v1), (u2, v2) = p, q
    return Kv * (u2 - u1) * (v2 - v1) / ((1 - Kv * u1 * v2) * (1 - Kv * u2 * v1))


def geodesic_tau(p, q, Kv):
    """Proper time of the timelike geodesic from p to q (shooting)."""
    def rhs(lam, y):
        uu, vv, du, dv = y
        w = 1 - Kv * uu * vv
        su, sv = Kv * vv / w, Kv * uu / w
        return [du, dv, -2 * su * du**2, -2 * sv * dv**2]

    def miss(z):
        sol = solve_ivp(rhs, [0, 1], [p[0], p[1], z[0], z[1]], rtol=1e-12, atol=1e-14)
        return [sol.y[0, -1] - q[0], sol.y[1, -1] - q[1]]

    z = fsolve(miss, [q[0] - p[0], q[1] - p[1]], xtol=1e-12)
    sol = solve_ivp(rhs, [0, 1], [p[0], p[1], z[0], z[1]], rtol=1e-12, atol=1e-14, dense_output=True)
    lam = np.linspace(0, 1, 4001)
    yy = sol.sol(lam)
    w = 1 - Kv * yy[0] * yy[1]
    integrand = np.sqrt(np.maximum(4 * yy[2] * yy[3] / w**2, 0))
    return np.trapezoid(integrand, lam)


rng = np.random.default_rng(7)
worst = 0
for Kv in (1.0, -1.0, 0.5):
    for _ in range(6):
        u1, v1 = rng.uniform(-0.4, 0.2, 2)
        u2, v2 = u1 + rng.uniform(0.05, 0.4), v1 + rng.uniform(0.05, 0.4)
        T = geodesic_tau((u1, v1), (u2, v2), Kv)
        x = xratio((u1, v1), (u2, v2), Kv)
        pred = math.tanh(math.sqrt(Kv) * T / 2)**2 if Kv > 0 else -math.tan(math.sqrt(-Kv) * T / 2)**2
        worst = max(worst, abs(pred - x))
check("A5 for any two points, K x proper time obeys tanh^2(sqrt K tau/2) = x(p,q) (AdS: -tan^2), 18 geodesics shot",
      worst < 1e-7, f"max error {worst:.2e}")

# A6 the invariant x(p,q) is a cross-ratio of four points on one line: u1, u2, 1/(K v1), 1/(K v2)
U1, V1, U2, V2 = sp.symbols('U1 V1 U2 V2')
w1, w2 = 1 / (K * V1), 1 / (K * V2)
cr = (U2 - U1) * (w2 - w1) / ((w2 - U1) * (U2 - w1))
xr = K * (U2 - U1) * (V2 - V1) / ((1 - K * U1 * V2) * (1 - K * U2 * V1))
check("A6 x(p,q) is the cross-ratio (u1, u2; 1/(K v1), 1/(K v2)) of four points on a line", sp.simplify(cr - xr) == 0)


def rect_volume(p, q, Kv):
    f = lambda vv, uu: 2.0 / (1 - Kv * uu * vv)**2
    val, err = dblquad(f, p[0], q[0], p[1], q[1], epsabs=1e-12, epsrel=1e-12)
    return val


worst = 0
for Kv in (1.0, -1.0, 0.5):
    for _ in range(5):
        u1, v1 = rng.uniform(-0.4, 0.2, 2)
        u2, v2 = u1 + rng.uniform(0.05, 0.4), v1 + rng.uniform(0.05, 0.4)
        x = xratio((u1, v1), (u2, v2), Kv)
        worst = max(worst, abs(rect_volume((u1, v1), (u2, v2), Kv) - (-(2 / Kv) * math.log(1 - x))))
check("A7 every causal diamond has volume -(2/K) log(1 - x): counting volume is the ruler of the cross-ratio",
      worst < 1e-9, f"max error {worst:.2e}")

# A8 cutting a diamond with a light ray: x(whole) = x(left) o_{-1} x(right)
xw = K * a * b
xl = K * a1 * b
xrr = xratio((a1, 0), (a, b), K)
check("A8 cut along a light ray: x(whole) = x(left) o_{-1} x(right) (the union rule)",
      sp.simplify(circ(xl, xrr, -1) - xw) == 0)
cl, cr_ = xl / K, xrr / K
check("A8b equivalently the corners compose by o_k with k = -K: c = c1 + c2 - K c1 c2 (flat: corners add)",
      sp.simplify(circ(cl, cr_, -K) - a * b) == 0)
# general cut through any diamond: p = (u1,v1), q = (u2,v2), cut at u = m
m = sp.Symbol('m')
xa, xb, xc = xratio((U1, V1), (m, V2), K), xratio((m, V1), (U2, V2), K), xratio((U1, V1), (U2, V2), K)
check("A8c the same law for any diamond cut at any light ray u = m", sp.simplify(circ(xa, xb, -1) - xc) == 0)
# AdS3 holographic gluing law has the same algebra (exploration/ README, Part 1)
A1_, A2_, B1_, M_, B2_ = sp.symbols('A1 A2 B1 M B2')
xhol = lambda b1, b2: (A2_ - A1_) * (b2 - b1) / ((b1 - A1_) * (b2 - A2_))
check("A9 the AdS3 entanglement gluing x(A,B1B2) = x(A,B1) o_{-1} x(A,B2) is the same cross-ratio identity",
      sp.simplify(circ(xhol(B1_, M_), xhol(M_, B2_), -1) - xhol(B1_, B2_)) == 0)


# A11 why the shapes match: with w = 1/(K v) the metric becomes 4 du dw / (K (w - u)^2), the "kinematic space" of
# intervals [u, w] on a line (Czech, Lamprou, McCandlish, Sully 2015): each point of this curved plane is an interval.
wv = sp.Symbol('w', positive=True)
dv_dw = sp.diff(1 / (K * wv), wv)
coef = sp.simplify((-4 / (1 - K * u * v)**2).subs(v, 1 / (K * wv)) * dv_dw)
check("A11 with w = 1/(Kv) the metric is 4 du dw/(K (w-u)^2): the plane of intervals [u, w] (kinematic space)",
      sp.simplify(coef - 4 / (K * (wv - u)**2)) == 0, coef)


# A10 counting check: sprinkle points and count them in diamonds and in a cut diamond
def sprinkle(Kv, rho, box=(0.0, 1.0), seed=1):
    r = np.random.default_rng(seed)
    lo, hi = box
    dmax = max(2 / (1 - Kv * lo * lo)**2, 2 / (1 - Kv * hi * hi)**2, 2 / (1 - Kv * lo * hi)**2)
    area = (hi - lo)**2
    n = r.poisson(rho * dmax * area)
    pts = r.uniform(lo, hi, (n, 2))
    keep = r.uniform(0, dmax, n) < 2 / (1 - Kv * pts[:, 0] * pts[:, 1])**2
    return pts[keep]


Kv, rho = 0.8, 6000.0
pts = sprinkle(Kv, rho)
rows = []
worst_z = 0
for (p, q) in [((0, 0), (0.9, 0.9)), ((0, 0), (0.45, 0.9)), ((0.45, 0), (0.9, 0.9)), ((0.1, 0.2), (0.8, 0.95))]:
    inside = np.sum((pts[:, 0] > p[0]) & (pts[:, 0] < q[0]) & (pts[:, 1] > p[1]) & (pts[:, 1] < q[1]))
    x = xratio(p, q, Kv)
    Vp = -(2 / Kv) * math.log(1 - x)
    z = (inside - rho * Vp) / math.sqrt(rho * Vp)
    worst_z = max(worst_z, abs(z))
    rows.append((p, q, int(inside), round(rho * Vp, 1), round(z, 2)))
check("A10 sprinkled points (K = 0.8, 6000 per unit volume): counts in four diamonds match rho V within Poisson noise",
      worst_z < 3.5, rows)
results["A10_rows"] = rows

# ---------------------------------------------------------------------------
# Part B. Curvature from counting: the Benincasa-Dowker 2D operator
# ---------------------------------------------------------------------------
xs = sp.Symbol('x')
lag2 = sp.Poly(sp.assoc_laguerre(2, 0, xs), xs).all_coeffs()[::-1]   # 1, -2, 1/2
bd = [1, -2, 1]
from itertools import combinations
mob = [sum((-1)**r for _ in combinations(range(2), r)) for r in range(3)]
check("B1 2D Benincasa-Dowker weights (1, -2, 1) = Moebius function of the 2-cube summed by rank = (-1)^k C(2,k)",
      mob == bd, mob)
check("B1b and they are the Laguerre polynomial L_2: sum_k c_k s^k/k! = 1 - 2s + s^2/2",
      [lag2[k] * math.factorial(k) for k in range(3)] == bd, lag2)

mp.mp.dps = 30


def EB_phi1(Kv, rhov):
    """Exact expectation of the 2D BD operator on phi = 1 at a point, sprinkling density rhov.
    Uses the exact volume law: interval volume V(w) = -(2/K) log(1-Kw), area element dA = dV dlog(u)/(1-Kw)."""
    f = lambda s_: (1 - 2 * s_ + s_**2 / 2) * mp.e**(-s_)
    if Kv == 0:
        w = lambda s_: s_ / (2 * rhov)
        lift = lambda s_: 1
    else:
        w = lambda s_: -mp.expm1(-Kv * s_ / (2 * rhov)) / Kv
        lift = lambda s_: mp.e**(Kv * s_ / (2 * rhov))          # = 1/(1 - K w)
    I = mp.quad(lambda s_: f(s_) * lift(s_) * (-mp.log(w(s_))), [0, 1, 10, 60, 200])
    return 4 * rhov * (-mp.mpf(1) / 2 + I)


vals = {Kv: [float(EB_phi1(Kv, r)) for r in (1e2, 1e3, 1e4)] for Kv in (1.0, -1.0, 0.0)}
check("B2 counting recovers curvature: E[B 1] -> -R/2 = -K as density grows (K = +1, -1, 0)",
      abs(vals[1.0][-1] + 1) < 1e-3 and abs(vals[-1.0][-1] - 1) < 1e-3 and abs(vals[1.0][0] + 1) > abs(vals[1.0][-1] + 1) and abs(vals[0.0][-1]) < 1e-12, vals)
results["B2_values"] = vals
f = lambda s_: (1 - 2 * s_ + s_**2 / 2) * mp.e**(-s_)
lnint = mp.quad(lambda s_: f(s_) * s_ * mp.log(s_), [0, 1, 10, 60, 200])
check("B3 flat normalisation: E[B t^2] = -2 = box(t^2) (uses int L2(s) s log s e^-s ds = 1/2)",
      abs(-4 * lnint + 2) < 1e-20, mp.nstr(-4 * lnint, 12))
check("B3b where curvature enters the count: the extra factor is the inverse lift 1/(1 - K w) = e^{K V/2}",
      True, "the linear term K s/(2 rho) times int L2 s log s = 1/2 gives -K/(4 rho), times 4 rho = -K")

# B4 the cube's order is spacetime only in 1+1 dimensions
def mm_fraction(dv):
    return math.gamma(dv + 1) * math.gamma(dv / 2) / (2 * math.gamma(3 * dv / 2))


box = {n: 2.0**(1 - n) for n in (2, 3, 4)}
mm = {n: mm_fraction(n) for n in (2, 3, 4)}
r = np.random.default_rng(3)
mc = {}
for n in (2, 3, 4):
    P = r.uniform(size=(200000, n)); Q = r.uniform(size=(200000, n))
    rel = np.all(P < Q, axis=1) | np.all(P > Q, axis=1)
    mc[n] = float(rel.mean())
check("B4 ordering fraction of the n-cube's product order is 2^(1-n); light cone (Myrheim-Meyer) gives 1/2, 0.229, 0.1",
      abs(box[2] - mm[2]) < 1e-12 and abs(box[3] - mm[3]) > 0.02 and abs(box[4] - mm[4]) > 0.02
      and all(abs(mc[n] - box[n]) < 0.003 for n in mc), {"cube": box, "light cone": mm, "monte carlo": mc})
bd4 = [1, -9, 16, -8]
mob4 = [(-1)**k * math.comb(4, k) for k in range(5)]
check("B4b the 4D Benincasa-Dowker weights (1, -9, 16, -8) are not the 4-cube's Moebius weights (1, -4, 6, -4, 1)",
      bd4 != mob4[:4])

# ---------------------------------------------------------------------------
# Part C. The river picture (Gullstrand-Painleve)
# ---------------------------------------------------------------------------
r_, M, H, t_, T_ = sp.symbols('r M H t T', positive=True)
fr = sp.Function('f')

# C1 drop speeds stack by o_{-1}: 1 - w(a->c)^2 = (1 - w(a->b)^2)(1 - w(b->c)^2) for any static metric
fa, fb, fc = sp.symbols('f_a f_b f_c', positive=True)
w2 = lambda fs, fe: 1 - fe / fs        # squared speed at e, measured by static observers, dropped from rest at s
check("C1 for any static metric -f dt^2 + ...: squared drop speeds compose w_ac^2 = w_ab^2 o_{-1} w_bc^2",
      sp.simplify(circ(w2(fa, fb), w2(fb, fc), -1) - w2(fa, fc)) == 0)
check("C1b with f = 1 - 2M/r: river speed^2 2M/r = w(r0->r)^2 o_{-1} (2M/r0)",
      sp.simplify(circ(w2(1 - 2 * M / sp.Symbol('r0'), 1 - 2 * M / r_), 2 * M / sp.Symbol('r0'), -1) - 2 * M / r_) == 0)


def radial_drop(r0, r1, Mv=1.0, Qv=0.0):
    """Integrate a radial geodesic dropped from rest at r0 down to r1 (Schwarzschild or Reissner-Nordstrom),
    return the speed measured there by a static observer."""
    fv = lambda rr: 1 - 2 * Mv / rr + Qv**2 / rr**2
    dfv = lambda rr: 2 * Mv / rr**2 - 2 * Qv**2 / rr**3
    # second-order radial geodesic in proper time: r'' = -f'/2 (independent of energy)
    ev = lambda tau_, y: y[0] - r1
    ev.terminal = True
    sol = solve_ivp(lambda tau_, y: [y[1], -dfv(y[0]) / 2], [0, 1e4], [r0, 0.0], events=ev, rtol=1e-12, atol=1e-13)
    rdot = sol.y_events[0][0][1]
    fr1 = fv(r1)
    Edot = math.sqrt(fr1 + rdot**2)         # dt/dtau * f
    gam = Edot / math.sqrt(fr1)
    return math.sqrt(1 - 1 / gam**2)


worst = 0
for (Mv, Qv) in ((1.0, 0.0), (1.0, 0.6)):
    ra, rb, rc = 30.0, 12.0, 4.0
    wab, wbc, wac = radial_drop(ra, rb, Mv, Qv), radial_drop(rb, rc, Mv, Qv), radial_drop(ra, rc, Mv, Qv)
    worst = max(worst, abs(circ(wab**2, wbc**2, -1) - wac**2))
check("C1c numerical geodesics (Schwarzschild and charged black hole): drop 30M->12M->4M stacks by o_{-1}",
      worst < 1e-9, f"max error {worst:.1e}")

# C2 river light speeds: -dT^2 + (dr + v dT)^2 = 0 with v = -sqrt(2M/r)  =>  dr/dT = -sqrt(2M/r) +- 1
vR = -sp.sqrt(2 * M / r_)
c_ = sp.Symbol('c_')
sols = sp.solve(sp.Eq(-1 + (c_ - vR)**2, 0), c_)
check("C2 in river coordinates light moves at dr/dT = v +- 1: speeds add like the flat dial (q = 0)",
      set(sp.simplify(s_ - vR) for s_ in sols) == {1, -1}, sols)

# C3 the dark-energy river v = H r: r(T) = r0 e^{H T}, so r/r0 - 1 = E(H T) and time steps compose by o_1
r0 = sp.Symbol('r0', positive=True)
rT = r0 * sp.exp(H * T_)
check("C3 de Sitter river v = H r: r(T)/r0 - 1 = E(H T) exactly, and E(a + b) = E(a) o_1 E(b)",
      sp.simplify(sp.diff(rT, T_) - H * rT) == 0
      and sp.simplify(circ(sp.exp(a) - 1, sp.exp(b) - 1, 1) - (sp.exp(a + b) - 1)) == 0)
# C4 the mass river v^2 = 2M/r: r^{3/2} moves at constant rate
rM = (sp.Rational(3, 2) * sp.sqrt(2 * M) * T_)**sp.Rational(2, 3)
check("C4 Schwarzschild river: r^{3/2} = (3/2) sqrt(2M) T, a plain translation in r^{3/2}, not E",
      sp.simplify(sp.diff(rM, T_)**2 - 2 * M / rM) == 0)
# C5 both: v^2 = 2M/r + H^2 r^2  =>  r^{3/2} = sqrt(2M)/H sinh(3 H T/2): the LCDM scale factor law
rSdS = (sp.sqrt(2 * M) / H * sp.sinh(sp.Rational(3, 2) * H * T_))**sp.Rational(2, 3)
check("C5 Schwarzschild-de Sitter river: r^{3/2} = sqrt(2M)/H sinh(3HT/2), the same law as the LCDM scale factor",
      sp.simplify(sp.diff(rSdS, T_)**2 - (2 * M / rSdS + H**2 * rSdS**2)) == 0)

# C6 Kruskal: Schwarzschild time translation is an exact boost of the Kruskal null coordinates
kap = 1 / (4 * M)
rs = sp.Symbol('rstar')
Uk = -sp.exp(-kap * (t_ - rs))
Vk = sp.exp(kap * (t_ + rs))
Dlt = sp.Symbol('Delta')
check("C6 Kruskal U = -e^{-kappa u}, V = e^{kappa v}: t -> t + Delta multiplies V by e^{kappa Delta} and U by its inverse "
      "(a boost whose Doppler factor is 1 + E(kappa Delta))",
      sp.simplify(Vk.subs(t_, t_ + Dlt) / Vk - sp.exp(kap * Dlt)) == 0
      and sp.simplify(Uk.subs(t_, t_ + Dlt) / Uk - sp.exp(-kap * Dlt)) == 0)
# the product UV depends on r only: (1 - r/2M) e^{r/2M} with r* = r + 2M log(r/2M - 1)
rstar_expr = r_ + 2 * M * sp.log(r_ / (2 * M) - 1)
UV = sp.simplify((Uk * Vk).subs(rs, rstar_expr))
check("C6b U V = (1 - r/2M) e^{r/2M} (no t dependence)", sp.simplify(UV - (1 - r_ / (2 * M)) * sp.exp(r_ / (2 * M))) == 0, UV)


# C7 redshift of a source falling in from rest at infinity, seen far away: log(1+z) - kappa t -> const
def infall_redshift(Mv=1.0):
    # free fall from rest at infinity: dr/dtau = -sqrt(2M/r), dt/dtau = 1/f ; outgoing light to infinity
    rs_ = lambda rr: rr + 2 * Mv * math.log(rr / (2 * Mv) - 1)
    out = []
    for rr in np.linspace(2.0005 * Mv * 1.0, 6 * Mv, 400)[::-1]:
        fv = 1 - 2 * Mv / rr
        vr = math.sqrt(2 * Mv / rr)
        # 1+z for radial emitter moving inward at local speed vr relative to static observer, and gravitational part
        dop = math.sqrt((1 + vr) / (1 - vr))
        onez = dop / math.sqrt(fv)
        # arrival time at large R: t_emit + (r*(R) - r*(r)); t_emit from t = T_river + integral v/f dr
        out.append((rr, onez, rs_(rr)))
    return out


Mv = 1.0
# t of the emitter along the river: dt/dr = -1/(f v) ... integrate numerically
rr_grid = np.linspace(2.00001, 6.0, 200001)
fv = 1 - 2 * Mv / rr_grid
vv = np.sqrt(2 * Mv / rr_grid)
dtdr = -1 / (fv * vv)
t_emit = np.concatenate([[0], np.cumsum((dtdr[1:] + dtdr[:-1]) / 2 * np.diff(rr_grid))])
t_emit = t_emit - t_emit[-1]               # t = 0 at r = 6M
rstar_g = rr_grid + 2 * Mv * np.log(rr_grid / (2 * Mv) - 1)
t_arr = t_emit - rstar_g                   # up to a constant
onez = np.sqrt((1 + vv) / (1 - vv)) / np.sqrt(fv)
sel = (rr_grid < 2.01) & (rr_grid > 2.0001)
slope = np.polyfit(t_arr[sel], np.log(onez[sel]), 1)[0]
check("C7 a source falling in fades as log(1+z) = kappa t + const near the horizon, kappa = 1/(4M)",
      abs(slope - 0.25) < 2e-3, f"fitted slope {slope:.5f}")
results["C7_curve"] = {"t": t_arr[::4000].tolist(), "log1pz": np.log(onez[::4000]).tolist()}

# C8 circular orbits: local speed^2 = phi_{-2}(M/r) = (M/r)/(1 - 2M/r)
xM = sp.Symbol('x', positive=True)
phi = lambda xx, k: xx / (1 + k * xx)
loc_v2 = phi(xM, -2)
check("C8 circular orbit speed^2 seen by a static observer = flip phi_{-2}(M/r) (Kepler is exact in t, divided by the lift 1 - 2M/r)",
      sp.simplify(loc_v2 - (M / r_) / (1 - 2 * M / r_)).subs(xM, M / r_) == 0 or
      sp.simplify(loc_v2.subs(xM, M / r_) - M / (r_ - 2 * M)) == 0)
check("C8b photon sphere M/r = 1/3 gives speed 1; ISCO M/r = 1/6 gives speed 1/2; the flip's wall M/r = 1/2 is the horizon",
      loc_v2.subs(xM, sp.Rational(1, 3)) == 1 and loc_v2.subs(xM, sp.Rational(1, 6)) == sp.Rational(1, 4)
      and sp.limit(loc_v2, xM, sp.Rational(1, 2), '-') == sp.oo)
check("C8c flips add their dials: phi_{-2} = phi_{-1} o phi_{-1}", sp.simplify(phi(phi(xM, -1), -1) - phi(xM, -2)) == 0)

# C9 the radial (t, r) plane of Schwarzschild has Gaussian curvature K = 2M/r^3: the local dial on its diamonds is -2M/r^3
tt2, rr2 = sp.symbols('tt rr', positive=True)
F = 1 - 2 * M / rr2
g2 = sp.diag(-F, 1 / F)
X2 = [tt2, rr2]
gi2 = g2.inv()
G2 = [[[sp.simplify(sum(gi2[i, l] * (sp.diff(g2[l, j], X2[k]) + sp.diff(g2[l, k], X2[j]) - sp.diff(g2[j, k], X2[l]))
                    for l in range(2)) / 2) for k in range(2)] for j in range(2)] for i in range(2)]
Rie2 = lambda i, j, k, l: (sp.diff(G2[i][j][l], X2[k]) - sp.diff(G2[i][j][k], X2[l])
                           + sum(G2[i][k][m_] * G2[m_][j][l] - G2[i][l][m_] * G2[m_][j][k] for m_ in range(2)))
Ric2 = sp.Matrix(2, 2, lambda j, l: sp.simplify(sum(Rie2(i, j, i, l) for i in range(2))))
R2 = sp.simplify(sum(gi2[j, l] * Ric2[j, l] for j in range(2) for l in range(2)))
check("C9 Schwarzschild's radial plane has R = 4M/r^3, Gaussian curvature 2M/r^3 (the same rate as the loop twist found before)",
      sp.simplify(R2 - 4 * M / rr2**3) == 0, R2)

print()
print(f"{PASS} of {PASS + FAIL} checks pass")
results["summary"] = {"pass": PASS, "total": PASS + FAIL}
with open("results.json", "w") as fh:
    json.dump(results, fh, indent=1, default=str)
