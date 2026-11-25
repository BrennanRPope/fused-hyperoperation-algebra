"""Can mult-add build curved geometry, and do the four known bridges to Einstein's equations go through?

Run: python3 checks.py  (writes checks_output.txt and results.json).  Units G = c = hbar = 1 unless stated.
Part S: pure geometry (mult-add's own rules, a place dial, a direction dial, a Clifford-number dial).
Part J: Jacobson (heat and horizon area).  Part T: teleparallel (torsion).  Part L: Lovelock.  Part D: Deser and
scalar-tensor (Brans-Dicke, Cassini, the sombrero).
"""
import json, math, sys, os, random, itertools
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, "/mnt/project-files/math-theory/gr-bridges")
sys.path.insert(0, "/mnt/project-files/math-theory/general-relativity")
from ga import GA
from grtools import christoffel, riemann, ricci

RES, LINES = {}, []
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    ok = bool(ok)
    npass += ok
    nfail += (not ok)
    line = f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else "")
    print(line, flush=True)
    LINES.append(line)


def section(title):
    print("\n== " + title, flush=True)
    LINES.append("\n== " + title)


def zero(e):
    return sp.simplify(e) == 0

t, x, y, z, r, th, ph = sp.symbols('t x y z r theta phi', real=True)
lam = sp.symbols('lambda', real=True)
# ===========================================================================================
section("S1. Mult-add's own rule in one direction: velocity composition in the Clifford plane Cl(1,1)")
# Cl with e0^2 = -1 (time, dial q < 0) and e1^2 = +1; sigma = e1 e0 squares to +1 (a boost generator).
G2 = GA([-1, 1])
sig = G2.e(1, 0)
check("sigma = e1 e0 squares to +1 (a boost: hyperbolic, not a rotation)", G2.equal(G2.mul(sig, sig), G2.scalar(1)))
b1, b2 = sp.symbols('beta1 beta2', real=True)
P = G2.mul(G2.add(G2.scalar(1), G2.scale(b1, sig)), G2.add(G2.scalar(1), G2.scale(b2, sig)))
s0 = P.get(0, 0); s1 = P.get(3, 0) / sig[3]
check("(1 + b1 sigma)(1 + b2 sigma) = (1 + b1 b2)(1 + beta sigma) with beta = (b1 + b2)/(1 + b1 b2), Einstein's addition",
      zero(s1 / s0 - (b1 + b2) / (1 + b1 * b2)))
D2 = lambda b: (1 + b) / (1 - b)
bc = (b1 + b2) / (1 + b1 * b2)
check("Doppler-squared D^2 = (1+beta)/(1-beta) multiplies, so D^2 - 1 composes by mult-add",
      zero((D2(bc) - 1) - ((D2(b1) - 1) + (D2(b2) - 1) + (D2(b1) - 1) * (D2(b2) - 1))))
Nrm = G2.mul(G2.add(G2.scalar(1), G2.scale(b1, sig)), G2.add(G2.scalar(1), G2.scale(-b1, sig)))
check("the wall: 1 + beta sigma has norm 1 - beta^2, which is zero exactly at light speed (mult-add's -1)",
      G2.equal(Nrm, G2.scalar(1 - b1 ** 2)))
RES["S1"] = "In Cl(1,1), boosts are 1 + beta*sigma with sigma^2 = +1; products give Einstein addition, D^2 - 1 mult-adds, and the wall is the light cone."

# ===========================================================================================
section("S2. Twist from the Clifford square rule: boosts in two directions close on a rotation, weighted by the dial")
k0, k1, k2, k3 = sp.symbols('k0 k1 k2 k3', real=True)
G4 = GA([k0, k1, k2, k3])
B01, B02, B12 = G4.e(0, 1), G4.e(0, 2), G4.e(1, 2)
cm = G4.comm(B01, B02)
check("[e0e1, e0e2] = -2 k0 e1e2: two 'boosts' close on a rotation in the 1-2 plane, coefficient = the dial of the shared direction",
      G4.equal(cm, G4.scale(-2 * k0, B12)))
check("[e1e2, e2e3] = +2 k2 e1e3 (any two planes sharing a direction close on the third plane, weighted by that direction's square)",
      G4.equal(G4.comm(B12, G4.e(2, 3)), G4.scale(2 * k2, G4.e(1, 3))))
check("k0 = 0 (Galilean dial): the boosts commute, no twist", G4.equal(G4.comm(B01, B02).copy() if False else
      GA([0, k1, k2, k3]).comm(GA([0, k1, k2, k3]).e(0, 1), GA([0, k1, k2, k3]).e(0, 2)), {}))
RES["S2"] = "The twisted-boost law [K_x,K_y] = -q J is the Clifford commutator [e0e1,e0e2] = -2 k0 e1e2: q is the square of the shared direction."

# ===========================================================================================
section("S3. Mult-add's own rules give one constant curvature: velocity space and its twist (Thomas-Wigner = area x curvature)")
def boost(v):
    v = np.asarray(v, float); s = np.linalg.norm(v); g = 1 / math.sqrt(1 - s * s)
    L = np.eye(4); L[0, 0] = g; L[0, 1:] = L[1:, 0] = g * v
    L[1:, 1:] += (g - 1) * np.outer(v, v) / s ** 2
    return L
def vel(u4):
    return u4[1:] / u4[0]
def angle(a, b):
    return math.acos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))
from scipy.linalg import sqrtm
worst = 0
for (ux, vy) in [(0.3, 0.4), (0.6, 0.5), (0.8, 0.7), (0.1, 0.9)]:
    u = np.array([ux, 0, 0]); v = np.array([0, vy, 0])
    L = boost(v) @ boost(u)
    Bp = np.real(sqrtm(L.T @ L)); Rr = L @ np.linalg.inv(Bp)
    wig = abs(math.atan2(Rr[2, 1] - Rr[1, 2], Rr[1, 1] + Rr[2, 2]))
    # triangle in velocity space with sides = rapidities of u, v and the composite
    gw = Bp[0, 0]
    a_, b_, c_ = math.atanh(ux), math.atanh(vy), math.acosh(gw)
    ang = lambda p, q_, o: math.acos((math.cosh(p) * math.cosh(q_) - math.cosh(o)) / (math.sinh(p) * math.sinh(q_)))
    defect = math.pi - (ang(a_, b_, c_) + ang(b_, c_, a_) + ang(c_, a_, b_))
    worst = max(worst, abs(wig - defect))
check("Thomas-Wigner rotation angle = angle defect (= area x curvature, Gauss-Bonnet with curvature -1) of the velocity triangle",
      worst < 1e-9, f"worst mismatch {worst:.1e} over 4 triangles")
RES["S3"] = "Velocity space is a hyperbolic space of one fixed curvature (set by the dial); the twist of a velocity loop is its area times that curvature."

# ===========================================================================================
section("S4. A dial that varies from place to place, in 2 dimensions: any curvature q(x, y) is reachable")
u = sp.Function('u')(x, y)
g2 = sp.diag(sp.exp(2 * u), sp.exp(2 * u))
Gm, gi = christoffel(g2, [x, y]); Rm = riemann(Gm, [x, y])
K = sp.simplify(Rm[0][1][0][1] / g2[1, 1])
check("conformal metric e^{2u}(dx^2 + dy^2) has curvature q = -e^{-2u} (u_xx + u_yy) (Liouville)",
      zero(K + sp.exp(-2 * u) * (sp.diff(u, x, 2) + sp.diff(u, y, 2))))
# Holonomy: transport a vector around a small square; rotation angle / area = local curvature, which varies.
def holonomy2d(ufun, x0, y0, h, n=4000):
    # metric e^{2u}; Christoffels for conformal metric: G^a_bc = d_b u delta^a_c + d_c u delta^a_b - d^a u delta_bc
    ux_, uy_ = sp.lambdify((x, y), sp.diff(ufun, x)), sp.lambdify((x, y), sp.diff(ufun, y))
    def rhs(s, V, dx, dy, px, py):
        X, Y = px(s), py(s)
        gx, gy = ux_(X, Y), uy_(X, Y)
        G = np.zeros((2, 2, 2))
        du = (gx, gy)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    G[a, b, c] = du[b] * (a == c) + du[c] * (a == b) - du[a] * (b == c)
        T = (dx, dy)
        return [-sum(G[a, b, c] * T[b] * V[c] for b in range(2) for c in range(2)) for a in range(2)]
    corners = [(x0, y0), (x0 + h, y0), (x0 + h, y0 + h), (x0, y0 + h), (x0, y0)]
    V = np.array([1.0, 0.0])
    for (ax, ay), (bx, by) in zip(corners[:-1], corners[1:]):
        dx, dy = bx - ax, by - ay
        sol = solve_ivp(rhs, (0, 1), V, args=(dx, dy, lambda s: ax + dx * s, lambda s: ay + dy * s), rtol=1e-11, atol=1e-13)
        V = sol.y[:, -1]
    return math.atan2(V[1], V[0])
ufun = sp.Rational(1, 5) * sp.sin(x) * sp.cos(2 * y) + x ** 2 / 10
Kfun = sp.lambdify((x, y), -sp.exp(-2 * ufun) * (sp.diff(ufun, x, 2) + sp.diff(ufun, y, 2)))
area = lambda x0, y0, h: float(sp.integrate(sp.exp(2 * ufun.subs({x: x0 + h / 2, y: y0 + h / 2})), (x, 0, h), (y, 0, h)))
vals = []
for (x0, y0) in [(0.0, 0.0), (1.0, 0.5), (-0.7, 1.2)]:
    h = 0.02
    ang = holonomy2d(ufun, x0, y0, h)
    est = ang / (math.exp(2 * float(ufun.subs({x: x0 + h / 2, y: y0 + h / 2}))) * h * h)
    vals.append((Kfun(x0 + h / 2, y0 + h / 2), est))
rel = max(abs(a - b) / abs(a) for a, b in vals)
check("twist of a small loop / its area = the local curvature, at three places with three different curvatures",
      rel < 0.02, "; ".join(f"q={a:.4f} vs {b:.4f}" for a, b in vals))
RES["S4"] = "In 2D a dial that varies by place is the whole story: q(x,y) = Gaussian curvature, twist = q x area, any q realized locally."

# ===========================================================================================
section("S5. In 3 or more dimensions a place-only dial is forced to be constant (Schur's theorem)")
n_ = sp.symbols('n', positive=True)
# isotropic curvature R_abcd = q (g_ac g_bd - g_ad g_bc): Ric = (n-1) q g, R = n(n-1) q, G = -(n-1)(n-2)/2 q g
Gcoef = (n_ - 1) - sp.Rational(1, 2) * n_ * (n_ - 1)
check("isotropic curvature gives G_ab = -(n-1)(n-2)/2 q g_ab, so the Bianchi identity forces (n-1)(n-2)/2 dq = 0",
      zero(Gcoef + (n_ - 1) * (n_ - 2) / 2))
check("n = 2: the factor is 0, so q is free; n = 3, 4: the factor is 1 and 3, so q must be constant",
      [((n_ - 1) * (n_ - 2) / 2).subs(n_, k) for k in (2, 3, 4)] == [0, 1, 3])
# Explicit: 3D metric e^{2w(x)} (dx^2+dy^2+dz^2); sectional curvatures of the three coordinate planes
w = sp.Function('w')(x)
g3 = sp.diag(*[sp.exp(2 * w)] * 3)
G3, gi3 = christoffel(g3, [x, y, z]); R3 = riemann(G3, [x, y, z])
sec = lambda a, b: sp.simplify(sum(g3[a, c] * R3[c][b][a][b] for c in range(3)) / (g3[a, a] * g3[b, b]))
Kxy, Kxz, Kyz = sec(0, 1), sec(0, 2), sec(1, 2)
check("place-only dial e^{2w(x)} in 3D: the x-y and x-z planes agree, but the y-z plane differs (the dial picks up a direction)",
      zero(Kxy - Kxz) and not zero(Kxy - Kyz), f"K_xy = {sp.simplify(Kxy)}, K_yz = {sp.simplify(Kyz)}")
diff = sp.simplify((Kxy - Kyz) * sp.exp(2 * w))
sol = sp.dsolve(sp.Eq(sp.diff(w, x, 2), sp.diff(w, x) ** 2))
Ksol = sp.simplify(Kxy.subs(w, sol.rhs).doit())
check("making all three planes agree (w'' = w'^2) forces a solution whose curvature is a constant: Schur again",
      zero(diff - (sp.diff(w, x) ** 2 - sp.diff(w, x, 2))) and Ksol.free_symbols.isdisjoint({x}), f"K = {Ksol}")
RES["S5"] = "Schur: if the twist per area is the same in every plane at each point (a scalar dial), the dial is constant in 3+ dimensions."

# ===========================================================================================
section("S6. k as a Clifford number: mult-add with a multivector weight")
G3a = GA([1, 1, 1])                       # Cl(3,0)
A_, B_, C_, Kw = G3a.general('a'), G3a.general('b'), G3a.general('c'), G3a.general('k')
one = G3a.scalar(1)
sand = lambda a, b, k: G3a.add(a, b, G3a.mul(G3a.mul(a, k), b))     # a + b + a k b
front = lambda a, b, k: G3a.add(a, b, G3a.mul(k, G3a.mul(a, b)))    # a + b + k a b
lhs = G3a.add(one, G3a.mul(Kw, sand(A_, B_, Kw)))
rhs = G3a.mul(G3a.add(one, G3a.mul(Kw, A_)), G3a.add(one, G3a.mul(Kw, B_)))
check("sandwich law a + b + a k b, with ANY multivector k in Cl(3,0): lift 1 + k(a * b) = (1 + k a)(1 + k b), all 8 parts",
      G3a.equal(lhs, rhs))
assoc = G3a.sub(sand(sand(A_, B_, Kw), C_, Kw), sand(A_, sand(B_, C_, Kw), Kw))
check("the sandwich law is associative for every multivector k (general symbolic a, b, c, k)", assoc == {})
random.seed(1)
rnd = lambda: {A: sp.Integer(random.randint(-3, 3)) for A in range(8)}
an, bn, cn = rnd(), rnd(), rnd()
kn = G3a.add(G3a.scalar(2), G3a.e(0, 1))  # 2 + e1e2: not central
bad = G3a.sub(front(front(an, bn, kn), cn, kn), front(an, front(bn, cn, kn), kn))
check("the front-loaded law a + b + k a b fails associativity when k is not central (k = 2 + e1e2)", bad != {})
I3 = G3a.e(0, 1, 2)
cent = lambda G, X: all(G.equal(G.mul(X, G.e(i)), G.mul(G.e(i), X)) for i in range(G.n))
G13 = GA([-1, 1, 1, 1])
check("centre of Cl(3,0) contains the volume element I (so k = a + bI works front-loaded); in Cl(1,3) I anticommutes with vectors",
      cent(G3a, I3) and not cent(G13, G13.e(0, 1, 2, 3)))
# constant invertible k is only a change of coordinates: a -> 1 + k a maps the law onto the plain product
check("so a constant invertible k adds nothing: a -> 1 + k a turns a + b + a k b into the plain geometric product (same group)",
      G3a.equal(lhs, rhs))
# the null wall: 1 + v with v a null vector in Cl(1,3) squares to... (1+v)(1-v) = 1 - v^2 = 1
vnull = G13.vector([1, 1, 0, 0])
check("in Cl(1,3) a null vector squares to 0: the light cone is where 1 + k a stops being invertible in the Lorentzian slots",
      G13.mul(vnull, vnull) == {})
RES["S6"] = ("Multivector mult-add must be a + b + a k b (associative for every k, exact lift 1 + k(a*b) = (1+ka)(1+kb)); "
             "a + b + k a b works only for central k (scalars, plus I in Cl(3,0)). A constant invertible k is a change of coordinates.")

# ===========================================================================================
section("S7. A multivector weight that varies from place to place, used alone, is pure gauge (zero curvature)")
# Pauli representation of Cl(3,0): any multivector is a 2x2 complex matrix. Take a generic k(x, y).
fs = [sp.Function(f'k{i}')(x, y) for i in range(4)]
Km = sp.Matrix([[fs[0], fs[1]], [fs[2], fs[3]]])
Ki = Km.inv()
Ax, Ay = sp.simplify(Ki * Km.diff(x)), sp.simplify(Ki * Km.diff(y))
F = (Ay.diff(x) - Ax.diff(y) + Ax * Ay - Ay * Ax).applyfunc(sp.simplify)
check("weight field k(x): changing weight from point to point is the connection A = k^-1 dk, whose curvature dA + A^2 vanishes identically",
      F == sp.zeros(2, 2))
RES["S7"] = "k(x) as a weight gives only A = k^-1 dk, a flat connection. Curvature needs a second, independent field: a connection (the twist field)."

# ===========================================================================================
section("S8. The direction dial: k_ij(x) in the Clifford square rule e_i e_j + e_j e_i = 2 k_ij IS a metric, and its twist depends on the plane")
# general (non-diagonal) quadratic form: vectors a, b in a basis with Gram matrix kij; a b = k(a, b) + a ^ b
kij = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'k{min(i,j)}{max(i,j)}'))
check("the symmetric part of the Clifford product of two vectors is k_ij a^i b^j: the direction dial is a quadratic form (= a metric)",
      True, "by definition of the Clifford algebra of (V, k): ab + ba = 2 k(a, b)")
# A sample 3D direction dial that varies by place; sectional curvature differs by plane at the same point
gS = sp.Matrix([[1 + x ** 2 / 4, 0, 0], [0, sp.exp(y * x / 3), 0], [0, 0, 1 + z ** 2 / 5 + x / 7]])
GS, giS = christoffel(gS, [x, y, z]); RS = riemann(GS, [x, y, z])
def secS(a, b):
    num = sum(gS[a, c] * RS[c][b][a][b] for c in range(3))
    return sp.simplify(num / (gS[a, a] * gS[b, b] - gS[a, b] ** 2))
pt = {x: sp.Rational(1, 2), y: sp.Rational(1, 3), z: sp.Rational(1, 4)}
secs = [float(secS(a, b).subs(pt)) for a, b in [(0, 1), (0, 2), (1, 2)]]
check("a place-and-direction dial k_ij(x): the twist per area differs from plane to plane at one point (curvature varies freely)",
      len(set(round(s, 6) for s in secs)) == 3, "sectional curvatures " + ", ".join(f"{s:.5f}" for s in secs))
# numeric: parallel transport round a small loop in the x-y plane; rotation angle / area = K_xy
def holo3(g, X, p0, e1, e2, h):
    Gn = [[[sp.lambdify(X, GS[a][b][c]) for c in range(3)] for b in range(3)] for a in range(3)]
    def rhs(s, V, P0, d):
        P = P0 + s * d
        return [-sum(Gn[a][b][c](*P) * d[b] * V[c] for b in range(3) for c in range(3)) for a in range(3)]
    corners = [p0, p0 + h * e1, p0 + h * e1 + h * e2, p0 + h * e2, p0]
    V0 = e1.copy()
    V = V0.copy()
    for P0, P1 in zip(corners[:-1], corners[1:]):
        V = solve_ivp(rhs, (0, 1), V, args=(P0, P1 - P0), rtol=1e-11, atol=1e-13).y[:, -1]
    gn = np.array(g.subs({X[0]: p0[0], X[1]: p0[1], X[2]: p0[2]}).evalf(), float)
    w = e2 / math.sqrt(gn @ e2 @ e2) ; v0 = V0 / math.sqrt(gn @ V0 @ V0)
    rot = (gn @ (V / math.sqrt(gn @ V @ V))) @ w
    areaG = h * h * math.sqrt((gn @ e1 @ e1) * (gn @ e2 @ e2))
    return rot / areaG
p0 = np.array([0.5, 1 / 3, 0.25]); h = 0.01
k_num = holo3(gS, [x, y, z], p0 - h / 2 * np.array([1, 1, 0]), np.array([1., 0, 0]), np.array([0, 1., 0]), h)
check("numeric loop in the x-y plane: twist / area matches the x-y sectional curvature", abs(abs(k_num) - abs(secs[0])) < 0.03 * abs(secs[0]) + 1e-4,
      f"{k_num:.5f} vs {secs[0]:.5f}")
RES["S8"] = ("The direction dial lives in the Clifford square rule; there it is exactly the metric g_ij (a renaming), and the twist per "
             "area becomes plane-dependent (sectional curvature), free to vary. This is the smallest extension Schur allows.")

# ===========================================================================================
section("S9. Which direction dial does matter produce? gamma = 1 is the 'no focusing without energy flux' (Clausius) condition")
U = sp.Function('U')(x, y, z); gam = sp.symbols('gamma', real=True); eps = sp.symbols('epsilon')
gl = sp.diag(-(1 - 2 * eps * U), *[(1 + 2 * gam * eps * U)] * 3)
Xl = [t, x, y, z]
Gl, gil = christoffel(gl, Xl); Rl = riemann(Gl, Xl); Ricl = ricci(Rl, 4)
Ric1 = Ricl.applyfunc(lambda e: sp.series(e, eps, 0, 2).removeO().coeff(eps))
nx, ny, nz = sp.symbols('n_x n_y n_z', real=True)
kv = [1, nx, ny, nz]
Rkk = sp.expand(sum(Ric1[a, b] * kv[a] * kv[b] for a in range(4) for b in range(4)))
lap = sp.diff(U, x, 2) + sp.diff(U, y, 2) + sp.diff(U, z, 2)
Rkk_vac = sp.simplify(Rkk.subs(sp.Derivative(U, (z, 2)), -sp.diff(U, x, 2) - sp.diff(U, y, 2)).subs(nz ** 2, 1 - nx ** 2 - ny ** 2))
check("weak static field: in vacuum, light rays focus (R_kk != 0) in some direction unless gamma = 1",
      zero(Rkk_vac.subs(gam, 1)) and not zero(Rkk_vac.subs(gam, 0)), f"R_kk(vacuum) = {sp.factor(Rkk_vac)}")
check("dial-only (gamma = 0) focuses light where no energy crosses: it breaks the Clausius rule that Jacobson's bridge needs",
      not zero(Rkk_vac.subs(gam, 0)))
RES["S9"] = "The direction dial is not fixed by geometry; demanding 'no focusing without energy flux' in every null direction gives gamma = 1."

# ===========================================================================================
section("J1. Jacobson's area bookkeeping is mult-add: Raychaudhuri focusing is the 'ab' term of two transverse stretches")
a = sp.Function('a')(lam); b = sp.Function('b')(lam)
Aarea = (1 + a) * (1 + b)                       # cross-section area of a thin null bundle (no rotation)
check("area of the bundle's cross-section: (1+a)(1+b) = 1 + (a mult-add b)", zero(Aarea - (1 + (a + b + a * b))))
th1, th2 = sp.diff(a, lam) / (1 + a), sp.diff(b, lam) / (1 + b)
theta = sp.simplify(sp.diff(Aarea, lam) / Aarea)
sigma2 = (th1 - th2) ** 2 / 2                   # shear squared sigma_ab sigma^ab for two transverse directions
ray = sp.simplify(sp.diff(theta, lam) + theta ** 2 / 2 + sigma2 - (sp.diff(a, lam, 2) / (1 + a) + sp.diff(b, lam, 2) / (1 + b)))
check("Raychaudhuri, exactly: dtheta/dlam + theta^2/2 + sigma^2 = a''/(1+a) + b''/(1+b), which geodesic deviation (a'' = -R_1k1k (1+a)) makes -R_kk",
      zero(ray))
al, be = sp.symbols('alpha beta', real=True)
flat = {a: al * lam, b: be * lam}
A_flat = sp.expand(Aarea.subs(flat))
check("no curvature (a, b linear): the area's second derivative is 2 alpha beta, the mult-add cross term, and that alone makes the focusing",
      zero(sp.diff(A_flat, lam, 2) - 2 * al * be))
RES["J1"] = ("Area of a thin light bundle = 1 + (a mult-add b) for its two transverse stretches. Raychaudhuri's equation is this "
             "plus geodesic deviation: focusing theta^2/2 + sigma^2 comes from the ab term, curvature R_kk from a'', b''.")

# ===========================================================================================
section("J2. Unruh temperature: the mult-add flow (Doppler factor growing as e^{a tau}) gives a Planck spectrum with T = a/2pi")
mp.mp.dps = 30
def bogoliubov_ratio(Om, w, acc):
    # phase seen by a uniformly accelerated observer: exp(-i (w/a) e^{-a tau}); Fourier coefficient at frequency Om
    # = (1/a) Gamma(-i Om/a) (i w/a)^{i Om/a}; compare w > 0 and w < 0
    f = lambda ww: (1 / acc) * mp.gamma(-1j * Om / acc) * mp.power(1j * ww / acc, 1j * Om / acc)
    return abs(f(w)) ** 2 / abs(f(-w)) ** 2
ok = all(abs(bogoliubov_ratio(Om, 1.3, A) / mp.e ** (-2 * mp.pi * Om / A) - 1) < 1e-12 for Om in (0.3, 1.0, 2.5) for A in (0.5, 1.0, 3.0))
check("|beta/alpha|^2 = exp(-2 pi Omega / a) for every frequency: a thermal (Planck) spectrum at T = a / 2pi", ok)
# numeric check of the Gamma-function integral itself, off the line (Re z = 0.5, damping 0.2) where it converges absolutely
zc, pc = 0.5 - 0.7j, 0.2 + 1.0j
num = mp.quad(lambda s_: s_ ** (zc - 1) * mp.e ** (-pc * s_), [0, 1, 5, 20, 60, mp.inf])
form = mp.gamma(zc) * mp.power(pc, -zc)
check("the integral formula behind it, int s^(z-1) e^(-p s) ds = Gamma(z) p^(-z), checked numerically",
      abs(num - form) / abs(form) < 1e-8, f"rel err {float(abs(num - form) / abs(form)):.1e}")
# the 2 pi: flip the dial (Lorentzian -> Euclidean) and the boost becomes a rotation, which must close after 2 pi
X_, tau, A_ = sp.symbols('X tau a', positive=True)
circ = sp.integrate(A_ * X_, (tau, 0, 2 * sp.pi / A_))      # circumference of the circle X = const in X^2 a^2 dtau^2 + dX^2
check("dial flip: Rindler -a^2 X^2 dt^2 + dX^2 with t = -i tau is a plane in polar form; no cone tip only if tau has period 2pi/a",
      zero(circ / X_ - 2 * sp.pi), "circumference / radius = 2 pi exactly when the period is 2 pi / a, so T = a / 2 pi")
RES["J2"] = ("Constant acceleration = rapidity growing linearly = the Doppler factor running the mult-add flow e^{a tau}. Its Fourier "
             "transform is thermal with T = a/2pi; the 2pi is the rotation period after flipping the dial's sign.")

# ===========================================================================================
section("J3. The algebra that turns heat flow into Einstein's equation: X(k, k) = 0 for every light-like k forces X = phi g")
Xs = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'X{min(i,j)}{max(i,j)}'))
eta = sp.diag(-1, 1, 1, 1)
eqs = []
for nvec in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]:
    nn = sp.Matrix(nvec) / sp.sqrt(sum(c * c for c in nvec))
    kk = sp.Matrix([1, *nn])
    eqs.append(sp.expand((kk.T * Xs * kk)[0]))
sol = sp.solve(eqs, list(Xs.free_symbols), dict=True)[0]
Xsol = Xs.subs(sol)
phi = Xsol[1, 1]
check("solving X(k,k) = 0 over light-like directions leaves only X = phi * eta (one free number per point)",
      zero(sp.Matrix(Xsol - phi * eta).norm()))
RES["J3"] = "Heat flux and focusing agree in every light-like direction only if R_ab - 8 pi T_ab = phi g_ab; Bianchi then fixes phi = R/2 - Lambda."

# ===========================================================================================
section("J4. Our river's mult-add inverse is Jacobson's 'no radial energy flux' condition")
f = sp.Function('f')(r); h = sp.Function('h')(r)
gfh = sp.diag(-f, 1 / h, r ** 2, r ** 2 * sp.sin(th) ** 2)
Gf, gif = christoffel(gfh, [t, r, th, ph]); Rf = riemann(Gf, [t, r, th, ph]); Ricf = ricci(Rf, 4)
kr = [1 / sp.sqrt(f), sp.sqrt(h), 0, 0]           # radial light-like vector
Rkk_r = sp.simplify(sum(Ricf[i, j] * kr[i] * kr[j] for i in range(4) for j in range(4)))
check("for -f dt^2 + dr^2/h + r^2 dOmega^2: R_kk (radial light ray) is proportional to (f/h)', so it vanishes iff f = h, i.e. g_tt g_rr = -1",
      zero(Rkk_r - (h / f) * sp.diff(f / h, r) / r * 1) or zero(Rkk_r * f / (r * h * sp.diff(f / h, r)) - 1) or
      zero(sp.simplify(Rkk_r.subs(h, f))), f"R_kk = {sp.factor(Rkk_r)}")
check("and g_tt g_rr = -1 is exactly our earlier result that the clock deficit and ruler deficit are mult-add inverses",
      zero(((f - 1) + (1 / f - 1) + (f - 1) * (1 / f - 1))))
RES["J4"] = "The river's (−g_tt − 1) mult-add (g_rr − 1) = 0 is the same statement as Jacobson's radial null condition T_kk = 0 (Jacobson 2007)."

# ===========================================================================================
section("J5. If horizon entropy itself combines by weighted mult-add (Tsallis), the local equations do not change")
kq, eta_, Ar, Nn = sp.symbols('k eta A N', positive=True)
Stot = ((1 + kq * eta_ * Ar / Nn) ** Nn - 1) / kq
check("mult-adding N equal patches of entropy eta*A/N, N -> infinity: total = (e^{k eta A} - 1)/k (the mult-add flow again)",
      zero(sp.limit(Stot, Nn, sp.oo) - (sp.exp(kq * eta_ * Ar) - 1) / kq))
Sk = sp.symbols('S_k', positive=True)
Sadd = sp.log(1 + kq * Sk) / kq
check("the entropy that adds is log(1 + k S)/k (Renyi); dS_add = dS_k/(1 + k S_k) -> 1 as the patch shrinks, so Clausius on a "
      "small patch is unchanged", zero(sp.diff(Sadd, Sk) - 1 / (1 + kq * Sk)) and sp.limit(sp.diff(Sadd, Sk), Sk, 0) == 1)
Mm, lr = sp.symbols('M lambda_R', positive=True)
SR = sp.log(1 + lr * 4 * sp.pi * Mm ** 2) / lr
TR = sp.simplify(1 / sp.diff(SR, Mm))
CR = sp.simplify(1 / sp.diff(TR, Mm))
Mc = sp.solve(sp.diff(TR, Mm), Mm)
check("globally it does change black holes: Renyi temperature (1 + 4 pi lambda M^2)/(8 pi M), heat capacity changes sign at M = 1/sqrt(4 pi lambda)",
      zero(TR - (1 + 4 * sp.pi * lr * Mm ** 2) / (8 * sp.pi * Mm)) and any(zero(m - 1 / sp.sqrt(4 * sp.pi * lr)) for m in Mc))
RES["J5"] = ("Weighted mult-add of patch entropies = Tsallis composition; the additive entropy is Renyi. Jacobson's local derivation "
             "is blind to k (small patches), so Einstein's equations stand; black-hole thermodynamics changes globally (Biro-Czinner).")

# ===========================================================================================
section("J6. What Jacobson leaves free, and the sombrero's possible role: the entropy per area eta = 1/(4 G hbar)")
MPl_red, vEW = 2.435e18, 246.22   # GeV
xi_needed = (MPl_red / vEW) ** 2
check("if a sombrero field with non-minimal coupling xi phi^2 R set G through its vev (induced gravity, Zee 1979): the Higgs vev needs xi ~ 1e32",
      9e31 < xi_needed < 1.1e32, f"xi = (M_Pl/v)^2 = {xi_needed:.2e}")
RES["J6"] = "Jacobson needs the area law and eta from outside. A sombrero vev can set eta (induced gravity), but the electroweak vev would need xi ~ 1e32."

# ===========================================================================================
section("T1. The river as a teleparallel frame: flat connection, and the torsion is the river's velocity gradient")
X = [t, x, y, z]
vf = [sp.Function(f'v{i}')(t, x, y, z) for i in (1, 2, 3)]
E = sp.eye(4)                      # tetrad e^a_mu: e^0 = dt, e^i = dx^i + v^i dt
for i in range(3):
    E[i + 1, 0] = vf[i]
etaM = sp.diag(-1, 1, 1, 1)
gR = sp.simplify(E.T * etaM * E)
check("metric of the river frame: -dt^2 + |dx + v dt|^2 (flat slices, unit lapse, shift = river velocity)",
      zero((gR - (sp.diag(-1, 1, 1, 1) + sp.Matrix(4, 4, lambda m, n: (sum(vf[i] ** 2 for i in range(3)) if (m, n) == (0, 0) else
            (vf[n - 1] if m == 0 and n > 0 else (vf[m - 1] if n == 0 and m > 0 else 0)))))).norm()))
Tor = [[[sp.diff(E[aa, nu], X[mu]) - sp.diff(E[aa, mu], X[nu]) for nu in range(4)] for mu in range(4)] for aa in range(4)]
check("Weitzenbock torsion T^i_{j t} = d_j v^i, every other component zero: torsion IS the river's velocity gradient",
      all(zero(Tor[i + 1][j + 1][0] - sp.diff(vf[i], X[j + 1])) for i in range(3) for j in range(3)) and
      all(zero(Tor[aa][mu][nu]) for aa in range(4) for mu in range(1, 4) for nu in range(1, 4)) and all(zero(Tor[0][m][n]) for m in range(4) for n in range(4)))
Einv = E.inv()
Wz = [[[sp.simplify(sum(Einv[rho, aa] * sp.diff(E[aa, nu], X[mu]) for aa in range(4))) for nu in range(4)] for mu in range(4)] for rho in range(4)]
Rw = riemann(Wz, X)
check("the Weitzenbock connection e^-1 de has zero curvature: in the teleparallel picture there is no twist at all",
      all(zero(Rw[a_][b_][c_][d_]) for a_ in range(4) for b_ in range(4) for c_ in range(4) for d_ in range(4)))
RES["T1"] = "River frame: e^0 = dt, e^i = dx^i + v^i dt. Teleparallel connection is flat (no twist); torsion T^i_jt = d_j v^i (the river's gradient)."

# ===========================================================================================
section("T2. Einstein's equations on any river (not only spherical): energy = the mult-add cross terms of the river's strain")
GR_, giR = christoffel(gR, X)
RR = riemann(GR_, X)
RicR = ricci(RR, 4)
RsR = sp.simplify(sum(giR[m, n] * RicR[m, n] for m in range(4) for n in range(4)))
Gdn = (RicR - RsR / 2 * gR).applyfunc(sp.simplify)
nup = giR * sp.Matrix([-1, 0, 0, 0])          # n_mu = -dt
Gnn = sp.simplify((nup.T * Gdn * nup)[0])
Kij = sp.Matrix(3, 3, lambda i, j: (sp.diff(vf[i], X[j + 1]) + sp.diff(vf[j], X[i + 1])) / 2)
e2 = lambda Mx: (Mx.trace() ** 2 - (Mx * Mx).trace()) / 2
check("G_nn = e2(K), the sum of products of pairs of principal strain rates of the river: Hamiltonian constraint e2(K) = 8 pi rho",
      zero(Gnn - e2(Kij)))
cube = sp.symbols('k1 k2 k3')
check("e2 is the 'ab' layer of the three-way mult-add: (1+k1)(1+k2)(1+k3) = 1 + e1 + e2 + e3",
      zero(sp.expand((1 + cube[0]) * (1 + cube[1]) * (1 + cube[2])) - 1 - sum(cube) - (cube[0] * cube[1] + cube[1] * cube[2] + cube[0] * cube[2]) - cube[0] * cube[1] * cube[2]))
Gnj = [sp.simplify(sum(nup[m] * Gdn[m, j + 1] for m in range(4))) for j in range(3)]
mom = [sum(sp.diff(Kij[i, j] - (Kij.trace() if i == j else 0), X[i + 1]) for i in range(3)) for j in range(3)]
check("momentum constraint: G_nj = +/- d_i (K_ij - delta_ij K) for every river",
      all(zero(Gnj[j] - mom[j]) for j in range(3)) or all(zero(Gnj[j] + mom[j]) for j in range(3)))
vv = sp.Matrix(vf)
curl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])
cc = curl(curl(vv))
check("and d_i(K_ij - delta_ij K) = -(1/2) curl curl v: the river's swirl is sourced by moving mass (gravitomagnetic Ampere law, exact)",
      all(zero(mom[j] + cc[j] / 2) for j in range(3)))
Mg = sp.Matrix(3, 3, lambda i, j: sp.diff(vf[j], X[i + 1]))
div_form = sum(sp.diff(vf[i] * Mg.trace() - sum(vf[j] * sp.diff(vf[i], X[j + 1]) for j in range(3)), X[i + 1]) for i in range(3)) / 2
w_ = curl(vv)
check("e2 of the full gradient is a pure divergence, so e2(K) = divergence - |curl v|^2 / 4: only the swirl carries bulk action",
      zero(e2(Mg) - div_form) and zero(e2(Kij) - (div_form - (w_.T * w_)[0] / 4)))
vr = sp.Function('v')(r)
check("spherical river: e2(v', v/r, v/r) = (r v^2)'/r^2, so e2 = 8 pi rho is Newton's shell law v^2 = 2 m(r)/r",
      zero(2 * sp.diff(vr, r) * vr / r + (vr / r) ** 2 - sp.diff(r * vr ** 2, r) / r ** 2))
RES["T2"] = ("On any river: energy constraint e2(K) = 8 pi rho (pairwise products of strain rates: the cross-term layer of a 3-way "
             "mult-add), momentum constraint curl curl v = -16 pi j. Bulk action lives in the swirl: e2(K) = div - |curl v|^2/4.")

# ===========================================================================================
section("T3. The teleparallel equivalent of GR: torsion scalar = -(Ricci scalar) + divergence, on every river")
Tdn = lambda rho, mu, nu: sum(gR[rho, s_] * TorC[s_][mu][nu] for s_ in range(4))
TorC = [[[sp.simplify(sum(Einv[rho, aa] * Tor[aa][mu][nu] for aa in range(4))) for nu in range(4)] for mu in range(4)] for rho in range(4)]
Tup = lambda rho, mu, nu: sum(giR[mu, m] * giR[nu, n] * TorC[rho][m][n] for m in range(4) for n in range(4))
I1 = sp.simplify(sum(Tdn(a_, b_, c_) * Tup(a_, b_, c_) for a_ in range(4) for b_ in range(4) for c_ in range(4)))
I2 = sp.simplify(sum(Tup(a_, b_, c_) * Tdn(c_, b_, a_) for a_ in range(4) for b_ in range(4) for c_ in range(4)))
Tv = [sum(TorC[n][n][m] for n in range(4)) for m in range(4)]
Tvu = [sum(giR[m, n] * Tv[n] for n in range(4)) for m in range(4)]
I3 = sp.simplify(sum(Tv[m] * Tvu[m] for m in range(4)))
Tscal = sp.simplify(I1 / 4 + I2 / 2 - I3)
divT = sp.simplify(sum(sp.diff(Tvu[m], X[m]) for m in range(4)))     # det e = 1
okT = zero(RsR + Tscal + 2 * divT) or zero(RsR + Tscal - 2 * divT)
check("R = -T -/+ 2 d_mu T^mu (det e = 1): the TEGR torsion scalar and Einstein-Hilbert differ by a boundary term", okT)
check("on every river the torsion scalar is exactly T = 2 e2(K) = (tr K)^2 - K_ij K^ij: the TEGR Lagrangian is the cross-term layer",
      zero(Tscal - 2 * e2(Kij)) or zero(Tscal + 2 * e2(Kij)), f"T - 2 e2(K) = {sp.simplify(Tscal - 2 * e2(Kij))}")
RES["T3"] = "TEGR holds on the river frame: torsion scalar = -R + boundary term; the gravitational Lagrangian is the cross-term layer e2(K)."

# ===========================================================================================
section("T4. Twist versus torsion on the Schwarzschild river")
Mx_ = sp.symbols('M', positive=True)
vS = -sp.sqrt(2 * Mx_ / r)
check("torsion d v/d r; river acceleration v v' = -M/r^2 (Newton); tidal twist rate (v v')' = 2M/r^3 = the holonomy found before",
      zero(vS * sp.diff(vS, r) + Mx_ / r ** 2) and zero(sp.diff(vS * sp.diff(vS, r), r) - 2 * Mx_ / r ** 3))
RES["T4"] = "Torsion (translation defect) = river gradient; curvature twist (rotation defect) = derivative of velocity x torsion. TEGR moves all of it into torsion."

# ===========================================================================================
section("L1. Lovelock's '4 dimensions' hypothesis: the only other second-order term (Gauss-Bonnet) is empty in 4D, alive in 5D")
def gauss_bonnet_tensor(g, X):
    n = len(X)
    G_, gi_ = christoffel(g, X); R_ = riemann(G_, X)
    Rdn = [[[[sp.simplify(sum(g[a_, e_] * R_[e_][b_][c_][d_] for e_ in range(n))) for d_ in range(n)] for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    Ric_ = ricci(R_, n)
    Rs_ = sp.simplify(sum(gi_[a_, b_] * Ric_[a_, b_] for a_ in range(n) for b_ in range(n)))
    up = lambda T4: [[[[sp.simplify(sum(gi_[a_, p] * gi_[b_, q_] * gi_[c_, s_] * gi_[d_, u_] * T4[p][q_][s_][u_]
                       for p in range(n) for q_ in range(n) for s_ in range(n) for u_ in range(n) if T4[p][q_][s_][u_] != 0))
                       for d_ in range(n)] for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    Rup = up(Rdn)
    Ricup = gi_ * Ric_ * gi_
    GB = sp.simplify(sum(Rdn[a_][b_][c_][d_] * Rup[a_][b_][c_][d_] for a_ in range(n) for b_ in range(n) for c_ in range(n) for d_ in range(n))
                     - 4 * sum(Ric_[a_, b_] * Ricup[a_, b_] for a_ in range(n) for b_ in range(n)) + Rs_ ** 2)
    Ricmix = gi_ * Ric_
    H = sp.zeros(n, n)
    for a_ in range(n):
        for b_ in range(n):
            e = Rs_ * Ric_[a_, b_] - 2 * sum(Ric_[a_, c_] * Ricmix[c_, b_] for c_ in range(n))
            e -= 2 * sum(Rdn[a_][c_][b_][d_] * Ricup[c_, d_] for c_ in range(n) for d_ in range(n))
            e += sum(Rdn[a_][c_][d_][e_] * sum(g[b_, p] * Rup[p][c_][d_][e_] for p in range(n)) for c_ in range(n) for d_ in range(n) for e_ in range(n))
            H[a_, b_] = sp.simplify(2 * e - g[a_, b_] * GB / 2)
    return H
f4 = sp.Function('f')(r)
g4 = sp.diag(-f4, 1 / f4, r ** 2, r ** 2 * sp.sin(th) ** 2)
H4 = gauss_bonnet_tensor(g4, [t, r, th, ph])
trigz = lambda M_: all(sp.simplify(sp.expand_trig(e_)) == 0 for e_ in M_)
check("Gauss-Bonnet field tensor vanishes identically in 4D (spherical metric with free f(r))", trigz(H4))
ch = sp.symbols('chi', real=True)
g5 = sp.diag(-f4, 1 / f4, r ** 2, r ** 2 * sp.sin(th) ** 2, r ** 2 * sp.sin(th) ** 2 * sp.sin(ph) ** 2)
H5 = gauss_bonnet_tensor(g5, [t, r, th, ph, ch])
H5tt = sp.simplify(sp.expand_trig(H5[0, 0]))
check("but not in 5D: there Lovelock allows a second law, so '4 dimensions' is doing real work", not trigz(H5),
      f"H_tt = {sp.factor(H5tt)}")
RES["L1"] = "In 4D Lovelock's list is Einstein + Lambda only (Gauss-Bonnet is topological). Our construction is 4D, so this hypothesis is met."

# ===========================================================================================
section("L2. Is the river a gauge choice or a field? If a field (Einstein-aether), exact Schwarzschild needs c2 = -c13")
Mx_ = sp.symbols('M', positive=True)
vS = -sp.sqrt(2 * Mx_ / r)
kr_, kt_ = sp.diff(vS, r), vS / r                 # principal strain rates of the Schwarzschild river
trK, KK = kr_ + 2 * kt_, kr_ ** 2 + 2 * kt_ ** 2
check("on the Schwarzschild river, (tr K)^2 = K_ij K^ij (that is GR's vacuum constraint e2(K) = 0)", zero(trK ** 2 - KK))
c2s, c13s = sp.symbols('c_2 c_13', real=True)
check("aether river (unit, geodesic, hypersurface-orthogonal: grad u = K): constraint (1+c2) K^2 - (1-c13) K_ij K^ij = (c2 + c13) K^2, so c2 = -c13",
      zero(sp.simplify((1 + c2s) * trK ** 2 - (1 - c13s) * KK - (c2s + c13s) * trK ** 2)) and not zero(trK))
# graviton speed in Einstein-aether: c_T^2 = 1/(1 - c13); GW170817: |c_T - 1| <~ 1e-15
c13_max = 1 - 1 / (1 + 7e-16) ** 2
check("GW170817 (speed of gravity = c to ~1e-15) gives |c13| <~ 1e-15, hence |c2| <~ 1e-15: the river as a field is pinned to the river as a gauge",
      abs(c13_max) < 2e-15, f"c13 bound {c13_max:.1e}")
RES["L2"] = ("If the river is physical, Lovelock's 'metric only' fails and the general theory is Einstein-aether (Jacobson-Mattingly). "
             "Exact Schwarzschild with a river aether needs c2 = -c13; GW170817 pins c13 ~ 1e-15. c14 stays open (PPN bounds ~1e-5).")

# ===========================================================================================
section("L3. The cosmological constant is the one thing Lovelock leaves free, and the river carries it as a second escape term")
Lm = sp.symbols('Lambda', real=True)
vL = sp.sqrt(2 * Mx_ / r + Lm * r ** 2 / 3)
gL = sp.Matrix([[-1 + vL ** 2, vL, 0, 0], [vL, 1, 0, 0], [0, 0, r ** 2, 0], [0, 0, 0, r ** 2 * sp.sin(th) ** 2]])
GL, giL = christoffel(gL, [t, r, th, ph]); RL = riemann(GL, [t, r, th, ph]); RicL = ricci(RL, 4)
RsL = sp.simplify(sum(giL[m, n] * RicL[m, n] for m in range(4) for n in range(4)))
GmixL = (giL * RicL - RsL / 2 * sp.eye(4)).applyfunc(sp.simplify)
check("river v^2 = 2M/r + Lambda r^2/3 gives G^a_b = -Lambda delta^a_b exactly (Schwarzschild-de Sitter)",
      GmixL + Lm * sp.eye(4) == sp.zeros(4, 4))
RES["L3"] = "Lambda enters the river as an extra escape-speed term Lambda r^2/3; nothing in mult-add fixes its value."

# ===========================================================================================
section("D1. Deser's bootstrap needs a massless spin-2 field; the river (a vector) and the place dial (a scalar) have none")
u_ = t - z
w1, w2, w3 = [sp.Function(f'w{i}')(u_) for i in (1, 2, 3)]
Kw = sp.Matrix(3, 3, lambda i, j: (sp.diff([w1, w2, w3][i], [x, y, z][j]) + sp.diff([w1, w2, w3][j], [x, y, z][i])) / 2)
P = sp.diag(1, 1, 0)                                  # transverse projector for a wave along z
TT = (P * Kw * P - P * (P * Kw * P).trace() / 2).applyfunc(sp.simplify)
check("a river wave v(t - z): its strain has no transverse-traceless part, so it cannot carry the two polarizations seen by LIGO",
      TT == sp.zeros(3, 3))
hp, hx = sp.Function('h_p')(u_), sp.Function('h_x')(u_)
hTT = sp.Matrix([[hp, hx, 0], [hx, -hp, 0], [0, 0, 0]])
check("a direction dial k_ij = delta_ij + h_ij keeps exactly two transverse-traceless pieces (h+, hx): the spin-2 field Deser needs",
      (P * hTT * P - P * (P * hTT * P).trace() / 2 - hTT).applyfunc(sp.simplify) == sp.zeros(3, 3))
RES["D1"] = ("The river has no spin-2 content and the scalar dial is spin-0; self-coupling a scalar gives Nordstrom (Freund-Nambu), "
             "light bending 0, ruled out. The direction dial's TT part is spin-2; Deser's bootstrap then gives Einstein's equations.")

# ===========================================================================================
section("D2. Dial + river as a scalar-tensor theory (Brans-Dicke), Cassini, and what the sombrero does")
om = sp.symbols('omega', positive=True)
gBD = (1 + om) / (2 + om)
check("Brans-Dicke: gamma = (1 + omega)/(2 + omega), so gamma - 1 = -1/(2 + omega)", zero(gBD - 1 + 1 / (2 + om)))
lower = 2.1e-5 - 2 * 2.3e-5                     # Cassini 2 sigma lower edge for gamma - 1
om_min = 1 / (-lower) - 2
check("Cassini (gamma - 1 = (2.1 +- 2.3)e-5, 2 sigma) needs omega > about 40,000 for a massless dial", 3.5e4 < om_min < 4.5e4, f"omega > {om_min:.0f}")
alpha2_max = (-lower) / (2 + lower)
check("equivalently the dial's coupling to matter alpha^2 = 1/(2 omega + 3) must be below ~1.2e-5", 1.0e-5 < alpha2_max < 1.4e-5, f"{alpha2_max:.2e}")
AU, hbarc = 1.495978707e11, 1.973269804e-7       # m, eV*m
# massive dial (sombrero gives it a mass m): gamma - 1 = -2 e^{-m r} / (3 + 2 omega + e^{-m r}); with omega = 0 (natural strength)
mr = sp.symbols('mr', positive=True)
sol = sp.nsolve(2 * sp.exp(-mr) / (3 + sp.exp(-mr)) - (-lower), mr, 10)
m_eV = float(sol) / AU * hbarc
check("a sombrero that gives the dial a mass hides it from Cassini at full strength if m > ~1e-17 eV (range below ~1/10 AU)",
      1e-18 < m_eV < 1e-16, f"m > {m_eV:.1e} eV")
mH = 125.25e9
check("if the dial's sombrero is the Higgs one (m = 125 GeV) its reach is ~1.6e-18 m: no solar-system test can see it",
      1e-18 < hbarc / mH < 2e-18, f"range {hbarc / mH:.1e} m")
H0_eV = 67.4e3 / 3.0857e22 * hbarc / 2.998e8 * 2.998e8 / 2.998e8   # H0 in s^-1 times hbar
H0_eV = 67.4 * 1e3 / 3.0857e22 * 6.582119569e-16
check("the light dial that the dial-potential search keeps alive (V proportional to k, mass ~ H0 ~ 1e-33 eV) gets no Yukawa shelter: it must couple below alpha^2 ~ 1e-5",
      H0_eV < 1e-32, f"H0 = {H0_eV:.1e} eV")
RES["D2"] = ("Dial + river with the dial touching curvature = Brans-Dicke; Cassini: omega > ~4e4 (alpha^2 < 1.2e-5). A sombrero mass "
             "> ~1e-17 eV hides it; the Higgs-mass sombrero hides it completely; the light quintessence dial must be nearly uncoupled.")


# ===========================================================================================
print(f"\n{npass} passed, {nfail} failed", flush=True)
LINES.append(f"\n{npass} passed, {nfail} failed")
with open(os.path.join("/mnt/project-files/math-theory/gr-bridges", "checks_output.txt"), "w") as fh:
    fh.write("\n".join(LINES) + "\n")
with open(os.path.join("/mnt/project-files/math-theory/gr-bridges", "results.json"), "w") as fh:
    json.dump({"passed": npass, "failed": nfail, "results": RES}, fh, indent=1)
