"""Electromagnetism through the mult-add lens.

Run: python3 checks.py  (writes checks_output.txt and results.json).  Natural units c = 1 unless stated.
Order follows history: electromagnetism set the dial (A), its fields ride the boosts (B), its gauge phase is a dial (C),
its field is the corner of a Clifford product and Maxwell is one equation (D), where the lens stops (E), and the
bridge to gravity (F).

Mult-add: a * b = a + b + ab, lift 1 + (a * b) = (1 + a)(1 + b).  Weighted: a + b + k ab.  Landscape law with dial q:
(u + v) / (1 - q u v): q = -1/c^2 is Einstein's velocity addition, q = 0 Galilean, q = +1 the tangent (circle) law.
"""
import json, math, sys, os, random, cmath
import sympy as sp
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/mnt/project-files/math-theory/gr-bridges")
from ga import GA

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
    return sp.simplify(sp.expand(e)) == 0


def ma(a, b):
    return a + b + a * b


t, x, y, z = sp.symbols('t x y z', real=True)
X4 = [t, x, y, z]

# =====================================================================================================
section("A. Electromagnetism sets the dial")
c_SI = 299792458.0
mu0 = 1.25663706212e-6          # CODATA 2018
eps0 = 8.8541878128e-12         # CODATA 2018
prod = c_SI ** 2 * mu0 * eps0
check("A1. c^2 eps0 mu0 = 1 (CODATA 2018): the speed in Maxwell's waves is 1/sqrt(eps0 mu0)",
      abs(prod - 1) < 1e-9, f"c^2 eps0 mu0 = {prod:.12f}")
RES["A1_c2_eps0_mu0"] = prod
qdial = -mu0 * eps0
check("A1b. so the dial of velocity addition is q = -1/c^2 = -eps0 mu0", abs(qdial * c_SI ** 2 + 1) < 1e-9,
      f"q = {qdial:.6e} s^2/m^2")
RES["A1b_q_SI"] = qdial

# A2: a dial boost x' = g(x - v t), t' = g(t + q v x).  Which q leaves the wave operator d_t^2 - c^2 d_x^2 alone?
v, q, c, g = sp.symbols('v q c g', real=True)
f = sp.Function('f')
Xp = g * (x - v * t)
Tp = g * (t + q * v * x)
F_ = f(Xp, Tp)
wave = sp.diff(F_, t, 2) - c ** 2 * sp.diff(F_, x, 2)
Xs, Ts = sp.symbols('Xs Ts')
fxx = sp.Derivative(f(Xs, Ts), Xs, 2)
ftt = sp.Derivative(f(Xs, Ts), Ts, 2)
fxt = sp.Derivative(f(Xs, Ts), Xs, Ts)
w2 = wave.subs({Xp: Xs, Tp: Ts})
w2 = sp.expand(w2.doit())
# collect coefficients of the three second derivatives (sympy writes them via Subs; use a generic polynomial approach)
a1, a2, a3 = sp.symbols('a1 a2 a3')
W = sp.expand(sp.diff((a1 * Xp ** 2 / 2 + a2 * Xp * Tp + a3 * Tp ** 2 / 2), t, 2)
              - c ** 2 * sp.diff((a1 * Xp ** 2 / 2 + a2 * Xp * Tp + a3 * Tp ** 2 / 2), x, 2))
cXX, cXT, cTT = [sp.expand(W).coeff(s) for s in (a1, a2, a3)]
sol_q = sp.solve(sp.Eq(cXT, 0), q)
check("A2. the cross term d_X d_T vanishes only for q = -1/c^2 (Galilean q = 0 leaves it)",
      sol_q == [-1 / c ** 2] and not zero(cXT.subs(q, 0)), f"q = {sol_q}; cross term at q=0: {sp.factor(cXT.subs(q, 0))}")
ratio = sp.simplify((cTT / cXX).subs({q: -1 / c ** 2, g: 1 / sp.sqrt(1 - v ** 2 / c ** 2)}))
check("A2b. at q = -1/c^2 the operator keeps its form exactly (same wave speed c in the new frame)",
      zero(ratio + 1 / c ** 2), f"coefficient ratio TT/XX = {ratio}")
RES["A2"] = "Maxwell's wave equation is invariant only under dial boosts with q = -1/c^2"

u1, u2 = sp.symbols('u1 u2', real=True)
law = lambda a, b, qq: (a + b) / (1 - qq * a * b)
check("A3. composing two dial boosts gives the landscape law (u+v)/(1 - q u v); q = -1 is Einstein's",
      zero(law(u1, u2, -1) - (u1 + u2) / (1 + u1 * u2)))
dop = lambda u: 2 * u / (1 - u)
check("A3b. Doppler-squared minus one composes by mult-add: d(u (+) v) = d(u) * d(v)",
      zero(dop(law(u1, u2, -1)) - ma(dop(u1), dop(u2))))

# =====================================================================================================
section("B. The fields ride the boosts: E and B on the light cone")
D, D1, D2 = sp.symbols('D D1 D2', positive=True)
Ex, Ey, Ez, Bx, By, Bz = sp.symbols('E_x E_y E_z B_x B_y B_z', real=True)


def vD(Dd):
    return (Dd ** 2 - 1) / (Dd ** 2 + 1)


def boost_std(fields, vv):
    ex, ey, ez, bx, by, bz = fields
    gg = 1 / sp.sqrt(1 - vv ** 2)
    return [ex, gg * (ey - vv * bz), gg * (ez + vv * by), bx, gg * (by + vv * ez), gg * (bz - vv * ey)]


def boost_D(fields, Dd):
    gg = (Dd ** 2 + 1) / (2 * Dd)
    gv = (Dd ** 2 - 1) / (2 * Dd)
    ex, ey, ez, bx, by, bz = fields
    return [ex, gg * ey - gv * bz, gg * ez + gv * by, bx, gg * by + gv * ez, gg * bz - gv * ey]


FLD = [Ex, Ey, Ez, Bx, By, Bz]
bs = boost_std(FLD, vD(D))
bd = boost_D(FLD, D)
check("B0. textbook boost of E, B along x equals the rational form in the Doppler factor D (gamma = (D^2+1)/2D)",
      all(zero(sp.simplify(a - b)) for a, b in zip(bs, bd)))
check("B1. light-cone parts: (E_y + B_z) -> (E_y + B_z)/D and (E_y - B_z) -> D (E_y - B_z)",
      zero(bd[1] + bd[5] - (Ey + Bz) / D) and zero(bd[1] - bd[5] - D * (Ey - Bz)))
check("B1b. the other pair: (E_z - B_y) -> /D and (E_z + B_y) -> x D",
      zero(bd[2] - bd[4] - (Ez - By) / D) and zero(bd[2] + bd[4] - D * (Ez + By)))
check("B2. the ratio (E_y - B_z)/(E_y + B_z) is multiplied by D^2 = (1+v)/(1-v) = 1 + d(v), mult-add's lift",
      zero(D ** 2 - (1 + vD(D)) / (1 - vD(D))) and zero(D ** 2 - (1 + dop(vD(D)))))
comp = boost_D(boost_D(FLD, D1), D2)
check("B3. two boosts compose as one with D = D1 D2 (all six components)",
      all(zero(a - b) for a, b in zip(comp, boost_D(FLD, D1 * D2))))
check("B3b. and v(D1 D2) is Einstein's sum of v(D1), v(D2); so d composes by mult-add",
      zero(vD(D1 * D2) - law(vD(D1), vD(D2), -1)) and zero(dop(vD(D1 * D2)) - ma(dop(vD(D1)), dop(vD(D2)))))
inv1 = lambda F: F[0] ** 2 + F[1] ** 2 + F[2] ** 2 - F[3] ** 2 - F[4] ** 2 - F[5] ** 2
inv2 = lambda F: F[0] * F[3] + F[1] * F[4] + F[2] * F[5]
check("B4. E^2 - B^2 and E.B are unchanged by every boost (the field's two invariants)",
      zero(inv1(bd) - inv1(FLD)) and zero(inv2(bd) - inv2(FLD)))
check("B4b. each light-cone pair keeps its product: (E_y+B_z)(E_y-B_z) fixed (a squeeze, like T^2 - X^2)",
      zero((bd[1] + bd[5]) * (bd[1] - bd[5]) - (Ey + Bz) * (Ey - Bz)))

# B5: the same boost as a Clifford sandwich F' = R F R~ with the half-velocity rotor
G = GA([1, -1, -1, -1])
sig = [G.e(k, 0) for k in (1, 2, 3)]
I4 = G.e(0, 1, 2, 3)


def Fmv(E, B):
    Ev = G.add(*[G.scale(E[k], sig[k]) for k in range(3)])
    Bv = G.add(*[G.scale(B[k], sig[k]) for k in range(3)])
    return G.add(Ev, G.mul(I4, Bv))


def rev(Xm):
    out = {}
    for A, cf in Xm.items():
        gr = bin(A).count("1")
        out[A] = cf * (-1 if gr * (gr - 1) // 2 % 2 else 1)
    return out


def read_EB(Fm):
    """read E_k (coefficient of sigma_k = e_k e_0) and B_k (from I B) back off a bivector"""
    E = [sp.simplify(-Fm.get((1 << k) | 1, 0)) for k in (1, 2, 3)]  # e_k e_0 = -e_0 e_k (bitmask stored as e_0 e_k)
    # I sigma_1 = e0123 e10 ; compute basis to read B robustly
    B = []
    for k in range(3):
        Ib = G.mul(I4, sig[k])
        (A, cf), = Ib.items()
        B.append(sp.simplify(Fm.get(A, 0) / cf))
    return E, B


Fm = Fmv([Ex, Ey, Ez], [Bx, By, Bz])
Eback, Bback = read_EB(Fm)
check("B5a. reading E and B back off the bivector F = E + I B works", Eback == [Ex, Ey, Ez] and Bback == [Bx, By, Bz])
s = sp.symbols('s', positive=True)  # s = sqrt(D)
ok_rotor = False
for sign in (1, -1):
    R = G.add(G.scalar((s + 1 / s) / 2), G.scale(sign * (s - 1 / s) / 2, sig[0]))
    Fp = G.mul(G.mul(R, Fm), rev(R))
    E2, B2 = read_EB(Fp)
    target = boost_D(FLD, s ** 2)
    if all(zero(a - b) for a, b in zip(E2 + B2, target)):
        ok_rotor = sign
check("B5. the boost is the sandwich F' = R F R~ with rotor R = cosh(phi/2) +/- sinh(phi/2) sigma_x (D = e^phi): "
      "half-velocities, as on the Clifford page", ok_rotor != False, f"sign = {ok_rotor}")
RES["B5_rotor_sign"] = ok_rotor

# B6: a plane wave is a null field on one light-cone axis; boosting rescales it by the Doppler factor
fw = sp.Function('w')
wave_F = [0, fw(t - x), 0, 0, 0, fw(t - x)]
check("B6a. a plane wave E_y = B_z = w(t - x) is null: E^2 - B^2 = 0 and E.B = 0 (it sits on the wall)",
      zero(inv1(wave_F)) and zero(inv2(wave_F)))
bw = boost_D(wave_F, D)
check("B6b. an observer moving along the wave sees amplitude x 1/D and energy density x 1/D^2 = 1/(1 + d(v))",
      zero(bw[1] - fw(t - x) / D) and zero(bw[5] - fw(t - x) / D) and zero((bw[1] ** 2 + bw[5] ** 2) / (2 * fw(t - x) ** 2) - 1 / D ** 2))
# Purcell: a pure electric field, boosted, acquires a magnetic field
pureE = [0, Ey, 0, 0, 0, 0]
pb = boost_D(pureE, D)
check("B7. magnetism from electricity: a pure E_y, boosted, gets B_z = -gamma v E_y (Purcell's picture)",
      zero(pb[5] + (D ** 2 - 1) / (2 * D) * Ey))

# =====================================================================================================
section("C. The gauge phase as a dial")
th1, th2, th = sp.symbols('theta1 theta2 theta', real=True)
a = sp.exp(sp.I * th) - 1
check("C1. a unit phase written as a = e^{i theta} - 1 satisfies a * conj(a) = 0: its conjugate is its mult-add inverse",
      zero(sp.simplify(ma(a, sp.conjugate(a)))))
ar, ai = sp.symbols('a_r a_i', real=True)
A_ = ar + sp.I * ai
check("C1b. in general a * conj(a) = |1 + a|^2 - 1, so 'phase' means exactly 'a * conj(a) = 0'",
      zero(sp.expand(ma(A_, sp.conjugate(A_)) - (sp.Abs(1 + A_) ** 2 - 1)).rewrite(sp.re)) or
      zero(sp.expand(ma(A_, sp.conjugate(A_)) - ((1 + ar) ** 2 + ai ** 2 - 1))))
check("C2. phases compose by mult-add: (e^{i t1} - 1) * (e^{i t2} - 1) = e^{i(t1+t2)} - 1 (a renaming of multiplication)",
      zero(sp.simplify(ma(sp.exp(sp.I * th1) - 1, sp.exp(sp.I * th2) - 1) - (sp.exp(sp.I * (th1 + th2)) - 1))))
check("C2b. the wall a = -1 fails the phase rule: (-1) * conj(-1) = -1, not 0, so U(1) never reaches the wall",
      ma(-1, -1) == -1)
tt1, tt2 = sp.tan(th1 / 2), sp.tan(th2 / 2)
check("C3. in the half-angle chart t = tan(theta/2) phases compose by the landscape law at q = +1: (t1+t2)/(1 - t1 t2)",
      zero(sp.simplify(sp.tan((th1 + th2) / 2) - law(tt1, tt2, 1)).rewrite(sp.tan)) or
      zero(sp.simplify(sp.expand_trig(sp.tan(th1 / 2 + th2 / 2)) - law(tt1, tt2, 1))))
p1, p2 = sp.symbols('phi1 phi2', real=True)
check("C3b. boosts in the half-rapidity chart s = tanh(phi/2) compose at the mirror dial q = -1: (s1+s2)/(1 + s1 s2)",
      zero(sp.simplify(sp.expand_trig(sp.tanh(p1 / 2 + p2 / 2)) - law(sp.tanh(p1 / 2), sp.tanh(p2 / 2), -1))))
# split-complex check: boosts as 1 + a with j^2 = +1, a * abar = 0
ph = sp.symbols('phi', real=True)
a0, a1j = sp.cosh(ph) - 1, sp.sinh(ph)          # a = a0 + j a1
# (a0 + j a1) * (a0 - j a1) mult-add, with j^2 = +1: a + abar + a abar = 2 a0 + (a0^2 - a1^2)
check("C4. a Lorentz boost as a split-complex number a = (cosh phi - 1) + j sinh phi also has a * abar = 0 "
      "(j^2 = +1): the same rule, mirror sign", zero(sp.simplify(2 * a0 + a0 ** 2 - a1j ** 2)))
check("C4b. the phase version with i^2 = -1: 2 a0 + a0^2 + a1^2 = 0 for a0 = cos - 1, a1 = sin",
      zero(sp.simplify(2 * (sp.cos(ph) - 1) + (sp.cos(ph) - 1) ** 2 + sp.sin(ph) ** 2)))

# C5 gauge invariance in 4D
phi_ = sp.Function('phi')(*X4)
Ak = [sp.Function(f'A{i}')(*X4) for i in 'xyz']
chi = sp.Function('chi')(*X4)


def EB_from_pot(ph0, Avec):
    grad = lambda fn: [sp.diff(fn, s_) for s_ in (x, y, z)]
    E = [-gp - sp.diff(Ai, t) for gp, Ai in zip(grad(ph0), Avec)]
    B = [sp.diff(Avec[2], y) - sp.diff(Avec[1], z), sp.diff(Avec[0], z) - sp.diff(Avec[2], x),
         sp.diff(Avec[1], x) - sp.diff(Avec[0], y)]
    return E, B


E0, B0 = EB_from_pot(phi_, Ak)
E1, B1 = EB_from_pot(phi_ - sp.diff(chi, t), [Ai + sp.diff(chi, s_) for Ai, s_ in zip(Ak, (x, y, z))])
check("C5. gauge invariance: A -> A + grad chi, phi -> phi - d_t chi leaves E and B unchanged",
      all(zero(a_ - b_) for a_, b_ in zip(E0 + B0, E1 + B1)))

# C6 holonomy = charge x flux: uniform B, A = B(-y, x)/2, exact on polygons
random.seed(7)
Bsym = sp.symbols('B', real=True)
okC6 = True
for trial in range(20):
    nv = random.randint(3, 7)
    angs = sorted(random.uniform(0, 2 * math.pi) for _ in range(nv))
    pts = [(sp.Rational(round(100 * (2 + random.random()) * math.cos(a_)), 37),
            sp.Rational(round(100 * (2 + random.random()) * math.sin(a_)), 41)) for a_ in angs]
    circ = 0
    area = 0
    for i in range(nv):
        (x0, y0), (x1, y1) = pts[i], pts[(i + 1) % nv]
        # straight segment, A linear: line integral is exact at the midpoint
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        circ += Bsym / 2 * (-ym * (x1 - x0) + xm * (y1 - y0))
        area += (x0 * y1 - x1 * y0) / 2
    okC6 &= zero(circ - Bsym * area)
check("C6. holonomy: the loop integral of A around 20 random rational polygons is exactly B x enclosed area "
      "(phase = charge x flux)", okC6)


def ab_loop(cx, cy, R, Phi=1.0):
    fx = lambda s_: cx + R * math.cos(s_)
    fy = lambda s_: cy + R * math.sin(s_)
    integrand = lambda s_: Phi / (2 * math.pi) * ((-fy(s_)) * (-R * math.sin(s_)) + fx(s_) * (R * math.cos(s_))) / (fx(s_) ** 2 + fy(s_) ** 2)
    return quad(integrand, 0, 2 * math.pi, limit=200)[0]


inside = [ab_loop(0.3, -0.2, 1.0), ab_loop(-0.5, 0.1, 2.0), ab_loop(0.0, 0.0, 0.1)]
outside = [ab_loop(3.0, 0.0, 1.0), ab_loop(-2.0, 2.5, 1.5)]
check("C6b. Aharonov-Bohm: around a thin solenoid (B = 0 on the path) the phase is the flux for any loop around it, "
      "0 for loops that miss it", all(abs(v_ - 1) < 1e-9 for v_ in inside) and all(abs(v_) < 1e-9 for v_ in outside),
      f"inside {[round(v_, 12) for v_ in inside]}, outside {[round(v_, 12) for v_ in outside]}")

# C7 lattice plaquette as mult-add of link values, and discrete Stokes
eps = sp.symbols('epsilon', positive=True)
cA = sp.symbols('c0:6', real=True)
Ax_lin = cA[0] + cA[1] * x + cA[2] * y
Ay_lin = cA[3] + cA[4] * x + cA[5] * y
Fxy = sp.diff(Ay_lin, x) - sp.diff(Ax_lin, y)


def link(Afun, p0, p1):
    xm, ym = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    dx_, dy_ = p1[0] - p0[0], p1[1] - p0[1]
    return sp.exp(sp.I * (Ax_lin.subs({x: xm, y: ym}) * dx_ + Ay_lin.subs({x: xm, y: ym}) * dy_)) - 1


inv = lambda u: -u / (1 + u)
corners = [(0, 0), (eps, 0), (eps, eps), (0, eps)]
l1 = link(None, corners[0], corners[1]); l2 = link(None, corners[1], corners[2])
l3 = link(None, corners[3], corners[2]); l4 = link(None, corners[0], corners[3])
plaq = ma(ma(ma(l1, l2), inv(l3)), inv(l4))
check("C7. a lattice plaquette is the mult-add of its four links (two run backwards, as mult-add inverses), "
      "and equals e^{i eps^2 F} - 1 exactly for a linear potential",
      zero(sp.simplify(plaq - (sp.exp(sp.I * eps ** 2 * Fxy) - 1))))
# discrete Stokes: two plaquettes side by side compose by mult-add into the 2x1 rectangle (shared link cancels)
cR = [(eps, 0), (2 * eps, 0), (2 * eps, eps), (eps, eps)]
m1 = link(None, cR[0], cR[1]); m2 = link(None, cR[1], cR[2]); m3 = link(None, cR[3], cR[2]); m4 = link(None, cR[0], cR[3])
plaq2 = ma(ma(ma(m1, m2), inv(m3)), inv(m4))
big = ma(ma(ma(ma(ma(l1, m1), m2), inv(m3)), inv(l3)), inv(l4))
check("C7b. discrete Stokes: plaquette * plaquette = the big loop (fluxes mult-add, the inner link cancels)",
      zero(sp.simplify(ma(plaq, plaq2) - big)))

# C8 non-commuting contrast: constant matrix potentials give a field only through the commutator
Pm = np.array([[0, 1], [1, 0]], dtype=complex) * 0.7j
Qm = np.array([[1, 0], [0, -1]], dtype=complex) * 0.4j
hol = []
for e_ in (1e-1, 5e-2, 2.5e-2):
    U = expm(e_ * Pm) @ expm(e_ * Qm) @ expm(-e_ * Pm) @ expm(-e_ * Qm)
    hol.append(np.linalg.norm((U - np.eye(2)) / e_ ** 2 - (Pm @ Qm - Qm @ Pm)))
check("C8. contrast (Yang-Mills): with constant but non-commuting potentials the loop still turns, by eps^2 x [A_x, A_y]; "
      "U(1) phases commute, so constant A gives no field", hol[2] < hol[1] < hol[0] and hol[2] < 0.05,
      f"error ratio {hol[0] / hol[1]:.2f}, {hol[1] / hol[2]:.2f} (about 2: first order)")
U1 = cmath.exp(0.1j * 0.7) * cmath.exp(0.1j * 0.4) * cmath.exp(-0.1j * 0.7) * cmath.exp(-0.1j * 0.4)
check("C8b. the U(1) version of the same loop closes exactly (phases commute)", abs(U1 - 1) < 1e-15)

# =====================================================================================================
section("D. F is the corner of a Clifford product; Maxwell is one equation")
va = sp.symbols('p0:4', real=True); vb = sp.symbols('r0:4', real=True)
Va = G.vector(list(va)); Vb = G.vector(list(vb))
ab = G.mul(Va, Vb); ba = G.mul(Vb, Va)
dotp = G.grade(ab, 0); wedge = G.grade(ab, 2)
check("D1. in spacetime algebra a b = a.b + a^b, and the order gap a b - b a = 2 a^b (the Clifford page's x*y - y*x)",
      G.equal(G.add(dotp, wedge), ab) and G.equal(G.comm(Va, Vb), G.scale(2, wedge)))
lift = G.mul(G.add(G.scalar(1), Va), G.add(G.scalar(1), Vb))
check("D1b. (1 + a)(1 + b) = 1 + a + b + a.b + a^b: the corner splits into a dot part and a wedge part",
      G.equal(lift, G.add(G.scalar(1), Va, Vb, dotp, wedge)))


def dmv(Xm, var):
    return G.clean({A: sp.diff(cf, var) for A, cf in Xm.items()})


def nabla(Xm):
    return G.add(*[G.mul({1 << m: sp.Integer(1) / G.k[m]}, dmv(Xm, X4[m])) for m in range(4)])


Apot = G.vector([phi_, Ak[0], Ak[1], Ak[2]])   # A = phi g0 + A^k gk
NA = nabla(Apot)
FA = Fmv(E0, B0)
check("D2. nabla A = nabla.A + nabla^A, and nabla^A is exactly F built from E = -grad phi - dA/dt, B = curl A",
      G.equal(G.grade(NA, 2), FA))
divA = sp.diff(phi_, t) + sum(sp.diff(Ak[i], s_) for i, s_ in enumerate((x, y, z)))
check("D2b. the scalar part nabla.A = d_t phi + div A (the Lorenz-gauge quantity)", zero(G.grade(NA, 0).get(0, 0) - divA))

Ef = [sp.Function(f'E{i}')(*X4) for i in 'xyz']
Bf = [sp.Function(f'B{i}')(*X4) for i in 'xyz']
Ff = Fmv(Ef, Bf)
NF = nabla(Ff)
divE = sum(sp.diff(Ef[i], s_) for i, s_ in enumerate((x, y, z)))
divB = sum(sp.diff(Bf[i], s_) for i, s_ in enumerate((x, y, z)))
curl = lambda V: [sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)]
cB, cE = curl(Bf), curl(Ef)
vec = G.grade(NF, 1)
tri = G.grade(NF, 3)
ok_vec = zero(vec.get(1, 0) - divE) and all(zero(vec.get(1 << (k + 1), 0) - (cB[k] - sp.diff(Ef[k], t))) for k in range(3))
# trivector components: e123 -> div B ; e0 e_i e_j -> -(curl E + dB/dt) on the missing axis
faraday = [cE[k] + sp.diff(Bf[k], t) for k in range(3)]
ok_tri = len(tri) == 4 and zero(tri.get(0b1110, 0) - divB)
for m in (1, 2, 3):
    blade = 0b1111 ^ (1 << m)           # e0 with the two axes other than m
    val = tri.get(blade, 0)
    ok_tri = ok_tri and (zero(val - faraday[m - 1]) or zero(val + faraday[m - 1]))
check("D3. nabla F = J holds iff Maxwell's four equations hold: vector part = (div E, curl B - dE/dt) = (rho, J); "
      "trivector part = (div B, curl E + dB/dt) = 0", ok_vec and ok_tri)
NNF = nabla(NF)
check("D4. charge conservation is automatic: the scalar part of nabla(nabla F) is identically 0, so nabla.J = 0",
      zero(G.grade(NNF, 0).get(0, 0)) and G.grade(NNF, 4) == {} and G.grade(NNF, 0) == {})
box = lambda fn: sp.diff(fn, t, 2) - sum(sp.diff(fn, s_, 2) for s_ in (x, y, z))
NNA = nabla(nabla(Apot))
check("D5. nabla nabla = the wave operator: nabla(nabla A) = (d_t^2 - laplacian) A, so in Lorenz gauge Maxwell is box A = J",
      G.equal(NNA, G.vector([box(phi_), box(Ak[0]), box(Ak[1]), box(Ak[2])])))
# D6 momentum-space form
ks = sp.symbols('k0:4', real=True)
Fs = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'F{min(i, j)}{max(i, j)}') * (1 if i < j else (-1 if i > j else 0)))
check("D6. for any antisymmetric F, k_mu k_nu F^{mu nu} = 0: in momentum space charge conservation k.J = 0 is forced",
      zero(sum(ks[i] * ks[j] * Fs[i, j] for i in range(4) for j in range(4))))
# D7: F alone as a multivector: F^2 = (E^2 - B^2) + 2 (E.B) I (the two invariants in one product)
F2 = G.mul(Fm, Fm)
sc = F2.get(0, 0); ps = F2.get(15, 0)
check("D7. F F = (E^2 - B^2) + 2 (E.B) I: the field's two invariants are the scalar and pseudoscalar parts of its square",
      zero(sc - inv1(FLD)) and (zero(ps - 2 * inv2(FLD)) or zero(ps + 2 * inv2(FLD))) and set(F2) <= {0, 15},
      f"pseudoscalar coefficient = {sp.factor(ps)}")
RES["D7_pseudo"] = str(sp.factor(ps))
# D8: plane wave solves nabla F = 0
Fw = Fmv([0, fw(t - x), 0], [0, 0, fw(t - x)])
check("D8. the plane wave solves nabla F = 0 (light in empty space), and F^2 = 0 (null)",
      nabla(Fw) == {} and G.mul(Fw, Fw) == {})

# =====================================================================================================
section("E. Where the lens stops")
F1w = Fmv([0, fw(t - x), 0], [0, 0, fw(t - x)])
gw = sp.Function('h')
F2w = Fmv([0, 0, gw(t + x)], [0, gw(t + x), 0])
check("E1. Maxwell in vacuum is linear: two solutions add by plain addition (k = 0, no corner term)",
      nabla(G.add(F1w, F2w)) == {} and nabla(F2w) == {})
Ee, bb = sp.symbols('E b', positive=True)
Lbi = bb ** 2 * (1 - sp.sqrt(1 - Ee ** 2 / bb ** 2))
Dbi = sp.diff(Lbi, Ee)
Dd_ = sp.symbols('D_e', positive=True)
Einv = sp.solve(sp.Eq(Dbi, Dd_), Ee)
check("E2. Born-Infeld (1934): D = E / sqrt(1 - E^2/b^2), so E = D / sqrt(1 + D^2/b^2) <= b: "
      "the same wall shape as speed from momentum, v = p / sqrt(1 + p^2)",
      zero(Dbi - Ee / sp.sqrt(1 - Ee ** 2 / bb ** 2)) and any(zero(e_ - Dd_ / sp.sqrt(1 + Dd_ ** 2 / bb ** 2)) for e_ in Einv))
Efun = lambda Dv: Dv / math.sqrt(1 + Dv * Dv)
d1_, d2_ = 0.8, 1.3
lhs = Efun(d1_ + d2_); rhs = (Efun(d1_) + Efun(d2_)) / (1 + Efun(d1_) * Efun(d2_))
check("E2b. but Born-Infeld fields do not compose by velocity addition (D adds, by Gauss's law): only the wall shape "
      "is shared, not the mult-add law", abs(lhs - rhs) > 1e-3, f"E(D1+D2) = {lhs:.4f} vs E1 (+) E2 = {rhs:.4f}")
me, e_ch, hbar = 9.1093837015e-31, 1.602176634e-19, 1.054571817e-34
Es = me ** 2 * c_SI ** 3 / (e_ch * hbar)
check("E3. the Schwinger field m^2 c^3/(e hbar) ~ 1.3e18 V/m: where QED makes the vacuum nonlinear "
      "(Maxwell's plain addition ends)", 1.30e18 < Es < 1.33e18, f"{Es:.4e} V/m")
RES["E3_schwinger_Vpm"] = Es

# =====================================================================================================
section("F. Bridge to gravity (side by side)")
# F1: the exterior derivative squares to zero: d(dA) = 0 is nabla^(nabla^A) = 0, the homogeneous half of Maxwell
check("F1. nabla^(nabla^A) = 0 for every potential (the 'Bianchi identity' of electromagnetism: no magnetic charge)",
      G.grade(nabla(G.grade(NA, 2)), 3) == {})
LINES.append("[NOTE] side by side, cited not re-checked: EM holonomy = charge x flux (C6); gravity's = rapidity (2M/r^3) x area "
             "(general-relativity page, route b). EM conservation (D4); gravity's Bianchi identity nabla.G = 0 (GR page).")

LINES.append(f"\n{npass} passed, {nfail} failed")
print(f"\n{npass} passed, {nfail} failed")
with open(os.path.join(HERE, "checks_output.txt"), "w") as fh:
    fh.write("\n".join(LINES) + "\n")
RES["passed"], RES["failed"] = npass, nfail
with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(RES, fh, indent=1, default=str)
