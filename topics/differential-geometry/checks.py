# Differential geometry through the circle operation: every claim on the page, checked.
import sympy as sp, numpy as np, math, random
from mpmath import mp, mpf, quad, log as mlog, atanh as matanh, atan as matan, sqrt as msqrt, pi as mpi
from scipy import integrate
mp.dps = 30
random.seed(3)
ok = 0; tot = 0
def check(name, cond):
    global ok, tot
    tot += 1; ok += bool(cond); print(('PASS ' if cond else 'FAIL ') + name)
z = lambda e: sp.simplify(sp.expand(e)) == 0
a, b, c, x, y, r, t, K, p, q, k = sp.symbols('a b c x y r t K p q k', real=True)
th1, th2 = sp.symbols('theta1 theta2', real=True)
F = lambda p, q, a, b: (a + b + p*a*b)/(1 - q*a*b)
V = lambda p, q, x: 1 + p*x + q*x**2

# --- 1. the three conics and their central charts ---
c1, s1, c2, s2 = sp.symbols('c1 s1 c2 s2', real=True)
C12, S12 = c1*c2 - K*s1*s2, s1*c2 + c1*s2
check("rotations of t² + Kx² compose by (c1c2 - K s1s2, s1c2 + c1s2) and keep the form",
      z(C12**2 + K*S12**2 - (c1**2 + K*s1**2)*(c2**2 + K*s2**2)))
r1, r2 = sp.symbols('r1 r2', real=True)
check("central chart r = x/t: the rotated point's chart value is (r1 + r2)/(1 - K r1 r2), the law F(0, K)",
      z((S12/C12).subs({s1: r1*c1, s2: r2*c2}) - F(0, K, r1, r2)))
check("K = 0: t = 1 stays put and chart values just add (flat, Galilean)", z(S12.subs({K: 0, c1: 1, c2: 1}) - (s1 + s2)) and C12.subs({K: 0, c1: 1, c2: 1}) == 1)
Kp = sp.symbols('Kp', positive=True)
cs = lambda th: (sp.cos(sp.sqrt(Kp)*th), sp.sin(sp.sqrt(Kp)*th)/sp.sqrt(Kp))
check("K > 0: (cos √K θ, sin √K θ / √K) is on the ellipse and turning angles add (rotation)",
      z(cs(th1)[0]**2 + Kp*cs(th1)[1]**2 - 1) and
      z(sp.expand_trig(cs(th1+th2)[0] - (cs(th1)[0]*cs(th2)[0] - Kp*cs(th1)[1]*cs(th2)[1]))))
ch = lambda th: (sp.cosh(sp.sqrt(Kp)*th), sp.sinh(sp.sqrt(Kp)*th)/sp.sqrt(Kp))
check("K < 0: (cosh √|K| θ, sinh √|K| θ / √|K|) is on the hyperbola and boost rapidities add",
      z(ch(th1)[0]**2 - Kp*ch(th1)[1]**2 - 1) and
      z(sp.expand_trig(ch(th1+th2)[0] - (ch(th1)[0]*ch(th2)[0] + Kp*ch(th1)[1]*ch(th2)[1]))))
check("arc length from the centre of the chart: θ = arctan(√K r)/√K, i.e. dθ = dr/(1 + K r²) = dr/V(r)",
      z(sp.diff(sp.atan(sp.sqrt(Kp)*r)/sp.sqrt(Kp), r) - 1/V(0, Kp, r)) and z(sp.diff(sp.atanh(sp.sqrt(Kp)*r)/sp.sqrt(Kp), r) - 1/V(0, -Kp, r)))

# --- 2. Gaussian curvature (Brioschi) ---
def gauss_curv(E, Fm, G, u, v):
    Eu, Ev, Fu, Fv, Gu, Gv = [sp.diff(f, w) for f in (E, Fm, G) for w in (u, v)]
    Evv, Guu, Fuv = sp.diff(E, v, 2), sp.diff(G, u, 2), sp.diff(Fm, u, v)
    M1 = sp.Matrix([[-Evv/2 + Fuv - Guu/2, Eu/2, Fu - Ev/2], [Fv - Gu/2, E, Fm], [Gv/2, Fm, G]])
    M2 = sp.Matrix([[0, Ev/2, Gu/2], [Ev/2, E, Fm], [Gu/2, Fm, G]])
    return (M1.det() - M2.det())/(E*G - Fm**2)**2
u_, v_ = sp.symbols('u v', real=True)
den = 1 + K*(u_**2 + v_**2)
E_b = ((1 + K*v_**2))/den**2; F_b = -K*u_*v_/den**2; G_b = (1 + K*u_**2)/den**2
check("central (Beltrami-Klein / gnomonic) chart metric has constant Gaussian curvature K",
      z(sp.simplify(gauss_curv(E_b, F_b, G_b, u_, v_) - K)))
check("its radial part is dr²/(1 + K r²)²·... : along a line through the centre ds = dr/(1 + K r²)",
      z(sp.simplify((E_b*1).subs(v_, 0) - 1/(1 + K*u_**2)**2)))
lam = 4/(1 + K*(u_**2 + v_**2))**2
check("conformal (Poincaré / stereographic) metric 4|dz|²/(1 + K|z|²)² has curvature K",
      z(sp.simplify(gauss_curv(lam, 0, lam, u_, v_) - K)))
check("conformal radial distance 2·arctan(√K r)/√K: the same tan law in r, with distances doubled",
      z(sp.diff(2*sp.atan(sp.sqrt(Kp)*r)/sp.sqrt(Kp), r) - 2/(1 + Kp*r**2)))
check("upper half-plane (dx² + dy²)/y² has curvature -1", z(sp.simplify(gauss_curv(1/v_**2, 0, 1/v_**2, u_, v_) + 1)))

# --- 3. doubling links the two kinds of chart ---
dbl = lambda q, s: 2*s/(1 - q*s**2)
check("doubling 2a/(1 - qa²) is F(a, a)", z(dbl(q, a) - F(0, q, a, a)))
check("doubling respects the law: dbl(F(a,b)) = F(dbl a, dbl b)", z(sp.simplify(dbl(q, F(0, q, a, b)) - F(0, q, dbl(q, a), dbl(q, b)))))
rP = sp.symbols('rP', positive=True)
check("K = -1: Klein radius = 2 r_P/(1 + r_P²) of the Poincaré radius (both give distance artanh r_K = 2 artanh r_P)",
      abs(matanh(2*mpf('0.3')/(1 + mpf('0.09'))) - 2*matanh(mpf('0.3'))) < 1e-25)
check("K = +1: gnomonic radius = 2 r_s/(1 - r_s²) of the stereographic radius (tan θ = 2 tan(θ/2)/(1 - tan²(θ/2)))",
      z(sp.expand_trig(sp.tan(2*th1)) - dbl(1, sp.tan(th1))))

# --- 4. curvature read off any landscape point ---
pp, qq = p + 2*c, q + p*c + c**2
check("re-charting x ↦ x/(1-cx) keeps 4q - p² (so K_proj = q - p²/4 and K_conf = 4q - p² are chart-free)",
      z(4*qq - pp**2 - (4*q - p**2)))
Lmu = lambda s: sp.log(1 + s)
check("The circle operation re-charted (c = -1/2) is F(0, -1/4), and its ruler log(1+x) = 2·artanh(φ(x)/2), φ(x) = x/(1 + x/2)",
      (pp.subs({p: 1, q: 0, c: sp.Rational(-1, 2)}), qq.subs({p: 1, q: 0, c: sp.Rational(-1, 2)})) == (0, sp.Rational(-1, 4))
      and abs(mlog(1 + mpf('0.7')) - 2*matanh((mpf('0.7')/(1 + mpf('0.35')))/2)) < 1e-25)
check("so the circle operation is a hyperbolic line: K_proj = q - p²/4 = -1/4 (Klein chart), K_conf = 4q - p² = -1, matching the half plane's curvature -1",
      (q - p**2/4).subs({p: 1, q: 0}) == sp.Rational(-1, 4) and (4*q - p**2).subs({p: 1, q: 0}) == sp.simplify(gauss_curv(1/v_**2, 0, 1/v_**2, u_, v_)))
# half-plane picture of circle-operation
check("half plane: lift y = 1 + x; the dilation y ↦ (1+a)y is an isometry and composes as the circle operation",
      z(((1 + a)*(1 + x)) - (1 + (a + x + a*x))) and z(sp.simplify(((1 + a)**2*(1 + 0))/((1 + a)*v_)**2 - 1/v_**2)))
check("distance up the vertical geodesic from y = 1 to y = 1 + x is log(1 + x), the circle-operation ruler; the wall x = -1 is the horizon y = 0, infinitely far",
      z(sp.integrate(1/v_, (v_, 1, 1 + x)).subs(x, 2) - sp.log(3)) and sp.limit(sp.log(1 + x), x, -1, '+') == -sp.oo)
check("dial line (p = k, q = 0): K_conf = -k² (radius 1/|k|), K_proj = -k²/4 (radius 2/|k|); flat only at k = 0",
      z((4*q - p**2).subs({p: k, q: 0}) + k**2) and z((q - p**2/4).subs({p: k, q: 0}) + k**2/4) and sp.solve(sp.Eq(-k**2, 0), k) == [0])
cl = sp.symbols('cl', positive=True)
check("velocity addition (u+v)/(1+uv/c²) is F(0, -1/c²): K_proj = -1/c², velocity space is hyperbolic with radius c",
      z(F(0, -1/cl**2, a, b) - (a + b)/(1 + a*b/cl**2)) and sp.sqrt(-1/((q - p**2/4).subs({p: 0, q: -1/cl**2}))) == cl)
check("failure union (k = -1) and parity (k = -2): K_conf = -1 and -4, both hyperbolic", [(4*q - p**2).subs({p: kv, q: 0}) for kv in (-1, -2)] == [-1, -4])

# --- 5. Cayley-Klein: the walls are the absolute ---
def CR(a_, b_, P, Q): return ((a_ - P)*(b_ - Q))/((a_ - Q)*(b_ - P))
good = True
for (pv, qv) in [(0, -1), (1, -2), (3, 2), (mpf('0.5'), mpf('-0.3')), (2, mpf('0.75'))]:
    pv, qv = mpf(pv), mpf(qv); D = pv*pv - 4*qv; rD = msqrt(D)
    u = (-pv + rD)/2; v = (-pv - rD)/2
    for xv in [mpf('0.1'), mpf('0.2'), mpf('-0.15')]:
        Lnum = quad(lambda s: 1/(1 + pv*s + qv*s*s), [0, xv])
        Lcr = mlog(CR(0, xv, 1/u, 1/v))/rD
        if abs(Lnum - Lcr) > 1e-15: good = False
check("ruler = Cayley-Klein distance: L(x) = log CR(0, x; walls)/√(p² - 4q) (5 line-type points)", good)
W_ = sp.symbols('W', positive=True)
check("The circle operation: walls -1 and ∞, and log CR(0, x; ∞, -1) = log(1 + x)", z(sp.limit(CR(0, x, W_, -1), W_, sp.oo) - (1 + x)))
check("Klein disk (K = -1): distance between two points on a chord = ½|log cross-ratio| with the chord's ends on the circle, = length of the straight segment in the Beltrami metric",
      (lambda P1, P2: (lambda d_seg, d_cr: abs(d_seg - d_cr) < 1e-9)(
          integrate.quad(lambda s: math.sqrt(((lambda X, Y, dX, dY: ((1 - (X*X + Y*Y))*(dX*dX + dY*dY) + (X*dX + Y*dY)**2)/(1 - (X*X + Y*Y))**2)(
              P1[0] + s*(P2[0]-P1[0]), P1[1] + s*(P2[1]-P1[1]), P2[0]-P1[0], P2[1]-P1[1]))), 0, 1, epsabs=1e-12)[0],
          (lambda ends: 0.5*abs(math.log(abs(np.linalg.norm(P1-ends[1])*np.linalg.norm(P2-ends[0])/(np.linalg.norm(P1-ends[0])*np.linalg.norm(P2-ends[1]))))))(
              (lambda dvec: [P1 + tt*dvec for tt in np.roots([dvec@dvec, 2*P1@dvec, P1@P1 - 1])])(P2 - P1))))(np.array([0.2, -0.3]), np.array([-0.5, 0.4])))
good = True
for pv, qv in [(0, 1), (1, 1), (-2, 3), (0, 4)]:
    Kp_ = qv - pv*pv/4
    good = good and abs(quad(lambda s: 1/(1 + pv*s + qv*s*s), [-mp.inf, mp.inf]) - mpi/msqrt(Kp_)) < 1e-20
check("circle-type law: whole line has length π/√K_proj, half a great circle (the elliptic plane, antipodes glued)", good)
check("tan law: the point ∞ is a quarter turn from 0, and F(a, ∞) = -1/a is the perpendicular direction", sp.limit(F(0, 1, a, x), x, sp.oo) == -1/a)

# --- 6. Einstein and Möbius addition (Ungar), collinear part ---
def einstein(u, v):
    uv = u @ v; g = 1/math.sqrt(1 - u @ u)
    return (u + v/g + (g/(1 + g))*uv*u)/(1 + uv)
def boost_klein(u, v):
    # Lorentz boost by velocity u applied to the event (1, v); return its Klein (velocity) coordinate
    g = 1/math.sqrt(1 - u @ u); n = u/np.linalg.norm(u)
    t0, xv = 1.0, v.copy()
    xpar = xv @ n
    t1 = g*(t0 + np.linalg.norm(u)*xpar)
    x1 = xv + ((g - 1)*xpar + g*np.linalg.norm(u)*t0)*n
    return x1/t1
U = np.array([0.5, 0.2]); Vv = np.array([-0.3, 0.6])
check("Einstein addition is a translation of the Klein disk: boost(u) sends velocity v to u ⊕_E v (2-D example)",
      np.allclose(einstein(U, Vv), boost_klein(U, Vv)))
mob_add = lambda A, B: (A + B)/(1 + np.conj(A)*B)
check("collinear: Einstein and Möbius addition both reduce to (a + b)/(1 + ab)",
      abs(einstein(np.array([0.4, 0.]), np.array([0.5, 0.]))[0] - 0.9/1.2) < 1e-12 and abs(mob_add(0.4+0j, 0.5+0j) - 0.9/1.2) < 1e-12)
check("Möbius addition z ↦ (a + z)/(1 + ā z) keeps the Poincaré disk metric (isometry, numeric)",
      (lambda A, Z, dZ: abs((abs((mob_add(A, Z + 1e-7*dZ) - mob_add(A, Z))/1e-7)/(1 - abs(mob_add(A, Z))**2)) - abs(dZ)/(1 - abs(Z)**2)) < 1e-5)(0.3+0.4j, -0.2+0.1j, 1+0.5j))
dbl2 = lambda Z: 2*Z/(1 + abs(Z)**2)
half = lambda W: W/(1 + math.sqrt(1 - abs(W)**2))
A_, B_ = 0.3 + 0.4j, -0.5 + 0.2j
ein_c = lambda A, B: (lambda w: w[0] + 1j*w[1])(einstein(np.array([A.real, A.imag]), np.array([B.real, B.imag])))
check("Ungar: doubling carries Möbius addition to Einstein addition in every direction: dbl(½a ⊕_M ½b) = a ⊕_E b (numeric)",
      abs(dbl2(mob_add(half(A_), half(B_))) - ein_c(A_, B_)) < 1e-12)

# --- 7. Gauss-Bonnet (only mentioned on the page) ---
P = [np.array([0.1, 0.1]), np.array([0.7, -0.2]), np.array([-0.3, 0.6])]
def hdist(A, B):
    return math.acosh((1 - A @ B)/math.sqrt((1 - A @ A)*(1 - B @ B)))
sides = [hdist(P[1], P[2]), hdist(P[0], P[2]), hdist(P[0], P[1])]
def ang(i):
    a_, b_, c_ = sides[i], sides[(i+1) % 3], sides[(i+2) % 3]
    return math.acos((math.cosh(b_)*math.cosh(c_) - math.cosh(a_))/(math.sinh(b_)*math.sinh(c_)))
angsum = sum(ang(i) for i in range(3))
def area_klein(P):
    A, B, C = P
    def integrand(s, tt):
        X = A + s*(B - A) + tt*(C - A)
        return abs((B-A)[0]*(C-A)[1]-(B-A)[1]*(C-A)[0])/(1 - X @ X)**1.5
    return integrate.dblquad(lambda tt, s: integrand(s, tt), 0, 1, 0, lambda s: 1 - s, epsabs=1e-11)[0]
check("Gauss-Bonnet, hyperbolic (K = -1): π - (angle sum) = area of a geodesic triangle (Klein-disk example)",
      abs((math.pi - angsum) - area_klein(P)) < 1e-7)
check("Gauss-Bonnet, sphere (K = +1): the octant triangle has three right angles and area π/2 = excess",
      abs(3*math.pi/2 - math.pi - 4*math.pi/8) < 1e-15)

print(f"\n{ok} of {tot} checks pass")
