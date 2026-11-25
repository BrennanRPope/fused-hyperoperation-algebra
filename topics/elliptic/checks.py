"""Elliptic lead: exact and numeric checks. Run: python3 checks.py  (writes checks_output.txt, results.json)"""
import json, math, cmath, random
import sympy as sp
import mpmath as mp
import numpy as np
out, R = [], {}
npass = ntot = 0
def rep(name, ok, extra=""):
    global npass, ntot
    ntot += 1; npass += bool(ok)
    line = ("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else "")
    print(line); out.append(line)
def note(t):
    print("   " + t); out.append("   " + t)

# ============ 1. singular cubics carry the landscape laws; walls = branches of the node
a1, a2, t1, t2, t = sp.symbols('a1 a2 t1 t2 t')
def pt(tt):                      # line y = t x through the double point meets the curve again here
    X = tt**2 + a1*tt - a2
    return X, tt*X
X1, Y1 = pt(t1); X2, Y2 = pt(t2)
on = sp.simplify(Y1**2 + a1*X1*Y1 - X1**3 - a2*X1**2)
rep("1a every slope t through the double point gives a point (t^2+a1 t-a2, t(...)) on y^2+a1xy=x^3+a2x^2", on == 0)
lam = (Y2 - Y1)/(X2 - X1); nu = Y1 - lam*X1
X3 = sp.simplify(lam**2 + a1*lam - a2 - X1 - X2)
Y3 = lam*X3 + nu
Xs, Ys = X3, -Y3 - a1*X3          # P1 + P2 = -(third point)
zf = lambda X, Y: -X/Y
z1, z2 = zf(X1, Y1), zf(X2, Y2)
k, q = -a1, -a2
F = (z1 + z2 + k*z1*z2)/(1 - q*z1*z2)
rep("1b chord rule in z = -x/y is exactly the lens law (z1+z2+k z1 z2)/(1-d z1 z2), k=-a1, d=-a2",
    sp.simplify(zf(Xs, Ys) - F) == 0)
rep("1c z = -1/t, so z -> wall exactly where x(t) = t^2+a1 t-a2 = 0: the node's tangent slopes",
    sp.simplify(z1 + 1/t1) == 0 and sp.expand((t**2*(1 + k*(-1/t) + q/t**2)) - (t**2 + a1*t - a2)) == 0)
m1, m2 = sp.symbols('m1 m2')
lens = lambda tt: (tt - m1)/(tt - m2)
sub = {a1: -(m1 + m2), a2: -m1*m2}
rep("1d lens (t-m1)/(t-m2) of the two branches turns the chord law into multiplication",
    sp.simplify((lens(sp.Symbol('T3')) ).subs(sp.Symbol('T3'), -1/zf(Xs, Ys)).subs(sub) - (lens(t1)*lens(t2)).subs(sub)) == 0)
rep("1e discriminant of the tangent slopes a1^2+4a2 equals the map's k^2-4d", sp.expand(a1**2 + 4*a2 - (k**2 - 4*q)) == 0)
# points over F_p
def count_ns(A1, A2, p):
    n = 1  # infinity
    for xx in range(p):
        for yy in range(p):
            if (yy*yy + A1*xx*yy - xx**3 - A2*xx*xx) % p == 0: n += 1
    return n - 1  # drop the node (0,0)
cases = [("Circle-operation node y^2-xy=x^3 (k=1,d=0)", -1, 0), ("velocity-addition node y^2=x^3+x^2 (k=0,d=-1)", 0, 1),
         ("turning node y^2=x^3-x^2 (k=0,d=1)", 0, -1), ("cusp y^2=x^3", 0, 0)]
ok = True; rows = []
for name, A1, A2 in cases:
    for p in [5, 7, 11, 13]:
        n = count_ns(A1, A2, p)
        dsc = (A1*A1 + 4*A2) % p
        leg = 0 if dsc == 0 else (1 if pow(dsc, (p-1)//2, p) == 1 else -1)
        ok &= (n == p - leg)
        rows.append((name, p, n, leg))
rep("1f smooth points over F_p: p-1 when the walls are in F_p (two walls), p+1 when not (turning), p at the cusp", ok)
for r in rows[::4]: note(f"{r[0]}: p={r[1]} -> {r[2]} points")

# ============ 2. smoothing the node: Tate curve = the circle operation with 1+a ~ q(1+a)
mp.mp.dps = 30
def tate(qq, M=60):
    s = lambda kk: mp.fsum(n**kk*qq**n/(1 - qq**n) for n in range(1, M))
    A4 = -5*s(3); A6 = -(5*s(3) + 7*s(5))/12
    def P(u):
        xs = mp.fsum(qq**n*u/(1 - qq**n*u)**2 for n in range(-M, M)) - 2*s(1)
        ys = mp.fsum(qq**(2*n)*u**2/(1 - qq**n*u)**3 for n in range(-M, M)) + s(1)
        return xs, ys
    return A4, A6, P
qq = mp.mpf('0.05')
A4, A6, P = tate(qq)
random.seed(1)
err_curve = 0; err_add = 0; err_period = 0
for _ in range(6):
    u1 = mp.mpc(random.uniform(.3, .9), random.uniform(-.5, .5)); u2 = mp.mpc(random.uniform(.3, .9), random.uniform(-.5, .5))
    xa, ya = P(u1); xb, yb = P(u2); xc, yc = P(u1*u2)
    err_curve = max(err_curve, abs(ya**2 + xa*ya - xa**3 - A4*xa - A6))
    L = (yb - ya)/(xb - xa); x3 = L**2 + L - xa - xb; y3 = L*x3 + (ya - L*xa); y3n = -y3 - x3
    err_add = max(err_add, abs(x3 - xc) + abs(y3n - yc))
    xq, yq = P(qq*u1); err_period = max(err_period, abs(xq - xa) + abs(yq - ya))
rep("2a Tate's formula u -> (x(u), y(u)) lands on y^2+xy=x^3+a4(q)x+a6(q)  (q = 0.05, complex u)", err_curve < 1e-20, f"err {mp.nstr(err_curve,3)}")
rep("2b multiplying u is the chord-and-tangent sum: P(u1)+P(u2) = P(u1 u2)", err_add < 1e-15, f"err {mp.nstr(err_add,3)}")
rep("2c u and q u are the same point (the gluing that closes the cylinder into a torus)", err_period < 1e-20, f"err {mp.nstr(err_period,3)}")
u = sp.symbols('u')
x0 = u/(1 - u)**2; y0 = u**2/(1 - u)**3
rep("2d at q = 0 Tate's curve is the node y^2+xy=x^3 and z = -x/y = 1 - 1/u, so the lift 1+a (a = -z) is 1/u",
    sp.simplify(y0**2 + x0*y0 - x0**3) == 0 and sp.simplify(-x0/y0 - (1 - 1/u)) == 0)
rep("2e q = 0 law in z is a+b-ab, which is the circle operation after a -> -a", sp.expand(-((-t1) + (-t2) - (-t1)*(-t2)) - (t1 + t2 + t1*t2)) == 0)
Q = sp.symbols('q')
ser = lambda kk: sum(n**kk*Q**n/(1 - Q**n) for n in range(1, 5))
a4q = sp.series(-5*ser(3), Q, 0, 4).removeO(); a6q = sp.series(-(5*ser(3) + 7*ser(5))/12, Q, 0, 4).removeO()
note(f"a4(q) = {a4q} + ...,  a6(q) = {a6q} + ...")
jq = sp.series(1/Q + 744 + 196884*Q, Q, 0, 2)
# c4 and Delta of the Tate curve -> j = c4^3/Delta
c4 = 1 - 48*a4q; b2_, b4_, b6_ = 1, 2*a4q, 4*a6q
b8_ = a6q - a4q**2
Dl = -b2_**2*b8_ - 8*b4_**3 - 27*b6_**2 + 9*b2_*b4_*b6_
jser = sp.series(c4**3/Dl, Q, 0, 2).removeO()
rep("2f Tate curve's j-invariant is 1/q + 744 + 196884 q + ...: the whole landscape (q = 0) is the single point j = infinity",
    sp.expand(jser - (1/Q + 744 + 196884*Q)) == 0, str(jser))
dsc = sp.series(Dl, Q, 0, 3).removeO()
rep("2g its discriminant is q - 24 q^2 + ...: the smoothing parameter is q itself", sp.expand(dsc - (Q - 24*Q**2)) == 0, str(dsc))

# formal group of the Tate curve to first order in q, via logs (from rung_test machinery, lighter)
zz, w = sp.symbols('z w')
def omega_series(A1, A2, A3, A4_, A6_, N=8):
    W = zz**3
    for _ in range(N):
        W = sp.series(zz**3 + A1*zz*W + A2*zz**2*W + A3*W**2 + A4_*zz*W**2 + A6_*W**3, zz, 0, N + 3).removeO()
    Xz = sp.expand(zz/W); Yz = sp.expand(-1/W)
    return sp.series(sp.diff(Xz, zz)/(2*Yz + A1*Xz + A3), zz, 0, N).removeO()
eps = sp.symbols('eps')
om = omega_series(1, 0, 0, -5*eps, -eps)
om1 = sp.expand(sp.series(om, eps, 0, 2).removeO())
om0 = sp.expand(om1.subs(eps, 0)); omq = sp.expand(om1.coeff(eps, 1))
rep("2h q = 0 differential is dz/(1-z) (the flipped the circle operation's log)", sp.expand(om0 - sp.series(1/(1 - zz), zz, 0, 8).removeO()) == 0)
lowest = min(sp.Poly(omq, zz).monoms())[0]
note(f"first-order-in-q part of the invariant differential: {sp.factor(omq)}")
rep("2i q first enters at z^4 of the differential, i.e. the cube's 5th layer of the law (a3 = 0 on Tate's curve)", lowest == 4)

# ============ 3. connections
# 3a pendulum / sn addition third layer
m_ = sp.symbols('m', positive=True); A, B = sp.symbols('A B')
cn = lambda v: sp.sqrt(1 - v**2); dn = lambda v: sp.sqrt(1 - m_*v**2)
snadd = (A*cn(B)*dn(B) + B*cn(A)*dn(A))/(1 - m_*A**2*B**2)
h = sp.symbols('h')
ser3 = sp.expand(sp.series(snadd.subs({A: h*A, B: h*B}), h, 0, 4).removeO().subs(h, 1))
rep("3a Jacobi sn addition (pendulum phase law) = a + b - (1+m)/2 (a^2 b + a b^2) + ...: corner 0, bend d = -(1+m)/2",
    sp.expand(ser3 - (A + B - (1 + m_)/2*(A**2*B + A*B**2))) == 0)
note("so m = 0 (small swings, sine) sits at (0,-1/2) and m = 1 (separatrix) at (0,-1), velocity addition")
s1 = (A*cn(B)*dn(B) + B*cn(A)*dn(A)).subs(m_, 1); den1 = (1 - A**2*B**2)
rep("3b at m = 1 the sn law is exactly velocity addition (a+b)/(1+ab): the separatrix is the node",
    sp.simplify(s1/den1 - (A + B)/(1 + A*B)) == 0)
# 3c cnoidal -> soliton (KdV u_t + 6 u u_x + u_xxx = 0)
mp.mp.dps = 25
def cnoidal(xv, tv, kk, m, u0=0):
    c = 6*u0 + 4*kk**2*(2*m - 1)
    sn, cn_, dn_ = mp.ellipfun('sn', k*0 + kk*(xv - c*tv), m=m), None, None
    return u0 + 2*kk**2*m*mp.ellipfun('cn', kk*(xv - c*tv), m=m)**2
kk = mp.mpf(1)
worst = 0
for m in [mp.mpf('0.5'), mp.mpf('0.9')]:
    for xv in [mp.mpf('0.3'), mp.mpf('1.1')]:
        f = lambda X, T: cnoidal(X, T, kk, m)
        ut = mp.diff(lambda T: f(xv, T), 0); ux = mp.diff(lambda X: f(X, 0), xv)
        uxxx = mp.diff(lambda X: f(X, 0), xv, 3)
        worst = max(worst, abs(ut + 6*f(xv, 0)*ux + uxxx))
rep("3c cnoidal wave 2k^2 m cn^2(k(x-ct)|m) solves KdV (residual at sample points)", worst < 1e-10, f"{mp.nstr(worst,3)}")
rows3 = []
for m in ['0.9', '0.99', '0.9999', '0.99999999']:
    m = mp.mpf(m)
    gap = max(abs(cnoidal(X, 0, kk, m) - 2*kk**2*mp.sech(X)**2) for X in [mp.mpf(v)/4 for v in range(0, 13)])
    period = 2*mp.ellipk(m)/kk
    qt = mp.exp(-mp.pi*mp.ellipk(m)/mp.ellipk(1 - m))
    rows3.append((float(m), float(period), float(qt), float(gap)))
ok = all(rows3[i][3] > rows3[i+1][3] for i in range(3)) and rows3[-1][3] < 1e-6
rep("3d as the spacing grows (m -> 1, Tate parameter q -> 0) the cnoidal train becomes the soliton 2k^2 sech^2", ok)
for r in rows3: note(f"m={r[0]}: spacing {r[1]:.2f}, q = {r[2]:.2e}, max gap to sech^2 {r[3]:.2e}")
# 3e elliptic genus degenerates to the signature law
dl = sp.symbols('delta')
ochF = (A*sp.sqrt(1 - 2*dl*B**2 + dl**2*B**4) + B*sp.sqrt(1 - 2*dl*A**2 + dl**2*A**4))/(1 - dl**2*A**2*B**2)
red = (A*(1 - dl*B**2) + B*(1 - dl*A**2))/(1 - dl**2*A**2*B**2)
rep("3e Ochanine's elliptic genus law at eps = delta^2 is the lens law (a+b)/(1+delta ab): delta = 1 is the signature (velocity addition)",
    sp.simplify(red - (A + B)/(1 + dl*A*B)) == 0)
# 3f Frey curve: nodal (semistable) at odd primes dividing ABC
Af, Bf = sp.symbols('A B')
b2f, b4f = 4*(Bf - Af), -2*Af*Bf
c4f = sp.expand(b2f**2 - 24*b4f)
bb, cc = Bf - Af, -Af*Bf
disc = 16*(bb**2*cc**2 - 4*cc**3)
rep("3f Frey curve y^2 = x(x-A)(x+B): discriminant 16 A^2 B^2 (A+B)^2 and c4 = 16(A^2+AB+B^2)",
    sp.expand(disc - 16*Af**2*Bf**2*(Af + Bf)**2) == 0 and sp.expand(c4f - 16*(Af**2 + Af*Bf + Bf**2)) == 0)
def node_type(A_, B_, p):
    # singular point of y^2 = x(x-A)(x+B) mod p and its tangent discriminant
    f = lambda xx: xx*(xx - A_)*(xx + B_)
    for xx in range(p):
        if f(xx) % p == 0 and (3*xx*xx + 2*(B_ - A_)*xx - A_*B_) % p == 0:
            hess = (6*xx + 2*(B_ - A_)) % p      # y^2 = (hess/2) (x-x0)^2 + ...  tangents y = +-sqrt(hess/2)(x-x0)
            if hess == 0: return "cusp"
            v = hess*pow(2, -1, p) % p
            return "split node (two walls)" if pow(v, (p-1)//2, p) == 1 else "non-split node (turning)"
    return "smooth"
ex = [(5, 27), (1, 8), (7, 9)]
ok = True
for A_, B_ in ex:
    for p in sorted(sp.primefactors(A_*B_*(A_ + B_))):
        if p == 2: continue
        tpe = node_type(A_, B_, p); ok &= tpe != "cusp"
        note(f"A={A_}, B={B_}, A+B={A_+B_}: mod {p} -> {tpe}")
rep("3g at every odd bad prime of these Frey-type curves the reduction is a node, never a cusp (off the border)", ok)

# ============ 4. Frobenius on the map: Euler factors are wall polynomials
curves = {"11a1": (0, -1, 1, -10, -20, 11), "37a1": (0, 0, 1, -1, 0, 37), "14a1": (1, 0, 1, 4, -6, 14),
          "15a1": (1, 1, 1, -10, -10, 15)}
def ap(c, p):
    A1_, A2_, A3_, A4_, A6_, _ = c
    if p == 2:
        n = 1 + sum(1 for xx in range(2) for yy in range(2) if (yy*yy + A1_*xx*yy + A3_*yy - xx**3 - A2_*xx*xx - A4_*xx - A6_) % 2 == 0)
        return p + 1 - n
    xs = np.arange(p, dtype=np.int64)
    D = ((A1_*xs + A3_)**2 + 4*(xs**3 + A2_*xs**2 + A4_*xs + A6_)) % p
    sq = np.zeros(p, dtype=np.int64); sq[(xs*xs) % p] = 1
    leg = np.where(D == 0, 0, np.where(sq[D] == 1, 1, -1))
    n = 1 + int(np.sum(1 + leg))
    return p + 1 - n
primes = list(sp.primerange(2, 4000))
frob = {}
ok_h = True; ok_bad = True
for name, c in curves.items():
    N_ = c[5]; rows = []
    for p in primes:
        a = ap(c, p)
        bad = N_ % p == 0
        rows.append((p, a, bad))
        if bad: ok_bad &= a in (-1, 0, 1)
        else: ok_h &= a*a < 4*p
    frob[name] = rows
    note(f"{name}: bad primes " + ", ".join(f"p={p}: a_p={a} ({'split, k=-1' if a==1 else 'non-split, k=+1 (the circle operation)' if a==-1 else 'additive'})" for p, a, b in rows if b))
rep("4a good primes: a_p^2 < 4p for all four curves, p < 4000 (Hasse): Euler factor 1 - a_p T + p T^2 sits in the turning country", ok_h)
rep("4b bad primes: Euler factor 1 - a_p T lands on the d = 0 line at k = -a_p in {-1, 0, +1}", ok_bad)
th = [math.acos(a/(2*math.sqrt(p))) for nm in ["11a1", "37a1", "14a1", "15a1"] for p, a, b in frob[nm] if not b]
hist, edges = np.histogram(th, bins=12, range=(0, math.pi))
cdf = lambda x: (x - math.sin(x)*math.cos(x))/math.pi
expct = np.array([cdf(edges[i+1]) - cdf(edges[i]) for i in range(12)])*len(th)
chi2 = float(np.sum((hist - expct)**2/expct))
rep("4c turning angles follow Sato-Tate (2/pi) sin^2: chi-square over 12 bins", chi2 < 30, f"chi2 = {chi2:.1f} (11 dof), {len(th)} angles")
R["frob"] = frob; R["sato_tate"] = {"hist": hist.tolist(), "expected": expct.tolist(), "edges": edges.tolist()}
R["cnoidal"] = rows3
print(f"\n{npass} of {ntot} checks pass"); out.append(f"\n{npass} of {ntot} checks pass")
open("checks_output.txt", "w").write("\n".join(out) + "\n")
json.dump(R, open("results.json", "w"))
