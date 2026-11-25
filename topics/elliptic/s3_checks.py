"""Why the good primes' turning angles live on the 3-sphere, in circle-operation terms. Run: python3 s3_checks.py"""
import json, math, itertools, random
import numpy as np, sympy as sp
from checks import ap as ap_count   # point counting from checks.py (imports run its checks once)
out, R = [], {}
npass = ntot = 0
def rep(name, ok, extra=""):
    global npass, ntot
    ntot += 1; npass += bool(ok)
    line = ("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""); print(line); out.append(line)
def note(t): print("   " + t); out.append("   " + t)

# 1. every law on the map is the circle operation's lift with 1 + a promoted to the matrix I + aG
a, b, k, q = sp.symbols('a b k q')
G = sp.Matrix([[k, 1], [-q, 0]]); I = sp.eye(2)
F = (a + b + k*a*b)/(1 - q*a*b)
lhs = (I + a*G)*(I + b*G); rhs = (1 - q*a*b)*(I + F*G)
rep("1a (I + aG)(I + bG) = (1 - d ab)(I + F(a,b) G) with G = [[k,1],[-d,0]]: the lens law is the circle operation's lift in 2x2 matrices",
    sp.simplify(lhs - rhs) == sp.zeros(2))
rep("1b Cayley-Hamilton G^2 = kG - dI: trace k, determinant d", sp.simplify(G*G - (k*G - q*I)) == sp.zeros(2))
rep("1c circle operation is the 1x1 case G = [1] (k=1, d=0): (1+a)(1+b) = 1 + (a+b+ab)", sp.expand((1 + a)*(1 + b) - (1 + (a + b + a*b))) == 0)

# 2. primes: G_p is (minus) Frobenius on the torus's two loops
curves = {"11a1": (0, -1, 1, -10, -20, 11), "37a1": (0, 0, 1, -1, 0, 37), "14a1": (1, 0, 1, 4, -6, 14),
          "15a1": (1, 1, 1, -10, -10, 15)}
cm = {"32a (y^2=x^3-x, CM by i)": (0, 0, 0, -1, 0, 32), "27a (y^2+y=x^3, CM by cube roots of 1)": (0, 0, 1, 0, 0, 27)}
primes = list(sp.primerange(5, 6000))
def data(c):
    return [(p, ap_count(c, p)) for p in primes if c[5] % p]
noncm = [x for c in curves.values() for x in data(c)]
cmd = [x for c in cm.values() for x in data(c)]
ok = all(abs(np.linalg.eigvals(np.array([[-a_, 1], [-p, 0]], float)/math.sqrt(p))).round(9).tolist() == [1.0, 1.0] for p, a_ in noncm[:3000])
rep("2a for every good prime, G_p/sqrt(p) (trace -a_p/sqrt p, det 1) has both eigenvalues on the unit circle: a turn", ok)
# 3. bad primes: Frobenius on the lift 1+a is the p-th power: identity (split) or circle inverse (non-split)
ok = True
for p in [7, 11, 37, 101]:
    ok &= all(pow(x, p, p) == x for x in range(1, p))                         # split torus F_p*
    d = next(n for n in range(2, p) if pow(n, (p-1)//2, p) == p - 1)            # F_p2 = F_p[sqrt d]
    mul = lambda u, v: ((u[0]*v[0] + d*u[1]*v[1]) % p, (u[0]*v[1] + u[1]*v[0]) % p)
    def pw(u, n):
        r = (1, 0)
        while n:
            if n & 1: r = mul(r, u)
            u = mul(u, u); n >>= 1
        return r
    norm1 = [(x, y) for x in range(p) for y in range(p) if (x*x - d*y*y) % p == 1]
    ok &= len(norm1) == p + 1
    ok &= all(mul(pw(u, p), u) == (1, 0) for u in norm1)
rep("3a on the split torus u^p = u (Frobenius does nothing, a_p = +1); on the turning torus u^p = 1/u (it flips to the circle inverse, a_p = -1)", ok)

# 4. the 3-sphere: unit quaternions = SU(2); a turn by theta is a latitude 2-sphere of radius sin(theta)
w, x, y, z = sp.symbols('w x y z', real=True)
U = sp.Matrix([[w + sp.I*x, y + sp.I*z], [-y + sp.I*z, w - sp.I*x]])
rep("4a the SU(2) matrix of a quaternion w+xi+yj+zk has trace 2w and determinant w^2+x^2+y^2+z^2",
    sp.expand(U.trace() - 2*w) == 0 and sp.expand(U.det() - (w**2 + x**2 + y**2 + z**2)) == 0)
rep("4b its wall discriminant (2w)^2 - 4 det = -4(x^2+y^2+z^2): depth into the turning country = size of its latitude sphere",
    sp.expand((2*w)**2 - 4*(w**2 + x**2 + y**2 + z**2) + 4*(x**2 + y**2 + z**2)) == 0)
th = sp.symbols('theta')
vol = sp.integrate(4*sp.pi*sp.sin(th)**2, (th, 0, sp.pi))
rep("4c latitude spheres 4 pi sin^2(theta) add up to the 3-sphere's volume 2 pi^2, so a random rotation's angle has density (2/pi) sin^2", sp.simplify(vol - 2*sp.pi**2) == 0)
rng = np.random.default_rng(7); v = rng.normal(size=(400000, 4)); v /= np.linalg.norm(v, axis=1)[:, None]
ang = np.arccos(np.clip(v[:, 0], -1, 1))
h, e = np.histogram(ang, bins=24, range=(0, np.pi)); cdf = lambda t: (t - np.sin(t)*np.cos(t))/np.pi
exp_ = (cdf(e[1:]) - cdf(e[:-1]))*len(ang); chi = float(np.sum((h - exp_)**2/exp_))
rep("4d simulated: 400,000 random points of the 3-sphere give the sin^2 angle law", chi < 45, f"chi2 {chi:.1f} on 23 dof")
# 5. primes vs the 3-sphere and the circle
def angles(dat): return np.array([math.acos(max(-1, min(1, a_/(2*math.sqrt(p))))) for p, a_ in dat])
th_n, th_c = angles(noncm), angles(cmd)
def chi_st(t, bins=12):
    h, e = np.histogram(t, bins=bins, range=(0, math.pi)); ex = (cdf(e[1:]) - cdf(e[:-1]))*len(t)
    return float(np.sum((h - ex)**2/ex)), h.tolist(), ex.tolist(), e.tolist()
cn = chi_st(th_n)
rep("5a 4 non-CM curves, primes to 6000: turning angles fit the 3-sphere's sin^2 law", cn[0] < 30, f"chi2 {cn[0]:.1f} on 11 dof, {len(th_n)} primes")
half = np.mean(np.isclose(th_c, math.pi/2, atol=1e-12))
rest = th_c[~np.isclose(th_c, math.pi/2, atol=1e-12)]
hh, ee = np.histogram(rest, bins=12, range=(0, math.pi)); exu = np.full(12, len(rest)/12); chiu = float(np.sum((hh - exu)**2/exu))
rep("5b 2 CM curves: half the primes are exact quarter turns (a_p = 0), the rest spread evenly on a circle (the 1-sphere), not sin^2",
    abs(half - 0.5) < 0.03 and chiu < 30, f"{half:.3f} quarter turns; uniform chi2 {chiu:.1f} on 11 dof")
# 6. moments through the cube
xn = np.array([a_/math.sqrt(p) for p, a_ in noncm]); xc = np.array([a_/math.sqrt(p) for p, a_ in cmd])
cat = [1, 2, 5, 14, 42]; cen = [2, 6, 20, 70, 252]   # middle Pascal entries C(2n,n)
mn = [float(np.mean(xn**(2*n))) for n in range(1, 6)]; mc = [float(np.mean(xc**(2*n))) for n in range(1, 6)]
note("non-CM moments E[(a_p/sqrt p)^2n], n=1..5: " + ", ".join(f"{m:.2f}" for m in mn) + "   Catalan: 1, 2, 5, 14, 42")
note("CM moments: " + ", ".join(f"{m:.2f}" for m in mc) + "   half the middle Pascal entry: 1, 3, 10, 35, 126")
rep("6a non-CM moments sit near the Catalan numbers (within 10% for n <= 3, from a finite set of primes) and far from the circle's", all(abs(m/c - 1) < 0.10 for m, c in zip(mn[:3], cat)) and all(abs(m - c) < abs(m - d) for m, c, d in zip(mn, cat, cen)))
rep("6b CM moments sit near half the middle layer of the 2n-cube (within 10% for n <= 3) and far from Catalan (from n = 2; both laws give 1 at n = 1)", all(abs(m/(c/2) - 1) < 0.10 for m, c in zip(mc[:3], cen)) and all(abs(m - c/2) < abs(m - d) for m, c, d in zip(mc[1:], cen[1:], cat[1:])))
# exact: cube layers carry frequencies; averaging picks layers
ok = True
for n in range(1, 9):
    m = 2*n
    layers = [sp.binomial(m, j) for j in range(m + 1)]           # corners of the m-cube with j ones; weight e^{i(m-2j)theta}
    circle = sum(L for j, L in enumerate(layers) if m - 2*j == 0)
    sphere = sum(L*(1 if m - 2*j == 0 else (-sp.Rational(1, 2) if abs(m - 2*j) == 2 else 0)) for j, L in enumerate(layers))
    integ = sp.integrate((2*sp.cos(th))**m*2/sp.pi*sp.sin(th)**2, (th, 0, sp.pi))
    dyck = sum(1 for c in itertools.product([1, -1], repeat=m) if sum(c) == 0 and all(s >= 0 for s in itertools.accumulate(c)))
    ok &= circle == sp.binomial(m, n) and sphere == sp.catalan(n) and sp.simplify(integ - sphere) == 0 and dyck == sp.catalan(n)
rep("6c exact, n = 1..8: the circle keeps the middle layer of the 2n-cube, the 3-sphere keeps middle minus neighbour = Catalan = corner walks that never dip below the start", ok)
rep("6d the sin^2 weight is (2 - e^{2i theta} - e^{-2i theta})/4 = |alpha - beta|^2/4: the squared distance between the two walls",
    sp.simplify(sp.sin(th)**2 - (2 - sp.exp(2*sp.I*th) - sp.exp(-2*sp.I*th))/4) == 0)
R.update({"noncm_hist": cn[1], "noncm_exp": cn[2], "edges": cn[3], "cm_angles": th_c.tolist()[:4000], "noncm_x": xn.tolist(), "cm_x": xc.tolist(),
          "moments_noncm": mn, "moments_cm": mc})
print(f"\n{npass} of {ntot} checks pass"); out.append(f"\n{npass} of {ntot} checks pass")
open("s3_checks_output.txt", "w").write("\n".join(out) + "\n"); json.dump(R, open("s3_results.json", "w"))
