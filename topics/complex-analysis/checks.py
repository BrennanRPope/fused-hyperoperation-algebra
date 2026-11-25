"""Complex analysis through the circle operation: every claim on the page, checked.
Run: python3 checks.py > checks_output.txt
"""
import cmath, math, json
import numpy as np
import mpmath as mp
import sympy as sp

rng = np.random.default_rng(7)
results = []
def check(name, ok, detail=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("   " + detail) if detail else ""))

def ma(z, w): return z + w + z * w
def rc(n=1): return rng.normal(size=n) + 1j * rng.normal(size=n)

# ---------------- 1. The circle operation is a spin and a stretch about the wall
z, w = rc(200), rc(200)
check("lift: 1 + z*w = (1+z)(1+w) on 200 random complex pairs",
      np.max(abs(1 + ma(z, w) - (1 + z) * (1 + w))) < 1e-12)
check("the wall -1 is fixed by every map z -> z*w", np.max(abs(ma(-1, w) + 1)) < 1e-12)
ratio = (ma(z, w) + 1) / (z + 1)
check("z -> z*w multiplies distance from the wall by |1+w| and turns by arg(1+w)",
      np.max(abs(ratio - (1 + w))) < 1e-10)
# similar triangles (-1, 0, w) and (-1, z, z*w)
t1 = (w - (-1)) / (0 - (-1)); t2 = (ma(z, w) - (-1)) / (z - (-1))
check("triangle (-1, 0, w) is similar to triangle (-1, z, z*w) (same complex shape ratio)",
      np.max(abs(t1 - t2)) < 1e-10)
check("the step z -> z*w is the vector w turned and stretched by 1+z: z*w - z = (1+z)w",
      np.max(abs(ma(z, w) - z - (1 + z) * w)) < 1e-12)
# affine maps fixing -1 are exactly the circle operation by something
al, be = rc(50), None
check("an affine map z -> a z + b fixes -1 iff b = a - 1, i.e. it is z -> z*(a-1)",
      np.max(abs((al * z[:50] + (al - 1)) - ma(z[:50], al - 1))) < 1e-12)
u, v = rc(50), rc(50)
check("composing two such maps combines their parameters by the circle operation",
      np.max(abs(ma(ma(z[:50], u), v) - ma(z[:50], ma(u, v)))) < 1e-10)
k = rc(50)
check("dial: z *_k w = z + w + k z w spins about its own wall -1/k by 1 + k w",
      np.max(abs((z[:50] + w[:50] + k * z[:50] * w[:50] + 1 / k) / (z[:50] + 1 / k) - (1 + k * w[:50]))) < 1e-9)

# ---------------- 2. the circle-operation circle and roots of zero
th1, th2 = rng.uniform(0, 2 * np.pi, 100), rng.uniform(0, 2 * np.pi, 100)
c1, c2 = np.exp(1j * th1) - 1, np.exp(1j * th2) - 1
check("circle |1+z| = 1 (through 0, centred on the wall) is closed under the circle operation",
      np.max(abs(abs(1 + ma(c1, c2)) - 1)) < 1e-12)
check("on that circle the circle operation adds angles seen from the wall",
      np.max(abs(ma(c1, c2) - (np.exp(1j * (th1 + th2)) - 1))) < 1e-12)
def mpow(z, n):
    r = 0
    for _ in range(n): r = ma(r, z)
    return r
ok = True
for n in range(1, 13):
    roots = [cmath.exp(2j * math.pi * j / n) - 1 for j in range(n)]
    ok &= all(abs(mpow(r, n)) < 1e-10 for r in roots)
    ok &= all(abs(abs(1 + r) - 1) < 1e-12 for r in roots)
    # no other solutions: (1+z)^n - 1 has exactly n roots
    ok &= len(set(np.round(np.roots([1] + [0] * (n - 1) + [-1]) - 1, 8))) == n
check("the n circle n-th roots of 0 are zeta - 1, a regular n-gon on |1+z| = 1 (n = 1..12)", ok)
ok = True
for n in range(1, 13):
    tot = 0
    for j in range(n): tot = ma(tot, cmath.exp(2j * math.pi * j / n) - 1)
    ok &= abs(tot - ((-1) ** (n + 1) - 1)) < 1e-9
check("combining all n roots by the circle operation gives 0 for odd n and -2 for even n", ok)

# ---------------- 3. the log of the lift, the branch cut, winding about the wall
L = lambda z: np.log(1 + z)
d = (L(ma(z, w)) - L(z) - L(w)) / (2j * np.pi)
check("Log(1+z*w) = Log(1+z) + Log(1+w) + 2 pi i m with m in {-1, 0, 1}",
      np.max(abs(d - np.round(d))) < 1e-10 and set(np.round(d.real).astype(int)) <= {-1, 0, 1})
def wind(center, r, about, N=4000):
    t = np.linspace(0, 2 * np.pi, N + 1)
    path = center + r * np.exp(1j * t)
    lg = np.unwrap(np.angle(1 + path))
    return (lg[-1] - lg[0]) / (2 * np.pi)
check("following log(1+z) once round a loop about the wall adds 2 pi i", abs(wind(-1, 0.5, -1) - 1) < 1e-9)
check("a loop round 0 that misses the wall adds nothing (the branch point is the wall, not 0)",
      abs(wind(0, 0.5, -1)) < 1e-9)
check("a loop that encloses both 0 and the wall adds 2 pi i", abs(wind(-0.5, 1.2, -1) - 1) < 1e-9)
E = lambda z: np.expm1(z)
a, b = rc(100), rc(100)
check("E(w) = e^w - 1 turns + into the circle operation: E(a+b) = E(a) * E(b)",
      np.max(abs(E(a + b) - ma(E(a), E(b)))) < 1e-9)
check("E(2 pi i) = 0: one trip round the wall is the period", abs(E(2j * np.pi)) < 1e-14)
y = np.linspace(-3, 3, 50)
check("E sends each vertical line Re w = s to the circle of radius e^s about the wall",
      all(np.max(abs(abs(E(s + 1j * y) + 1) - math.exp(s))) < 1e-12 for s in (-1, 0, 0.7)))
xs = np.linspace(-2, 1, 50)
check("E sends each horizontal line Im w = t to the ray from the wall at angle t",
      all(np.max(abs(np.angle(E(xs + 1j * t) + 1) - t)) < 1e-12 for t in (0.3, 1.0, 2.5)))

# ---------------- 4. Moebius: landscape laws, walls are the fixed points
A_, B_, P, Q = sp.symbols('a b p q')
F = (A_ + B_ + P * A_ * B_) / (1 - Q * A_ * B_)
num = sp.factor(sp.numer(sp.together(F - B_)))
check("fixed points of b -> F_pq(a,b): numerator of F - b is a(1 + p b + q b^2), the walls, for every a",
      sp.simplify(num - A_ * (1 + P * B_ + Q * B_ ** 2)) == 0)
M = lambda x: sp.Matrix([[1 + P * x, x], [-Q * x, 1]])
Fab = sp.together(F)
check("matrices M(a) = [[1+pa, a], [-qa, 1]] compose: M(a)M(b) = (1 - qab) M(F(a,b))",
      sp.simplify(M(A_) * M(B_) - (1 - Q * A_ * B_) * M(Fab)) == sp.zeros(2, 2))
check("det M(a) = 1 + pa + qa^2: the translation by a degenerates exactly at a wall",
      sp.expand(M(A_).det() - (1 + P * A_ + Q * A_ ** 2)) == 0)
check("tr^2 - 4 det = a^2 (p^2 - 4q): the discriminant alone sets the type of every translation",
      sp.expand(M(A_).trace() ** 2 - 4 * M(A_).det() - A_ ** 2 * (P ** 2 - 4 * Q)) == 0)
W1, W2 = sp.symbols('w1 w2')
zeta = lambda x: (1 - x / W1) / (1 - x / W2)
Fw = F.subs({P: -(W1 + W2) / (W1 * W2), Q: 1 / (W1 * W2)})
check("two walls w1, w2: zeta(x) = (1 - x/w1)/(1 - x/w2) turns the law into multiplication",
      sp.simplify(zeta(Fw) - zeta(A_) * zeta(B_)) == 0)
ok = True
for _ in range(50):
    p, q = rc(1)[0], rc(1)[0]
    w1, w2 = np.roots([q, p, 1])
    aa, bb = rc(1)[0] * 0.3, rc(1)[0] * 0.3
    FF = (aa + bb + p * aa * bb) / (1 - q * aa * bb)
    zt = lambda x: (1 - x / w1) / (1 - x / w2)
    ok &= abs(zt(FF) - zt(aa) * zt(bb)) < 1e-8 * (1 + abs(zt(FF)))
check("the same holds for 50 random complex points (p, q) of the landscape", ok)
aa, bb = 0.37, -1.2
check("tan law (0,1): walls +-i, zeta = (1+ix)/(1-ix) = e^{2i arctan x} (Cayley)",
      abs((1 + 1j * aa) / (1 - 1j * aa) - cmath.exp(2j * math.atan(aa))) < 1e-14)
check("tanh law (0,-1): walls +-1, zeta = (1+x)/(1-x) = e^{2 artanh x} (Doppler factor squared)",
      abs((1 + aa) / (1 - aa) - math.exp(2 * math.atanh(aa))) < 1e-14)
cc = 0.8
Fd = lambda a, b: (a + b + 2 * cc * a * b) / (1 - cc * cc * a * b)
psi = lambda x: x / (1 + cc * x)
check("double wall (p^2 = 4q, p = 2c): psi(x) = x/(1+cx) turns the law into plain addition",
      abs(psi(Fd(aa, 0.21)) - psi(aa) - psi(0.21)) < 1e-12)

# ---------------- 5. the rung step e^z - 1: a parabolic fixed point
zz = sp.symbols('z')
f = sp.exp(zz) - 1
check("e^z - 1 fixes 0 with derivative 1 (parabolic) and second derivative 1",
      f.subs(zz, 0) == 0 and sp.diff(f, zz).subs(zz, 0) == 1 and sp.diff(f, zz, 2).subs(zz, 0) == 1)
cs = sp.symbols('c1:8')
Aser = lambda x: -2 / x + sp.Rational(1, 3) * sp.log(x) + sum(c * x ** (i + 1) for i, c in enumerate(cs))
ser = sp.expand(sp.series(Aser(f) - Aser(zz) - 1, zz, 0, 9).removeO())
eqs = [ser.coeff(zz, i) for i in range(0, 9)]
sol = sp.solve(eqs[2:], cs, dict=True)[0]
cvals = [sol[c] for c in cs]
check("Abel ruler alpha(z) = -2/z + (1/3) log z + ...: the 1/z and log terms are forced (orders z^0, z^1 vanish)",
      eqs[0] == 0 and eqs[1] == 0)
check("next coefficients are rational: -1/36, 1/540, 1/7776, -71/435456",
      cvals[:4] == [sp.Rational(-1, 36), sp.Rational(1, 540), sp.Rational(1, 7776), sp.Rational(-71, 435456)])
mp.mp.dps = 40
cm = [mp.mpf(sp.Rational(c).p) / sp.Rational(c).q for c in cvals]
def Aal(z): return -2 / z + mp.log(z) / 3 + sum(c * z ** (k + 1) for k, c in enumerate(cm))
zt = mp.mpc(-0.01, 0.004)
check("alpha(f(z)) - alpha(z) = 1 to 1e-17 at |z| ~ 0.01", abs(Aal(mp.expm1(zt)) - Aal(zt) - 1) < 1e-17)
def orbit(z, n, g):
    for _ in range(n): z = g(z)
    return z
zl = -0.8 + 0.5j
for _ in range(20000): zl = cmath.exp(zl) - 1
zr = 0.8 + 0.5j
for _ in range(20000): zr = cmath.log(1 + zr)
check("attracting petal: forward orbits of e^z - 1 from the left creep into 0 like -2/n",
      abs(zl * 20000 + 2) < 0.02, f"n*z_n = {zl*20000:.4f}")
check("repelling petal: orbits of the inverse log(1+z) creep into 0 from the right like 2/n",
      abs(zr * 20000 - 2) < 0.02, f"n*z_n = {zr*20000:.4f}")
def att(z, n=3000):
    for _ in range(n): z = mp.expm1(z)
    return Aal(z) - n
def rep(z, n=3000):
    for _ in range(n): z = mp.log1p(z)
    return Aal(z) + n
est = []
for z0 in [mp.mpc(0, 1), mp.mpc(-0.3, 1), mp.mpc(0.3, 0.8), mp.mpc(0, 0.5), mp.mpc(0.1, 0.3)]:
    wv = rep(z0); D = att(z0) - wv
    est.append((complex(D * mp.exp(-2j * mp.pi * wv)), float(wv.imag), float(abs(D))))
a1 = est[-1][0]
check("the two Abel rulers (left petal, right petal) disagree on the upper overlap: D = a1 e^{2 pi i w}, a1 != 0",
      abs(a1) > 60, f"a1 = {a1.real:.6f} + {a1.imag:.6f}i, |a1| = {abs(a1):.4f}")
check("a1 is the same from five starting points at heights Im w from 2.4 to 6.4 (rel. spread < 1e-4)",
      max(abs(e[0] - a1) for e in est) / abs(a1) < 1e-4)
rate = math.log(est[0][2] / est[3][2]) / (est[3][1] - est[0][1])
check("the ripple dies like e^{-2 pi Im w} (fitted rate vs 2 pi)", abs(rate / (2 * math.pi) - 1) < 1e-3, f"rate = {rate:.4f}")
z0 = mp.mpc(0, -1); wv = rep(z0); D = att(z0) - wv
check("lower overlap carries the mirror ripple conj(a1) e^{-2 pi i w} (e^z - 1 is real)",
      abs(complex(D * mp.exp(2j * mp.pi * wv)) - a1.conjugate()) < 1e-2)

# ---------------- 6. products: circle-operation folds and their wall hits
def sinfold(z, N):
    t = 0
    for n in range(1, N + 1): t = ma(t, -z * z / n ** 2)
    return t
zs = [0.3 + 0.4j, 1.7 - 0.2j, -2.5 + 1.1j]
ok = all(abs((1 + sinfold(z, 200000)) - cmath.sin(math.pi * z) / (math.pi * z)) < 5e-4 * abs(z) ** 2 for z in zs)
check("sin(pi z)/(pi z) - 1 is the circle-operation fold of the pieces -z^2/n^2 (Euler), checked at 3 complex z", ok)
errs = [abs((1 + sinfold(1.3, N)) - math.sin(math.pi * 1.3) / (math.pi * 1.3)) for N in (100, 1000, 10000)]
check("the fold converges like 1/N because sum |z^2/n^2| converges", errs[0] / errs[1] > 8 and errs[1] / errs[2] > 8,
      " ".join(f"{e:.2e}" for e in errs))
check("its zeros are exactly the wall hits: piece -z^2/n^2 = -1 iff z = +-n",
      all(abs(1 + sinfold(n, 50)) < 1e-15 for n in (1, 2, 3)) and abs(1 + sinfold(2.5, 50)) > 0.01)
def gfold(z, N):
    t = 0
    for n in range(1, N + 1): t = ma(t, (1 + z / n) * cmath.exp(-z / n) - 1)
    return t
ok = all(abs(z * cmath.exp(0.5772156649015329 * z) * (1 + gfold(z, 200000)) - 1 / complex(mp.gamma(z))) < 1e-4 for z in (0.5 + 0.5j, 2.2 - 1j, -1.5 + 0.3j))
check("1/Gamma(z) = z e^{gamma z} times the fold of Weierstrass pieces (1+z/n)e^{-z/n} - 1 (complex z)", ok)
check("each Weierstrass piece is O(z^2/n^2) (the corner is all that is left), so the fold converges",
      all(abs(((1 + 0.7 / n) * math.exp(-0.7 / n) - 1) / (0.7 / n) ** 2 + 0.5) < 0.6 / n for n in (10, 100, 1000)))
def fold(seq):
    t = 0
    for x in seq: t = ma(t, x)
    return t
p4 = 1 + fold([(-1) ** n / math.sqrt(n) for n in range(2, 10 ** 4)])
p6 = 1 + fold([(-1) ** n / math.sqrt(n) for n in range(2, 10 ** 6)])
check("sum a_n converging is not enough: a_n = (-1)^n/sqrt n has prod(1+a_n) -> 0 (the a^2/2 corners pile up)",
      p6 < p4 / 5, f"N=1e4: {p4:.4f}  N=1e6: {p6:.5f}")
q6 = 1 + fold([(-1) ** n / n for n in range(2, 10 ** 6 + 2)])
check("a_n = (-1)^n/n: sum of corners finite, fold converges to 0 (product exactly 1: the pairs cancel)", abs(q6 - 1) < 1e-5)

# argument principle: wall hits
def wall_count(fun, dfun, c, r, N=20000):
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    zc = c + r * np.exp(1j * t); dz = 1j * r * np.exp(1j * t) * (2 * np.pi / N)
    return np.sum(dfun(zc) / (1 + fun(zc)) * dz) / (2j * np.pi)
fs = lambda z: np.sin(np.pi * z) / (np.pi * z) - 1
dfs = lambda z: (np.pi * z * np.pi * np.cos(np.pi * z) - np.pi * np.sin(np.pi * z)) / (np.pi * z) ** 2
cnt = wall_count(fs, dfs, 0, 3.5)
check("(1/2 pi i) closed integral of f'/(1+f) counts wall hits: sine fold on |z| = 3.5 gives 6",
      abs(cnt - 6) < 1e-6, f"{cnt.real:.6f}")
fe = lambda z: np.exp(z) - 1
cnt2 = wall_count(lambda z: 0 * z + fe(z), np.exp, 0, 5)
check("the rung step e^z - 1 never hits the wall (count 0 on |z| = 5)", abs(cnt2) < 1e-9)
fp = lambda z: (1 + z) ** 5 - 1
check("z^(*5) = (1+z)^5 - 1 hits the wall 5 times, all at z = -1 (count 5 on |z+1| = 1)",
      abs(wall_count(fp, lambda z: 5 * (1 + z) ** 4, -1, 1) - 5) < 1e-9)

# ---------------- 7. Tate: glue (1+a) ~ q(1+a)
qq = 0.2 * cmath.exp(0.7j)
a, b = rc(20), rc(20)
a2 = qq * (1 + a) - 1
check("the gluing (1+a) ~ q(1+a) respects the circle operation: a ~ a' implies a*b ~ a'*b (torus of the elliptic page)",
      np.max(abs((1 + ma(a2, b)) - qq * (1 + ma(a, b)))) < 1e-12)

n = len(results); p = sum(results)
print(f"\n{p} of {n} checks pass")
json.dump({"a1": [a1.real, a1.imag], "cvals": [str(c) for c in cvals],
           "est": [[e[0].real, e[0].imag, e[1], e[2]] for e in est]}, open("results.json", "w"), indent=1)
