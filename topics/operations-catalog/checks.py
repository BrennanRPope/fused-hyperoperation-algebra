"""Numerical checks for the operations catalog (math-theory/operations-catalog).
Each check prints PASS/FAIL. Plain Python + mpmath only."""
import random, math, itertools
from mpmath import mp, mpf, log, exp, quad, matrix, det, eye
mp.dps = 30
random.seed(7)
R = lambda lo=-0.4, hi=0.6: mpf(random.uniform(lo, hi))
res = []
def check(name, ok, info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ") + name + (("   " + info) if info else ""))
tol = mpf("1e-20")

def circ(k, a, b): return a + b + k*a*b
def lens(k, d, a, b): return (a + b + k*a*b) / (1 - d*a*b)

# 1 lift: 1+k(a o_k b) = (1+ka)(1+kb)
worst = max(abs((1 + k*circ(k, a, b)) - (1 + k*a)*(1 + k*b)) for k, a, b in [(R(-2, 2), R(), R()) for _ in range(200)])
check("lift turns o_k into multiplication", worst < tol, f"worst {float(worst):.1e}")
# 2 wall absorbing
worst = max(abs(circ(k, -1/k, b) - (-1/k)) for k, b in [(R(0.2, 2), R(-3, 3)) for _ in range(200)])
check("wall -1/k is the absorbing element: wall o_k b = wall", worst < tol)
# 3 quasi-inverse
qi = lambda k, a: -a/(1 + k*a)
worst = max(abs(circ(k, a, qi(k, a))) for k, a in [(R(-2, 2), R()) for _ in range(200)])
check("way back -a/(1+ka) is the quasi-inverse", worst < tol)
# 4 circle power
def cpow(k, a, n):
    x = mpf(0)
    for _ in range(n): x = circ(k, x, a)
    return x
worst = max(abs(cpow(k, a, n) - ((1 + k*a)**n - 1)/k) for k, a, n in [(R(0.1, 2), R(), random.randint(1, 9)) for _ in range(200)])
check("circle power a^(o n) = ((1+ka)^n - 1)/k", worst < mpf("1e-18"))
k, a = mpf("0.7"), mpf("0.3")
half = ((1 + k*a)**mpf("0.5") - 1)/k
check("real circle powers: half power twice = one", abs(circ(k, half, half) - a) < tol)
# 5 cube combine = sum k^(m-1) e_m ; det(I+kA) = 1 + k D_k(A)
def esym(xs, m): return sum(math.prod(c) for c in itertools.combinations(xs, m)) if m else 1
worst = mpf(0)
for _ in range(100):
    k = R(-2, 2); xs = [R() for _ in range(random.randint(1, 6))]
    f = mpf(0)
    for x in xs: f = circ(k, f, x)
    s = sum(k**(m - 1)*esym(xs, m) for m in range(1, len(xs) + 1))
    worst = max(worst, abs(f - s))
check("n-fold combine = sum over cube corners k^(m-1) e_m", worst < mpf("1e-18"))
A = matrix([[R() for _ in range(3)] for _ in range(3)]); k = mpf("0.8")
e1 = A[0, 0] + A[1, 1] + A[2, 2]
e2 = sum(A[i, i]*A[j, j] - A[i, j]*A[j, i] for i, j in [(0, 1), (0, 2), (1, 2)])
Dk = e1 + k*e2 + k*k*det(A)
check("matrix form: det(I + kA) = 1 + k(e1 + k e2 + k^2 e3)", abs(det(eye(3) + k*A) - (1 + k*Dk)) < tol)
# 6 Mobius inversion on the cube (n=3): recover a single corner weight from combines of subsets
# f(S) = combine of the set S; g(S) = sum_{T subset S} (-1)^{|S|-|T|} f(T) gives the top-corner term k^(|S|-1) prod S for |S|>=2
k = mpf("0.9"); xs = [R() for _ in range(3)]
def comb(sub):
    f = mpf(0)
    for x in sub: f = circ(k, f, x)
    return f
idx = range(3)
g = sum((-1)**(3 - r)*comb([xs[i] for i in T]) for r in range(4) for T in itertools.combinations(idx, r))
check("Mobius inversion over the cube isolates the top corner k^2 abc", abs(g - k**2*math.prod(xs)) < tol)
# 7 flip x/(1+kx): carries o_k to o_{-k}; flip = -(way back)
fl = lambda k, x: x/(1 + k*x)
worst = max(abs(fl(k, circ(k, a, b)) - circ(-k, fl(k, a), fl(k, b))) for k, a, b in [(R(-2, 2), R(0, 0.6), R(0, 0.6)) for _ in range(200)])
check("flip x/(1+kx) carries o_k to o_-k", worst < tol)
worst = max(abs(fl(k, a) + qi(k, a)) for k, a in [(R(-2, 2), R()) for _ in range(100)])
check("flip = minus the way back", worst < tol)
check("odds/probability is the flip at k = 1 / k = -1", abs(fl(1, mpf(3)) - mpf("0.75")) < tol and abs(fl(-1, mpf("0.75")) - 3) < tol)
# 8 mirror
worst = max(abs(circ(-k, a, b) + circ(k, -a, -b)) for k, a, b in [(R(-2, 2), R(), R()) for _ in range(100)])
check("mirror: a o_-k b = -((-a) o_k (-b))", worst < tol)
# 9 logarithm (bend): L_k = log(1+kx)/k turns o_k into +; general L_{k,d} = int dt/(1+kt+dt^2) turns the two-dial law into +
L = lambda k, x: log(1 + k*x)/k
worst = max(abs(L(k, circ(k, a, b)) - L(k, a) - L(k, b)) for k, a, b in [(R(0.1, 2), R(), R()) for _ in range(200)])
check("L_k(a o_k b) = L_k(a) + L_k(b)", worst < tol)
Lkd = lambda k, d, x: quad(lambda t: 1/(1 + k*t + d*t*t), [0, x])
worst = mpf(0)
for k, d in [(1, 0), (0, 1), (0, -1), (mpf("0.5"), mpf("0.3")), (mpf("1.2"), mpf("-0.2")), (2, 1)]:
    a, b = R(0, 0.3), R(0, 0.3)
    worst = max(worst, abs(Lkd(k, d, lens(k, d, a, b)) - Lkd(k, d, a) - Lkd(k, d, b)))
check("L_{k,d} = int_0^x dt/(1+kt+dt^2) turns the two-dial law into + (all three classes)", worst < mpf("1e-15"), f"worst {float(worst):.1e}")
# 10 k-log (Box-Cox) and k-exp (Tsallis): BoxCox_k(xy) = BoxCox_k(x) o_k BoxCox_k(y); L_k = log o exp_k
bc = lambda k, x: (x**k - 1)/k
worst = max(abs(bc(k, x*y) - circ(k, bc(k, x), bc(k, y))) for k, x, y in [(R(-2, 2), R(0.2, 3), R(0.2, 3)) for _ in range(200)])
check("Box-Cox k-log carries x to o_k", worst < tol)
expk = lambda k, x: (1 + k*x)**(1/k)
worst = max(abs(L(k, x) - log(expk(k, x))) for k, x in [(R(0.1, 2), R()) for _ in range(100)])
check("L_k = log of the Tsallis k-exponential", worst < tol)
# 11 transport: f_*(*)(a,b) = f(f^-1 a * f^-1 b)
def transport(f, finv, op): return lambda a, b: f(op(finv(a), finv(b)))
# slide psi_c(y) = y/(1+cy) ; check F_{k,d}(psi a, psi b) = psi(F_{k+2c, d+kc+c^2}(a, b))
worst = mpf(0)
for _ in range(200):
    k, d, c = R(-1, 1), R(-1, 1), R(-1, 1); a, b = R(-0.2, 0.2), R(-0.2, 0.2)
    psi = lambda y: y/(1 + c*y)
    worst = max(worst, abs(lens(k, d, psi(a), psi(b)) - psi(lens(k + 2*c, d + k*c + c*c, a, b))))
check("slide y/(1+cy): (k,d) -> (k+2c, d+kc+c^2)", worst < tol)
check("slide keeps the discriminant k^2 - 4d", all(abs((k + 2*c)**2 - 4*(d + k*c + c*c) - (k*k - 4*d)) < tol for k, d, c in [(R(-1, 1), R(-1, 1), R(-1, 1)) for _ in range(50)]))
# rescale x -> lam x : lam^-1 F(lam a, lam b) = F_{lam k, lam^2 d}(a,b)
worst = max(abs(lens(k, d, l*a, l*b)/l - lens(l*k, l*l*d, a, b)) for k, d, l, a, b in [(R(-1, 1), R(-1, 1), R(0.2, 3), R(-0.2, 0.2), R(-0.2, 0.2)) for _ in range(200)])
check("rescale x -> lam x: (k,d) -> (lam k, lam^2 d)", worst < tol)
check("o_k is o_1 rescaled: a o_k b = (ka o_1 kb)/k", abs(circ(mpf("0.37"), a, b) - circ(1, mpf("0.37")*a, mpf("0.37")*b)/mpf("0.37")) < tol)
# three classes: slide c=-k/2 then rescale lands on (0,-1), (0,0) or (0,1)
def normal_form(k, d):
    c = -k/2; d1 = d + k*c + c*c   # = d - k^2/4 = -disc/4
    if abs(d1) < mpf("1e-25"): return (0, 0), c, 1
    l = 1/mp.sqrt(abs(d1)); return (0, 1 if d1 > 0 else -1), c, l
ok = True
for k, d in [(1, 0), (2, 0), (mpf("0.5"), mpf("0.3")), (2, 1), (0, mpf("-0.25")), (1, mpf("0.25"))]:
    (k2, d2), c, l = normal_form(mpf(k), mpf(d))
    a, b = mpf("0.11"), mpf("-0.07")
    # law after slide then rescale should equal lens(k2,d2)
    psi = lambda y: y/(1 + c*y)
    # lens(k,d) at psi-chart equals lens(k+2c, d1); then rescale by l: (0, l^2 d1)
    lhs = lens(k + 2*c, d + k*c + c*c, l*a, l*b)/l
    ok &= abs(lhs - lens(k2, d2, a, b)) < tol
    print(f"   (k,d)=({float(k)},{float(d)}) disc={float(k*k-4*d):+.2f} -> class {('velocity add (a+b)/(1+ab)' if d2==-1 else 'plain +' if d2==0 else 'tangent add (a+b)/(1-ab)')}")
check("every two-dial law slides+rescales to exactly one of +, (a+b)/(1+ab), (a+b)/(1-ab)", ok)
# velocity addition is itself o_k in another chart: (0,-1) disc 4 ~ (2,0) disc 4
c = mpf(-1)   # slide (2,0) by c=-1 -> (0, 0-2+1) = (0,-1)
psi = lambda y: y/(1 - y)
a, b = mpf("0.2"), mpf("0.1")
check("velocity addition is o_2 seen through the slide y/(1-y)", abs(circ(2, psi(a), psi(b)) - psi(lens(0, -1, a, b))) < tol)
# 12 Maslov dequantization: h log(e^{a/h} + e^{b/h}) -> max(a,b)
vals = [float(h*log(exp(mpf(2)/h) + exp(mpf(1)/h))) for h in (1, mpf("0.1"), mpf("0.01"))]
check("rescaling limit h*log(e^(a/h)+e^(b/h)) -> max (Maslov dequantization)", abs(vals[-1] - 2) < 1e-6, str([round(v, 6) for v in vals]))
# 13 climb a rung = transport by E = e^x - 1, at whole heights
E = lambda x: exp(x) - 1; Linv = lambda x: log(1 + x)
B0 = lambda a, b: a + b
B1 = transport(E, Linv, B0)
B2 = transport(E, Linv, B1)
a, b = mpf("1.3"), mpf("0.7")
check("climb: E_*(+) = a+b+ab", abs(B1(a, b) - circ(1, a, b)) < tol)
check("climb twice: E_*(E_*(+)) = (1+a)(1+b)e^{L(a)L(b)} - 1", abs(B2(a, b) - ((1 + a)*(1 + b)*exp(Linv(a)*Linv(b)) - 1)) < tol)
Bm1 = transport(Linv, E, B0)
check("climb down: L_*(+) = log(e^a + e^b - 1)", abs(Bm1(a, b) - log(exp(a) + exp(b) - 1)) < tol)
# 14 fuse = transport by the unit shift S^-1 (x -> x-1): + -> a+b+1, x -> o_1, a^b -> (1+a)^(1+b) - 1
S, Sinv = (lambda x: x + 1), (lambda x: x - 1)
check("fuse +: a+b+1", abs(transport(Sinv, S, lambda x, y: x + y)(a, b) - (a + b + 1)) < tol)
check("fuse x: a o_1 b", abs(transport(Sinv, S, lambda x, y: x*y)(a, b) - circ(1, a, b)) < tol)
check("fuse power: (1+a)^(1+b) - 1 (Rungs.fexp)", abs(transport(Sinv, S, lambda x, y: x**y)(a, b) - ((1 + a)**(1 + b) - 1)) < tol)
check("E = S^-1 o exp, so a climb is Bennett's exp-climb followed by a fuse", abs(E(a) - Sinv(exp(a))) < tol)
# 15 translation T_a(x) = a o_k x is affine, fixes the wall, composes by o_k
k = mpf("0.6"); T = lambda a: (lambda x: circ(k, a, x))
x = mpf("0.25")
check("T_a o T_b = T_(a o_k b)", abs(T(a)(T(b)(x)) - T(circ(k, a, b))(x)) < tol)
check("every T_a fixes the wall -1/k", abs(T(a)(-1/k) + 1/k) < tol)
# 16 order gap in the stretch-and-shift group: (s,c): x -> (1+s)x + c
def app(m, x): return (1 + m[0])*x + m[1]
m1, m2 = (mpf("0.3"), mpf("0.5")), (mpf("-0.2"), mpf("0.9"))
gap = app(m1, app(m2, x)) - app(m2, app(m1, x))
check("order gap of two stretch-and-shifts = s1*c2 - s2*c1 (a wedge)", abs(gap - (m1[0]*m2[1] - m2[0]*m1[1])) < tol)
print(f"\n{sum(res)} of {len(res)} checks pass")
