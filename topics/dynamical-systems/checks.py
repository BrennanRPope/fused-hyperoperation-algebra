# Dynamical systems through mult-add: every claim on the page, checked.
import sympy as sp, numpy as np, math, cmath, mpmath as mp, json
ok = 0; tot = 0
def check(name, cond):
    global ok, tot
    tot += 1; ok += bool(cond); print(('PASS ' if cond else 'FAIL ') + name)
z = lambda e: sp.simplify(e) == 0
x, y, w, h, k, a, b, t, s, u, zz = sp.symbols('x y w h k a b t s u z')
ma = lambda p, q: p + q + p*q
wa = lambda kk, p, q: p + q + kk*p*q
msq = lambda v: ma(v, v)
inv = lambda v: -v/(1 + v)

print("--- 1. The lift is a change of coordinates: mult-add powers are z^n ---")
check("x ✳ x = (1+x)² − 1: the lift z = 1 + x conjugates it to z ↦ z²", z(1 + msq(x) - (1 + x)**2))
check("x ↦ x^(✳n) = (1+x)^n − 1 is conjugate to z ↦ z^n, n = 2..6, and composes: (x^(✳m))^(✳n) = x^(✳mn)",
      all(z(((1 + ((1 + x)**m - 1))**n - 1) - ((1 + x)**(m*n) - 1)) for m in range(2, 5) for n in range(2, 5)))
f = msq(x); itn = x
for n in range(1, 5): itn = sp.expand(msq(itn))
check("4 steps of x ✳ x: 1 + f⁴(x) = (1+x)^16", z(1 + itn - (1 + x)**16))
check("every dial square x ✳_k x is mult-add square in other units: k·(x ✳_k x) = (kx) ✳ (kx)", z(k*wa(k, x, x) - msq(k*x)))
check("so x ✳_k x (k ≠ 0) are all one map up to rescaling; k = 0 gives the linear map x ↦ 2x", z(wa(0, x, x) - 2*x))
check("k = −2 (parity): x ✳ x = 2x(1 − x), the logistic map at r = 2", z(wa(-2, x, x) - 2*x*(1 - x)))
check("k = −1 (failure union): x ✳ x = 1 − (1 − x)²", z(wa(-1, x, x) - (1 - (1 - x)**2)))

print("--- 2. Fixed points: the wall is superattracting, the identity repels ---")
check("wall: (−1 + h) ✳ (−1 + h) = −1 + h²  (multiplier 0: superattracting)", z(msq(-1 + h) - (-1 + h**2)))
check("identity: h ✳ h = 2h + h²  (multiplier 2: repelling)", z(msq(h) - (2*h + h**2)))
check("dial wall −1/k: (−1/k + h) ✳_k (−1/k + h) = −1/k + k·h²", z(wa(k, -1/k + h, -1/k + h) - (-1/k + k*h**2)))
check("logistic r = 2 has its superstable point at 1/2 = −1/k for k = −2", z(sp.diff(2*x*(1 - x), x).subs(x, sp.Rational(1, 2))))
check("x^(✳n): multiplier n at 0, multiplier 0 at the wall (n = 2..6)",
      all(sp.diff((1 + x)**n - 1, x).subs(x, 0) == n and sp.diff((1 + x)**n - 1, x).subs(x, -1) == 0 for n in range(2, 7)))
cnt = []
for n in range(1, 6):
    roots = [cmath.exp(2j*math.pi*j/(2**n - 1)) - 1 for j in range(2**n - 1)]
    good = all(abs((1 + r)**(2**n) - (1 + r)) < 1e-9 for r in roots)
    cnt.append(good and len(set(round(r.real, 9) + 1j*round(r.imag, 9) for r in roots)) == 2**n - 1)
check("points of period dividing n (besides the wall): 1 + x a (2ⁿ − 1)-th root of unity, 1, 3, 7, 15, 31 of them on |1+x| = 1", all(cnt))
rng = np.random.default_rng(3); good = True
for _ in range(300):
    r0 = rng.uniform(0.2, 1.8); th = rng.uniform(0, 2*math.pi); v = r0*cmath.exp(1j*th) - 1
    for _ in range(12):
        v = v + v + v*v
        if abs(v) > 1e6: break
    good &= (abs(v + 1) < 1e-6) if r0 < 0.95 else (abs(v) > 1e3) if r0 > 1.05 else True
check("inside the circle |1+x| = 1 orbits fall into the wall −1; outside they escape (300 random starts)", good)
v = cmath.exp(0.7j) - 1; mx = 0
for _ in range(40): v = v + v + v*v; v = (v + 1)/abs(v + 1) - 1; mx = max(mx, abs(abs(v + 1) - 1))
check("on the circle (the Julia set) the map is angle doubling: arg(1 + x✳x) = 2·arg(1 + x)",
      all(abs(cmath.phase((1 + msq(cmath.exp(1j*q) - 1))/(cmath.exp(1j*q))**2)) < 1e-12 for q in np.linspace(0, 6, 25)))

print("--- 3. Two rulers: Koenigs (multiplier 2) and Fatou (multiplier 1) ---")
check("Koenigs coordinate of x ✳ x at 0 is the log ruler: log(1 + x✳x) = 2·log(1 + x)", z(sp.expand_log(sp.log((1 + x)**2), force=True) - 2*sp.log(1 + x)))
mp.mp.dps = 40
def koenigs_limit(x0, n=60):
    v = mp.mpf(x0)
    for _ in range(n): v = mp.sqrt(1 + v) - 1          # inverse branch of x ✳ x
    return v*2**n
check("Koenigs's limit 2ⁿ·f⁻ⁿ(x) → log(1+x) (x = 0.5, 3, −0.6)", all(abs(koenigs_limit(x0) - mp.log(1 + x0)) < 1e-15 for x0 in (0.5, 3, -0.6)))
check("x^(✳n) for every n shares the same ruler: log(1 + x^(✳n)) = n·log(1+x)", all(z(sp.expand_log(sp.log((1 + x)**n), force=True) - n*sp.log(1 + x)) for n in range(2, 7)))
check("E(x) = e^x − 1 has multiplier 1 at 0 (parabolic)", sp.diff(sp.exp(x) - 1, x).subs(x, 0) == 1)
check("no Koenigs ruler for E: φ(E(x)) = φ(x) with φ = x + c₂x² + … fails at order x² (coefficient 1/2)",
      sp.expand(sp.series((lambda F: F.subs(x, sp.exp(x) - 1) - F)(x + sp.Symbol('c2')*x**2 + sp.Symbol('c3')*x**3), x, 0, 3).removeO()).coeff(x, 2) == sp.Rational(1, 2))
check("E turns + into ✳: E(a + b) = E(a) ✳ E(b), so x ✳ x = E(2·E⁻¹(x)): mult-add squaring is doubling seen through E",
      z(sp.exp(a + b) - 1 - ma(sp.exp(a) - 1, sp.exp(b) - 1)) and z(sp.exp(2*sp.log(1 + x)) - 1 - msq(x)))

print("--- 4. The rung flow: a vector field, a parabolic point, the Abel ruler ---")
N = 10
E = sp.series(sp.exp(x) - 1, x, 0, N + 1).removeO()
its = [x]
for n in range(1, N + 1): its.append(sp.series(E.subs(x, its[-1]), x, 0, N + 1).removeO())
coef = {i: sp.expand(sp.interpolate([(n, sp.expand(its[n]).coeff(x, i)) for n in range(i)], s)) for i in range(2, N + 1)}
vser = sum(sp.diff(coef[i], s).subs(s, 0)*x**i for i in range(2, N + 1))
check("vector field of the rung flow (d/ds of E^s at s=0): v = x²/2 − x³/12 + x⁴/48 − x⁵/180 + 11x⁶/8640 − …",
      z(sp.expand(vser).coeff(x, 2) - sp.Rational(1, 2)) and [sp.expand(vser).coeff(x, i) for i in range(2, 7)] ==
      [sp.Rational(1, 2), -sp.Rational(1, 12), sp.Rational(1, 48), -sp.Rational(1, 180), sp.Rational(11, 8640)])
check("E is the time-1 map of v: v(E(x)) = E′(x)·v(x) to order x¹⁰", sp.series(vser.subs(x, E) - sp.diff(E, x)*vser, x, 0, N + 1).removeO() == 0)
check("v ≥ 0 near 0: the flow moves right on both sides (negatives creep up to 0, positives leave)", sp.expand(vser).coeff(x, 2) > 0)
iv = sp.series(1/vser, x, 0, 4).removeO()
check("Abel ruler: α′ = 1/v = 2/x² + 1/(3x) − 1/36 + x/270 + …, so α = −2/x + ⅓·log|x| − x/36 + x²/540 + …",
      z(iv - (2/x**2 + 1/(3*x) - sp.Rational(1, 36) + x/270 + x**2/2592 - sp.Rational(71, 108864)*x**3)))
ac = [sp.Rational(-1, 36), sp.Rational(1, 540), sp.Rational(1, 7776), sp.Rational(-71, 435456)]
def alpha_asym(y):
    return -2/y + mp.log(abs(y))/3 + sum(float(c)*y**(j + 1) for j, c in enumerate(ac))
def alpha(x0, n=4000):
    y = mp.mpf(x0)
    if y > 0:
        for _ in range(n): y = mp.log(1 + y)
        return alpha_asym(y) + n
    for _ in range(n): y = mp.exp(y) - 1
    return alpha_asym(y) - n
mp.mp.dps = 30
good = all(abs(alpha(mp.exp(x0) - 1) - alpha(x0) - 1) < 1e-10 for x0 in (0.3, 1.0, 2.5, -0.2, -0.7, -2.0))
check("α(E(x)) = α(x) + 1 on both sides of 0 (x = 0.3, 1, 2.5, −0.2, −0.7, −2): one rung is one unit of the ruler", good)
def alpha_inv(A0, lo, hi):
    return mp.findroot(lambda yy: alpha(yy) - A0, (lo, hi), solver='bisect' if False else 'anderson')
half = alpha_inv(alpha(1.0) + 0.5, 1.05, 1.5)
check("half a rung by the ruler: E^(1/2)(1) = α⁻¹(α(1) + ½) = 1.2710274139 (fractional-rungs page)", abs(half - mp.mpf('1.2710274138906978')) < 1e-9)
mid = lambda xx: (alpha(xx + mp.mpf('1e-8')) - alpha(xx - mp.mpf('1e-8')))/mp.mpf('2e-8')
check("v = 1/α′ agrees with the series near 0 (x = ±0.05)", all(abs(1/mid(xx) - sum(float(sp.expand(vser).coeff(x, i))*xx**i for i in range(2, N + 1))) < 1e-9 for xx in (0.05, -0.05)))
fl = lambda tt, xx: xx/(1 - tt*xx/2)
check("near 0 the rung flow is the Möbius flow x/(1 − t·x/2) of x′ = x²/2 (first term of v), ruler −2/x adds t",
      z(fl(s, fl(t, x)) - fl(s + t, x)) and z(-2/fl(t, x) - (-2/x + t)))

print("--- 5. The logistic map is mult-add squaring folded by the mult-add inverse ---")
fold = lambda v: v*inv(v)/4
check("fold: X(w) = w·ι(w)/4 with ι(w) = −w/(1+w) the mult-add inverse; X(w) = −w²/(4(1+w))", z(fold(w) + w**2/(4*(1 + w))))
check("w ✳ ι(w) = 0, so w + ι(w) = −w·ι(w): X = −(w + ι(w))/4", z(ma(w, inv(w))) and z(fold(w) + (w + inv(w))/4))
check("2-to-1: X(ι(w)) = X(w)", z(fold(inv(w)) - fold(w)))
check("conjugacy: 4X(1 − X) = X(w ✳ w)  (logistic r = 4 is mult-add squaring folded)", z(4*fold(w)*(1 - fold(w)) - fold(msq(w))))
check("the Julia circle w = e^(iθ) − 1 folds onto [0, 1]: X = sin²(θ/2)",
      all(abs(complex(fold(sp.exp(sp.I*q) - 1).evalf()) - math.sin(q/2)**2) < 1e-12 for q in np.linspace(0.1, 6.2, 9)))
check("points off the circle and off the real line fold off the real line (sampled inside and outside)",
      all(abs(complex(fold(cmath.rect(r_, q) - 1)).imag) > 1e-9 for r_ in (0.5, 0.8, 1.3, 1.7) for q in (0.4, 1.3, 2.9, 4.0)))
check("the tent map is angle doubling folded: with X = sin²(πθ/2), θ ↦ 2θ (mod 2, folded) gives the logistic step",
      all(abs(4*math.sin(math.pi*q/2)**2*(1 - math.sin(math.pi*q/2)**2) - math.sin(math.pi*(1 - abs(1 - 2*q)) / 2)**2) < 1e-12 for q in np.linspace(0, 1, 21)))
rng = np.random.default_rng(7); xs = rng.uniform(0.01, 0.99, 2000); lyap = 0.0; xv = 0.3141; n_ = 200000
for _ in range(n_): lyap += math.log(abs(4 - 8*xv)); xv = 4*xv*(1 - xv)
check("Lyapunov exponent of the logistic map at r = 4 is log 2 (orbit average, 200 000 steps)", abs(lyap/n_ - math.log(2)) < 0.01)
check("Lyapunov exponent of x ✳ x on its Julia circle is log 2 exactly: |f′(x)| = 2|1 + x| = 2", all(abs(abs(2*(1 + (cmath.exp(1j*q) - 1))) - 2) < 1e-12 for q in np.linspace(0, 6, 13)))
hist, edges = np.histogram([4*math.sin(q)**2*0 + math.sin(q)**2 for q in rng.uniform(0, math.pi, 400000)], bins=20, range=(0, 1), density=True)
dens = [1/(math.pi*math.sqrt(m*(1 - m))) for m in (edges[:-1] + edges[1:])/2]
check("uniform angles on the circle fold to the logistic invariant density 1/(π√(x(1−x)))", max(abs(hh/dd - 1) for hh, dd in zip(hist[2:-2], dens[2:-2])) < 0.03)

print("--- 6. Chebyshev: z^n through the Joukowski map ---")
T = [sp.Integer(1), u]
for n in range(2, 9): T.append(sp.expand(2*u*T[-1] - T[-2]))
check("T_n(cos θ) = cos nθ, n ≤ 8", all(abs(float(T[n].subs(u, math.cos(0.37))) - math.cos(0.37*n)) < 1e-12 for n in range(9)))
check("2·T_n((z + 1/z)/2) = zⁿ + z⁻ⁿ  (Chebyshev is zⁿ seen through Joukowski)", all(z(2*T[n].subs(u, (zz + 1/zz)/2) - zz**n - zz**-n) for n in range(9)))
check("T_m ∘ T_n = T_mn (m, n ≤ 4)", all(z(T[m].subs(u, T[n]) - T[m*n]) for m in range(1, 5) for n in range(1, 5) if m*n <= 8))
check("logistic r = 4 in u = 1 − 2x is T₂: 1 − 2·4x(1−x) = T₂(1 − 2x)", z(1 - 2*4*x*(1 - x) - T[2].subs(u, 1 - 2*x)))
check("Joukowski in mult-add terms: (z + 1/z)/2 = 1 + (w + ι(w))/2 with z = 1 + w", z((1 + w + 1/(1 + w))/2 - (1 + (w + inv(w))/2)))
check("chain rule as mult-add: (f∘g)′ − 1 = (f′(g) − 1) ✳ (g′ − 1)  (a relabelling of the product rule)",
      z((sp.diff(sp.sin(sp.exp(x)), x) - 1) - ma(sp.cos(sp.exp(x)) - 1, sp.exp(x) - 1)))

print("--- 7. Flows: mult-add is 'run for time a, then time b' of y′ = 1 + y ---")
ysol = (1 + a)*sp.exp(t) - 1
check("y′ = 1 + y from y(0) = a gives y(t) = a ✳ E(t); the wall −1 is the equilibrium", z(sp.diff(ysol, t) - (1 + ysol)) and z(ysol - ma(a, sp.exp(t) - 1)))
yk = a + (sp.exp(k*t) - 1)/k + k*a*(sp.exp(k*t) - 1)/k
check("dial: y′ = 1 + k·y gives y(t) = a ✳_k (e^(kt) − 1)/k, equilibrium at the wall −1/k", z(sp.diff(yk, t) - (1 + k*yk)) and z(yk.subs(t, 0) - a))
# Riccati y' = 1 + p y + q y^2 from 0, numerically: compose times
def riccati(t_end, y0, p, q):
    from scipy.integrate import solve_ivp
    r = solve_ivp(lambda tt, yy: [1 + p*yy[0] + q*yy[0]**2], (0, t_end), [y0], rtol=1e-12, atol=1e-14); return r.y[0, -1]
pf, qf = 0.4, -0.3
Ya, Yb, Yab = riccati(0.7, 0, pf, qf), riccati(0.5, 0, pf, qf), riccati(1.2, 0, pf, qf)
check("landscape law F_{p,q} is the time law of y′ = 1 + p·y + q·y² (p = 0.4, q = −0.3: Y(0.7 + 0.5) = F(Y(0.7), Y(0.5)))",
      abs(Yab - (Ya + Yb + pf*Ya*Yb)/(1 - qf*Ya*Yb)) < 1e-9)

print("--- 8. The two half-rulers disagree by a periodic wobble (Ecalle's horn map) ---")
mp.mp.dps = 34
acm = [mp.mpf(1)*c for c in ac]
def aa(yv): return -2/yv + mp.log(yv)/3 + sum(c*yv**(j + 1) for j, c in enumerate(acm))
def aplus(x0, n=6000):
    yv = mp.mpc(x0)
    for _ in range(n): yv = mp.log(1 + yv)
    return aa(yv) + n
def aminus(x0, n=6000):
    yv = mp.mpc(x0)
    for _ in range(n): yv = mp.exp(yv) - 1
    return aa(yv) - n
cs = []; horn = []
for x0 in (0.45j, 0.5j, 0.1 + 0.5j, -0.1 + 0.5j):
    Ap = aplus(x0); D = Ap - aminus(x0); cs.append(D*mp.exp(-2j*mp.pi*Ap))
check("in the upper half-plane, α₊ − α₋ = c₁·e^(2πiα) with one constant c₁ at four points (spread < 1e-7)",
      max(abs(c_ - cs[0]) for c_ in cs) < 1e-7)
check("c₁ ≈ −10.748807 − 63.246298i, |c₁| ≈ 64.15318 (ruler normalized with no constant term)",
      abs(cs[0] - mp.mpc(-10.748807093, -63.2462981528)) < 1e-6 and abs(abs(cs[0]) - 64.1531845) < 1e-6)
Dl = aplus(-0.5j) - aminus(-0.5j)
check("lower half-plane is the mirror: c₋₁ = conj(c₁)", abs(Dl*mp.exp(2j*mp.pi*aplus(-0.5j)) - mp.conj(cs[0])) < 1e-6)
for sv in [0.32 + 0.06*i for i in range(15)]:
    Ap = aplus(sv*1j, 3000); D = Ap - aminus(sv*1j, 3000)
    if abs(D) > 1e-25: horn.append([float(Ap.imag), float(mp.log10(abs(D)))])
slope = (horn[-1][1] - horn[0][1])/(horn[-1][0] - horn[0][0])
check("|α₊ − α₋| falls like e^(−2π·Im α): slope of log₁₀|D| against Im α is −2π/ln 10 = −2.729", abs(slope + 2*math.pi/math.log(10)) < 0.01)
check("c₁ ≠ 0, so E has no analytic flow through 0 (Écalle–Voronin), matching Baker 1958's divergence", abs(cs[0]) > 1)

print(f"\n{ok}/{tot} checks pass")
json.dump({'horn': horn, 'c1': [float(cs[0].real), float(cs[0].imag)], 'alpha_coef': [float(c) for c in ac], 'v': [float(sp.expand(vser).coeff(x, i)) for i in range(2, N + 1)]}, open('series.json', 'w'))
