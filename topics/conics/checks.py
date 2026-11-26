"""Conics in Hypercube Composition Theory: checks. Run: python3 checks.py > checks_output.txt"""
import numpy as np, sympy as sp
rng = np.random.default_rng(7)
R = []
def rec(tag, ok, note):
    R.append((tag, bool(ok), note)); print(f"{tag}: {'PASS' if ok else 'FAIL'}  {note}")

s, t, k, d, a, b = sp.symbols('s t k d a b')
F = (a + b + k*a*b)/(1 - d*a*b)
q = lambda u: 1 + k*u + d*u**2
P = lambda u: (u/q(u), u**2/q(u))

# ---- A. the chord conic x^2 + kxy + dy^2 = y
x, y = P(s)
rec("A1", sp.simplify(x**2 + k*x*y + d*y**2 - y) == 0,
    "P(s) = (s/q, s^2/q), q = 1+ks+ds^2, lies on x^2+kxy+dy^2 = y (the line y = s x meets it again at P(s))")
xP, yP = P(s); xQ, yQ = P(t)
slope = sp.simplify((yP - yQ)/(xP - xQ))
rec("A2", sp.simplify(slope - F.subs({a: s, b: t})) == 0,
    "chord rule: slope of chord PQ = (s+t+kst)/(1-dst), so the line through O parallel to PQ hits P(s ∘ t). Exact (sympy)")
M = sp.Matrix([[1, k/2, 0], [k/2, d, -sp.Rational(1, 2)], [0, -sp.Rational(1, 2), 0]])
rec("A3", sp.simplify(M.det()) == -sp.Rational(1, 4),
    "3x3 matrix of the chord conic has det -1/4 for every (k,d): it is never degenerate")
ok = True
for _ in range(400):
    kk, dd = rng.normal(0, 2, 2)
    ev = np.linalg.eigvalsh([[1, kk/2], [kk/2, dd]])
    typ = 'ellipse' if ev.min() > 0 else 'hyperbola'
    ok &= (typ == 'ellipse') == (kk**2 - 4*dd < 0)
kk = 1.3; dd = kk**2/4
ev = np.linalg.eigvalsh([[1, kk/2], [kk/2, dd]])
ok &= abs(ev.min()) < 1e-12
rec("A4", ok, "affine type of the chord conic = sign of the plane dial's discriminant k^2-4d (400 random dials): "
    "inside the parabola d = k^2/4 ellipse, on it parabola (quadratic part rank 1, conic still smooth), outside hyperbola")
ok = True
for kk, dd in [(0, -1), (1.5, 0.2), (-1, 0), (2, 1)]:
    disc = kk**2 - 4*dd
    if disc <= 0: continue
    w = np.roots([dd, kk, 1]) if dd != 0 else np.array([-1/kk])
    for wi in w:
        u = wi + 1e-7
        pt = np.array([u/(1 + kk*u + dd*u*u), u*u/(1 + kk*u + dd*u*u)])
        ok &= np.linalg.norm(pt) > 1e5
rec("A5", ok, "the walls (roots of 1+ka+da^2) are exactly the directions where the chord conic runs off to infinity "
    "(its points at infinity); d=0 also has the wall at a=∞")
rec("A6", sp.simplify(F.subs({k: 0, d: 0})) == a + b and slope.subs({k: 0, d: 0}).simplify() == s + t,
    "dial origin (0,0): chord conic is the parabola y = x^2 and the chord rule is plain addition")

# ---- B. the unit conic of the plane number eps^2 = k eps - d
x1, y1, x2, y2 = sp.symbols('x1 y1 x2 y2')
N = lambda X, Y: X**2 + k*X*Y + d*Y**2
X3 = x1*x2 - d*y1*y2; Y3 = x1*y2 + x2*y1 + k*y1*y2
rec("B1", sp.expand(N(X3, Y3) - N(x1, y1)*N(x2, y2)) == 0,
    "norm N = x^2+kxy+dy^2 of x+y·eps is multiplicative (eps^2 = k eps - d)")
rec("B2", sp.simplify(Y3/X3 - F.subs({a: y1/x1, b: y2/x2})) == 0,
    "multiplying plane numbers composes their directions a = y/x by the plane law: the unit conic N = 1 carries the law")
T1, T2 = sp.symbols('T1 T2'); kap = d - k**2/4
amap = lambda T: T/(1 - (k/2)*T)
lhs = amap((T1 + T2)/(1 - kap*T1*T2)); rhs = F.subs({a: amap(T1), b: amap(T2)})
rec("B3", sp.simplify(lhs - rhs) == 0 and
    sp.expand((x1 + k*y1/2)**2 + kap*y1**2 - N(x1, y1)) == 0,
    "completing the square (X = x + ky/2) is the Cayley map a = T/(1-(k/2)T): N = X^2 + κY^2, law -> (T1+T2)/(1-κT1T2), κ = d - k^2/4")
X_, Y_ = sp.symbols('X Y')
f_e = sp.factor_list(X_**2 + Y_**2); f_p = sp.factor_list(X_**2); f_h = sp.factor_list(X_**2 - Y_**2)
rec("B4", len(f_e[1]) == 1 and f_p[1][0][1] == 2 and len(f_h[1]) == 2,
    "radius-0 level N = 0 over the reals: κ>0 only the point 0 (X^2+Y^2 irreducible); κ=0 a double line X^2; κ<0 two crossing lines (X-Y)(X+Y)")
ok = True
for kk, dd in [(0, -1), (-1, 0), (1, -2), (3, 1)]:
    w = np.roots([dd, kk, 1]) if dd != 0 else np.array([-1/kk])
    for wi in w: ok &= abs(1 + kk*wi + dd*wi*wi) < 1e-9   # N(1,w) = 0: ray is a null line = asymptote of N=1
rec("B5", ok, "the walls are the asymptote directions of the unit hyperbola N = 1 (N(1,w) = 0), i.e. its light cone")
rec("B6", sp.expand(N(1 + 0*y1, y1).subs({k: 0, d: 0})) == 1,
    "flat dial: N(1,a) = 1 for every a at (0,0): the unit 'circle' x^2 = 1 is two parallel lines x = ±1 and multiplication slides along x = 1")

# ---- C. the ruler is Cayley's log cross-ratio with the walls as the absolute
def ruler(aa, kk, dd, n=20001):
    u = np.linspace(0, aa, n); fu = 1/(1 + kk*u + dd*u*u)
    return np.trapezoid(fu, u)
ok = True
for kk, dd, aa in [(0, -1, 0.6), (1.5, 0.2, 0.3), (2, -1, 0.4), (-1, 0, 0.7)]:
    L = ruler(aa, kk, dd)
    if dd == 0:
        cay = np.log(1 + kk*aa)/kk
    else:
        w1, w2 = np.roots([dd, kk, 1])
        cay = np.log(((aa - w1)*(0 - w2))/((aa - w2)*(0 - w1)))/(dd*(w1 - w2))
    ok &= abs(L - np.real(cay)) < 1e-6
for kk, dd, aa in [(0, 1, 0.8), (1, 1, 0.5)]:
    kap = dd - kk*kk/4
    T = aa/(1 + kk*aa/2); cay = np.arctan(np.sqrt(kap)*T)/np.sqrt(kap)
    ok &= abs(ruler(aa, kk, dd) - cay) < 1e-6
kk = 2.0; dd = 1.0; aa = 0.7
ok &= abs(ruler(aa, kk, dd) - aa/(1 + kk*aa/2)) < 1e-6
rec("C1", ok, "the law's ruler ∫da/(1+ka+da^2) = Cayley's distance: (1/√disc)·log cross-ratio with the two walls (hyperbolic; "
    "d=0 gives log(1+ka)/k), arctan with no real walls (elliptic), a/(1+ka/2) with a double wall (parabolic)")

def hilbert(p, qq, inside, tmax=50):
    p = np.asarray(p, float); qq = np.asarray(qq, float); v = qq - p
    def edge(sign):
        lo, hi = 0.0, 1.0
        while inside(p + sign*hi*v) and hi < 1e6: hi *= 2
        for _ in range(200):
            m = (lo + hi)/2
            if inside(p + sign*m*v): lo = m
            else: hi = m
        return sign*lo
    t1 = edge(1.0); t0 = edge(-1.0)   # boundary at p + t v
    # points at params 0 (p), 1 (q), t0 (<0), t1 (>1)
    return 0.5*np.log(((t1 - 0)*(1 - t0))/((t1 - 1)*(0 - t0)))
disc_in = lambda z: z[0]**2 + z[1]**2 < 1
par_in = lambda z: z[1] > z[0]**2
to_disc = lambda z: np.array([2*z[0]/(z[1] + 1), (z[1] - 1)/(z[1] + 1)])
ok = True; ok2 = True
for _ in range(30):
    xs = rng.uniform(-1, 1, 2); p = [xs[0], xs[0]**2 + rng.uniform(0.05, 2)]
    xs = rng.uniform(-1, 1, 2); qq = [xs[1], xs[1]**2 + rng.uniform(0.05, 2)]
    dpar = hilbert(p, qq, par_in); P1, Q1 = to_disc(p), to_disc(qq)
    ddisc = hilbert(P1, Q1, disc_in)
    hyp = np.arccosh((1 - P1 @ Q1)/np.sqrt((1 - P1 @ P1)*(1 - Q1 @ Q1)))
    ok &= abs(dpar - ddisc) < 1e-6; ok2 &= abs(ddisc - hyp) < 1e-6
rec("C2", ok2, "Hilbert's distance inside a circle (ellipse) = the hyperbolic distance of the Klein model (30 random pairs)")
rec("C3", ok, "Hilbert's distance inside the parabola y > x^2 equals that inside the disc after the projective map "
    "(x,y) -> (2x/(y+1), (y-1)/(y+1)): the parabola gives the same hyperbolic plane")
half_in = lambda z: z[0]**2 + z[1]**2 < 1 and z[1] > 0
p, qq = np.array([0.0, 0.1]), np.array([0.0, 0.6])
dh, dd_ = hilbert(p, qq, half_in), hilbert(p, qq, disc_in)
br_to_disc = lambda z: np.array([z[0]/z[1], 1/z[1]])
inbranch = br_to_disc(np.array([0.7, 1.5])); ok = inbranch[1] > 0 and inbranch @ inbranch < 1
rec("C4", ok and abs(dh - dd_) > 0.05,
    f"inside one branch of y^2-x^2 = 1 maps to a half-disc (the line at infinity becomes a diameter); its Hilbert distance "
    f"differs from the disc's ({dh:.4f} vs {dd_:.4f}), so one branch alone is not the hyperbolic plane; both branches together are")

# ---- D. curvature: the dial is the conic's affine curvature, not its eccentricity
u = sp.symbols('u', real=True)
def affine_curv(kk, dd, uu):
    qq_ = 1 + kk*u + dd*u**2; r = sp.Matrix([u/qq_, u**2/qq_])
    r1 = r.diff(u); r2 = r1.diff(u)
    w = sp.Matrix.hstack(r1, r2).det()
    D = lambda v: v.diff(u)/sp.cbrt(w) if True else None
    rs = r1/sp.cbrt(w); rss = rs.diff(u)/sp.cbrt(w); rsss = rss.diff(u)/sp.cbrt(w)
    return complex(sp.Matrix.hstack(rss, rsss).det().subs(u, uu).evalf()).real
ok = True; rows = []
for kk, dd in [(0, 1), (1, 1), (0, 0.25), (2, 1), (0, -1), (1, -0.5), (-1, 0), (3, 1)]:
    kap = dd - kk**2/4
    vals = [affine_curv(sp.Rational(kk).limit_denominator(1000), sp.Rational(dd).limit_denominator(1000), uu)
            for uu in (0.13, 0.41)]
    # w may be negative for some param ranges: cube root handled via real cube root through orientation
    pred = 2**(4/3)*kap
    rows.append((kk, dd, kap, vals[0], vals[1], pred))
    ok &= all(abs(v - pred) < 1e-8 for v in vals)
rec("D1", ok, "equi-affine curvature of the chord conic is constant and equals 2^(4/3)·κ, κ = d - k^2/4 "
    + "; ".join(f"(k,d)=({r[0]},{r[1]}): {r[3]:.6f} vs {r[5]:.6f}" for r in rows))
def ecc(kk, dd):
    ev = np.linalg.eigvalsh([[1, kk/2], [kk/2, dd]])
    lp, lm = ev.max(), ev.min()
    if lm > 0: return np.sqrt(1 - lm/lp)
    return np.sqrt(1 + abs(lm)/lp) if lp > 0 else np.nan
e1, e2 = ecc(0.5, 1.0), ecc(0.5*2, 1.0*4)
iso = sp.simplify(F.subs({k: 2*k, d: 4*d}) - F.subs({a: 2*a, b: 2*b})/2)
rec("D2", iso == 0 and abs(e1 - e2) > 0.05,
    f"(k,d) and (λk, λ^2 d) are the same law (a -> λa), but the chord conics' eccentricities differ (λ=2: {e1:.4f} vs {e2:.4f}): eccentricity is a ruler choice")
ok = True
for _ in range(300):
    kk, dd = rng.normal(0, 2, 2)
    e = ecc(kk, dd)
    if np.isnan(e): continue
    ok &= np.sign(e*e - 1) == np.sign(kk*kk - 4*dd)
rec("D3", ok, "sign(e^2 - 1) = sign(k^2 - 4d) on 300 random dials: the dial fixes ellipse/parabola/hyperbola, not e itself")
ok = True
for _ in range(50):
    r = rng.normal(0, 1, 2); v = rng.normal(0, 1, 2); mu = 1.0
    E = v @ v/2 - mu/np.linalg.norm(r); L = r[0]*v[1] - r[1]*v[0]
    evec = (np.array([v[1]*L, -v[0]*L]) / mu) - r/np.linalg.norm(r)
    ok &= abs(evec @ evec - 1 - 2*E*L*L/mu**2) < 1e-9
rec("D4", ok, "Kepler: e^2 - 1 = 2EL^2/mu^2 (50 random states): energy (the plane dial, per gravity-from-theory) sets the type; angular momentum sets e")

# ---- E. the other dials touch the conics once
def lse_level(T, c=1.0, n=720):
    th = np.linspace(0, 2*np.pi, n, endpoint=False); pts = []
    e1 = np.array([1, -1, 0])/np.sqrt(2); e2 = np.array([1, 1, -2])/np.sqrt(6)
    for a_ in th:
        dvec = np.cos(a_)*e1 + np.sin(a_)*e2
        f = lambda rr: T*np.log(np.sum(np.exp(rr*dvec/T))) - T*np.log(3) - c
        lo, hi = 0.0, 50.0
        for _ in range(100):
            m = (lo + hi)/2
            lo, hi = (m, hi) if f(m) < 0 else (lo, m)
        pts.append(lo)
    return np.array(pts)
spreads = {T: (lambda r: (r.max() - r.min())/r.mean())(lse_level(T)) for T in (0.1, 1, 10, 100)}
rec("E1", spreads[100] < spreads[10] < spreads[1] < spreads[0.1] and abs(spreads[100]/spreads[10] - 10**-0.5) < 0.05,
    "roundness (Maslov temperature) dial: free-energy level set in the plane becomes a circle only in the hot limit (out-of-round spread falls like T^-1/2) "
    + ", ".join(f"T={T}: spread {v:.4f}" for T, v in spreads.items()) + "; cold end is the triangle (three lines)")
def conic_resid(p):
    th = np.linspace(0, np.pi/2, 200); xs = np.cos(th)**(2/p); ys = np.sin(th)**(2/p)
    A = np.column_stack([xs**2, xs*ys, ys**2, xs, ys, np.ones_like(xs)])
    return np.linalg.svd(A, compute_uv=False)[-1]
res = {p: conic_resid(p) for p in (1.5, 2, 3, 6)}
rec("E2", res[2] < 1e-12 and min(res[1.5], res[3], res[6]) > 1e-6,
    "ℓ^p roundness dial |x|^p+|y|^p = 1 is a conic only at p = 2 (smallest conic-fit residual: "
    + ", ".join(f"p={p}: {v:.1e}" for p, v in res.items()) + "); its ends, diamond and square, are pairs of parallel-line pairs")

print(f"\n{sum(r[1] for r in R)} of {len(R)} pass")
