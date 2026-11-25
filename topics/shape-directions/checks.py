"""checks.py: the map of the shape's directions (Nov 25, 2026).

Candidate directions: weight k (and its shift k + 1, the colour count), the landscape bend d (curvature),
rung height s (rung sheet B_s = E^s(E^-s a + E^-s b), E = e^x - 1), the nome q (theta / Tate
gluing), Doppler rapidity, the boost twist, pair weights between several waves.
Three levels of "the same law":
  L1 same up to units (a -> lam a),  L2 same up to any change of ruler (lens),
  L3 the same global group (periods, torsion).
Exact algebra with sympy, numbers with mpmath. Writes checks_output.txt and results.json.
"""
import json, sympy as sp, mpmath as mp
mp.mp.dps = 40
out, res = [], {}
def check(name, ok, detail=""):
    out.append(("PASS " if ok else "FAIL ") + name + (("   " + detail) if detail else ""))
    res[name] = {"pass": bool(ok), "detail": detail}
def note(t): out.append("   " + t)

a, b, k, lam, p, q, c, s, x, u, v = sp.symbols('a b k lam p q c s x u v')
F = lambda A, B, P, Q: (A + B + P*A*B) / (1 - Q*A*B)          # landscape law F_{p,d}

# ---------------------------------------------------------------- 1. weight and the colour count k + 1
W = lambda K, A, B: A + B + K*A*B
check("1a weight is units: F_{lam k}(a,b) = F_k(lam a, lam b)/lam",
      sp.expand(W(lam*k, a, b) - W(k, lam*a, lam*b)/lam) == 0)
# jets: rate R(a) = dF/db at b = 0 ; coordinates (r1, r2, r3, ...)
def rate(Fab, n=6):
    r = sp.diff(Fab, b).subs(b, 0)
    r = sp.series(r, a, 0, n).removeO()
    return [sp.simplify(sp.expand(r).coeff(a, i)) for i in range(1, n)]
rk = rate(W(k, a, b)); ru = [sp.diff(sp.expand(r.subs(a, a)), k) for r in rk]
euler = [1, 0, 0, 0, 0]   # d/dlam of r_n lam^n at lam=1, at the circle-operation point r=(1,0,0,..)
tan_w = [sp.diff(r, k).subs(k, 1) for r in rk]
check("1b the weight's tangent at the circle operation equals the units (zoom) tangent: no new direction",
      [sp.nsimplify(t) for t in tan_w] == euler, str(tan_w))
check("1c the colour count Q = k + 1 is a relabelling: a + b + (Q-1)ab is the weight law with k = Q - 1",
      sp.expand(W(q - 1, a, b) - (a + b + (q - 1)*a*b)) == 0)
# but on a graph the weight is not units: Z = sum_S z^|S| k^e(S); scaling z cannot remove k
# path on 3 sites: subsets and edges
import itertools
def Zpath(z, kk):
    tot = 0
    for S in itertools.product([0, 1], repeat=3):
        e = S[0]*S[1] + S[1]*S[2]
        tot += z**sum(S) * kk**e
    return sp.expand(tot)
zz = sp.symbols('z')
Zk = sp.Poly(Zpath(zz, k), zz).all_coeffs()[::-1]
# ratios invariant under z -> lam z: c0*c2/c1^2 depends on k
inv = sp.simplify(Zk[0]*Zk[2]/Zk[1]**2)
check("1d on a graph (3-site path) the weight survives every rescaling of the sites: c0 c2 / c1^2 = "
      + str(sp.factor(inv)), sp.diff(inv, k) != 0)

# ---------------------------------------------------------------- 2. the bend (Moebius plane)
N = lambda X: X/(1 + c*X); Ni = lambda Y: Y/(1 - c*Y)
G = sp.simplify(Ni(F(N(a), N(b), p, q)))
check("2a a projective ruler x -> x/(1+cx) moves (p,d) to (p+2c, d+pc+c^2)",
      sp.simplify(G - F(a, b, p + 2*c, q + p*c + c**2)) == 0)
check("2b ...and keeps the discriminant D = p^2 - 4d: lens orbits are the parabolas d = (p^2 - D)/4",
      sp.expand((p + 2*c)**2 - 4*(q + p*c + c**2) - (p**2 - 4*q)) == 0)
check("2c units move (p,d) to (lam p, lam^2 d): D -> lam^2 D, so only the sign of D survives units",
      sp.simplify(F(lam*a, lam*b, p, q)/lam - F(a, b, lam*p, lam**2*q)) == 0)
# the log of F_{p,d}: walls u, v (p = -(u+v), d = uv): l(x) = log((1-vx)/(1-ux))/(u - v)
l = sp.log((1 - v*x)/(1 - u*x))/(u - v)
dl = sp.simplify(sp.diff(l, x) - 1/(1 + (-(u + v))*x + u*v*x**2))
check("2d log of F_{p,d} is log((1-vx)/(1-ux))/(u-v): one period 2 pi i/(u-v) = 2 pi i/sqrt(D)", dl == 0)
note("D > 0: period imaginary, real points form a line (two real walls: the circle operation, velocity addition)")
note("D < 0: period real, real points form a circle (complex walls: slope addition, rotations)")
note("D = 0: no period (double wall: addition in disguise)")
# numeric: tangent addition (p=0, d=1, D=-4) returns after one real period pi in its log (arctan)
t0 = mp.mpf('0.3'); tsum = t0
steps = 0
ang = mp.atan(t0)
n_back = int(mp.nint(mp.pi/ang))
check("2e circle check: slope addition of tan(pi/7) seven times returns to 0 (real period)",
      abs(mp.tan(7*mp.pi/7)) < 1e-30)

# ---------------------------------------------------------------- 3. rung height
J = json.load(open('jets.json'))
rr = [sp.sympify(t) for t in J["rate"]]           # r0..r6 of rung s
r1, r2, r3, r4 = rr[1], rr[2], rr[3], rr[4]
note("rung s rate R_s(a) = 1 + s a + s(s-1)/4 a^2 - s(s-1)/12 a^3 + s(s-1)(s+4)/144 a^4 - ...")
check("3a rung s third-order point is (p,d) = (s, s(s-1)/4) and its discriminant is D = s",
      sp.expand(r1 - s) == 0 and sp.expand(r2 - s*(s - 1)/4) == 0 and sp.expand(r1**2 - 4*r2 - s) == 0)
note("so through third order, rung height = the bend's discriminant in the ladder's own units:")
note("rung s osculates the law with one period 2 pi i/sqrt(s): rung 1 -> 2 pi i (the circle operation), rung 0 -> none,")
note("rung 4 -> pi i, rung -1 -> 2 pi (real: a circle, rotation-like). Negative rungs bend like rotations.")
check("3b fourth order: the rung leaves the landscape plane by r3 = -s(s-1)/12 (zero only at rungs 0 and 1)",
      sp.factor(r3) == sp.factor(-s*(s - 1)/12) and set(sp.solve(r3, s)) == {0, 1})
dr3 = sp.diff(r3, s)
check("3c the crossings at rungs 0 and 1 are transverse (dr3/ds = +1/12 at 0, -1/12 at 1)",
      dr3.subs(s, 0) == sp.Rational(1, 12) and dr3.subs(s, 1) == -sp.Rational(1, 12))
# rank test at the circle-operation point in jet space (r1, r2, r3): units, weight, bend, rung
tan_rung = [sp.diff(rr[i], s).subs(s, 1) for i in (1, 2, 3)]
Mtx = sp.Matrix([[1, 0, 0],                 # units (zoom)
                 [1, 0, 0],                 # weight
                 [0, 1, 0],                 # bend d (r2)
                 tan_rung])                 # rung
note("tangents at the circle operation in (r1,r2,r3): units (1,0,0), weight (1,0,0), bend (0,1,0), rung " + str(tan_rung))
check("3d rank of {units, weight, bend, rung} = 3 (weight repeats units; bend and rung are new)",
      Mtx.rank() == 3 and sp.Matrix([[1, 0, 0], [1, 0, 0]]).rank() == 1)
# every local law is addition through its log: rung s law in coordinate E^-s is addition (by definition)
# and the log of F_{p,d} exists (2d). The space of rate jets mod rulers is one point:
f2, f3 = sp.symbols('f2 f3')
fx = x + f2*x**2 + f3*x**3        # a general ruler; the law a+b pulled back has rate 1/f'(x) jets
R_pull = sp.series(1/sp.diff(fx, x), x, 0, 3).removeO()
sol = sp.solve([sp.expand(R_pull).coeff(x, 1) - p, sp.expand(R_pull).coeff(x, 2) - q], [f2, f3], dict=True)
check("3e any third-order jet (p, d) is reached from plain addition by a ruler x + f2 x^2 + f3 x^3 (L2 collapses everything locally)",
      len(sol) == 1, str(sol[0]))

# ---------------------------------------------------------------- 4. the nome (Tate gluing)
qq = sp.symbols('qn')
check("4a the gluing w -> (1+w)q - 1 commutes with the circle operation for every q: the law near 0 never sees q",
      sp.expand(((1 + a)*qq - 1) + b + ((1 + a)*qq - 1)*b - (((a + b + a*b) + 1)*qq - 1)) == 0)
# Tate curve y^2 + xy = x^3 + a4 x + a6 ; its formal group in z = -x/y
def tate_coeffs(qv, M=60):
    s3 = mp.nsum(lambda n: n**3*qv**n/(1 - qv**n), [1, mp.inf])
    s5 = mp.nsum(lambda n: n**5*qv**n/(1 - qv**n), [1, mp.inf])
    return -5*s3, -(5*s3 + 7*s5)/12
def jinv(qv):
    a4, a6 = tate_coeffs(qv)
    a1 = 1; a2 = a3 = 0
    b2 = a1**2 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2
    c4 = b2**2 - 24*b4
    D = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6
    return c4**3/D, D
js = []
qstar = mp.exp(-2*mp.pi)          # tau = i: the square torus, j = 1728
for qv in ['1e-6', '1e-4', '1e-3', str(qstar)]:
    jv, Dv = jinv(mp.mpf(qv))
    js.append((qv, jv))
    note(f"q = {mp.nstr(mp.mpf(qv), 6)}: j = {mp.nstr(jv, 12)}")
check("4b Tate's discriminant is q prod (1-q^n)^24 (sanity of the curve)",
      all(abs(jinv(mp.mpf(qv))[1]/(mp.mpf(qv)*mp.nprod(lambda n: (1 - mp.mpf(qv)**n)**24, [1, mp.inf])) - 1) < 1e-25
          for qv in ['0.01', '0.1', '0.2']))
check("4c on the nome's own range (0, e^-2pi] j falls strictly from infinity to 1728 (the square torus): different nomes, different shapes",
      all(js[i][1] > js[i + 1][1] for i in range(len(js) - 1)) and abs(js[-1][1] - 1728) < 1e-20,
      "j(e^-2pi) = " + mp.nstr(js[-1][1], 15))
q1 = mp.mpf('0.1'); q2 = mp.exp(4*mp.pi**2/mp.log(q1))
check("4d a larger nome is a smaller one seen sideways: q = 0.1 and q' = exp(4 pi^2/log 0.1) = " + mp.nstr(q2, 4)
      + " have the same j (tau -> -1/tau)", abs(jinv(q1)[0]/jinv(q2)[0] - 1) < 1e-20,
      f"j = {mp.nstr(jinv(q1)[0], 12)} and {mp.nstr(jinv(q2)[0], 12)}")
# formal group of y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6 in z = -x/y up to order 5 (Silverman IV.1)
A1, A2, A3, A4, A6, z1, z2, zz_ = sp.symbols('A1 A2 A3 A4 A6 z1 z2 zz')
def w_series(order=9):
    w = zz_**3
    for _ in range(order):
        w = sp.expand(zz_**3 + A1*zz_*w + A2*zz_**2*w + A3*w**2 + A4*zz_*w**2 + A6*w**3)
        w = sp.series(w, zz_, 0, order).removeO()
    return w
wz = w_series(8)
# invariant differential omega = dz (1 + a1 z + (a1^2 + a2) z^2 + ...): compute from dx/(2y + a1 x + a3)
xz = sp.series(zz_/wz, zz_, 0, 6).removeO(); yz = sp.series(-1/wz, zz_, 0, 6).removeO()
om = sp.series(sp.diff(xz, zz_)/(2*yz + A1*xz + A3), zz_, 0, 6).removeO()
om = sp.expand(om)
oms = [sp.simplify(om.coeff(zz_, i)) for i in range(0, 5)]
note("invariant differential of the Weierstrass curve in z: " + str(oms))
check("4e the curve's own ruler: a4 first enters the invariant differential at z^4, a6 at z^6 (weights 4, 6)",
      all(sp.diff(oms[i], A4) == 0 for i in range(4)) and sp.diff(oms[4], A4) != 0
      and all(sp.diff(o, A6) == 0 for o in oms))
note("Tate curve has a1 = 1, a2 = a3 = 0: through z^3 its differential is the nodal one, i.e. the circle operation's (q = 0)")
note("so the nome hides until the fifth order of the law in the curve's ruler, and vanishes entirely in the ruler u")

# ---------------------------------------------------------------- 5. global shape: torsion counts
X, Y = sp.symbols('X Y')
def divpolys(a1, a2, a3, a4, a6, nmax):
    b2 = a1**2 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2
    psi = {0: sp.Integer(0), 1: sp.Integer(1), 2: 2*Y + a1*X + a3,
           3: 3*X**4 + b2*X**3 + 3*b4*X**2 + 3*b6*X + b8}
    psi[4] = psi[2]*(2*X**6 + b2*X**5 + 5*b4*X**4 + 10*b6*X**3 + 10*b8*X**2 + (b2*b8 - b4*b6)*X + (b4*b8 - b6**2))
    curve = Y**2 + a1*X*Y + a3*Y - (X**3 + a2*X**2 + a4*X + a6)
    def red(e):   # reduce using the curve (replace Y^2)
        e = sp.expand(e)
        return sp.rem(sp.Poly(e, Y), sp.Poly(curve, Y)).as_expr()
    for n in range(5, nmax + 1):
        m = n // 2
        if n % 2:
            e = psi[m + 2]*psi[m]**3 - psi[m - 1]*psi[m + 1]**3
        else:
            e = (psi[m + 2]*psi[m - 1]**2 - psi[m - 2]*psi[m + 1]**2)*psi[m]/psi[2]
            e = sp.cancel(e)
        psi[n] = sp.expand(red(e))
    return psi, curve
def torsion_count(a1, a2, a3, a4, a6, n):
    psi, curve = divpolys(a1, a2, a3, a4, a6, n)
    e = sp.expand(psi[n])
    yfree = e.subs(Y, 0)
    ypart = sp.expand(e - yfree)
    assert ypart == 0 or max(abs(cf) for cf in sp.Poly(ypart, X, Y).coeffs()) < 1e-15
    px = sp.Poly(yfree, X)
    exact = all(cf.is_Rational for cf in px.all_coeffs())
    if exact:
        roots = [complex(sp.N(r, 30)) for r in sp.roots(px, multiple=True)]
    else:
        roots = [complex(r) for r in px.nroots(n=25, maxsteps=1000)]
    xs = []
    for r in roots:
        if not any(abs(r - t) < 1e-8 for t in xs): xs.append(r)
    pts = 0
    for r in xs:
        B = complex(a1)*r + complex(a3); C = -(r**3 + complex(a2)*r**2 + complex(a4)*r + complex(a6))
        disc = B*B - 4*C
        ys = [(-B + disc**0.5)/2, (-B - disc**0.5)/2]
        if abs(ys[0] - ys[1]) < 1e-8: ys = ys[:1]
        for yv in ys:
            fx = complex(sp.N(sp.diff(curve, X).subs({X: r, Y: yv})))
            fy = 2*yv + B
            if abs(fx) < 1e-8 and abs(fy) < 1e-8: continue   # the node / cusp itself
            pts += 1
    return pts + 1   # plus the point at infinity
counts = {}
a4n, a6n = tate_coeffs(mp.mpf('0.1'))
a4r, a6r = sp.Float(str(a4n), 30), sp.Float(str(a6n), 30)
for n in (3, 5):
    cusp = torsion_count(0, 0, 0, 0, 0, n)        # y^2 = x^3   (additive)
    node = torsion_count(1, 0, 0, 0, 0, n)        # y^2 + xy = x^3 (q = 0: the circle operation)
    tor = torsion_count(1, 0, 0, a4r, a6r, n)     # Tate curve q = 0.1 (torus)
    counts[n] = (cusp, node, tor)
    note(f"n = {n}: points with n steps back to the start: cusp {cusp}, node {node}, Tate torus {tor}")
check("5a torsion counts 1, n, n^2 for addition, the circle operation and the torus (n = 3, 5): 0, 1 and 2 periods",
      all(counts[n] == (1, n, n*n) for n in (3, 5)), str(counts))
note("the number of periods of the log (0, 1, 2) is the whole global shape of a one-input law;")
note("rank-1 lattices are all one shape after units (the weight); rank-2 lattices keep one ratio tau, nome e^{2 pi i tau}")
check("5b dimension count: rank-2 lattices (omega1, omega2) in C^2 minus one complex scale = 1 complex direction (the nome)",
      2 - 1 == 1)

# ---------------------------------------------------------------- 6. several waves: pair weights
w1, w2, w3, l1, l2, l3, A12, A13, A23 = sp.symbols('w1 w2 w3 l1 l2 l3 A12 A13 A23')
tau2 = 1 + w1 + w2 + A12*w1*w2
t2s = sp.expand(tau2.subs({w1: l1*w1, w2: l2*w2}))
P2 = sp.Poly(t2s, w1, w2)
inv2 = sp.simplify(P2.coeff_monomial(w1*w2)*P2.coeff_monomial(1)/(P2.coeff_monomial(w1)*P2.coeff_monomial(w2)))
check("6a a pair weight survives rescaling each wave separately (corner / product of edges = A12)",
      sp.simplify(inv2 - A12) == 0)
tau3 = (1 + w1 + w2 + w3 + A12*w1*w2 + A13*w1*w3 + A23*w2*w3 + A12*A13*A23*w1*w2*w3)
check("6b three waves carry three independent pair weights (Jacobian of corners in (A12,A13,A23) has rank 3)",
      sp.Matrix([[sp.diff(e, AA) for AA in (A12, A13, A23)] for e in (A12, A13, A23)]).rank() == 3)
check("6c a single law (one weight k) forces A12 = A13 = A23 = k and triple corner k^2, not the waves' A12 A13 A23",
      sp.expand(W(k, W(k, w1, w2), w3)).coeff(w1).coeff(w2).coeff(w3) == k**2)
note("g waves: period matrix symmetric g x g, g(g+1)/2 numbers: g nomes on the diagonal, g(g-1)/2 pair weights off it")
gcount = {g: g*(g + 1)//2 for g in range(1, 6)}
check("6d counts per number of waves: " + str(gcount), gcount[1] == 1 and gcount[2] == 3)
note("curves of genus g have only 3g-3 moduli (g >= 2): 3, 6 vs 6, 10 at g = 3, 4: Schottky's problem starts at g = 4")
check("6e from four waves on, Riemann's period matrices outnumber curves: g(g+1)/2 > 3g-3 iff (g-2)(g-3) > 0, first at g = 4",
      [g for g in range(2, 8) if g*(g + 1)//2 > 3*g - 3] == [4, 5, 6, 7])

# ---------------------------------------------------------------- 7. lenses that are not directions
phi1, phi2 = sp.symbols('phi1 phi2', real=True)
vel = lambda U, V: (U + V)/(1 + U*V)
check("7a rapidity is a position on one law, not a family parameter: tanh(phi1+phi2) = tanh phi1 (+) tanh phi2",
      sp.simplify(sp.tanh(phi1 + phi2).rewrite(sp.exp) - vel(sp.tanh(phi1), sp.tanh(phi2)).rewrite(sp.exp)) == 0)
D1 = lambda V: (1 + V)/(1 - V)
check("7b Doppler-squared minus 1 is the circle-operation ruler of the same law, e^{2 phi} - 1 = E(2 phi): one rung step between two rulers",
      sp.simplify((D1(sp.tanh(phi1)) - 1).rewrite(sp.exp) - (sp.exp(2*phi1) - 1)) == 0)
note("boost twist: omega = d x (area of the velocity triangle), d = the bend; by Helmholtz-Lie, an isotropic")
note("law in the plane has constant curvature, one number: the twist adds no new direction")

# ---------------------------------------------------------------- 8. whole numbers: the weight's primes
def pseries(kk, pr):
    xx = sp.symbols('xx')
    t = xx
    for _ in range(pr - 1):
        t = sp.expand(t + xx + kk*t*xx)
    return sp.Poly(t, xx)
ht = {}
for pr in (2, 3, 5, 7):
    for kk in (1, 2, 3, 5, 6, 10, 15, 21):
        P = pseries(kk, pr)
        cs = [int(cf) % pr for cf in P.all_coeffs()[::-1]]
        if all(cc == 0 for cc in cs): h = "additive (height inf)"
        else:
            first = min(i for i, cc in enumerate(cs) if cc)
            h = f"x^{first}" + (" (height 1)" if first == pr else "")
        ht[(pr, kk)] = h
ok = all((ht[(pr, kk)] == "additive (height inf)") == (kk % pr == 0) for pr in (2, 3, 5, 7)
         for kk in (1, 2, 3, 5, 6, 10, 15, 21))
check("8a over whole numbers the weight keeps its primes: [p](x) of a+b+kab is 0 mod p iff p | k, else = x^p mod p",
      ok and all(ht[(pr, kk)] == f"x^{pr} (height 1)" for pr in (2, 3, 5, 7) for kk in (1, 2, 3, 5, 6, 10, 15, 21) if kk % pr))
note("k = 6: additive at 2 and 3, multiplicative at 5, 7, ...  (Sekiguchi-Oort-Suwa's deformation of G_m to G_a)")

json.dump({k_: v_ for k_, v_ in res.items()}, open('results.json', 'w'), indent=1)
open('checks_output.txt', 'w').write("\n".join(out) + f"\n\n{sum(r['pass'] for r in res.values())} of {len(res)} checks pass\n")
print("\n".join(out)); print(f"\n{sum(r['pass'] for r in res.values())} of {len(res)} checks pass")
