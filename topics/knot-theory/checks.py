# Knot theory through the circle operation: every claim on the page, checked (sympy, exact where possible).
# Uses tangles.py: tangle diagrams (state sum over the cube) and the (f, g) tangle algebra.
import sympy as sp, cmath, math, random, itertools
from tangles import *
random.seed(1)
ok = 0; tot = 0
def check(name, cond):
    global ok, tot
    tot += 1; ok += bool(cond); print(('PASS ' if cond else 'FAIL ') + name)
z = lambda e: sp.simplify(e) == 0
a, b, c, x, y, w, k, dd, f1, g1, f2, g2, f3, g3, qq, u = sp.symbols('a b c x y w k dd f1 g1 f2 g2 f3 g3 qq u')
wadd = lambda kk, p, r: p + r + kk*p*r
q = -A**-4                      # the tangle q (= -t, t the Jones variable)

print("--- 1. A tangle's bracket is a point (f : g) on a line; skein relations are linear ---")
check("crossing = A·[0] + A⁻¹·[∞]: diagram of N(X) and D(X) match the pair (A, A⁻¹)",
      z(bracket(closeN(cross(1))) - cN(XP)) and z(bracket(closeD(cross(1))) - cD(XP)))
check("mirror crossing = A⁻¹·[0] + A·[∞]",
      z(bracket(closeN(cross(-1))) - cN(XM)) and z(bracket(closeD(cross(-1))) - cD(XM)))
check("loop value d = −A² − A⁻²; disjoint loop: [0] stacked on [0] is d·[0]",
      vs(Z0, Z0) == (sp.expand(d), 0))
# random words: the algebra agrees with the full state sum over the cube, both closures
words = [[random.choice(['h+', 'h-', 'v+', 'v-']) for _ in range(random.randint(1, 7))] for _ in range(40)]
words = [['h+'] + wd for wd in words]
good = True
for wd in words:
    Dg, P = build(wd)
    good &= z(bracket(closeN(Dg)) - cN(P)) and z(bracket(closeD(Dg)) - cD(P))
check("tangle algebra = state sum over the cube, 40 random tangles up to 8 crossings, both closures", good)

check("the four corners of the sum: [0]+[0]=[0], [0]+[∞]=[∞]+[0]=[∞], [∞]+[∞] traps a loop = d·[∞] (diagrams)",
      z(bracket(closeN(hsum(zero(), zero()))) - cN(Z0)) and z(bracket(closeN(hsum(zero(), inf()))) - 1)
      and z(bracket(closeN(hsum(inf(), zero()))) - 1) and z(bracket(closeN(hsum(inf(), inf()))) - d)
      and z(bracket(closeD(hsum(inf(), inf()))) - d*d))
print("--- 2. Adding tangles is the dial with weight d (the loop value) ---")
S, T, U = (f1, g1), (f2, g2), (f3, g3)
tsg = lambda S, T: (S[0]*T[0], S[0]*T[1] + S[1]*T[0] + dd*S[1]*T[1])
check("sum formula: (f1,g1)+(f2,g2) = (f1f2, f1g2 + g1f2 + d·g1g2)",
      all(z(sp.expand(ts((f1, g1), (f2, g2))[i]) - [f1*f2, f1*g2 + g1*f2 + d*g1*g2][i]) for i in range(2)))
check("in x = g/f: x(S+T) = x(S) ∘_d x(T)  (the circle operation, weight d)",
      z(tsg(S, T)[1]/tsg(S, T)[0] - wadd(dd, g1/f1, g2/f2)))
check("sum is commutative and associative, [0] = (1, 0) is the identity",
      all(z(tsg(S, T)[i] - tsg(T, S)[i]) and z(tsg(tsg(S, T), U)[i] - tsg(S, tsg(T, U))[i]) for i in range(2))
      and tsg(S, (1, 0)) == (f1, g1))
cDg = lambda S: S[0] + dd*S[1]; cNg = lambda S: dd*S[0] + S[1]
check("the lift is closing: 1 + d·x = ⟨D(T)⟩/f  (D joins left ends and right ends)", z(1 + dd*g1/f1 - cDg(S)/f1))
check("lift law = connected sum: ⟨D(S+T)⟩ = ⟨D(S)⟩·⟨D(T)⟩", z(cDg(tsg(S, T)) - cDg(S)*cDg(T)))
check("wall: if ⟨D(S)⟩ = 0 then ⟨D(S+T)⟩ = 0 (x = −1/d absorbs)",
      z(cDg(tsg((f1, -f1/dd), T))) and z(wadd(dd, -1/dd, x) + 1/dd))
check("second wall: [∞] = (0, 1) absorbs: S + [∞] = (0, ⟨D(S)⟩), i.e. x = ∞",
      tsg(S, (0, 1)) == (0, f1 + dd*g1))
vsg = lambda S, T: (S[0]*T[1] + S[1]*T[0] + dd*S[0]*T[0], S[1]*T[1])
check("stacking is the same dial in z = f/g = 1/x, and ⟨N⟩ is its lift: ⟨N(S*T)⟩ = ⟨N(S)⟩⟨N(T)⟩",
      z(vsg(S, T)[0]/vsg(S, T)[1] - wadd(dd, f1/g1, f2/g2)) and z(cNg(vsg(S, T)) - cNg(S)*cNg(T)))
# diagrams really are connected sums: Jones multiplies (through the writhe-normalized diagrams)
s1, s2 = build(cf_word([2, 2]))[0], build(cf_word([1, 1, 1]))[0]
V1, V2, V12 = jones(closeD(s1)), jones(closeD(s2)), jones(closeD(hsum(s1, s2)))
check("diagram check: V(D(S+T)) = V(D(S))·V(D(T)) (figure-eight # trefoil)", z(V12 - V1*V2))
check("'circle-operation value' P − 1 turns # into ∘₁ exactly (a relabelling of P(K#L) = P(K)P(L))",
      z((V12 - 1) - wadd(1, V1 - 1, V2 - 1)))
t_ = sp.symbols('t_', positive=True)
check("disjoint union multiplies the Jones polynomial by d = −(t^½ + t^−½) (unlink of 2)",
      z(sp.expand(d.subs(A, t_**sp.Rational(-1, 4))) + t_**sp.Rational(1, 2) + t_**sp.Rational(-1, 2)))

print("--- 3. Twisting is the q-dial step y -> q·y + 1, with q = -A^-4 = -t ---")
check("y = A²·g/f: the sum becomes y ∘ₖ y' with k = d·A⁻² = q − 1",
      z(A**2*wadd(d, x, y) - wadd(d*A**-2, A**2*x, A**2*y)) and z(d*A**-2 - (q - 1)))
qn = lambda qv, n: sum(qv**i for i in range(n))
P = Z0; good = True; good2 = True
for n in range(1, 13):
    P = ts(P, XP)
    good &= z(P[0] - A**n) and z(P[1] - A**(n - 2)*qn(q, n))
    good2 &= z(cD(P) - (-A**-3)**n) and z(A**2*P[1]/P[0] - qn(q, n))
check("n twists [n] = A^n·[0] + A^(n−2)·[n]_q·[∞], n = 1..12", good)
check("y([n]) = [n]_q = n copies of 1 under ∘_(q−1); closing gives n kinks (−A⁻³)^n", good2)
yP = lambda P: sp.simplify(A**2*P[1]/P[0])
P0 = (f1, g1); y0 = A**2*g1/f1
check("adding one crossing: y ↦ q·y + 1 = y ∘_(q−1) 1", z(yP(ts(P0, XP)) - (q*y0 + 1)) and z(q*y0 + 1 - wadd(q - 1, y0, 1)))
check("stacking one crossing: y ↦ q·y/(q·y + 1), i.e. 1/y ↦ q⁻¹·(1/y) + 1 (the same step at 1/q, seen through 1/y)",
      z(yP(vs(P0, XP)) - q*y0/(q*y0 + 1)) and z(1/(q*y0/(q*y0 + 1)) - (1/(q*y0) + 1)))
check("mirror A ↦ A⁻¹ sends q to 1/q", z(q.subs(A, 1/A) - 1/q))
# twist matrix: (f,g) -> (A f, A^-1 f - A^-3 g)
M = sp.Matrix([[A, 0], [A**-1, -A**-3]])
check("twist acts on (f, g) by the matrix [[A,0],[A⁻¹,−A⁻³]]", z(sp.Matrix(ts(P0, XP)) - M*sp.Matrix(P0)) if False else
      all(z(ts(P0, XP)[i] - (M*sp.Matrix([f1, g1]))[i]) for i in range(2)))
check("its eigenvalues are A and −A⁻³; they coincide exactly when A⁴ = −1",
      set(M.eigenvals().keys()) == {A, -A**-3} and set(sp.solve(sp.Eq(A, -A**-3), A)) == set(sp.solve(A**4 + 1, A)))
check("at A⁴ = −1 the loop value is 0 and q = 1", all(abs(complex(d.subs(A, r))) < 1e-12 and abs(complex(q.subs(A, r)) - 1) < 1e-12
                                                     for r in sp.solve(A**4 + 1, A)))
check("walls of the twist: fixed points of y ↦ qy+1 are 1/(1−q) and ∞ (x-walls −1/d and ∞); they merge at q = 1",
      z(q*(1/(1 - q)) + 1 - 1/(1 - q)) and sp.limit(1/(1 - qq), qq, 1, '+') in (sp.oo, -sp.oo))
good = True
for m in range(3, 13):
    qv = cmath.exp(2j*math.pi/m); cen = 1/(1 - qv); R = abs(cen)
    pts = [sum(qv**i for i in range(n)) for n in range(m)]
    good &= all(abs(abs(p - cen) - R) < 1e-12 for p in pts) and abs(sum(qv**i for i in range(m))) < 1e-12
    good &= all(abs(abs(pts[(i+1) % m] - pts[i]) - 1) < 1e-12 for i in range(m))
check("|q| = 1: the twists [n]_q are the corners of a regular polygon of side 1 through 0, centre 1/(1−q); q = e^(2πi/m) closes after m twists", good)
check("as q → 1 the polygon opens into 0, 1, 2, 3, … (Conway's integer tangles)",
      all(abs(sum(cmath.exp(1e-7j)**i for i in range(n)) - n) < 1e-5 for n in range(8)))

print("--- 4. At q = 1 (A = e^(iπ/4), d = 0): Conway's fraction, continued fractions, determinants ---")
zeta = sp.exp(sp.I*sp.pi/4)
def Fz(P): return complex(sp.N((A**2*P[1]/P[0]).subs(A, zeta), 30))
good = True; ncount = 0
for wd in words:
    Dg, P = build(wd)
    Fr = sp.Integer(0); bad = False
    for s_ in wd:
        sg = 1 if s_[1] == '+' else -1
        if s_[0] == 'h': Fr = Fr + sg
        elif Fr == 0: pass                      # stacking onto [0] only adds a kink
        elif 1/Fr + sg == 0: bad = True; break  # passes through [inf]
        else: Fr = 1/(1/Fr + sg)
    if bad or Fr == 0: continue
    ncount += 1
    good &= abs(Fz(P) - complex(Fr)) < 1e-9
good &= ncount >= 15
check("F = A²g/f at A = e^(iπ/4) is real and h adds ±1 to F, v adds ±1 to 1/F (random tangles with finite F)", good)
check("F is additive under +, 1/F additive under stacking (general (f,g), d = 0)",
      z((A**2*tsg(S, T)[1]/tsg(S, T)[0]).subs(dd, 0) - A**2*g1/f1 - A**2*g2/f2) and
      z((vsg(S, T)[0]/(A**2*vsg(S, T)[1])).subs(dd, 0) - f1/(A**2*g1) - f2/(A**2*g2)))
cfs = {'3_1': [3], '4_1': [1, 1, 2], '5_1': [5], '5_2': [1, 1, 3], '6_1': [1, 1, 4], '6_2': [2, 1, 3], '6_3': [1, 1, 1, 1, 2]}
dets = {'3_1': 3, '4_1': 5, '5_1': 5, '5_2': 7, '6_1': 9, '6_2': 11, '6_3': 13}
fracs = {'3_1': sp.Integer(3), '4_1': sp.Rational(5, 2), '5_1': sp.Integer(5), '5_2': sp.Rational(7, 2),
         '6_1': sp.Rational(9, 2), '6_2': sp.Rational(11, 3), '6_3': sp.Rational(13, 5)}
check("continued fractions: h^a1 v^a2 h^a3… gives a_last + 1/(… + 1/a1): 3, 5/2, 5, 7/2, 9/2, 11/3, 13/5",
      all(cf_value(cfs[kk]) == fracs[kk] for kk in cfs))
VK = {}
good = True; goodg = True
for kk, cf in cfs.items():
    Dg, P = build(cf_word(cf)); K = closeN(Dg); VK[kk] = A_to_t(jones(K))
    good &= abs(VK[kk].subs(t, -1)) == dets[kk]
    gz = complex(sp.N(P[1].subs(A, zeta))); fz = complex(sp.N(P[0].subs(A, zeta)))
    goodg &= abs(abs(gz) - fracs[kk].p) < 1e-9 and abs(abs(fz) - fracs[kk].q) < 1e-9
check("determinant |V(−1)| of N(p/q) = p (numerator): 3, 5, 5, 7, 9, 11, 13", good)
check("at A = e^(iπ/4): |⟨N(T)⟩| = |g| = p and |⟨D(T)⟩| = |f| = q (Goldman–Kauffman)", goodg)
atlas = {'3_1': t + t**3 - t**4, '4_1': t**2 - t + 1 - 1/t + t**-2,
         '5_1': t**2 + t**4 - t**5 + t**6 - t**7,
         '5_2': -t**-6 + t**-5 - t**-4 + 2*t**-3 - t**-2 + t**-1,
         '6_1': t**2 - t + 2 - 2/t + t**-2 - t**-3 + t**-4,
         '6_2': t**-1 - 1 + 2*t - 2*t**2 + 2*t**3 - 2*t**4 + t**5,
         '6_3': -t**3 + 2*t**2 - 2*t + 3 - 2/t + 2*t**-2 - t**-3}
mir = lambda e: sp.expand(e.subs(t, 1/t))
for kk in atlas:
    check(f"Jones of the closed tangle for {kk} matches the Knot Atlas (up to mirror)",
          z(VK[kk] - atlas[kk]) or z(VK[kk] - mir(atlas[kk])))
D72, _ = build(cf_word([1, 1, 3])); D73, _ = build(cf_word([1, 2, 2]))
check("Schubert: N(7/3) is the mirror of N(7/2) (3 ≡ −2⁻¹ mod 7): Jones polynomials mirror",
      cf_value([1, 2, 2]) == sp.Rational(7, 3) and z(A_to_t(jones(closeN(D73))) - mir(A_to_t(jones(closeN(D72))))))
Rm = sp.Matrix([[1, 1], [0, 1]]); Lm = sp.Matrix([[1, 0], [1, 1]])
check("h and v act on F = p/q by [[1,1],[0,1]] and [[1,0],[1,1]] (det 1): they generate SL(2,Z), with S = R·L⁻¹·R",
      Rm.det() == 1 and Lm.det() == 1 and Rm*Lm.inv()*Rm == sp.Matrix([[0, 1], [-1, 0]]))
qy = sp.symbols('q')
Rq = lambda v: qy*v + 1; Lq = lambda v: qy*v/(qy*v + 1)
yfig = Rq(Rq(Lq(Rq(0))))
check("the q-step pair R_q: y ↦ qy+1, L_q: y ↦ qy/(qy+1) gives [5/2]_q = (1+2q+q²+q³)/(1+q) (Morier-Genoud–Ovsienko)",
      z(yfig - (1 + 2*qy + qy**2 + qy**3)/(1 + qy)))
_, Pf = build(cf_word([1, 1, 2]))
check("and the tangle's own y = A²g/f for 5/2 is exactly that q-rational at q = −A⁻⁴", z(yP(Pf) - yfig.subs(qy, q)))
check("[r/s]_q from the tangle → r/s at q = 1 (A = e^(iπ/4)) for all seven",
      all(sp.nsimplify(complex(sp.N(yP(build(cf_word(cf))[1]).subs(A, zeta))).real) == fracs[kk] for kk, cf in cfs.items()))

print("--- 5. The landscape: Möbius steps, but no bend ---")
F = lambda p, qb, aa, bb: (aa + bb + p*aa*bb)/(1 - qb*aa*bb)
mob = sp.Matrix([[1 + a*k, a], [u*a, 1]])  # (1+pa)b + a over (-q a) b + 1, with p=k, -q=u
check("every landscape law, 'add a' is a Möbius map of b: F(a,b) = ((1+pa)b + a)/(−qa·b + 1)",
      z(F(k, -u, a, b) - (mob[0, 0]*b + mob[0, 1])/(mob[1, 0]*b + mob[1, 1])))
check("tangle sums sit on the dial line (bend 0): x(S+T) = F_{d,0}(x(S), x(T))", z(F(dd, 0, x, y) - wadd(dd, x, y)))
wv = lambda y1, y2: 1/wadd(1/qy - 1, 1/y1, 1/y2)
check("the stacked law in y has identity ∞ (y ⋆ ∞ = y), so it is the dial seen through 1/y, not a landscape point",
      z(sp.limit(wv(x, w), w, sp.oo) - x))
check("vertical sum in y: 1/y = 1/y1 ∘_(1/q−1) 1/y2 (general pair, q = −A⁻⁴)",
      z(1/yP(vs(S, T)) - wadd(1/q - 1, 1/(A**2*g1/f1), 1/(A**2*g2/f2))))

print("--- 6. Torus knots: Alexander polynomials are ratios of q-numbers with q = t ---")
tt = sp.symbols('t')
qnt = lambda qv, n: sum(qv**i for i in range(n))
good = True; good1 = True; goodsym = True
for p_ in range(2, 8):
    for r_ in range(p_ + 1, 9):
        if math.gcd(p_, r_) != 1: continue
        Dl = sp.cancel((tt**(p_*r_) - 1)*(tt - 1)/((tt**p_ - 1)*(tt**r_ - 1)))
        good &= sp.Poly(Dl, tt).is_polynomial if hasattr(sp.Poly(Dl, tt), 'is_polynomial') else sp.denom(Dl) == 1
        good &= z(Dl - qnt(tt**p_, r_)/qnt(tt, r_)) and z(Dl - qnt(tt**r_, p_)/qnt(tt, p_))
        good1 &= Dl.subs(tt, 1) == 1
        deg = sp.degree(Dl, tt); goodsym &= z(sp.expand(tt**deg*Dl.subs(tt, 1/tt)) - Dl)
check("Δ_T(p,q) = (t^pq−1)(t−1)/((t^p−1)(t^q−1)) = [q]_(t^p)/[q]_t = [p]_(t^q)/[p]_t is a polynomial (coprime p<q≤8)", good)
check("Δ(1) = 1 and Δ is palindromic", good1 and goodsym)
check("[pq]_t = [p]_t·[q]_(t^p) for p, q ≤ 7", all(z(qnt(tt, p_*r_) - qnt(tt, p_)*qnt(tt**p_, r_)) for p_ in range(1, 8) for r_ in range(1, 8)))
good = True
for p_, r_ in [(2, 3), (2, 5), (3, 4), (3, 5), (2, 7), (4, 5)]:
    Dl = sp.cancel((tt**(p_*r_) - 1)*(tt - 1)/((tt**p_ - 1)*(tt**r_ - 1)))
    rts = [complex(rr) for rr in sp.Poly(Dl, tt).nroots()]
    want = [cmath.exp(2j*math.pi*j/(p_*r_)) for j in range(p_*r_) if j % p_ != 0 and j % r_ != 0]
    good &= len(rts) == (p_ - 1)*(r_ - 1) == len(want) and all(min(abs(r0 - w0) for w0 in want) < 1e-8 for r0 in rts)
check("zeros of Δ_T(p,q): the pq-th roots of unity that are neither p-th nor q-th roots, (p−1)(q−1) of them", good)
check("T(2,n), n odd: Δ = [n]_(−t) = 1 − t + t² − … + t^(n−1) (n = 3..15)",
      all(z(sp.cancel((tt**(2*n) - 1)*(tt - 1)/((tt**2 - 1)*(tt**n - 1))) - qnt(-tt, n)) for n in range(3, 16, 2)))
check("known values: Δ(3_1) = 1 − t + t², Δ(5_1) = 1−t+t²−t³+t⁴, Δ(T(3,4) = 8_19) = 1 − t + t³ − t⁵ + t⁶",
      z(qnt(-tt, 3) - (1 - tt + tt**2)) and z(qnt(-tt, 5) - (1 - tt + tt**2 - tt**3 + tt**4)) and
      z(sp.cancel((tt**12 - 1)*(tt - 1)/((tt**3 - 1)*(tt**4 - 1))) - (1 - tt + tt**3 - tt**5 + tt**6)))
check("|Δ(−1)| = determinant: 3 (3_1), 5 (5_1), 3 (8_19), and [n]_(−t) at t = −1 is n",
      qnt(-tt, 3).subs(tt, -1) == 3 and qnt(-tt, 5).subs(tt, -1) == 5 and
      abs(sp.cancel((tt**12 - 1)*(tt - 1)/((tt**3 - 1)*(tt**4 - 1))).subs(tt, -1)) == 3)
check("the twist tangle and the torus knot share the q: y([n]) = [n]_q with q = −t, Δ(T(2,n)) = [n]_(−t)", True and
      z(qn(q, 5).subs(A, tt**sp.Rational(-1, 4)) - qnt(-tt, 5)))

print("--- 7. Dial points of the Jones polynomial: loop value d, Potts Q = d² ---")
dsq = lambda tv: 2 + tv + 1/tv
pts = {'t=-1 (d=0, q=1)': (-1, 0), 't=e^(2πi/3) (d²=1)': (cmath.exp(2j*math.pi/3), 1), 't=i (d²=2)': (1j, 2),
       't=e^(iπ/3) (d²=3)': (cmath.exp(1j*math.pi/3), 3), 't=1 (d²=4)': (1, 4)}
check("Potts number Q = d² = 2 + t + 1/t: 0, 1, 2, 3, 4 at t = −1, e^(2πi/3), i, e^(iπ/3), 1",
      all(abs(dsq(v) - Q) < 1e-12 for v, Q in pts.values()))
ev = lambda e, tv: complex(sp.lambdify(t, e, 'numpy')(tv))
check("V(e^(2πi/3)) = 1 for all seven knots (Q = 1)", all(abs(ev(VK[kk], cmath.exp(2j*math.pi/3)) - 1) < 1e-9 for kk in VK))
check("V(i) = ±1 for all seven (Q = 2; the sign is the Arf invariant)", all(abs(abs(ev(VK[kk], 1j)) - 1) < 1e-9 and abs(ev(VK[kk], 1j).imag) < 1e-9 for kk in VK))
check("|V(e^(iπ/3))| is a power of √3 (Q = 3): 3_1 √3, 4_1 1, 6_1 √3", all(min(abs(abs(ev(VK[kk], cmath.exp(1j*math.pi/3))) - math.sqrt(3)**m) for m in range(4)) < 1e-9 for kk in VK)
      and abs(abs(ev(VK['3_1'], cmath.exp(1j*math.pi/3))) - math.sqrt(3)) < 1e-9 and abs(abs(ev(VK['4_1'], cmath.exp(1j*math.pi/3))) - 1) < 1e-9)
check("V(1) = 1 for every knot (Q = 4)", all(VK[kk].subs(t, 1) == 1 for kk in VK))
check("the twist's q at those points: q = −t = 1, −e^(2πi/3), −i, −e^(iπ/3), −1: orders 1, 6, 4, 3, 2",
      [next(m for m in range(1, 13) if abs((-v)**m - 1) < 1e-9) for v, _ in pts.values()] == [1, 6, 4, 3, 2])

import json
json.dump({kk: {str(e): int(cc) for e, cc in zip(*[[m[0] for m in sp.Poly(sp.expand(VK[kk]*t**20), t).monoms()],
           sp.Poly(sp.expand(VK[kk]*t**20), t).coeffs()])} for kk in VK}, open('jones.json', 'w'))
print(f"\n{ok}/{tot} checks pass")
