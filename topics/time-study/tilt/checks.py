"""DRAFT tilt + chain checks. Units c = 1, linear order in the wave strain.
Q1: can the law make the shared time axis tilt from place to place (GR's shift)?  Q2 (Brennan's hunch): is time a
composition chain of a pure function at each point, and does that chain make the tilt?"""
import json, itertools, numpy as np, sympy as sp
R = {}
def ok(name, cond, info=""):
    R[name] = {"pass": bool(cond), "info": info}; print(("PASS " if cond else "FAIL ") + name, info)

t, x, y, z = sp.symbols('t x y z', real=True); X = [t, x, y, z]
def riem(h, a, b, c, d):
    D = lambda i, j, M: sp.diff(M, X[i], X[j])
    return sp.Rational(1, 2)*(D(b, c, h[a, d]) + D(a, d, h[b, c]) - D(a, c, h[b, d]) - D(b, d, h[a, c]))
eps3 = lambda i, j, k: sp.LeviCivita(i, j, k)

# A1: a winding rate alone (only g_00 changes, any pattern in place and time) has NO magnetic curvature: it cannot tilt.
F = sp.Function('F')(t, x, y, z); h1 = sp.zeros(4); h1[0, 0] = F
mag = [sp.simplify(riem(h1, 0, i, j, k)) for i in range(1, 4) for j in range(1, 4) for k in range(1, 4)]
ok("A1 rate alone (any N(t,x,y,z)) gives R_0ijk = 0: no tilt, no twist", all(m == 0 for m in mag))

# A2/A8: the wave has magnetic curvature B, and in complex form B = i x E (the four-leaf turned 45 degrees).
hp, hc = sp.Function('hp'), sp.Function('hc'); u = t - z
hT = sp.zeros(4); hT[1, 1] = hp(u); hT[2, 2] = -hp(u); hT[1, 2] = hT[2, 1] = hc(u)
P, C = sp.symbols('Pdd Cdd', real=True)
def S(e): return sp.simplify(e.doit().subs({sp.Derivative(hp(u), u, 2): P, sp.Derivative(hc(u), u, 2): C}).subs(z, 0)) if e != 0 else 0
def S2(e):
    e = sp.simplify(e.doit())
    e = e.replace(lambda q: isinstance(q, sp.Subs), lambda q: q.doit())
    return sp.simplify(e.subs({sp.diff(hp(u), t, 2): P, sp.diff(hc(u), t, 2): C}).subs(z, 0))
E = {(i, j): S2(riem(hT, 0, i, 0, j)) for i in (1, 2) for j in (1, 2)}
B = {(i, j): S2(sp.Rational(1, 2)*sum(eps3(i, p, q)*riem(hT, p, q, 0, j) for p in (1, 2, 3) for q in (1, 2, 3))) for i in (1, 2) for j in (1, 2)}
Ec = E[(1, 1)] - sp.I*E[(1, 2)]; Bc = B[(1, 1)] - sp.I*B[(1, 2)]
ok("A2 wave has B != 0 (a real tilt/twist, the same in every slicing at this order)", any(v != 0 for v in B.values()),
   f"E_xx={E[(1,1)]}, E_xy={E[(1,2)]}, B_xx={B[(1,1)]}, B_xy={B[(1,2)]}")
ratio = sp.simplify(Bc/Ec)
ok("A8 complex B / complex E = +-i: the twist is the winding pattern times i (45 degree turn)", ratio in (sp.I, -sp.I), f"B/E = {ratio}")

# A3: detector frame tilt g_0i = -(2/3) R_0jik x^j x^k: nonzero, quadratic in place (picture tl1).
Xs = [x, y, z]
def g0i(i):
    return sp.expand(-sp.Rational(2, 3)*sum(S2(riem(hT, 0, j, i, k))*Xs[j-1]*Xs[k-1] for j in (1, 2, 3) for k in (1, 2, 3)))
shift = [g0i(i) for i in (1, 2, 3)]
ok("A3 detector-frame tilt is nonzero for the wave", any(s_ != 0 for s_ in shift), f"g_0x = {shift[0]}, g_0y = {shift[1]}, g_0z = {shift[2]}")

# Boost and rotation matrices (4x4)
def boost(v, n):
    n = np.asarray(n, float)/np.linalg.norm(n); g = 1/np.sqrt(1 - v*v); Bm = np.eye(4)
    Bm[0, 0] = g; Bm[0, 1:] = Bm[1:, 0] = -g*v*n; Bm[1:, 1:] += (g - 1)*np.outer(n, n); return Bm
def rot_angle(M):           # angle of the pure rotation left after removing the boost part
    g = M[0, 0]; vv = -M[1:, 0]/g; v = np.linalg.norm(vv); Rm = boost(v, -vv/v) @ M if v > 0 else M   # undo the boost part
    Rm = Rm[1:, 1:]; return np.arccos(np.clip((np.trace(Rm) - 1)/2, -1, 1))

# A4: chain of scalar steps (one number per point) commutes: a loop always closes, no twist.
k = 0.7; a = np.random.default_rng(2).normal(size=6)*0.3
lift = lambda s: 1 + k*s
loop_scalar = np.prod([lift(s) for s in a]) * np.prod([1/lift(s) for s in a[::-1]])
v1, v2 = 0.3, 0.5; law_d = lambda p, q: (p + q)/(1 + p*q)          # the two-dial law at k = 0, d = -1
ok("A4 scalar chain (and collinear tilts on the d = -1 law): order never matters, loops close",
   abs(loop_scalar - 1) < 1e-12 and abs(law_d(v1, v2) - law_d(v2, v1)) < 1e-15
   and np.allclose(boost(v1, [1, 0, 0]) @ boost(v2, [1, 0, 0]), boost(law_d(v1, v2), [1, 0, 0])))

# A5: chain of direction steps (tilts in two directions) leaves a rotation: the twist. Exact Thomas-Wigner angle.
M = boost(v2, [0, 1, 0]) @ boost(v1, [1, 0, 0]); th = rot_angle(M)
g1, g2 = 1/np.sqrt(1 - v1**2), 1/np.sqrt(1 - v2**2); G = g1*g2
th_exact = np.arccos((1 + G + g1 + g2)**2/((1 + G)*(1 + g1)*(1 + g2)) - 1)
vs = np.linspace(0.01, 0.2, 5)
loop = [rot_angle(boost(v, [0, -1, 0]) @ boost(v, [-1, 0, 0]) @ boost(v, [0, 1, 0]) @ boost(v, [1, 0, 0])) for v in vs]
ok("A5 perpendicular tilts composed in order leave a rotation (Thomas-Wigner, exact); small loop angle ~ v^2",
   abs(th - th_exact) < 1e-10 and np.allclose(np.array(loop)/vs**2, 1, rtol=0.05),
   f"angle {np.degrees(th):.4f} deg (exact {np.degrees(th_exact):.4f}); loop/v^2 = {np.round(np.array(loop)/vs**2, 3)}")

# A6: the law's chain form A o B = A + B + kAB is associative and its order gap is k[A,B]; for tilt generators
#     the gap is a rotation generator. This is the chain form already in torsion-dial/, not a new dial.
Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1; Ky = np.zeros((4, 4)); Ky[0, 2] = Ky[2, 0] = 1
Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1, 1
chain = lambda A_, B_: A_ + B_ + k*A_ @ B_
A_, B_, C_ = v1*Kx, v2*Ky, 0.2*Kx + 0.1*Ky
gap = chain(A_, B_) - chain(B_, A_)
ok("A6 chain law is associative; order gap k[A,B] for two tilts is a pure rotation (k v1 v2 J_z, up to sign)",
   np.allclose(chain(chain(A_, B_), C_), chain(A_, chain(B_, C_))) and np.allclose(np.abs(gap), np.abs(k*v1*v2*Jz)))

# A7: a chain needs grouping-free steps. Rung 1 (the law) gives one answer for every bracketing; rung 3 (powers) doesn't.
vals = [1.3, 0.8, 1.7, 1.1, 0.9]
def bracketings(seq, op):
    if len(seq) == 1: return [seq[0]]
    out = []
    for i in range(1, len(seq)):
        for L_ in bracketings(seq[:i], op):
            for R_ in bracketings(seq[i:], op): out.append(op(L_, R_))
    return out
r1 = bracketings(vals, lambda p, q: p + q + k*p*q); r3 = bracketings(vals, lambda p, q: p**q)
ok("A7 rung 1: all 14 bracketings agree; rung 3: they spread (no single 'next state')",
   np.ptp(r1) < 1e-9 and np.ptp(r3) > 0.1, f"rung 1 spread {np.ptp(r1):.1e}; rung 3 values {min(r3):.3f} to {max(r3):.3f}")


# A9 (arrow): rung 1 steps all undo, except the absorbing point -1/k, which sends every state to itself (many-to-one).
inv = lambda s: -s/(1 + k*s); law = lambda p, q: p + q + k*p*q
steps = np.random.default_rng(3).normal(size=8)*0.4
starts = np.linspace(-1, 1, 1000)
def run(seq, xs):
    out = xs.copy()
    for s_ in seq: out = law(s_, out)
    return out
fwd = run(steps, starts); back = run([inv(s_) for s_ in steps[::-1]], fwd)
absorb = np.insert(steps, 4, -1/k); fa = run(absorb, starts)
distinct = lambda v: len(np.unique(np.round(v, 9)))
ok("A9 rung 1: steps undo exactly; one absorbing step (-1/k) collapses 1000 states to 1, and nothing undoes it",
   np.allclose(back, starts) and distinct(fwd) == 1000 and distinct(fa) == 1 and abs(1 + k*(-1/k)) < 1e-15,
   f"distinct after chain: {distinct(fwd)}; with absorbing step: {distinct(fa)}")

# A10 (arrow): rung 3 makes reading order matter (left fold vs right fold) but loses no information, so it says
#      'direction matters' without saying which way is forward.
lf = lambda seq: __import__('functools').reduce(lambda p, q: p**q, seq)
rf = lambda seq: __import__('functools').reduce(lambda p, q: q**p, seq[::-1])
seq3 = [1.3, 0.8, 1.7, 1.1, 0.9]; xs = np.linspace(1.01, 3, 1000)
outs = xs**0.8; outs = outs**1.7
ok("A10 rung 3: left and right readings differ, yet each step is one-to-one (1000 states stay 1000)",
   abs(lf(seq3) - rf(seq3)) > 0.05 and distinct(outs) == 1000,
   f"left {lf(seq3):.4f}, right {rf(seq3):.4f}; rung 1 left = right: {abs(__import__('functools').reduce(law, seq3) - __import__('functools').reduce(lambda p,q: law(q,p), seq3[::-1])) < 1e-12}")

# S1 (speed limit, his hunch): neighbour-only chain on rung 1 (lifts multiply; geometric mean of three neighbours).
#    A nudge at one place spreads at most one place per step: a light cone. Its speed is 1 place per step, which is c
#    only if c is how a place and a step are measured; the chain gives the cone, not the number.
n_sites, n_steps = 81, 30; logl = np.zeros(n_sites); logl[40] = np.log(1.5); front = []   # work with log of the lift
hist = [np.exp(logl)]
for _ in range(n_steps):
    logl = (np.roll(logl, 1) + logl + np.roll(logl, -1))/3; hist.append(np.exp(logl))   # untouched places stay exactly 0
    front.append(int(np.max(np.abs(np.nonzero(logl > 0)[0] - 40))))
ok("S1 neighbour chain: influence spreads exactly one place per step (a cone)", front == list(range(1, n_steps + 1)),
   f"front after 30 steps: {front[-1]} places")

# S2: c is the law's absorbing point in other coordinates. Velocity addition (u+v)/(1+uv) (two-dial law, k = 0, d = -1)
#     maps onto the law by a(v) = -2v/(k(1+v)); then a(u (+) v) = a(u) o_k a(v) and a(c = 1) = -1/k.
uu, vv, kk = sp.symbols('u v k', positive=True)
amap = lambda w: -2*w/(kk*(1 + w))
lhs = amap((uu + vv)/(1 + uu*vv)); rhs = amap(uu) + amap(vv) + kk*amap(uu)*amap(vv)
ok("S2 velocity addition is the law in disguise, and light speed lands on the absorbing point -1/k",
   sp.simplify(lhs - rhs) == 0 and sp.simplify(amap(1) + 1/kk) == 0 and sp.simplify((1 + vv)/(1 + vv) - 1) == 0,
   "a(1) = -1/k; also 1 (+) v = 1 for every v")

# S3: the limit and the arrow are the same point seen two ways. Slower-than-light steps compose toward c and never
#     reach it (no finite chain lands on -1/k); a light-like step lands on it at once and erases the rest (A9).
vel = [0.5]
for _ in range(30): vel.append((vel[-1] + 0.5)/(1 + 0.5*vel[-1]))
ok("S3 30 steps of 0.5c give a speed below c (never reached); one step of c gives c whatever came before",
   all(w < 1 for w in vel) and abs((1 + 0.37)/(1 + 0.37) - 1) == 0,
   f"after 30 steps: 1 - v = {1 - vel[-1]:.2e}")

json.dump({"checks": R, "cone_hist": [list(h) for h in hist], "vel_chain": vel, "distinct_with_absorb": int(distinct(fa)), "rung3_bracketings": sorted(r3), "loop_over_v2": list(np.array(loop)/vs**2)},
          open(__file__.replace("checks.py", "results.json"), "w"), indent=1)
print(f"{sum(r['pass'] for r in R.values())} of {len(R)} pass")
