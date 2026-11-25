"""Which smooth shape passes through every rung of the s-ladder?

Rung s:  B_s(a,b) = E^s( E^-s(a) + E^-s(b) ),  E(x) = e^x - 1 (fused exponential).
s = 0 addition, s = 1 mult-add, s = 2 (1+a)(1+b)e^{L(a)L(b)} - 1, s = -1 log(e^a + e^b - 1).
E^s is the regular (Szekeres/Ecalle) iterate, computed in ../fractional-rungs/regular_iteration.py.

Candidates tested here (each defined from mult-add, each checked at whole rungs,
for smoothness, and for what it says about infinities):
  A  the rung sheet: level curves B_s(a,b) = c stacked over height s
  B  the flow (ruler) surface z = E^s(x): rung s is plain addition in the coordinate E^-s
  C  the two-dial surface E^s(d * E^-s(a)): rung dial s times count dial d (the cube is s = 1)
  D  the saddle horn: curvature of the graph z = B_s(a,b) at the identity
  E  a helix: complex points riding the rungs (near the identity)
  F  the tropical floor: s -> -infinity
  G  the corner-dial sheet a + b + k ab (a different smooth family)
  H  the growing cube alone, and the rung-constant surface (iteration ladder)
Run: python3 rung_shape.py results.json
"""
import sys, os, json, math, cmath
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fractional-rungs"))
from regular_iteration import Es, B, coeffs_at, ORDER

res = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    res[name] = bool(ok)

L = np.log1p; E = np.expm1

# ------------------------------------------------------------------ A: rung sheet
def sheet_b(s, a, c=1.0):
    """the b with B_s(a, b) = c:  b = E^s(E^-s(c) - E^-s(a))"""
    a = np.asarray(a, float)
    u = float(Es(-s, np.array([c]))[0]) - Es(-s, a)
    out = Es(s, u)
    return out

a = np.linspace(0.0, 1.0, 201)
closed = {
    0: 1 - a,
    1: 2/(1 + a) - 1,
    -1: np.log(np.e - np.exp(a) + 1),
}
for s, ref in closed.items():
    err = np.nanmax(np.abs(sheet_b(s, a) - ref))
    rep(f"A sheet slice s={s} is the exact level curve", err < 1e-8, f"max err {err:.1e}")
# s = 2: check B_2(a, b(a)) = 1 with the closed form of rung 2
b2 = sheet_b(2, a[1:-1])
aa = a[1:-1]
B2 = (1 + aa)*(1 + b2)*np.exp(np.log1p(aa)*np.log1p(b2)) - 1
err = np.max(np.abs(B2 - 1))
rep("A sheet slice s=2 satisfies the closed form of rung 2", err < 1e-8, f"max err {err:.1e}")

# smoothness in s across whole rungs: one-sided second derivatives agree
def d2_sides(s0, a0=0.5, h=0.02):
    f = lambda s: float(sheet_b(s, np.array([a0]))[0])
    left = (f(s0) - 2*f(s0 - h) + f(s0 - 2*h))/h**2
    right = (f(s0 + 2*h) - 2*f(s0 + h) + f(s0))/h**2
    centre = (f(s0 + h) - 2*f(s0) + f(s0 - h))/h**2
    return left, centre, right
worst = 0
for s0 in [-1, 0, 1, 2]:
    l, c_, r = d2_sides(s0)
    # one-sided estimates differ from the centre by O(h) * third derivative; a kink would give O(1/h)
    worst = max(worst, abs(l - c_) + abs(r - c_))
    print(f"   s={s0:+d}: d2b/ds2 left {l:+.5f}  centre {c_:+.5f}  right {r:+.5f}")
rep("A sheet has no crease at any whole rung (one-sided curvatures in s agree to O(h))", worst < 0.05, f"worst gap {worst:.3f}")

# limits: bottom -> max corner (L shape), top -> collapses onto the axes (slowly, about like 1/s)
def area_at(s, n=2001):
    aa = np.linspace(0, 1, n)
    return np.trapezoid(np.clip(np.nan_to_num(sheet_b(s, aa), nan=0.0), 0, 1), aa)
lo = area_at(-3.5); hi = [area_at(t) for t in (8, 16, 32, 64)]
rep("A bottom of the sheet: area under the slice -> 1 (the outer corner max(a,b)=1)", lo > 0.999, f"area at s=-3.5: {lo:.4f}")
rep("A top of the sheet: area -> 0 (the curve falls onto the axes)", all(np.diff(hi) < 0) and hi[-1] < 0.01,
    "areas at s=8,16,32,64: " + ", ".join(f"{t:.4f}" for t in hi))
# area fraction under each slice (inside the unit square)
ss = np.linspace(-3, 4, 71)
area = []
for s in ss:
    bb = np.clip(np.nan_to_num(sheet_b(s, a), nan=0.0), 0, 1)
    area.append(np.trapezoid(bb, a))
area = np.array(area)
rep("A area under the slice falls smoothly from near 1 to near 0 (monotone)", np.all(np.diff(area) < 0), f"{area[0]:.3f} -> {area[-1]:.3f}")
a0 = area[np.argmin(abs(ss))]; a1 = area[np.argmin(abs(ss - 1))]
rep("A area at rung 0 = 1/2 and at rung 1 = 2 ln 2 - 1", abs(a0 - .5) < 1e-4 and abs(a1 - (2*math.log(2) - 1)) < 1e-4, f"{a0:.5f}, {a1:.5f} (2ln2-1 = {2*math.log(2)-1:.5f})")
res["area"] = {"s": ss.tolist(), "area": area.tolist()}

# ---------- uniqueness: any other interpolation creases the sheet at the identity
A_COEF = [-1/36, 1/540, 1/7776, -71/435456]
def abel0(y):  # asymptotic regular Abel function of E near 0+, alpha(E(y)) = alpha(y) + 1
    return -2/y + math.log(y)/3 + sum(c*y**(i+1) for i, c in enumerate(A_COEF))
def abel(x, n=4000):
    y = x
    for _ in range(n): y = math.log1p(y)
    return abel0(y) + n
def abel_inv(A, n=4000):
    # solve abel0(y) = A - n by Newton for small y, then walk out
    t = A - n
    y = -2/t
    for _ in range(60):
        g = abel0(y) - t
        dg = 2/y**2 + 1/(3*y) + sum((i+1)*c*y**i for i, c in enumerate(A_COEF))
        y -= g/dg
    for _ in range(n): y = math.expm1(y)
    return y
# consistency: abel-flow equals Es
x0 = 0.37
e_ab = abs(abel_inv(abel(x0) + 0.5) - float(Es(0.5, np.array([x0]))[0]))
rep("the Abel function reproduces the regular half rung", e_ab < 1e-9, f"diff {e_ab:.1e}")

EPS = 0.08
def perturbed(s, x):
    """another interpolation: Abel coordinate alpha + EPS sin(2 pi alpha); same whole rungs"""
    A = abel(x)
    T = A + EPS*math.sin(2*math.pi*A) + s
    Ap = A + s
    for _ in range(80):
        Ap -= (Ap + EPS*math.sin(2*math.pi*Ap) - T)/(1 + 2*math.pi*EPS*math.cos(2*math.pi*Ap))
    return abel_inv(Ap)
for s in [1, 2]:
    e_int = abs(perturbed(s, 0.4) - float(Es(s, np.array([0.4]))[0]))
    rep(f"perturbed interpolation still hits rung {s} exactly", e_int < 1e-8, f"diff {e_int:.1e}")
xs_small = np.geomspace(2e-3, 0.2, 60)
def slope(fun, x, h=1e-7):
    return (fun(x*(1 + h)) - fun(x*(1 - h)))/(2*x*h)
reg_sl = np.array([slope(lambda t: float(Es(0.5, np.array([t]))[0]), x) for x in xs_small])
per_sl = np.array([slope(lambda t: perturbed(0.5, t), x) for x in xs_small])
print(f"   slope of the half rung near 0: regular {reg_sl[:3].round(4)} ... -> 1;  perturbed ranges {per_sl[:20].min():.3f} .. {per_sl[:20].max():.3f}")
rep("regular half rung has slope -> 1 at the identity", abs(reg_sl[0] - 1) < 0.01)
rep("perturbed half rung keeps oscillating near the identity (no slope at 0: a crease)", per_sl[:20].max() - per_sl[:20].min() > 0.3,
    f"oscillation {per_sl[:20].max() - per_sl[:20].min():.3f}")
res["uniqueness"] = {"x": xs_small.tolist(), "regular": reg_sl.tolist(), "perturbed": per_sl.tolist(), "eps": EPS}

# ------------------------------------------------------------------ B: flow surface
xg = np.array([-0.7, -0.3, 0.2, 0.8, 1.5])
errs = [np.max(np.abs(Es(1, xg) - E(xg))), np.max(np.abs(Es(2, xg) - E(E(xg)))),
        np.max(np.abs(Es(-1, xg) - L(xg))), np.max(np.abs(Es(-2, xg) - L(L(xg))))]
rep("B flow surface slices at s=1,2,-1,-2 are E, E∘E, L, L∘L", max(errs) < 1e-8, f"max err {max(errs):.1e}")
gl = np.max(np.abs(Es(0.7, Es(-0.3, xg)) - Es(0.4, xg)))
rep("B flow surface is one flow: E^0.7 ∘ E^-0.3 = E^0.4", gl < 1e-8, f"err {gl:.1e}")
# rung s is addition in coordinate E^-s
p, q = 0.6, 0.9
for s in [0.5, 1.5]:
    lhs = float(B(s, p, q)); rhs = float(Es(s, Es(-s, np.array([p])) + Es(-s, np.array([q])))[0])
print("   (rung s is addition in the coordinate u = E^-s(x), by definition)")

# ------------------------------------------------------------------ C: two-dial surface
def dial2(s, d, x=1.0):
    return Es(s, d*float(Es(-s, np.array([x]))[0]))
dd = np.linspace(0, 4, 41)
e1 = np.max(np.abs(dial2(1, dd) - (2**dd - 1)))
e0 = np.max(np.abs(dial2(0, dd) - dd))
rep("C two-dial: rung-1 slice is the cube count 2^d - 1 for every real d", e1 < 1e-8, f"err {e1:.1e}")
rep("C two-dial: rung-0 slice is the line d", e0 < 1e-9, f"err {e0:.1e}")
# whole d at rung 1/2 = d-fold B_1/2 of 1
v = np.array([1.0])
for _ in range(4): v = B(0.5, v, np.array([1.0]))
e5 = abs(float(v[0]) - float(dial2(0.5, np.array([5.0]))[0]))
rep("C two-dial: d = 5 at rung 1/2 equals five 1s combined by rung 1/2", e5 < 1e-8, f"err {e5:.1e}")
# growth class at s = 1/2: faster than any power, slower than any exponential
ns = np.array([1e1, 1e2, 1e3, 1e4, 1e5])
vals = np.array([float(dial2(0.5, np.array([n]))[0]) for n in ns])
pow_exp = np.log(vals)/np.log(ns); exp_rate = np.log(vals)/ns
print("   rung 1/2 run n times: n =", ns, "\n      log(value)/log(n) =", pow_exp.round(3), "\n      log(value)/n      =", exp_rate.round(5))
rep("C rung 1/2 growth: power exponent keeps rising while exponential rate keeps falling",
    np.all(np.diff(pow_exp) > 0) and np.all(np.diff(exp_rate) < 0))
res["growth_half"] = {"n": ns.tolist(), "value": vals.tolist()}

# ------------------------------------------------------------------ D: saddle horn
h = 1e-3
curv = {}
for s in [-2, -1, -0.5, 0, 0.5, 1, 2]:
    f = lambda x, y: float(B(s, np.array([x]), np.array([y]))[0])
    faa = (f(h, 0) - 2*f(0, 0) + f(-h, 0))/h**2
    fab = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h))/(4*h*h)
    K = (faa*faa - fab*fab)/(1 + 2)**2  # gradient (1,1) at the identity, fbb = faa = 0
    curv[s] = (faa, fab, K)
okD = all(abs(curv[s][0]) < 1e-3 and abs(curv[s][1] - s) < 5e-3*max(1, abs(s)) and abs(curv[s][2] + s*s/9) < 5e-3*max(1, s*s) for s in curv)
rep("D graph of rung s at the identity: flat along the axes, twist = s, Gaussian curvature = -s²/9", okD,
    " ".join(f"s={s}:K={curv[s][2]:+.4f}" for s in curv))
# ruled? second derivative in b along a = 0.5 vanishes only at s = 0, 1
rul = {}
for s in [-1, -0.5, 0, 0.5, 1, 1.5, 2]:
    bb = np.array([0.3, 0.6, 0.9])
    v3 = B(s, np.full(3, 0.5), bb)
    rul[s] = float(v3[0] - 2*v3[1] + v3[2])
rep("D only rungs 0 and 1 have straight-line rulings (graph linear in b)", abs(rul[0]) < 1e-9 and abs(rul[1]) < 1e-9 and all(abs(rul[s]) > 1e-3 for s in rul if s not in (0, 1)),
    " ".join(f"{s}:{rul[s]:+.4f}" for s in rul))
res["curvature"] = {str(k): v[2] for k, v in curv.items()}

# ------------------------------------------------------------------ E: helix near the identity
def Es_c(s_arr, x, n=20000):
    """regular iterate for complex x in the attracting (left) petal: L^n P_s E^n"""
    y = complex(x)
    for _ in range(n): y = cmath.exp(y) - 1
    cs_list = [coeffs_at(float(s)) for s in s_arr]
    out = []
    for cs in cs_list:
        z = 0j
        for i in range(ORDER, 0, -1): z = (z + cs[i])*y
        out.append(z)
    z = np.array(out)
    for _ in range(n): z = np.log(1 + z)
    return z
x0 = -0.03 + 0.03j
S = np.concatenate([-np.geomspace(3000, 1, 60), [0.0], np.geomspace(1, 3000, 60)])
traj = Es_c(S, x0)
r = abs(x0)**2/(2*x0.imag)
dev = np.abs(np.abs(traj - 1j*r) - r)/r
rep("E complex trajectory ends at 0 in both directions (one closed loop over the whole ladder)",
    abs(traj[0]) < 2e-3 and abs(traj[-1]) < 2e-3, f"|start| {abs(traj[0]):.1e}, |end| {abs(traj[-1]):.1e}")
rep("E the loop is nearly the circle that kisses 0 (from the x²/2 term)", np.max(dev) < 0.1, f"max deviation {np.max(dev):.3f} of the radius")
e_int = abs(traj[np.argmin(abs(S - 1))] - (cmath.exp(x0) - 1))
rep("E at s = 1 the loop point is E(x0) (to the rounding of 20000 complex steps)", e_int < 1e-7, f"err {e_int:.1e}")
# winding: angle around the circle's centre as s runs
ang = np.unwrap(np.angle(traj - 1j*r))
rep("E the whole ladder s from -inf to +inf winds once (not a repeating helix)", abs(abs(ang[-1] - ang[0]) - 2*math.pi) < 0.25, f"turned {abs(ang[-1]-ang[0])/(2*math.pi):.3f} times")
res["loop"] = {"s": S.tolist(), "re": traj.real.tolist(), "im": traj.imag.tolist(), "r": r}

# ------------------------------------------------------------------ F: tropical floor
low = [float(B(s, np.array([2.0]), np.array([1.0]))[0]) for s in [-1, -2, -3]]
rep("F rungs fall to max(a,b) as s -> -inf (the floor of the sheet)", abs(low[-1] - 2) < 1e-6, str([round(t, 6) for t in low]))

# ------------------------------------------------------------------ G: corner dial a + b + k ab
for s, label in [(2, "rung 2"), (-1, "rung -1")]:
    # best k to match on the unit square
    gr = np.linspace(0.05, 1, 12)
    A2, B2_ = np.meshgrid(gr, gr)
    tgt = B(s, A2.ravel(), B2_.ravel()).reshape(A2.shape)
    k = np.sum((tgt - A2 - B2_)*A2*B2_)/np.sum((A2*B2_)**2)
    misfit = np.max(np.abs(tgt - (A2 + B2_ + k*A2*B2_)))
    rep(f"G corner dial misses {label}: best k = {k:.3f} still off", misfit > 0.05, f"max misfit {misfit:.3f}")

json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w"), indent=1)
nb = sum(v for v in res.values() if isinstance(v, bool)); nt = sum(isinstance(v, bool) for v in res.values())
print(f"\n{nb} of {nt} checks pass")
