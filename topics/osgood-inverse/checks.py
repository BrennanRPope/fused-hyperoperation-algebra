"""Numerical checks for "The inverse of the Osgood test" (Nov 25, 2026).
Each line prints PASS/FAIL. Mult-add notation: E(x)=e^x-1, L(y)=log(1+y)=log1p (E's inverse)."""
import sys, json, math
import numpy as np
from scipy.integrate import solve_ivp, quad
sys.path.insert(0, '/mnt/project-files/math-theory/fractional-rungs')
from regular_iteration import Es

out = []
def rep(name, ok, extra=''):
    out.append(dict(name=name, ok=bool(ok), extra=extra))
    print(('PASS ' if ok else 'FAIL ') + name + (('  ' + extra) if extra else ''))

E, Lg = np.expm1, np.log1p
def Ln(y, n):
    for _ in range(n): y = np.log1p(y)
    return y
def fn(y, n):  # ladder f_n(y) = prod_{k<n} (1 + L^k y)
    p = 1.0
    for k in range(n): p = p * (1 + Ln(y, k))
    return p
def En(t, n):
    for _ in range(n): t = np.expm1(t)
    return t

# 1. Blow-up time is the Osgood integral, exactly (scalar autonomous ODE)
for p, Tstar in [(2.0, 1.0), (1.5, 2.0), (3.0, 0.5)]:
    I = quad(lambda s: (1 + s) ** -p, 0, np.inf)[0]
    ev = lambda t, y: y[0] - 1e6
    ev.terminal = True
    sol = solve_ivp(lambda t, y: [(1 + y[0]) ** p], [0, 3*Tstar], [0.0], events=ev, max_step=Tstar/50, rtol=1e-12, atol=1e-12)
    tb = sol.t_events[0][0]
    # time to reach 1e8 should equal T* minus tail integral from 1e8
    tail = (1 + 1e6) ** (1 - p) / (p - 1)
    rep(f"y'=(1+y)^{p}: blow-up time = integral = {Tstar}", abs(I - Tstar) < 1e-9 and abs(tb + tail - Tstar) < 1e-6,
        f"integral {I:.10f}, ODE reaches 1e6 at {tb:.8f} (+ tail {tail:.2e})")

# 2. The global ladder: E^n(t) solves y' = f_n(y); time to reach y is L^n(y)
for n in [1, 2, 3]:
    t = np.linspace(0.1, 1.2 if n == 3 else 2.0, 40); h = 1e-6
    d = (En(t + h, n) - En(t - h, n)) / (2 * h)
    ok = np.allclose(d, fn(En(t, n), n), rtol=1e-6)
    y = np.array([10.0, 1e3, 1e6])
    I = [quad(lambda s: 1 / fn(s, n), 0, yy, limit=400)[0] for yy in y]
    rep(f"ladder n={n}: E^{n}(t) solves y' = f_{n}(y) and the clock is L^{n}(y)", ok and np.allclose(I, Ln(y, n), rtol=1e-7))

# 3. Rung-shift invariance: int ds/((1+s)g(s)) = int dx/g(E(x))
for name, g in [('g=1+s', lambda s: 1 + s), ('g=(1+L(s))^1.3', lambda s: (1 + Lg(s)) ** 1.3), ('g=1+L(s)', lambda s: 1 + Lg(s))]:
    S = 1e6
    a = quad(lambda s: 1 / ((1 + s) * g(s)), 0, S, limit=400)[0]
    b = quad(lambda x: 1 / g(E(x)), 0, Lg(S), limit=400)[0]
    rep(f"rung shift: s=E(x) carries the test for f=(1+s)*{name} one rung down", abs(a - b) < 1e-8 * max(1, a), f"{a:.10f} vs {b:.10f}")

# 4. Edge at each level: f_n passes, f_n*(1+L^{n-1})^eps fails (n=2, eps=0.3): tail integrals
def clock(f, a, b):  # integral in the variable u = L^2(s) to avoid huge s: ds = (1+s)(1+L s) du
    return quad(lambda u: (1 + E(E(u))) * (1 + E(u)) / f(E(E(u))), a, b, limit=400)[0]
I_pass = [clock(lambda s: fn(s, 2), 0, U) for U in [3, 5, 6.5]]
I_fail = [clock(lambda s: fn(s, 2) * (1 + Lg(s)) ** 0.3, 0, U) for U in [3, 5, 6.5]]
exact_fail = [(1 - (1 + E(U)) ** -0.3) / 0.3 for U in [3, 5, 6.5]]
rep("n=2: f_2=(1+y)(1+Ly) keeps the clock running (= L^2 y); f_2*(1+Ly)^0.3 has a finite clock, 1/0.3",
    abs(I_pass[-1] - 6.5) < 1e-6 and np.allclose(I_fail, exact_fail, rtol=1e-7),
    f"pass {I_pass}, fail {[round(v,6) for v in I_fail]} -> limit {1/0.3:.6f}")

# 5. Uniqueness: non-Osgood modulus gives a fan of solutions out of 0; Osgood modulus does not
w = lambda r: 2 * np.sqrt(r)
Om = quad(lambda r: 1 / w(r), 0, 1)[0]
rep("x'=2*sqrt(x): leaving time Omega(1) = 1 is finite, so x=(t-c)^2 (t>=c) all solve from 0", abs(Om - 1) < 1e-9)
for c in [0.0, 0.5, 1.3]:
    t = np.linspace(c + 0.01, c + 2, 50); x = (t - c) ** 2
    d = 2 * (t - c)
    assert np.allclose(d, w(x))
wlog = lambda r: r * np.log(1 / r)
# substitute r = e^{-v}: dr/(r log(1/r)) = dv/v
tails = [quad(lambda v: np.exp(-v) / wlog(np.exp(-v)), math.log(2), k * math.log(10), limit=400)[0] for k in [5, 20, 80]]
rep("x'=x log(1/x) (Yudovich's modulus): the leaving time from r to 1/2 grows without bound as r -> 0",
    all(np.diff(tails) > 0.9) and abs(tails[-1] - (math.log(80 * math.log(10)) - math.log(math.log(2)))) < 1e-6,
    f"{[round(v,4) for v in tails]} (= loglog(1/r) - loglog 2)")

# 6. The mirror: s = 1/x turns x' = w(x) into s' = -s^2 w(1/s); leaving 0 = arriving from infinity
x = np.linspace(0.01, 2, 30); sol_x = x  # check the vector-field identity pointwise
s = 1 / x
rep("mirror: d/dt(1/x) = -s^2 w(1/s) when x' = w(x)", np.allclose(-w(x) / x ** 2, -s ** 2 * w(1 / s)))
Ttop = quad(lambda s: 1 / (s ** 2 * w(1 / s)), 1, np.inf)[0]
rep("mirror of 2*sqrt(x): s'=-2 s^{3/2} comes in from infinity in finite time, the same time Omega(1)=1", abs(Ttop - 1) < 1e-9)
# mult-add form: u = multAddInv(v) = -v/(1+v) swaps v -> infinity with u -> wall -1
v = np.array([1e2, 1e6, 1e12]); u = -v / (1 + v)
rep("mult-add inverse sends infinity to the wall: -v/(1+v) -> -1", np.allclose(u, -1, atol=1e-2) and abs(u[-1] + 1) < 1e-11)

# 7. Systems: same size, opposite outcomes (direction decides)
def rot(t, z): r2 = z[0] ** 2 + z[1] ** 2; return [-r2 * z[1], r2 * z[0]]
def rad(t, z): r2 = z[0] ** 2 + z[1] ** 2; return [r2 * z[0], r2 * z[1]]
sr = solve_ivp(rot, [0, 50], [1.0, 0.0], rtol=1e-11, atol=1e-12)
ev = lambda t, z: z[0] ** 2 + z[1] ** 2 - 1e10; ev.terminal = True
sd = solve_ivp(rad, [0, 5], [1.0, 0.0], events=ev, rtol=1e-11, atol=1e-12)
nr = np.hypot(sr.y[0], sr.y[1])
rep("rotation field |x|^2 Jx (speed r^3, fails Osgood) keeps |x| = 1 forever", np.allclose(nr, 1, atol=1e-7), f"max drift {abs(nr-1).max():.1e}")
rep("radial field |x|^2 x (same speed) blows up at T = 1/(2 r0^2) = 0.5", abs(sd.t_events[0][0] - 0.5) < 1e-4, f"reached r^2=1e10 at {sd.t_events[0][0]:.6f}")

# 8. Rung height: test int dx/E^t(x), computed in the rung flow's own coordinate x = E^a(x0)
#    (dx = V(E^a x0) da), which reaches x ~ 1e10 by a = 2.5 from x0 = 2.
x0, hh0 = 2.0, 1e-4
a = np.linspace(0, 2.6, 400)
Vg = np.array([(Es(ai + hh0, np.array([x0]))[0] - Es(ai - hh0, np.array([x0]))[0]) / (2 * hh0) for ai in a])
res = {}
for t in [0.0, 0.25, 0.5, 1.0]:
    with np.errstate(over='ignore', invalid='ignore'):
        Pt = np.array([Es(ai + t, np.array([x0]))[0] for ai in a])
        g = np.nan_to_num(Vg / Pt, nan=0.0, posinf=0.0)
    I = np.concatenate([[0], np.cumsum((g[1:] + g[:-1]) / 2 * np.diff(a))])
    res[t] = [I[np.searchsorted(a, v)] for v in [2.0, 2.5]]
rep("test int dx/E^t(x): rung 0 keeps growing (it is log x)", res[0.0][1] - res[0.0][0] > 10, f"{res[0.0][0]:.3f} -> {res[0.0][1]:.3f}")
rep("rungs 0.25, 0.5, 1: the integral has stopped growing by x ~ 1e10",
    all(res[t][1] - res[t][0] < 0.01 * res[t][1] for t in [0.25, 0.5, 1.0]),
    ', '.join(f't={t}: {res[t][0]:.5f} -> {res[t][1]:.5f}' for t in [0.25, 0.5, 1.0]))

# 9. The generator of the rung flow sits at the Osgood edge: V = d/dt E^t |_{t=0}
y = np.array([3.0, 10.0, 1e2, 1e4, 1e8, 1e16, 1e40, 1e100])
hh = 1e-4
V = (Es(hh, y) - Es(-hh, y)) / (2 * hh)
VE = (Es(hh, E(np.array([1.0, 2.0, 3.0, 5.0]))) - Es(-hh, E(np.array([1.0, 2.0, 3.0, 5.0])))) / (2 * hh)
V0 = (Es(hh, np.array([1.0, 2.0, 3.0, 5.0])) - Es(-hh, np.array([1.0, 2.0, 3.0, 5.0]))) / (2 * hh)
rep("generator law V(E(x)) = (1 + E(x)) V(x)", np.allclose(VE, (1 + E(np.array([1.0, 2.0, 3.0, 5.0]))) * V0, rtol=1e-5))
ratio = V / ((1 + y) * Lg(y))
rep("V(y)/((1+y)L(y)) keeps rising: V outgrows the n=2 rung of the ladder", np.all(np.diff(ratio) > 0), ' '.join(f'{r:.3g}' for r in ratio))
# Abel clock: int_x^{E(x)} dy/V = 1 exactly
x0 = 2.0
A = quad(lambda yy: 1 / ((Es(hh, np.array([yy])) - Es(-hh, np.array([yy]))) / (2 * hh))[0], x0, float(E(x0)), limit=200)[0]
rep("the rung flow's clock: int from x to E(x) of dy/V = 1 (one rung = one unit of time)", abs(A - 1) < 1e-4, f"{A:.6f}")

json.dump(out, open('/mnt/project-files/math-theory/osgood-inverse/checks.json', 'w'), indent=1)
print(sum(o['ok'] for o in out), '/', len(out), 'pass')
