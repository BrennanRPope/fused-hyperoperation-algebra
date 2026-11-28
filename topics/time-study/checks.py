"""DRAFT time-study checks. Units c = 1. Weak gravitational wave along z, linear order in h.
Question: does one shared time axis with a winding rate per place (the lapse N) describe a passing wave?"""
import json, numpy as np, sympy as sp
R = {}; ok = lambda name, cond, info="": (R.__setitem__(name, {"pass": bool(cond), "info": info}), print(("PASS " if cond else "FAIL ")+name, info))

t, x, y, z, eps = sp.symbols('t x y z epsilon', real=True)
hp, hc = sp.Function('hp'), sp.Function('hc')       # plus and cross strains, functions of u = t - z
u = t - z
X = [t, x, y, z]

def lin_riemann(h, a, b, c, d):
    """Linearised Riemann R_abcd from perturbation h (4x4 sympy matrix)."""
    D = lambda i, j, M: sp.diff(M, X[i], X[j])
    return sp.Rational(1, 2)*(D(b, c, h[a, d]) + D(a, d, h[b, c]) - D(a, c, h[b, d]) - D(b, d, h[a, c]))

# --- T1: TT gauge. Lapse is 1 everywhere, shift 0: one shared axis, every free clock winds at the same rate.
hTT = sp.zeros(4); hTT[1, 1] = hp(u); hTT[2, 2] = -hp(u); hTT[1, 2] = hTT[2, 1] = hc(u)
eta = sp.diag(-1, 1, 1, 1); g = eta + eps*hTT
ginv = g.inv()
Gam_i00 = [sp.simplify(sp.series(sum(ginv[i, l]*(2*sp.diff(g[l, 0], t) - sp.diff(g[0, 0], X[l]))/2 for l in range(4)), eps, 0, 2).removeO()) for i in range(1, 4)]
ok("T1 TT gauge: lapse N = 1, shift 0, test masses at rest stay at rest (Gamma^i_00 = 0)",
   g[0, 0] == -1 and all(g[0, i] == 0 for i in range(1, 4)) and all(G == 0 for G in Gam_i00))

# --- T2: yet the measured light round trip along an x-arm changes: the wave lives in the spatial stretch, not in N.
w, L, A = 2*np.pi, 0.1, 1e-3            # wave period 1, arm 0.1 wavelengths
hfun = lambda tt: A*np.cos(w*tt)
def round_trip(t0, n=20000):            # integrate dx/dt = +-1/sqrt(1+h) at z = 0
    tt, xx, dt = t0, 0.0, 2.2*L/n; out = True
    while True:
        v = 1/np.sqrt(1 + hfun(tt)); xn = xx + (v if out else -v)*dt
        if out and xn >= L: tt += (L - xx)/v; xx = L; out = False; continue
        if not out and xn <= 0: tt += xx/v; return tt - t0
        xx = xn; tt += dt
t0s = np.linspace(0, 1, 9)[:-1]
num = np.array([round_trip(t0) - 2*L for t0 in t0s])
formula = np.array([(A/(2*w))*(np.sin(w*(t0 + 2*L)) - np.sin(w*t0)) for t0 in t0s])   # +1/2 * integral of h
ok("T2 light round trip changes although N = 1 (numeric ray vs +1/2 integral of h)",
   np.max(np.abs(num - formula)) < 2e-3*A*L, f"max signal {np.max(np.abs(num)):.2e}, max error {np.max(np.abs(num-formula)):.1e}")

# --- T3: detector frame. Winding rate per place AND moment, N^2 = 1 + 2 Phi, Phi = -1/4 hdd_ij x^i x^j,
#         reproduces the wave's tidal pull exactly (geodesic deviation).
hdp, hdc = sp.symbols('hddot_p hddot_c', real=True)
Phi = -sp.Rational(1, 4)*(hdp*(x**2 - y**2) + 2*hdc*x*y)
force = [-sp.diff(Phi, v) for v in (x, y)]
RTT = {(i, j): lin_riemann(hTT, 0, i, 0, j) for i in (1, 2) for j in (1, 2)}
sub = {sp.diff(hp(u), t, 2): hdp, sp.diff(hc(u), t, 2): hdc}
def s(e): return sp.simplify(e.subs({sp.Derivative(hp(u), t, t): hdp}).doit().subs(sub)).subs({z: 0})
dev = [-(s(RTT[(1, 1)])*x + s(RTT[(1, 2)])*y), -(s(RTT[(2, 1)])*x + s(RTT[(2, 2)])*y)]
R00 = [s(RTT[(1, 1)]), s(RTT[(1, 2)])]
ok("T3 tidal force from the winding-rate field = geodesic deviation from the wave",
   all(sp.simplify(force[k] - dev[k]) == 0 for k in range(2)), f"R_0x0x = {R00[0]}, R_0x0y = {R00[1]}")

# --- T4: a winding rate fixed in time cannot do it: best static pattern over a cycle is zero, residual 100%.
ts = np.linspace(0, 1, 400, endpoint=False); hdd = -w**2*A*np.cos(w*ts)
best_static = hdd.mean(); resid = np.sqrt(np.mean((hdd - best_static)**2))/np.sqrt(np.mean(hdd**2))
ok("T4 static per-place winding rate fails for the wave (best fit 0, residual 100%)",
   abs(best_static) < 1e-12*w**2*A and abs(resid - 1) < 1e-9, f"residual fraction {resid:.3f}")

# --- T5: how good is 'shared axis + winding field' alone? Exact arm response / winding-only response = sinc(wL).
def ratio(wL): return np.sin(wL)/wL
cases = {"LIGO 4 km, 100 Hz": 2*np.pi*100*4e3/3e8, "LIGO 4 km, 1 kHz": 2*np.pi*1e3*4e3/3e8,
         "LISA 2.5 Gm, 10 mHz": 2*np.pi*0.01*2.5e9/3e8, "LISA 2.5 Gm, 0.1 Hz": 2*np.pi*0.1*2.5e9/3e8}
errs = {k: 1 - ratio(v) for k, v in cases.items()}
wl = 0.3; ok("T5 winding-only error is (wL)^2/6 for small arms", abs((1 - ratio(wl)) - wl**2/6) < wl**4/100,
   "; ".join(f"{k}: {100*e:.3g}%" for k, e in errs.items()))

# --- T6: spin 2. The winding pattern flips sign under a 90 degree turn and returns after 180; 45 degrees turns plus into cross.
def pattern(px, cx, xx, yy): return -0.25*(px*(xx**2 - yy**2) + 2*cx*xx*yy)
rng = np.random.default_rng(1); P = rng.normal(size=(2, 50))
rot = lambda a: (np.cos(a)*P[0] - np.sin(a)*P[1], np.sin(a)*P[0] + np.cos(a)*P[1])
base = pattern(1, 0, *P)
ok("T6 spin 2: 90 deg flips, 180 deg returns, 45 deg turns plus into cross",
   np.allclose(pattern(1, 0, *rot(np.pi/2)), -base) and np.allclose(pattern(1, 0, *rot(np.pi)), base)
   and np.allclose(pattern(1, 0, *rot(np.pi/4)), pattern(0, -1, *P)))

# --- T7: the law's helix reading per place. x = i v, 1 + x o x = (1 + i v)^2 (k = 1): real part 1 - v^2 is the clock
#         rate squared relative to the fastest clock, imag/real is turns per loop. In the wave both swing in time.
k = 1; comp = lambda a, b: a + b + k*a*b
gx = np.linspace(-1, 1, 41); GX, GY = np.meshgrid(gx, gx); hA = 0.05
N2 = 1 + 2*pattern(hA, 0, GX, GY); Nref2 = N2.max(); v = np.sqrt(1 - N2/Nref2)
lift = 1 + comp(1j*v, 1j*v)
ok("T7 law per place: Re(1 + iv o iv) = clock rate squared, works with a reference clock (the fastest one)",
   np.allclose(lift.real, N2/Nref2) and np.allclose(lift.imag, 2*v), f"reference clock rate^2 {Nref2:.4f}")

# --- T8: the slicing is a choice, the curvature is not. TT (N = 1) and detector frame (N varies) give the same R_0i0j.
hDet = sp.zeros(4); hDet[0, 0] = -2*Phi.subs({hdp: sp.diff(hp(t), t, 2), hdc: sp.diff(hc(t), t, 2)})
Rdet = sp.simplify(lin_riemann(hDet, 0, 1, 0, 1).subs({x: 0, y: 0}))
RTTo = sp.simplify(lin_riemann(hTT, 0, 1, 0, 1).subs(z, 0).doit())
target = -sp.Rational(1, 2)*sp.diff(hp(t), t, 2)
ok("T8 same tidal curvature in both slicings (frame-free content is how the winding rate bends, not its value)",
   sp.simplify(Rdet - target) == 0 and sp.simplify(RTTo - target) == 0, f"detector {Rdet}, TT {RTTo}")

json.dump({"checks": R, "lapse_only_errors": errs}, open(__file__.replace("checks.py", "results.json"), "w"), indent=1)
print(f"{sum(r['pass'] for r in R.values())} of {len(R)} pass")
