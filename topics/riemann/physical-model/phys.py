"""A physical model for the zeros' heat flow dz_j/dt = 2 sum 1/(z_j - z_k) (Rodgers-Tao normalisation).
Tests: (1) Burgers poles via Cole-Hopf (exact, on and off the line); (2) Calogero-Moser second-order form and its
zero energy; (3) the 2D dipole potential on the line; (4) why no Euclidean-gradient system matches off the line;
(5) collapse laws of 1/r (line charges, superconductor vortices) vs 1/r^2 (point charges on a wire) beads;
(6) the threshold as a Hill radius.
"""
import json, math
import numpy as np
from numpy.polynomial import polynomial as Pn
from scipy.integrate import solve_ivp
out = {}
rng = np.random.default_rng(11)

def F(z):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf); return (1 / D).sum(1)
def flow(rot=1.0):
    def g(t, v):
        n = len(v) // 2; z = v[:n] + 1j * v[n:]; w = rot * 2 * F(z); return np.concatenate([w.real, w.imag])
    return g
pack = lambda z: np.concatenate([z.real, z.imag]); unpack = lambda v: v[:len(v)//2] + 1j * v[len(v)//2:]

# ---- 1. Burgers: phi_s = phi_xx (s = -t); u = -2 phi_x/phi; poles of u = zeros of phi move dz/ds = -2 sum 1/(z-w)
z0 = np.array([-3.0, -1.2, 0.4, 2.5, 0.9 + 0.8j, 0.9 - 0.8j])
c0 = np.real(Pn.polyfromroots(z0))
def heat_poly(c, s):          # exact: e^{s D^2} on a polynomial = sum s^k/k! D^{2k}
    outc = np.zeros_like(c); term = c.copy(); k = 0
    while np.any(term != 0):
        outc[:len(term)] += term * s**k / math.factorial(k)
        term = Pn.polyder(term, 2) if len(term) > 2 else np.zeros(1); k += 1
        if len(term) == 1 and term[0] == 0: break
    return outc
s1 = 0.15
roots_exact = np.sort_complex(Pn.polyroots(heat_poly(c0, s1)))
sol = solve_ivp(flow(-1.0), [0, s1], pack(z0.astype(complex)), rtol=1e-12, atol=1e-14, method="DOP853")
roots_ode = np.sort_complex(unpack(sol.y[:, -1]))
out["burgers_poles"] = {"max_err": float(np.max(np.abs(roots_exact - roots_ode))), "s": s1}
# Burgers velocity profile as a conjugate pole pair approaches the axis (de Bruijn forward = Burgers backward)
xs = np.linspace(-1.5, 3.5, 1200)
prof = {}
for y in [0.8, 0.4, 0.2, 0.08]:
    zz = np.array([-3.0, 4.5, 0.9 + 1j * y, 0.9 - 1j * y])
    u = -2 * np.real(np.sum(1 / (xs[:, None] - zz[None, :]), 1))
    prof[str(y)] = u.tolist()
    # front steepness at the pair: max |u_x| near x = 0.9 ~ 1/y^2
near = (xs > 0.3) & (xs < 1.5)
out["burgers_front_slope"] = {y: float(np.max(np.abs(np.gradient(np.array(prof[y]), xs))[near])) for y in prof}
np.savez("burgers.npz", xs=xs, **{f"u_{k}": np.array(v) for k, v in prof.items()})

# ---- 2. Calogero-Moser: along the heat flow, z'' = -8 sum 1/(z-w)^3 and E = 1/2 sum z'^2 - 4 sum_{j<k} 1/(z_j-z_k)^2 = 0
z = np.sort(rng.normal(0, 3, 9)).astype(complex)
h = 2e-3
def step(z, dt):
    s = solve_ivp(flow(), [0, dt], pack(z), rtol=1e-13, atol=1e-15, method="DOP853"); return unpack(s.y[:, -1])
zp = step(z, h)
s = solve_ivp(flow(), [0, -h], pack(z), rtol=1e-13, atol=1e-15, method="DOP853"); zm = unpack(s.y[:, -1])
D = z[:, None] - z[None, :]; np.fill_diagonal(D, 1.0)
I3 = 1 / D**3; np.fill_diagonal(I3, 0)
acc_cm = -8 * I3.sum(1)
I1 = 1 / D; np.fill_diagonal(I1, 0); Fz = I1.sum(1)
acc_num = -4 * ((Fz[:, None] - Fz[None, :]) * I1**2).sum(1)   # chain rule on z' = 2F
v = 2 * F(z); iu = np.triu_indices(len(z), 1)
E = 0.5 * np.sum(v**2) - 4 * np.sum(1 / (z[:, None] - z[None, :])[iu]**2)
out["calogero_moser"] = {"max_rel_err_accel": float(np.max(np.abs(acc_num - acc_cm)) / np.max(np.abs(acc_cm))), "energy": float(abs(E)),
                         "identity_sumF2": float(abs(np.sum(F(z)**2) - 2 * np.sum((1 / (z[:, None] - z[None, :]))[iu]**2)))}
# with an off-line pair the identity and zero energy still hold (holomorphic), but the force is not a physical gradient
z2 = np.concatenate([np.sort(rng.normal(0, 3, 6)), [0.5 + 0.7j, 0.5 - 0.7j]]).astype(complex)
v2 = 2 * F(z2); iu2 = np.triu_indices(len(z2), 1)
out["calogero_moser"]["energy_offline"] = float(abs(0.5 * np.sum(v2**2) - 4 * np.sum((1 / (z2[:, None] - z2[None, :]))[iu2]**2)))

# ---- 3. 2D parallel dipoles pointing along the wire: energy -cos(2 phi)/r^2 = -Re 1/(z-w)^2. On the line = CM potential.
phis = np.linspace(0, 2 * math.pi, 13); r = 1.7
dip = np.array([-(1 - 2 * math.cos(p)**2) * -1 / r**2 for p in phis])     # 2D: p1.p2 - 2(p1.rhat)(p2.rhat) over r^2, p along x
cm = np.array([-np.real(1 / (r * np.exp(1j * p))**2) for p in phis])
out["dipole_match"] = float(np.max(np.abs(dip - cm)))

# ---- 4. any real potential's gradient is ANTI-holomorphic: force = -conj(g'(z)) for V = Re g. Check numerically.
g = lambda z: -4 / z**2; gp = lambda z: 8 / z**3
pts = np.array([1.2 + 0.7j, -0.4 + 1.1j, 2.0 - 0.3j])
eps = 1e-6
Vx = (np.real(g(pts + eps)) - np.real(g(pts - eps))) / (2 * eps)
Vy = (np.real(g(pts + 1j * eps)) - np.real(g(pts - 1j * eps))) / (2 * eps)
force = -(Vx + 1j * Vy)
out["gradient_is_antiholomorphic"] = float(np.max(np.abs(force - (-np.conj(gp(pts))))))
# holomorphic flow f = 2/z is the gradient of Phi = 2 log r in the metric dx^2 - dy^2: xdot = Phi_x, ydot = -Phi_y
Phi = lambda z: 2 * np.log(np.abs(z))
Px = (Phi(pts + eps) - Phi(pts - eps)) / (2 * eps); Py = (Phi(pts + 1j * eps) - Phi(pts - 1j * eps)) / (2 * eps)
out["lorentzian_gradient"] = float(np.max(np.abs((Px - 1j * Py) - 2 / pts)))

# ---- 5. collapse on rewind: 1/r beads (ours, line charges, SC vortices) vs 1/r^2 beads (point charges on a wire)
def rewind_pair(power, d0=0.5, N=30):
    x = np.arange(N, dtype=float); x[N // 2] = x[N // 2 - 1] + d0
    x = np.sort(x)
    def rhs(t, x):
        D = x[:, None] - x[None, :]; np.fill_diagonal(D, np.inf)
        return -2 * (np.sign(D) / np.abs(D)**power).sum(1)        # rewind: attraction
    def ev(t, x): return np.min(np.diff(np.sort(x))) - 1e-5
    ev.terminal = True
    s = solve_ivp(rhs, [0, 5], x, rtol=1e-12, atol=1e-14, events=ev, dense_output=True, method="DOP853")
    T = s.t_events[0][0]
    left = np.logspace(-8, -2.5, 30)
    gaps = np.array([np.min(np.diff(np.sort(s.sol(T - l)))) for l in left])
    return float(T), float(np.polyfit(np.log(left), np.log(gaps), 1)[0]), left, gaps
T1, a1, l1, g1 = rewind_pair(1); T2, a2, l2, g2 = rewind_pair(2)
out["collapse"] = {"one_over_r": {"T": T1, "exponent": a1}, "one_over_r2": {"T": T2, "exponent": a2}}
np.savez("collapse.npz", l1=l1, g1=g1, l2=l2, g2=g2)

# ---- 6. threshold as a Hill radius: pair push 4/d vs crowd's tidal spring 2 S d  -> d* = sqrt(2/S)
snaps = np.load("../known-tools/rewind_snaps.npz")
S = snaps["S"]; d0 = snaps["d0"]; mid = snaps["mid"]
dstar = np.sqrt(2 / S)
out["hill"] = {"median_dstar_over_gap": float(np.median(dstar / d0)), "share_inside": float(np.mean(d0 < dstar)),
               "median_dstar_over_meangap_at_gamma_2700": None}
json.dump(out, open("phys.json", "w"), indent=1)
print(json.dumps(out, indent=1))
