"""Three maps onto the zeros' heat flow (Rodgers-Tao normalisation z = 2*gamma, dz_j/dt = 2 sum 1/(z_j - z_k)).
A. Blow-up: does a backward collision look like a forward finite-time blow-up (Keller-Segel, u_t = u_xx + u^p,
   bubble collapse)?  Collapse law, threshold ("critical mass") S d0^2 = 2.
B. Pressure: add an outside force -a z (harmonic pressure) or a general V'(z). Equilibria, equivalence to zoom +
   clock change, what pressure would hold zeta's zeros still at t = 0.
C. Twist: heat time t -> tau e^{i theta}. theta = pi/2 is the Schrodinger-type flow. Do zeros stay on the line?
   What is conserved? Does anything map back?
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
out = {}
rng = np.random.default_rng(3)
snaps = np.load("../known-tools/rewind_snaps.npz")

def ring_rhs(P, extra=None, rot=1.0):
    def f(t, v):
        n = len(v) // 2; z = v[:n] + 1j * v[n:]
        D = z[:, None] - z[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            C = (math.pi / P) / np.tan(math.pi * D / P)
        np.fill_diagonal(C, 0)
        w = rot * 2 * C.sum(1)
        if extra is not None: w = w + extra(z)
        return np.concatenate([w.real, w.imag])
    return f
def line_rhs(extra=None, rot=1.0):
    def f(t, v):
        n = len(v) // 2; z = v[:n] + 1j * v[n:]
        D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf)
        w = rot * 2 * (1 / D).sum(1)
        if extra is not None: w = w + extra(z)
        return np.concatenate([w.real, w.imag])
    return f
pack = lambda z: np.concatenate([np.real(z), np.imag(z)])
unpack = lambda v: v[:len(v) // 2] + 1j * v[len(v) // 2:]

# ===================== A. blow-up ====================================================================
# A1. collapse law of a close pair rewound on a ring of 40 (backward = run with rot = -1)
P = 40.0
base = np.arange(P) + 0.5
z = base.copy(); z[20] = z[19] + 0.35                    # a close pair (gap 0.35) among unit spacing
z = np.sort(z)
def gap_ev(t, v): zz = unpack(v).real; return np.min(np.diff(np.sort(zz))) - 1e-4
gap_ev.terminal = True
sol = solve_ivp(ring_rhs(P, rot=-1.0), [0, 1], pack(z.astype(complex)), rtol=1e-12, atol=1e-14, events=gap_ev,
                dense_output=True, method="DOP853")
Tc = float(sol.t_events[0][0])
tt = Tc - np.logspace(-7, -1.5, 40)
gaps = np.array([np.min(np.diff(np.sort(unpack(sol.sol(t)).real))) for t in tt])
alpha = np.polyfit(np.log(Tc - tt), np.log(gaps), 1)[0]
out["A1_collapse"] = {"rewind_time_to_collision": Tc, "exponent_gap_vs_time_left": float(alpha),
                      "prefactor_gap2_over_timeleft": float(np.median(gaps**2 / (Tc - tt)))}
print("A1", out["A1_collapse"])
np.savez("collapse.npz", left=Tc - tt, gaps=gaps)

# A2. threshold ("critical mass"): frozen-crowd rule S d0^2 < 2 collapses on rewind; test with the full ring
res = []
for d0 in np.linspace(0.3, 1.4, 23):
    z = base.copy(); z[19] = 19.5 - (d0 - 1) / 2 * 0 ; z[19] = 20.0 - d0 / 2; z[20] = 20.0 + d0 / 2
    others = np.delete(z, [19, 20]); m = 20.0
    S = np.sum((math.pi / P)**2 / np.sin(math.pi * (m - others) / P)**2)
    s2 = solve_ivp(ring_rhs(P, rot=-1.0), [0, 6.0], pack(z.astype(complex)), rtol=1e-10, atol=1e-12,
                   events=gap_ev, method="DOP853")
    col = len(s2.t_events[0]) > 0
    which = None
    if col:
        zz = np.sort(unpack(s2.y_events[0][0]).real); which = int(np.argmin(np.diff(zz)))
    x_ = S * d0**2 / 2
    tform = float(-math.log(1 / (1 - x_)) / (4 * S)) if x_ < 1 else None
    res.append({"d0": float(d0), "S_d0sq": float(S * d0**2), "collapses_within_6": col,
                "t": float(s2.t_events[0][0]) if col else None, "own_pair": which == 19,
                "t_formula_frozen": (-tform if tform is not None else None)})
out["A2_threshold"] = res
print("A2", [(round(r["S_d0sq"], 2), r["own_pair"], round(r["t"], 3), r["t_formula_frozen"] and round(r["t_formula_frozen"], 3)) for r in res])
# share of our zeros below the threshold
x = snaps["S"] * snaps["d0"]**2
out["A2_our_share_below"] = float(np.mean(x < 2))

# ===================== B. pressure =====================================================================
# B1. harmonic pressure -a z : equilibrium = scaled Hermite zeros (Stieltjes 1885)
n = 12; a = 1.0
z0 = np.sort(rng.normal(0, 2, n)).astype(complex)
sB = solve_ivp(line_rhs(extra=lambda z: -a * z), [0, 40], pack(z0), rtol=1e-12, atol=1e-14, method="DOP853")
zeq = np.sort(unpack(sB.y[:, -1]).real)
from numpy.polynomial.hermite import hermroots
herm = np.sort(hermroots([0] * n + [1]))                  # physicists' H_n
# stationarity: sum 1/(x_j - x_k) = (a/2) x_j  <=> y = x sqrt(a/2) are zeros of H_n (sum 1/(y_j-y_k) = y_j)
out["B1_hermite"] = {"max_err": float(np.max(np.abs(zeq * math.sqrt(a / 2) - herm)))}
print("B1", out["B1_hermite"])
# B2. pressure = zoom + clock: z(s) = e^{-a s} w(t(s)), t(s) = (e^{2 a s} - 1)/(2a), w the plain flow
z0 = np.sort(rng.normal(0, 2, n)).astype(complex)
s_end = 0.6
sp = solve_ivp(line_rhs(extra=lambda z: -a * z), [0, s_end], pack(z0), rtol=1e-12, atol=1e-14, method="DOP853")
t_end = (math.exp(2 * a * s_end) - 1) / (2 * a)
wp = solve_ivp(line_rhs(), [0, t_end], pack(z0), rtol=1e-12, atol=1e-14, method="DOP853")
out["B2_zoom_clock"] = {"s_end": s_end, "t_end": t_end,
                        "max_err": float(np.max(np.abs(unpack(sp.y[:, -1]) - math.exp(-a * s_end) * unpack(wp.y[:, -1]))))}
print("B2", out["B2_zoom_clock"])
# B3. "stabilised" backward flow: with pressure a, running s -> -inf only reaches t = -1/(2a).
#     Ring pair from A1 collides at rewind Tc. Pressure a > 1/(2 Tc) keeps the collision out of reach forever.
out["B3_clock"] = {"Tc": Tc, "a_needed_to_never_reach_collision": 1 / (2 * Tc)}
# B4. the pressure that would hold zeta's zeros still at t = 0 is minus their velocity field
z0s, z1s = snaps["0.0"], snaps["0.01"]
v = (z1s - z0s) / 0.01
lo, hi = 1000, 4500
ker = np.ones(101) / 101
vs = np.convolve(v, ker, mode="same")                   # smooth drift (the crowd spreading as density grows)
res_v = (v - vs)[lo:hi]
gap = np.diff(z0s)[lo:hi]
out["B4_velocity"] = {"rms_residual_velocity": float(np.sqrt(np.mean(res_v**2))),
                      "rms_over_smooth_drift": float(np.sqrt(np.mean(res_v**2)) / np.mean(np.abs(vs[lo:hi]))),
                      "mean_gap": float(gap.mean()),
                      "frac_sign_changes": float(np.mean(np.sign(res_v[1:]) != np.sign(res_v[:-1])))}
# same at later times (finite differences between snapshots)
for (ta, tb) in [("0.01", "0.05"), ("0.05", "0.2")]:
    vv = (snaps[tb] - snaps[ta]) / (float(tb) - float(ta))
    rr = (vv - np.convolve(vv, ker, mode="same"))[lo:hi]
    out["B4_velocity"][f"rms_residual_{ta}_{tb}"] = float(np.sqrt(np.mean(rr**2)))
print("B4", out["B4_velocity"])
np.savez("velocity.npz", z=z0s[lo:hi], res=res_v, vs=vs[lo:hi])

# ===================== C. twist ========================================================================
# 300 consecutive zeros put on a ring (keeps their local statistics), flowed in twisted time tau e^{i theta}
w = snaps["0.0"][2000:2300]
Pz = (w[-1] - w[0]) * 300 / 299
zr = (w - w[0]).astype(complex)
twist = {}
for th in [0.0, math.pi / 8, math.pi / 4, math.pi / 2]:
    rot = complex(math.cos(th), math.sin(th))
    s = solve_ivp(ring_rhs(Pz, rot=rot), [0, 0.05], pack(zr), rtol=1e-9, atol=1e-11, t_eval=np.linspace(0, 0.05, 26),
                  method="DOP853")
    Z = np.array([unpack(s.y[:, k]) for k in range(s.y.shape[1])])
    twist[f"{th:.4f}"] = {"tau": s.t.tolist(), "max_im": np.max(np.abs(Z.imag), 1).tolist(),
                          "rms_im": np.sqrt(np.mean(Z.imag**2, 1)).tolist(),
                          "sum_z2_change_re": (np.sum(Z**2, 1) - np.sum(Z[0]**2)).real.tolist(),
                          "sum_z2_change_im": (np.sum(Z**2, 1) - np.sum(Z[0]**2)).imag.tolist()}
    if abs(th - math.pi / 2) < 1e-9:
        np.savez("twist_traj.npz", tau=s.t, Z=Z)
out["C_twist"] = twist
# conserved / linear quantities on the line (polynomial case): sum z (conserved), sum z^2 = sum z0^2 + 2N(N-1) t (any complex t)
N = 10; z0 = np.sort(rng.normal(0, 3, N)).astype(complex); tc_ = 0.03 * complex(0, 1)
sC = solve_ivp(line_rhs(rot=1j), [0, 0.03], pack(z0), rtol=1e-12, atol=1e-14, method="DOP853")
zz = unpack(sC.y[:, -1])
out["C_moments"] = {"sum_z_change": float(abs(zz.sum() - z0.sum())),
                    "sum_z2_err": float(abs((zz**2).sum() - (z0**2).sum() - 2 * N * (N - 1) * tc_)),
                    "max_im_after": float(np.max(np.abs(zz.imag)))}
print("C", out["C_moments"], {k: v["max_im"][-1] for k, v in twist.items()})
# initial speed off the line at theta = pi/2: Im z_j' = 2 F_j (the real heat-flow velocity, turned 90 degrees)
json.dump(out, open("pbt.json", "w"), indent=1)
