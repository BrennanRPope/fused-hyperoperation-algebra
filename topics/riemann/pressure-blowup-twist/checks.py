"""Checks for pressure-blowup-twist. Run after pbt.py."""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
O = json.load(open("pbt.json")); res = []
def ck(n, ok, d=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + n + ("  " + d if d else ""))
A = O["A1_collapse"]
ck("K1 backward collision: gap ∝ (time left)^½, gap² = 8 × time left",
   abs(A["exponent_gap_vs_time_left"] - 0.5) < 0.01 and abs(A["prefactor_gap2_over_timeleft"] / 8 - 1) < 0.01,
   f"exponent {A['exponent_gap_vs_time_left']:.4f}, gap²/time {A['prefactor_gap2_over_timeleft']:.3f}")
r = O["A2_threshold"]
below = [q for q in r if q["S_d0sq"] < 1.6]
ck("K2 under the threshold the chosen pair collapses, at the frozen-crowd time within 10%",
   all(q["own_pair"] and abs(q["t"] / q["t_formula_frozen"] - 1) < 0.10 for q in below), f"{len(below)} runs")
above = [q for q in r if q["S_d0sq"] > 2]
ck("K3 over the threshold the chosen pair survives, but some other pair still collapses (every run)",
   all(not q["own_pair"] for q in above) and all(q["collapses_within_6"] for q in r), f"{len(above)} runs above, {len(r)} total")
ck("K4 about half of our zeros' neighbour pairs sit under the threshold", 0.4 < O["A2_our_share_below"] < 0.6,
   f"{O['A2_our_share_below']:.3f}")
ck("K5 inward pressure −a·z: the settled zeros are the scaled Hermite zeros (Stieltjes)", O["B1_hermite"]["max_err"] < 1e-10,
   f"{O['B1_hermite']['max_err']:.1e}")
ck("K6 pressure = zoom e^(−as) + clock t = (e^(2as) − 1)/(2a), exact", O["B2_zoom_clock"]["max_err"] < 1e-10,
   f"{O['B2_zoom_clock']['max_err']:.1e}")
a = 1.0; s = -30.0
ck("K7 with pressure a, running backward forever only reaches t = −1/(2a)", abs((math.exp(2 * a * s) - 1) / (2 * a) + 1 / (2 * a)) < 1e-12)
B = O["B4_velocity"]
ck("K8 the force that would hold ζ's zeros still at t = 0 flips sign often and exceeds the smooth drift",
   B["frac_sign_changes"] > 0.4 and B["rms_over_smooth_drift"] > 1, f"{B['frac_sign_changes']:.2f}, {B['rms_over_smooth_drift']:.2f}")
ck("K9 that force shrinks as heat time goes on (relaxation toward the crystal)",
   B["rms_residual_velocity"] > B["rms_residual_0.01_0.05"] > B["rms_residual_0.05_0.2"],
   f"{B['rms_residual_velocity']:.2f} > {B['rms_residual_0.01_0.05']:.2f} > {B['rms_residual_0.05_0.2']:.2f}")
T = O["C_twist"]
ck("K10 untwisted heat flow keeps every zero on the line", max(T["0.0000"]["max_im"]) < 1e-9)
sp = {k: T[k]["max_im"][1] / T[k]["tau"][1] for k in ("0.3927", "0.7854", "1.5708")}
rat = [sp[k] / math.sin(float(k)) for k in sp]
ck("K11 any twist θ > 0 moves zeros off the line at once, starting speed ∝ sin θ (within 15%)",
   all(v > 0 for v in sp.values()) and (max(rat) / min(rat) < 1.15), str([round(v, 2) for v in rat]))
M = O["C_moments"]
ck("K12 under complex heat time: Σz conserved, Σz² = Σz₀² + 2N(N−1)t exactly", M["sum_z_change"] < 1e-10 and M["sum_z2_err"] < 1e-10)
# K13 untwisting: run θ = π/2 forward then back, recover the start
def f(rot):
    def g(t, v):
        n = len(v) // 2; z = v[:n] + 1j * v[n:]; D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf)
        w = rot * 2 * (1 / D).sum(1); return np.concatenate([w.real, w.imag])
    return g
z0 = np.sort(np.random.default_rng(5).normal(0, 3, 10)).astype(complex); v0 = np.concatenate([z0.real, z0.imag])
v1 = solve_ivp(f(1j), [0, 0.05], v0, rtol=1e-12, atol=1e-14, method="DOP853").y[:, -1]
v2 = solve_ivp(f(-1j), [0, 0.05], v1, rtol=1e-12, atol=1e-14, method="DOP853").y[:, -1]
ck("K13 removing the twist returns the zeros exactly to the line (it is the same family, run back)",
   np.max(np.abs(v2 - v0)) < 1e-9 and np.max(np.abs(v1[10:])) > 1e-2, f"{np.max(np.abs(v2 - v0)):.1e}")
print(f"\n{sum(res)} of {len(res)} checks pass")
