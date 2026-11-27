import json, math
d = json.load(open("dialwick.json"))
res = []
def chk(name, ok, val):
    res.append((name, ok)); print(f"{'PASS' if ok else 'FAIL'}  {name}: {val}")

q = d["quarter_turn"]
chk("W1 quarter turn exact: G_t(w) = H_t(iw)", q["G_equals_H_of_iw_relerr"] < 1e-12, q["G_equals_H_of_iw_relerr"])
chk("W2 turned frame runs FORWARD heat dG/dt = +G_ww", q["forward_heat_relerr"] < 1e-5, q["forward_heat_relerr"])
p = d["physical_every_t"]
chk("W3 phi = G_t(x) > 0 on t in [-1,1], |x| <= 60 (physical Burgers, RH or not)", p["min_G_over_grid"] > 0, p["min_G_over_grid"])
b = d["burgers_residual"]
chk("W4 u = -2G'/G solves u_t + u u_x = u_xx (viscosity +1)", b["max_abs"] < 1e-5 * b["scale"] * 10 or b["max_abs"] < 1e-9, b)
pc = d["poles_on_centerline"]
chk("W5 poles of u on the centreline at 2*gamma_1..5", pc["max_err"] < 1e-3, pc["max_err"])
chk("W6 magnetisation m(h) = G'/G increasing in h (variance > 0) for t = -0.3, 0, 0.3", d["magnetisation"]["increasing_in_h_all_t"], d["magnetisation"]["m_at_h40"])
t = d["toy"]
chk("W7 toy c + e^t cos z: off the line before ln c, real after", t["off_line_before"] and t["real_after"], t["Lambda"])
chk("W8 toy in the turned frame: pole columns merge on the centreline at t = ln c; Lee-Yang iff t >= ln c",
    t["columns_merge_at_x0"] and t["lee_yang_at"] == {"0.0": False, "0.3": False, str(t["Lambda"]): True, "0.6": True}, t["lee_yang_at"])
chk("W9 turned-frame pole law dw/dt = -2 sum 1/(w_j - w_k) matches exact heat flow", d["rotated_pole_ode"]["max_err"] < 1e-9, d["rotated_pole_ode"]["max_err"])
chk("W10 side pole columns merge in finite time (forward heat, turned frame)", d["rotated_merger"]["merge_time"] is not None, d["rotated_merger"]["merge_time"])
chk("W11 metric half turn: v = (1-d)/2 Hol + (1+d)/2 conj(Hol) (d = metric value)", d["kappa_split_identity_err"] < 1e-12, d["kappa_split_identity_err"])
chk("W12 ripple growth = d * 2(pi q - q^2/2) on a crystal (exact ring sums)", d["kappa_dispersion"]["max_relerr"] < 1e-9, d["kappa_dispersion"]["max_relerr"])
kp = d["kappa_pair"]
ok = kp["-1.0"]["y_end"] < 2e-3 and kp["-0.5"]["y_end"] < 2e-3 and abs(kp["0.0"]["y_end"] - 0.3) < 1e-9 and kp["0.5"]["y_end"] > 0.3 and kp["1.0"]["y_end"] > kp["0.5"]["y_end"]
r = kp["-0.5"]["t_end"] / kp["-1.0"]["t_end"]
chk("W13 off-line pair: captured for d < 0 (time ~ 1/|d|), frozen at d = 0, escapes for d > 0", ok and abs(r - 2) < 0.05, {k: kp[k] for k in ["-1.0", "-0.5", "0.0", "1.0"]})
vp = d["viscosity_phase"]
okv = abs(vp["0"][0] + 1) < 1e-12 and abs(vp["0.7853981633974483"][1] - 1) < 1e-12 and abs(vp["1.5707963267948966"][0] - 1) < 1e-12
chk("W14 viscosity(theta) = -e^{-2i theta}: -1, i, +1 at 0, pi/4, pi/2; PDE check at pi/4", okv and d["viscosity_phase_check_relerr"] < 1e-4, d["viscosity_phase_check_relerr"])
# W15 the quarter turn conjugates the plane law's d to -d: law_d(ix, ia) = i * law_{-d}(x, a)
import random
law = lambda x, a, dd: (x + a) / (1 - dd * a * x)
random.seed(1); errs = []
for _ in range(50):
    x, a, dd = random.uniform(-2, 2), random.uniform(-2, 2), random.choice([1.0, 0.5, -1.0, 2.0])
    errs.append(abs(law(1j * x, 1j * a, dd) - 1j * law(x, a, -dd)))
chk("W15 x -> ix turns the plane law at d into the law at -d (tan addition <-> tanh addition)", max(errs) < 1e-12, max(errs))
n = sum(o for _, o in res); print(f"\n{n}/{len(res)} checks pass")
