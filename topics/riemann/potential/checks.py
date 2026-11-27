"""Checks for the potential-energy page. Run after potential.py."""
import json, math
import numpy as np
O = json.load(open("potential.json"))
R = np.load("ring_traj.npz")
res = []
def ck(name, ok, detail=""):
    res.append((name, bool(ok))); print(("PASS" if ok else "FAIL"), name, detail)

# P1 Jensen: window excess >= 0, = 0 only for the crystal (periodic, exact version)
rng = np.random.default_rng(7)
def ex_periodic(z, P, K=5):
    n = len(z); tot = 0
    for m in range(1, K + 1):
        s = (np.roll(z, -m) - z) % P
        tot += -np.log(s / (m * P / n)).sum()
    return tot / (n * K)
vals = [ex_periodic(np.sort(rng.uniform(0, 100, 100)), 100.0) for _ in range(200)]
ck("P1 Jensen: excess log energy >= 0 for 200 random rings, 0 for the crystal",
   min(vals) > 0 and abs(ex_periodic(np.arange(100.0), 100.0)) < 1e-12, f"min {min(vals):.3f}")

# P2 heat flow is gradient descent of W on the line
rt = O["real_toy"]
ck("P2 real configuration: W falls monotonically, dW/dt = -2 sum F^2", rt["monotone"] and rt["max_rel_err_dWdt"] < 1e-2,
   f"rel err {rt['max_rel_err_dWdt']:.1e}")

# P3 general formula dW/dt = -2 Re sum F_j^2 (holds off the line too), checked at an off-line configuration
z = np.array([-3, -2, -1, 1.3, 2, 3, 0.2 + 0.5j, 0.2 - 0.5j], complex)
def F(z):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf); return (1 / D).sum(1)
def W(z):
    D = np.abs(z[:, None] - z[None, :]); iu = np.triu_indices(len(z), 1); return -np.log(D[iu]).sum()
h = 1e-6
num = (W(z + h * 2 * F(z)) - W(z - h * 2 * F(z))) / (2 * h)
th = -2 * np.real(np.sum(F(z)**2))
ck("P3 dW/dt = -2 Re sum F_j^2 holds off the line too (so its sign can flip); at the ring start it is POSITIVE",
   abs(num - th) < 1e-6 * abs(th) and O["ring_toy"]["dW_dt_start_numeric"] > 0,
   f"numeric {num:.6f}, formula {th:.6f}; ring start dW/dt = {O['ring_toy']['dW_dt_start_numeric']:.2f}")

# P4 ring toy: off-line start below the crystal, W rises to capture, falls after toward the crystal
r = O["ring_toy"]
ck("P4 ring: start below crystal, rises to capture, falls after, ends at crystal",
   r["W_start"] < r["W_crystal"] and r["W_rises_to_capture"] and r["W_falls_after"]
   and abs(r["W_after_0.5"] - r["W_crystal"]) < 1e-2,
   f"start {r['W_start']:.3f}, crystal {r['W_crystal']:.3f}, end {r['W_after_0.5']:.3f}")

# P5 the log energy has no floor off the line: lifting the pair lowers W without limit
base = np.arange(-10, 11, dtype=float); base = base[base != 0]
Ws = [W(np.concatenate([base, [1j * y, -1j * y]]).astype(complex)) for y in (0.1, 0.5, 1, 2, 5, 20, 100)]
ck("P5 lifting an off-line pair higher keeps lowering W (no minimum off the line)", np.all(np.diff(Ws) < 0),
   f"W: {Ws[0]:.2f} -> {Ws[-1]:.2f}")

# P6 ring capture time vs the pair formula t = ln(1 + 2 S y^2)/(4S), S from the other zeros (ring sum)
z0 = R["z0"]; P = 40.0; x0 = z0[-2].real; y0 = z0[-2].imag
others = z0[:-2].real
S = np.sum((math.pi / P)**2 / np.sin(math.pi * (x0 - others) / P)**2)
tf = math.log1p(2 * S * y0**2) / (4 * S)
ck("P6 ring capture time within 10% of the frozen-crowd pair formula", abs(tf - r["t_capture"]) / r["t_capture"] < 0.10,
   f"run {r['t_capture']:.4f}, formula {tf:.4f}")

# P7 our zeros: excess energy falls at every snapshot; t = 0 sits near GUE, far from Poisson
e = O["excess"]; ours = [e["ours"][k] for k in ("0.0", "0.01", "0.05", "0.2")]
ck("P7 our zeros: excess energy falls 0 -> .01 -> .05 -> .2", np.all(np.diff(ours) < 0), str([round(v, 4) for v in ours]))
ck("P8 at t = 0 our zeros are GUE-like (between crystal and GUE, 10x below Poisson)",
   0 < ours[0] < e["GUE"] and ours[0] < e["Poisson"] / 10, f"ours {ours[0]:.4f}, GUE {e['GUE']:.4f}, Poisson {e['Poisson']:.3f}")

# P9 the Lehmer pair: pair formula u(t_c) = 0 (the spike)
L = O["lehmer_pair"]
ck("P9 closest pair: gap^2 hits 0 at t_c (energy spike), t_c < 0", abs(L["u_at_tc"]) < 1e-8 and L["t_c"] < 0,
   f"t_c {L['t_c']:.3e}")
# P10 every latest collision is before 0
ck("P10 every height band's latest collision is at t < 0 (no spike after 0 among our zeros)",
   all(t < 0 for _, t in O["latest_collision_by_height"]))
# P11 Li numbers: negative before capture, all positive after
li = O["li_toy"]
ck("P11 toy Li numbers: a negative one exists before capture, none after",
   all(v < 0 for v in li["min_lambda"]) and all(v > 0 for v in li["min_lambda_after"]))
# P12 capture is uphill in the toy too (finite, mirrored): pair term grows
T = np.load("toy_traj.npz")
ck("P12 toy: the off-line height y shrinks to 0 at capture (de Bruijn's height² falls)", np.all(np.diff(T["y2"]) < 0) and T["y2"][-1] < 2e-3)
n_ok = sum(o for _, o in res); print(f"\n{n_ok}/{len(res)} pass")
