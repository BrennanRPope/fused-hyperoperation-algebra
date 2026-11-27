"""Spin-orbit resonance (Schnittman 2004) applied to RH, two ways.

I. Primes as spins. Each Euler factor rotates: p^{-it} has phase -t log p. At the zeros gamma:
   - each prime's phase leans toward 180 deg: mean cos(gamma log p) ~ -(log p / sqrt p) / log(T / 2pi)  (Landau 1912,
     Gonek 1993; the lower-order term of Ford-Zaharescu 2005 in the u.d. of alpha*gamma, Hlawka 1975)
   - two primes never lock to each other: the phase difference gamma log(p/q) has no lean (no Landau spike at p/q)
II. Zeros as spins under heat flow. An off-line pair x +- iy is pulled onto the line by forward heat time:
   y' = -1/y - 2 S y  =>  capture at t = ln(1 + 2 S y0^2) / (4S) <= y0^2/2 (de Bruijn 1950's bound).
   Real pairs split into two families under rewind by S d0^2 < 2 (collide, become a conjugate pair) or not.
"""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = np.array(json.load(open("../round9/zeros.json"))["gamma"])
N = len(g); T = g[-1]
out = {}
# ---------- I. primes as spins ----------
P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
lean = {}
for p in P:
    c = np.cos(g * math.log(p))
    pred = -np.mean((math.log(p) / math.sqrt(p)) / np.log(g / (2 * math.pi)))
    lean[p] = {"mean_cos": float(c.mean()), "predicted": float(pred), "se": float(c.std() / math.sqrt(N))}
out["lean"] = lean
pairs = {}
for (p, q) in [(3, 2), (5, 3), (7, 5), (5, 2), (7, 3), (11, 7)]:
    c = np.cos(g * math.log(p / q))
    pairs[f"{p}/{q}"] = {"mean_cos": float(c.mean()), "se": float(c.std() / math.sqrt(N))}
out["pairs"] = pairs
# control: the same statistics at random heights with the zeros' density
rng = np.random.default_rng(3)
tr = np.sort(rng.uniform(g[0], T, N))
out["control_random_heights_mean_cos2"] = float(np.cos(tr * math.log(2)).mean())
# lean of p = 2 in height windows
edges = np.linspace(g[0], T, 12)
win = []
for a, b in zip(edges[:-1], edges[1:]):
    m = (g >= a) & (g < b); c = np.cos(g[m] * math.log(2))
    win.append({"mid": float((a + b) / 2), "mean": float(c.mean()), "se": float(c.std() / math.sqrt(m.sum())),
                "pred": float(-(math.log(2) / math.sqrt(2)) / math.log((a + b) / 2 / (2 * math.pi)))})
out["lean2_windows"] = win

# ---------- II. zeros as spins under heat flow ----------
z = 2 * g
allz = np.concatenate([z, -z])
def S_at(x, skip):
    d = x - allz
    d = d[np.abs(d) > 1e-9]
    s = (1 / d**2).sum()
    return s
# insert an artificial off-line pair at x0 +- i y0 between real zeros, run forward heat flow on the pair with the
# real zeros moving too (local window), compare capture time with the formula and de Bruijn's bound
def capture(j, y0, half=120):
    x0 = (z[j] + z[j + 1]) / 2
    idx = np.arange(j - half, j + half + 2)
    w = np.concatenate([z[idx].astype(complex), [x0 + 1j * y0, x0 - 1j * y0]])
    frozen = np.concatenate([np.delete(z, idx), -z])
    S = (1 / (x0 - np.concatenate([z[idx], frozen]))**2).sum()
    t = 0.0; ys = [(0.0, y0)]
    while True:
        D = w[:, None] - w[None, :]; np.fill_diagonal(D, np.inf)
        v = 2 * (1 / D).sum(1) + 2 * (1 / (w[:, None] - frozen[None, :])).sum(1)
        y = w[-2].imag
        dt = min(2e-4, 0.01 * y**2)
        w = w + dt * v; t += dt
        ys.append((t, w[-2].imag))
        if w[-2].imag < 2e-3 or t > 5:
            break
    tf = math.log(1 + 2 * S * y0**2) / (4 * S)
    return {"j": int(j), "x0": float(x0), "y0": y0, "S": float(S), "t_capture_nbody": t, "t_formula": tf,
            "debruijn_bound": y0**2 / 2}, np.array(ys)
caps = []; traj = {}
for j, y0 in [(3000, 0.3), (3000, 0.6), (3000, 1.0), (1200, 0.6)]:
    r, ys = capture(j, y0); caps.append(r); traj[f"{j}_{y0}"] = ys
    print(r, flush=True)
out["capture"] = caps
json.dump(out, open("spinorbit.json", "w"), indent=1)
np.savez("spinorbit_traj.npz", **traj)
print(json.dumps({k: out[k] for k in ("pairs", "control_random_heights_mean_cos2")}, indent=1))
print({p: (round(v["mean_cos"], 4), round(v["predicted"], 4), round(v["se"], 4)) for p, v in lean.items()})
