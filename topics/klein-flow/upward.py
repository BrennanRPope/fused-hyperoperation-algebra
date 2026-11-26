"""Experiment B: 'everything flows upward in the rectangle'.

Same local picture in the fundamental rectangle [0,pi] x [0,2pi]: fluid moving up (+y) at speed 1.
  Trect (straight gluing): that is the uniform flow V = 1, a harmonic flow. It never slows down.
  K (flipped gluing): going through the flipped edge turns 'up' into 'down', so upward flow in the rectangle
     is glued to downward flow: a shear layer sits on the flipped edge. v(x) = tanh(sin(x)/delta) on the cover.
Both get the same tiny random nudge (projected onto each surface). Usage: python3 upward.py N nu T
"""
import sys
import numpy as np
from flow import Flow, Grid, project

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
nu = float(sys.argv[2]) if len(sys.argv) > 2 else 5e-4
T = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
delta = 0.15
g = Grid(N)
rng = np.random.default_rng(2)
kk = np.sqrt(g.K2)
noise = g.ifft(np.where((kk > 0) & (kk < 12), 1.0, 0.0) * np.exp(2j * np.pi * rng.random(g.K2.shape)))
noise *= 1e-2 / np.std(noise)

# K: v = tanh(sin x / delta), u = 0  ->  vorticity w = dv/dx
vK = np.tanh(np.sin(g.X) / delta)
wK = np.cos(g.X) / delta / np.cosh(np.sin(g.X) / delta) ** 2
runs = {
    "K": Flow(N, nu, "K", wK + noise, U=0.0, V=0.0),
    "Trect": Flow(N, nu, "Trect", noise, U=0.0, V=1.0),
}
times = np.linspace(0, T, 81)
snap_t = [0, 6, 12, 20, T]
out = {}
for name, f in runs.items():
    E, Emean = [], []
    for t in times:
        f.advance_to(t)
        E.append(f.energy()); Emean.append(0.5 * (f.U**2 + f.V**2))
        if any(abs(t - s) < 1e-9 for s in snap_t):
            u, v = f.velocity()
            out[f"{name}_w_{t:g}"] = f.w().astype(np.float32)
            out[f"{name}_v_{t:g}"] = (v + f.V).astype(np.float32)
            out[f"{name}_u_{t:g}"] = (u + f.U).astype(np.float32)
    out[f"{name}_E"] = np.array(E); out[f"{name}_Emean"] = np.array(Emean)
    print(name, "E start", E[0], "E end", E[-1], flush=True)
out["t"] = times
np.savez_compressed("upward.npz", **out)
