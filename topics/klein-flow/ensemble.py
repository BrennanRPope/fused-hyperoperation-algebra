"""Experiment A: matched decaying turbulence on the square torus T2, the rectangle torus Trect and the Klein bottle K.

Same random seed on the cover, projected onto each surface's allowed fields, then rescaled to the same energy
(0.5 per unit area, no mean flow). Same viscosity, same grid. Records energy, enstrophy, spectra, vorticity
skewness, the lowest-shell mode amplitudes, and snapshots.

Usage: python3 ensemble.py N nu T nseeds [outdir]
"""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool
from flow import Flow, Grid, project, random_vorticity, normalize_energy

SURFACES = os.environ.get("SURFACES", "T2,Trect,K").split(",")


def lowest_shell(flow):
    """Amplitudes of the |k|=1 modes of the streamfunction (cover Fourier coefficients, normalised so that a
    unit-amplitude mode cos(x) has amplitude 1)."""
    g = flow.g
    N = g.N
    ph = flow.wh * g.K2inv / (N * N)
    # rfft layout: [kx index, ky index]; kx=+1 -> index 1, ky=+1 -> index 1, kx=0 ky=1 -> [0,1]
    ax = 2 * np.abs(ph[1, 0])  # cos x / sin x pair amplitude
    ay = 2 * np.abs(ph[0, 1])  # cos y / sin y pair amplitude
    # split the y pair into its cos y and sin y parts (on K only sin y is allowed)
    cy = 2 * ph[0, 1].real
    sy = -2 * ph[0, 1].imag
    return dict(ax=float(ax), ay=float(ay), cosy=float(cy), siny=float(sy))


def skewness(w):
    return float(np.mean(w**3) / np.mean(w**2) ** 1.5)


def run(args):
    surface, seed, N, nu, T, outdir = args
    w0 = random_vorticity(N, k0=float(os.environ.get("K0", "6")), seed=seed)
    w0 = normalize_energy(N, project(Grid(N), w0, surface), 0.5)
    f = Flow(N, nu, surface, w0)
    times = np.unique(np.concatenate([np.linspace(0, 20, 41), np.linspace(20, T, 61), [60.0, 150.0]]))
    snap_times = [0, 5, 20, 60, 150, T]
    rec = dict(t=[], E=[], Z=[], skew=[], shell_frac=[], ax=[], ay=[], cosy=[], siny=[])
    spectra = {}
    snaps = {}
    t0 = time.time()
    for t in times:
        f.advance_to(t)
        w = f.w()
        E, Z = f.energy(), f.enstrophy()
        sp = f.spectrum()
        ls = lowest_shell(f)
        rec["t"].append(float(f.t)); rec["E"].append(float(E)); rec["Z"].append(float(Z))
        rec["skew"].append(skewness(w)); rec["shell_frac"].append(float(sp[1] / sp.sum()))
        for kk in ("ax", "ay", "cosy", "siny"):
            rec[kk].append(ls[kk])
        if any(abs(t - s) < 1e-9 for s in [0, 5, 20, 60, 150, T]):
            spectra[str(round(t, 3))] = sp[: N // 2].tolist()
        if seed < 2 and any(abs(t - s) < 1e-9 for s in snap_times):
            snaps[f"w_{round(t,3)}"] = w.astype(np.float32)
            snaps[f"psi_{round(t,3)}"] = f.psi().astype(np.float32)
    rec["spectra"] = spectra
    rec["wall"] = time.time() - t0
    name = f"{surface}_s{seed}"
    with open(os.path.join(outdir, name + ".json"), "w") as fh:
        json.dump(rec, fh)
    if snaps:
        np.savez_compressed(os.path.join(outdir, name + "_snaps.npz"), **snaps)
    return name, rec["wall"], rec["shell_frac"][-1]


if __name__ == "__main__":
    N = int(sys.argv[1]); nu = float(sys.argv[2]); T = float(sys.argv[3]); ns = int(sys.argv[4])
    outdir = sys.argv[5] if len(sys.argv) > 5 else f"runs_N{N}"
    os.makedirs(outdir, exist_ok=True)
    first = int(os.environ.get("FIRST_SEED", "0"))
    jobs = [(s, seed, N, nu, T, outdir) for seed in range(first, first + ns) for s in SURFACES]
    with Pool(4) as p:
        for name, wall, frac in p.imap_unordered(run, jobs):
            print(name, f"{wall:.0f}s", f"lowest-shell fraction {frac:.3f}", flush=True)
