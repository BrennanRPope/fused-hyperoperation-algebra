"""Checks for the Klein-bottle flow solver and the facts used in the write-up. Run: python3 checks.py"""
import json, time
import numpy as np
from flow import Flow, Grid, project, project_mean, random_vorticity, normalize_energy

results = []


def check(name, ok, detail):
    results.append(dict(name=name, ok=bool(ok), detail=detail))
    print(("PASS " if ok else "FAIL ") + name + " :: " + detail, flush=True)


t_start = time.time()
N = 128
g = Grid(N)
X, Y = g.X, g.Y

# C1 symmetry is kept by the equations themselves (no re-projection), to round-off
for s in ["Trect", "K"]:
    w0 = normalize_energy(N, project(g, random_vorticity(N, 6, 3), s))
    f = Flow(N, 5e-4, s, w0, reproject=False)
    f.advance_to(10.0)
    w = f.w()
    err = np.max(np.abs(w - (g.shift(w) if s == "Trect" else -g.glide(w)))) / np.max(np.abs(w))
    check(f"C1 {s}: symmetry kept by Navier-Stokes without re-projection (t=10)", err < 1e-10, f"max rel. error {err:.1e}")

# C2 lone vortex psi = cos x + sin y is allowed on K, is an exact steady Euler flow, decays as exp(-nu t)
psi = np.cos(X) + np.sin(Y)
allowed = np.max(np.abs(psi + g.glide(psi)))
w = psi  # -Lap psi = psi for this shell
f = Flow(N, 1e-2, "K", w, reproject=False)
nl = np.max(np.abs(g.ifft(f.rhs(f.wh))))
f.advance_to(20.0)
amp = np.max(np.abs(f.psi())) / 2.0
check("C2 lone vortex cos x + sin y: allowed on K, steady (nonlinear term 0), decays exactly like e^(-nu t)",
      allowed < 1e-12 and nl < 1e-10 and abs(amp - np.exp(-0.2)) < 1e-6,
      f"glide-odd error {allowed:.1e}; nonlinear term {nl:.1e}; amplitude at t=20 {amp:.8f} vs e^-0.2 = {np.exp(-0.2):.8f}")

# C3 energy budget dE/dt = -2 nu Z on K (decaying turbulence)
w0 = normalize_energy(N, project(g, random_vorticity(N, 6, 4), "K"))
f = Flow(N, 1e-3, "K", w0)
Es, Zs, ts = [], [], []
for T in np.linspace(0, 4, 81):
    f.advance_to(T); Es.append(f.energy()); Zs.append(f.enstrophy()); ts.append(f.t)
Es, Zs, ts = map(np.array, (Es, Zs, ts))
dE = np.gradient(Es, ts)
rel = np.max(np.abs(dE[2:-2] + 2 * 1e-3 * Zs[2:-2])) / np.max(2 * 1e-3 * Zs)
check("C3 K: energy budget dE/dt = -2 nu * enstrophy", rel < 1e-3, f"max rel. mismatch {rel:.1e}")

# C4 harmonic (uniform, steady, never-decaying) flows: torus 2, Klein bottle 1
dims = {}
for s in ["T2", "Trect", "K"]:
    basis = [project_mean(1.0, 0.0, s), project_mean(0.0, 1.0, s)]
    dims[s] = int(np.linalg.matrix_rank(np.array(basis)))
check("C4 harmonic flows (first Betti number): T2 2, Trect 2, K 1", dims == {"T2": 2, "Trect": 2, "K": 1}, str(dims))

# C5 spectrum of the (twisted) Laplacian: multiplicities of the first eigenvalues on each surface
def multiplicities(s, nmax=4):
    out = {}
    for lam in range(1, nmax + 1):
        modes = []
        for kx in range(-3, 4):
            for ky in range(-3, 4):
                if kx * kx + ky * ky != lam:
                    continue
                for ph in (0.0, np.pi / 2):
                    m = np.cos(kx * X + ky * Y + ph)
                    modes.append(project(g, m, s).ravel())
        out[lam] = int(np.linalg.matrix_rank(np.array(modes), tol=1e-8)) if modes else 0
    return out
mult = {s: multiplicities(s) for s in ["T2", "Trect", "K"]}
ok = mult["T2"][1] == 4 and mult["Trect"][1] == 2 and mult["K"][1] == 3
check("C5 lowest streamfunction shell (|k|=1): T2 4 modes, Trect 2, K 3 (cos x, sin x, sin y)", ok, json.dumps(mult))

# C6 odd moments of vorticity vanish identically on K, not on Trect
wK = project(g, random_vorticity(N, 6, 7), "K")
wR = project(g, random_vorticity(N, 6, 7), "Trect")
m1K, m3K = np.mean(wK), np.mean(wK**3) / np.mean(wK**2) ** 1.5
m3R = np.mean(wR**3) / np.mean(wR**2) ** 1.5
check("C6 K: mean and skewness of vorticity are 0 for every allowed field; Trect skewness free",
      abs(m1K) < 1e-12 and abs(m3K) < 1e-12 and abs(m3R) > 1e-3,
      f"K mean {m1K:.1e}, K skew {m3K:.1e}, Trect skew {m3R:.3f}")

# C7 'total spin' depends on where you cut the Klein bottle open; on the torus it is always 0
def cut_integral(w, a_index):
    n = N // 2
    idx = (np.arange(n) + a_index) % N
    return np.sum(w[idx, :]) * (2 * np.pi / N) ** 2
cK = np.array([cut_integral(wK, a) for a in range(N)])
cR = np.array([cut_integral(wR, a) for a in range(N)])
antip = np.max(np.abs(cK[: N // 2] + cK[N // 2:]))
check("C7 integral of vorticity over a fundamental rectangle: Trect 0 for every cut; K depends on the cut and flips sign after half a turn",
      np.max(np.abs(cR)) < 1e-10 and np.ptp(cK) > 1e-2 and antip < 1e-10,
      f"Trect max |.| {np.max(np.abs(cR)):.1e}; K range {cK.min():.3f}..{cK.max():.3f}; K(a)+K(a+pi) max {antip:.1e}")

# C8 Taylor-Green (allowed on T2) decays exactly: solver sanity
psi = np.sin(X) * np.sin(Y)
f = Flow(N, 1e-2, "T2", 2 * psi)
f.advance_to(10.0)
amp = np.max(np.abs(f.psi()))
check("C8 Taylor-Green vortex on T2 decays exactly like e^(-2 nu t)", abs(amp - np.exp(-0.2)) < 1e-6, f"{amp:.8f} vs {np.exp(-0.2):.8f}")

# C9 resolution: N=128 vs N=256 agree on energy and enstrophy for a K run to t=10
vals = {}
# same initial field on both grids: build on 256 from a band-limited field defined analytically
rng = np.random.default_rng(11)
modes = [(kx, ky, rng.normal(), rng.uniform(0, 2 * np.pi)) for kx in range(-9, 10) for ky in range(-9, 10) if 0 < kx * kx + ky * ky <= 81]
for n in (128, 256):
    gg = Grid(n)
    w0 = sum(a * np.cos(kx * gg.X + ky * gg.Y + ph) for kx, ky, a, ph in modes)
    w0 = normalize_energy(n, project(gg, w0, "K"))
    f = Flow(n, 1e-3, "K", w0)
    f.advance_to(10.0)
    vals[n] = (f.energy(), f.enstrophy())
eE = abs(vals[128][0] - vals[256][0]) / vals[256][0]
eZ = abs(vals[128][1] - vals[256][1]) / vals[256][1]
check("C9 resolution: N=128 and N=256 agree (K, nu=1e-3, t=10)", eE < 1e-3 and eZ < 1e-2, f"energy rel. diff {eE:.1e}, enstrophy rel. diff {eZ:.1e}")

# C10 a uniform 'upward' flow is impossible on K: its best K-version is zero; upward in the rectangle = shear layer at the flipped edge
v = np.ones_like(X)
# v is a velocity component: on K it must satisfy v(x+pi,-y) = -v(x,y); projection of v=1:
vK = 0.5 * (v - g.glide(v))
check("C10 uniform flow along the flipped edge (v = 1) is not allowed on K (projection is 0)", np.max(np.abs(vK)) < 1e-14, f"max |projection| {np.max(np.abs(vK)):.1e}")

npass = sum(r["ok"] for r in results)
print(f"\n{npass} of {len(results)} pass ({time.time()-t_start:.0f} s)")
json.dump(dict(results=results, passed=npass, total=len(results)), open("checks_results.json", "w"), indent=1)
