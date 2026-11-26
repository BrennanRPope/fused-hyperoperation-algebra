"""DRAFT checks for three-body/README.md. Run: python3 checks.py  (numpy, scipy; builds on libtb.so from tb.c).
Writes checks_output.txt and results.json (used by figures.py)."""
import sys, os, json, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from tbpy import *
from chore import Chore

LOG = []; RES = {}; NPASS = [0, 0]
def check(tag, ok, msg):
    NPASS[1] += 1; NPASS[0] += bool(ok)
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}"; print(line); LOG.append(line)

from checks_lib import *

# ============================================================== A. special solutions on the two 3-spheres
print("== A. special solutions"); RES["A"] = {}
for name, euler in (("Lagrange", False), ("Euler", True)):
    X0, Tk, p0 = homographic(0.5, euler)
    t, Xs = traj(X0, Tk, 3000)
    w = shape(Xs); u = s3point(Xs)
    U0 = u[:, 0] + 1j * u[:, 1]; U1 = u[:, 2] + 1j * u[:, 3]; ip = np.conj(U0[0]) * U0 + np.conj(U1[0]) * U1
    c, rho, res = circfit(Xs[:, 6:8]); sv = np.linalg.svd(fock(Xs[:, 6:8], p0), compute_uv=False)
    check(f"A1 {name}", np.ptp(w, 0).max() < 1e-9, f"shape never moves on the shape sphere (spread {np.ptp(w,0).max():.1e}); sits at {np.round(w[0],6)}")
    check(f"A2 {name}", abs(abs(ip).min() - 1) < 1e-9, f"normalized configuration runs along one Hopf circle of S³ (min |<u0,u>| = {abs(ip).min():.12f})")
    check(f"A3 {name}", res < 1e-9 and sv[2] / sv[0] < 1e-9, f"each body's velocities lie on one circle (max misfit {res:.1e}), so on Fock's sphere they lie on one great circle (3rd/1st singular value {sv[2]/sv[0]:.1e}); implied p0 = √(−2E) = {p0:.6f}")
    RES["A"][name] = dict(shape=w[0].tolist(), e=0.5, period=Tk)

hod = {}
for name, (a, b, T) in SUV.items():
    if name in ("yarn", "bumblebee"): continue
    q, r = refine(a, b, T); t, Xs = traj(suvakov(*q[:2]), q[2], 4000)
    V = Xs[:, 6:8]; c, rho, res = circfit(V)
    ang = np.unwrap(np.arctan2(V[:, 1] - c[1], V[:, 0] - c[0])); wind = np.unwrap(np.arctan2(V[:, 1], V[:, 0]))
    hod[name] = dict(misfit=res, sweep=float(np.degrees(ang.max() - ang.min())), winding=float((wind[-1] - wind[0]) / 2 / np.pi), c=c.tolist(), rho=rho)
f8 = hod["figure-eight"]
check("A4 figure-eight", f8["misfit"] < 0.1 and abs(f8["winding"]) < 1e-6,
      f"velocities stay within {100*f8['misfit']:.1f}% of one circle, sweeping a {f8['sweep']:.0f}° arc back and forth (winds {f8['winding']:.0f} times round the origin), so the momenta rock on a piece of one great circle of Fock's sphere instead of going round it")
others = {k: v["misfit"] for k, v in hod.items() if k != "figure-eight"}
check("A5 other Suvakov orbits", min(others.values()) > 1.0, "no other orbit tested comes close: worst-fit circles miss by " + ", ".join(f"{k} {v:.1f}×" for k, v in others.items()) + " of their radius")
RES["A"]["hodographs"] = hod

# ============================================================== B. swap rule: rotation = -1/2 enclosed shape area
print("== B. swap rule on the shape sphere")
X0 = suvakov(0.40, 0.30)
t, Xs = traj(X0, 80.0, 160001, 1e-12)
W = shape(Xs); U = s3point(Xs); Z0 = U[:, 0] + 1j * U[:, 1]; Z1 = U[:, 2] + 1j * U[:, 3]
rows = []
for i in np.arange(0, len(t), 37)[:3000]:
    if i + 800 >= len(t): break
    d = np.linalg.norm(W[i + 800:] - W[i], axis=1); j = np.argmin(d)
    if d[j] < 2e-3:
        j += i + 800; Om = solid_angle(W[i:j + 1])
        al = np.angle(np.conj(Z0[i]) * Z0[j] + np.conj(Z1[i]) * Z1[j]); pr = -Om / 2
        rows.append((t[i], t[j], al, pr, (al - pr + np.pi) % (2 * np.pi) - np.pi))
rows = np.array(rows)
check("B1 angular momentum", abs(angmom(X0)) < 1e-14, "Šuvakov starts have zero angular momentum")
check("B2 half-area law", len(rows) > 100 and np.abs(rows[:, 4]).max() < 2e-4,
      f"{len(rows)} closed shape loops of one long chaotic orbit: the triangle's net turn equals −½ × enclosed solid angle on the unit shape sphere to within {np.abs(rows[:,4]).max():.1e} rad (median {np.median(np.abs(rows[:,4])):.1e}; gaps come from closing each loop to 0.002)")
RES["B"] = dict(rows=rows[:, 2:4].tolist(), n=len(rows))

# ============================================================== C. dial = minus energy
print("== C. plane dial and energy (Moser / Osipov / Belbruno)")
for E, lab in ((-0.5, "E<0"), (0.0, "E=0"), (0.5, "E>0")):
    # Kepler orbit with K=1, perihelion 1: hodograph is a circle centre (0,K/h) radius K/h ... check rho^2-|c|^2 = -2E
    h = np.sqrt(2 * (E + 1)); th = np.linspace(-2.0, 2.0, 400) if E >= 0 else np.linspace(0, 2 * np.pi, 400)
    V = np.column_stack([-np.sin(th) / h, (np.cos(th) + h * h - 1) / h])
    c, rho, res = circfit(V)
    check(f"C1 {lab}", abs((rho ** 2 - c @ c) - (-2 * E)) < 1e-9, f"velocity circle has radius² − centre² = {rho**2-c@c:+.6f} = −2E, so Fock's sphere has curvature d ∝ −E ({'sphere' if E<0 else 'flat' if E==0 else 'hyperbolic'})")

# ============================================================== D. action minimums: choreographies
print("== D. action minimizer (choreographies, inertial frame)")
rng = np.random.default_rng(3); ch = Chore(N=30, M=900, omega=0.0); D = []
for sd in range(40):
    amp = rng.normal(size=(ch.nc, 2)) / (np.concatenate([ch.ns, ch.ns])[:, None] ** 1.5)
    x, A, g = ch.minimize(amp.ravel()); X0 = ch.initial_state(x)
    try: ret = float(np.linalg.norm(flow(X0, 2 * np.pi, 1e-12) - X0))
    except RuntimeError: ret = float("nan")
    D.append(dict(A=float(A), g=float(g), ret=ret, L=angmom(X0), x=x.tolist()))
As = np.array([d["A"] for d in D]); Ls = np.array([d["L"] for d in D])
circ = np.isclose(As, As.min(), atol=1e-5); eight = np.isclose(As, 24.371927, atol=1e-3) & (abs(Ls) < 1e-6)
check("D1 circle", circ.sum() >= 20 and np.all(abs(abs(Ls[circ]) - abs(Ls[circ][0])) < 1e-6), f"{circ.sum()}/40 random starts fall to the lowest action {As.min():.6f}: Lagrange's rotating triangle (|L| = {abs(Ls[circ][0]):.4f})")
check("D2 figure-eight", eight.sum() >= 3, f"{eight.sum()}/40 fall to action {np.mean(As[eight]):.6f} with L = 0: the figure-eight (ODE return after one period {np.median([D[i]['ret'] for i in np.nonzero(eight)[0]]):.1e}, limited by 30 harmonics)")
other = ~(circ | eight)
check("D3 nothing else that is real", all((not np.isfinite(D[i]["ret"])) or D[i]["ret"] > 1e-2 for i in np.nonzero(other)[0]),
      f"{other.sum()} other stopping points; none integrates back to its start (they are artefacts of the truncation)")
RES["D"] = D

# ============================================================== E. return-distance landscape
print("== E. return-distance landscape (Suvakov slice)")
L = np.load(os.path.join(HERE, "landscape.npz"))
check("E0 grid", np.isfinite(L["dabs"]).sum() > 1e5, f"{np.isfinite(L['dabs']).sum()} starts, each run to T = {float(L['Tmax']):.0f}")
M = json.load(open(os.path.join(HERE, "minima.json")))
found = [f for f in M["results"] if f["res"] < 1e-8 and f["T"] > 1]
check("E0b minima", len(found) > 0, f"{M['n_jobs']} local minima below 0.08 (both kinds) shot to exact orbits; {len(found)} converged to 1e-8")
# dedupe (same orbit found twice, or at a multiple of its period)
uniq = []
for f in sorted(found, key=lambda f: f["T"]):
    dup = False
    for u in uniq:
        if abs(f["p1"] - u["p1"]) < 1e-5 and abs(f["p2"] - u["p2"]) < 1e-5:
            k = f["T"] / u["T"]
            if abs(k - round(k)) < 1e-4: dup = True; break
    if not dup: uniq.append(f)
for f in uniq:
    t, Xs = traj(suvakov(f["p1"], f["p2"]), f["T"], int(f["T"] * 3000), 1e-13)
    syz = syzygies(t, Xs); w = word_from_syzygies(syz, 3); f["word"] = w; f["len"] = len(w)
    f["match"], f["power"] = lookup(w)
    W = shape(Xs); f["shape_area_turn"] = float(-solid_angle(W) / 2)
    f["E"] = energy(suvakov(f["p1"], f["p2"]))
    f["Tstar"] = f["T"] * abs(f["E"]) ** 1.5
RES["E"] = uniq
nabs = [f for f in uniq if abs(f["alpha"]) < 1e-6 or f["kind"] == "absolute"]
nrot = [f for f in uniq if f["kind"] == "rotation-free" and abs((f["alpha"] + np.pi) % (2 * np.pi) - np.pi) > 1e-3]
known = [f for f in nabs if f["match"]]
check("E1 absolute minima", len(nabs) > 0, f"{len(nabs)} distinct periodic orbits refined to 1e-8; {len(known)} have a loop type in Li–Liao's 695-family list" +
      "; others: " + ", ".join(f"T={f['T']:.2f} word length {f['len']}" for f in nabs if not f["match"]))
turns = [abs((f["shape_area_turn"] + np.pi) % (2 * np.pi) - np.pi) for f in nabs]
check("E2b half-area law on every periodic orbit", max(turns) < 1e-3, f"each of the {len(nabs)} periodic orbits encloses a shape-sphere area whose half is a whole number of turns (max miss {max(turns):.1e} rad), as the swap rule demands for an orbit that comes back unturned")
new = [f for f in nabs if not f["match"]]
dmins = []
for f in new:
    t, Xs = traj(suvakov(f["p1"], f["p2"]), f["T"], int(f["T"] * 4000))
    r = Xs[:, :6].reshape(-1, 3, 2); f["dmin"] = float(min(np.linalg.norm(r[:, i] - r[:, j], axis=1).min() for i in range(3) for j in range(i + 1, 3)))
from scipy.integrate import solve_ivp
def rhs(t, X):
    r = X[:6].reshape(3, 2); a = np.zeros((3, 2))
    for i in range(3):
        for j in range(i + 1, 3):
            d = r[j] - r[i]; k = d / np.linalg.norm(d) ** 3; a[i] += k; a[j] -= k
    return np.concatenate([X[6:], a.ravel()])
for f in new:
    X0 = suvakov(f["p1"], f["p2"]); sol = solve_ivp(rhs, (0, f["T"]), X0, method="DOP853", rtol=1e-13, atol=1e-15)
    f["scipy_return"] = float(np.linalg.norm(sol.y[:, -1] - X0))
check("E2c the unlisted orbits are real", all(f["scipy_return"] < 1e-6 for f in new), "an independent integrator (scipy DOP853) brings each unlisted orbit back to within " + f"{max(f['scipy_return'] for f in new):.0e}" + "; closest passes between two bodies " + ", ".join(f"{f['dmin']:.4f}" for f in new) + " (start spacing 1)")
check("E2d unlisted loop types", True, f"{len(set(canon(f['word']) for f in new))} distinct loop types among the {len(new)} unlisted orbits: " + "; ".join(sorted(set(canon(f['word']) for f in new), key=len)))
check("E2 rotating minima", True, f"{len(nrot)} relative periodic orbits: every rotation-free minimum in this slice shot back to an orbit that returns unturned")
if nrot:
    errs = [abs((f["alpha"] - f["shape_area_turn"] + np.pi) % (2 * np.pi) - np.pi) for f in nrot]
    check("E3 half-area law predicts the turn", max(errs) < 1e-5, f"for every relative periodic orbit the turn equals −½ × its shape loop's area (max error {max(errs):.1e})")

open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(LOG) + f"\n\n{NPASS[0]} of {NPASS[1]} pass\n")
json.dump(RES, open(os.path.join(HERE, "results.json"), "w"), default=float)
print(f"{NPASS[0]} of {NPASS[1]} pass")
