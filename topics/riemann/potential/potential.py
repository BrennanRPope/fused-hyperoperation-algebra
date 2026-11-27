"""Is there a potential energy for the zeros under heat flow whose minimum or balance point sits at t = 0?
Brennan's hypothesis (Nov 27 2026 call): "I think it's action minimization." / "maybe it's potential energy."

Setting (Rodgers-Tao normalisation, z = 2*gamma): real zeros of H_t move by
    dz_j/dt = 2 sum_{k != j} 1/(z_j - z_k)                       (de Bruijn; Csordas-Smith-Varga 1994; Rodgers-Tao 2018)
The log-gas (Coulomb gas) potential energy is W = - sum_{j<k} log|z_j - z_k|.
On the real line dW/dt = -2 sum_j F_j^2 <= 0, F_j = sum_k 1/(z_j - z_k): the heat flow is gradient descent of W.

Parts
  1. The snag, measured on our 5500 zeros: excess energy over the evenly spaced 'crystal' at t = 0, .01, .05, .2.
  2. Off the line the flow is NOT downhill: a conjugate pair being captured climbs W, which spikes to +inf at capture.
  3. Lambda is the time of the last spike of W (last collision), not a minimum.
  4. Balance of 'trend vs wiggle' (Li / Bombieri-Lagarias numbers) sits at Lambda wherever Lambda is: toy demo.
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp

out = {}
snaps = np.load("../known-tools/rewind_snaps.npz")
R = json.load(open("../known-tools/rewind.json"))

# ---------- 1. excess energy over the crystal ----------------------------------------------------------
def excess(z, lo=1000, hi=4500, K=10, win=101):
    """mean over j in [lo,hi), m = 1..K of -log(s_jm / <s_m>_local); s_jm = z_{j+m} - z_j.
    By Jensen (log is concave), the window average is >= 0, with equality only for even spacing."""
    tot = 0.0; cnt = 0
    ker = np.ones(win) / win
    for m in range(1, K + 1):
        s = z[m:] - z[:-m]
        loc = np.convolve(s, ker, mode="same")
        r = s[lo:hi] / loc[lo:hi]
        tot += -np.log(r).sum(); cnt += len(r)
    return tot / cnt

ts = [0.0, 0.01, 0.05, 0.2]
E = {t: excess(snaps[str(t)]) for t in ts}
print("excess energy per zero-neighbour:", E)
# references at t = 0: GUE (sine process) and Poisson, same statistic
rng = np.random.default_rng(1)
def gue_unfolded(n=3000):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); H = (A + A.conj().T) / 2
    ev = np.linalg.eigvalsh(H) / math.sqrt(n)                         # semicircle on [-2, 2]... radius sqrt(2)
    ev = ev / math.sqrt(2) * 2                                        # radius 2
    Nsc = lambda x: n * (0.5 + (x * np.sqrt(4 - x**2) / 4 + 2 * np.arcsin(x / 2)) / (2 * math.pi))
    return Nsc(np.clip(ev, -2, 2))
g = gue_unfolded()
Egue = excess(g, lo=600, hi=2400)
pois = np.cumsum(rng.exponential(size=5500))
Epois = excess(pois)
lat = np.arange(5500.0)
Elat = excess(lat)
print("GUE", Egue, "Poisson", Epois, "crystal", Elat)
out["excess"] = {"ours": {str(t): E[t] for t in ts}, "GUE": Egue, "Poisson": Epois, "crystal": Elat}

# ---------- toy dynamics in the complex plane -------------------------------------------------------------
def rhs(t, v, n):
    z = v[:n] + 1j * v[n:]
    D = z[:, None] - z[None, :]
    np.fill_diagonal(D, np.inf)
    f = 2 * (1 / D).sum(1)
    return np.concatenate([f.real, f.imag])
def W(z):
    D = np.abs(z[:, None] - z[None, :]); iu = np.triu_indices(len(z), 1)
    return -np.log(D[iu]).sum()
def Fvec(z):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf); return (1 / D).sum(1)

# 2a. real configuration: W only falls, dW/dt = -2 sum F^2
zr = np.concatenate([np.arange(-15, 16) + rng.normal(0, 0.25, 31)]).astype(complex)
zr = np.sort(zr.real).astype(complex)
n = len(zr)
sol = solve_ivp(rhs, [0, 0.3], np.concatenate([zr.real, zr.imag]), args=(n,), rtol=1e-10, atol=1e-12,
                dense_output=True)
tt = np.linspace(0, 0.3, 61)
Wr = np.array([W(sol.sol(t)[:n] + 1j * sol.sol(t)[n:]) for t in tt])
dW_num = np.gradient(Wr, tt)
dW_th = np.array([-2 * np.sum(np.abs(Fvec(sol.sol(t)[:n] + 1j * sol.sol(t)[n:]))**2) for t in tt])
out["real_toy"] = {"monotone": bool(np.all(np.diff(Wr) < 0)),
                   "max_rel_err_dWdt": float(np.max(np.abs(dW_num[2:-2] - dW_th[2:-2]) / np.abs(dW_th[2:-2])))}
print("real toy", out["real_toy"])

# 2b. one off-line conjugate pair (and its mirror -z, since H_t is even): capture climbs W
def offline_config(y0=0.6, x0=5.0, L=20):
    base = np.array([k + 0.5 for k in range(L)], float)
    base = base[np.abs(base - x0) > 0.9]                              # leave room for the pair
    pts = list(base) + list(-base)
    pts = np.array(pts, complex)
    pair = np.array([x0 + 1j * y0, x0 - 1j * y0, -x0 + 1j * y0, -x0 - 1j * y0])
    return np.concatenate([pts, pair])
z0 = offline_config()
n = len(z0)
def capture_event(t, v, n): return v[2 * n - 4] - 1e-3          # imag part of pair's upper member
capture_event.terminal = True
sol2 = solve_ivp(rhs, [0, 1.0], np.concatenate([z0.real, z0.imag]), args=(n,), rtol=1e-10, atol=1e-12,
                 events=capture_event, dense_output=True)
tcap = float(sol2.t_events[0][0])
tt2 = np.linspace(0, tcap, 200)
traj2 = np.array([sol2.sol(t)[:n] + 1j * sol2.sol(t)[n:] for t in tt2])
W2 = np.array([W(z) for z in traj2])
y2 = traj2[:, n - 4].imag
out["offline_toy"] = {"t_capture": tcap, "W_start": float(W2[0]), "W_end": float(W2[-1]),
                      "W_rises": bool(np.all(np.diff(W2) > 0))}
print("off-line toy", out["offline_toy"])
# after capture: pair is real, separate; continue from a just-split real pair to show W falls again
zs = traj2[-1].copy()
upper, lower = n - 4, n - 3
xm = zs[upper].real
zs[upper] = xm + 2e-3; zs[lower] = xm - 2e-3
zs[n - 2] = -xm + 2e-3; zs[n - 1] = -xm - 2e-3
zs = zs.real.astype(complex)
sol3 = solve_ivp(rhs, [0, 0.4], np.concatenate([zs.real, zs.imag]), args=(n,), rtol=1e-10, atol=1e-12,
                 dense_output=True)
tt3 = np.linspace(0, 0.4, 200)
W3 = np.array([W(sol3.sol(t)[:n] + 1j * sol3.sol(t)[n:]) for t in tt3])
out["offline_toy"]["W_falls_after"] = bool(np.all(np.diff(W3) < 0))
np.savez("toy_traj.npz", tt2=tt2, W2=W2, y2=y2, tt3=tt3 + tcap, W3=W3, tt=tt, Wr=Wr)


# 2c. same test on a ring (period P), so the crowd cannot spread out and the even crystal is stationary.
#     kernel: sum over images of 1/(z + kP) = (pi/P) cot(pi z / P);  W_ring = -sum log|sin(pi (z_j - z_k)/P)|
P = 40.0
def rhs_ring(t, v, n):
    z = v[:n] + 1j * v[n:]
    D = z[:, None] - z[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        C = (math.pi / P) / np.tan(math.pi * D / P)
    np.fill_diagonal(C, 0)
    f = 2 * C.sum(1)
    return np.concatenate([f.real, f.imag])
def W_ring(z):
    D = np.abs(np.sin(math.pi * (z[:, None] - z[None, :]) / P)); iu = np.triu_indices(len(z), 1)
    return -np.log(D[iu]).sum()
base = np.arange(P) + 0.5
base = base[np.abs(base - 20.0) > 0.9]
zring = np.concatenate([base.astype(complex), [20.0 + 0.6j, 20.0 - 0.6j]])
nr = len(zring)
def cap_r(t, v, n): return v[2 * n - 2] - 1e-3
cap_r.terminal = True
solr = solve_ivp(rhs_ring, [0, 1.0], np.concatenate([zring.real, zring.imag]), args=(nr,), rtol=1e-10,
                 atol=1e-12, events=cap_r, dense_output=True)
tcr = float(solr.t_events[0][0])
ttr = np.linspace(0, tcr, 300)
trr = np.array([solr.sol(t)[:nr] + 1j * solr.sol(t)[nr:] for t in ttr])
Wrr = np.array([W_ring(z) for z in trr])
Wpair = -np.log(np.abs(np.sin(math.pi * (trr[:, -2] - trr[:, -1]) / P)))
yr = trr[:, -2].imag
zs = trr[-1].copy(); xm = zs[-2].real
zs[-2] = xm + 2e-3; zs[-1] = xm - 2e-3; zs = zs.real.astype(complex)
solr2 = solve_ivp(rhs_ring, [0, 0.5], np.concatenate([zs.real, zs.imag]), args=(nr,), rtol=1e-10, atol=1e-12,
                  dense_output=True)
ttr2 = np.linspace(0, 0.5, 300)
Wrr2 = np.array([W_ring(solr2.sol(t)[:nr] + 1j * solr2.sol(t)[nr:]) for t in ttr2])
crystal = np.arange(P) + 0.5
out["ring_toy"] = {"t_capture": tcr, "capture_formula": float(np.log1p(2 * (math.pi**2 / 3 - 2 * 0) * 0) if False else 0),
                   "W_rises_to_capture": bool(np.all(np.diff(Wrr) > 0)), "W_falls_after": bool(np.all(np.diff(Wrr2) < 0)),
                   "W_start": float(Wrr[0]), "W_at_capture": float(Wrr[-1]), "W_after_0.5": float(Wrr2[-1]),
                   "W_crystal": float(W_ring(crystal.astype(complex))),
                   "dW_dt_start_numeric": float((Wrr[1] - Wrr[0]) / (ttr[1] - ttr[0]))}
del out["ring_toy"]["capture_formula"]
print("ring toy", out["ring_toy"])
np.savez("ring_traj.npz", ttr=ttr, Wrr=Wrr, Wpair=Wpair, yr=yr, ttr2=ttr2 + tcr, Wrr2=Wrr2,
         z0=zring, zcap=trr[-1], zend=solr2.sol(0.5)[:nr])

# ---------- 3. Lambda = time of the last spike: pair energy near a collision ----------------------------
# pair formula (frozen background S): d^2(t) = 2/S + (d0^2 - 2/S) exp(-4 S t). Pair energy -1/2 log|d^2|.
tc = snaps["tc"]; mid = snaps["mid"]; S = snaps["S"]; d0 = snaps["d0"]
j = int(np.nanargmax(tc))
u = lambda t, S_, d_: 2 / S_ + (d_**2 - 2 / S_) * np.exp(-4 * S_ * t)
out["lehmer_pair"] = {"gamma_mid": float(mid[j] / 2), "t_c": float(tc[j]), "S": float(S[j]), "d0": float(d0[j]),
                      "u_at_tc": float(u(tc[j], S[j], d0[j]))}
# the latest collision per height band moves toward 0 from below as height grows
ok = np.isfinite(tc)
bands = np.linspace(mid[0], mid[-1], 12)
late = []
for a, b in zip(bands[:-1], bands[1:]):
    m = ok & (mid >= a) & (mid < b)
    if m.any(): late.append([float((a + b) / 4), float(np.max(tc[m]))])
out["latest_collision_by_height"] = late
print("latest collisions by height", late)
# off-line world: capture of a beta = 1 zero above 3e12, from squeeze.json
Q = json.load(open("../spin-orbit/squeeze.json"))
out["offline_world_tcap_3e12"] = Q["tcap_at_3e12_beta1"]

# ---------- 4. Li / Bombieri-Lagarias numbers on the toy: the trend-vs-wiggle balance flips at capture ----
def li_min(z, nmax=3000):
    rho = 0.5 + 1j * z / 2                                            # z real  <=>  rho on the critical line
    w = 1 - 1 / rho
    ns = np.arange(1, nmax + 1)
    lam = np.array([np.sum(1 - w**k).real for k in ns])
    return float(lam.min()), int(ns[np.argmin(lam)])
tli = np.linspace(0, tcap * 0.999, 25)
lis = [li_min(sol2.sol(t)[:n] + 1j * sol2.sol(t)[n:]) for t in tli]
tli_after = np.linspace(0.002, 0.4, 10)
lis_after = [li_min(sol3.sol(t)[:n] + 1j * sol3.sol(t)[n:]) for t in tli_after]
out["li_toy"] = {"t": tli.tolist(), "min_lambda": [a for a, _ in lis], "argmin_n": [b for _, b in lis],
                 "t_after": (tli_after + tcap).tolist(), "min_lambda_after": [a for a, _ in lis_after]}
print("Li min before capture:", [round(a, 3) for a, _ in lis][:5], "... after:", [round(a, 3) for a, _ in lis_after][:3])
json.dump(out, open("potential.json", "w"), indent=1)
