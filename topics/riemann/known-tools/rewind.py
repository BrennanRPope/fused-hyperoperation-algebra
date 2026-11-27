"""Newman's Lambda applied to our 5500 zeros (z = 2*gamma, Rodgers-Tao normalisation).

Zero dynamics of H_t (de Bruijn, Csordas-Smith-Varga, Rodgers-Tao): dz_j/dt = 2 sum_{k != j} 1/(z_j - z_k).
Rewind (t < 0): close pairs attract and collide; a collision at t_c shows Lambda >= t_c (heuristic, not rigorous).
Pair formula with the rest as a frozen background S = sum_k 1/(m - z_k)^2:
   (d^2)' = 8 - 4 S d^2  =>  t_c = -ln(1/(1 - S d0^2/2)) / (4S)   when S d0^2 / 2 < 1.
Forward (t > 0): gaps relax toward sqrt(2/S), the evenly spaced state Rodgers-Tao prove at large height.
"""
import json, math, sys, time
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, "../unified-hypothesis")
import kern

g = np.array(json.load(open("../round9/zeros.json"))["gamma"])
z = 2 * g
N = len(z)
Z = z[-1]
dens = lambda w: np.log(np.asarray(w, float) / (4 * math.pi)) / (4 * math.pi)       # zeros per unit z

def tail_S(x):
    f = lambda w: dens(w) * (1 / (x - w)**2 + 1 / (x + w)**2)
    return quad(f, Z + 1.0, np.inf, limit=200)[0]
def tail_F(x):   # force 2*sum over the missing tail of 1/(x - w) + 1/(x + w)
    f = lambda w: dens(w) * (1 / (x - w) + 1 / (x + w))
    return 2 * quad(f, Z + 1.0, np.inf, limit=200)[0]

out = {}
# ---- pair collision times ----------------------------------------------------------------
allz = np.concatenate([z, -z])
mid = (z[:-1] + z[1:]) / 2; d0 = np.diff(z)
S = np.zeros(N - 1)
for i in range(0, N - 1, 500):
    m = mid[i:i + 500, None]
    D = m - allz[None, :]
    s = (1 / D**2).sum(1)
    # remove the pair itself
    n_ = len(s)
    s -= 1 / (mid[i:i+n_] - z[i:i+n_])**2 + 1 / (mid[i:i+n_] - z[i+1:i+1+n_])**2
    S[i:i + 500] = s
tgrid = np.linspace(z[0], Z, 60)
St = np.array([tail_S(x) for x in tgrid])
S += np.interp(mid, tgrid, St)
x = S * d0**2 / 2
tc = np.where(x < 1, -np.log(1 / np.maximum(1 - x, 1e-300)) / (4 * S), np.nan)
order = np.argsort(-np.nan_to_num(tc, nan=-1e9))
best = order[:10]
out["pairs"] = {"n_pairs_isolated": int(np.isfinite(tc).sum()),
                "top10": [{"j": int(j) + 1, "gamma": [float(g[j]), float(g[j+1])], "gap_z": float(d0[j]),
                           "S": float(S[j]), "x": float(x[j]), "t_c": float(tc[j])} for j in best]}
print("closest-to-0 collision", out["pairs"]["top10"][:3])

# ---- check 1: the pair formula against a direct rewind of H_t at low height (float64 ok to z ~ 75)
K = kern.riemann(500)
low = {}
for t in np.arange(-0.5, -12.01, -0.25):
    cz = kern.complex_zeros(K, 72.0, 14.0, t=float(t))
    nonreal = [c for c in cz if c.imag > 1e-4]
    if nonreal:
        low = {"t_first_nonreal": float(t), "zeros": [[float(c.real), float(c.imag)] for c in nonreal]}
        break
out["low_direct"] = low
# formula for low pairs, with the full background
lowpairs = []
for j in range(5):
    lowpairs.append({"gamma": [float(g[j]), float(g[j+1])], "x": float(x[j]), "t_c": float(tc[j]) if np.isfinite(tc[j]) else None})
# the first zero also pairs with its mirror -z_1
m0 = 0.0; d00 = 2 * z[0]
S0 = (1 / (m0 - allz)**2).sum() - 2 / z[0]**2 + tail_S(1e-9)
x0 = S0 * d00**2 / 2
out["low_formula"] = {"pairs": lowpairs, "mirror_pair": {"x": float(x0),
     "t_c": float(-math.log(1 / (1 - x0)) / (4 * S0)) if x0 < 1 else None}}
print("low direct", low, "\nlow formula", out["low_formula"])

# ---- check 2: full N-body rewind of the closest pair (local window moves, the rest frozen)
def nbody_rewind(j, half=150, dt0=1e-7):
    idx = np.arange(max(0, j - half), min(N, j + half + 2))
    w = z[idx].astype(complex).copy()
    frozen = np.concatenate([np.delete(z, idx), -z])
    tf = np.interp(w.real, tgrid, [tail_F(xx) for xx in tgrid]) if False else 0.0
    t = 0.0; a = list(idx).index(j)
    while True:
        D = w[:, None] - w[None, :]; np.fill_diagonal(D, np.inf)
        v = 2 * (1 / D).sum(1) + 2 * (1 / (w[:, None] - frozen[None, :])).sum(1)
        gap = abs(w[a + 1] - w[a])
        dt = min(1e-4, 0.002 * gap**2)
        w = w - dt * v; t -= dt
        if gap < 1e-4 or t < -0.05:
            return t, gap
j0 = int(best[0])
tN, gapN = nbody_rewind(j0)
out["nbody_closest"] = {"j": j0 + 1, "t_collision": tN, "formula": float(tc[j0]), "final_gap": float(gapN)}
print("nbody", out["nbody_closest"])

# ---- forward flow to t = 0.2: gaps even out --------------------------------------------
tgF = np.linspace(1.0, Z, 80)
TF = np.array([tail_F(xx) for xx in tgF])
w = z.copy(); t = 0.0
snaps = {}
targets = [0.0, 0.01, 0.05, 0.2]
def vel(w):
    v = np.zeros_like(w)
    allw = np.concatenate([w, -w])
    for i in range(0, N, 1000):
        D = w[i:i+1000, None] - allw[None, :]
        D[np.arange(D.shape[0]), np.arange(i, i + D.shape[0])] = np.inf
        v[i:i+1000] = 2 * (1 / D).sum(1)
    return v + np.interp(w, tgF, TF)
ti = 0
t0 = time.time()
while t < 0.2 + 1e-12:
    if ti < len(targets) and t >= targets[ti] - 1e-12:
        snaps[str(targets[ti])] = w.copy(); ti += 1
        if ti == len(targets): break
    gmin = np.diff(w).min()
    dt = min(1e-3, 0.02 * gmin**2, targets[ti] - t)
    k1 = vel(w); k2 = vel(w + dt * k1)
    w = w + dt * (k1 + k2) / 2; t += dt
print("forward done", time.time() - t0)
def norm_gaps(w, lo=1000, hi=4500):
    ww = w[lo:hi]; m = (ww[:-1] + ww[1:]) / 2
    return np.diff(ww) * dens(m)          # in units of the mean gap
fw = {}
for k_, ws in snaps.items():
    ng = norm_gaps(ws)
    fw[k_] = {"var": float(ng.var()), "min": float(ng.min()), "frac_below_half": float((ng < 0.5).mean())}
out["forward"] = fw
np.savez("rewind_snaps.npz", **{k_: v for k_, v in snaps.items()}, tc=tc, mid=mid, S=S, d0=d0)
json.dump(out, open("rewind.json", "w"), indent=1)
print(json.dumps(fw, indent=1))
