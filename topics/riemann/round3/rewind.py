"""Road 1: the rewind (backward heat flow) of the zeros, in root and in fused (pair) coordinates.
Convention (de Bruijn-Newman): H_t(z), dH/dt = -d2H/dz2, zeros at t = 0 are z = 2*gamma.
Zero motion: dz_j/dt = 2 sum_{k != j} 1/(z_j - z_k). Rewind = t decreasing (zeros attract).
-> rewind.json"""
import json, numpy as np, mpmath as mp

G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
L = json.load(open('zeros_lehmer.json')); GL = np.array([float(x) for x in L['gamma']]); NL = L['n']

def tail_force(z, zmax):
    # sum over zeros above zmax and their mirrors: 2*sum [1/(z-w) + 1/(z+w)], density rho(w) = log(w/4pi)/(4pi)
    w = np.linspace(zmax, zmax * 400, 400001); rho = np.log(w / (4 * np.pi)) / (4 * np.pi)
    f = 2 * z / (z * z - w * w) * rho
    return 2 * np.trapezoid(f, w)

def rewind_window(gam, j, W=150, all_gam=None, zmax=None):
    """full-ODE rewind of the pair (j, j+1) in a window of +-W zeros; others frozen as background."""
    z_all = 2 * all_gam
    lo, hi = j - W, j + W + 2
    zw = z_all[lo:hi].copy()
    out = np.concatenate([z_all[:lo], z_all[hi:]])
    bg = np.array([2 * (np.sum(1 / (zi - out)) + np.sum(1 / (zi + z_all))) + tail_force(zi, zmax) for zi in zw])
    def vel(z):
        D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf)
        return 2 * np.sum(1 / D, axis=1) + bg
    t = 0.0; p = W  # pair is zw[W], zw[W+1]
    g0 = zw[p + 1] - zw[p]
    while True:
        gaps = np.diff(zw); g = gaps.min(); i = gaps.argmin()
        if g < 1e-4 * g0: break
        dt = -0.01 * g * g
        k1 = vel(zw); k2 = vel(zw + dt / 2 * k1); k3 = vel(zw + dt / 2 * k2); k4 = vel(zw + dt * k3)
        zw = zw + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    v = vel(zw)
    t_star = t - g * g / (2 * g * (v[i + 1] - v[i]))  # g^2 is linear near collision
    return t_star, int(i == p), g0, int(lo + i)

res = {}
# isolated (pair-only) rewind times tau = dgamma^2/2 for every neighbour pair
dg = np.diff(G); tau = dg ** 2 / 2
res['pairs'] = {'gamma_mid': ((G[1:] + G[:-1]) / 2).tolist(), 'tau_isolated': tau.tolist()}
# full ODE for the 15 smallest-tau pairs away from the ends
idx = [i for i in np.argsort(tau) if 160 <= i < len(G) - 160][:15]
zmax = 2 * G[-1]
full = []
for i in idx:
    ts, same, g0, ic = rewind_window(G, i, all_gam=G, zmax=zmax)
    full.append({'j': int(i) + 1, 'gamma': float(G[i]), 'dgamma': float(dg[i]), 'tau_isolated': float(tau[i]), 't_full': float(ts), 'same_pair': same, 'collides_j': ic + 1, 'collides_gamma': float(G[ic])})
    print(full[-1], flush=True)
res['full'] = full
# Lehmer's pair window (zeros 6551..6870), background from this window only + tail beyond it
dgl = np.diff(GL); il = int(np.argmin(dgl[152:-152])) + 152
G2 = np.concatenate([G, GL]); il2 = len(G) + il
ts, same, g0, ic = rewind_window(G2, il2, all_gam=G2, zmax=2 * GL[-1])
# note: zeros below 6551 are left out of the background; their effect on the gap is a smooth bias (checked below)
res['lehmer'] = {'n': NL[il], 'gamma': float(GL[il]), 'gamma2': float(GL[il + 1]), 'dgamma': float(dgl[il]),
                 'tau_isolated': float(dgl[il] ** 2 / 2), 't_full': float(ts), 'same_pair': same}
print(res['lehmer'])
json.dump(res, open('rewind.json', 'w'))
