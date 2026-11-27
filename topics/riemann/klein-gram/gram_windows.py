"""Gram statistics by height window, against Riemann zeta at the height with the same zero spacing,
and the mean Gram value (-1)^n Re[e^{i theta} F(1/2 + i g_n)] for each zeta.
Matched height: zeta's zero density (1/2pi) log(H/2pi) equals the torus density (1/pi) log(t/pi)
when H = 2 t^2 / pi."""
import json, numpy as np
from zetas import parts, theta_T, theta_R, gram_points
from riemann_siegel import Z_rs
from gram import stats

G = json.load(open('gram_2000.json'))
gt = np.array(G['gram_t']); gn = np.array(G['gram_n'])
mid = G['strips']['middle']['counts']; rt = G['strips']['right']['counts']
cnt = {k: np.array(mid[k]) + np.array(rt[k]) for k in ('torus', 'klein', 'odd')}
WIN = [(10, 500), (500, 1000), (1000, 1500), (1500, 2000)]

def zeta_window(H, n_int=3000):
    """Gram stats for zeta over n_int consecutive Gram intervals starting near height H."""
    L = 2*np.pi/np.log(H/(2*np.pi))
    ns, gs = gram_points(theta_R, H, H + L*(n_int + 2))
    gs = gs[:n_int + 1]
    counts = []
    for g0, g1 in zip(gs[:-1], gs[1:]):
        t = np.linspace(g0, g1, 81)
        z = Z_rs(t)
        n = int(np.sum(np.sign(z[1:]) != np.sign(z[:-1])))
        a = np.abs(z)   # a close pair can hide between samples: refine round local minima of |Z|
        for i in np.nonzero((a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]) & (np.sign(z[:-2]) == np.sign(z[2:])) & (np.sign(z[1:-1]) == np.sign(z[:-2])))[0] + 1:
            zz = Z_rs(np.linspace(t[i-1], t[i+1], 401))
            n += int(np.sum(np.sign(zz[1:]) != np.sign(zz[:-1])))
        counts.append(n)
    vals = ((-1.0)**ns[:n_int])*Z_rs(gs[:n_int])
    return np.array(counts), float(np.mean(vals))

out = {'windows': []}
FT = parts(0.5 + 1j*gt)
vals = {'torus': FT[0], 'klein': FT[0]/2 - FT[1], 'odd': FT[0]/2 + FT[1]}
gv = {k: ((-1.0)**gn)*np.real(np.exp(1j*theta_T(gt))*v) for k, v in vals.items()}
out['mean_gram_value'] = {k: float(np.mean(v)) for k, v in gv.items()}
out['gram_sign_ok'] = {k: float(np.mean(v > 0)) for k, v in gv.items()}
out['gram_values'] = {k: v.tolist() for k, v in gv.items()}
for a, b in WIN:
    sel = (gt[:-1] >= a) & (gt[:-1] < b)
    tm = 0.5*(a + b); H = 2*tm*tm/np.pi
    zc, zmean = zeta_window(H)
    row = {'torus_window': [a, b], 'zeta_height': H, 'zeta': stats(zc), 'zeta_mean_gram_value': zmean}
    for k in cnt: row[k] = stats(cnt[k][sel])
    out['windows'].append(row)
    print(f"[{a},{b}] zeta@{H:.0f}: one {row['zeta']['one']:.3f} empty {row['zeta']['empty']:.3f} | " +
          " | ".join(f"{k} one {row[k]['one']:.3f} empty {row[k]['empty']:.3f} mean {row[k]['mean']:.3f}" for k in cnt), flush=True)
print('mean Gram value', out['mean_gram_value']); print('weak Gram law (sign) holds', out['gram_sign_ok'])
json.dump(out, open('gram_windows.json', 'w'))
