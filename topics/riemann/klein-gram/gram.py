"""Gram-interval statistics for the torus, Klein-bottle and odd-half zetas, and Riemann zeta.

Torus and Riemann: sign changes of the real Hardy-type function Z(t) = e^{i theta(t)} F(1/2+it)
between consecutive Gram points (theta(g_n) = n pi).
Klein and odd half have no functional equation, so no real Z(t). We use the torus Gram points
(same Gamma factor pi^-s Gamma(s), same smooth zero count) and count zeros of any real part
in each Gram interval by the argument principle on rectangles [a, b] x [g_n, g_{n+1}],
for the strips  far-left (SL, -1),  middle (-1, 1),  right (1, SR).
Usage: python3 gram.py T   (writes gram_T.json)
"""
import json, sys, time, numpy as np
from zetas import parts, theta_T, theta_R, gram_points, zeta

SL, CUTS, SR = -3.5, (-1.0, 1.0), 2.5
XS = [SL, *CUTS, SR]

def combos(zt, e):
    return {'torus': zt, 'klein': zt/2 - e, 'odd': zt/2 + e}

def darg_path(pts, h_ok=0.5, depth=0):
    """Total change of arg along the polyline pts (complex array), for all three kinds at once,
    refining any step whose phase jump exceeds h_ok radians."""
    F = combos(*parts(pts))
    tot = {}
    for k, f in F.items():
        d = np.angle(f[1:]/f[:-1])
        bad = np.nonzero(np.abs(d) > h_ok)[0]
        tot[k] = d
    bad = set()
    for k in F:
        bad |= set(np.nonzero(np.abs(tot[k]) > h_ok)[0].tolist())
    if bad and depth < 12:
        for i in sorted(bad):
            sub = np.linspace(pts[i], pts[i+1], 9)
            r = darg_path(sub, h_ok, depth+1)
            for k in F: tot[k][i] = r[k]
    return {k: float(np.sum(v)) for k, v in tot.items()} if depth else {k: v for k, v in tot.items()}

def horizontal(g, h=0.01):
    """arg change from SR to SL along Im = g, split at the cuts: list per segment (right to left)."""
    out = []
    xs = XS[::-1]
    for a, b in zip(xs[:-1], xs[1:]):
        n = int(round((a - b)/h))
        pts = np.linspace(a, b, n + 1) + 1j*g
        d = darg_path(pts)
        out.append({k: float(np.sum(v)) for k, v in d.items()})
    return out

def vertical(x, gs, h=0.005):
    """arg change along Re = x from g_n to g_{n+1}, for each n."""
    res = {k: [] for k in ('torus', 'klein', 'odd')}
    for g0, g1 in zip(gs[:-1], gs[1:]):
        n = max(8, int(np.ceil((g1 - g0)/h)))
        pts = x + 1j*np.linspace(g0, g1, n + 1)
        d = darg_path(pts)
        for k in res: res[k].append(float(np.sum(d[k])))
    return {k: np.array(v) for k, v in res.items()}

def hardy_signs(theta, F, gs):
    return np.real(np.exp(1j*theta(gs))*F(0.5 + 1j*gs))

def sign_change_counts(theta, F, gs, h=0.004):
    """number of sign changes of Z(t) in each Gram interval (refined grid)."""
    counts = []
    for g0, g1 in zip(gs[:-1], gs[1:]):
        t = np.linspace(g0, g1, max(16, int((g1 - g0)/h)) + 1)
        z = np.real(np.exp(1j*theta(t))*F(0.5 + 1j*t))
        n = int(np.sum(np.sign(z[1:]) != np.sign(z[:-1])))
        # a close pair can hide between samples: refine round each local minimum of |Z| with no sign change
        a = np.abs(z)
        for i in np.nonzero((a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]) & (np.sign(z[:-2]) == np.sign(z[2:])) & (np.sign(z[1:-1]) == np.sign(z[:-2])))[0] + 1:
            tt = np.linspace(t[i-1], t[i+1], 2001)
            zz = np.real(np.exp(1j*theta(tt))*F(0.5 + 1j*tt))
            n += int(np.sum(np.sign(zz[1:]) != np.sign(zz[:-1])))
        counts.append(n)
    return np.array(counts)

def stats(c):
    c = np.asarray(c); n = len(c)
    return {'intervals': n, 'one': float(np.mean(c == 1)), 'empty': float(np.mean(c == 0)),
            'two': float(np.mean(c == 2)), 'three_plus': float(np.mean(c >= 3)), 'mean': float(np.mean(c))}

if __name__ == '__main__':
    T = float(sys.argv[1]); t0 = 10.0; tic = time.time()
    ns, gs = gram_points(theta_T, t0, T)
    print(f'torus Gram points in [{t0},{T}]: {len(gs)}', flush=True)
    H = [horizontal(g) for g in gs]
    print('horizontals done', round(time.time()-tic), flush=True)
    V = {x: vertical(x, gs) for x in XS}
    print('verticals done', round(time.time()-tic), flush=True)
    out = {'T': T, 'gram_n': ns.tolist(), 'gram_t': gs.tolist(), 'strips': {}}
    xs = XS[::-1]  # segment j of H runs xs[j] -> xs[j+1] (right to left)
    for j in range(3):
        a, b = xs[j+1], xs[j]          # strip a < Re < b
        name = ['right', 'middle', 'far_left'][j]
        cnt = {}
        for k in ('torus', 'klein', 'odd'):
            bottom = -np.array([H[i][j][k] for i in range(len(gs)-1)])     # a -> b along g_n
            top = np.array([H[i+1][j][k] for i in range(len(gs)-1)])       # b -> a along g_{n+1}
            w = (bottom + V[b][k] + top - V[a][k])/(2*np.pi)
            cnt[k] = np.rint(w).astype(int).tolist()
            assert np.max(np.abs(w - np.rint(w))) < 0.05, (name, k, np.max(np.abs(w - np.rint(w))))
        out['strips'][name] = {'a': a, 'b': b, 'counts': cnt}
    # torus and Riemann sign changes
    tor_sc = sign_change_counts(theta_T, lambda s: parts(s)[0], gs)
    out['torus_sign_changes'] = tor_sc.tolist()
    out['torus_gram_sign_ok'] = float(np.mean(((-1.0)**ns)*hardy_signs(theta_T, lambda s: parts(s)[0], gs) > 0))
    nR, gR = gram_points(theta_R, t0, None or float(sys.argv[2]) if len(sys.argv) > 2 else T)
    rz_sc = sign_change_counts(theta_R, zeta, gR)
    out['riemann'] = {'T': float(gR[-1]), 'gram_t': gR.tolist(), 'sign_changes': rz_sc.tolist(),
                      'gram_sign_ok': float(np.mean(((-1.0)**nR)*hardy_signs(theta_R, zeta, gR) > 0))}
    print('done', round(time.time()-tic), flush=True)
    json.dump(out, open(f'gram_{int(T)}.json', 'w'))
    for k in ('torus', 'klein', 'odd'):
        mid = np.array(out['strips']['middle']['counts'][k]); rt = np.array(out['strips']['right']['counts'][k])
        fl = np.array(out['strips']['far_left']['counts'][k])
        print(k, 'middle', stats(mid), 'middle+right', stats(mid+rt), 'far-left per interval', fl.mean())
    print('torus sign', stats(tor_sc), 'riemann sign', stats(rz_sc))
