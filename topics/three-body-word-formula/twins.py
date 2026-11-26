"""DRAFT. Are the same-word 'twin' orbits new, or the catalog orbit seen through a symmetry?
Symmetries of the problem: relabel the bodies (6), reflect, reverse time, rotate, scale, shift in time.
Fingerprints that none of these change: the path on the shape sphere up to the 12 relabel/reflect images,
the sorted swells and dips of the size R (scaled by |E|), and the Floquet multipliers. (The action is no test:
for any periodic orbit the virial theorem makes it 3 T*, so equal periods force equal actions.)"""
import json, os, sys, ctypes, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from scipy.spatial import cKDTree
from tbw import ic, flow, polish, raw
if not os.path.exists(os.path.join(H, "libhyp.so")):  # the compiled file is not kept in the repo; build it from hyp.c
    import subprocess
    subprocess.run(["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-o", os.path.join(H, "libhyp.so"), os.path.join(H, "hyp.c"), "-lm"], check=True)
lib = ctypes.CDLL(os.path.join(H, 'libhyp.so')); Dd = ctypes.c_double
lib.sample.argtypes = [Dd, Dd, Dd, ctypes.c_int, Dd, ctypes.POINTER(Dd)]
def pos(p1, p2, T, n):
    out = (Dd * (6 * (n + 1)))(); assert lib.sample(p1, p2, T, n, 1e-12, out) == 0
    return np.frombuffer(out, float).reshape(-1, 3, 2).copy()
def shape(r):
    z = r[..., 0] + 1j * r[..., 1]; z1 = (z[:, 1] - z[:, 0]) / np.sqrt(2); z2 = np.sqrt(2 / 3) * (z[:, 2] - .5 * (z[:, 0] + z[:, 1]))
    R2 = abs(z1)**2 + abs(z2)**2; w = np.conj(z1) * z2
    return np.stack([(abs(z1)**2 - abs(z2)**2) / R2, 2 * w.real / R2, 2 * w.imag / R2], 1)
def images(S):
    out = []
    for k in range(3):
        th = 2 * np.pi * k / 3; M = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        for a in (1, -1):
            for b in (1, -1): out.append((S * np.array([1, a, b])) @ M.T)
    return out
def gap(A, B): return max(cKDTree(A).query(B)[0].max(), cKDTree(B).query(A)[0].max())
def Rext(r, E):
    R = np.sqrt((r**2).sum((1, 2))) * abs(E); i = np.arange(1, len(R) - 1)
    mx = R[i][(R[i] > R[i - 1]) & (R[i] >= R[i + 1])]; mn = R[i][(R[i] < R[i - 1]) & (R[i] <= R[i + 1])]
    return np.sort(mx), np.sort(mn)
def floq(p1, p2, T, eps=1e-7):
    X0 = ic(p1, p2); M = np.zeros((12, 12))
    for j in range(12):
        e = np.zeros(12); e[j] = eps; M[:, j] = (flow(X0 + e, T) - flow(X0 - e, T)) / (2 * eps)
    return np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
def fp(p1, p2, T, n):
    E = -2.5 + 3 * (p1**2 + p2**2); r = pos(p1, p2, T, n)
    return dict(S=shape(r), R=Rext(r, E), F=floq(p1, p2, T), Ts=T * abs(E)**1.5, E=E)
def compare(a, b, n=40000):
    A = fp(*a, n); B = fp(*b, n); A2 = fp(*a, 2 * n)
    g = min(gap(A['S'], X) for X in images(B['S'])); floor = gap(A['S'], A2['S'])
    same_R = len(A['R'][0]) == len(B['R'][0]) and len(A['R'][1]) == len(B['R'][1])
    dR = float(max(np.max(abs(A['R'][0] - B['R'][0])), np.max(abs(A['R'][1] - B['R'][1])))) if same_R else None
    return dict(Ts=[A['Ts'], B['Ts']], shape_gap=float(g), sampling_floor=float(floor), R_counts=[[len(x) for x in A['R']], [len(x) for x in B['R']]], R_diff=dR,
                maxlam=[float(A['F'][0]), float(B['F'][0])], mods=[A['F'].tolist(), B['F'].tolist()])
if __name__ == '__main__':
    sys.argv = sys.argv[:1]; import figures as F
    d = {o['cls']: o for o in F.d}; T = F.tally()
    R = json.load(open(os.path.join(H, 'predict2_sample.json'))); Q = json.load(open(os.path.join(H, 'retry.json')))
    pairs = [('I.A1', (0.3471168881, 0.5327249454, 6.3259139829), (0.3394074119878319, 0.5361849781525929, 6.290818765065431))]
    for c, v in T.items():
        if v['v'] != 'twin': continue
        fs = [f for rr in R[c]['rotations'] for f in rr['finals']] + Q.get(c, {}).get('finals', [])
        f = [f for f in fs if f['res'] < 1e-7 and f.get('word_ok')][0]; o = d[c]
        q, res = polish(f['v1'], f['v2'], f['T']); pairs.append((c, (o['v1'], o['v2'], o['T']), (float(q[0]), float(q[1]), float(q[2]))))
    out = {}
    for c, a, b in pairs:
        qa, ra = polish(*a); a = tuple(float(x) for x in qa)
        r = compare(a, b); r['cat'] = a; r['twin'] = b; r['res_cat'] = ra
        r['word_cat'] = raw(a[0], a[1], a[2] * 0.999999, maxev=600)[0]; r['word_twin'] = raw(b[0], b[1], b[2] * 0.999999, maxev=600)[0]
        # any other catalog orbit with the same scale-free period (catalog gives 3 decimals)
        r['same_period_in_catalog'] = [x['cls'] for x in F.d if abs(x['Tstar'] - r['Ts'][1]) < 2e-3]
        out[c] = r; json.dump(out, open(os.path.join(H, 'twins.json'), 'w'), indent=1)
        print(c, 'T* %.7f %.7f' % tuple(r['Ts']), 'gap %.3f floor %.3f' % (r['shape_gap'], r['sampling_floor']), 'R', r['R_counts'], r['R_diff'], 'lam %.4g %.4g' % tuple(r['maxlam']), r['same_period_in_catalog'], r['word_cat'] == r['word_twin'], flush=True)
