"""Round 6, lattice-side pairs not covered by the zero runs:
 s+f: the height dial on the face-power lattice zeta (a-points F = X_s near the zeros), F normalised to -> 1;
 flux pass: lattice-point discrepancy N(R) - area R^2 for the face-power and roundness gauges (its growth exponent)."""
import json, numpy as np
from lattice_engine import *
res = {}
Xs = {'0.8': -0.14189434629030684, '0.95': -0.031641082794077846, '1.0': 0.0}
d = json.load(open('lattice_zeros.json'))
sf = {}
for name, (B, N), p in [('square', square_setup(), 2.0), ('hexagon', hex_setup(), 3.0), ('hexagon', hex_setup(), 4.0)]:
    Z = LatticeZeta(B, N, p)
    gmin2 = Z.gmin2; w = int(np.sum(np.abs(Z.g2 - gmin2) < 1e-9))
    F = lambda s: Z(np.atleast_1d(s))*np.exp(np.atleast_1d(s)*np.log(gmin2))/w
    zs = [complex(*z) for z in d[name][str(p)]['zeros'][:6]]
    out = {}
    for sh, X in Xs.items():
        pts = []
        for z0 in zs:
            z = z0 + 0.02
            for _ in range(60):
                f = F(z)[0] - X; h = 1e-6
                dz = f/((F(z + h)[0] - F(z - h)[0])/(2*h)); z -= dz
                if abs(dz) < 1e-11: break
            if abs(F(z)[0] - X) < 1e-8: pts.append([z.real, z.imag, abs(z - z0)])
        out[sh] = dict(a_points=pts, max_off_line=max(abs(q[0] - 0.5) for q in pts) if pts else None)
    sf[f'{name}_p{p}'] = out
    print(name, p, {k: (len(v['a_points']), v['max_off_line']) for k, v in out.items()}, flush=True)
res['sf_height_a_points'] = sf

flux = {}
Rs = np.linspace(40, 900, 400)
for name, (B, N), p in [('square', square_setup(), 2.0), ('square', square_setup(), 3.0), ('square', square_setup(), 4.0),
                        ('square', square_setup(), np.inf), ('hexagon', hex_setup(), 3.0), ('hexagon', hex_setup(), 4.0),
                        ('hexagon', hex_setup(), 6.0), ('hexagon', hex_setup(), np.inf)]:
    M = 1400
    a, c = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1))
    pts = np.stack([a.ravel(), c.ravel()], 1) @ B
    g = np.sort(gauge(pts, N, p)); g = g[g > 0]
    th = np.linspace(0, 2*np.pi, 200001)[:-1]
    G = gauge(np.stack([np.cos(th), np.sin(th)], 1), N, p)
    area = 0.5*np.mean(G**-2.0)*2*np.pi/abs(np.linalg.det(B))
    E = np.searchsorted(g, Rs, side='right') - area*Rs**2
    # growth exponent of the root-mean-square discrepancy over dyadic windows
    rms = []; mids = []
    for lo in [40, 80, 160, 320, 640]:
        m = (Rs >= lo) & (Rs < 2*lo) if lo < 640 else (Rs >= lo)
        rms.append(np.sqrt(np.mean(E[m]**2))); mids.append(np.sqrt(lo*min(2*lo, 900)))
    slope = np.polyfit(np.log(mids), np.log(rms), 1)[0]
    key = f"{name}_p{'inf' if np.isinf(p) else p}"
    flux[key] = dict(rms_exponent=float(slope), rms=list(map(float, rms)))
    print(key, round(slope, 3), flush=True)
res['flux_discrepancy'] = flux
json.dump(res, open('lattice_pairs.json', 'w'), indent=1)
print('done')
