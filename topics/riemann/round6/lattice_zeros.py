"""Round 6: zeros of the lattice zeta as the roundness dial (square) and the face-power dial (hexagon) run from 2 to infinity,
plus the curvature runs (lattice shape) and a scale check for the ratio dial."""
import json, time, numpy as np
from lattice_engine import *
res = {}
PS = [2.0, 2.1, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0, 20.0, 40.0, np.inf]
def run(name, B, N, plist):
    out = {}
    for p in plist:
        t0 = time.time()
        Z = LatticeZeta(B, N, p)
        zs, (S, T, V) = find_zeros(Z)
        out['inf' if np.isinf(p) else str(p)] = dict(zeros=[[z.real, z.imag] for z in zs], A=Z.A, fit_resid=Z.fit_resid)
        print(name, p, len(zs), round(time.time() - t0, 1), [(round(z.real, 4), round(z.imag, 3)) for z in zs[:6]], flush=True)
        json.dump(res, open('lattice_zeros.json', 'w'), indent=1)
    return out
res['square'] = run('square', *square_setup(), PS)
res['hexagon'] = run('hexagon', *hex_setup(), PS)
# curvature: rectangle lattice Z x sqrt5 Z (x^2 + 5y^2, discriminant -20, class number 2) with the square's normals
B5 = np.array([[1.0, 0.0], [0.0, np.sqrt(5.0)]])
res['disc-20_square'] = run('disc-20', B5, square_setup()[1], [2.0, 3.0, 4.0, np.inf])
# curvature with face power: lattice x^2 + xy + 6y^2 (discriminant -23, class number 3) with its own Voronoi hexagon
b1 = np.array([1.0, 0.0]); b2 = np.array([0.5, np.sqrt(23)/2])
V = [b1, b2, b2 - b1]; Nv = np.array([2*v/np.dot(v, v) for v in V] + [-2*v/np.dot(v, v) for v in V])
res['disc-23_voronoi'] = run('disc-23', np.array([b1, b2]), Nv, [2.0, 3.0, 4.0, np.inf])
# ratio dial: the hexagon at p = 3 on a lattice scaled by 1.7 has the same zeros
res['hexagon_scaled_1.7'] = run('hex-scaled', 1.7*hex_setup()[0], hex_setup()[1]/1.7, [3.0])
json.dump(res, open('lattice_zeros.json', 'w'), indent=1)
print('done')
