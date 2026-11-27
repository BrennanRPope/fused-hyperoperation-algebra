"""Which 3D disturbances grow on the shared 2D base? Split disturbances on the cover by their behaviour under the
straight shift tau and the glide g (both keep the base): four classes. K keeps (tau+,g+) and (tau-,g+);
TR keeps (tau+,g+) and (tau+,g-). usage: python3 classes.py seed out.json"""
import sys, json, numpy as np
from flow3d import Box, NS3, glide, shift, symmetrize
from wall import coeffs, field
def P(u, op, sign, N):
    if op == 'g':
        s = np.array([1.0, -1.0, 1.0]); img = np.array([s[i]*glide(u[i], N) for i in range(3)])
    else:
        img = np.array([shift(c, N) for c in u])
    return 0.5*(u + sign*img)
def run(seed, ts, gs, h=1.0, N=64, nu=0.01, T=5.0, a=0.001, E2=0.5):
    b = Box(N, h); ns = NS3(b, nu, 'T3')
    a0, a1 = coeffs(N, seed)
    u2 = symmetrize(symmetrize(field(b, a0, 0, 2, 4), 'K', N), 'TR', N)
    u3 = P(P(field(b, a1, 1, 1, 3), 't', ts, N), 'g', gs, N)
    U2 = b.leray(np.array([b.fwd(c) for c in u2]))*b.mask; U3 = b.leray(np.array([b.fwd(c) for c in u3]))*b.mask
    U2 *= np.sqrt(E2/ns.energy(U2)); U3 *= np.sqrt(a*E2/ns.energy(U3)); U = U2 + U3
    e0 = ns.energy3d(U); t = 0.0; tt, ee = [0.0], [1.0]
    while t < T - 1e-12:
        dt = min(ns.dt_cfl(U), T - t, 0.1); U = ns.step(U, dt); t += dt
        tt.append(t); ee.append(ns.energy3d(U)/e0)
    return dict(seed=seed, tau=ts, g=gs, h=h, growth=ee[-1], t=tt[::4], e=ee[::4])
if __name__ == '__main__':
    seed = int(sys.argv[1]); res = []
    for ts in (1, -1):
        for gs in (1, -1):
            r = run(seed, ts, gs); res.append(r); print(seed, ts, gs, round(r['growth'], 2), flush=True)
            json.dump(res, open(sys.argv[2], 'w'))
