"""Height-dial scan: a strong 2D base flow plus a 3D part of relative energy a, on a box of height h.
For each surface (K x S1, TR x S1, optionally T3), height h and amplitude a, record the enstrophy
gain G = max_t Z(t)/Z(0), the 3D energy fraction and max|omega|. In 2D (a = 0) Z can only fall,
so G > 1 means the third direction has taken over from viscosity: control in the enstrophy sense is lost.
usage: python3 wall.py N nu T seed surfaces heights amps out.json
"""
import sys, json, time
import numpy as np
from flow3d import Box, NS3, symmetrize, glide_defect

def coeffs(N, seed, N0=64):
    """Random coefficients drawn on a fixed 64 x 64 wavenumber grid, so the same seed gives the same
    initial flow at every resolution (only |kx|, |ky| <= 4 are used)."""
    rng = np.random.default_rng(seed)
    a0 = rng.standard_normal((3, N0, N0)) + 1j*rng.standard_normal((3, N0, N0))
    a1 = rng.standard_normal((3, N0, N0)) + 1j*rng.standard_normal((3, N0, N0))
    if N == N0:
        return a0, a1
    k0 = np.fft.fftfreq(N0, 1.0/N0).astype(int); k = np.fft.fftfreq(N, 1.0/N).astype(int)
    out = []
    for a in (a0, a1):
        b = np.zeros((3, N, N), complex)
        idx = {kk: i for i, kk in enumerate(k)}
        for i, kx in enumerate(k0):
            for j, ky in enumerate(k0):
                if abs(kx) <= 8 and abs(ky) <= 8 and kx in idx and ky in idx:
                    b[:, idx[kx], idx[ky]] = a[:, i, j]
        out.append(b)
    return tuple(out)

def field(box, a, kzi, kmin, kmax):
    b = box
    U = np.zeros((3,) + b.KX.shape, complex)
    U[..., kzi] = a
    kh = np.sqrt(b.KX**2 + b.KY**2)
    U *= (kh >= kmin) & (kh <= kmax)
    if kzi == 0:
        U[2] = 0
    U = b.leray(U)*b.mask
    return np.array([b.inv(c) for c in U])

def initial(box, surf, seed, a, E2=0.5, joint=False):
    N = box.N
    a0, a1 = coeffs(N, seed)
    ns = NS3(box, 0.0, surf)
    parts = []
    for coef, kzi, (kmin, kmax) in ((a0, 0, (2, 4)), (a1, 1, (1, 3))):
        u = symmetrize(field(box, coef, kzi, kmin, kmax), surf, N)
        if joint and kzi == 0:   # base flow kept by both gluings: identical on K and TR
            u = symmetrize(symmetrize(field(box, coef, kzi, kmin, kmax), 'K', N), 'TR', N)
        U = box.leray(np.array([box.fwd(c) for c in u]))*box.mask
        parts.append(U)
    U2, U3 = parts
    U2 *= np.sqrt(E2/ns.energy(U2))
    U3 *= np.sqrt(a*E2/ns.energy(U3)) if a > 0 else 0.0
    return U2 + U3

def run(N, nu, T, seed, surf, h, a, joint=False):
    box = Box(N, h)
    ns = NS3(box, nu, surf)
    U = initial(box, surf, seed, a, joint=joint)
    Z0 = ns.enstrophy(U); E0 = ns.energy(U)
    t = 0.0; Zmax = Z0; tZ = 0.0; wm0 = ns.wmax(U); wmax = wm0
    ts, Zs, E3s, E3a = [0.0], [Z0], [ns.energy3d(U)/E0], [ns.energy3d(U)]
    while t < T - 1e-12:
        dt = min(ns.dt_cfl(U), T - t, 0.1)
        U = ns.step(U, dt); t += dt
        Z = ns.enstrophy(U)
        if Z > Zmax: Zmax, tZ = Z, t
        ts.append(t); Zs.append(Z); e3 = ns.energy3d(U); E3a.append(e3); E3s.append(e3/ns.energy(U))
        if len(ts) % 5 == 0:
            wmax = max(wmax, ns.wmax(U))
    out = dict(N=N, Nz=box.Nz, nu=nu, T=T, seed=seed, surf=surf, h=h, a=a, E0=E0, Z0=Z0,
               G=Zmax/Z0, tG=tZ, wmax0=wm0, wmax=wmax, Ef=ns.energy(U),
               t=ts[::max(1, len(ts)//80)], Z=[z/Z0 for z in Zs][::max(1, len(ts)//80)],
               E3=E3s[::max(1, len(ts)//80)], E3abs=E3a[::max(1, len(ts)//80)], joint=joint)
    if surf == 'K':
        u = np.array([box.inv(c) for c in U])
        out['glide_defect'] = glide_defect(u, N)
    return out

if __name__ == '__main__':
    N, nu, T, seed = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    surfs = sys.argv[5].split(','); hs = [float(x) for x in sys.argv[6].split(',')]
    amps = [float(x) for x in sys.argv[7].split(',')]; out = sys.argv[8]
    joint = len(sys.argv) > 9 and sys.argv[9] == 'joint'
    res = []
    for h in hs:
        for a in amps:
            for s in surfs:
                t0 = time.time()
                r = run(N, nu, T, seed, s, h, a, joint)
                r['secs'] = time.time() - t0
                res.append(r)
                print(f"h={h} a={a} {s}: G={r['G']:.3f} tG={r['tG']:.2f} E3f={r['E3'][-1]:.3f} "
                      f"wmax={r['wmax']:.2f} ({r['secs']:.0f}s)" + (f" gd={r.get('glide_defect',0):.1e}" if s=='K' else ''), flush=True)
                json.dump(res, open(out, 'w'))
