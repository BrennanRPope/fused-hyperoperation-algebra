"""Blow-up hunt inside the twisted class: an antiparallel pair of x-directed vortex tubes (Kerr-type).
On the cover [0,2pi)^3 the glide g(x,y,z)=(x+pi,-y,z) maps an x-tube at height y to an ANTIPARALLEL x-tube
at -y, so the pair is built as tube A plus its glide image. With a symmetric (Kerr) bend of even wavenumber
the pair is also tau-periodic (tau: x -> x+pi), so it lies in both the Klein class K x S1 and the torus class
TR x S1. A small 3D seed is then added in one of the four (tau +-, g +-) classes:
  (t-,g+) = K only   (period-doubled, mirror symmetric)    -> flow lives on K x S1
  (t+,g-) = TR only  (period pi, mirror broken)            -> flow lives on TR x S1
  (t+,g+) = shared, (t-,g-) = neither (T3 only), 'none' = base alone.
usage: python3 tubes.py N nu T cls seed out.json [eps] [a] [kbend]
"""
import sys, json, time, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from flow3d import Box, NS3, glide, shift
from wall import coeffs
from scipy import fft as sfft
_W = int(os.environ.get("WORKERS", 1))
if _W > 1:
    Box.fwd = lambda self, f: sfft.rfftn(f, axes=(0, 1, 2), workers=_W)
    Box.inv = lambda self, F: sfft.irfftn(F, s=self.shape, axes=(0, 1, 2), workers=_W)

SIG = float(os.environ.get("SIG", 0.30)); DELTA = float(os.environ.get("DELTA", 0.55)); GAMMA = float(os.environ.get("GAMMA", 2.0))

def P(u, op, sign, N):
    if op == 'g':
        s = np.array([1.0, -1.0, 1.0]); img = np.array([s[i]*glide(u[i], N) for i in range(3)])
    else:
        img = np.array([shift(c, N) for c in u])
    return 0.5*(u + sign*img)

def tube_A(box, a, kbend):
    """vorticity of one x-tube: centre y = DELTA - a cos(k x), z = pi + a cos(k x)."""
    N, Nz = box.N, box.Nz
    x = np.arange(N)*2*np.pi/N; z = np.arange(Nz)*2*np.pi*box.h/Nz
    X, Y, Z = np.meshgrid(x, x, z, indexing='ij')
    yc = DELTA - a*np.cos(kbend*X); zc = np.pi + a*np.cos(kbend*X)
    dy = a*kbend*np.sin(kbend*X); dz = -a*kbend*np.sin(kbend*X)
    w = np.zeros((3,) + X.shape)
    for iy in (-1, 0, 1):
        for iz in (-1, 0, 1):
            r2 = (Y - yc + 2*np.pi*iy)**2 + (Z - zc + 2*np.pi*iz*box.h)**2
            w[0] += np.exp(-r2/SIG**2)
    amp = GAMMA/(np.pi*SIG**2)
    nrm = np.sqrt(1 + dy**2 + dz**2)
    return amp*np.array([w[0]/nrm, w[0]*dy/nrm, w[0]*dz/nrm])

def vel_from_vort(box, w):
    W = np.array([box.fwd(c) for c in w])
    KX, KY, KZ = box.KX, box.KY, box.KZ
    U = 1j*np.array([KY*W[2] - KZ*W[1], KZ*W[0] - KX*W[2], KX*W[1] - KY*W[0]])*box.K2i
    return box.leray(U)*box.mask

def initial(box, cls, seed, eps, a, kbend):
    N = box.N; ns = NS3(box, 0.0, 'T3')
    wA = tube_A(box, a, kbend)
    s = np.array([1.0, -1.0, 1.0])
    # vorticity is a pseudovector: image vorticity = det(Dg) Dg w(g p) = -Dg w(g p)
    wB = np.array([-s[i]*glide(wA[i], N) for i in range(3)])
    U = vel_from_vort(box, wA + wB)
    if cls == 'none' or eps == 0:
        return U
    ts, gs = {'t-g+': (-1, 1), 't+g-': (1, -1), 't+g+': (1, 1), 't-g-': (-1, -1)}[cls]
    rng = np.random.default_rng(1000 + seed)
    R = (rng.standard_normal((3,) + box.KX.shape) + 1j*rng.standard_normal((3,) + box.KX.shape))
    kk = np.sqrt(box.K2)
    R *= (kk >= 1) & (kk <= 4)
    R = box.leray(R)*box.mask
    r = np.array([box.inv(c) for c in R])
    r = P(P(r, 't', ts, N), 'g', gs, N)
    Rs = box.leray(np.array([box.fwd(c) for c in r]))*box.mask
    Rs *= np.sqrt(eps*ns.energy(U)/ns.energy(Rs))
    return U + Rs

def tail(ns, U):
    """fraction of energy in the outer shell 0.8*kmax < |k| <= kmax (resolution monitor)."""
    b = ns.b; kk = np.sqrt(b.K2); kmax = b.N/3
    V = U*((kk > 0.8*kmax) & (b.mask > 0))
    return sum(ns._wsum(c) for c in V)/max(sum(ns._wsum(c) for c in U), 1e-300)

def defects(box, U):
    N = box.N; u = np.array([box.inv(c) for c in U]); m = np.abs(u).max()
    s = np.array([1.0, -1.0, 1.0])
    dg = np.abs(u - np.array([s[i]*glide(u[i], N) for i in range(3)])).max()/m
    dt = np.abs(u - np.array([shift(c, N) for c in u])).max()/m
    return dg, dt

def run(N, nu, T, cls, seed, eps=0.01, a=0.15, kbend=2, every=float(os.environ.get("EVERY", 0.25))):
    box = Box(N, 1.0, Nz=N); ns = NS3(box, nu, 'T3')
    U = initial(box, cls, seed, eps, a, kbend)
    E0, Z0, W0 = ns.energy(U), ns.enstrophy(U), ns.wmax(U)
    dg0, dt0 = defects(box, U)
    t = 0.0; nxt = 0.0
    rec = dict(t=[], wmax=[], Z=[], E=[], tail=[])
    t0 = time.time()
    while t < T - 1e-12:
        if t >= nxt - 1e-12:
            rec['t'].append(t); rec['wmax'].append(ns.wmax(U)); rec['Z'].append(ns.enstrophy(U))
            rec['E'].append(ns.energy(U)); rec['tail'].append(tail(ns, U)); nxt += every
        dt = min(ns.dt_cfl(U, 0.8), T - t, 0.05)
        U = ns.step(U, dt); t += dt
    rec['t'].append(t); rec['wmax'].append(ns.wmax(U)); rec['Z'].append(ns.enstrophy(U))
    rec['E'].append(ns.energy(U)); rec['tail'].append(tail(ns, U))
    dg1, dt1 = defects(box, U)
    i = int(np.argmax(rec['wmax']))
    return dict(N=N, nu=nu, T=T, cls=cls, seed=seed, eps=eps, a=a, kbend=kbend, E0=E0, Z0=Z0, W0=W0,
                wpeak=rec['wmax'][i], tpeak=rec['t'][i], Zpeak=max(rec['Z']), tailmax=max(rec['tail']),
                defects0=[dg0, dt0], defects1=[dg1, dt1], secs=time.time() - t0, **rec)

if __name__ == '__main__':
    N, nu, T, cls, seed, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6]
    eps = float(sys.argv[7]) if len(sys.argv) > 7 else 0.01
    a = float(sys.argv[8]) if len(sys.argv) > 8 else 0.15
    kb = int(sys.argv[9]) if len(sys.argv) > 9 else 2
    r = run(N, nu, T, cls, seed, eps, a, kb)
    json.dump(r, open(out, 'w'))
    print(cls, seed, 'wpeak %.3f/%.3f at t=%.2f  Zpeak/Z0 %.3f  tail %.2e  defects0 %s defects1 %s  %.0fs' % (
        r['wpeak'], r['W0'], r['tpeak'], r['Zpeak']/r['Z0'], r['tailmax'], np.round(r['defects0'], 3), np.round(r['defects1'], 3), r['secs']))
