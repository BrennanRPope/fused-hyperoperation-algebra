"""DRAFT. Read the Klein-class antiparallel pair (klein-flow/periodic-box/blowup) through the two-plane picture.
Same initial data and solver as tubes.py (GAMMA=8 DELTA=0.5 SIG=0.3, a=0.25, kbend=2, no seed).
Tracks, every 0.1:
  Z        total enstrophy (total spin)
  plane    share of enstrophy in the kx=0 Fourier plane = the straight pair, a flat flow on the yz-plane
  P_cross  enstrophy production <w.(w.grad)u> (all of it needs the bend: a flat flow alone has none)
  bend_self / bend_plane: production split by whether the stretched spin and the stretching gradient
           come from the flat part (F) or the bend (D):  P = P(F,F,F)=0 + P(F,F,D) + ... (8 corners)
usage: GAMMA=8 DELTA=0.5 SIG=0.3 WORKERS=4 python3 klein_planes.py N nu T out.json"""
import sys, os, json, time, itertools, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'klein-flow', 'periodic-box', 'blowup'))
import tubes
from flow3d import Box, NS3

N, nu, T, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
box = Box(N, 1.0, Nz=N); ns = NS3(box, nu, 'T3')
U = tubes.initial(box, 'none', 0, 0.0, 0.25, 2)
KX, KY, KZ = box.KX, box.KY, box.KZ
flat = (KX == 0)
wz = np.where(KZ == 0, 1.0, 2.0)
def curl(V): return 1j*np.array([KY*V[2]-KZ*V[1], KZ*V[0]-KX*V[2], KX*V[1]-KY*V[0]])
def inv(V): return np.array([box.inv(c) for c in V])
def diag(U, t):
    W = curl(U); tot = (wz*(np.abs(W)**2).sum(0)).sum(); pl = (wz*(np.abs(W)**2).sum(0)*flat).sum()
    parts = {'F': U*flat, 'D': U*(~flat)}
    w = {p: inv(curl(parts[p])) for p in parts}
    gu = {p: np.array([[box.inv(1j*[KX, KY, KZ][j]*parts[p][i]) for j in range(3)] for i in range(3)]) for p in parts}
    P = {}
    for a, b, c in itertools.product('FD', repeat=3):
        st = np.einsum('jxyz,ijxyz->ixyz', w[b], gu[c])
        P[a+b+c] = float(np.mean((w[a]*st).sum(0)))
    return dict(t=t, Z=float(ns.enstrophy(U)), wmax=float(ns.wmax(U)), plane=float(pl/tot), P=P)
rows = [diag(U, 0.0)]; t = 0.0; nxt = 0.1; t0 = time.time()
while t < T - 1e-12:
    dt = min(ns.dt_cfl(U, 0.8), T - t, 0.05); U = ns.step(U, dt); t += dt
    if t >= nxt - 1e-9:
        rows.append(diag(U, t)); nxt += 0.1
        r = rows[-1]; print(f"t={t:.2f} Z={r['Z']:.2f} wmax={r['wmax']:.1f} plane={r['plane']:.3f} P={sum(r['P'].values()):.2f}", flush=True)
json.dump(dict(N=N, nu=nu, T=T, rows=rows, secs=time.time()-t0), open(out, 'w'))
