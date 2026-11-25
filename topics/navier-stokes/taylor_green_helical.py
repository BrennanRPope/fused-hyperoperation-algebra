"""3D Taylor-Green vortex, pseudo-spectral NS at N^3, with each velocity mode
split into right- and left-handed helical waves (curl eigenmodes, Waleffe).
Tracks max vorticity and splits enstrophy production  P = <w . (w.grad) u>
into its 8 handedness corners (sign of w, sign of w, sign of u).
Run: python3 taylor_green_helical.py  -> prints a table, saves tg_helical.npz"""
import numpy as np, itertools, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 64
nu = 1/400; dt = 0.01; T = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
k1 = np.fft.fftfreq(N, 1/N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
K = np.array([KX, KY, KZ]); K2 = (K**2).sum(0); K2[0,0,0] = 1; KA = np.sqrt(K2)
dealias = (np.abs(KX) < N/3) & (np.abs(KY) < N/3) & (np.abs(KZ) < N/3)
xs = np.arange(N)*2*np.pi/N
X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
u = np.array([np.sin(X)*np.cos(Y)*np.cos(Z), -np.cos(X)*np.sin(Y)*np.cos(Z), 0*X])
uh = np.fft.fftn(u, axes=(1,2,3))
fft = lambda f: np.fft.fftn(f, axes=(1,2,3)); ifft = lambda f: np.real(np.fft.ifftn(f, axes=(1,2,3)))
def proj(fh): return fh - K*(K*fh).sum(0)/K2
def curlh(fh): return 1j*np.cross(K, fh, axis=0)
def rhs(uh):
    u = ifft(uh); w = ifft(curlh(uh))
    nl = fft(np.cross(u, w, axis=0))*dealias
    return proj(nl) - nu*K2*uh
def helical(uh):
    c = curlh(uh)/KA
    return 0.5*(uh + c), 0.5*(uh - c)   # right-handed (+), left-handed (-)
def corners(uh):
    up, um = helical(uh); us = {'+': up, '-': um}
    w = {s: ifft(curlh(us[s])) for s in us}; uu = {s: us[s] for s in us}
    out = {}
    for a, b, c in itertools.product('+-', repeat=3):
        grad_u = np.array([[np.real(np.fft.ifftn(1j*K[j]*uu[c][i])) for j in range(3)] for i in range(3)])
        stretch = np.einsum('jxyz,ijxyz->ixyz', w[b], grad_u)
        out[a+b+c] = np.mean((w[a]*stretch).sum(0))
    return out
times, wmax, rows = [], [], []
t = 0.0; step = 0
while t <= T + 1e-9:
    if step % 50 == 0:
        w = ifft(curlh(uh)); times.append(t); wmax.append(np.sqrt((w**2).sum(0)).max())
        c = corners(uh); rows.append(c)
        same = c['+++'] + c['---']; tot = sum(c.values())
        print(f"t={t:4.1f} max|w|={wmax[-1]:7.3f}  P_total={tot:8.4f}  single-handed={same:8.4f}  mixed={tot-same:8.4f}", flush=True)
    k1_ = rhs(uh); k2_ = rhs(uh+0.5*dt*k1_); k3_ = rhs(uh+0.5*dt*k2_); k4_ = rhs(uh+dt*k3_)
    uh = uh + dt/6*(k1_+2*k2_+2*k3_+k4_); t += dt; step += 1
keys = sorted(rows[0])
np.savez('tg_helical.npz', times=times, wmax=wmax, keys=keys, P=np.array([[r[k] for k in keys] for r in rows]))
