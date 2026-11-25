"""Full Navier-Stokes vs the single-handed (homochiral) version of Biferale et al.
Same right-handed random start; in the single-handed run the nonlinear term is
projected back onto right-handed waves every step. Tracks enstrophy and max|w|.
Run: python3 homochiral.py"""
import numpy as np
N = 48; nu = 1/300; dt = 0.005; T = 6.0
k1 = np.fft.fftfreq(N, 1/N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
K = np.array([KX, KY, KZ]); K2 = (K**2).sum(0); K2[0,0,0] = 1; KA = np.sqrt(K2)
dealias = (np.abs(KX) < N/3) & (np.abs(KY) < N/3) & (np.abs(KZ) < N/3)
F = lambda f: np.fft.fftn(f, axes=(1,2,3)); Fi = lambda f: np.real(np.fft.ifftn(f, axes=(1,2,3)))
def proj(fh): return fh - K*(K*fh).sum(0)/K2
def curlh(fh): return 1j*np.cross(K, fh, axis=0)
def right(fh): return 0.5*(fh + curlh(fh)/KA)
rng = np.random.default_rng(1)
spec = KA**4*np.exp(-2*(KA/3)**2)
uh = proj(F(rng.standard_normal((3, N, N, N))))
uh *= np.sqrt(spec)/(np.sqrt((np.abs(uh)**2).sum(0))+1e-30)
uh = F(Fi(right(uh)))       # keep only right-handed waves, real field
uh = right(uh); uh = F(Fi(uh))
E0 = np.mean((Fi(uh)**2).sum(0)); uh *= np.sqrt(0.5/E0)   # energy 1/4 per unit
def rhs(uh, single):
    u = Fi(uh); w = Fi(curlh(uh))
    nl = proj(F(np.cross(u, w, axis=0))*dealias)
    if single: nl = right(nl)
    return nl - nu*K2*uh
def stats(uh):
    w = Fi(curlh(uh)); u = Fi(uh)
    hl = np.mean((u*w).sum(0)); ens = np.mean((w**2).sum(0))
    return ens, np.sqrt((w**2).sum(0)).max(), hl
for single in (False, True):
    v = uh.copy(); t = 0; step = 0; out = []
    while t <= T+1e-9:
        if step % 100 == 0:
            ens, wm, hl = stats(v); out.append((t, ens, wm, hl))
        k1_ = rhs(v, single); k2_ = rhs(v+0.5*dt*k1_, single); k3_ = rhs(v+0.5*dt*k2_, single); k4_ = rhs(v+dt*k3_, single)
        v = v + dt/6*(k1_+2*k2_+2*k3_+k4_); t += dt; step += 1
    name = 'single-handed' if single else 'full NS'
    for r in out: print(f"{name:13s} t={r[0]:4.1f} enstrophy={r[1]:8.3f} max|w|={r[2]:7.3f} helicity={r[3]:+.3f}", flush=True)
    np.save(f"homochiral_{'single' if single else 'full'}.npy", np.array(out))
