"""DRAFT. Faster blow-up hunt solver: same equations, diagnostics, stopping rules and output as hunt.py.
du/dt = -P_theta[(u.grad)u] + nu Lap u,  P_theta = I - theta k k^T/|k|^2  (theta = 1: Navier-Stokes / Euler; 0: pressureless dust).
Speed-ups over hunt.py (asked for by Brennan: "Please refactor the approach to improve performance"):
  - single precision (float32 / complex64): half the memory traffic, faster FFTs;
  - wavenumbers kept as 1D vectors that broadcast, never as full 3D arrays;
  - every array operation done per slab in a thread pool (numpy releases the lock on large arrays), with no np.cross temporaries;
  - at theta = 1 the |u|^2/2 gradient is skipped (the projection removes it exactly): 9 FFTs per stage instead of 10;
  - the integrating factor is skipped at nu = 0;
  - the spectrum fit for delta uses shells above 1e-10 of the peak (single-precision round-off floor; hunt.py used 1e-28),
    and the stopping rules apply only once those shells reach at least halfway to the cutoff; peak spin and the spectrum reuse the first RK stage's fields.
usage: [WORKERS=4] [DIAG=0.025] python3 hunt2.py ic N nu theta T out.json      ic in tg, rings, klein"""
import sys, os, json, time, numpy as np, scipy.fft as sf
from concurrent.futures import ThreadPoolExecutor
W = int(os.environ.get('WORKERS', 4)); DIAG = float(os.environ.get('DIAG', 0.025))
ic, N, nu, theta, T, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
f4, c8 = np.float32, np.complex64
Nh = N//2 + 1; pool = ThreadPoolExecutor(W)
k1 = np.fft.fftfreq(N, 1/N).astype(f4); kzv = np.fft.rfftfreq(N, 1/N).astype(f4)
kx, ky, kz = k1[:, None, None], k1[None, :, None], kzv[None, None, :]
K2 = (kx**2 + ky**2 + kz**2).astype(f4); iK2 = np.where(K2 == 0, 0, 1/np.where(K2 == 0, 1, K2)).astype(f4)
mx = (np.abs(k1) < N/3).astype(f4); mz = (np.abs(kzv) < N/3).astype(f4)
dea = mx[:, None, None]*mx[None, :, None]*mz[None, None, :]
wz = np.where(kzv == 0, 1.0, 2.0).astype(f4); wz[-1] = 1.0 if N % 2 == 0 else 2.0
SL = [slice(i*N//W, (i+1)*N//W) for i in range(W)]
def par(fn): list(pool.map(fn, SL))
fwd = lambda a: sf.rfftn(a, axes=(-3, -2, -1), workers=W)
inv = lambda a: sf.irfftn(a, s=(N, N, N), axes=(-3, -2, -1), workers=W)
def ksl(s): return kx[s], ky, kz
def mean_sq(ah):                                   # mean of the field squared, from Fourier coefficients (ah: (..., N, N, Nh))
    return float(sum((wz*(np.abs(a)**2)).sum(dtype=np.float64) for a in ah.reshape(-1, N, N, Nh))/N**6)

def proj_inplace(ah, th):
    def f(s):
        a, (qx, qy, qz) = ah[:, s], ksl(s)
        d = (qx*a[0] + qy*a[1] + qz*a[2])*iK2[s]*th
        a[0] -= qx*d; a[1] -= qy*d; a[2] -= qz*d
    par(f)
def initial():
    x = (np.arange(N)*2*np.pi/N).astype(f4); X, Y, Z = x[:, None, None], x[None, :, None], x[None, None, :]
    tg = np.zeros((3, N, N, N), f4); tg[0] = np.sin(X)*np.cos(Y)*np.cos(Z); tg[1] = -np.cos(X)*np.sin(Y)*np.cos(Z)
    tgh = fwd(tg).astype(c8); Zref = sum(mean_sq((q*1j*tgh[None])) for q in []) if False else None
    def gradZ(uh): return sum(mean_sq((1j*q)*uh) for q in (kx, ky, kz))
    if ic == 'tg': return tgh
    if ic == 'rings':
        ug = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs/drain/rings_d3_n32_nu0.05_u.npz'))['u'].astype(float)
        M = ug.shape[1]; src = sf.fftn(ug, axes=(1, 2, 3)); full = np.zeros((3, N, N, N), complex); h = M//2; idx = np.r_[0:h, -h:0]
        full[np.ix_(range(3), idx % N, idx % N, idx % N)] = src[np.ix_(range(3), idx % M, idx % M, idx % M)]
        uh = fwd(np.real(sf.ifftn(full, axes=(1, 2, 3)))*(N/M)**3).astype(c8)
    elif ic == 'klein':
        sig, a, d0 = 0.3, 0.25, 0.5; zc = a*np.cos(2*Y)
        g = np.exp(-((X - np.pi + d0)**2 + (Z - np.pi - zc)**2)/sig**2) - np.exp(-((X - np.pi - d0)**2 + (Z - np.pi - zc)**2)/sig**2)
        gh = fwd(g.astype(f4)).astype(c8)                                   # spin along y only: u = curl^-1 (0, g, 0)
        uh = np.zeros((3, N, N, Nh), c8); uh[0] = 1j*(-kz*gh)*iK2; uh[2] = 1j*(kx*gh)*iK2
    proj_inplace(uh, 1.0)
    return (uh*np.sqrt(gradZ(tgh)/gradZ(uh))).astype(c8)

def curl_h(uh):
    wh = np.empty_like(uh)
    def f(s):
        a, (qx, qy, qz) = uh[:, s], ksl(s)
        wh[0, s] = 1j*(qy*a[2] - qz*a[1]); wh[1, s] = 1j*(qz*a[0] - qx*a[2]); wh[2, s] = 1j*(qx*a[1] - qy*a[0])
    par(f); return wh
def rhs(uh, keep=None):
    u = inv(uh); w = inv(curl_h(uh))
    if keep is not None: keep['u'], keep['w'] = u, w
    c = np.empty_like(u); q = np.empty_like(u[0]) if theta != 1 else None
    def f(s):
        c[0, s] = u[1, s]*w[2, s] - u[2, s]*w[1, s]; c[1, s] = u[2, s]*w[0, s] - u[0, s]*w[2, s]; c[2, s] = u[0, s]*w[1, s] - u[1, s]*w[0, s]
        if q is not None: q[s] = 0.5*(u[0, s]**2 + u[1, s]**2 + u[2, s]**2)
    par(f)
    nl = fwd(c).astype(c8, copy=False)                                       # -(u.grad)u = u x w - grad |u|^2/2
    if q is not None:
        qh = fwd(q).astype(c8, copy=False)
        def g(s):
            qx, qy, qz = ksl(s); nl[0, s] -= 1j*qx*qh[s]; nl[1, s] -= 1j*qy*qh[s]; nl[2, s] -= 1j*qz*qh[s]
        par(g)
    proj_inplace(nl, theta)
    def m(s): nl[:, s] *= dea[s]
    par(m); return nl

kk = np.arange(1, N//2); shell = np.rint(np.sqrt(K2)).astype(np.int32).ravel()
def diag(uh, t, u=None, w=None):
    G = [[inv((1j*q)*uh[i]) for q in (kx, ky, kz)] for i in range(3)]
    gmax = float(np.sqrt(sum(G[i][j].astype(np.float64)**2 for i in range(3) for j in range(3))).max())
    wv = [G[2][1] - G[1][2], G[0][2] - G[2][0], G[1][0] - G[0][1]]; wmax = float(np.sqrt(sum(x.astype(np.float64)**2 for x in wv)).max())
    sq = float(max(0.0, -(G[0][0] + G[1][1] + G[2][2]).min())); del G, wv
    E = np.bincount(shell, weights=(0.5*wz*(np.abs(uh)**2).sum(0)).ravel().astype(np.float64), minlength=N)[:N//2]/N**6
    km = int(N/3) - 1; sel = (kk >= 3) & (kk <= km) & (E[kk] > 1e-10*E.max()); kf, Ef = kk[sel], E[kk][sel]
    nfit = int(len(kf))
    if len(kf) > 4:
        Am = np.vstack([np.ones_like(kf, float), -np.log(kf), -2*kf]).T; cf = np.linalg.lstsq(Am, np.log(Ef), rcond=None)[0]; delta, n = float(cf[2]), float(cf[1])
    else: delta, n = float('nan'), float('nan')
    tail = float(E[int(0.8*N/3):int(N/3)].sum()/max(E.sum(), 1e-300))
    return dict(t=t, nfit=nfit, tail=tail, E=0.5*mean_sq(uh), Z=sum(mean_sq((1j*q)*uh) for q in (kx, ky, kz)), gmax=gmax, wmax=wmax, squeeze=sq, delta=delta, n=n,
                Ek=E[1:N//2].tolist() if abs((t*10) - round(t*10)) < 1e-6 else None)

uh = initial(); def_m = None
def mdea(s): uh[:, s] *= dea[s]
par(mdea)
rows = [diag(uh, 0.0)]; t = 0.0; t0 = time.time(); kmax = N/3; nxt = DIAG; stop = ''
while t < T - 1e-12:
    r = rows[-1]; keep = {}
    k1_ = rhs(uh, keep); umax = float(np.abs(keep['u']).max())
    dt = min(0.5*(2*np.pi/N)/max(umax, 1e-9), 0.25/max(r['gmax'], 1e-9), 0.02, T - t); del keep
    if nu > 0:
        E1 = np.exp(-nu*K2*dt/2).astype(f4); E2 = E1*E1
        k2_ = rhs(E1*(uh + (dt/2)*k1_)); k3_ = rhs(E1*uh + (dt/2)*k2_); k4_ = rhs(E2*uh + dt*E1*k3_)
        uh = E2*uh + (dt/6)*(E2*k1_ + 2*E1*(k2_ + k3_) + k4_)
    else:
        k2_ = rhs(uh + (dt/2)*k1_); k3_ = rhs(uh + (dt/2)*k2_); k4_ = rhs(uh + dt*k3_)
        def upd(s): uh[:, s] += (dt/6)*(k1_[:, s] + 2*(k2_[:, s] + k3_[:, s]) + k4_[:, s])
        par(upd)
    del k1_, k2_, k3_, k4_
    t += dt
    if t >= nxt - 1e-9:
        rows.append(diag(uh, t)); nxt = t + DIAG; r = rows[-1]
        if len(rows) % 10 == 0: print(f"t={t:.3f} E={r['E']:.4f} Z={r['Z']:.3f} gmax={r['gmax']:.2f} wmax={r['wmax']:.2f} squeeze={r['squeeze']:.2f} delta={r['delta']:.4f} n={r['n']:.2f} {time.time()-t0:.0f}s", flush=True)
        if not np.isfinite(r['E']) or r['E'] > 10*rows[0]['E']: stop = 'diverged'; break
        if r['nfit'] < 0.5*(kmax - 3): continue                       # spectrum not yet filled toward the cutoff: fit not meaningful
        if r['delta']*kmax < 1: stop = 'unresolved (delta*kmax < 1)'; break
        if r['n'] < 1 and r['delta']*kmax < 3: stop = 'unresolved (spectrum flat to the cutoff, n < 1)'; break
json.dump(dict(ic=ic, N=N, nu=nu, theta=theta, T=T, stop=stop, t_end=t, rows=rows, secs=time.time() - t0, solver='hunt2', diag_every=DIAG),
          open(out, 'w'))
print('done', ic, N, nu, theta, 't_end', round(t, 3), stop, f'{time.time()-t0:.0f}s')
