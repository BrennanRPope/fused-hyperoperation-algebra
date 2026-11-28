"""DRAFT. Blow-up hunt: evolve a candidate flow and watch for a blow-up.
Pseudo-spectral, 2/3 dealiasing, RK4 with an integrating factor for viscosity, 2*pi periodic box.
The push-back dial theta (Brennan's gravity idea): the pressure removes theta of the squeezing part of advection,
  du/dt = -P_theta[(u.grad)u] + nu Lap u,   P_theta = I - theta k k^T/|k|^2
theta = 1 is Navier-Stokes (Euler at nu = 0); theta = 0 is pressureless dust (Burgers; Zel'dovich's model of
gravitational collapse in cosmology), which folds in finite time.
Diagnostics every DIAG in time (default 0.01, and then also every 5 steps): energy, mean|grad u|^2, max|grad u|, max spin, max squeeze (-div u),
and delta, the width of the analyticity strip, fitted from E(k) ~ C k^-n exp(-2 delta k) (Sulem, Sulem, Frisch 1983):
a blow-up at time T is delta -> 0 at T. The run stops once delta*kmax < 1, or the spectrum is flat to the cutoff (no longer resolved).
'tail' records the share of energy in the top fifth of kept wavenumbers.
usage: python3 hunt.py ic N nu theta T out.json      ic in tg, rings, klein"""
import sys, os, json, time, numpy as np, scipy.fft as sf
W = int(os.environ.get('WORKERS', 4)); DIAG = float(os.environ.get('DIAG', 0.01))
ic, N, nu, theta, T, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
x = np.arange(N)*2*np.pi/N; X = np.meshgrid(x, x, x, indexing='ij')
k1 = np.fft.fftfreq(N, 1/N); kz = np.fft.rfftfreq(N, 1/N)
K = np.array(np.meshgrid(k1, k1, kz, indexing='ij')); K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
dea = np.all(np.abs(K) < N/3, axis=0)
wz = np.where(K[2] == 0, 1.0, 2.0); wz[..., -1] = 1.0 if N % 2 == 0 else 2.0
fwd = lambda a: sf.rfftn(a, axes=(-3, -2, -1), workers=W)
inv = lambda a: sf.irfftn(a, s=(N, N, N), axes=(-3, -2, -1), workers=W)
def proj(uh, th=1.0): return uh - th*K*(K*uh).sum(0)/K2s
def mean_sq(uh): return float((wz*np.abs(uh)**2).sum()/N**6)        # mean of the field squared (sum over components)

def initial():
    if ic == 'tg':                                                   # Taylor-Green: a single layered flat flow (u_z = 0)
        u = np.array([np.sin(X[0])*np.cos(X[1])*np.cos(X[2]), -np.cos(X[0])*np.sin(X[1])*np.cos(X[2]), 0*X[0]])
        return fwd(u)
    if ic == 'rings':                                                # round 4's fastest-growing 3D flow (colliding rings), 32^3, visc 0.05
        ug = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs/drain/rings_d3_n32_nu0.05_u.npz'))['u'].astype(float)
        M = ug.shape[1]; src = sf.fftn(ug, axes=(1, 2, 3)); full = np.zeros((3, N, N, N), complex); h = M//2
        idx = np.r_[0:h, -h:0]                                       # Fourier interpolation M^3 -> N^3 (kept modes are |k| <= M/4)
        full[np.ix_(range(3), idx % N, idx % N, idx % N)] = src[np.ix_(range(3), idx % M, idx % M, idx % M)]
        uh = proj(fwd(np.real(sf.ifftn(full, axes=(1, 2, 3)))*(N/M)**3))
        ref = mean_sq(1j*K[:, None]*initial_tg()[None])               # same mean |grad u|^2 as Taylor-Green
        return uh*np.sqrt(ref/mean_sq(1j*K[:, None]*uh[None]))
    if ic == 'klein':                                                # antiparallel vortex tubes with a bend (the Klein pair's shape)
        sig, a, d0 = 0.3, 0.25, 0.5
        zc = a*np.cos(2*X[1])                                        # bent tubes along y, bend toward each other
        r1 = (X[0] - np.pi + d0)**2 + (X[2] - np.pi - zc)**2; r2 = (X[0] - np.pi - d0)**2 + (X[2] - np.pi - zc)**2
        g = np.exp(-r1/sig**2) - np.exp(-r2/sig**2)
        w = np.array([0*g, g, 0*g]); wh = fwd(w)
        uh = proj(1j*np.cross(K, wh, axis=0)/K2s*1.0)
        ref = mean_sq(1j*K[:, None]*initial_tg()[None])
        return uh*np.sqrt(ref/mean_sq(1j*K[:, None]*uh[None]))
def initial_tg():
    return fwd(np.array([np.sin(X[0])*np.cos(X[1])*np.cos(X[2]), -np.cos(X[0])*np.sin(X[1])*np.cos(X[2]), 0*X[0]]))

def rhs(uh):
    u = inv(uh); w = inv(1j*np.cross(K, uh, axis=0))
    nl = fwd(np.cross(u, w, axis=0)) - 1j*K*fwd(0.5*(u**2).sum(0))[None]   # -(u.grad)u = u x w - grad |u|^2/2
    return proj(nl, theta)*dea

kk = np.arange(1, N//2); shell = np.rint(np.sqrt(K2)).astype(int)
def diag(uh, t):
    G = np.array([[inv(1j*K[j]*uh[i]) for j in range(3)] for i in range(3)])
    gmax = float(np.sqrt((G**2).sum((0, 1))).max())
    w = np.array([G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]]); wmax = float(np.sqrt((w**2).sum(0)).max())
    div = G[0, 0] + G[1, 1] + G[2, 2]; sq = float(max(0.0, -div.min()))
    E = np.bincount(shell.ravel(), weights=(0.5*wz*(np.abs(uh)**2).sum(0)).ravel(), minlength=N)[:N//2]/N**6
    km = int(N/3) - 1; sel = (kk >= 3) & (kk <= km) & (E[kk] > 1e-28)
    kf, Ef = kk[sel], E[kk][sel]
    if len(kf) > 4:
        Am = np.vstack([np.ones_like(kf, float), -np.log(kf), -2*kf]).T; c = np.linalg.lstsq(Am, np.log(Ef), rcond=None)[0]; delta, n = float(c[2]), float(c[1])
    else: delta, n = float('nan'), float('nan')
    tail = float(E[int(0.8*N/3):int(N/3)].sum()/max(E.sum(), 1e-300))
    return dict(t=t, tail=tail, E=0.5*mean_sq(uh), Z=mean_sq(1j*K[:, None]*uh[None]), gmax=gmax, wmax=wmax, squeeze=sq, delta=delta, n=n,
                Ek=E[1:N//2].tolist() if abs((t*10) - round(t*10)) < 1e-6 else None)

uh = initial()*dea; uh = proj(uh, 1.0)
rows = [diag(uh, 0.0)]; t = 0.0; step = 0; t0 = time.time(); kmax = N/3; nxt = DIAG; stop = ''
while t < T - 1e-12:
    r = rows[-1]
    umax = float(np.abs(inv(uh)).max()); dt = min(0.5*(2*np.pi/N)/max(umax, 1e-9), 0.25/max(r["gmax"], 1e-9), 0.02, T - t)
    E1 = np.exp(-nu*K2*dt/2)[None]; E2 = E1**2
    k1_ = rhs(uh); k2_ = rhs(E1*(uh + dt/2*k1_)); k3_ = rhs(E1*uh + dt/2*k2_); k4_ = rhs(E2*uh + dt*E1*k3_)
    uh = E2*uh + dt/6*(E2*k1_ + 2*E1*(k2_ + k3_) + k4_); t += dt; step += 1
    if t >= nxt - 1e-9 or (step % 5 == 0 and DIAG <= 0.01):
        rows.append(diag(uh, t)); nxt = t + DIAG
        r = rows[-1]
        if len(rows) % 10 == 0: print(f"t={t:.3f} E={r['E']:.4f} Z={r['Z']:.3f} gmax={r['gmax']:.2f} wmax={r['wmax']:.2f} squeeze={r['squeeze']:.2f} delta={r['delta']:.4f} n={r['n']:.2f}", flush=True)
        if not np.isfinite(r['E']) or r['E'] > 10*rows[0]['E']: stop = 'diverged'; break
        if r['delta'] * kmax < 1: stop = 'unresolved (delta*kmax < 1)'; break
        if r['n'] < 1 and r['delta']*kmax < 3: stop = 'unresolved (spectrum flat to the cutoff, n < 1)'; break
json.dump(dict(ic=ic, N=N, nu=nu, theta=theta, T=T, stop=stop, t_end=t, rows=rows, secs=time.time() - t0), open(out, 'w'))
print('done', ic, N, nu, theta, 't_end', round(t, 3), stop, f'{time.time()-t0:.0f}s')
