"""DRAFT. Round 2: three perpendicular planes. (Solver and flow builders copied from twoplanes.py.)
Modes: two-tw (A + gA), two-un (A + gC), three-tw (A + gA + D, D a flat xz flow with the quarter-turn symmetry),
three-un (A + gC + D2, D2 a flat xz flow with no symmetry). Total starting energy 0.5 in every mode.
Diagnostics every 0.2: stretching split into 27 terms T[i,j,m] = <w_i w_j d_j u^(m)_i>, u^(m) = Biot-Savart of w_m e_m
(plane m = coordinate plane with normal e_m), plus each plane's own spin.
Usage: python3 threeplanes.py N nu T mode seed out.json

Original header follows.
P1 = xy-plane flow A(x,y); P2 = yz-plane flow B(y,z). Shared line = y-axis.
  twisted   : B = g A, g = quarter turn about the y-axis + shift pi/2 along y (the quarter-turn flat
              3-manifold's twist). The field A + gA is exactly invariant (A carries the half-turn symmetry),
              so the run stays inside that twisted space for all time.
  untwisted : B = g C with C an independent flow of the same statistics (no symmetry).
k dial on the crossing (a o_k b = a + b + k ab applied to flows):
  N_k(u) = N(uA) + N(uB) + k [N(u) - N(uA) - N(uB)],  uA = modes with kz=0, uB = modes with kx=0.
  k=0: each plane only talks to itself (and through the shared line's modes); k=1: true Navier-Stokes.
  N_k conserves energy for every k (each piece does).
Pseudo-spectral, rotational form, Leray projection, 2/3 dealiasing, integrating-factor RK4.
Usage: python3 twoplanes.py N nu T kdial mode seed out.json   (mode: twisted|untwisted|A|B)"""
import sys, json, time, numpy as np
from scipy import fft as sfft

N = int(sys.argv[1]); nu = float(sys.argv[2]); T = float(sys.argv[3]); kd = 1.0
mode = sys.argv[4]; seed = int(sys.argv[5]); out = sys.argv[6]
W = 4
k1 = sfft.fftfreq(N, 1/N); kr = sfft.rfftfreq(N, 1/N)
KX, KY, KZ = np.meshgrid(k1, k1, kr, indexing='ij'); K = np.array([KX, KY, KZ])
K2 = (K**2).sum(0); K2i = np.where(K2 > 0, 1/np.where(K2 > 0, K2, 1), 0)
mask = ((np.abs(KX) < N/3) & (np.abs(KY) < N/3) & (np.abs(KZ) < N/3)).astype(float)
mA = (KZ == 0); mB = (KX == 0); axis = mA & mB
# rfft halves z: modes with kz=0 appear once; weights for Parseval sums
wgt = np.where(KZ == 0, 1.0, 2.0); wgt[..., -1] = np.where(N % 2 == 0, 1.0, 2.0)
fwd = lambda f: sfft.rfftn(f, axes=(1, 2, 3), workers=W)
inv = lambda F: sfft.irfftn(F, s=(N, N, N), axes=(1, 2, 3), workers=W)
def leray(U): return U - K*(K*U).sum(0)*K2i
def curl(U): return 1j*np.array([KY*U[2]-KZ*U[1], KZ*U[0]-KX*U[2], KX*U[1]-KY*U[0]])
def Nl(U):
    u = inv(U); w = inv(curl(U))
    return leray(fwd(np.cross(u, w, axis=0)))*mask        # = -P[(u.grad)u]
def rhs_nl(U):
    if kd == 1.0: return Nl(U)
    UA = U*mA; UB = U*mB
    n0 = Nl(UA) + Nl(UB)
    return n0 if kd == 0.0 else n0 + kd*(Nl(U) - n0)

xs = np.arange(N)*2*np.pi/N
def flat2d(rng, k0=3.0, Ng=64):
    """random 2D stream function psi(a,b) on [0,2pi)^2, energy spectrum peaked at k0. Drawn on a fixed
    64x64 grid (so every resolution starts from the same flow) and spectrally padded to N."""
    a = sfft.fftfreq(Ng, 1/Ng); A1, B1 = np.meshgrid(a, a, indexing='ij'); kk = np.sqrt(A1**2+B1**2)
    amp = np.where(kk > 0, kk**1.5*np.exp(-(kk/k0)**2)/np.where(kk > 0, kk, 1)**2, 0)  # psi amplitude
    ph = np.exp(2j*np.pi*rng.random((Ng, Ng)))
    P = amp*ph*(np.abs(A1) < Ng/3)*(np.abs(B1) < Ng/3)
    if N != Ng:
        Q = np.zeros((N, N), complex); h = Ng//2
        idx = np.r_[0:h, N-h:N]; Q[np.ix_(idx, idx)] = P[np.ix_(np.r_[0:h, Ng-h:Ng], np.r_[0:h, Ng-h:Ng])]*(N/Ng)**2
        P = Q
    return np.real(sfft.ifftn(P))
def half_turn_sym(psi):
    """psi(a,b) -> make A from psi carry the half-turn g^2: A(x,y) -> (-Ax(-x,y-pi), Ay(-x,y-pi)).
    In stream-function terms psi(x,y) -> -psi(-x, y-pi) keeps the velocity map; average the two."""
    flipped = -np.roll(psi[(-np.arange(N)) % N, :], N//2, axis=1)
    return 0.5*(psi + flipped)
def A_field(psi):  # P1 flow from psi(x,y): u = (d_y psi, -d_x psi, 0), constant in z
    a = sfft.fftfreq(N, 1/N); A1, B1 = np.meshgrid(a, a, indexing='ij'); P = sfft.fftn(psi)
    ux = np.real(sfft.ifftn(1j*B1*P)); uy = np.real(sfft.ifftn(-1j*A1*P))
    u = np.zeros((3, N, N, N)); u[0] = ux[:, :, None]; u[1] = uy[:, :, None]; return u
def g_image(u):
    """(g u)(x,y,z) = R u(z, y - pi/2, -x), R(x,y,z) = (-z, y, x): quarter turn about y plus a quarter shift."""
    idx = np.arange(N); Ii, Jj, Ll = np.meshgrid(idx, idx, idx, indexing='ij')
    s = u[:, Ll, (Jj - N//4) % N, (-Ii) % N]
    return np.array([-s[2], s[1], s[0]])

import itertools, os
def xz_flow(psi):   # plane y (xz-plane) flow from psi(x,z): u = (d_z psi, 0, -d_x psi), constant in y
    a = sfft.fftfreq(N, 1/N); A1, B1 = np.meshgrid(a, a, indexing='ij'); P = sfft.fftn(psi)
    ux = np.real(sfft.ifftn(1j*B1*P)); uz = np.real(sfft.ifftn(-1j*A1*P))
    u = np.zeros((3, N, N, N)); u[0] = ux[:, None, :]; u[2] = uz[:, None, :]; return u
def g_orbit_avg(u):
    v = [u]
    for _ in range(3): v.append(g_image(v[-1]))
    return sum(v)/4
rng = np.random.default_rng(seed)
psiA = half_turn_sym(flat2d(rng)); uA = A_field(psiA)
if mode.endswith('tw'): uB = g_image(uA)
else: uB = g_image(A_field(half_turn_sym(flat2d(np.random.default_rng(seed + 1000)))))
u0 = uA + uB
if mode.startswith('three'):
    D = xz_flow(flat2d(np.random.default_rng(seed + 2000)))
    if mode.endswith('tw'): D = g_orbit_avg(D)
    D *= np.sqrt(((uA**2).sum(0).mean())/max((D**2).sum(0).mean(), 1e-30))   # same energy as plane A's flow
    u0 = u0 + D
u0 *= np.sqrt(1.0/(u0**2).sum(0).mean())     # total energy 0.5
if mode.endswith('tw'): print('twist invariance error', np.abs(g_image(u0) - u0).max()/np.abs(u0).max())
U = leray(fwd(u0))*mask
Eye = np.eye(3)
def diag(U, t):
    Wh = curl(U); w = inv(Wh)
    Z = 0.5*(w**2).sum(0).mean()
    up = [curl(np.array([Eye[m][i]*Wh[m]*K2i for i in range(3)])) for m in range(3)]
    Tt = np.zeros((3, 3, 3))
    for m in range(3):
        for i in range(3):
            if i == m: continue
            gi = inv(np.array([1j*K[j]*up[m][i] for j in range(3)]))
            for j in range(3): Tt[i, j, m] = np.mean(w[i]*w[j]*gi[j])
    return dict(t=t, Z=float(Z), wmax=float(np.sqrt((w**2).sum(0)).max()),
                comp_Z=[float(0.5*(w[i]**2).mean()) for i in range(3)], comp_max=[float(np.abs(w[i]).max()) for i in range(3)],
                T=Tt.tolist())
rows = [diag(U, 0.0)]; t = 0.0; nxt = 0.2; t0 = time.time()
Ef = lambda h: np.exp(-nu*K2*h)
while t < T - 1e-12:
    umax = np.abs(inv(U)).max(); dt = min(0.4*(2*np.pi/N)/umax, 0.05, T - t)
    e1, e2 = Ef(dt/2), Ef(dt)
    a = rhs_nl(U); b = rhs_nl(e1*(U + dt/2*a)); c = rhs_nl(e1*U + dt/2*b); d = rhs_nl(e2*U + dt*e1*c)
    U = e2*U + dt/6*(e2*a + 2*e1*(b + c) + d); t += dt
    if t >= nxt - 1e-9 or t >= T - 1e-12:
        rows.append(diag(U, t)); nxt += 0.2
json.dump(dict(N=N, nu=nu, T=T, mode=mode, seed=seed, rows=rows, secs=time.time()-t0), open(out, 'w'))
r0 = rows[0]; zm = max(r['Z'] for r in rows)
print(f"{mode} s{seed}: Z0={r0['Z']:.3f} Zpeak/Z0={zm/r0['Z']:.3f} wmax_peak={max(r['wmax'] for r in rows):.1f} {time.time()-t0:.0f}s")
