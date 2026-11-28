"""DRAFT (round 3). How fast can two planes feed each other's spin?

Instantaneous growth of total spin (enstrophy Z = mean |w|^2) under Navier-Stokes:
    dZ/dt = 2 R(u),   R(u) = P(u) - nu * mean |grad w|^2,   P(u) = mean( w_i w_j d_j u_i )   (stretching)
We maximize R over flows with fixed Z (the Lu-Doering problem, J. Math. Phys. / Indiana 2008) in two classes:
    'any'  : any incompressible flow in the 2*pi box
    'two'  : flows made of two perpendicular planes only (spin along x and z, none along y):
             w_y = 0  <=>  per Fourier mode u is along  e_y|k|^2 - k_y k   (k not along y)
    'twod' : the same, but the two planes share the box diagonal: no spin along (1,1,1)/sqrt(3)
    'one'  : one layered flat flow only (no velocity along x, may vary in x). Its spin is not purely along x,
             so it can stretch; a strict flat flow (also constant in x) cannot.
Modes are kept with |k_i| <= N/4, so every product below is computed exactly on the N-grid (no aliasing).
Ascent: nonlinear conjugate gradient in the enstrophy metric (gradient divided by |k|^2), on the sphere Z = Z0.

usage: python3 feed.py N nu cls seed out.json [Z0=1] [iters=600]
"""
import sys, json, time, numpy as np, scipy.fft as sf
W = 4
def setup(N):
    k1 = np.fft.fftfreq(N, 1/N)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
    K = np.array([KX, KY, KZ]); K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
    keep = (np.abs(KX) <= N//4) & (np.abs(KY) <= N//4) & (np.abs(KZ) <= N//4) & (K2 > 0)
    return K, K2, K2s, keep
def fft(a): return sf.fftn(a, axes=(-3, -2, -1), workers=W)
def ifft(a): return sf.ifftn(a, axes=(-3, -2, -1), workers=W).real
class Prob:
    def __init__(s, N, nu, cls):
        s.N, s.nu, s.cls = N, nu, cls
        s.K, s.K2, s.K2s, s.keep = setup(N)
        K = s.K
        if cls == 'two':
            ey = np.zeros_like(K); ey[1] = 1
            v = ey*s.K2 - K[1]*K                              # direction allowed when k not along y
            vn = np.sqrt((v**2).sum(0)); along_y = (K[0] == 0) & (K[2] == 0)
            s.v = np.where(vn > 0, v/np.where(vn > 0, vn, 1), 0); s.along_y = along_y
        if cls == 'twod':                                     # two planes sharing the box diagonal: no spin along n = (1,1,1)/sqrt3
            n = np.array([1.0, 1.0, 1.0])/np.sqrt(3); nk = np.cross(n[:, None, None, None]*np.ones_like(K), K, axis=0)
            v = np.cross(K, nk, axis=0); vn = np.sqrt((v**2).sum(0))
            s.v = np.where(vn > 0, v/np.where(vn > 0, vn, 1), 0); s.along_y = np.sqrt((nk**2).sum(0)) < 1e-9
        if cls == 'one':
            ex = np.zeros_like(K); ex[0] = 1                  # u_x = 0 and div u = 0: u along k x e_x
            v = np.cross(K, ex, axis=0); vn = np.sqrt((v**2).sum(0))
            s.v = np.where(vn > 0, v/np.where(vn > 0, vn, 1), 0)
    def proj(s, uh):
        """project Fourier velocity onto the class (and onto kept modes)"""
        K = s.K
        uh = uh - K*(K*uh).sum(0)/s.K2s                       # divergence free
        if s.cls in ('two', 'twod'):
            pv = s.v*(s.v*uh).sum(0)
            uh = np.where(s.along_y[None], uh, pv)
        elif s.cls == 'one':
            uh = s.v*(s.v*uh).sum(0)
        return uh*s.keep
    def fields(s, uh):
        K = s.K; ik = 1j*K
        w_h = np.array([ik[1]*uh[2]-ik[2]*uh[1], ik[2]*uh[0]-ik[0]*uh[2], ik[0]*uh[1]-ik[1]*uh[0]])
        du = np.array([[ifft(ik[j]*uh[i]) for j in range(3)] for i in range(3)])   # du[i,j] = d_j u_i
        return w_h, ifft(w_h), du
    def value(s, uh):
        w_h, w, du = s.fields(uh)
        P = np.einsum('iabc,jabc,ijabc->', w, w, du)/s.N**3
        D = (s.K2[None]*np.abs(w_h)**2).sum()/s.N**6
        Z = (np.abs(w_h)**2).sum()/s.N**6
        return P - s.nu*D, P, D, Z, w_h, w, du
    def grad(s, uh, w, du):
        """L2 gradient of R (as Fourier array), mean inner product"""
        K = s.K; ik = 1j*K
        a = np.einsum('jabc,ijabc->iabc', w, du) + np.einsum('iabc,ijabc->jabc', w, du)   # (w.grad)u + (grad u)^T w
        a_h = fft(a)
        curl_a = np.array([ik[1]*a_h[2]-ik[2]*a_h[1], ik[2]*a_h[0]-ik[0]*a_h[2], ik[0]*a_h[1]-ik[1]*a_h[0]])
        wgw = np.array([sum(w[j]*ifft(ik[j]*fft(w[i])) for j in range(3)) for i in range(3)])
        g = curl_a - fft(wgw) - s.nu*2*(s.K2**2)[None]*uh
        return s.proj(g)
def Zof(p, uh): return (p.K2[None]*np.abs(uh)**2).sum()/p.N**6
def normalize(p, uh, Z0): return uh*np.sqrt(Z0/Zof(p, uh))
def init(p, seed, Z0):
    rng = np.random.default_rng(seed)
    uh = fft(rng.standard_normal((3, p.N, p.N, p.N)))
    uh = uh*np.exp(-p.K2/8.0)[None]                         # start at large scales
    return normalize(p, p.proj(uh), Z0)
def ascend(p, uh, Z0, iters=600, tol=1e-9):
    """Riemannian nonlinear conjugate gradient (Polak-Ribiere+) on the sphere Z = Z0, enstrophy metric."""
    ip = lambda a, b: (p.K2[None]*(a.conj()*b).real).sum()/p.N**6
    def rgrad(uh, w, du):
        g = p.proj(p.grad(uh, w, du)/p.K2s[None])
        return g - ip(uh, g)/ip(uh, uh)*uh
    R, P, D, Z, w_h, w, du = p.value(uh); g = rgrad(uh, w, du); d = g; step = 0.05; hist = []
    for it in range(iters):
        if ip(d, g) <= 0: d = g                              # restart if not an ascent direction
        dn = np.sqrt(ip(d, d)/Z0)
        if np.sqrt(ip(g, g)/Z0) < tol: break
        while True:
            trial = normalize(p, p.proj(uh + step*d/dn*np.sqrt(Z0)), Z0)
            Rt = p.value(trial)
            if Rt[0] > R: break
            step *= 0.5
            if step < 1e-10: return uh, R, P, D, hist
        uh = trial; R, P, D, Z, w_h, w, du = Rt
        gn = rgrad(uh, w, du)
        d = d - ip(uh, d)/ip(uh, uh)*uh                       # transport old direction to new tangent space
        beta = max(0.0, ip(gn, gn - g)/ip(g, g))
        d = gn + beta*d; g = gn; step = min(step*1.6, 1.0)
        hist.append((it, R, P, D))
        if len(hist) > 40 and abs(hist[-1][1]-hist[-40][1]) < 1e-6*abs(hist[-1][1]): break
    return uh, R, P, D, hist
def gradcheck(p, Z0=1.0, seed=0):
    uh = init(p, seed, Z0); R0 = p.value(uh); g = p.grad(uh, R0[5], R0[6])
    rng = np.random.default_rng(seed+99); dh = p.proj(fft(rng.standard_normal(uh.shape))); dh = dh/np.sqrt(Zof(p, dh))*1e-6
    num = (p.value(uh+dh)[0] - p.value(uh-dh)[0])/2
    ana = (g.conj()*dh).real.sum()/p.N**6
    return num, ana
def diag(p, uh):
    """how the stretching at the optimum splits: by plane pair, and spin shares"""
    K = p.K; ik = 1j*K
    w_h, w, du = p.fields(uh)
    comp = [(np.abs(w_h[i])**2).sum()/p.N**6 for i in range(3)]
    M = np.einsum('iabc,jabc->ij', w, w)/p.N**3                 # spin covariance: its smallest eigenvalue says how far from two planes in its best frame
    um = []
    for m in range(3):                                       # Biot-Savart of w_m e_m
        wm = np.zeros_like(w_h); wm[m] = w_h[m]
        um.append(np.array([ik[1]*wm[2]-ik[2]*wm[1], ik[2]*wm[0]-ik[0]*wm[2], ik[0]*wm[1]-ik[1]*wm[0]])/p.K2s[None])
    T = np.zeros((3, 3, 3))
    for m in range(3):
        dum = np.array([[ifft(ik[j]*um[m][i]) for j in range(3)] for i in range(3)])
        for i in range(3):
            for j in range(3): T[i, j, m] = (w[i]*w[j]*dum[i, j]).mean()
    kmean = np.sqrt((p.K2[None]**2*np.abs(uh)**2).sum()/(p.K2[None]*np.abs(uh)**2).sum())
    spec = np.abs(uh)**2*p.K2[None]
    kmax = int(np.sqrt(p.K2).max()); Kn = np.sqrt(p.K2).round().astype(int)
    shell = np.bincount(Kn.ravel(), spec.sum(0).ravel(), minlength=kmax+1)
    return dict(comp_Z=[float(c) for c in comp], M=M.tolist(), M_eig=sorted(np.linalg.eigvalsh(M).tolist()), T=T.tolist(), k_rms=float(kmean), shell=shell.tolist())
if __name__ == '__main__':
    N, nu, cls, seed, out = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
    Z0 = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    iters = int(sys.argv[7]) if len(sys.argv) > 7 else 600
    p = Prob(N, nu, cls); t0 = time.time()
    uh = init(p, seed, Z0); uh, R, P, D, hist = ascend(p, uh, Z0, iters)
    d = diag(p, uh)
    res = dict(N=N, nu=nu, cls=cls, seed=seed, Z0=Z0, R=float(R), P=float(P), D=float(D), iters=len(hist),
               hist=[(int(a), float(b), float(c), float(e)) for a, b, c, e in hist[::5]], secs=time.time()-t0, **d)
    json.dump(res, open(out, 'w'))
    import os
    if os.environ.get('SAVE'): np.savez_compressed(out.replace('.json', '_field.npz'), w=p.fields(uh)[1].astype(np.float32))
    print(f'{cls} N={N} nu={nu} seed={seed}: R={R:.5g} P={P:.5g} D={D:.5g} iters={len(hist)} k_rms={d["k_rms"]:.2f} {time.time()-t0:.0f}s')
