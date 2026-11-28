"""DRAFT (round 4). Fastest instantaneous growth of total spin in any number of dimensions d (3 or 4).

Spin is measured as Z = mean |grad u|^2 (in 3D this equals mean |w|^2 for incompressible flow; in 4D it is the
natural spin measure, half the squared spin tensor). Under Navier-Stokes, in any dimension,
    (1/2) dZ/dt = R(u) = -mean( d_j u_i  d_j u_k  d_k u_i ) - nu * mean |Lap u|^2 .
We maximize R at fixed Z (the Lu-Doering problem) on the d-torus of side 2*pi. Modes are kept with |k_i| <= N/4 so the
cubic term is computed exactly. Gradient of the cubic term: div_b C_ab with C = A^2 + A^T A + A A^T, A_ij = d_j u_i.
Ascent: nonlinear conjugate gradient on the sphere Z = Z0 in the metric of Z itself.

usage: python3 feedd.py d N nu seed out.json [Z0=1] [iters=500]
"""
import sys, json, time, numpy as np, scipy.fft as sf
W = 4
class Prob:
    def __init__(s, d, N, nu):
        s.d, s.N, s.nu = d, N, nu
        k1 = np.fft.fftfreq(N, 1/N)
        s.K = np.array(np.meshgrid(*([k1]*d), indexing='ij'))
        s.K2 = (s.K**2).sum(0); s.K2s = np.where(s.K2 == 0, 1, s.K2)
        s.keep = np.all(np.abs(s.K) <= N//4, axis=0) & (s.K2 > 0)
        s.ax = tuple(range(-d, 0))
    def fft(s, a): return sf.fftn(a, axes=s.ax, workers=W)
    def ifft(s, a): return sf.ifftn(a, axes=s.ax, workers=W).real
    def proj(s, uh):
        return (uh - s.K*(s.K*uh).sum(0)/s.K2s)*s.keep
    def A(s, uh):
        return np.array([[s.ifft(1j*s.K[j]*uh[i]) for j in range(s.d)] for i in range(s.d)])   # A[i,j] = d_j u_i
    def value(s, uh):
        A = s.A(uh)
        C3 = np.einsum('jiX,jkX,kiX->', *(A.reshape(s.d, s.d, -1),)*3)/A[0, 0].size    # sum d_j u_i d_j u_k d_k u_i
        D = (s.K2[None]**2*np.abs(uh)**2).sum()/s.N**(2*s.d)
        Z = (s.K2[None]*np.abs(uh)**2).sum()/s.N**(2*s.d)
        return -C3 - s.nu*D, -C3, D, Z, A
    def grad(s, uh, A):
        d = s.d; Af = A.reshape(d, d, -1)
        # F = -mean(A_ij A_kj A_ik) with A_ij = d_j u_i ;  dF = mean(dA_ab * (-B_ab)) ; B = A A + A^T A + A A^T (index form below)
        B = np.einsum('ikX,kjX->ijX', Af, Af) + np.einsum('ijX,ikX->kjX', Af, Af) + np.einsum('ijX,kjX->ikX', Af, Af)
        B = B.reshape(A.shape)
        # dF = -mean(d_b du_a B_ab) = mean(du_a d_b B_ab)
        g = np.array([sum(1j*s.K[b]*s.fft(B[a, b]) for b in range(d)) for a in range(d)])
        g = g - s.nu*2*(s.K2**2)[None]*uh
        return s.proj(g)
def Zof(p, uh): return (p.K2[None]*np.abs(uh)**2).sum()/p.N**(2*p.d)
def normalize(p, uh, Z0): return uh*np.sqrt(Z0/Zof(p, uh))
def init(p, seed, Z0):
    rng = np.random.default_rng(seed)
    uh = p.fft(rng.standard_normal((p.d,) + (p.N,)*p.d))*np.exp(-p.K2/8.0)[None]
    return normalize(p, p.proj(uh), Z0)
def ascend(p, uh, Z0, iters=500):
    ip = lambda a, b: (p.K2[None]*(a.conj()*b).real).sum()/p.N**(2*p.d)
    def rgrad(uh, A):
        g = p.proj(p.grad(uh, A)/p.K2s[None]); return g - ip(uh, g)/ip(uh, uh)*uh
    R, P, D, Z, A = p.value(uh); g = rgrad(uh, A); dd = g; step = 0.05; hist = []
    for it in range(iters):
        if ip(dd, g) <= 0: dd = g
        dn = np.sqrt(ip(dd, dd)/Z0)
        while True:
            trial = normalize(p, p.proj(uh + step*dd/dn*np.sqrt(Z0)), Z0); Rt = p.value(trial)
            if Rt[0] > R: break
            step *= 0.5
            if step < 1e-10: return uh, R, P, D, hist
        uh = trial; R, P, D, Z, A = Rt; gn = rgrad(uh, A)
        dd = dd - ip(uh, dd)/ip(uh, uh)*uh
        beta = max(0.0, ip(gn, gn - g)/ip(g, g)); dd = gn + beta*dd; g = gn; step = min(step*1.6, 1.0)
        hist.append((it, R, P, D))
        if len(hist) > 40 and abs(hist[-1][1] - hist[-40][1]) < 1e-6*abs(hist[-1][1]): break
    return uh, R, P, D, hist
def gradcheck(p, seed=0):
    uh = init(p, seed, 1.0); v = p.value(uh); g = p.grad(uh, v[4])
    rng = np.random.default_rng(seed + 7); dh = p.proj(p.fft(rng.standard_normal(uh.shape))); dh = dh/np.sqrt(Zof(p, dh))*1e-6
    num = (p.value(uh + dh)[0] - p.value(uh - dh)[0])/2; ana = (g.conj()*dh).real.sum()/p.N**(2*p.d)
    return num, ana
def diag(p, uh):
    """how many dimensions the winner really uses: share of Z in modes that do not vary along each axis"""
    e = p.K2[None]*np.abs(uh)**2; tot = e.sum()
    flat = [float(e[:, p.K[a] == 0].sum()/tot) for a in range(p.d)]
    kr = float(np.sqrt((p.K2[None]**2*np.abs(uh)**2).sum()/tot))
    # spin-direction covariance of the velocity gradient: eigenvalues of mean A^T A tell which directions the flow varies along
    A = p.A(uh); G = np.einsum('ijX,ikX->jk', A.reshape(p.d, p.d, -1), A.reshape(p.d, p.d, -1))/A[0, 0].size
    return dict(flat_share=flat, k_rms=kr, vary_eig=sorted(np.linalg.eigvalsh(G).tolist()))
if __name__ == '__main__':
    d, N, nu, seed, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    Z0 = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    iters = int(sys.argv[7]) if len(sys.argv) > 7 else 500
    p = Prob(d, N, nu); t0 = time.time()
    uh, R, P, D, hist = ascend(p, init(p, seed, Z0), Z0, iters)
    res = dict(d=d, N=N, nu=nu, seed=seed, Z0=Z0, R=float(R), P=float(P), D=float(D), iters=len(hist),
               hist=[(int(a), float(b), float(c), float(e)) for a, b, c, e in hist[::5]], secs=time.time() - t0, **diag(p, uh))
    json.dump(res, open(out, 'w'))
    import os
    if os.environ.get('SAVE'): np.savez_compressed(out.replace('.json', '_u.npz'), u=np.array([p.ifft(c) for c in uh]).astype(np.float32))
    print(f'd={d} N={N} nu={nu} seed={seed}: R={R:.5g} P={P:.5g} D={D:.5g} iters={len(hist)} k_rms={res["k_rms"]:.2f} flat={np.round(res["flat_share"],3)} {time.time()-t0:.0f}s')
