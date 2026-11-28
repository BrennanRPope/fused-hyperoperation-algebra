"""DRAFT (round 4). Can motion hidden in a 4th direction slow the spin growth that a 3D view shows?

Take a 3D flow U(x) (here the fastest-growing 3D flow: colliding vortex rings) and add hidden motion u' that varies
along a 4th direction w as cos/sin(kw w), so it averages to zero over w. The 3D view is the w-average, which is U.
Exact (derived here, checked numerically below): the hidden motion changes the growth of the 3D view's spin by
    Q(u') = + mean( u'_i u'_j d_j V_i ),   V = -Lap U   (Q < 0 is a drain),
so it is quadratic in u' and only the symmetric part S of grad V matters. Write u' = Re[(a(x), b(x)) e^{i kw w}],
with b = i div(a)/kw (incompressible in 4D). Then Q = (1/2) mean(a.S.a) and the hidden spin is
    Z' = mean|grad u'|^2 = (1/2) sum_k (|k|^2 + kw^2)(|a_k|^2 + |k.a_k|^2/kw^2).
We find the hidden motion with the strongest drain per unit hidden spin (largest lambda in -(1/2) S a = lambda M a; lambda > 0 is a drain), then
measure, in a full 4D grid, what that motion does to the 3D view, to itself, and to the total.

usage: python3 drain.py rings_u.npz nu out.json
"""
import sys, json, numpy as np, scipy.fft as sf
from scipy.sparse.linalg import LinearOperator, eigsh
W = 4
fft = lambda a, ax: sf.fftn(a, axes=ax, workers=W); ifft = lambda a, ax: sf.ifftn(a, axes=ax, workers=W).real
def load(f):
    U = np.load(f)['u'].astype(float); N = U.shape[1]
    k1 = np.fft.fftfreq(N, 1/N); K = np.array(np.meshgrid(k1, k1, k1, indexing='ij')); return U, N, K
def strain_of_V(U, N, K):
    Uh = fft(U, (1, 2, 3)); K2 = (K**2).sum(0); Vh = K2[None]*Uh           # V = -Lap U
    G = np.array([[ifft(1j*K[j]*Vh[i], (0, 1, 2)) for j in range(3)] for i in range(3)])   # G[i,j] = d_j V_i
    return 0.5*(G + G.transpose(1, 0, 2, 3, 4))
def best_drain(S, N, K, kw, nev=3):
    keep = np.all(np.abs(K) <= N//4, axis=0)
    K2 = (K**2).sum(0); n = 3*N**3
    def Aop(x):                                                    # -S a, restricted to kept modes (symmetric)
        a = x.reshape(3, N, N, N); ah = fft(a, (1, 2, 3))*keep; a = ifft(ah, (1, 2, 3))
        y = -np.einsum('ijX,jX->iX', S.reshape(3, 3, -1), a.reshape(3, -1)).reshape(3, N, N, N)
        return (ifft(fft(y, (1, 2, 3))*keep, (1, 2, 3))*0.5).ravel()
    def Mop(x, inv=False):                                         # hidden spin as a quadratic form in a
        ah = fft(x.reshape(3, N, N, N), (1, 2, 3))
        c = (K2 + kw**2)[None]
        if not inv:
            yh = c*(ah + K*(K*ah).sum(0)/kw**2)*0.5
        else:                                                      # (I + k k^T/kw^2)^-1 = I - k k^T/(kw^2 + |k|^2)
            yh = (ah - K*(K*ah).sum(0)/(kw**2 + K2))/c*2
        return ifft(yh, (1, 2, 3)).ravel()
    A = LinearOperator((n, n), matvec=Aop, dtype=float)
    M = LinearOperator((n, n), matvec=Mop, dtype=float)
    Minv = LinearOperator((n, n), matvec=lambda x: Mop(x, True), dtype=float)
    vals, vecs = eigsh(A, k=nev, M=M, Minv=Minv, which='LA', tol=1e-6, maxiter=3000)
    i = np.argmax(vals); a = vecs[:, i].reshape(3, N, N, N)
    a = ifft(fft(a, (1, 2, 3))*keep, (1, 2, 3))
    return float(vals[i]), a
def rate4(u, N4, nu):
    """(1/2) d/dt mean|grad u|^2 for a 4D field (full grid), and the same for its w-average and for the rest"""
    ax = (1, 2, 3, 4); k1 = np.fft.fftfreq(N4[0], 1/N4[0]); kw1 = np.fft.fftfreq(N4[3], 1/N4[3])
    K = np.array(np.meshgrid(k1, k1, k1, kw1, indexing='ij')); K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
    uh = fft(u, ax)
    sax = (0, 1, 2, 3)
    A = np.array([[ifft(1j*K[j]*uh[i], sax) for j in range(4)] for i in range(4)])
    adv = np.einsum('jX,ijX->iX', u.reshape(4, -1), A.reshape(4, 4, -1)).reshape(u.shape)    # (u.grad)u
    Nh = -fft(adv, ax); Nh = Nh - K*(K*Nh).sum(0)/K2s                                         # Leray-projected
    dudt = Nh - nu*K2[None]*uh
    tot = (K2[None]*(uh.conj()*dudt).real).sum()/uh[0].size**2
    m0 = (K[3] == 0)
    vis = (K2[None]*(uh.conj()*dudt).real)[:, m0].sum()/uh[0].size**2
    Z = (K2[None]*np.abs(uh)**2).sum()/uh[0].size**2; Zv = (K2[None]*np.abs(uh)**2)[:, m0].sum()/uh[0].size**2
    return dict(total=float(tot), view=float(vis), hidden=float(tot - vis), Z=float(Z), Z_view=float(Zv), Z_hidden=float(Z - Zv))
if __name__ == '__main__':
    f, nu, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    U, N, K = load(f); S = strain_of_V(U, N, K)
    res = dict(nu=nu, N=N)
    for kw in [1.0, 0.5, 0.25]:
        lam, a = best_drain(S, N, K, kw)
        res[f'lambda_kw{kw}'] = lam
        if kw == 1.0: a1 = a
        print(f'kw={kw}: strongest drain per unit hidden spin lambda = {lam:.4g}', flush=True)
    # full 4D accounting for kw = 1 (w-period 2*pi), Nw points along w
    Nw = 8; w = np.arange(Nw)*2*np.pi/Nw
    ah = fft(a1, (1, 2, 3)); b = ifft(1j*1j*(K*ah).sum(0), (0, 1, 2))   # b = i div(a); as real field: u'_w = Re(i div a e^{iw}) = -div a sin w
    base = np.zeros((4, N, N, N, Nw)); base[:3] = U[..., None]
    hid = np.zeros_like(base)
    hid[:3] = a1[..., None]*np.cos(w)[None, None, None, None]
    diva = ifft(1j*(K*ah).sum(0), (0, 1, 2))
    hid[3] = -diva[..., None]*np.sin(w)                                  # d_w(-div a sin w) = -div a cos w : 4D divergence 0
    r0 = rate4(base, (N, N, N, Nw), nu)
    Zh_unit = rate4(hid, (N, N, N, Nw), nu)['Z']
    lam = res['lambda_kw1.0']; R3 = r0['view']
    eps_cancel = R3/lam                                                  # hidden spin that cancels the view's growth (to first order)
    res.update(R3_view_alone=R3, Z_view=r0['Z_view'], hidden_spin_to_cancel=eps_cancel)
    scan = []
    for frac in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 4.0]:
        eps = frac*eps_cancel; u = base + hid*np.sqrt(eps/Zh_unit)
        r = rate4(u, (N, N, N, Nw), nu); r['hidden_spin'] = eps; r['frac_of_cancel'] = frac; scan.append(r)
        print(f"hidden spin {eps:.4f} ({frac}x): view growth {r['view']:+.4f}  hidden growth {r['hidden']:+.4f}  total {r['total']:+.4f}", flush=True)
    res['scan'] = scan
    json.dump(res, open(out, 'w'), indent=1)
