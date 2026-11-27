"""Fig 29: the pendulum picture against the Lee-Yang work. Run: python3 figures_pendulum_leeyang.py"""
import os, pickle, numpy as np, mpmath as mp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipj
HERE = os.path.dirname(os.path.abspath(__file__)); LY = os.path.join(HERE, '..', 'lee-yang')
kc = 3+2*np.sqrt(2); kO = lambda k: 4*k/(k-1)**2
fig, ax = plt.subplots(2, 3, figsize=(16, 10))

# A: energy dial along the dial gas temperature
a = ax[0, 0]; k = np.geomspace(1.02, 200, 800); m = kO(k)**2
a.semilogx(k, np.log10(m), color='k', lw=2)
a.axhline(0, color='tab:orange', lw=1.5); a.axvline(kc, color='tab:orange', ls='--', lw=1)
a.fill_between(k, -4, 4, where=k < kc, color='tab:purple', alpha=.08)
a.fill_between(k, -4, 4, where=k > kc, color='tab:blue', alpha=.08)
kk = k[k > kc]; a.semilogx(kk, (1-kO(kk)**2)**(1/8)*2-2, color='tab:green', lw=1.5)
a.text(1.9, 3.0, 'over the top\n(disordered, m > 1)', color='tab:purple', fontsize=10)
a.text(15, 2.9, 'swinging\n(ordered, m < 1)', color='tab:blue', fontsize=10)
a.text(kc*1.08, -3.3, r'$k_c=3+2\sqrt{2}$: separatrix, $m=1$', color='tab:orange', fontsize=9)
a.text(30, -0.55, r"Yang's $M=(1-m)^{1/8}$ (scaled)", color='tab:green', fontsize=9)
a.set_ylim(-4, 4); a.set_xlabel('dial k = e^{4K}'); a.set_ylabel(r'$\log_{10} m$,  $m = k_O^2$')
a.set_title('A. The dial gas transition is the separatrix\n(Kramers–Wannier swaps m and 1/m)')

# B: Onsager via the pendulum period against the dial gas torus
b = ax[0, 1]; C = pickle.load(open(os.path.join(LY, 'torus_counts.pkl'), 'rb'))
def onsager(k):
    K = np.log(k)/4; ko = kO(k); kap = np.minimum(ko, 1/ko)
    return 1/np.tanh(2*K)/2*(1+2/np.pi*(2*np.tanh(2*K)**2-1)*(1+kap)*ellipk(kap**2))
kg = np.geomspace(1.15, 40, 600); kg = kg[np.abs(kg-kc) > 1e-3]
b.semilogx(kg, onsager(kg), 'k', lw=2, label='Onsager, as (1+k_O) × pendulum period')
def torus(L, k):
    c = C[L].astype(float); N = L*L; mm = np.arange(c.shape[0])[:, None]; e = np.arange(c.shape[1])[None, :]
    lw = np.where(c > 0, np.log(np.where(c > 0, c, 1)) - 2*mm*np.log(k) + e*np.log(k), -np.inf)
    w = np.exp(lw - lw.max()); return 2*(w*e).sum()/w.sum()/N - 1
for L, col in zip((4, 5, 6), ('tab:red', 'tab:orange', 'tab:green')):
    ks = np.geomspace(1.2, 40, 30); b.semilogx(ks, [torus(L, x) for x in ks], 'o', ms=3.5, color=col, label=f'dial gas torus L={L}')
b.axvline(kc, color='tab:orange', ls='--', lw=1); b.legend(fontsize=8)
b.set_xlabel('dial k'); b.set_ylabel(r'$\langle\sigma\sigma\rangle$ (nearest neighbour)')
b.set_title('B. The dial gas on the torus follows Onsager,\nwhose integral is the pendulum quarter period')

# C: Fisher zeros of the dial gas against |k_O| = 1
c = ax[0, 2]; F = pickle.load(open(os.path.join(LY, 'fisher_zeros.pkl'), 'rb'))
for L, col, mk in zip((4, 5, 6), ('tab:red', 'tab:orange', 'tab:green'), ('s', '^', 'o')):
    z = np.array(F[L]); c.plot(z.real, z.imag, mk, ms=5, color=col, mfc='none', label=f'Fisher zeros L={L}')
phi = np.linspace(0, 2*np.pi, 2000); P = []
for f in phi:
    w = np.exp(1j*f); P.append(np.roots([w, -(2*w+4), w]))
P = np.array(P)
for j in range(2): c.plot(P[:, j].real, P[:, j].imag, '.', ms=1, color='k')
c.plot([kc, 1/kc], [0, 0], '*', color='tab:orange', ms=12, label=r'$k_c$, $1/k_c$')
c.set_xlim(-8, 10); c.set_ylim(-9, 9); c.set_aspect('equal'); c.legend(fontsize=8, loc='upper left')
c.set_xlabel('Re k'); c.set_ylabel('Im k')
c.set_title(r'C. Complex-temperature zeros close in on $|k_O|=1$:' + '\nthe separatrix |m| = 1 in complex energy (black)')

# D: pendulum reading in complex time; |G| = 1 on the halfway lines; zeros of G^{2N}+1 on them
d = ax[1, 0]; mm = 0.7; K = ellipk(mm); Kp = ellipk(1-mm)
X, Y = np.meshgrid(np.linspace(-K, K, 400), np.linspace(0, 2*Kp, 400)); U = X+1j*Y
Gf = np.vectorize(lambda u: complex(mp.ellipfun('sn', u, m=mm)*mp.ellipfun('dn', u, m=mm)/mp.ellipfun('cn', u, m=mm)))
mp.mp.dps = 15
Gv = Gf(U[::4, ::4]); lg = np.log10(np.abs(Gv)+1e-30)
d.contourf(X[::4, ::4]/K, Y[::4, ::4]/Kp, np.clip(lg, -2, 2), 40, cmap='RdBu_r')
d.contour(X[::4, ::4]/K, Y[::4, ::4]/Kp, lg, [0], colors='k', linewidths=1.5)
N = 3; pts = []
for j in range(2*N):
    t = mp.exp(1j*mp.pi*(2*j+1)/(2*N))
    for x in np.linspace(-.9, .9, 7):
        for y in np.linspace(.05, 1.95, 7):
            try:
                r = mp.findroot(lambda u: mp.ellipfun('sn', u, m=mm)*mp.ellipfun('dn', u, m=mm)/mp.ellipfun('cn', u, m=mm)-t, K*x+1j*Kp*y)
                pts.append(((float(mp.re(r)/K)+1) % 2-1, float(mp.im(r)/Kp) % 2))
            except Exception: pass
pts = np.unique(np.round(pts, 5), axis=0); d.plot(pts[:, 0], pts[:, 1], 'o', color='lime', mec='k', ms=7)
d.set_xlabel('real time / K'); d.set_ylabel("imaginary time / K'")
d.set_title('D. Pendulum reading |G| (red > 1, blue < 1), m = 0.7:\n|G| = 1 on the halfway lines; zeros of $G^6+1$ (green) sit on them')

# E: 1D Ising ring Lee-Yang zeros = equal-time samples of a swing
e = ax[1, 1]; th = np.linspace(0, 2*np.pi, 400); e.plot(np.cos(th), np.sin(th), color='0.8')
Nr = 24
for kk, col, r in ((2, 'tab:red', 1.0), (20, 'tab:blue', 1.0)):
    cc = np.sqrt(1-1/kk); psi = np.pi*(2*np.arange(Nr)+1)/(2*Nr); half = np.arccos(cc*np.cos(psi))
    thj = np.concatenate([2*half, -2*half]); s = 1.0 if kk == 2 else 0.86
    e.plot(s*np.cos(thj), s*np.sin(thj), 'o', ms=4, color=col, label=f'k={kk}: swing size c = {cc:.3f}')
e.set_aspect('equal'); e.legend(fontsize=8, loc='center')
e.text(1.05, 0.05, 'arc end\n(zeros crowd,\nlike a turning point)', fontsize=8)
e.set_title('E. 1D ring, N = 24: cos(θ/2) = c cos(ψ), ψ evenly spaced.\nColder k = bigger swing; bigger N = finer sampling')
e.set_xlim(-1.3, 1.6); e.set_ylim(-1.2, 1.2); e.axis('off')

# F: anisotropic couplings ride the pendulum's imaginary time; isotropic = halfway
f = ax[1, 2]
for ko, col in ((0.3, 'tab:blue'), (0.7, 'tab:green'), (0.97, 'tab:orange')):
    Kp = ellipk(1-ko**2); u = np.linspace(0.01, 0.99, 300)*Kp
    v = np.array([float(mp.re(-1j*mp.ellipfun('sn', 1j*x, m=ko**2))) for x in u])
    f.semilogy(u/Kp, np.sqrt(ko)*v, color=col, label=f'k_O = {ko}')
u = np.linspace(0.01, 0.99, 300); f.semilogy(u, np.tan(u*np.pi/2), 'k--', label='k_O = 1 (critical): tan')
f.axvline(0.5, color='0.5', ls=':'); f.axhline(1, color='0.5', ls=':')
f.set_xlabel("imaginary time u / K'"); f.set_ylabel(r'$\sqrt{k_O}\,\sinh 2K_1$')
f.legend(fontsize=8); f.set_title('F. Couplings in imaginary time: sinh2K₁ = −i sn(iu),\nsinh2K₂ at K′ − u; equal at halfway (isotropic)')
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'fig29-pendulum-lee-yang.png'), dpi=130)
print('saved fig29')
