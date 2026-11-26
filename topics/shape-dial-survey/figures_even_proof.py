"""Figure 28: even roundness settled from Polya's end, and the fused-integral test."""
import numpy as np, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
mp.mp.dps = 60

def coeffs(p, a, K=200):
    p = mp.mpf(p); a = mp.mpf(a)
    return [mp.gamma(a+1)*mp.gamma((2*k+1)/p)/(p*mp.factorial(2*k)*mp.gamma((2*k+1)/p+a+1)) for k in range(K)]
def G(c, s):
    w = -s*s; acc = mp.mpf(0)
    for ck in reversed(c): acc = acc*w + ck
    return acc
def phase_grid(f, X, Y):
    Z = np.zeros((len(Y), len(X)), dtype=complex)
    for i, y in enumerate(Y):
        for j, x in enumerate(X):
            v = f(mp.mpc(x, y)); Z[i, j] = complex(v/abs(v)) if v != 0 else 0
    return Z
def show_phase(ax, Z, X, Y, title):
    h = (np.angle(Z)/(2*np.pi)) % 1.0
    rgb = hsv_to_rgb(np.dstack([h, 0.75*np.ones_like(h), 0.95*np.ones_like(h)]))
    ax.imshow(rgb, extent=[X[0], X[-1], Y[0], Y[-1]], origin='lower', aspect='auto')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_title(title, fontsize=10); ax.set_xlabel('Re s'); ax.set_ylabel('Im s')

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.28)

# A: the scaling factor q(z) = 1/Gamma(rho z + beta)
ax = fig.add_subplot(gs[0, 0])
p, a = 10, 3; rho = 2/p; beta = a+1+1/p
z = np.linspace(-40, 14, 2000)
qv = np.array([float(1/mp.gamma(rho*t+beta)) for t in z])
ax.plot(z, np.sign(qv)*np.log10(1+np.abs(qv)*1e3), color='tab:blue')
zeros = [-(beta+n)/rho for n in range(0, 4)]
ax.plot(zeros, [0]*len(zeros), 'o', color='tab:red', label='zeros of q (all negative)')
ks = np.arange(0, 14); ax.plot(ks, [np.log10(1+float(1/mp.gamma(rho*k+beta))*1e3) for k in ks], 's', color='tab:green', ms=4, label='q(k), k = 0, 1, 2, ...')
ax.axvline(0, color='grey', lw=0.6); ax.axhline(0, color='grey', lw=0.6)
ax.set_title('A. The scaling factor q(z) = 1/Γ(2z/p + a + 1 + 1/p)\np = 10, a = 3 (signed log scale)', fontsize=10)
ax.set_xlabel('z'); ax.legend(fontsize=8, loc='lower left')

# B: phase portraits, even p past 8
X = np.linspace(0.2, 30, 220); Y = np.linspace(-5, 5, 70)
for col, (pp, aa, ttl) in enumerate([(10, 25, 'B. p = 10, a = 25 (above a = p): every zero on the line'),
                                      (2.5, 5, 'C. p = 2.5, a = 5 (not even, above a = p):\nzeros leave the line past the crossover')]):
    c = coeffs(pp, aa)
    Zp = phase_grid(lambda s: G(c, s), X, Y)
    show_phase(fig.add_subplot(gs[0, col+1]), Zp, X, Y, ttl)

# D: the map
ax = fig.add_subplot(gs[1, 0])
P = np.linspace(1, 16, 400)
P8 = P[P <= 8]; P16 = P[P >= 8]
ax.fill_between(P8, 0, P8, color='#9ecae1', alpha=0.8, label='a ≤ p, 1.1 ≤ p ≤ 8: proved (near + far)')
ax.fill_between(P16, 0, P16, color='#d9d9d9', alpha=0.8, label='a ≤ p, p > 8 not even: open')
ax.fill_between(P, P, 18, color='#fcbba1', alpha=0.7, label='a > p, p not even: breaks')
for pe in [2, 4, 6, 8, 10, 12, 14, 16]:
    ax.plot([pe, pe], [0, 18], color='#006d2c', lw=3 if pe > 2 else 3)
ax.plot([], [], color='#006d2c', lw=3, label='even p, every a: proved (this section)')
ax.set_xlim(1, 16.3); ax.set_ylim(0, 18); ax.set_xlabel('roundness p'); ax.set_ylabel('a = (n − 1)/p')
ax.set_title('D. The (p, a) map now: every even line is settled', fontsize=10); ax.legend(fontsize=7.5, loc='upper left')

# E: fused integral at the multadd level, (e^{kG} - 1)/k with k = 1
c = coeffs(10, 3)
X2 = np.linspace(0.2, 20, 200); Y2 = np.linspace(-6, 6, 90)
Zf = phase_grid(lambda s: (mp.exp(G(c, s)) - 1), X2, Y2)
show_phase(fig.add_subplot(gs[1, 1]), Zf, X2, Y2, 'E. Fused version (e^{kG} − 1)/k, k = 1, p = 10, a = 3:\nsame zeros on the line, new ones off it')
Zg = phase_grid(lambda s: G(c, s), X2, Y2)
show_phase(fig.add_subplot(gs[1, 2]), Zg, X2, Y2, 'F. Plain G at the same p, a for comparison')
fig.suptitle('Figure 28. Even roundness from Pólya’s end: G is Pólya’s function with its k-th coefficient scaled by q(k).\n'
             'q has zeros only on the negative side, so Laguerre’s theorem keeps every zero real. Zeros show where all colours meet.', fontsize=11)
fig.savefig('figures/fig28-lee-yang-even-proof.png', dpi=110, bbox_inches='tight')
print('ok')
