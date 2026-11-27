"""Fig 30: one function F whose imaginary part is s*G, and where its zeros sit. Run: python3 figures_hb.py"""
import os, math, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
HERE = os.path.dirname(os.path.abspath(__file__))
N, W = leggauss(3000); X = (N+1)/2; W = W/2
def make_F(p, a):
    rho = W*a*p*X**(p-1)*(1-X**p)**(a-1)
    def F(w):
        w = np.atleast_1d(np.asarray(w, complex)); out = np.empty(w.shape, complex)
        for i in range(0, len(w), 2000): out[i:i+2000] = np.exp(1j*np.outer(w[i:i+2000], X)) @ rho
        return out
    return F, (lambda w: np.exp(1j*w*X) @ (rho*1j*X))
height = lambda p, a, w: math.log(p**a*math.gamma(a+1)/(a*math.gamma(p+1))) + (p-a)*math.log(abs(w))
def zeros(p, a, R):
    F, dF = make_F(p, a); out = []
    for x in np.arange(6, R, 1.0):
        w = x + 1j*height(p, a, x)
        h0 = w.imag; conv = False
        for _ in range(60):
            d = dF(w)
            if not np.isfinite(d) or d == 0: break
            step = F(w)[0]/d; w = w - step
            if abs(w.imag - h0) > 40: break
            if abs(step) < 1e-4: conv = True; break
        if conv and 0 < w.real < R and abs(F(w)[0]) < 1e-3*np.min(np.abs(F(np.array([w+0.7, w-0.7, w+0.7j, w-0.7j])))) \
           and all(abs(w-z) > 1e-5 for z in out): out.append(w)
    return np.array(sorted(out, key=lambda z: z.real))
fig, ax = plt.subplots(2, 3, figsize=(16, 10))

# A: the densities behind F
a_ = ax[0, 0]; x = np.linspace(0, 1, 2000)
for p, col in ((1.5, 'tab:purple'), (3, 'tab:blue'), (10, 'tab:green'), (20.5, 'tab:orange')):
    a_.plot(x, p*p*x**(p-1)*(1-x**p)**(p-1), color=col, label=f'p = a = {p}')
a_.set_xlabel('x'); a_.set_ylabel('density'); a_.legend(fontsize=8)
a_.set_title('A. The weight seen from the rim: X = U^(1/p), U ~ Beta(1, a).\nF(w) = average of e^{iwX}, and s·G(s) = Im F(s)')

# B: |F| over the plane, p = 10, a = 10, with its zeros and the real line
b = ax[0, 1]; p = 10; F, _ = make_F(p, p)
U, V = np.meshgrid(np.linspace(0, 160, 321), np.linspace(-10, 30, 161)); Z = np.abs(F((U+1j*V).ravel())).reshape(U.shape)
cs = b.contourf(U, V, np.log10(Z+1e-30), 40, cmap='viridis'); plt.colorbar(cs, ax=b, label='log10 |F|')
z = zeros(p, p, 160); b.plot(z.real, z.imag, 'o', color='w', mec='k', ms=6)
b.axhline(0, color='r', lw=2); b.axhline((p-1)*math.log(p), color='w', ls='--', lw=1)
b.text(5, (p-1)*math.log(p)+1.2, '(p − 1) ln p', color='w'); b.text(85, -7, 'real line: G lives here', color='r')
b.set_xlabel('Re w'); b.set_ylabel('Im w'); b.set_title('B. p = a = 10: every zero of F (white) sits about\n20 above the real line; none below it')

# C: zero heights against distance, for a below, at and above p
c = ax[0, 2]
for (p, a, col, ls) in ((3, 2.4, 'tab:blue', ':'), (3, 3, 'tab:blue', '-'), (3, 3.6, 'tab:blue', '--'),
                        (10, 9.3, 'tab:green', ':'), (10, 10, 'tab:green', '-'), (10, 12, 'tab:green', '--')):
    z = zeros(p, a, 250); c.plot(z.real, z.imag, ls, marker='o', ms=3, color=col, label=f'p={p}, a={a}')
w = np.linspace(10, 250, 200)
for p, a in ((3, 3), (10, 10)): c.plot(w, [height(p, a, q) for q in w], color='k', lw=0.7)
c.axhline(0, color='r', lw=1.5); c.set_xscale('log'); c.legend(fontsize=7, ncol=2)
c.set_xlabel('Re w (log)'); c.set_ylabel('height of zero'); c.set_title('C. Below a = p the zeros climb, at a = p they level off,\nabove it they sink like −(a − p) ln|w| (thin black: formula)')

# D: the cushion grows with p
d = ax[1, 0]; ps = [1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9.5, 10.5, 12.5, 15.3]
low_p = []; low_pm = []; psA = []; psB = []
for p in ps:
    z1 = zeros(p, p, 250); z2 = zeros(p, max(p-0.75, p/2), 250)
    print(p, len(z1), len(z2), flush=True)
    if len(z1): psA.append(p); low_p.append(min(z1.imag))
    if len(z2): psB.append(p); low_pm.append(min(z2.imag))
pp = np.linspace(1.05, 21, 200)
d.plot(pp, (pp-1)*np.log(pp), 'k', label='(p − 1) ln p')
d.plot(psA, low_p, 'o', color='tab:green', label='lowest zero, a = p'); d.plot(psB, low_pm, 's', color='tab:blue', mfc='none', label='lowest zero, a = p − 0.75 (p/2 below p = 2)')
d.axvline(8, color='0.6', ls=':'); d.text(8.3, 2, 'past 8: still open', color='0.4')
d.set_xlabel('roundness p'); d.set_ylabel('distance from the real line'); d.legend(fontsize=8)
d.set_title('D. The cushion: how far F\'s nearest zero sits\nfrom the real line, across the top strip')

# E: even p above a = p: the criterion fails though G is real
e = ax[1, 1]; z = zeros(4, 6, 300); e.plot(z.real, z.imag, 'o-', ms=3, color='tab:red', label='p = 4, a = 6 (G real: proved)')
z2 = zeros(2.5, 5, 300); e.plot(z2.real, z2.imag, 's-', ms=3, color='tab:orange', label='p = 2.5, a = 5 (G breaks)')
e.axhline(0, color='k'); e.legend(fontsize=8); e.set_xscale('log')
e.set_xlabel('Re w (log)'); e.set_ylabel('height of zero')
e.set_title('E. Above a = p the criterion fails at every p,\neven where G stays real (p = 4): sufficient, not necessary')

# F: where it would prove the rule
f = ax[1, 2]; pg = np.linspace(1, 16, 400)
f.fill_between(pg, 0, pg, color='tab:blue', alpha=.25, label='criterion holds (all checks): a ≤ p')
for k in range(2, 17, 2): f.axvline(k, color='tab:green', lw=1)
f.fill_between(pg[pg <= 8], 0, pg[pg <= 8], color='none', hatch='//', edgecolor='0.35', label='already proved (p ≤ 8)')
f.plot(pg, pg, 'k'); f.set_xlim(1, 16); f.set_ylim(0, 16)
f.set_xlabel('roundness p'); f.set_ylabel('a = (n − 1)/p'); f.legend(fontsize=8, loc='upper left')
f.set_title('F. Proving "no zeros of F below the line" for a ≤ p\nwould settle the rule\'s holding side at every p (green: even p)')
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'fig30-lee-yang-hermite-biehler.png'), dpi=130); print('saved fig30')
