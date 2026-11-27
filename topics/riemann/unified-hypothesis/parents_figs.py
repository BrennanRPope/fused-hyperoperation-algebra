"""Pictures for the parent hunt: Φ parents scan, antiderivative zeros, multiplier tests."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
import kern

P = json.load(open('parents.json'))
base = len(P['H0_real_zeros'])
sc = P['scans']
# ---- 1. parent scan
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
def series(kind):
    xs, ys = [], []
    for k, v in sc.items():
        kk, par = k.split('|')
        if kk == kind:
            xs.append(float(par)); ys.append(len(v))
    o = np.argsort(xs); return np.array(xs)[o], np.array(ys)[o]
x, y = series('heat')
axs[0].plot(x, y, 'o-', color='tab:green'); axs[0].set_xlabel('rewind time s (parent kernel Φ·e^(−s u²))')
axs[0].set_title('Heat parents of Φ: rewind by s\nlow zeros all stay real (Rodgers–Tao: high ones do not, for every s > 0)', fontsize=9)
x, y = series('factor')
axs[1].plot(x, y, 'o-', color='tab:orange'); axs[1].axvline(np.pi/8, color='gray', ls=':')
axs[1].text(np.pi/8 + 0.03, 3, 'c = π/8: pole meets the edge of Φ’s strip', fontsize=8)
axs[1].set_xscale('log'); axs[1].set_xlabel('c (parent kernel Φ / (1 + u²/c²))')
axs[1].set_title('Factor parents: one Laguerre step back\nall 18 low zeros real only for c ≥ 1', fontsize=9)
x, y = series('bath')
axs[2].plot(x, y, 'o-', color='tab:purple'); axs[2].axvline(4, color='gray', ls=':')
axs[2].text(4.1, 3, 'J = 4: pole iπ/(2J) meets the strip', fontsize=8)
axs[2].set_xlabel('J (parent kernel Φ / cosh(J u), the two-level bath undone)')
axs[2].set_title('Bath parents: undo a half-filled two-level bath\nall 18 low zeros real only for J ≤ 1', fontsize=9)
for a in axs:
    a.axhline(base, color='k', lw=0.8); a.set_ylabel('real zeros with z < 150'); a.set_ylim(-1, base + 2)
    a.text(a.get_xlim()[0], base + 0.6, f' Ξ itself: {base}', fontsize=8)
fig.suptitle('Parents of Φ tried directly: count of real zeros below z = 150 (Ξ has 18 there). A parent that loses any is not real-rooted.')
fig.tight_layout(); fig.savefig('figures/fig_parents_of_phi.png', dpi=130); plt.close(fig)

# ---- 2. antiderivative parent
R = kern.riemann(600)
L = np.pi / 2 * kern.Phi(np.array([0.0]))[0]
x = np.arange(0.01, 74, 0.002)
F = np.sin(np.outer(x, R.u)) @ (R.w / R.u) - L
real = x[np.where(np.sign(F[:-1]) * np.sign(F[1:]) < 0)[0]]
off = np.array(json.load(open('anti_offaxis_plus.json'))['zeros'])
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
ax = axs[0]
ax.plot(real, 0 * real, 'o', color='tab:blue', label='zeros of the parent F = ∫H₀ − L (real, z > 0)')
g = kern.real_zeros(R, 74)
ax.plot(g, 0 * g, 'x', color='k', ms=9, label='zeros of Ξ (z = 2γ)')
ax.plot(-off[:, 0], off[:, 1], 's', color='tab:red', label='zeros of F off the line (left side)')
ax.plot(-off[:, 0], -off[:, 1], 's', color='tab:red')
ax.axhline(0, color='gray', lw=0.6); ax.axvline(0, color='gray', lw=0.6)
ax.set_xlabel('Re z'); ax.set_ylabel('Im z'); ax.legend(fontsize=8, loc='lower right')
ax.set_title('Antiderivative parent: real zeros interlace with Ξ’s on the right,\nbut a second string of zeros sits off the line on the left', fontsize=9)
ax = axs[1]
xx = np.linspace(0.5, 74, 2000)
FF = np.sin(np.outer(xx, R.u)) @ (R.w / R.u) - L
FL = -np.sin(np.outer(xx, R.u)) @ (R.w / R.u) - L
ax.plot(xx, np.sign(FF) * np.log10(1 + np.abs(FF) * 1e9), label='F(z), z > 0 (signed log scale)')
ax.plot(-xx, np.sign(FL) * np.log10(1 + np.abs(FL) * 1e9), label='F(z), z < 0')
ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel('z (real)'); ax.set_title('Along the real line: F → 0 on the right (crossing at each zero),\nF → −2L on the left (never crosses)', fontsize=9)
ax.legend(fontsize=8)
fig.suptitle('Rolle parent of Ξ: F′ = H₀. Jensen’s formula needs as many zeros as Ξ has, but the real ones are on one side only, so the rest are off the line')
fig.tight_layout(); fig.savefig('figures/fig_antiderivative_parent.png', dpi=130); plt.close(fig)

# ---- 3. multiplier tests: first Turán quantities
M1 = json.load(open('multiplier.json')); M2 = json.load(open('multiplier2.json'))
fig, ax = plt.subplots(figsize=(10, 4.8))
rows = [('Φ* (Pólya 1926)', M1['Phi_over_PhiStar']['turan_fail'], 40, M1['Phi_over_PhiStar']['jensen_fail'])]
for p in ['2', '4', '6', '8', '12']:
    rows.append((f'Pólya bump e^(−u^{p})' + (' = Gaussian: same as RH' if p == '2' else ''), M2[p]['turan_fail'], 120, M2[p]['jensen_fail']))
rows.append(('control: p-ball over e^(−u⁴), should pass', M1['control_pball_over_polya_p4_a2']['turan_fail'], 40, []))
rows.append(('control: 1/(k+1), should fail', M1['control_1_over_k_plus_1']['turan_fail'], 40, M1['control_1_over_k_plus_1']['jensen_fail']))
for i, (nm, fails, kmax, jf) in enumerate(rows):
    ks = np.arange(1, kmax)
    ok = [k for k in ks if k not in fails]
    ax.plot(ok, [i] * len(ok), '|', color='tab:green', ms=10)
    ax.plot(fails, [i] * len(fails), '|', color='tab:red', ms=10)
    ax.text(kmax + 2, i, f'Jensen degree 2 {"fails" if any(d == 2 for d, _ in jf) else "holds"}', va='center', fontsize=8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xlabel('k  (green: Turán r_k² ≥ r_(k−1) r_(k+1) holds; red: fails)')
ax.set_xlim(0, 150)
ax.set_title('Can Φ be a Laguerre-multiplier child of a known real-rooted parent?\nr_k = (moment 2k of Φ) / (moment 2k of the parent). A multiplier sequence must pass every Turán test.', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_multiplier_parents.png', dpi=130); plt.close(fig)
print('ok')
