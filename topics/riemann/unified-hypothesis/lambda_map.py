import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
d = json.load(open('lam_pball.json'))
fig, axs = plt.subplots(1, 2, figsize=(15, 5.5))
ax = axs[0]
pts = []
for k, v in d.items():
    p, a = map(float, k.split('|'))
    lam = v.get('Lambda')
    if lam is None and p == 1 and a == 1:
        lam = 0.0
    if lam is not None:
        pts.append((p, a, lam))
pts = np.array(pts)
sc = ax.scatter(pts[:, 0], pts[:, 1], c=pts[:, 2], cmap='viridis_r', s=150, edgecolor='k', vmin=-8.5, vmax=0)
for p, a, l in pts:
    ax.text(p + 0.08, a + 0.08, f'{l:.2f}', fontsize=7)
pp = np.linspace(1, 6.5, 50)
ax.plot(pp, pp, 'r--', lw=1); ax.text(5.3, 5.7, 'a = p (rule boundary)', color='r', fontsize=8, rotation=35)
ax.fill_between(pp, pp, 9, color='red', alpha=0.06)
ax.text(1.3, 7.6, 'above the line, non-even p: Λ = +∞\n(rim wave beats the centre kink at every heat time)', fontsize=8, color='darkred')
ax.set_xlabel('p (ball shape)'); ax.set_ylabel('a = (n−1)/p (weight power)')
plt.colorbar(sc, ax=ax, label='measured rewind constant Λ')
ax.set_title('p-ball kernels (1−|u|^p)^a: rewind constant Λ where the rule holds\n(even p = 2, 4, 6 hold above the line too; triangle p = a = 1 has Λ = 0 exactly)', fontsize=9)
ax.set_ylim(0, 8.8)
ax = axs[1]
ax.set_xlim(-4.2, 4.2); ax.set_ylim(0, 4); ax.axvline(0, color='k', lw=1)
rows = [
    (3.4, 'Pólya bumps e^(−u^p)', [(-4.0, 'p = 2, 4: Λ = −∞ (Newman 1976)'), (-1.6, 'even p ≥ 6: finite Λ ≤ 0'), (3.2, 'non-even p > 2: +∞')]),
    (2.4, 'p-balls (1−|u|^p)^a', [(-3.2, 'a ≤ p or p even: Λ < 0 (−8.5 to −0.07)'), (0.0, 'triangle: 0'), (3.2, 'a > p, p not even: +∞')]),
    (1.4, 'Riemann Φ', [(0.1, 'Λ in [0, 0.2]\nRodgers–Tao ≥ 0; Platt–Trudgian ≤ 0.2')]),
]
for y, nm, items in rows:
    ax.text(-4.1, y + 0.35, nm, fontsize=10, weight='bold')
    for x, t in items:
        ax.plot([x], [y], 'o', ms=10, color='tab:blue' if x < 0 else ('tab:orange' if x < 3 else 'tab:red'))
        ax.text(x, y - 0.35, t, fontsize=7.5, ha='center')
ax.text(-2, 0.3, 'real zeros (Λ ≤ 0 side)', ha='center', color='tab:blue')
ax.text(2.2, 0.3, 'not real (Λ > 0)', ha='center', color='tab:red')
ax.set_yticks([]); ax.set_xticks([-4, -2, 0, 2, 4]); ax.set_xticklabels(['−∞', '', '0', '', '+∞'])
ax.set_xlabel('rewind constant Λ (schematic axis)')
ax.set_title('One axis for all three families (Newman’s framework, 1974–76). RH says Φ sits exactly at Λ = 0.\nThe p-ball family jumps from negative to +∞; it never touches 0 except at the triangle.', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_lambda_map.png', dpi=130)
print(len(pts))
