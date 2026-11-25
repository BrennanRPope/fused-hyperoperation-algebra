import json, math, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BLUE, ORANGE, AQUA, VIOLET, INK, MUTE = '#257de0', '#f47230', '#00b894', '#6a4ce2', '#1d2330', '#8a93a3'
R = json.load(open('fillspill/rungs_speed.json')); B = json.load(open('fillspill/border.json'))
O = json.load(open('fillspill/border_orbit.json')); H = json.load(open('fillspill/heat.json'))
fig, ax = plt.subplots(2, 2, figsize=(10, 7.2), dpi=200)
def style(a, title):
    a.set_title(title, fontsize=10, color=INK, loc='left')
    for s in ('top', 'right'): a.spines[s].set_visible(False)
    for s in ('left', 'bottom'): a.spines[s].set_color(MUTE)
    a.tick_params(colors=MUTE, labelsize=8); a.grid(color='#e6e8ec', lw=0.6)
# 1 rungs
a = ax[0, 0]; S = np.array(R['S'])
a.semilogy(S, R['v1'], color=BLUE, lw=2, label='standard flow')
a.semilogy(S, np.maximum(R['vc'], 1e-3), color=MUTE, lw=1.2, ls=(0, (3, 2)), label='a staircase put in by hand')
for n in (-1, 0, 1): a.axvline(n, color=ORANGE, lw=0.8, alpha=.6)
a.set_xlabel('rung height s', fontsize=8, color=INK); a.set_ylabel('speed of E^s(1) in s', fontsize=8, color=INK)
a.legend(fontsize=7, frameon=False, loc='upper left')
style(a, 'Rungs: the flow never pauses at a rung')
# 2 border orbit
a = ax[0, 1]; arr = np.array(O['a'][:13500]); arr[~np.isfinite(arr)] = np.nan
arr = np.where(np.abs(arr) > 8, np.nan, arr)
a.plot(np.arange(len(arr)), arr, color=BLUE, lw=1.4, label='just past the border: lingers, spills')
k, q, x = 1.0, 0.25 - 1e-2/4, 0.05; b = 0.0; tw = []
for n in range(13500):
    den = 1 - q*b*x; b = (b + x + k*b*x)/den if den != 0 else np.inf; tw.append(b if abs(b) < 8 else np.nan)
a.plot(np.arange(len(tw)), tw, color=ORANGE, lw=1.4, label='just inside two walls: settles')
a.axhline(-2, color=MUTE, lw=0.8, ls=(0, (2, 2))); a.text(13400, -2.3, 'ghost of the double wall', fontsize=7, color=MUTE, ha='right', va='top')
a.set_ylim(-6.5, 6.5); a.set_xlabel('steps of a ↦ a ⊕ x', fontsize=8, color=INK); a.set_ylabel('a', fontsize=8, color=INK)
a.legend(fontsize=7, frameon=False, loc='upper right')
style(a, 'Landscape border: fill, linger, spill')
# 3 linger scaling
a = ax[1, 0]; d = np.array(B['d']); T = np.array(B['T'])
a.loglog(d, T, 'o', color=BLUE, ms=4, label='steps for one turn')
a.loglog(d, T[6] * (d / d[6]) ** -0.5, color=MUTE, lw=1, ls=(0, (3, 2)), label='slope −1/2 (any fold)')
a.set_xlabel('distance past the border, 4q − k²', fontsize=8, color=INK); a.set_ylabel('steps', fontsize=8, color=INK)
a.legend(fontsize=7, frameon=False)
style(a, 'The linger grows without limit at the border')
# 4 latent heat
a = ax[1, 1]; Q = np.array(H['Q']); L = np.array(H['L'])
a.plot(Q - 1, L, color=VIOLET, lw=2, label='latent heat, exact (Baxter)')
for qq, (eo, ed) in H['mc'].items():
    qq = int(qq); a.plot(qq - 1, ed - eo, 'o', color=ORANGE if qq > 4 else 'none', mec=ORANGE, ms=5)
a.plot([], [], 'o', color=ORANGE, ms=5, label='simulation, 64 × 64')
a.plot([], [], 'o', color='none', mec=ORANGE, ms=5, label='simulation, slow settling (no latent heat)')
a.axvline(3, color=MUTE, lw=0.8, ls=(0, (2, 2))); a.text(3.2, 1.12, 'k = 3', fontsize=7, color=MUTE)
a.set_xlim(0.5, 19.5); a.set_xlabel('dial k = colours − 1', fontsize=8, color=INK); a.set_ylabel('per site', fontsize=8, color=INK)
a.legend(fontsize=7, frameon=False, loc='lower right')
style(a, 'Colour counts: a real heat plateau, only past k = 3')
plt.tight_layout()
plt.savefig('figures/f5_fillspill.png'); plt.savefig('page/figures/f5_fillspill.png')
