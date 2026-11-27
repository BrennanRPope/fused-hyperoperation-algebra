"""figH: viscosity sweep for the antiparallel pair (base alone and Klein-only seed).
Left: peak spin. Right: Kerr's test sqrt(nu)*Z(t); curves for different nu crossing at one time is the
known non-singular (convergent) signature; a blow-up would make the peak run away as nu shrinks.
usage: python3 fig_nu.py out.png"""
import sys, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
runs = {}
for nu, tag in ((0.004, 'nu0.004'), (0.002, 'n128'), (0.001, 'nu0.001')):
    for c in ('none', 't-g+'):
        try: runs[(nu, c)] = json.load(open('runs/%s_%s_s0.json' % (tag, c)))
        except FileNotFoundError: pass
shade = {0.004: 0.35, 0.002: 0.65, 0.001: 1.0}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for (nu, c), r in sorted(runs.items()):
    col = {0.004: '#f5b041', 0.002: '#e67e22', 0.001: '#922b21'}[nu]
    ls = '-' if c == 't-g+' else ':'
    t = np.array(r['t']); m = t <= 2.6
    lab = ('with Klein-only seed' if c == 't-g+' else 'pair alone (dotted)') + ', ν=%g' % nu
    ax[0].plot(t[m], np.array(r['wmax'])[m], ls, color=col, lw=2, label=lab)
    ax[1].plot(t[m], np.sqrt(nu)*np.array(r['Z'])[m], ls, color=col, lw=2, label=lab)
try:
    r = json.load(open('runs/res192_none_s0.json')); t = np.array(r['t'])
    ax[0].plot(t, r['wmax'], color='#1b4f72', lw=2.5, label='pair alone, ν=0.001, finer grid 192³')
    ax[1].plot(t, np.sqrt(0.001)*np.array(r['Z']), color='#1b4f72', lw=2.5, label='finer grid')
except FileNotFoundError: pass
ax[0].set_xlabel('time'); ax[0].set_ylabel('peak spin on the grid  max|ω|'); ax[0].set_title('Peak spin: not settled at the lowest viscosity')
ax[1].axvline(1.47, color='#888', lw=1, ls='--'); ax[1].text(1.5, 0.08, 'all cross here', color='#555', fontsize=9)
ax[1].set_xlabel('time'); ax[1].set_ylabel('√ν × enstrophy'); ax[1].set_title("Kerr's crossing test")
ax[0].legend(frameon=False, fontsize=8)
for a in ax: a.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=140)
