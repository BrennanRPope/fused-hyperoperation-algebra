import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans'})
d = np.load('tg_helical.npz'); keys = list(d['keys']); P = d['P']; t = d['times']
single = P[:, keys.index('+++')] + P[:, keys.index('---')]; total = P.sum(1)
full = np.load('homochiral_full.npy'); one = np.load('homochiral_single.npy')
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax[0].fill_between(t, 0, single, color='#2e9e6a', alpha=0.8, label='same-handed corners (++, −−)')
ax[0].fill_between(t, single, total, color='#e8833a', alpha=0.7, label='mixed-handed corners (+−)')
ax[0].set_xlabel('time'); ax[0].set_ylabel('vortex stretching (enstrophy production)')
ax[0].set_title('3D Taylor–Green vortex, 64³, ν = 1/400\nabout 90% of the stretching comes from + meeting −', fontsize=10.5)
ax[0].legend(frameon=False, fontsize=9.5, loc='upper left'); ax[0].spines[['top','right']].set_visible(False)
ax[1].plot(full[:,0], full[:,2], color='#e8833a', lw=2.2, label='full Navier–Stokes')
ax[1].plot(one[:,0], one[:,2], color='#2e9e6a', lw=2.2, label='right-handed waves only')
ax[1].set_xlabel('time'); ax[1].set_ylabel('largest spin anywhere, max |ω|')
ax[1].set_title('Same right-handed start, 48³, ν = 1/300\none handedness alone never stretches spin up', fontsize=10.5)
ax[1].legend(frameon=False, fontsize=9.5); ax[1].set_ylim(0, None); ax[1].spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('figures/f5_handedness.png', dpi=150)
print('peak single share', (single/np.where(total>1e-9,total,np.nan))[2:].round(3))
