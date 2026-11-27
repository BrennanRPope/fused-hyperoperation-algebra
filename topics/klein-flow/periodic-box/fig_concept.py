import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# ---------- figure A: known dial paths to the periodic box ----------
fig, ax = plt.subplots(figsize=(11, 6.2))
paths = [
    ("Box height  (thin slab -> cube)", "thin slab: smooth for large flows\n(Raugel & Sell 1993; Montgomery-Smith 1999)", 0.30),
    ("Dissipation power  (|k|^2a, a: 5/4 -> 1)", "a >= 5/4: smooth (Lions 1969);\nlog below it (Tao 2009; Barbato-Morandin-Romito 2014)", 0.18),
    ("Rotation rate  (fast -> none)", "fast rotation: smooth\n(Babin, Mahalov & Nicolaenko 1997-99)", 0.25),
    ("Cross-hand weight  (0 -> 1)", "one hand only (weight 0): smooth\n(Biferale & Titi 2013); proved only at the start", 0.02),
    ("Dimension by decimation  (D: 2 -> 3)", "D <= 2 behaves 2D; simulations only\n(Frisch, Pomyalov, Procaccia & Ray 2012)", 0.20),
]
for i, (name, note, stall) in enumerate(paths):
    y = len(paths) - i
    ax.plot([0, 1], [y, y], color='0.85', lw=8, solid_capstyle='round', zorder=1)
    sim = 'simulations' in note
    ax.plot([0, stall], [y, y], color='tab:blue' if sim else 'tab:green', lw=8, solid_capstyle='round', zorder=2)
    ax.plot([stall, 1], [y, y], color='tab:red', lw=2, ls=':', zorder=2)
    ax.text(-0.02, y + 0.22, name, ha='left', fontsize=10.5, weight='bold')
    ax.text(0.0, y - 0.38, note, ha='left', fontsize=8.5, color='0.25')
    ax.text(max(stall, 0.42) + 0.01, y - 0.30, 'open beyond here for strong flows' if not sim else 'no proof anywhere on this dial', fontsize=7.5, color='tab:red', va='center')
ax.scatter([1]*len(paths), [len(paths)-i for i in range(len(paths))], s=160, color='black', zorder=3)
ax.text(1.02, 3, "Clay's\nperiodic box\n(3D, full\nstrength)", fontsize=11, va='center', weight='bold')
ax.text(0.0, 0.25, "green: proved smooth for strong flows   blue: seen in simulations only   dotted red: open   black dot: the Clay setting (options B and D)   bar lengths are schematic",
        fontsize=8.5, color='0.3')
ax.set_xlim(-0.05, 1.2); ax.set_ylim(0, len(paths) + 0.7); ax.axis('off')
ax.set_title("Figure A. Five known dials that start where Navier–Stokes is proved smooth and end at Clay's box.\n"
             "Each gets part of the way for strong flows; none reaches the end. The twist is not on this list: it is a pick (Figure B).",
             fontsize=10.5, loc='left')
fig.tight_layout(); fig.savefig('figures/figA-dial-paths-to-the-periodic-box.png', dpi=130); plt.close(fig)

# ---------- figure B: the twist unfolds into the box; tubes come in antiparallel pairs ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
ax = axs[0]
ax.add_patch(plt.Rectangle((0, 0), np.pi, 2*np.pi, fc='tab:blue', alpha=0.15, ec='tab:blue'))
ax.add_patch(plt.Rectangle((np.pi, 0), np.pi, 2*np.pi, fc='tab:orange', alpha=0.15, ec='tab:orange'))
ax.text(np.pi/2, 2*np.pi + 0.25, 'Klein bottle box\n(one copy)', ha='center', fontsize=9.5, color='tab:blue')
ax.text(3*np.pi/2, 2*np.pi + 0.25, 'its mirror image\n(shifted by half, flipped)', ha='center', fontsize=9.5, color='tab:orange')
ax.annotate('', xy=(1.0 + np.pi - 0.3, 2*np.pi - 1.0 - 0.3), xytext=(1.3, 1.35),
            arrowprops=dict(arrowstyle='->', color='0.3', lw=1.5, connectionstyle='arc3,rad=0.25'))
ax.text(0.25, 3.4, 'glide:\n(x, y) -> (x + pi, -y)', fontsize=9, color='0.3')
# a z-tube (out of page) and its image
x0, y0 = 1.0, 1.0
ax.add_patch(plt.Circle((x0, y0), 0.35, fc='tab:red', alpha=0.8)); ax.text(x0, y0, 'out', ha='center', va='center', fontsize=7, color='w')
ax.add_patch(plt.Circle((x0 + np.pi, 2*np.pi - y0), 0.35, fc='tab:blue', alpha=0.8)); ax.text(x0+np.pi, 2*np.pi-y0, 'in', ha='center', va='center', fontsize=7, color='w')
ax.set_xlim(-0.3, 2*np.pi + 0.3); ax.set_ylim(-0.3, 2*np.pi + 1.1); ax.set_aspect('equal')
ax.set_xticks([0, np.pi, 2*np.pi]); ax.set_xticklabels(['0', 'pi', '2pi']); ax.set_yticks([0, np.pi, 2*np.pi]); ax.set_yticklabels(['0', 'pi', '2pi'])
ax.set_title('(a) Two copies make the ordinary periodic box (x, y), times\nthe vertical circle z. A vortex tube pointing out of the page\nreturns from the flip pointing in: an antiparallel partner.', fontsize=9.5, loc='left')
ax = axs[1]
xs = np.linspace(0, 2*np.pi, 300)
for m, ls, col in ((1, '-', 'tab:purple'), (2, '-', 'tab:green')):
    off = 0 if m == 1 else 2.6
    yA = 0.6 + 0.25*np.cos(m*xs)          # tube A, running along x
    yB = -(0.6 + 0.25*np.cos(m*(xs - np.pi)))   # its image: shifted by pi, mirrored
    ax.plot(xs, yA + off, color=col, lw=3); ax.plot(xs, yB + off, color=col, lw=3, ls='--')
    lab = ('odd wavenumber (1): both tubes bend the same way,\nthe gap never closes: the "antisymmetric" Crow pattern' if m == 1 else
           'even wavenumber (2): bends meet,\nthe "symmetric" Crow pattern (reconnection studies)')
    ax.text(2*np.pi + 0.15, off, lab, fontsize=8.5, va='center', color=col)
ax.axhline(0, color='0.8', lw=0.8); ax.axhline(2.6, color='0.8', lw=0.8)
ax.set_xlim(0, 2*np.pi + 4.2); ax.set_ylim(-1.3, 3.6); ax.axis('off')
ax.set_title('(b) Tubes lying along the glide direction (solid) and their forced\npartners (dashed, opposite spin). Which classic pair you get\ndepends only on whether the wiggle has odd or even wavenumber.', fontsize=9.5, loc='left')
fig.suptitle('Figure B. The twist comes off by unfolding: every Klein bottle x circle flow is a periodic-box flow with one extra symmetry.', fontsize=10.5, x=0.01, ha='left')
fig.tight_layout(); fig.savefig('figures/figB-twist-unfolds-into-the-box.png', dpi=130); plt.close(fig)
print('ok')
