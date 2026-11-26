import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
sig, Rb = 3.0, 1.0
def f(r, s=sig, R=Rb): return (np.tanh(s*(r+R)) - np.tanh(s*(r-R)))/(2*np.tanh(s*R))
def fp(r, s=sig, R=Rb): return s*(1/np.cosh(s*(r+R))**2 - 1/np.cosh(s*(r-R))**2)/(2*np.tanh(s*R))
g = np.linspace(-2.2, 2.2, 400); X, Y = np.meshgrid(g, g); Rr = np.sqrt(X**2+Y**2)+1e-9
OUT = 'figures/'

# Fig 1: the bubble is a river
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
im = ax[0].imshow(f(Rr), extent=[-2.2, 2.2, -2.2, 2.2], origin='lower', cmap='Blues', vmin=0, vmax=1.1)
q = slice(None, None, 25)
ax[0].quiver(X[q, q], Y[q, q], f(Rr)[q, q], 0*X[q, q], color='navy', scale=12, width=0.004)
ax[0].add_patch(Circle((0, 0), 0.25, color='orange')); ax[0].text(0, -0.5, 'ship (at rest in\nits own patch)', ha='center', fontsize=8)
ax[0].text(0, 1.75, 'outside: space still', ha='center'); ax[0].text(1.25, 1.05, 'wall', color='crimson')
ax[0].add_patch(Circle((0, 0), Rb, fill=False, ls='--', color='crimson'))
ax[0].set_title('(a) Alcubierre bubble = a river of space\nmoving at v_s inside, still outside'); ax[0].set_xlabel('x (direction of travel)'); ax[0].set_ylabel('y')
plt.colorbar(im, ax=ax[0], label='river speed / v_s  (= f)')
r1 = np.linspace(-2.5, 2.5, 500)
ax[1].plot(r1, f(np.abs(r1)), 'b', lw=2, label='f: river speed / v_s')
ax[1].plot(r1, fp(np.abs(r1))*np.sign(r1)/sig, 'crimson', lw=1.5, label="f' / sigma: wall steepness")
ax[1].axvspan(-1.25, -0.75, color='crimson', alpha=0.08); ax[1].axvspan(0.75, 1.25, color='crimson', alpha=0.08)
ax[1].text(0, 1.05, 'carried along', ha='center'); ax[1].text(1.9, 0.12, 'still', ha='center'); ax[1].text(-1.9, 0.12, 'still', ha='center')
ax[1].set_xlabel('position along x through the centre'); ax[1].legend(fontsize=8, loc='lower center'); ax[1].set_title('(b) shape function f(r): flat inside, steep wall')
plt.tight_layout(); plt.savefig(OUT+'fig1_bubble_is_a_river.png', dpi=130); plt.close()

# Fig 2: energy = cross term of strain rates
vs = 1.0
gx = vs*fp(Rr)*X/Rr; gy = vs*fp(Rr)*Y/Rr
e2 = -(gy**2)/4; theta = gx
fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
im = ax[0].imshow(theta, extent=[-2.2, 2.2, -2.2, 2.2], origin='lower', cmap='RdBu_r', vmin=-1.5, vmax=1.5)
ax[0].text(1.55, 0, 'squeeze\n(ahead)', ha='center', fontsize=8); ax[0].text(-1.55, 0, 'stretch\n(behind)', ha='center', fontsize=8)
ax[0].set_title("(a) volume change e1 = expansion\n(Alcubierre's famous picture)"); plt.colorbar(im, ax=ax[0], fraction=0.046)
im = ax[1].imshow(e2/(8*np.pi), extent=[-2.2, 2.2, -2.2, 2.2], origin='lower', cmap='Reds_r')
ax[1].set_title('(b) energy density = e2 / 8pi\n(cross term of strain rates): all negative'); plt.colorbar(im, ax=ax[1], fraction=0.046)
ax[1].text(0, 1.55, 'ring of negative energy\n(above and below the path)', ha='center', fontsize=8)
for a in ax[:2]: a.set_xlabel('x (travel)'); a.set_ylabel('y')
# panel c: square stretched one way, squeezed the other
a = ax[2]; a.set_aspect('equal'); a.axis('off')
sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])*0.6
a.add_patch(Polygon(sq, fill=False, lw=1.5))
M = np.array([[1, 0.45], [0.45, 1]]); w_, V_ = np.linalg.eigh(M)
sq2 = sq @ M.T; a.add_patch(Polygon(sq2, fill=False, lw=2, color='crimson'))
for lam, vec in zip(w_, V_.T):
    col = 'green' if lam > 1 else 'purple'
    a.annotate('', xy=1.3*vec, xytext=-1.3*vec, arrowprops=dict(arrowstyle='<->', color=col, lw=2))
a.text(0.95, 1.05, 'stretch (+)', color='green'); a.text(0.95, -1.15, 'squeeze (-)', color='purple')
a.text(-1.5, -1.65, 'Sideways shear of the river = one stretch + one squeeze.\nProduct of the pair is negative, so the cross term\ne2 = l1 l2 + l2 l3 + l3 l1 < 0  ->  negative energy.', fontsize=8)
a.set_xlim(-1.6, 1.8); a.set_ylim(-1.8, 1.6); a.set_title('(c) why: a square of space in the wall')
plt.tight_layout(); plt.savefig(OUT+'fig2_energy_is_the_cross_term.png', dpi=130); plt.close()

# Fig 3: two composition rules and the horizon
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
r1 = np.linspace(0, 2.2, 400)
for v_, c_ in [(0.5, 'tab:blue'), (1.0, 'tab:orange'), (2.0, 'tab:red')]:
    ax[0].plot(r1, v_*f(r1), color=c_, lw=2, label=f'v_s = {v_}c')
ax[0].axhline(1, color='k', ls='--'); ax[0].text(1.5, 1.05, 'light speed: the wall of the velocity rule', fontsize=8)
ax[0].annotate('horizon: river\ncrosses light speed', xy=(1.06, 1.0), xytext=(1.3, 1.5), arrowprops=dict(arrowstyle='->'), fontsize=8)
ax[0].set_xlabel('distance from bubble centre'); ax[0].set_ylabel('river speed (units of c)'); ax[0].legend(fontsize=8)
ax[0].set_title('(a) a faster-than-light bubble makes horizons\n(same rule that makes r = 2M a horizon)')
a = ax[1]; a.axis('off'); a.set_xlim(0, 10); a.set_ylim(0, 6)
a.annotate('', xy=(0.3, 3.6), xytext=(0.3, 1.0), arrowprops=dict(arrowstyle='->')); a.text(0.45, 3.4, 'time', fontsize=8)
for xc, v_, lab in [(1.6, 0.0, 'outside\n(river 0)'), (4.6, 0.6, 'wall (river 0.6c)'), (7.4, 2.0, 'inside, river 2c:\nboth edges forward')]:
    a.plot([xc, xc+0.7*(v_+1)], [1, 3.2], 'b'); a.plot([xc, xc+0.7*(v_-1)], [1, 3.2], 'b'); a.plot([xc, xc], [1, 3.2], color='gray', ls=':')
    a.text(xc, 0.3, lab, ha='center', fontsize=8)
a.text(5, 5.3, 'Light cones tip over with the river: dx/dt = v + 1 and v - 1.', ha='center', fontsize=9)
a.text(5, 4.5, 'Across the river speeds add plainly (dial k = 0, no wall);\nat one point they compose by the relativity rule (k = -1, wall at c).', ha='center', fontsize=8)
a.set_title('(b) where the theory puts the trick')
plt.tight_layout(); plt.savefig(OUT+'fig3_two_rules_and_the_horizon.png', dpi=130); plt.close()

# Fig 4: energy budget
KG = 2.998e8**2/6.674e-11; Lp = 1.616e-35
thick = np.logspace(np.log10(100*Lp), 1, 300)
fig, ax = plt.subplots(figsize=(8, 4.8))
for Rm, c_ in [(10, 'tab:blue'), (100, 'tab:red')]:
    M = (Rm**2*(2/thick)/36)*KG
    ax.loglog(thick, M, color=c_, lw=2, label=f'bubble radius {Rm} m, v_s = c')
for val, lab in [(1.9e27, 'Jupiter'), (2e30, 'Sun'), (1.5e53, 'all visible matter (order)')]:
    ax.axhline(val, color='gray', ls=':'); ax.text(2e-33, val*1.6, lab, fontsize=8, color='gray')
ax.axvline(100*Lp, color='k', ls='--'); ax.text(1.5*100*Lp, 1e30, 'wall allowed by\nquantum inequalities\n(Pfenning-Ford)', fontsize=8)
ax.set_xlabel('wall thickness (m)'); ax.set_ylabel('negative energy needed, as mass (kg)')
ax.set_title('Fig 4. Alcubierre energy budget: E = -(c^2/G) v^2 R^2 / (18 x wall)  (scales with v squared)')
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig(OUT+'fig4_energy_budget.png', dpi=130); plt.close()

# Fig 5: irrotational warp, positive and negative cancel
n = 161; g1 = np.linspace(-2.2, 2.2, n); h = g1[1]-g1[0]
X3, Y3, Z3 = np.meshgrid(g1, g1, g1, indexing='ij'); R3 = np.sqrt(X3**2+Y3**2+Z3**2)+1e-9
phi = X3*f(R3)
U = np.gradient(phi, h, h, h); Hs = [np.gradient(u, h, h, h) for u in U]
tr = Hs[0][0]+Hs[1][1]+Hs[2][2]; trK2 = sum(Hs[i][j]*Hs[j][i] for i in range(3) for j in range(3))
rho2 = ((tr**2 - trK2)/(16*np.pi))[:, :, n//2]   # z = 0 slice of the full 3D energy density
Xg, Yg = X3[:, :, n//2], Y3[:, :, n//2]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
lim = np.percentile(np.abs(rho2), 99)
im = ax[0].imshow(rho2.T, extent=[-2.2, 2.2, -2.2, 2.2], origin='lower', cmap='RdBu', vmin=-lim, vmax=lim)
ax[0].set_title('(a) swirl-free bubble (river = slope of a hill):\nred = negative, blue = positive energy'); plt.colorbar(im, ax=ax[0], fraction=0.046)
ax[0].set_xlabel('x (travel)'); ax[0].set_ylabel('y')
import json; Rj = json.load(open('results.json'))
A5, A6 = Rj['A5_values'], Rj['A6_values']
ax[1].bar(['Alcubierre\n(swirling)', 'swirl-free:\nnegative part', 'swirl-free:\npositive part', 'swirl-free:\ntotal'],
          [A5['E_closed'], A6['Eneg'], A6['Epos'], A6['Etot']], color=['crimson', 'crimson', 'steelblue', 'gray'])
ax[1].axhline(0, color='k', lw=0.8)
ax[1].set_title('(b) total = -(1/32pi) x (total swirl squared)\nno swirl: totals cancel, but negative patches stay'); ax[1].set_ylabel('energy (geometric units, R = 1, v_s = c)')
plt.tight_layout(); plt.savefig(OUT+'fig5_swirl_free_bubble.png', dpi=130); plt.close()

# Fig 6: artificial gravity
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
zz = np.linspace(0, 0.4, 200)
hz = 0.45 - zz
ax[0].plot(np.sqrt(2*hz), zz, 'b', lw=2, label='w^2 linear in height: vacuum (flat, no tides)')
ax[0].plot(np.sqrt(2*hz + 3*hz**2), zz, 'r--', lw=1.5, label='w^2 bent: needs sideways stress')
ax[0].fill_between([0, 1.2], -0.03, 0, color='k', alpha=0.4); ax[0].text(0.6, -0.018, 'floor (must be pushed up at g)', ha='center', va='center', fontsize=8, color='w')
ax[0].set_xlabel('river speed falling toward the floor'); ax[0].set_ylabel('height'); ax[0].legend(fontsize=7, loc='upper right')
ax[0].set_title('(a) uniform gravity is free, but it is\njust an accelerating floor (equivalence principle)')
L = np.logspace(0, 4, 100)
ax[1].loglog(L, 9.81*L**2/6.674e-11, 'k', lw=2, label='ball of radius L')
ax[1].loglog(L, 9.81/(2*np.pi*6.674e-11)*L**2, 'g--', lw=2, label='plate of side L')
for val, lab in [(1e15, 'Mt Everest ~1e15 kg'), (7.3e22, 'the Moon')]:
    ax[1].axhline(val, color='gray', ls=':'); ax[1].text(1.2, val*1.8, lab, fontsize=8, color='gray')
ax[1].set_xlabel('size L (m)'); ax[1].set_ylabel('mass for 1 g at the surface (kg)'); ax[1].legend(fontsize=8)
ax[1].set_title('(b) real (tidal) gravity needs real mass')
rpm = np.linspace(0.5, 10, 200); r_ = 9.81/(rpm*2*np.pi/60)**2
ax[2].semilogy(rpm, r_, 'purple', lw=2); ax[2].axvspan(0.5, 2, color='green', alpha=0.1); ax[2].text(0.7, 30, 'comfortable\n(< 2 rpm)', fontsize=8)
ax[2].scatter([2, 4], [224, 56], color='k'); ax[2].text(2.2, 260, '224 m', fontsize=8); ax[2].text(4.2, 62, '56 m', fontsize=8)
ax[2].set_xlabel('spin rate (rpm)'); ax[2].set_ylabel('radius for 1 g (m)')
ax[2].set_title('(c) spin: the swap rule at d = +1 gives\nthe Coriolis side push (precession page)')
plt.tight_layout(); plt.savefig(OUT+'fig6_artificial_gravity.png', dpi=130); plt.close()

# Fig 7: wormhole throat
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
r0 = 1.0; rr = np.linspace(r0, 4, 200); zemb = 2*np.sqrt(r0*(rr-r0))   # Ellis/M-T b = r0^2/r embedding: z = r0 arccosh(r/r0)
zemb = r0*np.arccosh(rr/r0)
for s_ in (1, -1):
    ax[0].plot(rr, s_*zemb, 'k'); ax[0].plot(-rr, s_*zemb, 'k')
ax[0].annotate('throat: smallest sphere', xy=(1, 0), xytext=(1.6, -0.4), arrowprops=dict(arrowstyle='->'), fontsize=8)
ax[0].set_title('(a) a wormhole throat (embedding picture)'); ax[0].set_xlabel('r'); ax[0].set_ylabel('depth')
lam = np.linspace(-2, 2, 200); A = 1 + lam**2
ax[1].plot(lam, A, 'b', lw=2, label='bundle area 1 + (a o b)')
ax[1].plot(lam, 1 - 0.5*lam**2, 'r--', label='normal matter: area only focuses')
ax[1].axvline(0, color='gray', ls=':'); ax[1].text(0.1, 2.5, 'area minimum\n= throat', fontsize=8)
ax[1].set_xlabel('distance along the light rays'); ax[1].legend(fontsize=8)
ax[1].set_title('(b) light through a throat must spread again:\nfocusing rule then needs negative energy')
b_ = r0**2/rr; bp = -r0**2/rr**2
ax[2].plot(rr, (rr*bp - b_)/(8*np.pi*rr**3), 'crimson', lw=2, label='rho + p_r')
ax[2].plot(rr, bp/(8*np.pi*rr**2), 'b', lw=1, label='rho')
ax[2].axhline(0, color='k', lw=0.8); ax[2].set_xlabel('r'); ax[2].legend(fontsize=8)
ax[2].set_title('(c) Morris-Thorne (b = r0^2/r): null energy\nnegative at and near the throat')
plt.tight_layout(); plt.savefig(OUT+'fig7_wormhole_throat.png', dpi=130); plt.close()

# Fig 8: link probability through the thermofield-double wormhole
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
tt = np.linspace(0, 1.6, 300)
for ell, c_ in [(0.5, 'tab:blue'), (1.0, 'tab:green'), (2.0, 'tab:red')]:
    s = np.sinh(np.pi*ell); ch = np.cosh(2*np.pi*tt)
    ax[0].plot(tt, s**2/(s**2+ch**2), color=c_, lw=2, label=f'region width {ell} (in units of beta)')
ax[0].axhline(0.5, color='k', ls='--'); ax[0].text(1.0, 0.53, 'above 1/2: wormhole links A and B\n(mutual information > 0)', fontsize=8)
ax[0].set_xlabel('time t on both sides (units of beta)'); ax[0].set_ylabel('link probability x = two-sided cross-ratio')
ax[0].set_title('(a) ER=EPR in our reading: chance a random\ngeodesic links A (left) to B (right) through the throat'); ax[0].legend(fontsize=8)
C9 = Rj['C9']
ax[1].plot(-np.array(C9['a'])/abs(C9['a_crit']), C9['eta'], 'o-', color='purple')
ax[1].axhline(1, color='k', ls=':'); ax[1].text(-0.55, 1.005, 'link probability 1 = light-connected = traversable', fontsize=8)
ax[1].set_xlabel('<- ordinary matter (longer throat)      negative energy (shorter throat) ->\nshift, as a fraction of the shift that opens the throat'); ax[1].set_ylabel('link probability (one light-cone direction)')
ax[1].set_title('(b) Gao-Jafferis-Wall in this reading: negative\nenergy pushes the link probability up to 1')
plt.tight_layout(); plt.savefig(OUT+'fig8_link_probability_through_throat.png', dpi=130); plt.close()

# Fig 9: one wall
fig, ax = plt.subplots(figsize=(11, 4.6)); ax.axis('off'); ax.set_xlim(0, 11); ax.set_ylim(0, 5)
boxes = [(0.2, 'Warp bubble', 'energy = cross term e2\nof the river\'s strain rates;\nsideways shear makes\ne2 < 0', 'crimson'),
         (2.9, 'Wormhole', 'throat = light bundle area\n1 + (a o b) at a minimum;\nfocusing rule then makes\nnull energy < 0', 'crimson'),
         (5.6, 'Artificial gravity', 'uniform field = free but it\nis an accelerating floor;\nreal tides need real mass\n(or spin, or thrust)', 'gray'),
         (8.3, 'ER = EPR', 'two-sided cross-ratio is a\nlink probability; open iff > 1/2;\nnegative energy lifts it to 1\n(traversable)', 'steelblue')]
for x0, title, body, col in boxes:
    ax.add_patch(plt.Rectangle((x0, 1.4), 2.5, 2.9, fill=False, ec=col, lw=2))
    ax.text(x0+1.25, 3.95, title, ha='center', fontsize=11, weight='bold', color=col)
    ax.text(x0+1.25, 2.6, body, ha='center', va='center', fontsize=8.5)
ax.text(5.5, 0.7, 'Same wall in both faster-than-light routes: the corner (cross term) of the composition has to turn negative.\n'
        'The theory makes that visible; it does not lower the wall.', ha='center', fontsize=10)
ax.set_title('Fig 9. What the theory sees in each topic')
plt.savefig(OUT+'fig9_one_wall.png', dpi=130, bbox_inches='tight'); plt.close()
print('ok')
