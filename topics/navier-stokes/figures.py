import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
OUT = 'figures/'
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans'})
BLUE, ORANGE, GREEN, PURPLE, GREY = '#2a6fdb', '#e8833a', '#2e9e6a', '#8a4fd8', '#777777'

# F1: square (2D) vs cube (3D), layers labelled
fig = plt.figure(figsize=(11, 4.6))
ax = fig.add_subplot(1, 2, 1)
sq = np.array([[0,0],[1,0],[1,1],[0,1],[0,0]])
ax.fill(sq[:4,0], sq[:4,1], color=PURPLE, alpha=0.15)
ax.plot(sq[:,0], sq[:,1], color=BLUE, lw=3)
ax.plot(0, 0, 'o', color='k', ms=10)
ax.text(-0.08, -0.12, 'corner: pressure\n(1 number)', ha='right', va='top')
ax.text(0.5, -0.1, 'edges: velocity (2 arrows)', ha='center', va='top', color=BLUE)
ax.text(0.5, 0.5, 'face: vorticity\n(1 number)\nTOP LAYER', ha='center', va='center', color=PURPLE, fontweight='bold')
ax.set_title('2D flow = the square, layers 1, 2, 1\nspin sits on the top layer: it can only be carried\nand smeared, never stretched  →  proved smooth', fontsize=10.5)
ax.set_xlim(-0.6, 1.4); ax.set_ylim(-0.45, 1.25); ax.axis('off'); ax.set_aspect('equal')
ax = fig.add_subplot(1, 2, 2, projection='3d')
V = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
faces = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
ax.add_collection3d(Poly3DCollection([V[f] for f in faces], facecolor=PURPLE, alpha=0.10, edgecolor=GREY))
for e in ([0,1],[0,3],[0,4]): ax.plot(*V[e].T, color=BLUE, lw=3)
ax.scatter(*V[0], color='k', s=60)
ax.text(-0.15, -0.15, -0.15, 'pressure (1)', fontsize=9)
ax.text(0.5, -0.25, 0, 'velocity (3)', color=BLUE, fontsize=9)
ax.text(1.05, 0.5, 0.5, 'vorticity (3)\non the faces', color=PURPLE, fontsize=9)
ax.text(0.2, 0.6, 1.35, 'helicity (1): the whole\nvolume, top layer', color=GREEN, fontsize=9, fontweight='bold')
ax.set_title('3D flow = the cube, layers 1, 3, 3, 1\nspin sits on the faces, one layer below the top:\nfaces can shrink while volume stays fixed  →  open', fontsize=10.5)
ax.set_axis_off(); ax.view_init(22, -55)
plt.tight_layout(); plt.savefig(OUT+'f1_square_vs_cube.png', dpi=150); plt.close()

# F2: material cube squeezed at fixed volume
fig = plt.figure(figsize=(11, 5.0))
for i, L in enumerate([1, 2, 4]):
    ax = fig.add_subplot(1, 3, i+1, projection='3d')
    s = np.array([1/np.sqrt(L), 1/np.sqrt(L), L])
    W = V*s - np.array([s[0]/2, s[1]/2, 0])
    ax.add_collection3d(Poly3DCollection([W[f] for f in faces], facecolor=BLUE, alpha=0.12, edgecolor=GREY))
    top = W[[4,5,6,7]]
    ax.add_collection3d(Poly3DCollection([top], facecolor=PURPLE, alpha=0.55))
    ax.quiver(0, 0, L+0.1, 0, 0, 0.3*L, color=PURPLE, lw=2.5, arrow_length_ratio=0.25)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(0, 6.0); ax.set_box_aspect((1, 1, 2.2))
    ax.set_axis_off(); ax.view_init(18, -60)
    ax.set_title(f'stretch ×{L}\nvolume 1, face area 1/{L}\nspin ×{L}', fontsize=10.5)
fig.text(0.5, 0.03, 'Volume is locked at 1 and the spin flux through the purple face is locked too, so when the face shrinks the spin must grow.\nBlow-up would mean a face shrinking to zero in finite time.', ha='center', fontsize=10.5)
plt.tight_layout(rect=(0, 0.1, 1, 0.94)); plt.savefig(OUT+'f2_stretching_cube.png', dpi=150); plt.close()

# F3: right- and left-handed helical waves (curl eigenmodes)
fig = plt.figure(figsize=(11, 4.6))
zz = np.linspace(0, 4*np.pi, 33)
for i, (sgn, name, col) in enumerate([(1, 'right-handed: curl u = +k u', ORANGE), (-1, 'left-handed: curl u = −k u', GREEN)]):
    ax = fig.add_subplot(1, 2, i+1, projection='3d')
    ux, uy = np.sin(zz), sgn*np.cos(zz)
    ax.quiver(0*zz, 0*zz, zz, ux, uy, 0*zz, color=col, lw=1.6, arrow_length_ratio=0.2)
    zf = np.linspace(0, 4*np.pi, 400); ax.plot(np.sin(zf), sgn*np.cos(zf), zf, color='k', lw=0.8, alpha=0.6)
    ax.plot([0,0],[0,0],[0,4*np.pi], color=GREY, lw=1)
    ax.set_title(name + '\nthe arrow turns a quarter turn per quarter wavelength', fontsize=10.5)
    ax.set_axis_off(); ax.view_init(15, -50); ax.set_box_aspect((1,1,2.2))
fig.suptitle('Every incompressible flow splits into these two helices (Waleffe 1992). A pure single helix is an exact Beltrami flow:\nits velocity and spin point the same way, so the troublesome cross term vanishes.', fontsize=10.5)
plt.tight_layout(); plt.savefig(OUT+'f3_helical_waves.png', dpi=150); plt.close()

# F4: draining bathtub acoustic black hole
A, B, cs = 1.0, 2.0, 1.0
fig, ax = plt.subplots(figsize=(6.8, 6.8))
th = np.linspace(0, 2*np.pi, 400)
ax.fill(np.hypot(A,B)*np.cos(th), np.hypot(A,B)*np.sin(th), color=ORANGE, alpha=0.12)
ax.fill(A*np.cos(th), A*np.sin(th), color='k', alpha=0.85)
ax.plot(np.hypot(A,B)*np.cos(th), np.hypot(A,B)*np.sin(th), color=ORANGE, lw=2)
gx, gy = np.meshgrid(np.linspace(-4,4,23), np.linspace(-4,4,23)); R2 = gx**2+gy**2; R2[R2 < 0.3] = np.nan
ax.streamplot(np.linspace(-4,4,23), np.linspace(-4,4,23), (-A*gx - B*gy)/R2, (-A*gy + B*gx)/R2, color=BLUE, density=1.1, linewidth=0.7, arrowsize=0.8)
# sound fronts emitted at several points: circles of radius cs*dt carried by flow
for r0, ang in [(3.2, 0.3), (2.6, 2.2), (1.7, 4.0), (1.25, 5.3)]:
    px, py = r0*np.cos(ang), r0*np.sin(ang)
    vx, vy = (-A*px - B*py)/r0**2, (-A*py + B*px)/r0**2
    dt = 0.35
    ax.plot(px + vx*dt + cs*dt*np.cos(th), py + vy*dt + cs*dt*np.sin(th), color=PURPLE, lw=1.4)
    ax.plot(px, py, 'o', color=PURPLE, ms=4)
ax.text(0, 0, 'acoustic\nhorizon', color='w', ha='center', va='center', fontsize=9)
ax.text(1.0, 2.05, 'ergosphere edge:\nswirl alone outruns sound', color=ORANGE, fontsize=9)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Sound in a draining, spinning bathtub (Unruh, Visser)\nblue: the flow; purple: a sound pulse a moment after it is made.\nThe spin of the medium makes the ergosphere, like a rotating black hole.', fontsize=10.5)
plt.tight_layout(); plt.savefig(OUT+'f4_acoustic_black_hole.png', dpi=150); plt.close()
print('ok')
