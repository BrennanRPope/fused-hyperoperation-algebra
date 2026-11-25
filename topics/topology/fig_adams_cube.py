"""Picture of the Adams-Atiyah proof on the mult-add cube (CP^2 case)."""
import numpy as np, matplotlib, itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B='#1d2bd8'; R='#cf3519'; K='#15131a'; M='#b9b7c2'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'savefig.dpi':160,'savefig.bbox':'tight','savefig.facecolor':'white'})
fig,axs=plt.subplots(1,3,figsize=(15,5.4))
col={0:K,1:B,2:R,3:M}
def dot(ax,p,g,lab=None):
    ax.scatter(*p,s=230 if g else 120,color=col[g],zorder=4,edgecolor='white',lw=1.5)
    if lab: ax.text(p[0],p[1]-0.17,lab,ha='center',va='top',fontsize=10,color=col[g])
# panel 1: square = psi^2
ax=axs[0]; P={(0,0):(0,0),(1,0):(1,0),(0,1):(0,1),(1,1):(1,1)}
for a,b in [((0,0),(1,0)),((0,0),(0,1)),((1,0),(1,1)),((0,1),(1,1))]:
    ax.plot(*zip(P[a],P[b]),color=M,lw=2,zorder=1)
for c,p in P.items(): dot(ax,p,sum(c),{0:'1',1:'x',2:'x²'}[sum(c)])
ax.set_title('ψ²(x) = x ⊕ x = 2x + x²\nsquare: 2 edge corners, 1 diagonal corner',fontsize=12)
ax.text(0.5,-0.55,'mod 2 the two x corners cancel:\nψ²(x) ≡ x², so μ₂ = 1 is odd',ha='center',fontsize=11,color=R)
ax.set_xlim(-0.5,1.5); ax.set_ylim(-0.9,1.4); ax.set_aspect('equal'); ax.axis('off')
# panel 2: cube = psi^3
ax=axs[1]
def proj(c): return (c[0]+0.45*c[2], c[1]+0.32*c[2])
corners=list(itertools.product([0,1],repeat=3))
for a in corners:
    for b in corners:
        if sum(abs(i-j) for i,j in zip(a,b))==1 and a<b:
            ax.plot(*zip(proj(a),proj(b)),color=M,lw=2,zorder=1)
for c in corners:
    g=sum(c); dot(ax,proj(c),g if g<3 else 3)
ax.text(*np.add(proj((1,1,1)),(0.12,0.08)),'x³ = 0\n(no 6-cell in CP²)',fontsize=10,color='#7d7a86')
ax.set_title('ψ³(x) = x ⊕ x ⊕ x = 3x + 3x²\ncube: 3 edge corners, 3 face corners',fontsize=12)
ax.text(0.72,-0.55,'μ₃ = 3',ha='center',fontsize=11,color=R)
ax.set_xlim(-0.5,1.9); ax.set_ylim(-0.9,1.6); ax.set_aspect('equal'); ax.axis('off')
# panel 3: psi^2 psi^3 = psi^6 as three squares
ax=axs[2]
centers=[(0,0),(2.2,0),(1.1,1.9)]
blues=[]
for cx,cy in centers:
    sq={(0,0):(cx,cy),(1,0):(cx+0.8,cy),(0,1):(cx,cy+0.8),(1,1):(cx+0.8,cy+0.8)}
    for a,b in [((0,0),(1,0)),((0,0),(0,1)),((1,0),(1,1)),((0,1),(1,1))]:
        ax.plot(*zip(sq[a],sq[b]),color=M,lw=1.6,zorder=1)
    for c,p in sq.items():
        if sum(c)>0: dot(ax,p,sum(c))
    blues+= [sq[(1,0)],sq[(0,1)]]
for i,j in itertools.combinations(range(6),2):
    if i//2!=j//2: ax.plot(*zip(blues[i],blues[j]),color=R,lw=0.7,alpha=0.45,zorder=0)
ax.set_title('ψ²ψ³ = ψ⁶: a cube of squares\n6 x corners; x² corners 3 inside + 12 across = 15 = C(6,2)',fontsize=12)
ax.text(1.5,-0.75,'any order gives the same 6-cube, so ψ² and ψ³ commute',ha='center',fontsize=11,color=K)
ax.set_xlim(-0.5,3.5); ax.set_ylim(-1.1,3.0); ax.set_aspect('equal'); ax.axis('off')
plt.savefig('figures/t7_adams_cube.png'); print('ok')
