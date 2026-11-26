import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
BLUE='#2f6fb3'; ORANGE='#e08a2e'; GREY='#9aa0a6'

# Fig 1: coverage grid
names=['1. Affine Composition Theory','2. Conic Composition Theory','3. Affine Sphere Composition Theory','4. Unit Ball Composition Theory','(now) Affine Composition Theory']
parts=['core law\n∘ₖ','homothety\nk','curvature\ndial d','conics &\nwalls','round→pointy\nℓp shapes','simplex','cubic\nform C','affine\nsphere','rungs','order\n(on a dial)']
# 2 full, 1 partial, 0 none
M=np.array([
 [2,2,2,2,1,1,2,2,1,2],
 [2,2,2,2,0,0,0,1,0,0],
 [2,2,2,2,2,2,2,2,0,0],
 [1,0,0,0,2,1,0,0,0,0],
 [2,1,0,0,0,0,0,0,1,1]])
fig,ax=plt.subplots(figsize=(12,4.2))
col={2:BLUE,1:'#a9c6e8',0:'#eeeeee'}
for i in range(M.shape[0]):
  for j in range(M.shape[1]):
    ax.add_patch(Rectangle((j,-i),0.95,0.9,color=col[M[i,j]]))
    ax.text(j+0.47,-i+0.45,{2:'yes',1:'part',0:''}[M[i,j]],ha='center',va='center',fontsize=9,color='white' if M[i,j]==2 else '#333')
ax.set_xlim(-0.1,len(parts)); ax.set_ylim(-M.shape[0]+0.8,1.9)
for j,p in enumerate(parts): ax.text(j+0.47,1.0,p,ha='center',va='bottom',fontsize=9)
for i,n in enumerate(names): ax.text(-0.15,-i+0.45,n,ha='right',va='center',fontsize=10,weight='bold' if i==0 else 'normal')
ax.axis('off')
ax.set_title('Figure 1. What each name\'s words honestly cover (my judgement, not a measurement)',fontsize=11,loc='left',x=-0.33)
plt.savefig('fig1-coverage.png',dpi=150,bbox_inches='tight'); plt.close()

# Fig 2: composition, same centre vs different centres
fig,axs=plt.subplots(1,2,figsize=(11,3.6))
def H(c,r,y): return c+r*(y-c)
for ax,(c1,c2),title in [(axs[0],(-1,-1),'Same centre (our law, k = 1): either order lands in the same place'),
                          (axs[1],(-1,3.0),'Different centres: the two orders land apart')]:
  x0=0.0; r1,r2=1.6,1.4
  a=H(c1,r1,H(c2,r2,x0)); b=H(c2,r2,H(c1,r1,x0))
  ax.axhline(0,color=GREY,lw=1)
  for c,lab in [(c1,'centre A'),(c2,'centre B')]:
    ax.plot(c,0,'k^',ms=9); ax.text(c,-0.35,lab if c1!=c2 else 'shared centre\n(the wall −1/k)',ha='center',va='top',fontsize=9)
  ax.plot(x0,0,'o',color='k'); ax.text(x0,-0.35 if c1!=c2 else 0.25,'start',ha='center',fontsize=9)
  ax.annotate('',xy=(a,0.9),xytext=(x0,0.9),arrowprops=dict(arrowstyle='->',color=BLUE,lw=2)); ax.text((a+x0)/2,1.05,'B then A',color=BLUE,ha='center',fontsize=9)
  ax.annotate('',xy=(b,-0.9+0.0),xytext=(x0,-0.9),arrowprops=dict(arrowstyle='->',color=ORANGE,lw=2)); ax.text((b+x0)/2,-1.25,'A then B',color=ORANGE,ha='center',fontsize=9)
  ax.plot([a,a],[0.9,0],':',color=BLUE); ax.plot([b,b],[-0.9,0],':',color=ORANGE)
  if abs(a-b)>1e-9: ax.annotate('',xy=(a,0.15),xytext=(b,0.15),arrowprops=dict(arrowstyle='<->',color='k')); ax.text((a+b)/2,0.3,'gap = order (holonomy ω)',ha='center',fontsize=9)
  ax.set_ylim(-1.8,1.6); ax.set_xlim(-2,4.5); ax.axis('off'); ax.set_title(title,fontsize=10)
fig.suptitle('Figure 2. Order is off at the base and a dial turns it on: same centre, the orders agree; move the centres apart and order costs an oriented area',fontsize=10,y=1.04)
plt.savefig('fig2-composition.png',dpi=150,bbox_inches='tight'); plt.close()

# Fig 3: affine spheres round and pointy
fig=plt.figure(figsize=(11,4.6))
ax=fig.add_subplot(121,projection='3d')
u=np.linspace(0,1.3,40); v=np.linspace(0,2*np.pi,60); U,V=np.meshgrid(u,v)
ax.plot_surface(np.sinh(U)*np.cos(V),np.sinh(U)*np.sin(V),np.cosh(U),color=BLUE,alpha=0.7,lw=0)
R=np.linspace(0,1.8,10); ax.plot_wireframe(np.outer(np.cos(v),R).T*1,np.outer(np.sin(v),R).T,np.tile(R,(60,1)).T,color=GREY,lw=0.4,rstride=1,cstride=6)
ax.set_title('Round: one sheet of a hyperboloid\ninside the round cone (cubic form C = 0)',fontsize=10); ax.set_axis_off(); ax.view_init(18,30)
ax=fig.add_subplot(122,projection='3d')
s=np.linspace(-1.2,1.2,40); S,T=np.meshgrid(s,s)
X=np.exp(S); Y=np.exp(T); Z=1/(X*Y)
m=(Z<4)&(X<4)&(Y<4); Z=np.where(m,Z,np.nan)
ax.plot_surface(X,Y,Z,color=ORANGE,alpha=0.8,lw=0)
for e in [(0,0,1),(0,1,0),(1,0,0)]:
  ax.plot([0,4*e[0]],[0,4*e[1]],[0,4*e[2]],color='k',lw=1.2)
ax.set_title("Pointy: Țițeica's surface xyz = 1\ninside the corner of a cube (cone over a triangle)",fontsize=10); ax.set_axis_off(); ax.view_init(22,40)
fig.suptitle('Figure 3. Affine spheres come in round and pointy kinds (Cheng–Yau: one for every convex cone)',fontsize=11,y=1.06)
plt.savefig('fig3-affine-spheres.png',dpi=150,bbox_inches='tight'); plt.close()
print('ok')
