# DRAFT figure for round 10. Reuses the exact entropy functions from checks.py (top part only).
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('checks.py').read().split("g=np.linspace")[0])
d=json.load(open('results.json')); per={int(k):v for k,v in d['per_n'].items()}
fig=plt.figure(figsize=(13,7.2)); show=[3,4,5,8,13,30]
for i,n in enumerate(show):
    Rn=1.02/np.cos(np.pi/n); g=np.linspace(-Rn,Rn,91); GX=np.array([(x,y) for y in g for x in g])
    a=fig.add_subplot(2,6,i+1); Z=np.full(len(GX),np.nan); m=inside(n,GX)
    Z[m]=np.abs(mix_ent(n,GX[m])-meas_ent(n,GX[m]))
    im=a.imshow(Z.reshape(len(g),len(g)),extent=[g[0],g[-1],g[0],g[-1]],origin='lower',cmap='magma',vmin=0,vmax=1)
    V=verts(n); a.plot(*np.vstack([V,V[:1]]).T,color='gray',lw=0.8); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'{n} sides'+(' (classical)' if n==3 else ''),fontsize=9)
cax=fig.add_axes([0.92,0.56,0.012,0.32]); fig.colorbar(im,cax=cax,label='|mixing − measurement| (bits)')
b=fig.add_subplot(2,2,3); ns=sorted(per)
for par,mk,lab in [(0,'o','even sides'),(1,'s','odd sides')]:
    xs=[n for n in ns if n%2==par and n>=4]; b.plot(xs,[per[n]['max'] for n in xs],mk+'-',ms=4,label=lab)
b.set_ylim(0,1.1); b.set_xlabel('number of sides'); b.set_ylabel('worst disagreement (bits)')
b.set_title('The worst disagreement never shrinks: it sits on the edges',fontsize=10); b.legend(fontsize=8)
c=fig.add_subplot(2,2,4)
for par,mk,lab in [(0,'o','even sides'),(1,'s','odd sides')]:
    xs=[n for n in ns if n%2==par and n>=4]; c.plot(xs,[per[n]['mean'] for n in xs],mk+'-',ms=4,label=lab)
far={int(k):v for k,v in d['coarse_far'].items()}
c.plot([61,62],[far[61],far[62]],'k^',ms=6,label='61 and 62 sides (coarse grid)')
c.set_xlabel('number of sides'); c.set_ylabel('average disagreement (bits)')
c.set_title('On average it fades, so the inside rounds off',fontsize=10); c.legend(fontsize=8)
fig.suptitle('DRAFT: one entropy for mixing and measuring holds only for the triangle and the exact disk',fontsize=11)
fig.subplots_adjust(left=0.06,right=0.9,top=0.9,bottom=0.08,hspace=0.35,wspace=0.3); fig.savefig('pe1_polygons.png',dpi=120)
