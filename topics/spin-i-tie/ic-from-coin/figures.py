# DRAFT figure for round 9.
import json, math, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('results.json'))
def h(p): return 0.0 if p<=0 or p>=1 else -(p*math.log2(p)+(1-p)*math.log2(1-p))
def H(w): w=w[w>1e-15]; return float(-(w*np.log2(w)).sum())
def gap(x,y):
    best=9
    for a in np.linspace(0,1,401):
        w=np.array([a,(1+x)/2-a,(1+y)/2-a,0]); w[3]=1-w[:3].sum()
        if (w>=-1e-12).all(): best=min(best,H(np.clip(w,0,1)))
    return best-min(h((1+x)/2),h((1+y)/2))
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
a=ax[0]; E=np.array(d['E']); a.plot(E,d['ic_sum'],lw=2,color='C0')
a.axhline(1,color='k',ls='--',lw=1); a.text(0.02,1.05,'information causality: at most the 1 bit sent',fontsize=9)
for e,lab,c,off in [(0.5,'best classical','C2',(-0.27,0.25)),(1/math.sqrt(2),'quantum 2√2','C1',(0.05,-0.45)),(1.0,'PR box','C3',(-0.3,0.0))]:
    y=np.interp(e,E,d['ic_sum']); a.plot(e,y,'o',color=c,ms=8); a.annotate(f'{lab}\n{y:.2f} bits',xy=(e,y),xytext=(e+off[0],y+off[1]),fontsize=9,arrowprops=dict(arrowstyle='->'))
a.set_xlabel('box strength E   (CHSH = 2 + 2E)'); a.set_ylabel('information Bob gets, summed over both bits')
a.set_title('No-signalling plus Shannon records still lets 2 bits through 1',fontsize=10); a.set_ylim(0,2.4)
b=ax[1]; n=41; xs=np.linspace(-0.98,0.98,n); G=np.array([[gap(x,y) for x in xs] for y in xs])
im=b.imshow(G,extent=[-1,1,-1,1],origin='lower',cmap='magma'); fig.colorbar(im,ax=b,label='mixing entropy − measurement entropy (bits)')
b.set_title('A gbit (square): two entropies disagree, so IC is not guaranteed',fontsize=10)
b.set_xlabel('x   (a qubit, whose states fill a ball, has zero gap everywhere: check N5)'); b.set_ylabel('y')
fig.suptitle('DRAFT: information causality needs the states themselves to carry one entropy',fontsize=11)
fig.tight_layout(); fig.savefig('ic2_square.png',dpi=130)
