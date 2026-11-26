"""DRAFT figures for three-body-order. Orange = known, blue = ours, red = does not work, grey = neither."""
import json,sys,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'../three-body'); import tbpy
from words import rho, readings, walk_area
from physics import middle
OR,BL,RD,GR='#e8871e','#2b6cb0','#c53030','#888888'
plt.rcParams.update({'font.size':10,'figure.dpi':130})
R=json.load(open('words.json')); T=np.array([r['Tstar'] for r in R])
K=np.zeros((len(R),5))
for i,r in enumerate(R):
    for k in r['runs']: K[i,min(k,5)-1]+=1
dbl=np.array([sum(1 for j in range(len(r['word'])) if r['word'][j]==r['word'][(j+1)%len(r['word'])]) for r in R],float)
def loo(X,y):
    H=X@np.linalg.pinv(X); e=(y-H@y)/(1-np.diag(H))/y; return np.median(abs(e))*100
def save(f,name): f.tight_layout(); f.savefig('figures/'+name,bbox_inches='tight'); plt.close(f)

# 1 ladder
f,ax=plt.subplots(figsize=(10,6.2)); ax.axis('off')
rungs=[('1. The wedge (oriented area)', [('On the shape sphere: turn between line-ups = half the enclosed area (swap rule)',OR,'holds exactly, known'),
         ('On the word: signed area of the letter walk',RD,'time reversal flips it; adds nothing')]),
       ('2. Whole-number windings', [('How many times each pair circles per period',GR,'fixed by the word; twins identical')]),
       ('3. Letters as maps that compose', [('Stretch-and-shift (affine): letters as digits',OR,'already the word-to-start recipe'),
         ('Order dial t (Burau); t = -1 is the Mobius / modular picture',BL,'runs = cusp trips; size of the product adds little'),
         ('Order of the runs inside a word',BL,'real but small (0.12% to 0.03%)')]),
       ('4. More derived', [('Temperley-Lieb, 3-sphere, non-commuting rungs, path integrals',GR,'same as rung 3, done before, or not testable here')])]
y=0.97
for title,items in rungs:
    ax.text(0.0,y,title,fontsize=12,weight='bold',transform=ax.transAxes); y-=0.065
    for txt,col,ver in items:
        ax.add_patch(plt.Rectangle((0.02,y-0.035),0.6,0.055,color=col,alpha=0.18,transform=ax.transAxes))
        ax.text(0.03,y-0.008,txt,fontsize=9.5,transform=ax.transAxes,va='center')
        ax.text(0.64,y-0.008,ver,fontsize=9.5,color=col,transform=ax.transAxes,va='center',weight='bold'); y-=0.075
    y-=0.02
ax.annotate('',xy=(0.97,0.08),xytext=(0.97,0.95),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=2,color=GR))
ax.text(0.985,0.5,'more derived',rotation=270,transform=ax.transAxes,va='center',color=GR)
ax.set_title('Our order tools, most fundamental first, applied to three bodies',loc='left',fontsize=12)
save(f,'fig1-ladder.png')

# 2 swap rule
d={o['cls']:o for o in json.load(open('../three-body-word-formula/data/li-liao-695.json'))}
o=d['I.A2']; ts,X=tbpy.traj(tbpy.suvakov(o['v1'],o['v2']),o['T'],8000); S=tbpy.shape(X)
idx=np.where(np.sign(S[:-1,2])*np.sign(S[1:,2])<0)[0]; i=idx[0]; mi=middle(X[i])
j=[jj for jj in idx[1:] if middle(X[jj])==mi][0]
f=plt.figure(figsize=(11,4.8)); ax=f.add_subplot(121,projection='3d')
u=np.linspace(0,2*np.pi,60); v=np.linspace(0,np.pi,30)
ax.plot_wireframe(np.outer(np.cos(u),np.sin(v)),np.outer(np.sin(u),np.sin(v)),np.outer(np.ones(60),np.cos(v)),color='#ddd',lw=0.4)
ax.plot(*S.T,color=GR,lw=0.6,alpha=0.6,label='whole orbit I.A2')
ax.plot(*S[i:j+1].T,color=BL,lw=2.2,label='arc between two line-ups, same middle body')
g=np.linspace(0,2*np.pi,200); ax.plot(np.cos(g),np.sin(g),0*g,color='k',lw=0.8,label='equator = line-ups (area of triangle zero)')
for c,lab in [(np.array([-1,0,0]),'collision')]: pass
coll=tbpy.shape(np.array([[0,0,0,0,1,0]+[0]*6,[0,0,1,0,0,0]+[0]*6,[1,0,0,0,0,0]+[0]*6],float)+0.0)
ax.scatter(*np.array([[-1,0,0],[0.5,np.sqrt(3)/2,0],[0.5,-np.sqrt(3)/2,0]]).T,color=RD,s=30,label='collisions (removed)')
ax.set_box_aspect((1,1,1)); ax.set_axis_off(); ax.legend(loc='lower left',fontsize=7.5); ax.set_title('The shape sphere: each point is a triangle shape',fontsize=10)
ax2=f.add_subplot(122); P=json.load(open('swap_polished.json'))
for c,e in P.items():
    rows=np.array(e['rows']); ax2.scatter((-rows[:,0]/2)%np.pi,rows[:,1],s=14,label=c)
ax2.plot([0,np.pi],[0,np.pi],color=OR,lw=1.5,label='turn = half the enclosed area')
ax2.set_xlabel('minus half the enclosed area on the shape sphere (mod pi)'); ax2.set_ylabel('turn of the line-up (mod pi)')
ax2.set_title('Swap rule: 123 arcs from 6 orbits, worst miss 1e-4 rad',fontsize=10); ax2.legend(fontsize=7)
save(f,'fig2-swap-rule-on-the-shape-sphere.png')

# 3 word walk
f,(a1,a2)=plt.subplots(1,2,figsize=(11,4.6))
step={'a':(1,0),'A':(-1,0),'b':(0,1),'B':(0,-1)}
for w,col,lab in [('BabA',BL,'figure-eight BabA'),(d['I.A2']['word'],OR,'I.A2 '+d['I.A2']['word'])]:
    p=[(0,0)]
    for ch in w: p.append((p[-1][0]+step[ch][0],p[-1][1]+step[ch][1]))
    p=np.array(p,float)+(0.06 if col==OR else 0); a1.plot(p[:,0],p[:,1],'-o',color=col,label=lab+f'  (area {walk_area(w)[2]:+g})'); a1.fill(p[:,0],p[:,1],color=col,alpha=0.15)
a1.set_xlabel('net a steps'); a1.set_ylabel('net b steps'); a1.set_aspect('equal'); a1.grid(alpha=0.3); a1.legend(fontsize=8)
a1.set_title('A word as a walk: a = right, b = up, capitals go back',fontsize=10)
closed=np.array([r['netA']==0 and r['netB']==0 for r in R]); area=np.array([r['area'] for r in R])
c=np.linalg.lstsq(K,T,rcond=None)[0]; res=(T-K@c)/T*100
a2.scatter(np.abs(area[closed]),res[closed],s=10,color=RD,alpha=0.6)
a2.set_xlabel('size of the enclosed area |area| (closed walks, 311 words)'); a2.set_ylabel('period left over after run prices (%)')
a2.set_title('The word\'s wedge says nothing about the period (r = %.2f)'%np.corrcoef(np.abs(area[closed]),res[closed])[0,1],fontsize=10)
save(f,'fig3-the-word-as-a-walk.png')

# 4 windings
W=json.load(open('wind_cat.json')); TW=json.load(open('twin_windings.json'))
f,(a1,a2)=plt.subplots(1,2,figsize=(11,4.4))
n=np.array([e['netA'] for e in W]); w=np.array([[e['w'][k] for k in ['12','13','23']] for e in W])
jit=np.random.default_rng(0).normal(0,0.05,len(n))
a1.scatter(n+jit,w[:,1],color=BL,s=14,label='pair 1-3'); a1.scatter(n+jit,w[:,2],color=OR,s=14,label='pair 2-3'); a1.scatter(n+jit,w[:,0],color=GR,s=14,label='pair 1-2')
a1.set_xlabel('net count of a in the word'); a1.set_ylabel('times the pair circles in one period'); a1.legend(fontsize=8); a1.grid(alpha=0.3)
a1.set_title('90 shortest orbits: windings are (0, n, -n), set by the word',fontsize=10)
labs=list(TW.keys()); x=np.arange(len(labs))
a2.bar(x-0.18,[TW[c]['cat']['13'] for c in labs],0.36,color=OR,label='catalog orbit'); a2.bar(x+0.18,[TW[c]['twin']['13'] for c in labs],0.36,color=BL,label='its same-word twin')
a2.set_xticks(x); a2.set_xticklabels(labs,rotation=30,fontsize=8); a2.set_ylabel('pair 1-3 windings per period'); a2.legend(fontsize=8)
a2.set_title('Twins circle exactly alike: windings cannot tell them apart',fontsize=10)
save(f,'fig4-whole-number-windings.png')

# 5 runs are cusp trips
ell=np.array([r['ell'] for r in R]); L=np.array([r['L'] for r in R])
ce=np.linalg.lstsq(K,ell,rcond=None)[0]; ct=np.linalg.lstsq(K,T,rcond=None)[0]; ks=np.arange(1,6)
f,(a1,a2)=plt.subplots(1,2,figsize=(11,4.4))
a1.plot(ks,(ct/ks)/(ct[1]/2),'-o',color=OR,label='real time per letter (period law)')
a1.plot(ks,(ce/ks)/(ce[1]/2),'-s',color=BL,label='Mobius length per letter (t = -1)')
a1.axhline(1,color=GR,lw=0.6); a1.set_xlabel('letters in the run = close passes'); a1.set_ylabel('price per letter, runs of 2 = 1')
a1.set_title('Same pieces, opposite trend: time per pass grows, Mobius length per pass shrinks',fontsize=9.5); a1.legend(fontsize=8)
frac=(K[:,2]*3+K[:,3]*4+K[:,4]*5)/L
sc=a2.scatter(ell/L,T/L,c=frac,cmap='viridis',s=9); plt.colorbar(sc,ax=a2,label='share of letters in runs of 3 or more')
a2.set_xlabel('Mobius length per letter'); a2.set_ylabel('period per letter'); a2.set_title('Long runs: short on the modular surface, long in time (r = %.2f)'%np.corrcoef(ell/L,T/L)[0,1],fontsize=9.5)
save(f,'fig5-runs-are-cusp-trips.png')

# 6 burau dial
ths=np.linspace(0,1,41); m1=[];m2=[]
for th in ths:
    lr=np.array([np.log(max(rho(r['word'],np.exp(1j*np.pi*th)),1.0)) for r in R])
    m1.append(loo(np.column_stack([K,lr]),T)); m2.append(loo(np.column_stack([K,dbl,lr]),T))
base=loo(K,T); based=loo(np.column_stack([K,dbl]),T)
f,ax=plt.subplots(figsize=(8,4.4))
ax.plot(ths,m1,color=BL,lw=2,label='run prices + size of the composed letter matrices'); ax.axhline(base,color=GR,ls='--',label='run prices alone')
ax.plot(ths,m2,color=RD,lw=2,label='... once doubled letters (aa, bb) are counted'); ax.axhline(based,color=RD,ls=':',label='run prices + doubled letters')
ax.set_xlabel('order dial t = exp(i pi x): x = 0 forgets every pure braid, x = 1 is the Mobius picture'); ax.set_ylabel('typical period miss, held out (%)')
ax.set_title('Letters as matrices: the dial helps only by noticing doubled letters',fontsize=10); ax.legend(fontsize=8)
save(f,'fig6-order-dial-on-letters.png')
json.dump(dict(theta=ths.tolist(),runs_plus_rho=m1,runs_dbl_rho=m2,base=base,base_dbl=based),open('dial_sweep.json','w'))

# 7 run order control
C=json.load(open('pairs_ctrl.json'))
f,ax=plt.subplots(figsize=(7.5,4.2))
ax.bar([0,1],[C['counts'][0],C['real'][0]],color=[GR,BL]); ax.scatter(np.full(len(C['shuffled']),2)+np.random.default_rng(1).normal(0,0.05,len(C['shuffled'])),C['shuffled'],color=RD,s=14)
ax.set_xticks([0,1,2]); ax.set_xticklabels(['how many runs\nof each length','+ which runs sit\nnext to which','same, runs shuffled\n(20 shuffles)'])
ax.set_ylabel('typical period miss, held out (%)'); ax.set_title('Order of runs matters, but only a little: real order beats every shuffle',fontsize=10)
save(f,'fig7-order-of-runs.png')

# 8 stability
S=json.load(open('stability.json')); el=np.array([readings(s['word'])['ell']/2 for s in S]); ly=np.array([s['logmaxlam'] for s in S]); Ls=np.array([len(s['word']) for s in S])
f,ax=plt.subplots(figsize=(7.5,4.6)); sc=ax.scatter(el,ly,c=Ls,cmap='plasma',s=12); plt.colorbar(sc,label='letters in the word')
mx=max(el.max(),ly.max()); ax.plot([0,mx],[0,mx],color=OR,label='braid entropy = instability')
ax.set_xlabel('braid entropy of the word (log of the Mobius stretch, per period)'); ax.set_ylabel('log of the largest Floquet multiplier')
ax.set_title('239 shortest orbits: the word\'s braid entropy does not predict instability',fontsize=10); ax.legend(fontsize=8)
save(f,'fig8-stability-vs-braid-entropy.png')
print('figures done')
