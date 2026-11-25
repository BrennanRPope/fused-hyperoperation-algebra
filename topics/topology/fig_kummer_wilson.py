"""Pictures for Kummer's theorem (carries) and Wilson's theorem on the mult-add circle."""
import numpy as np, matplotlib
from math import comb, factorial
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
B='#1d2bd8'; R='#cf3519'; K='#15131a'; M='#b9b7c2'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'savefig.dpi':160,'savefig.bbox':'tight','savefig.facecolor':'white'})

def v(p,x):
    e=0
    while x%p==0: x//=p; e+=1
    return e
def carries(p,a,b):
    c=0; n=0
    while a or b or c:
        s=a%p+b%p+c; c=s//p; n+=c; a//=p; b//=p
    return n

# checks behind the pictures
for p in (2,3,5,7):
    for n in range(60):
        for k in range(n+1):
            assert v(p,comb(n,k))==carries(p,k,n-k)
for p in (2,3,5,7,11,13):
    assert factorial(p-1)%p==p-1
for n in range(2,60):
    isprime=all(n%d for d in range(2,n))
    assert (factorial(n-1)%n==n-1)==isprime
for n in range(1,9):
    assert sum((-1)**k*comb(n,k)*(n-k)**n for k in range(n+1))==factorial(n)
print('kummer/wilson checks pass')

# Figure 11: Kummer
fig,axs=plt.subplots(1,2,figsize=(13,5.6),gridspec_kw={'width_ratios':[1.35,1]})
ax=axs[0]; p=3; N=27
cm=['#f1f0f4','#9fa8f2',B,R]
for n in range(N):
    for k in range(n+1):
        e=v(p,comb(n,k)); assert e<=2
        ax.add_patch(plt.Rectangle((k-n/2-0.5,-n-0.5),1,1,facecolor=cm[e],edgecolor='white',lw=0.4))
ax.set_xlim(-N/2-1,N/2+2); ax.set_ylim(-N-6.5,1); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Layers of the n-cube (rows n = 0 to 26), coloured by\nhow many times 3 divides the count',fontsize=12)
for e,lab in enumerate(['no factor of 3 (no carries)','exactly 3 (1 carry)','exactly 3² (2 carries)']):
    x0=-N/2; y0=-N-2.5-e*1.7
    ax.add_patch(plt.Rectangle((x0,y0),1.2,1.2,facecolor=cm[e],edgecolor=K,lw=0.5))
    ax.text(x0+1.6,y0+0.6,lab,fontsize=10,va='center')
ax=axs[1]; ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(0,10)
ax.set_title('Why: 5 + 4 = 9 in base 3',fontsize=12)
rows=[('carry','1 1 ',R),('5 =','  1 2',K),('4 =','+ 1 1',K),('9 =','1 0 0',B)]
for i,(a,b,c) in enumerate(rows):
    ax.text(2.2,8.6-i*1.25,a,fontsize=15,ha='right',color=c,family='DejaVu Sans Mono')
    ax.text(2.6,8.6-i*1.25,b,fontsize=17,color=c,family='DejaVu Sans Mono')
ax.plot([2.5,5.6],[5.55,5.55],color=K,lw=1.2)
ax.text(0.2,3.6,'Two carries, so 3² is the exact power\nof 3 in C(9,4) = 126 = 9 × 14.',fontsize=12)
ax.text(0.2,1.3,'On the cube: corners of the 9-cube with\n4 switches on. Each carry is one more\nfactor of p in the count.',fontsize=12,color='#444')
plt.savefig('figures/t11_kummer_carries.png'); plt.close()

# Figure 12: Wilson on the mult-add circle, p = 13
p=13; g=2  # 2 generates the units mod 13
pos={}
x=1
for j in range(p-1):
    pos[(x-1)%p]=j  # point a = x - 1 sits at turn j/(p-1)
    x=x*g%p
fig,axs=plt.subplots(1,2,figsize=(13,6))
ax=axs[0]
th=lambda j: np.pi/2-2*np.pi*j/(p-1)
ax.add_patch(plt.Circle((0,0),1,fill=False,color=M,lw=1.2))
done=set()
for a,j in pos.items():
    inv=[b for b in pos if ((1+a)*(1+b)-1)%p==0][0]
    if a in done: continue
    done.update({a,inv})
    if inv!=a:
        t1,t2=th(j),th(pos[inv])
        ax.plot([np.cos(t1),np.cos(t2)],[np.sin(t1),np.sin(t2)],color=B,lw=1.4,alpha=0.75)
for a,j in pos.items():
    t=th(j); fixed=((1+a)*(1+a)-1)%p==0
    ax.add_patch(plt.Circle((np.cos(t),np.sin(t)),0.07,facecolor=R if fixed else 'white',edgecolor=K,lw=1,zorder=3))
    ax.text(1.2*np.cos(t),1.2*np.sin(t),str(a),ha='center',va='center',fontsize=11,color=R if fixed else K)
ax.text(0,-1.55,'Mult-add mod 13 turns the circle: each point is a turn.\nBlue chords join each point to its mult-add inverse.\nOnly 0 and 11 (= −2) are their own inverses.\nAdding every point by mult-add lands on 11, the half turn:\n(1+0)(1+1)···(1+11) = 12! ≡ −1, which is Wilson.',ha='center',va='top',fontsize=10.5)
ax.set_xlim(-1.45,1.45); ax.set_ylim(-2.75,1.4); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('The mult-add circle for p = 13',fontsize=12)
ax=axs[1]; p=7; n=p-1
terms=[(-1)**k*comb(n,k)*(n-k)**n for k in range(n+1)]
mods=[t%p for t in terms]
ax.bar(range(n+1),mods,color=[B if m==1 else M for m in mods],edgecolor=K,width=0.6)
for k,t in enumerate(terms):
    ax.text(k,max(mods[k],0)+0.06,f'{t:+,}',ha='center',fontsize=9,rotation=0)
ax.set_xticks(range(n+1)); ax.set_xlabel('layer k of the 6-cube (k values left out)')
ax.set_ylabel('contribution mod 7'); ax.set_ylim(0,1.6); ax.set_yticks([0,1])
ax.set_title('What Lean checks: inclusion–exclusion on the 6-cube\nsum = 720 = 6!, and mod 7 it is 1+1+1+1+1+1+0 = 6 ≡ −1',fontsize=11.5)
for s in ('top','right'): ax.spines[s].set_visible(False)
plt.savefig('figures/t12_wilson_circle.png'); plt.close()
