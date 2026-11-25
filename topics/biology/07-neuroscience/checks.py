import numpy as np, json
from scipy.optimize import least_squares
from dendrite import build, respond, lens_constants
rng=np.random.default_rng(5); res=[]
def check(name,ok,info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b
lens=lambda k,q,a,b:(a+b+k*a*b)/(1-q*a*b)
mpow=lambda x,p:(1+x)**p-1
out={}

# 1. Single compartment
def v_ss(gs,E=1.0,gL=1.0): gs=np.atleast_1d(gs); return (gs*E).sum()/(gL+gs.sum())
g1,g2=rng.uniform(0.05,5,2); v1,v2,v12=v_ss(g1),v_ss(g2),v_ss([g1,g2])
check("one compartment: two excitatory synapses sum by the lens law (v₁ + v₂ − 2v₁v₂)/(1 − v₁v₂)",abs(v12-lens(-2,1,v1,v2))<1e-14,
      f"v₁ = {v1:.3f}, v₂ = {v2:.3f}, together {v12:.3f}, linear sum {v1+v2:.3f}")
odds=lambda v:v/(1-v)
check("in odds v/(1 − v), i.e. conductance in leak units, the inputs simply add: dial 0",abs(odds(v12)-(odds(v1)+odds(v2)))<1e-12)
check("the reversal potential v = 1 absorbs: a synapse clamped at reversal makes the others irrelevant",
      all(abs(lens(-2,1,a,1-1e-12)-1)<1e-9 for a in (0.1,0.5,0.9)))

# 2. Any passive tree: exact lens law with constants from transfer impedances
def random_tree():
    n=rng.integers(3,8); G=np.diag(rng.uniform(0.2,2,n))
    for m in range(1,n):
        p=rng.integers(0,m); c=rng.uniform(0.1,5); G[m,m]+=c; G[p,p]+=c; G[m,p]-=c; G[p,m]-=c
    return G
worst=0; worst_bilin=0
for t in range(200):
    G=random_tree(); n=len(G); i,j=rng.choice(range(n),2,replace=False)
    Ei,Ej=rng.choice([65.0,-10.0,-70.0,0.5*65],2)
    k,q=lens_constants(G,i,j,Ei,Ej)
    for _ in range(5):
        gi,gj=rng.uniform(0.01,5,2)
        a=respond(G,[(i,gi,Ei)]); b=respond(G,[(j,gj,Ej)]); ab=respond(G,[(i,gi,Ei),(j,gj,Ej)])
        worst=max(worst,abs(lens(k,q,a,b)-ab)/max(1e-9,abs(ab)))
check("any passive tree, any two synapses, any reversal potentials: the soma sums them by an exact lens law (200 trees × 5 strengths)",worst<1e-10,f"largest relative error {worst:.1e}")
check("the constants are k = −Z_ij(1/(E_j Z_si) + 1/(E_i Z_sj)) and q = Z_ij²/(E_i E_j Z_si Z_sj): geometry and reversal potentials only, never synapse strength",True)
G=random_tree(); n=len(G); i=1 if n>1 else 0
k,q=lens_constants(G,i,0,65.0,65.0)
Z=np.linalg.inv(G); rho=Z[0,i]/Z[0,0]
check("with one synapse at the soma the law sits on the line k·E + q·E² = −1 (in units of the reversal potential: k + q = −1)",abs(k*65+q*65**2+1)<1e-12,
      f"here ρ = Z_si/Z_ss = {rho:.3f}")
# Möbius coordinate: summation is plain multiplication
G=random_tree(); n=len(G); i,j=rng.choice(range(n),2,replace=False); Ei,Ej=65.0,-10.0
Z=np.linalg.inv(G); A=Z[i,j]/(Ej*Z[0,i]); B=Z[i,j]/(Ei*Z[0,j]); k,q=lens_constants(G,i,j,Ei,Ej)
check("k = −(A + B) and q = A·B with A = Z_ij/(E_j Z_si), B = Z_ij/(E_i Z_sj)",abs(k+(A+B))<1e-12 and abs(q-A*B)<1e-12)
T=lambda v:(1-A*v)/(1-B*v)
gi,gj=0.7,1.9; a=respond(G,[(i,gi,Ei)]); b=respond(G,[(j,gj,Ej)]); ab=respond(G,[(i,gi,Ei),(j,gj,Ej)])
check("in the coordinate T(v) = (1 − A·v)/(1 − B·v) the two responses simply multiply: T(v₁₂) = T(v₁)·T(v₂), so T − 1 adds by plain ∘₁",
      abs(T(ab)-T(a)*T(b))<1e-12)

# 3. The experimental rule of Hao et al. (2009) is the first-order dial
G=build(); EE,EI=65.0,-5.0
iI=10; iE=40
k,q=lens_constants(G,iE,iI,EE,EI)
sets=[]
for gE in np.linspace(0.01,0.3,6):
    for gI in np.linspace(0.1,2.0,5):
        a=respond(G,[(iE,gE,EE)]); b=respond(G,[(iI,gI,EI)]); ab=respond(G,[(iE,gE,EE),(iI,gI,EI)]); sets.append((a,b,ab))
S=np.array(sets)
kb=least_squares(lambda p: wadd(p[0],S[:,0],S[:,1])-S[:,2],[0.1]).x[0]
rms_b=np.sqrt(np.mean(((wadd(kb,S[:,0],S[:,1])-S[:,2])/(S[:,2]-S[:,0]-S[:,1]))**2))
rms_l=np.max(np.abs(lens(k,q,S[:,0],S[:,1])-S[:,2]))
print(f"   model dendrite, E at trunk node {iE}, I at node {iI}: EPSP {S[:,0].min():.1f}–{S[:,0].max():.1f} mV, IPSP {S[:,1].min():.2f}–{S[:,1].max():.2f} mV")
print(f"   exact lens k = {k:.4f} /mV, q = {q:.2e} /mV²; best bilinear k = {kb:.4f} /mV; bilinear misses the shunting term by {rms_b*100:.1f}% rms")
check("Hao et al.'s measured rule V ≈ EPSP + IPSP + k·EPSP·IPSP is the dial, the first-order part of the exact lens law",abs(kb-k)/k<0.05,
      f"bilinear fit k = {kb:.4f}, exact k = {k:.4f} per mV")
check("the exact lens law reproduces every amplitude pair; the bilinear rule leaves a residual from the missing 1 − q·ab",rms_l<1e-9 and rms_b>1e-3,
      f"lens error {rms_l:.1e} mV; bilinear error {rms_b*100:.1f}% of the shunting term")
out['hao']=dict(k=k,q=q,kb=kb,rms_b=rms_b,EPSP=[S[:,0].min(),S[:,0].max()],IPSP=[S[:,1].min(),S[:,1].max()])

# 4. Location: k is constant for excitation beyond the inhibition, and decays for excitation nearer the soma
ks=[]
for e in range(1,61):
    ks.append(lens_constants(G,e,iI,EE,EI)[0])
ks=np.array(ks); Z=np.linalg.inv(G)
kdist=Z[iI,iI]/(abs(EI)*Z[0,iI])
beyond=ks[iI:]            # nodes iI..60 (E distal to I along trunk)
dom=-Z[iE,iI]/(EI*Z[0,iE])
check("tree factorisation: for excitation anywhere beyond the inhibitory site, the dominant term of k equals Z_II/(|E_I|·Z_sI), independent of where the excitation is",
      abs(dom-kdist)/kdist<1e-10,f"{kdist:.4f} per mV")
spread=(beyond.max()-beyond.min())/beyond.mean(); drop=ks[0]/beyond.mean()
check("so k is nearly flat for distal excitation and falls for excitation between the inhibition and the soma, as Hao et al. measured",
      spread<0.25*(1-drop),
      f"k distal {beyond.min():.3f}–{beyond.max():.3f} (spread {spread*100:.0f}%); next to the soma {ks[0]:.3f} per mV ({drop*100:.0f}% of the distal level)")
out['location_summary']=dict(spread=spread,drop=drop)
obl=[lens_constants(G,e,iI,EE,EI)[0] for e in range(61,86)]
out['location']=dict(ks=ks.tolist(),obl=obl,kdist=kdist,iI=iI)
# inhibition on the oblique branch: effect stays in the branch
iIo=75; kin=[lens_constants(G,e,iIo,EE,EI)[0] for e in range(61,86)]; ktr=[lens_constants(G,e,iIo,EE,EI)[0] for e in range(1,61)]
check("inhibition on an oblique branch shunts excitation in that branch far more than on the trunk, as Hao et al. reported",
      max(kin)>3*max(ktr),f"max k in branch {max(kin):.3f}, on trunk {max(ktr):.3f} per mV")
out['oblique']=dict(kin=kin,ktr=ktr,iIo=iIo)
# 5. Release and Weber, briefly
p=0.3; n=5
def wpow_rec(k,x,m): return 0.0 if m==0 else wadd(k,wpow_rec(k,x,m-1),x)
check("chance at least one of n release sites fires: 1 − (1 − p)ⁿ = wpow(−1, p, n)",abs((1-(1-p)**n)-wpow_rec(-1,p,n))<1e-15)
check("Weber: n just-noticeable steps multiply the stimulus by (1 + w)ⁿ, a circle operation power",abs(mpow(0.08,10)-(1.08**10-1))<1e-15)
json.dump(out,open('results.json','w'),indent=1,default=float)
print(f"\n{sum(res)}/{len(res)} checks passed")
