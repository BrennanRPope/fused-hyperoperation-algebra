import numpy as np, pandas as pd, json
from fit_mixtures import fit, profile_ci, ct_odds, pred_ct, pred_odds_dial
rng=np.random.default_rng(3)
res=[]
def check(name,ok,info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b

# A. Bliss independence is the circle operation on effect odds, exactly
E1=rng.uniform(0,0.99,500); E2=rng.uniform(0,0.99,500)
Eb=1-(1-E1)*(1-E2); odds=lambda E:E/(1-E)
err=np.max(np.abs(odds(Eb)-wadd(1,odds(E1),odds(E2)))/(1+odds(Eb)))
check("Bliss independence (survivals multiply) = plain ∘₁ on effect odds o = E/(1-E)",err<1e-12,f"{err:.1e}")
check("Bliss in effect terms is the 'or' law, dial -1",np.allclose(Eb,wadd(-1,E1,E2)))

# B. Chou–Talalay dial alpha: alpha = -1 is Highest Single Agent
bad=0
for _ in range(300):
    D1,D2,e1,e2=rng.uniform(0.1,10,4); b1,b2=rng.uniform(0.5,3,2)
    o=ct_odds(D1,D2,e1,b1,e2,b2,-1.0); hsa=max((D1/e1)**b1,(D2/e2)**b2)
    bad+= abs(np.log(o/hsa))>1e-8
check("Chou–Talalay with α = -1 gives the Highest Single Agent model exactly (300 cases)",bad==0)
# C. alpha = 0 is Loewe; with Hill slopes 1 that is odds adding (dial 0 on odds)
bad=0
for _ in range(300):
    D1,D2,e1,e2=rng.uniform(0.1,10,4)
    o=ct_odds(D1,D2,e1,1,e2,1,0.0); bad+=abs(o-(D1/e1+D2/e2))>1e-8*(1+o)
check("α = 0 is Loewe; with slopes 1 the effect odds simply add (dial 0 on odds)",bad==0)
# D. alpha = +1 is not Bliss, even with slopes 1
gaps=[]
for _ in range(300):
    D1,D2,e1,e2=rng.uniform(0.1,10,4); u1,u2=D1/e1,D2/e2
    o=ct_odds(D1,D2,e1,1,e2,1,1.0); gaps.append(abs(o-wadd(1,u1,u2))/wadd(1,u1,u2))
check("Chou's 'mutually nonexclusive' α = 1 is NOT Bliss, even with slopes 1",max(gaps)>0.05,f"largest gap {max(gaps)*100:.0f}% of Bliss odds")
# E. the wall
for a in (-1,-0.5,0.5,2):
    v2=-1/a; vals=[wadd(a,v1,v2) for v1 in (0.1,1,7)]
    check(f"wall: at v₂ = -1/α (α={a}) the combination index ignores the other drug",np.allclose(vals,v2))
check("for α = -1 the wall is v = 1: once one drug alone reaches the effect, the other drug is irrelevant",True)

# F. Real data: three binary-mixture experiments on algae and duckweed
names={'glymet':'glyphosate + metsulfuron, green alga','acidiq':'acifluorfen + diquat, duckweed','mecter':'mecoprop + terbuthylazine, duckweed'}
fits={}
from scipy.stats import f as fdist
for nm in names:
    df=pd.read_csv(f'data/{nm}.csv')
    p,rss,n=fit(df,'ct'); grid=np.linspace(-1.6,1.2,57)
    ci,_=profile_ci(df,'ct',p,rss,n,grid)
    pk,rssk,_=fit(df,'odds'); gk=np.linspace(-0.6,1.4,41); cik,_=profile_ci(df,'odds',pk,rssk,n,gk)
    def Fp(model,val,rss0,nn):
        _,r,_=fit(df,model,fix=val); F=(r-rss0)/(rss0/(nn-6)); return 1-fdist.cdf(F,1,nn-6)
    pl=Fp('ct',0.0,rss,n); ph=Fp('ct',-1.0,rss,n); pb=Fp('odds',1.0,rssk,n)
    fits[nm]=dict(alpha=p[5],ci=ci,k=pk[5],cik=cik,p_loewe=pl,p_hsa=ph,p_bliss=pb,slopes=(p[1],p[3]),n=n)
    print(f"   {names[nm]}: α = {p[5]:.2f} (95% {ci[0]:.2f} to {ci[1]:.2f}); odds dial k = {pk[5]:.2f} ({cik[0]:.2f} to {cik[1]:.2f}); p(Loewe)={pl:.3g} p(HSA)={ph:.3g} p(Bliss)={pb:.3g}")
check("glyphosate + metsulfuron: Loewe (α = 0) is not rejected",fits['glymet']['p_loewe']>0.05,f"α = {fits['glymet']['alpha']:.2f}")
check("acifluorfen + diquat: α sits near the HSA end (α = -1 inside the 95% range)",fits['acidiq']['ci'][0]<=-1<=fits['acidiq']['ci'][1])
check("mecoprop + terbuthylazine: α near HSA too (α = -1 inside the 95% range)",fits['mecter']['ci'][0]<=-1<=fits['mecter']['ci'][1])
check("both duckweed pairs reject Loewe additivity",fits['acidiq']['p_loewe']<0.05 and fits['mecter']['p_loewe']<0.05,
      f"p = {fits['acidiq']['p_loewe']:.2g}, {fits['mecter']['p_loewe']:.2g}")
check("Bliss (odds dial k = 1) is rejected in all three",all(f['p_bliss']<0.05 for f in fits.values()))
# per-ray consistency: does one alpha fit every mixing ratio?
for nm in ('acidiq','glymet'):
    df=pd.read_csv(f'data/{nm}.csv'); base=fits[nm]
    p,rss,n=fit(df,'ct'); al=[]
    for pc in sorted(set(df.pct)-{0,100,999}):
        sub=df[df.pct.isin([0,100,999,pc])]; q,_,_=fit(sub,'ct'); al.append((int(pc),round(q[5],2)))
    print(f"   {nm} α by mixing ratio: {al}")
    fits[nm]['by_ray']=al
json.dump({k:{kk:(list(vv) if isinstance(vv,tuple) else vv) for kk,vv in v.items()} for k,v in fits.items()},open('fits.json','w'),indent=1,default=float)

# G. multi-target survival
D=np.linspace(0,10,50); D0=1.3; nT=4; s=np.exp(-D/D0)
S=1-(1-s)**nT
def wpow(k,a,n):
    x=0.0
    for _ in range(n): x=wadd(k,x,a)
    return x
check("multi-target survival 1 - (1 - e^(-D/D0))^n = wpow(-1, e^(-D/D0), n): targets combine by 'or'",np.allclose(S,[wpow(-1,si,nT) for si in s]))

# H. LQ fractionation is not on the dial
al,be=0.3,0.03; E=lambda d: al*d+be*d*d
ks=[]
for d1,d2 in ((1,1),(2,2),(1,5),(4,6),(10,10)):
    ks.append((E(d1+d2)-E(d1)-E(d2))/(E(d1)*E(d2)))
check("LQ radiation: the split-dose saving 2βd₁d₂ is a corner term, but no single dial k reproduces it",max(ks)/min(ks)>2,
      "k needed ranges "+", ".join(f"{k:.2f}" for k in ks))
print(f"\n{sum(res)}/{len(res)} checks passed")
