"""Test: do the four 'quarter turns' of water oxidation predict the overpotential
better than the standard descriptor (dG_O - dG_OH)?  numpy only.
Real points: published DFT free energies quoted in the text of
 Romeo, Illas, Calle-Vallejo, Chem. Sci. 2023 (BaNiO2) and Calle-Vallejo, Acc. Chem. Res. 2025 (Sr2NiIrO6),
 and the IrO2/RuO2 values quoted in arXiv 2512.05938.
Simulated points: catalysts drawn around the published 3.2 +/- 0.2 eV OOH-OH scaling line."""
import numpy as np, json
E0=1.23; TOT=4.92
def steps(oh,o,ooh): return np.array([oh,o-oh,ooh-o,TOT-ooh])
def arcs(s): return 360*s/TOT
def eta(s): return s.max(axis=-1)-E0
def eta_orr(s): return E0-s.min(axis=-1)
def essi(s):
    out=[]
    for r in np.atleast_2d(s):
        up=r[r>=E0]; out.append(up.mean()-E0)
    return np.array(out)
def rms(s): return np.sqrt(((s-E0)**2).mean(axis=-1))
res={}
# 1 identity: largest arc minus 90 deg, times 1.23/90, equals the overpotential
real={'BaNiO2 (Romeo 2023)':(1.85,3.46,4.78),'Sr2NiIrO6 (Calle-Vallejo 2025)':(0.77,2.12,3.25),
      'Sr2NiIrO6 optimized':(0.93,2.44,3.41),'IrO2/RuO2 typical (arXiv 2512.05938)':(1.23,2.83,4.43)}
pub={'BaNiO2 (Romeo 2023)':dict(eta=0.62,essi=0.36),'Sr2NiIrO6 (Calle-Vallejo 2025)':dict(eta=0.44,essi=0.28),'Sr2NiIrO6 optimized':dict(eta=0.28),'IrO2/RuO2 typical (arXiv 2512.05938)':dict(eta=0.37)}
rows=[]
for k,(oh,o,ooh) in real.items():
    s=steps(oh,o,ooh); a=arcs(s)
    rows.append(dict(name=k,steps=s.round(3).tolist(),arcs=a.round(1).tolist(),
      eta_from_largest_arc=round((a.max()-90)*E0/90,3),eta=round(float(eta(s)),3),essi=round(float(essi(s)[0]),3),
      rms=round(float(rms(s)),3),descriptor=round(o-oh,3),round_trip=round(float(eta(s)+eta_orr(s)),3),published=pub[k]))
res['real_points']=rows
# 2 simulated catalysts around the published scaling line
rng=np.random.default_rng(7); N=20000
oh=rng.uniform(0.0,2.2,N); d=rng.uniform(0.6,2.6,N)          # dG_OH and the descriptor dG_O - dG_OH
o=oh+d; ooh=oh+3.2+rng.normal(0,0.2,N)                        # published 3.2 +/- 0.2 eV link
S=np.stack([oh,o-oh,ooh-o,TOT-ooh],axis=1)
ok=(S>-0.5).all(axis=1); S=S[ok]; d=d[ok]; oh=oh[ok]
true=eta(S)
volcano=np.maximum(d,3.2-d)-E0                                 # standard one-descriptor prediction
es=essi(S); rm=rms(S)
def fit(x,y):
    A=np.vstack([x,np.ones_like(x)]).T; c,*_=np.linalg.lstsq(A,y,rcond=None); p=A@c
    return dict(slope=round(float(c[0]),3),intercept=round(float(c[1]),3),mae=round(float(np.abs(p-y).mean()),3),r2=round(float(1-((p-y)**2).sum()/((y-y.mean())**2).sum()),3))
res['simulated_N']=int(len(true))
res['predict_eta']=dict(largest_arc='exact by definition (MAE 0)',
  descriptor_volcano=dict(mae=round(float(np.abs(volcano-true).mean()),3),r2=round(float(1-((volcano-true)**2).sum()/((true-true.mean())**2).sum()),3)),
  essi_fit=fit(es,true),rms_fit=fit(rm,true))
# among good catalysts only (eta<0.6)
g=true<0.6
res['predict_eta_good_only']=dict(n=int(g.sum()),descriptor_volcano_mae=round(float(np.abs(volcano[g]-true[g]).mean()),3),essi_fit=fit(es[g],true[g]),rms_fit=fit(rm[g],true[g]))
# 3 round trip (both directions): descriptor alone vs descriptor + dG_OH
rt=eta(S)+eta_orr(S)
res['round_trip']=dict(descriptor_only=fit(volcano,rt),
  spread_of_arcs='exact by definition',
  best_possible_on_scaling_line=None)
# analytic floor with exact 3.2 link: minimise max-min over oh,d
ohg,dg=np.meshgrid(np.linspace(0,2.2,441),np.linspace(0.6,2.6,401))
Sg=np.stack([ohg,dg,3.2-dg,TOT-ohg-3.2],axis=-1)
spread=Sg.max(-1)-Sg.min(-1); i=np.unravel_index(spread.argmin(),spread.shape)
res['round_trip']['best_possible_on_scaling_line']=dict(spread_V=round(float(spread[i]),3),dG_OH=round(float(ohg[i]),3),descriptor=round(float(dg[i]),3))
# 4 how much DFT noise would hide: add 0.2 and 0.5 eV noise to each adsorbate, see how often the ranking of two catalysts differing by 0.1 V flips
flip={}
for sig in (0.1,0.2,0.5):
    n=20000; base=rng.normal(0,sig*np.sqrt(2),n)  # difference of two max-step errors, rough
    flip[sig]=round(float((base>0.1).mean()),3)
res['ranking_flip_probability_for_0p1V_gap']=flip
json.dump(res,open('quarter_turns_output.json','w'),indent=1)
print(json.dumps(res,indent=1))
