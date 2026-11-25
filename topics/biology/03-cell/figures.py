import sys, json, numpy as np, pandas as pd
sys.path.insert(0,"../tools"); from style import *
setup()
from fit_mixtures import fit, pred_ct
fits=json.load(open('fits.json'))
fig,ax=plt.subplots(1,2,figsize=(11,4.4),gridspec_kw={"width_ratios":[1,1.15]})
v=np.linspace(0,1,400)
for a,c,lab in ((-1,S4,"α = −1: Highest Single Agent"),(-0.5,MUTE,"α = −0.5"),(0,S1,"α = 0: Loewe additivity"),(1,S2,"α = 1: Chou 'nonexclusive'")):
    if a==-1:
        ax[0].plot([0,1,1],[1,1,0],color=c,label=lab)
    else:
        v2=(1-v)/(1+a*v); ax[0].plot(v,v2,color=c,label=lab,lw=2 if a!=-0.5 else 1.2,ls="-" if a!=-0.5 else "--")
ax[0].set_xlim(0,1.05); ax[0].set_ylim(0,1.05); ax[0].set_aspect("equal")
ax[0].set_xlabel("v₁ = dose of drug 1 / dose giving the effect alone"); ax[0].set_ylabel("v₂ (same for drug 2)")
ax[0].set_title("Iso-effect lines: v₁ + v₂ + α·v₁v₂ = 1"); ax[0].legend(fontsize=8.5,loc="lower left")
lab={'glymet':'glyphosate + metsulfuron\n(green alga)','acidiq':'acifluorfen + diquat\n(duckweed)','mecter':'mecoprop + terbuthylazine\n(duckweed)'}
for i,(nm,col) in enumerate(zip(('glymet','mecter','acidiq'),(S1,S3,S4))):
    f=fits[nm]; lo,hi=f['ci']
    ax[1].errorbar(f['alpha'],i,xerr=[[f['alpha']-lo],[hi-f['alpha']]],fmt="o",ms=9,color=col,mec="white",mew=2,capsize=4,lw=2)
    ax[1].text(0.62,i,lab[nm],va="center",fontsize=9,color=INK)
for x,t in ((-1,"HSA"),(0,"Loewe"),(1,"nonexclusive")):
    ax[1].axvline(x,color=GRID if x else MUTE,lw=1.2,zorder=0); ax[1].text(x,2.55,t,ha="center",fontsize=9,color=MUTE)
ax[1].set_xlim(-1.4,1.9); ax[1].set_ylim(-0.6,2.8); ax[1].set_yticks([])
ax[1].set_xlabel("fitted interaction dial α  (95% range)"); ax[1].set_title("Three real mixture experiments on the dial")
save(fig,"figures/f1_combination_dial.png")
# f2 the duckweed data against Loewe and the fitted dial
df=pd.read_csv('data/acidiq.csv'); p,_,_=fit(df,'ct'); p0,_,_=fit(df,'ct',fix=0.0)
fig,ax=plt.subplots(1,2,figsize=(10.5,4.0),sharey=True)
for a_,pc in zip(ax,(33,67)):
    sub=df[df.pct==pc]; a_.semilogx(sub.dose,sub.rgr,"o",color=INK,ms=5,alpha=0.7,label="measured growth rate")
    D=np.logspace(np.log10(sub.dose.min()/1.5),np.log10(sub.dose.max()*1.5),120); grid=pd.DataFrame(dict(dose=D,pct=pc))
    a_.semilogx(D,pred_ct(p0,grid),"--",color=S1,label="Loewe (α = 0), best fit")
    a_.semilogx(D,pred_ct(p,grid),color=S4,label=f"fitted dial α = {p[5]:.2f}")
    a_.set_title(f"Acifluorfen + diquat, mixing ratio {pc}:{100-pc}"); a_.set_xlabel("total dose")
ax[0].set_ylabel("relative growth rate of duckweed"); ax[0].legend(fontsize=9,loc="lower left")
save(fig,"figures/f2_duckweed_fit.png")
