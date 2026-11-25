"""Twists of the photosynthesis oxygen clock in Fused Hyperoperation Algebra.
Run: python3 twists.py  (numpy only)"""
import numpy as np, math, json
a,b=0.084,0.028; g=1-a-b              # miss, double hit, single hit (Shinkarev 2005 fit)
def kok(alpha,beta,gamma,n=4,sign=+1):
    M=np.zeros((n,n))
    for s in range(n):
        M[s,s]+=alpha; M[(s+sign)%n,s]+=gamma; M[(s+2*sign)%n,s]+=beta
    return M
def yields(M,S,release_from,flashes=16):
    # O2 released when the clock passes from release_from back to 0
    out=[];S=S.copy()
    for f in range(flashes):
        y=S[release_from]*(g+b)+S[release_from-1]*b if release_from>0 else 0
        out.append(y); S=M@S
    return np.array(out)
S0=np.array([0.25,0.73,0.0,0.02])
res={}
# 1 forward clock
M=kok(a,b,g); lam=a-b+1j*g
Yf=yields(M,S0,3); res['forward']=dict(lam=[lam.real,lam.imag],abs=abs(lam),deg=math.degrees(math.atan2(lam.imag,lam.real)),yields=list((Yf/Yf.mean()).round(3)))
# 2 reversed direction: S0->S3->S2->S1, O2 released on S1->S0; start mirrored (S_k -> S_-k)
Mr=kok(a,b,g,sign=-1); Sr=S0[[0,3,2,1]]
Yr=[];S=Sr.copy()
for f in range(16):
    Yr.append(S[1]*(g+b)+S[2]*b); S=Mr@S
Yr=np.array(Yr); res['reversed']=dict(lam=[lam.real,-lam.imag],identical_yields=bool(np.allclose(Yr,Yf)),eig=sorted([round(complex(e).imag,4) for e in np.linalg.eigvals(Mr)]))
# 3 circle inverse = un-flashing (time reversal of the probabilities)
Mi=np.linalg.inv(M)
neg=int((Mi<-1e-12).sum()); ev=np.linalg.eigvals(Mi)
S=np.array([0.25,0.25,0.25,0.25]); traj=[]
S=np.array([0.0,1.0,0.0,0.0]); 
for f in range(6):
    traj.append(list(S.round(3))); S=Mi@S
res['inverse']=dict(negative_entries=neg,min_entry=float(Mi.min().round(3)),eig_abs=sorted(np.abs(ev).round(4).tolist()),lam_inv=[(1/lam).real,(1/lam).imag],abs_inv=abs(1/lam),undo_from_S1=traj)
# 4 complement: swap hit and miss
lc=g-b+1j*a; res['complement']=dict(lam=[lc.real,lc.imag],abs=abs(lc),deg=math.degrees(math.atan2(lc.imag,lc.real)),period=360/math.degrees(math.atan2(lc.imag,lc.real)))
# 5 hydroxylamine: two states below S0 -> two-flash delay. Model: start all in "S-1" of a 4-clock = S3 index? do 6-state chain: S-2,S-1 then clock
def nh2oh(start_idx,flashes=16):
    # states: -2,-1,0,1,2,3 ; below 0 they climb without releasing
    P=np.zeros(6); P[start_idx]=1; out=[]
    for f in range(flashes):
        y=P[5]*(g+b)+P[4]*b; out.append(y)
        Q=np.zeros(6)
        for s in range(6):
            Q[s]+=a*P[s]
            def step(k):
                t=s+k
                if t>5: t=2+(t-6)  # wrap to S0
                return t
            Q[step(1)]+=g*P[s]; Q[step(2)]+=b*P[s]
        P=Q
    return np.array(out)
Yn=nh2oh(1)  # start in S-1 (≈ S1 shifted by two)
Y1=nh2oh(3)
res['hydroxylamine']=dict(first_max_untreated=int(np.argmax(Y1[:6])+1),first_max_treated=int(np.argmax(Yn[:8])+1),treated=list((Yn/Yn[4:].mean()).round(3)),untreated=list((Y1/Y1[2:].mean()).round(3)))
# 6 polygon family with the same hit/miss
fam={}  # replaced in twists_output.json by a no-double-hit version (hit = 1 - miss)
for n,name in [(1,'rhodopsin pump: 1 photon per turn'),(2,'H2 or sulfide: 2 electrons'),(4,'O2 from water: 4 electrons'),(8,'N2 (nitrogenase): 8 electrons')]:
    z=a+g*np.exp(2j*np.pi/n)
    fam[n]=dict(name=name,lam=[round(z.real,4),round(z.imag,4)],abs=round(abs(z),4),deg=round(math.degrees(math.atan2(z.imag,z.real)),2))
res['polygons']=fam
# 7 energetics (standard E0' at pH 7, V; Thauer et al. 1977 values)
E=dict(O2_H2O=0.818,NO3_NO2=0.433,cytc=0.254,S_HS=-0.270,NAD=-0.320,NADP=-0.324,H2=-0.414,Fd=-0.420)
F=96.485; hc=1239.84
ph=lambda nm: hc/nm
en={}
en['oxygenic']=dict(lift_V=E['O2_H2O']-E['NADP'],photons_per_e=2,photon_eV=ph(680),stored_frac=(E['O2_H2O']-E['NADP'])/(2*ph(680)))
en['green_sulfur']=dict(lift_V=E['S_HS']-E['Fd'],photons_per_e=1,photon_eV=ph(840),stored_frac=(E['S_HS']-E['Fd'])/ph(840))
en['cyt_oxidase']=dict(drop_V=E['O2_H2O']-E['cytc'],kJ_per_O2=4*F*(E['O2_H2O']-E['cytc']))
en['respiration_NADH_O2']=dict(drop_V=E['O2_H2O']-E['NAD'],kJ_per_O2=4*F*(E['O2_H2O']-E['NAD']))
lift=E['NO3_NO2']-E['NAD']
en['reverse_flow_nitrite']=dict(lift_V=lift,kJ_per_NADH=2*F*lift,protons_per_e_at_0p15=lift/0.15,protons_per_e_at_0p20=lift/0.20)
en['photon_counts_per_CO2']=dict(oxygenic_electrons=4,oxygenic_photons=8,green_sulfur_electrons=4,green_sulfur_photons=4)
en['bacteriorhodopsin']=dict(photon_eV=ph(568),stored_eV_at_0p2V=0.2,frac=0.2/ph(568))
en['electrolysis']=dict(per_e_V=1.229,per_O2_eV=4*1.229)
# OER square: scaling dG(OOH)-dG(OH)=3.2 eV
tot=4*1.229; s23=3.2
steps=[(tot-s23)/2,s23/2,s23/2,(tot-s23)/2]
en['OER']=dict(ideal=[1.229]*4,scaling_best=[round(x,3) for x in steps],eta_min=max(steps)-1.229,arcs_deg_ideal=[90]*4,arcs_deg_real=[round(360*x/tot,1) for x in steps])
res['energy']=en
# 8 mirror sugar: thermodynamics identical
res['mirror']=dict(dG_L_glucose_equals_D=True,note='enantiomers have equal free energy in an achiral setting; parity-violating difference never measured in a molecule')
# leaf circular polarization (Patty et al. 2019): field reflection ~7e-4, lab transmission ~8e-3
res['leaf_CP']=dict(reflection_field=7e-4,transmission_lab=8e-3)
def conv(o):
    if isinstance(o,(np.floating,)): return float(o)
    if isinstance(o,(np.integer,)): return int(o)
    if isinstance(o,complex): return [o.real,o.imag]
    raise TypeError(type(o))
json.dump(res,open('twists_output.json','w'),indent=1,default=conv)
print(json.dumps(res,indent=1,default=conv))
