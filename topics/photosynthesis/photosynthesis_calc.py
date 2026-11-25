import numpy as np, itertools, math
# Hall & Rao chain (Wikipedia): fractional losses
steps=[("Outside 400-700 nm",0.47),("Not absorbed in leaf",0.30),("Blue/green cooled to red",0.24),("ATP+NADPH to sugar",0.68),("Respiration",0.375)]
eps=np.array([e for _,e in steps])
eta=np.prod(1-eps); print("eta",eta, "total loss",1-eta, "naive sum",eps.sum())
n=len(eps)
# inclusion-exclusion orders
for k in range(1,n+1):
    s=sum(np.prod(c) for c in itertools.combinations(eps,k))
    print("order",k,"sum",round(s,4))
# waterfall attribution in given order vs reversed
def waterfall(order):
    rem=1;out={}
    for i in order:
        out[steps[i][0]]=rem*eps[i]; rem*=1-eps[i]
    return out
print("forward",{k:round(v,3) for k,v in waterfall(range(n)).items()})
print("reverse",{k:round(v,3) for k,v in waterfall(reversed(range(n))).items()})
# log shares
L=-np.log(1-eps); print("nepers",L.round(3),"total",L.sum().round(3),"shares",(L/L.sum()).round(3))
# Shapley: average over orders
sh=np.zeros(n)
for perm in itertools.permutations(range(n)):
    w=waterfall(perm)
    for i in range(n): sh[i]+=w[steps[i][0]]
sh/=math.factorial(n); print("shapley",sh.round(4),sh.sum())
print("log-proportional",(L/L.sum()*(1-eta)).round(4))
# energy
h=6.62607015e-34;c=2.99792458e8;NA=6.02214076e23;eV=1.602176634e-19
E680=h*c/680e-9*NA/1000; print("photon 680 kJ/mol",E680, "eV",h*c/680e-9/eV)
dG=2870; perO2=dG/6
for p in (8,9.4,10): print(p,"photons/O2 eff",perO2/(p*E680))
R=8.314;T=298.15; print("lnK",-dG*1000/(R*T),"log10K",-dG*1000/(R*T)/math.log(10), "photon lnK per photon",E680*1000/(R*T))
# Kok model, Shinkarev classic fit
a,b=0.084,0.028; g=1-a-b
M=np.zeros((4,4))
for s in range(4):
    M[s,s]+=a; M[(s+1)%4,s]+=g; M[(s+2)%4,s]+=b
S=np.array([0.25,0.73,0.0,0.02])
Y=[]
for f in range(1,21):
    # O2 released on S3->S0 (hit) and S2->S0 or S3->S1 via double hit
    y=S[3]*(g+b)+S[2]*b
    Y.append(y); S=M@S
print("yields",np.round(np.array(Y)/np.mean(Y),3))
ev=np.linalg.eigvals(M); print("eig",np.round(ev,4), np.abs(ev).round(4))
lam=a-b+1j*g; print("lambda formula",lam,abs(lam),"angle deg",math.degrees(math.atan2(g,a-b)))
x=lam-1; print("x",x)
# FRET chain loss odds
R0=1.0
