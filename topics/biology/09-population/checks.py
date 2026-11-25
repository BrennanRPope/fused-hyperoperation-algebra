import numpy as np, json
from scipy.optimize import brentq
res=[]; out={}
def check(name,ok,info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
star=lambda a,b:a+b+a*b; wadd=lambda k,a,b:a+b+k*a*b; inv=lambda u:-u/(1+u)
def wpow(k,x,n):
    y=0.0
    for _ in range(n): y=wadd(k,y,x)
    return y
# 1. Reed–Frost
p,I=0.002,7
check("Reed–Frost: a susceptible meeting I infectives escapes with (1 − p)^I, so the infection chance is wpow(−1, p, I)",abs((1-(1-p)**I)-wpow(-1,p,I))<1e-15)
N=2000; R0=p*N; S=N-1.0; C=1.0; tot=1.0; curve=[C]
while C>1e-9:
    newC=S*(1-(1-p)**C); S-=newC; C=newC; tot+=C; curve.append(C)
z=brentq(lambda z:1-z-np.exp(-R0*z),1e-6,1)
check("the Reed–Frost epidemic's final size matches the final-size relation 1 − z = e^(−R₀z) for large populations",abs(tot/N-z)<0.01,
      f"R₀ = {R0}: chain {tot/N:.3f}, relation {z:.3f}")
out['rf']=dict(curve=curve,z=z,R0=R0)
# 2. Interventions
e1,e2=0.4,0.3; R=3.0
check("independent interventions that cut transmission by e₁ and e₂ cut it by the 'or' law e₁ ∘₋₁ e₂ together",abs(R*(1-e1)*(1-e2)-R*(1-wadd(-1,e1,e2)))<1e-14,
      f"{e1} and {e2} together: {wadd(-1,e1,e2):.2f}")
check("the herd-immunity threshold 1 − 1/R₀ is −inv(R₀ − 1), the circle operation's inverse of the excess reproduction number",abs((1-1/R)-(-inv(R-1)))<1e-15,
      f"R₀ = {R}: {1-1/R:.3f}")
need=next(n for n in range(1,50) if wpow(-1,0.25,n)>=1-1/R)
check("n independent 25% interventions stop spread at R₀ = 3 once wpow(−1, 0.25, n) ≥ 1 − 1/R₀",need==4,f"n = {need}")
# 3. Beverton–Holt
bh=lambda R,M:(lambda x:R*x/(1+x/M))
R1,M1,R2,M2=2.5,100.,1.6,40.
f=lambda x:bh(R2,M2)(bh(R1,M1)(x)); Rc=R1*R2; Mc=1/(1/M1+R1/M2)
xs=np.linspace(1,300,50)
check("two Beverton–Holt seasons in a row are one Beverton–Holt season: growth excesses compose by ∘₁, and 1/M = 1/M₁ + R₁/M₂",
      np.max(np.abs(f(xs)-bh(Rc,Mc)(xs)))<1e-10 and abs((Rc-1)-star(R1-1,R2-1))<1e-14,f"R = {Rc}, M = {Mc:.2f}")
Rr,Mr=2.0,50.; K=(Rr-1)*Mr; x=5.0; u0=K/x-1; T=6
for _ in range(T): x=bh(Rr,Mr)(x)
check("the room left below carrying capacity, K/N − 1, shrinks by the factor 1/R each generation exactly",abs((K/x-1)-u0*Rr**-T)<1e-10)
# 4. Hassell
hass=lambda lam,a,b:(lambda N:lam*N/(1+a*N)**b)
check("Hassell's N' = λN/(1 + aN)^b has denominator 1 + mpow(aN, b); b = 1 is Beverton–Holt",
      abs(hass(2,0.02,1)(30)-bh(2,50)(30))<1e-14)
Ns=np.linspace(0,200,50); c=0.02
gap=np.max(np.abs(hass(2,c/2000,2000)(Ns)-2*Ns*np.exp(-c*Ns)))
check("with a·b held fixed and b large, Hassell becomes Ricker: (1 + c·N/b)^(−b) → e^(−cN)",gap<0.1,f"largest gap {gap:.3f} at b = 2000")
# 5. Survivorship and colonisation
q=np.array([0.05,0.01,0.02,0.08])
d=0.0
for qi in q: d=wadd(-1,d,qi)
check("the chance of dying across several age classes is the 'or' fold of the per-class mortalities",abs(d-(1-np.prod(1-q)))<1e-15,f"{d:.4f}")
check("Levins colonisation: an empty patch with n occupied neighbours, each colonising with chance c, is colonised with wpow(−1, c, n)",abs(wpow(-1,0.1,5)-(1-0.9**5))<1e-15)
json.dump(out,open('results.json','w'),indent=1)
print(f"\n{sum(res)}/{len(res)} checks passed")
