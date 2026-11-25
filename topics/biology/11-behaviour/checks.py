import numpy as np, json
res=[]; out={}
def check(name,ok,info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
star=lambda a,b:a+b+a*b; wadd=lambda k,a,b:a+b+k*a*b
def wpow(k,x,n):
    y=0.0
    for _ in range(n): y=wadd(k,y,x)
    return y
# 1. Rescorla–Wagner
lam,a=1.0,0.15; V=0.0; Vs=[0.0]
for n in range(30): V+=a*(lam-V); Vs.append(V)
check("Rescorla–Wagner learning ΔV = αβ(λ − V) gives V_n = λ·wpow(−1, αβ, n): the learning curve is the 'or' power",
      max(abs(Vs[n]-lam*wpow(-1,a,n)) for n in range(31))<1e-14,f"after 10 trials V = {Vs[10]:.3f}")
check("two training phases of n₁ and n₂ trials compose by the 'or' law: wpow(−1, a, n₁ + n₂) = wpow(−1, a, n₁) ∘₋₁ wpow(−1, a, n₂)",
      abs(wpow(-1,a,17)-wadd(-1,wpow(-1,a,9),wpow(-1,a,8)))<1e-15)
a1,a2=0.1,0.2; VA=VB=0.0
for n in range(12):
    err=lam-VA-VB; VA+=a1*err; VB+=a2*err
check("a compound of two cues learns its summed strength as one cue with rate α₁ + α₂ (shared prediction error)",
      abs((VA+VB)-lam*wpow(-1,a1+a2,12))<1e-14)
# blocking
VA=0.0
for n in range(40): VA+=a1*(lam-VA)
VB=0.0
for n in range(20):
    err=lam-VA-VB; VA+=a1*err; VB+=a2*err
check("blocking: a cue pretrained near λ leaves almost nothing for a new cue added later",VB<0.02,f"new cue reaches {VB:.3f} of λ")
out['rw']=Vs
# 2. Discounting
D1,D2,k=3.0,5.0,0.4
u_h=lambda D:k*D
check("hyperbolic discounting d = 1/(1 + kD): the 'interest' u = 1/d − 1 = kD adds over split delays, dial 0",abs(u_h(D1+D2)-(u_h(D1)+u_h(D2)))<1e-15)
r=0.1; u_e=lambda D:(1+r)**D-1
check("exponential discounting: u composes by plain ∘₁ over split delays, dial 1, the only time-consistent case",abs(u_e(D1+D2)-star(u_e(D1),u_e(D2)))<1e-14)
dh=lambda D:1/(1+k*D)
gap=1/(dh(D1)*dh(D2))-1/dh(D1+D2)
check("subadditivity: discounting a hyperbolic delay in two pieces discounts more than in one, by exactly k²D₁D₂ in 1/d",abs(gap-k*k*D1*D2)<1e-14,
      f"k = {k}, delays {D1:.0f} + {D2:.0f}: 1/d grows by {gap:.2f}")
kap,c=0.3,0.12; u_k=lambda D:(np.exp(c*D)-1)/kap
check("the discount family u(D) = (e^(cD) − 1)/κ composes over split delays by the fixed dial κ: exponential measured from the wall −1/κ",
      abs(u_k(D1+D2)-wadd(kap,u_k(D1),u_k(D2)))<1e-12)
kk=0.4; Ds=np.linspace(0,20,50)
lim=max(abs((np.exp(1e-6*kk*D)-1)/1e-6-kk*D) for D in Ds)
check("…κ = 1 is exponential and κ → 0 with c/κ = k held fixed is hyperbolic",lim<1e-4)
al,be=0.5,0.2
gh=lambda D,al,be:(1+al*D)**(-be/al)
check("Loewenstein–Prelec's generalised hyperbola (1 + αD)^(−β/α) is 1/(1 + mpow(αD, β/α)); α = β is hyperbolic and α → 0 exponential",
      abs(gh(4,be,be)-1/(1+be*4))<1e-15 and abs(gh(4,1e-7,be)-np.exp(-be*4))<1e-6)
# preference reversal under hyperbolic
small,large,gapD=50,100,6; kpr=0.25
pref=lambda t:(small/(1+kpr*t))>(large/(1+kpr*(t+gapD)))
check("hyperbolic discounting reverses preferences as both options recede, exponential never does",
      pref(0) and not pref(30) and all(((small*(1+r)**-t)>(large*(1+r)**-(t+gapD)))==((small)>(large*(1+r)**-gapD)) for t in (0,30)))
# 3. Many eyes
p=0.2
check("many eyes: a group of n vigilant animals spots a predator with wpow(−1, p, n)",abs(wpow(-1,p,6)-(1-0.8**6))<1e-15,f"p = {p}, n = 6: {wpow(-1,p,6):.3f}")
# 4. Herrnstein
re=2.0; rr=np.array([1,5,20.]); inv=lambda u:-u/(1+u)
check("Herrnstein's hyperbola R = k·r/(r + r_e) is k times −inv(r/r_e), saturation again",np.allclose(rr/(rr+re),-inv(rr/re)))
json.dump(out,open('results.json','w'),indent=1)
print(f"\n{sum(res)}/{len(res)} checks passed")
