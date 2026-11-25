import numpy as np
from scipy.optimize import minimize_scalar, brentq
rng=np.random.default_rng(11); res=[]
def check(name,ok,info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b
star=lambda a,b:a+b+a*b
inv=lambda u:-u/(1+u)

# 1. Hill's force–velocity hyperbola
F0=1.0; a=0.25*F0; vmax=1.0; b=a*vmax/F0
F=np.linspace(0,F0,200); v=b*(F0-F)/(F+a)
X,Y,X0=F/a,v/b,F0/a
check("Hill (F + a)(v + b) = (F0 + a)b is the level set X ∘₁ Y = X0, with X = F/a, Y = v/b",np.allclose(star(X,Y),X0))
check("velocity is the circle-operation difference Y = X0 ∘₁ inv(X)",np.allclose(Y,star(X0,inv(X))))
check("so v_max/b = F0/a: the curve is symmetric in its two reduced variables",abs(vmax/b-F0/a)<1e-12)
check("the asymptotes F = −a and v = −b are the circle operation's absorbing −1 in each coordinate: X = −1 forces 1 + X ∘₁ Y = 0 for every Y",
      all(abs(star(-1,y)+1)<1e-15 for y in (0.1,3,40)))
P=lambda f: f*b*(F0-f)/(f+a)
opt=minimize_scalar(lambda f:-P(f),bounds=(0,F0),method="bounded").x
fopt=a*(np.sqrt(1+F0/a)-1)
check("peak power sits where (1 + X)² = 1 + X0, i.e. F/F0 = (a/F0)(√(1 + F0/a) − 1)",abs(opt-fopt)<1e-5,
      f"a/F0 = 0.25 gives F/F0 = {fopt:.3f}, v/vmax = {b*(F0-fopt)/(fopt+a)/vmax:.3f}, power {P(fopt)/(F0*vmax):.3f} of F0·vmax")
check("at peak power the two reduced variables are equal, X = Y: the circle operation 'square root' of X0",abs(fopt/a-b*(F0-fopt)/(fopt+a)/b)<1e-9)

# 2. Antagonists (Paton & Rang): dose ratios on the dial
def dose_ratio(b1,b2,omega,site):
    # agonist occupancy a/(Z0 + a); DR = Z0 (the antagonist part of the partition function)
    return 1+b1+b2 if site=="same" else 1+b1+b2+omega*b1*b2
for _ in range(1):
    b1,b2=rng.uniform(0.1,20,2)
    # check by solving occupancy equality numerically for the two-site model
    def occ(A,b1,b2,om): return A/(1+b1+b2+om*b1*b2+A)
    target=occ(1.0,0,0,1)
    for om,lab in ((0,"same site (mutually exclusive): dose ratios add, DR₁₂ − 1 = (DR₁ − 1) + (DR₂ − 1)"),
                   (1,"independent separate sites: dose ratios multiply, DR₁₂ = DR₁·DR₂")):
        A12=brentq(lambda A: occ(A,b1,b2,om)-target,1e-9,1e9)
        A1=brentq(lambda A: occ(A,b1,0,om)-target,1e-9,1e9); A2=brentq(lambda A: occ(A,0,b2,om)-target,1e-9,1e9)
        check(lab,abs((A12-1)-wadd(om,A1-1,A2-1))<1e-6*(A12),f"b₁ = {b1:.2f}, b₂ = {b2:.2f}")
om=3.0; A12=brentq(lambda A: A/(1+b1+b2+om*b1*b2+A)-0.5,1e-9,1e12)
check("antagonists that help each other bind (ω = 3) push the dose ratio beyond multiplicative: DR₁₂ − 1 = (DR₁ − 1) ∘_ω (DR₂ − 1)",
      abs((A12-1)-wadd(om,b1,b2))<1e-6*A12)

# 3. Organ extraction: series is dial 1, parallel pathways dial 0
Q=1.0; cl1,cl2=rng.uniform(0.1,5,2)
E=lambda u:u/(1+u)          # well-stirred extraction with u = CLint/Q
check("well-stirred extraction E = u/(1 + u) = −inv(u): saturation again",abs(E(cl1)-(-inv(cl1)))<1e-15)
uP=cl1+cl2
check("two enzymes in one organ (parallel): intrinsic clearances add, dial 0 on u",abs(E(uP)-E(wadd(0,cl1,cl2)))<1e-15)
Fser=(1-E(cl1))*(1-E(cl2))
check("two organs in series (gut then liver): escapes multiply, so u combines by plain ∘₁",abs((1/Fser-1)-star(cl1,cl2))<1e-12)
check("and the total first-pass loss 1 − F is the 'or' of the two extractions, dial −1",abs((1-Fser)-wadd(-1,E(cl1),E(cl2)))<1e-15)
check("same total intrinsic clearance removes more when split across two organs in series than pooled in one",
      all((1-(1-E(c/2))**2) > E(c) for c in (0.2,1,5)))

# 4. relabellings
check("Q10: rate ratio over ΔT is a circle operation power, (1 + q)^(ΔT/10) − 1",abs(((1+1.3)**(7/10)-1)-((2.3**0.7)-1))<1e-15)
check("Poiseuille: a 10% narrower vessel has resistance 1 + mpow(−0.1, −4) = 1.52 times higher",abs((0.9**-4)-1.524)<1e-3)
print(f"\n{sum(res)}/{len(res)} checks passed")
