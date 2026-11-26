"""Gravity's 6 through the classical double copy (Kerr-Schild form).
Schwarzschild in Kerr-Schild form: g = eta + phi k k, phi = 2M/r, k null.  Its inverse is exact,
g^-1 = eta - phi k k, so a probe's mass shell is  eta(p,p) - phi (k.p)^2 = -1.
The single copy (Coulomb, same k):                 eta(p,p) - 2 (Q/r) (k.p) = -1.
Double copy of the probe: the charge coupling Q (k.p) becomes M (k.p)^2  (Gonzo-Shi: charge -> momentum).
(k.p)^2 = E^2 + 2 E p_r + p_r^2.  Tag the three terms a, b, c and measure what each does to the orbit."""
import sympy as sp, mpmath as mp, json
res={}; ok_all=True
def check(name,ok,msg):
    global ok_all; ok_all&=bool(ok); res[name]=dict(ok=bool(ok),msg=msg); print(("PASS " if ok else "FAIL ")+name+": "+msg)

M,r,E,pr,L,a,b,c,x=sp.symbols('M r E p_r L a b c x',positive=True)
phi=2*M/r
g=sp.Matrix([[-(1-phi),phi],[phi,1+phi]]); ginv=sp.Matrix([[-1-phi,phi],[phi,1-phi]])
check("K1",sp.simplify(g*ginv)==sp.eye(2),"Kerr-Schild inverse is exact: (eta + phi kk)^-1 = eta - phi kk (k null)")
p=sp.Matrix([-E,pr]); kp=E+pr
shell=sp.expand((p.T*ginv*p)[0]+L**2/r**2+1)
check("K2",sp.simplify(shell-(-E**2+pr**2+L**2/r**2+1-phi*kp**2))==0,
      "probe mass shell = flat shell - (2M/r)(k.p)^2, with k.p = E + p_r: the single copy's (Q/r)(k.p) squared")
check("K3",sp.expand(kp**2)==E**2+2*E*pr+pr**2,"(k.p)^2 = E^2 [a: Newton's pull] + 2E p_r [b: cross term] + p_r^2 [c: space]")

# tagged shell -> complete the square in p_r (a canonical shift, orbits unchanged) -> static A, B
Aeff=1/(1+a*phi+b**2*phi**2/(1-c*phi)); Beff=1/(1-c*phi)
def T_of(A,B):
    Ap=sp.diff(A,r);App=sp.diff(A,r,2); N=3+r*App/Ap; v2=r*Ap/(2*A)
    return sp.simplify(((N-4*v2)/B).subs(M,x*r))
T100=T_of(Aeff.subs({a:1,b:0,c:0}),1); T101=T_of(Aeff.subs({a:1,b:0,c:1}),Beff.subs(c,1))
T110=T_of(Aeff.subs({a:1,b:1,c:0}),1); T111=T_of(Aeff.subs({a:1,b:1,c:1}),Beff.subs(c,1))
v2=x/(1-2*x)
check("K4",sp.simplify(T100-1)==0,"Newton's term alone (a): T = 1 exactly, at every field strength: no precession, orbits close")
check("K5",sp.simplify(T101-(1-2*x))==0,"Newton + space (a + c): T = 1 - 2M/r exactly: the space factor, same as the gyroscope's")
check("K6",sp.simplify(T111-(1-2*x)*(1-4*v2))==0 and sp.simplify(T111-(1-6*x))==0,
      "all three (a + b + c): T = 1 - 6M/r = (1 - 2M/r)(1 - (2v)^2) exactly")
check("K6b",sp.simplify(T110-1/(1+4*x))==0,"Newton + cross term, no space (a + b): T = 1/(1 + 4M/r) exactly")
Tg=T_of(Aeff,Beff); weak=sp.expand(sp.series(Tg,x,0,2).removeO())
check("K7",sp.simplify(weak-(1-(4*b**2/a+2*c)*x))==0,
      "first order: T = 1 - (4b^2/a + 2c) M/r, i.e. advance per lap = pi (4b^2 + 2ac) M^2/L^2; at a=b=c=1: 4 + 2 = 6")

# real orbit integrations
mp.mp.dps=40
def advance(W,L,u1,u2g,Eg):
    E0=mp.re(mp.findroot(lambda E: W(u1,E),Eg)); u2=mp.re(mp.findroot(lambda u: W(u,E0),u2g(E0) if callable(u2g) else u2g))
    lo,hi=min(u1,u2),max(u1,u2); mid=(lo+hi)/2; half=(hi-lo)/2
    def f(th):
        w=abs(W(mid-half*mp.cos(th),E0))
        return L*half*mp.sin(th)/mp.sqrt(w) if w>0 else L*half/mp.sqrt(L**2*half**2)  # exact-zero guard at an endpoint (limit value for a quadratic W)
    return float(2*mp.quad(f,[0,mp.pi])-2*mp.pi)
def grav(a_,b_,c_,M_=1,L_=200,e=0.3):
    pp=mp.mpf(L_)**2/(a_*M_)
    W=lambda u,E:(1-2*c_*M_*u)*((1+2*a_*M_*u)*E**2+4*b_**2*M_**2*u**2*E**2/(1-2*c_*M_*u)-1-L_**2*u**2)
    return advance(W,L_,(1+e)/pp,(1-e)/pp,1-a_*M_/(2*pp))
sims={}; worst=0
for t in [(1,1,1),(1,1,0),(1,0,1),(1,.5,1),(1,.5,0),(2,1,1)]:
    d=grav(*t); pred=float(mp.pi*(4*t[1]**2+2*t[0]*t[2])/200**2); sims[str(t)]=[d,pred]; worst=max(worst,abs(d/pred-1))
d0=grav(1,0,0); d0s=advance(lambda u,E:(1+2*u)*E**2-1-64*u**2,8,mp.mpf(1.5)/64,lambda E0:2*E0**2/64-mp.mpf(1.5)/64,0.99)   # strong field: r ~ 43-128M
sims["(1,0,0)"]=[d0,0.0]; sims["(1,0,0) strong"]=[d0s,0.0]
check("K8",worst<1e-3 and abs(d0)<1e-15 and abs(d0s)<1e-15,
      "simulated orbits (e = 0.3) match pi(4b^2+2ac)M^2/L^2 to %.1e for 6 tag sets; Newton's term alone precesses %.0e (weak) and %.0e (r ~ 43-128M)"%(worst,d0,d0s))

# single copy (electric, Kerr-Schild gauge): -E^2 - 2kuE*a_ - 2ku p_r*b_ + p_r^2 + L^2u^2 + 1 = 0
def elec(a_,b_,k=1,L_=200,e=0.3):
    pp=mp.mpf(L_)**2/(a_*k)
    W=lambda u,E: E**2+2*a_*k*u*E+b_**2*k**2*u**2-1-L_**2*u**2
    return advance(W,L_,(1+e)/pp,(1-e)/pp,1-a_*k/(2*pp))
e11=elec(1,1); e10=elec(1,0); som=float(mp.pi/200**2)
sims["electric (1,1)"]=[e11,som]; sims["electric (1,0)"]=[e10,0.0]
check("K9",abs(e11/som-1)<1e-3 and abs(e10)<1e-15,
      "single copy: Coulomb term alone precesses %.0e; the cross term gives Sommerfeld's pi k^2/L^2 (ratio %.5f): electric 1 = cross term at strength k"%(e10,e11/som))
check("K10",abs(sims[str((1,1,0))][0]/(4*som)-1)<1e-3,
      "gravity's cross term = the single copy's, at strength 2M instead of M (the 2 of 2E p_r): (2M)^2 -> 4, so 1 + 3 = 2^2 (ratio %.5f)"%(sims[str((1,1,0))][0]/(4*som)))
json.dump(dict(checks=res,sims=sims),open('results_dcopy.json','w'),indent=1)
print("ALL PASS" if ok_all else "SOME FAIL", "%d/%d"%(sum(v['ok'] for v in res.values()),len(res)))
