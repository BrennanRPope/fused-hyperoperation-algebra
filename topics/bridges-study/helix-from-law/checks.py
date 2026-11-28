"""DRAFT. Can the law give the Euclidean Kerr helix? Checks K1-K8."""
import sympy as sp, json
M,a,al,v,x,y,k,d,p,q,Q,u = sp.symbols('M a alpha v x y k d p q Q u', real=True)
res = {}
def chk(name, cond, msg):
    res[name] = bool(cond); print(name, "PASS" if cond else "FAIL", " ", msg)
Z = lambda e: sp.simplify(e) == 0

# Two-dial law: conjugate multiplication by f(x) = (1+p x)/(1+q x)
f = lambda z: (1+p*z)/(1+q*z)
law2 = (x+y+k*x*y)/(1-d*x*y)
chk("K1", Z((f(x)*f(y) - f(law2)).subs({k:p+q, d:p*q})),
    "f(x)f(y) = f(x∘y) with x∘y = (x+y+k xy)/(1-d xy), k=p+q, d=pq: associative, one family (Coleman–McGuinness 1991; the repo's d dial)")
z = sp.symbols('z')
lawz = lambda A,B: (A+B+k*A*B)/(1-d*A*B)
chk("K2", Z(lawz(lawz(x,y),z) - lawz(x,lawz(y,z))) and Z(law2.subs(d,0) - (x+y+k*x*y)),
    "associative for every k, d; d = 0 is the original law")
chk("K3", Z(law2.subs({k:0,d:-1}) - (x+y)/(1+x*y)) and Z(law2.subs({k:0,d:1}) - (x+y)/(1-x*y)),
    "k=0, d=-1 is relativistic velocity addition (tanh); k=0, d=+1 is tangent addition (rotation)")

# Kerr horizon read as a spun Schwarzschild horizon
rp = M + sp.sqrt(M**2-a**2)
R = sp.simplify((rp**2+a**2)/rp)
Om = a/(rp**2+a**2); kap = (rp-M)/(rp**2+a**2)
vr = sp.simplify(Om*R)
chk("K4", Z(R-2*M) and Z(vr - a/rp), "equatorial horizon radius is 2M (as Schwarzschild); rim speed v = Ω_H·2M = a/r+")
chk("K5", Z(kap - (1-vr**2)/(4*M)), "surface gravity κ = (1 - v²)/4M: Schwarzschild's 1/4M times the rim's 1/γ²")
helixL = sp.simplify(Om/kap)
chk("K6", Z(helixL - law2.subs({k:0,d:1,x:vr,y:vr})) and Z(helixL - a/sp.sqrt(M**2-a**2)),
    "real time: helix Ω_H/κ = v∘v with k=0, d=+1 (rim composed with itself)")
rpE = M + sp.sqrt(M**2+al**2); t = al/rpE
helixE = al/sp.sqrt(M**2+al**2)
tcheck = all(abs(sp.N((t - sp.tanh(u/2)).subs(al, M*sp.sinh(u)).subs({M:m0,u:u0}))) < 1e-12
             for m0 in (0.7,1.3,2.0) for u0 in (0.2,0.9,2.1))
chk("K7", Z(helixE - law2.subs({k:0,d:-1,x:t,y:t})) and tcheck,
    "imaginary time (a=iα): helix tanh u = t∘t with k=0, d=-1; t = α/r+ = tanh(u/2). The i flips d from +1 to -1")
# multiplicative coordinate read off the moved mass
XE = (1+t)/(1-t)
XL = (1+sp.I*vr)/(1-sp.I*vr)
okE = Z(XE - (rpE+al-M)/M)
okL = Z(sp.simplify(XL - (rp-M+sp.I*a)/M)) and Z(sp.simplify((rp-M)**2 + a**2 - M**2))
# Kerr-Newman: rotation and charge
rpKN = M + sp.sqrt(M**2-a**2-Q**2)
RKN = (rpKN**2+a**2)/rpKN; vKN = a/rpKN
kapKN = (rpKN - (a**2+Q**2)/rpKN)/(2*(rpKN**2+a**2))
okKN = Z(kapKN*2*RKN - (1 - vKN**2 - Q**2/rpKN**2))
chk("K8", okE and okL and okKN,
    "the law's multiplicative coordinate is (r+ + α - M)/M = e^u (imaginary time) and (r+ - M + ia)/M = e^{iχ} (real time); "
    "Kerr-Newman: κ·2R = 1 - v² - Q²/r+², spin and charge slowdowns add with no cross term")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump(res, open("results.json","w"), indent=1)
