"""Checks: the velocity-dependent second-order (2PN) terms of Gupta's Lagrangian that are linear in one coupling
(the charge-charge block at G = 0, and the G-linear gravity and scalar block) against the time-symmetric Fokker
exchange of a photon, a graviton and a scalar between two worldlines in flat space. At linear order in a coupling the
Lagrangian is fixed up to a total time derivative, so the test is that the Euler-Lagrange expression of the difference
vanishes identically (evaluated on random jets). Negative controls: a small change to each block must be detected.
Writes checks_kinetic_output.txt."""
import sympy as sp, random, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pn2
t,t1,t2,e=sp.symbols('t t1 t2 e')
D=3
X1=[sp.Function(f'x1{i}') for i in range(D)]; X2=[sp.Function(f'x2{i}') for i in range(D)]
def fok(F0,F2,F4):
    """order-4 part of time-symmetric exchange: sum_n (-1)^n/(2n)! d1^n d2^n [F r^(2n-1)]; F's given as functions of (v1,v2) lists"""
    x1=[f(t1) for f in X1]; x2=[f(t2) for f in X2]
    v1=[sp.diff(c,t1) for c in x1]; v2=[sp.diff(c,t2) for c in x2]
    r=sp.sqrt(sum((a-b)**2 for a,b in zip(x1,x2)))
    ex = F4(v1,v2)/r - sp.Rational(1,2)*sp.diff(F2(v1,v2)*r,t1,t2) + sp.Rational(1,24)*sp.diff(F0(v1,v2)*r**3,t1,2,t2,2)
    return ex.subs({t1:t,t2:t})
dot=lambda a,b: sum(x*y for x,y in zip(a,b))
# Gupta blocks from K symbols
x1=[f(t) for f in X1]; x2=[f(t) for f in X2]
v1=[sp.diff(c,t) for c in x1]; v2=[sp.diff(c,t) for c in x2]; a1=[sp.diff(c,t,2) for c in x1]; a2=[sp.diff(c,t,2) for c in x2]
R=[a-b for a,b in zip(x1,x2)]; r=sp.sqrt(dot(R,R)); n=[c/r for c in R]
Kv={'v1v1':dot(v1,v1),'v2v2':dot(v2,v2),'v1v2':dot(v1,v2),'nv1':dot(n,v1),'nv2':dot(n,v2),'na1':dot(n,a1),'na2':dot(n,a2),
    'a1a2':dot(a1,a2),'v1a2':dot(v1,a2),'a1v2':dot(a1,v2),'v1a1':dot(v1,a1),'v2a2':dot(v2,a2)}
l1,l2=pn2.halves(); A2,B2=pn2.image(l2); full2=sp.expand(A2+B2)
Gs,q1,q2,m1,m2,G12,gb=pn2.G,pn2.q1,pn2.q2,pn2.m1,pn2.m2,pn2.G12,pn2.gb
def block(expr):
    return expr.subs({pn2.K[k]:v for k,v in Kv.items()}).subs(pn2.r,r)
qq=full2.subs({Gs:0,G12:0}); qq=sp.expand(qq - qq.subs({q1:0}))   # q1 q2 block at G=0
qq=sp.expand(qq.coeff(q1).coeff(q2))
gl=sp.expand(full2.subs({q1:0,q2:0}).coeff(G12,1))      # G12-linear block
gl=sp.expand(gl.coeff(m1,1).coeff(m2,1))
b0=gl.subs(gb,0); bg=sp.expand(gl.coeff(gb,1)); assert sp.expand(gl-b0-gb*bg)==0
print('qq block:',qq); print('gl block free:',gl.free_symbols)
FEM=(lambda a,b:1, lambda a,b:-dot(a,b), lambda a,b:0)
FS=(lambda a,b:1, lambda a,b:-(dot(a,a)+dot(b,b))/2, lambda a,b:-(dot(a,a)**2+dot(b,b)**2)/8+dot(a,a)*dot(b,b)/4)
def FGR(k):
    def f(a,b):
        A,B,C=sp.symbols('A B C')  # A=v1^2 B=v2^2 C=v1v2
        F=(2*(1-e**2*C)**2-(1-e**2*A)*(1-e**2*B))/sp.sqrt((1-e**2*A)*(1-e**2*B))
        c=sp.series(F,e,0,6).removeO().coeff(e,k)
        return c.subs({A:dot(a,a),B:dot(b,b),C:dot(a,b)})
    return f
FG=(FGR(0),FGR(2),FGR(4))
tests={'EM (q1q2, G=0)':(block(qq), -fok(*FEM)),
       'graviton (G-linear, no scalar)':(block(b0), fok(*FG)),
       'scalar (G-linear, coefficient of alpha1 alpha2)':(block(b0-2*bg), fok(*FS))}
from sympy.calculus.euler import euler_equations
random.seed(3); lines=[]; ok_all=True
def elmax(Dl, ntr=2):
    eqs=euler_equations(Dl,[f(t) for f in X1+X2],t); vals=[]
    for trial in range(ntr):
        jet={}
        for f in X1+X2:
            for k in range(5,-1,-1):
                jet[sp.diff(f(t),t,k) if k else f(t)]=sp.Rational(random.randint(-9,9),random.randint(1,5))+(3 if (k==0 and f in X1[:1]) else 0)
        vals.append(max(abs(sp.nsimplify(eq.lhs.subs(jet))) for eq in eqs))
    return max(vals)
bump=Kv['v1v1']*Kv['v2v2']/r/100
tot=sp.diff(Kv['nv1']*Kv['v2v2']*r, t)          # a total derivative: must pass
for name,(Lg,Lf) in tests.items():
    d0=elmax(Lg-Lf); d1=elmax(Lg-Lf+bump); d2=elmax(Lg-Lf+tot)
    ok = d0==0 and d1!=0 and d2==0; ok_all&=ok
    lines.append(f"{'PASS' if ok else 'FAIL'}  {name}: Gupta's block minus Fokker exchange is a total derivative (EL max {d0}); "
                 f"control with v1^2 v2^2/(100 r) added detected (EL max {sp.N(d1,4)}); control with a total derivative added passes (EL max {d2})")
    print(lines[-1], flush=True)
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
open(os.path.join(HERE,'checks_kinetic_output.txt'),'w').write("\n".join(lines)+"\n")
