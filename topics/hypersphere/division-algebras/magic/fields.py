"""Verify exact static solutions of S = (1/16pi) int sqrt(-g)[R - 2(dphi)^2 - e^{-2a phi} F^2] (G=1)."""
import sympy as sp
t,r,th,ph=sp.symbols('t r theta phi_',positive=True)
X=[t,r,th,ph]
def einstein(g):
    gi=g.inv(); n=4
    Gam=[[[sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(n))/2 for c in range(n)] for b in range(n)] for a in range(n)]
    def Ric(b,c):
        return sp.simplify(sum(sp.diff(Gam[a][b][c],X[a])-sp.diff(Gam[a][b][a],X[c])+sum(Gam[a][a][d]*Gam[d][b][c]-Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
    R=sp.Matrix(4,4,lambda b,c:Ric(b,c)); Rs=sp.simplify(sum(gi[i,j]*R[i,j] for i in range(n) for j in range(n)))
    return sp.simplify(R-Rs*g/2), gi
def check(g,phi,At,a):
    G,gi=einstein(g)
    F=sp.zeros(4); F[0,1]=-sp.diff(At,r); F[1,0]=-F[0,1]   # F_tr = d_t A_r - d_r A_t
    Fup=gi*F*gi; F2=sp.simplify(sum(F[i,j]*Fup[i,j] for i in range(4) for j in range(4)))
    dphi=[sp.diff(phi,x) for x in X]; dphi2=sum(gi[i,j]*dphi[i]*dphi[j] for i in range(4) for j in range(4))
    T=sp.Matrix(4,4,lambda m,n: 2*(dphi[m]*dphi[n]-g[m,n]*dphi2/2)+2*sp.exp(-2*a*phi)*(sum(F[m,k]*F[n,l]*gi[k,l] for k in range(4) for l in range(4))-g[m,n]*F2/4))
    sq=sp.sqrt(-g.det())
    box=sp.simplify(sum(sp.diff(sq*gi[i,j]*dphi[j],X[i]) for i in range(4) for j in range(4))/sq)
    dil=sp.simplify(box+a*sp.exp(-2*a*phi)*F2/2)
    maxw=[sp.simplify(sum(sp.diff(sq*sp.exp(-2*a*phi)*Fup[m,nn],X[m]) for m in range(4))) for nn in range(4)]
    num=lambda e: max(abs(complex(sp.N(e.subs({r:rv,th:0.7}),30))) for rv in (sp.Rational(37,7),sp.Rational(91,5)))
    return [num(x) for x in (G-T)], num(dil), [num(x) for x in maxw]
def ghs(a,rp,rm):
    a=sp.S(a)
    s=a**2/(1+a**2); lam2=(1-rp/r)*(1-rm/r)**((1-a**2)/(1+a**2)); R2=r**2*(1-rm/r)**(2*s)
    phi=a/(1+a**2)*sp.log(1-rm/r)   # e^{2a phi}=(1-rm/r)^{2a^2/(1+a^2)}, so e^{2a phi}/R^2 = 1/r^2
    Q=sp.sqrt(rp*rm/(1+a**2)); M=rp/2+(1-a**2)/(1+a**2)*rm/2
    At=Q/r   # since F_tr = e^{2a phi} Q / R^2  and  e^{2a phi}/R^2 = 1/r^2
    g=sp.diag(-lam2,1/lam2,R2,R2*sp.sin(th)**2)
    return g,phi,At,Q,M
def jnw(gam,b):
    f=1-b/r; g=sp.diag(-f**gam,f**(-gam),f**(1-gam)*r**2,f**(1-gam)*r**2*sp.sin(th)**2)
    phi=sp.sqrt((1-gam**2)/4)*sp.log(f); return g,phi,sp.Integer(0),gam*b/2,sp.sqrt(1-gam**2)*b/2
if __name__=="__main__":
    rp,rm=sp.Rational(3),sp.Rational(1)
    for a in [sp.Integer(0),1/sp.sqrt(3),sp.Integer(1),sp.sqrt(3)]:
        g,phi,At,Q,M=ghs(a,rp,rm); E,d,mx=check(g,phi,At,a)
        print("GHS a=",a," Einstein ok:",max(E)<1e-12," dilaton ok:",d<1e-12," Maxwell ok:",max(mx)<1e-12)
    g,phi,At,M,S=jnw(1/sp.sqrt(2),sp.Integer(2)); E,d,mx=check(g,phi,At,0)
    print("JNW gamma=1/sqrt2 Einstein ok:",max(E)<1e-12," dilaton ok:",d<1e-12)
