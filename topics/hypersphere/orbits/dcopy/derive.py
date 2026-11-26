import sympy as sp
M,r,E,pr,L,a,b,c,x=sp.symbols('M r E p_r L a b c x',positive=True)
phi=2*M/r
# Kerr-Schild Schwarzschild: g^{mu nu}=eta - phi k k, k^mu=(-1, rhat); k.p = E + p_r
ginv = sp.Matrix([[-1-phi,phi],[phi,1-phi]])  # (t,r) block, inverse metric
g = sp.Matrix([[-(1-phi),phi],[phi,1+phi]])     # KS metric (t,r) block: eta + phi k k, k_mu=(1,1)
print("ginv*g =",sp.simplify(ginv*g))
p=sp.Matrix([-E,pr])
shell=sp.expand((p.T*ginv*p)[0]+L**2/r**2+1)
print("shell:",shell)
print("matches -E^2+pr^2+L^2/r^2+1-phi(E+pr)^2:",sp.simplify(shell-(-E**2+pr**2+L**2/r**2+1-phi*(E+pr)**2))==0)
# tagged version -> equivalent static metric after completing square
Aeff=1/(1+a*phi+b**2*phi**2/(1-c*phi)); Beff=1/(1-c*phi)
print("A(1,1,1)=",sp.simplify(Aeff.subs({a:1,b:1,c:1})))
# exact T for static metric: T=(N-4v^2)/B, N=3+rA''/A', v^2 = r A'/(2A - rA') (local circular speed^2 = rA/(2A))
def T_of(A,B):
    Ap=sp.diff(A,r);App=sp.diff(A,r,2)
    N=3+r*App/Ap
    v2=r*Ap/(2*A)
    return sp.simplify((N-4*v2)/B)
for tag in [(1,1,1),(1,1,0),(1,0,1),(1,0,0),(0,1,0),(0,1,1)]:
    T=T_of(Aeff.subs({a:tag[0],b:tag[1],c:tag[2]}),Beff.subs({c:tag[2]}))
    Tx=sp.simplify(T.subs(M,x*r))
    print(tag,"T=",sp.factor(Tx),"  weak:",sp.series(Tx,x,0,2).removeO())
Tg=sp.simplify(T_of(Aeff,Beff).subs(M,x*r))
print("general weak T:",sp.expand(sp.series(Tg,x,0,2).removeO()))
