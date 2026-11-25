import sympy as sp, random, json
random.seed(11)
A23,u,s=sp.symbols('A23 u s')
def E(t,rj,rk,A): return sp.expand((1-t)*(rk-rj*t)-A*(rk-t)*(1-rj*t))
def cr(ai,bi,aj,bj): return (ai-aj)*(bi-bj)/((ai-bj)*(bi-aj))
res=[]
for trial in range(5):
    # random points, normalized a0=0,b0=oo,a1=1
    a=[0,1]+[sp.Rational(random.randint(-40,40),random.randint(1,9)) for _ in range(2)]
    b=[None]+[sp.Rational(random.randint(-40,40),random.randint(1,9)) for _ in range(3)]
    r=[None]+[sp.Rational(a[j])/b[j] for j in (1,2,3)]
    A=lambda i,j: cr(a[i],b[i],a[j],b[j])
    E1=E(u,r[1],r[2],A(1,2)); E2=E(s*u,r[1],r[3],A(1,3)); E3=E(s,r[2],r[3],A23)
    R=sp.resultant(E3,E2,s); R2=sp.factor(sp.resultant(R,E1,u))
    P=sp.Poly(sp.numer(sp.together(R2)),A23)
    fl=sp.factor_list(P.as_expr())
    true=A(2,3)
    # which factor vanishes at the true value
    facs=[(sp.Poly(f,A23).degree(),f.subs(A23,true)==0) for f,m in fl[1] if f.has(A23)]
    roots=set(sp.Poly(f,A23).nroots().__iter__().__next__() for f,m in fl[1] if f.has(A23))
    distinct=[sp.N(x,8) for f,m in fl[1] if f.has(A23) for x in sp.Poly(f,A23).nroots()]
    res.append(dict(deg=P.degree(),factors=facs,true=str(true),candidates=[str(x) for x in distinct]))
    print(res[-1])
json.dump(res,open('sixth.json','w'),indent=1)
