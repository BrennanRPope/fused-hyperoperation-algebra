import sympy as sp, time, json
r1,r2,r3,A12,A13,A23,u,s=sp.symbols('r1 r2 r3 A12 A13 A23 u s')
def E(t,rj,rk,A): return sp.expand((1-t)*(rk-rj*t)-A*(rk-t)*(1-rj*t))
E1=E(u,r1,r2,A12); E2=E(s*u,r1,r3,A13); E3=E(s,r2,r3,A23)
t0=time.time()
R=sp.resultant(sp.Poly(E3,s),sp.Poly(E2,s))
R=R.as_expr() if hasattr(R,'as_expr') else R
print('first',time.time()-t0, sp.Poly(R,u).degree())
R2=sp.resultant(sp.Poly(R,u),sp.Poly(E1,u))
R2=R2.as_expr() if hasattr(R2,'as_expr') else R2
print('second',time.time()-t0)
fl=sp.factor_list(sp.expand(R2))
print('factored',time.time()-t0)
out=[]
for f,m in fl[1]:
    P=sp.Poly(f,r1,r2,r3,A12,A13,A23)
    print(m,P.total_degree(),len(P.terms()),P.degree_list(),str(f)[:150]); out.append(str(f))
json.dump(out,open('relation4_factors.json','w'))
