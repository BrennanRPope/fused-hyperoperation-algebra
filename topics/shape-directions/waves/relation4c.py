import sympy as sp, time, json
r1,r2,r3,A12,A13,A23,u,s=sp.symbols('r1 r2 r3 A12 A13 A23 u s')
def E(t,rj,rk,A): return sp.expand((1-t)*(rk-rj*t)-A*(rk-t)*(1-rj*t))
E1=E(u,r1,r2,A12); E2=E(s*u,r1,r3,A13); E3=E(s,r2,r3,A23)
g1=(s,u,r1,r2,r3,A12,A13,A23)
t0=time.time()
Rs=sp.Poly(E3,*g1).resultant(sp.Poly(E2,*g1))
print('first',time.time()-t0,Rs.gens,flush=True)
g2=(u,r1,r2,r3,A12,A13,A23)
R2=sp.Poly(Rs.as_expr(),*g2).resultant(sp.Poly(E1,*g2))
print('second',time.time()-t0,flush=True)
fl=R2.factor_list()
print('factor',time.time()-t0,flush=True)
out=[]
for f,m in fl[1]:
    print(m,f.total_degree(),len(f.terms()),f.degree_list(),str(f.as_expr())[:200],flush=True); out.append(str(f.as_expr()))
json.dump(out,open('relation4_factors.json','w'))
