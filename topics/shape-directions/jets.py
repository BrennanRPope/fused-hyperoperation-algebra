import sympy as sp
x,s,a,b=sp.symbols('x s a b')
N=8
E=sp.series(sp.exp(x)-1,x,0,N+1).removeO()
# Julia iterative log V: V(E(x)) = E'(x) V(x), V = sum v_n x^n, v_2=1/2
v=sp.symbols('v2:%d'%(N+1))
V=x**2/2+sum(v[i]*x**(i+2) for i in range(1,N-1))
eq=sp.expand(V.subs(x,E) - sp.diff(E,x)*V)
eq=sp.series(eq,x,0,N+1).removeO()
sol=sp.solve([eq.coeff(x,n) for n in range(3,N+1)],v[1:N-1],dict=True)[0]
V=sp.expand(V.subs(sol))
print("V =",V)
def trunc(e): return sp.series(sp.expand(e),x,0,N+1).removeO()
# E^s = exp(s V d/dx) x
term=x; Es=x
for m in range(1,N+1):
    term=trunc(V*sp.diff(term,x))
    Es+= s**m/sp.factorial(m)*term
Es=sp.expand(trunc(Es))
print("E^s =",[sp.factor(Es.coeff(x,n)) for n in range(1,N+1)])
# check s=1
print("check s=1:",sp.simplify(Es.subs(s,1)-E))
Ems=Es.subs(s,-s)
# rate R_s(a) = 1/ (E^{-s})'(a)
R=sp.series(1/sp.diff(Ems,x),x,0,N).removeO()
R=sp.expand(R)
rc=[sp.factor(R.coeff(x,n)) for n in range(0,N)]
print("rate coeffs r0..:",rc)
import json
json.dump({"Es":[str(sp.factor(Es.coeff(x,n))) for n in range(1,N+1)],"rate":[str(c) for c in rc],"V":str(V)},open("jets.json","w"),indent=1)
