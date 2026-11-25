"""Where the q-dial sits relative to the rung flow B_s(a,b) = E^s(E^-s(a) + E^-s(b)), E(x) = e^x - 1.
Exact symbolic (sympy). Result: B_s = a + b + s·ab + c3(s)(a²b + ab²) + ..., while the q-dial law is
a + b + k·ab exactly. So k = s to second order (the dial is the flow's tangent direction), and the
two agree again at s = 1 (k = 1): the dial is a straight chord, the flow a curved arc."""
import sympy as sp
s,x,a,b,k=sp.symbols('s x a b k')
ORD=5
E=sp.series(sp.exp(x)-1,x,0,ORD+1).removeO()
def comp(f,g):
    return sp.expand(sp.series(f.subs(x,g),x,0,ORD+1).removeO())
its=[x]
for n in range(1,ORD+1): its.append(comp(E,its[-1]))
Es=0
for i in range(1,ORD+1):
    pts=list(range(i)); vals=[its[n].coeff(x,i) for n in pts]
    Es+=sp.expand(sp.interpolate(list(zip(pts,vals)),s))*x**i if i>1 else x
Es=sp.expand(Es)
Ems=sp.expand(Es.subs(s,-s))
print('E^s =',sp.collect(Es,x))
inner=sp.expand(Ems.subs(x,a)+Ems.subs(x,b))
t=sp.Symbol('t')
B=sp.expand(Es.subs(x,inner.subs({a:a*t,b:b*t})))
B=sum(B.coeff(t,d) for d in range(1,4))
P=sp.Poly(B,a,b)
for (i,j),c in sorted(P.terms()):
    if i+j<=3: print(f'a^{i} b^{j}:',sp.factor(c))
c2=P.coeff_monomial(a*b); c3=P.coeff_monomial(a**2*b)
assert sp.simplify(c2-s)==0
assert sp.simplify(c3.subs(s,0))==0 and sp.simplify(c3.subs(s,1))==0
print('c2 = s: the dial k matches flow time s to second order')
print('c3 =',sp.factor(c3),' (zero at s=0 and s=1 only; the q-dial has none)')
# chord vs arc: the gap at s=1/2
print('c3(1/2) =',c3.subs(s,sp.Rational(1,2)))
# q-numbers along the flow: E^s(n E^-s(x))/x, series in x
for n in [2,3]:
    qn=sp.expand(Es.subs(x,n*Ems))
    qn=sum(qn.coeff(x,d)*x**(d-1) for d in range(1,4))
    print(f'[{n}] along the flow:',sp.factor(qn.coeff(x,0)),'+',sp.factor(qn.coeff(x,1)),'x +',sp.factor(qn.coeff(x,2)),'x^2')
    print('   at s=1 vs [n]_q with q=1+x:',sp.expand(((1+x)**n-1)/x))
