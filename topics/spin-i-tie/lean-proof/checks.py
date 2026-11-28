# DRAFT. Round 15: symbolic check (sympy) of each step in CoinRule.lean, since the Lean file could not be compiled here.
import json, sympy as sp
k,a,b,c,q=sp.symbols('k a b c q',real=True)
law=lambda k,a,b: a+b+k*a*b; res={}
res['L1 lift']=sp.simplify(1+k*law(k,a,b)-(1+k*a)*(1+k*b))==0
res['L2 comm+assoc']=sp.simplify(law(k,a,b)-law(k,b,a))==0 and sp.simplify(law(k,law(k,a,b),c)-law(k,a,law(k,b,c)))==0
p=sp.symbols('p0:2',positive=True); r=sp.symbols('r0:3',positive=True)
S=lambda xs: (1-sum(x**q for x in xs))/(q-1)
joint=[x*y for x in p for y in r]
lhs=sp.expand_power_base(S(joint),force=True)
res['L3 Tsallis composes by the law, k = 1 - q']=sp.simplify(lhs-law(1-q,S(p),S(r)))==0
sA,sB=sp.symbols('sA sB',nonzero=True)
res['L4 coin rule forces k = 0']=sp.solve(sp.Eq(law(k,sA,sB)-sA,sB),k)==[0]
coin=S([sp.Rational(1,2)]*2)
res['L5 a fair coin has nonzero Tsallis size for q = 2, 3, 1/2']=all(coin.subs(q,v)!=0 for v in (2,3,sp.Rational(1,2)))
for n,v in res.items(): print(n,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({n:bool(v) for n,v in res.items()},open('results.json','w'),indent=1)
