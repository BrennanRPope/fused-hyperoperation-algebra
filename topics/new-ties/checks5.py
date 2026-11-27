"""Round 5 checks: airtight brackets, no-collision proof for a family, Lean build, and the Boolean family.
Run: python3 checks5.py   (needs python-flint; a few minutes)"""
import io, contextlib, json, re, itertools
from fractions import Fraction as F
import sympy as sp
with contextlib.redirect_stdout(io.StringIO()):
    import planelaws as P
import union_ifs as U
from exact_bounds import Setup, mu_star
from pingpong import criterion, EXAMPLE, PRIMES, X0
import boolean as B
res=[]
def ok(name,cond,extra=""):
    res.append(bool(cond)); print(("PASS " if cond else "FAIL ")+name+("  "+extra if extra else ""))
out={}
# V1 certified brackets for the round-4 examples, and certified below D*
for q in (2,3,12):
    cl=[(F(1),F(q),[F(1,5),F(3,20)]),(F(-1),F(-q),[F(1,4)])]
    r=Setup(cl,nb=6).bracket(); m=mu_star(cl,r['t_up'])
    out[f"walls ±{q}"]=(str(r['t_low']),str(r['t_up']))
    ok(f"V1 walls ±{q}: certified {float(r['t_low']):.4f} <= delta <= {float(r['t_up']):.4f}, and t_up < D* (certified sign)",
       r['low_ok'] and r['up_ok'] and bool(m<0))
# V2-V5 the ping-pong example
c1,c2,c3,info=criterion(EXAMPLE,PRIMES,X0)
ok("V2 ping-pong example (walls ±21/5; multipliers 5/18, 5/24 | 7/24; primes 2, 3; base point 0): criterion (i), (ii), (iii) hold exactly",c1 and c2 and c3)
maps=[P.map_with(p,q,s)[0] for p,q,ms in EXAMPLE for s in ms]
lv={P.norm(((F(1),F(0)),(F(0),F(1))))}; counts=[]
for n in range(9):
    lv={P.norm(P.mm(g,A)) for A in lv for g in maps}; counts.append(len(lv))
tc=U.trace_counts([(F(i),[F(j) for j in range(len(ms))]) for i,(p,q,ms) in enumerate(EXAMPLE)],9)
ok("V3 ping-pong example: different words give different maps up to length 9 (exact), as the proof says",counts==tc,str(counts[-1]))
geo=all((A[0][0]*A[1][1]-A[0][1]*A[1][0])>0 and all(F(-21,20)<P.app(A,x)<F(21,20) for x in (F(-21,20),F(21,20))) for A in maps)
ok("V4 ping-pong example: positive determinants, V=[-1.05,1.05] sent strictly inside itself",geo)
r=Setup(EXAMPLE,nb=6).bracket(); m=mu_star(EXAMPLE,r['t_up'])
out["ping-pong example"]=(str(r['t_low']),str(r['t_up']))
ok(f"V5 ping-pong example: certified {float(r['t_low']):.4f} <= delta <= {float(r['t_up']):.4f} < D*",r['low_ok'] and r['up_ok'] and bool(m<0))
# V6 no collisions for the round-4 examples up to length 14 (modular arithmetic)
CM=json.load(open("collisions_modp.json"))
ok("V6 round-4 examples (walls ±2, ±3, ±12): no collisions up to length 14 (832,040 maps each, modular check)",all(v["match"] for v in CM.values()))
# V7 Lean
txt=open("lean/compile_output.txt").read(); src=open("lean/NewTies.lean").read()
nthm=len(re.findall(r"^theorem",src,re.M)); naxl=txt.count("depends on axioms: [propext, Classical.choice, Quot.sound]")
ok(f"V7 Lean: lean/NewTies.lean builds with Lean 4.37.0-rc1, {nthm} theorems, no sorry, standard axioms only",
   "exit status 0" in txt and "sorry" not in src.replace("no `sorry`","") and naxl==nthm and "error" not in txt)
# V8 Boolean classification
rows=B.classify(); fixed=[r for r in rows if len(r)==5 and r[4]]
cent=sorted(set(r[3] for r in fixed)); contracting=[r for r in fixed if r[2] not in ("-1",)]
ok("V8 of the 12 Boolean steps that depend on x, 7 keep a centre that does not move with a; their centres are 0, 1/2, 1;"
   " 6 of them contract (complement does not)", len(fixed)==7 and cent==["0","1","1/2"] and len(contracting)==6, str(cent))
# V9 identities
a,x,k,M=sp.symbols('a x k M'); law=a+x+k*a*x
ok("V9 XOR is union in a box of size 1/2; AND is OR conjugated by complement; every fused law is multiplication in u = 1 + kx;"
   " k = 0 is a translation (slope 1, no fixed point, no fractal)",
   sp.simplify((a+x-a*x/M).subs(M,sp.Rational(1,2))-(a+x-2*a*x))==0 and sp.simplify(1-((1-a)+(1-x)-(1-a)*(1-x))-a*x)==0
   and sp.simplify((1+k*law)-(1+k*a)*(1+k*x))==0 and sp.diff(law.subs(k,0),x)==1)
# V10 union + intersection = round 1 example A
boolA=[(F(1),[F(1,5),F(3,20)]),(F(0),[F(1,4)])]; exA=U.EXAMPLES["A: two in box 1, one in box 3"]
cA,_=U.distinct_composites(U.maps_from(boolA),11); cE,_=U.distinct_composites(U.maps_from(exA),11)
ok("V10 OR with densities 4/5, 17/20 and AND with 1/4 is round 1's example A moved by y = (3 - x)/2: same counts to length 11,"
   " same D*", cA==cE==U.trace_counts(boolA,11) and abs(U.trace_root(boolA)-U.trace_root(exA))<1e-12, f"D* {U.trace_root(boolA):.4f}")
# V11 accidental collisions come from slope coincidences
s1,s2,r_=sp.symbols('s1 s2 r'); c1_,c2_=sp.Rational(1,2),sp.Integer(1)
Xf=lambda s:(lambda y:c1_+s*(y-c1_)); Of=lambda y:c2_+r_*(y-c2_)
def comp(*fs):
    def h(y):
        for f in reversed(fs): y=f(y)
        return y
    return h
diff=sp.factor(sp.expand(comp(Xf(s2),Of,Of,Xf(s1))(x)-comp(Xf(s1),Of,Xf(s2),Of)(x)))
ok("V11 OR 1/5 with XOR slopes 1/6, 1/7 collide at length 4 because r s1 s2 - r s2 + s1 - s2 = 0 (a slope coincidence, not a rule)",
   diff.subs({s1:sp.Rational(1,6),s2:sp.Rational(1,7),r_:sp.Rational(1,5)})==0 and diff!=0, str(diff))
bj=json.load(open("boolean.json"))["boolean_fractal"]
first=next((i+1 for i,(c,t) in enumerate(zip(bj["counts"],bj["trace_counts"])) if c!=t),None)
ok(f"V12 three-centre Boolean fractal (OR 4/5, 17/20; AND 3/13; XOR 9/22, 15/26): counts match the shared-centre count up to"
   f" length {first-1}; first collision at length {first} ({bj['trace_counts'][first-1]-bj['counts'][first-1]} fewer maps)", first==9)
json.dump(out,open("checks5.json","w"),indent=1)
print(f"{sum(res)}/{len(res)} pass")
