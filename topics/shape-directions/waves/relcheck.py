import sympy as sp, random, itertools, json
random.seed(9)
A01,A02,A03,A12,A13,A23=X=sp.symbols('A01 A02 A03 A12 A13 A23')
R=sp.sympify(open('relation4.txt').read())
idx={(0,1):A01,(0,2):A02,(0,3):A03,(1,2):A12,(1,3):A13,(2,3):A23}
def cr(a,b,i,j): return (a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j]))
ok=0
for _ in range(30):
    while True:
        a=[sp.Rational(random.randint(-99,99),random.randint(1,13)) for _ in range(4)]; b=[sp.Rational(random.randint(-99,99),random.randint(1,13)) for _ in range(4)]
        if len(set(a+b))==8: break
    sub={idx[k]:cr(a,b,*k) for k in idx}
    ok+= R.subs(sub)==0
print('exact zero at points',ok,'/30')
# symbolic exact check
a=sp.symbols('a0:4'); b=sp.symbols('b0:4')
sub={idx[k]:cr(a,b,*k) for k in idx}
num=sp.numer(sp.together(R.subs(sub)))
print('symbolic', sp.expand(num)==0)
# random arrays off the relation
print('random array value', R.subs({x:sp.Rational(random.randint(1,9),random.randint(1,9)) for x in X}))
# S4 symmetry
sym=True
for perm in itertools.permutations(range(4)):
    m={}
    for (i,j),v in idx.items():
        pi,pj=sorted((perm[i],perm[j])); m[v]=idx[(pi,pj)]
    Rp=R.subs(m,simultaneous=True)
    if sp.expand(Rp-R)!=0 and sp.expand(Rp+R)!=0: sym=False;break
print('S4 symmetric',sym)
U=sp.symbols('u01 u02 u03 u12 u13 u23')
Ru=sp.expand(R.subs({x:1-u for x,u in zip(X,U)},simultaneous=True))
print('terms in u=1-A',len(sp.Add.make_args(Ru)))
print('degrees',sp.Poly(R,*X).degree_list(),'total',sp.Poly(R,*X).total_degree())
# rational second root: given 5, solve quadratic in A23 at a sample
