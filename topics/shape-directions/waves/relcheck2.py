import sympy as sp, itertools, json, random
X=sp.symbols('A01 A02 A03 A12 A13 A23'); A01,A02,A03,A12,A13,A23=X
R=sp.Poly(sp.sympify(open('relation4.txt').read()),*X)
b1,a2,b2,a3,b3=sp.symbols('b1 a2 b2 a3 b3')
a={1:1,2:a2,3:a3}; b={1:b1,2:b2,3:b3}
def cr(i,j): return (a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j]))
sub={A01:sp.Integer(1)/b1,A02:a2/b2,A03:a3/b3,A12:cr(1,2),A13:cr(1,3),A23:cr(2,3)}  # wave 0 = (0, infinity)
# evaluate with common denominators via Poly composition: multiply out
num=sp.numer(sp.together(R.as_expr().subs(sub)))
print('normalized symbolic identity:', sp.expand(num)==0)
res={'normalized_identity':bool(sp.expand(num)==0)}
# S4 symmetry up to sign
idx={(0,1):A01,(0,2):A02,(0,3):A03,(1,2):A12,(1,3):A13,(2,3):A23}
Re=R.as_expr(); sym=True
for perm in itertools.permutations(range(4)):
    m={v:idx[tuple(sorted((perm[i],perm[j])))] for (i,j),v in idx.items()}
    Rp=sp.expand(Re.subs(m,simultaneous=True))
    if Rp!=sp.expand(Re) and Rp!=sp.expand(-Re): sym=False;break
print('S4 symmetric',sym); res['S4']=sym
print('degrees',R.degree_list(),'total',R.total_degree(),'terms',len(R.terms())); res['degrees']=R.degree_list(); res['terms']=len(R.terms())
print('irreducible',len(sp.factor_list(Re)[1])==1)
json.dump(res,open('relcheck.json','w'))
