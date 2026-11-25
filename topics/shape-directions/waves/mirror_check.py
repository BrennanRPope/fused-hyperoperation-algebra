import sympy as sp
p1,p2,q1,q2,ai,bi,aj,bj,l,m,g=sp.symbols('p1 p2 q1 q2 ai bi aj bj l m g')
P=lambda kx,ky,kt: kx**4-4*kx*kt+3*ky**2
k=lambda p,q:(p-q,p**2-q**2,p**3-q**3)
k1,k2=k(p1,q1),k(p2,q2)
N=(p1-p2)*(q1-q2); D=(p1-q2)*(q1-p2)
chk={}
chk['dispersion']=sp.expand(P(*k(p1,q1)))==0
chk['two_soliton']=sp.expand(P(*[a+b for a,b in zip(k1,k2)])*N+P(*[a-b for a,b in zip(k1,k2)])*D)==0
chk['free_iff_trivial']=sp.expand((ai-aj)*(bi-bj)-(ai-bj)*(bi-aj)-(ai-bi)*(aj-bj))==0
f=lambda x:l*x+m
chk['affine_num']=sp.expand((f(ai)-f(aj))*(f(bi)-f(bj))-l**2*(ai-aj)*(bi-bj))==0
chk['affine_den']=sp.expand((f(ai)-f(bj))*(f(bi)-f(aj))-l**2*(ai-bj)*(bi-aj))==0
# inversion: cr(1/a..) = cr(a..): cleared cross-multiplication
inv=lambda x:1/x
lhs=(inv(ai)-inv(aj))*(inv(bi)-inv(bj))*((ai-bj)*(bi-aj)); rhs=((ai-aj)*(bi-bj))*(inv(ai)-inv(bj))*(inv(bi)-inv(aj))
chk['inversion']=sp.simplify(lhs-rhs)==0
chk['kdv']=sp.expand((p1-p2)*(-p1+p2)+(p1-p2)**2)==0 and sp.expand((p1+p2)*(-p1-p2)+(p1+p2)**2)==0
chk['excess_points']=sp.expand(g*(g-1)-2*(2*g-3)-(g-2)*(g-3))==0
chk['excess_curves']=sp.expand(g*(g+1)-2*(3*g-3)-(g-2)*(g-3))==0
chk['excess_kdv']=sp.expand(g*(g-1)-2*(g-1)-(g-1)*(g-2))==0
for k_,v in chk.items(): print(k_,v)
A12,A13,A14,A23,A24,A34,a,b,c,d=sp.symbols('A12 A13 A14 A23 A24 A34 a b c d')
def tau4(A12,A13,A14,A23,A24,A34,a,b,c,d):
    return (1+a+b+c+d+A12*a*b+A13*a*c+A14*a*d+A23*b*c+A24*b*d+A34*c*d
      +A12*A13*A23*a*b*c+A12*A14*A24*a*b*d+A13*A14*A34*a*c*d+A23*A24*A34*b*c*d
      +A12*A13*A14*A23*A24*A34*a*b*c*d)
import itertools
W={(0,1):A12,(0,2):A13,(0,3):A14,(1,2):A23,(1,3):A24,(2,3):A34}; X=[a,b,c,d]
gen=0
for S in itertools.product([0,1],repeat=4):
    idx=[i for i in range(4) if S[i]]; w=1
    for i in idx: w*=X[i]
    for i,j in itertools.combinations(idx,2): w*=W[(i,j)]
    gen+=w
print('tau4_def', sp.expand(gen-tau4(A12,A13,A14,A23,A24,A34,a,b,c,d))==0)
print('tau4_free', sp.expand(tau4(1,1,1,1,1,1,a,b,c,d)-(1+a)*(1+b)*(1+c)*(1+d))==0)
print('tau4_splits', sp.expand(tau4(A12,A13,1,A23,1,1,a,b,c,d)-(1+a+b+c+A12*a*b+A13*a*c+A23*b*c+A12*A13*A23*a*b*c)*(1+d))==0)
print('tau4_exclusion', sp.expand(tau4(0,A13,A14,A23,A24,A34,a,b,c,d)-(1+a+b+c+d+A13*a*c+A14*a*d+A23*b*c+A24*b*d+A34*c*d+A13*A14*A34*a*c*d+A23*A24*A34*b*c*d))==0)
print('one_wave', sp.expand((1+a)*(1+b)-(1+a+b+a*b))==0)
print('excess_from_four', all(((max(g-2,0))*(max(g-3,0))>0)==(g>=4) for g in range(2,40)))
print('kdv_from_three', all(((g-1)*max(g-2,0)>0)==(g>=3) for g in range(2,40)))
X=sp.symbols('A01 A02 A03 A12 A13 A23'); R=sp.sympify(open('waves/relation4.txt').read())
def crq(ai,bi,aj,bj): return sp.Rational((ai-aj)*(bi-bj),(ai-bj)*(bi-aj))
c1=[crq(0,2,1,5),crq(0,2,3,-1),crq(0,2,7,4),crq(1,5,3,-1),crq(1,5,7,4),crq(3,-1,7,4)]
c2=[crq(1,-3,2,6),crq(1,-3,4,0),crq(1,-3,9,5),crq(2,6,4,0),crq(2,6,9,5),crq(4,0,9,5)]
print('config1', c1==[sp.Rational(-3,5),9,sp.Rational(7,10),-3,-1,sp.Rational(5,2)], R.subs(dict(zip(X,c1)))==0)
print('config2', c2==[sp.Rational(9,25),sp.Rational(-9,7),sp.Rational(4,3),-3,sp.Rational(-7,9),sp.Rational(25,9)], R.subs(dict(zip(X,c2)))==0)
print('all2', R.subs({x:2 for x in X})==-32)
print('second', R.subs(dict(zip(X,c1[:5]+[sp.Rational(92,77)])))==0)
