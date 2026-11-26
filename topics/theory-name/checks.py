import sympy as sp
x,a,b,k,d=sp.symbols('x a b k d')
law=lambda u,v:(u+v+k*u*v)/(1-d*u*v)
# 1. left translation M_a(x)=law(a,x) composes: M_a(M_b(x)) == M_{a∘b}(x)
lhs=law(a,law(b,x)); rhs=law(law(a,b),x)
print("1 composition of maps closes (all k,d):", sp.simplify(lhs-rhs)==0)
# 2. at d=0, M_a is a homothety about -1/k with ratio 1+ka
h=law(a,x).subs(d,0)
print("2 d=0 map is affine, ratio 1+ka:", sp.simplify(sp.diff(h,x)-(1+k*a))==0, " fixes -1/k:", sp.simplify(h.subs(x,-1/k)+1/k)==0)
# 3. M_a is Mobius for d!=0 (degree 1/1 in x) with fixed points the walls d x^2 + k x + 1 = 0
num,den=sp.fraction(sp.together(law(a,x)-x))
print("3 fixed points of M_a:", sp.factor(num))
# 4. homotheties with different centres don't commute: gap is a translation
c1,c2,r1,r2=sp.symbols('c1 c2 r1 r2')
H=lambda c,r,y: c+r*(y-c)
gap=sp.simplify(H(c1,r1,H(c2,r2,x))-H(c2,r2,H(c1,r1,x)))
print("4 order gap for different centres:", sp.factor(gap))
# Titeica's xyz=1 as the affine sphere over the octant: cited from Loftin, Survey on Affine Spheres (arXiv 0809.1186), not re-derived here.
