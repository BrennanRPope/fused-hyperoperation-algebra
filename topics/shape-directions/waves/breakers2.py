import sympy as sp, numpy as np, json, itertools
X=sp.symbols('A01 A02 A03 A12 A13 A23'); x,y,z,w=sp.symbols('x y z w')
R=sp.sympify(open('relation4.txt').read())
pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
def general(edges_w):
    return {X[n]:edges_w.get(p,1) for n,p in enumerate(pairs)}
res={}
for name,E in {'path, three weights':{(0,1):x,(1,2):y,(2,3):z},'star, three weights':{(0,1):x,(0,2):y,(0,3):z},
               'triangle + lone, three weights':{(0,1):x,(0,2):y,(1,2):z},'cycle, four weights':{(0,1):x,(1,2):y,(2,3):z,(0,3):w}}.items():
    f=sp.factor(R.subs(general(E))); res[name]=str(f); print(name,':',f)
# tetrahedral configuration
v=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
proj=lambda p:(p[0]+1j*p[1])/(1-p[2])
a=[proj(p) for p in v]; b=[proj(-p) for p in v]
cr=lambda i,j:(a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j]))
tet=[cr(i,j) for i,j in pairs]; print('tetrahedron + antipodes weights',np.round(tet,10)); res['tetrahedral']=[str(np.round(t,10)) for t in tet]

# exact: uniform k, normalized waves (0,inf),(1,1/k),(t,t/k)
t,kk=sp.symbols('t k')
f=lambda T,rj,rk:(1-T)*(rk-rj*T)/((rk-T)*(1-rj*T))
q=sp.factor(sp.numer(sp.together(f(t,kk,kk)-kk)))
print('three-wave condition:',q); res['three_wave_condition']=str(q)
w=sp.Rational(-1,2)+sp.sqrt(3)*sp.I/2
A=[0,1,w,w**2]; B=[sp.oo,sp.Rational(-1,2),-w/2,-w**2/2]
def crx(i,j):
    if i==0: return sp.simplify(A[j]/B[j])
    return sp.simplify((A[i]-A[j])*(B[i]-B[j])/((A[i]-B[j])*(B[i]-A[j])))
cfg=[sp.nsimplify(sp.simplify(sp.expand(crx(i,j)))) for i,j in pairs]
print('k=-2 configuration weights',cfg); res['k_minus_two_config']=[str(c) for c in cfg]
json.dump(res,open('breakers2.json','w'),indent=1)
