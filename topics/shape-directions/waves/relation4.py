import sympy as sp, random, json, time
A01,A02,A03,A12,A13,A23,a2,a3=sp.symbols('A01 A02 A03 A12 A13 A23 a2 a3')
a={0:None,1:sp.Integer(1),2:a2,3:a3}
b={1:1/A01,2:a2/A02,3:a3/A03}
def cr(i,j): return (a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j]))
e12=sp.numer(sp.together(cr(1,2)-A12)); e13=sp.numer(sp.together(cr(1,3)-A13)); e23=sp.numer(sp.together(cr(2,3)-A23))
print('deg in a2',sp.degree(e12,a2),'deg a3',sp.degree(e13,a3))
t=time.time()
r1=sp.resultant(e23,e12,a2); r1=sp.factor(r1)
print('r1 done',time.time()-t)
r2=sp.resultant(r1,e13,a3); 
F=sp.factor_list(r2)
print('r2 done',time.time()-t)
facs=[(f,m) for f,m in F[1]]
for f,m in facs: print(m, sp.Poly(f).total_degree(), len(sp.Poly(f).terms()), str(f)[:120])
json.dump([str(f) for f,m in facs],open('relation4_factors.json','w'))
