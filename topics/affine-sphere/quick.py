import sympy as sp, time
from blaschke import evaluate
u,v=sp.symbols('u v',real=True)
k,d=sp.Rational(3,5),sp.Rational(7,10)
q=1+k*u+d*(u**2+v**2)
pos=[u/q,(u**2+v**2)/q,v/q]
t=time.time()
X,xi,h,K2=evaluate(pos,[u,v],{u:sp.Rational(1,3),v:sp.Rational(-1,4)})
c=d-k**2/4
center=sp.Matrix([-k,2,0])/(4*d-k**2)
print(time.time()-t, list(xi), K2)
print('xi+H(x-c)', list((xi+2**sp.Rational(3,2)*c*(X-center)).evalf(20)))
X,xi,h,K2=evaluate(pos,[u,v],{u:0,v:0})
print(list(xi), sp.sqrt(2)*(-k/2), sp.sqrt(2).evalf())
