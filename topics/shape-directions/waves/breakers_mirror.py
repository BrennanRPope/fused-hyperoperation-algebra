# sympy mirror of every statement in lean/RuleBreakers.lean (run from shape-directions/)
import sympy as sp
from fractions import Fraction as F
X=sp.symbols('A01 A02 A03 A12 A13 A23'); R=sp.sympify(open('waves/relation4.txt').read())
x,y,z,w,k,t=sp.symbols('x y z w k t')
ev=lambda *v: sp.expand(R.subs(dict(zip(X,v)),simultaneous=True))
chk={'path':ev(x,1,1,y,1,z)==0,'star':ev(x,y,z,1,1,1)==0,'triangle':ev(x,y,1,z,1,1)==0,'paw':ev(x,y,1,z,1,w)==0,
'cycle':sp.expand(ev(x,1,w,y,1,z)-(x-1)**2*(y-1)**2*(z-1)**2*(w-1)**2)==0,'diamond':sp.expand(ev(k,k,k,k,k,1)-k*(k-1)**8)==0,
'cycle_two':ev(2,1,2,2,1,2)==1,'uniform':sp.expand(ev(k,k,k,k,k,k)+k**3*(k-1)**8*(k+2))==0,
'uniform_three':ev(3,3,3,3,3,3)==-34560,'uniform_minus_two':ev(-2,-2,-2,-2,-2,-2)==0,
'third_wave':sp.expand((1-t)*(k-k*t)-k*((k-t)*(1-k*t))+k*(k-1)*(t*t-(k+1)*t+1))==0,
'fourth_wave':sp.expand((t*t-(k+1)*t+1)-(t*t+t+1)+(k+2)*t)==0}
def mul(x,y): a,b=x;c,d=y; return (a*c-b*d, a*d+b*c-b*d)
sub=lambda x,y:(x[0]-y[0],x[1]-y[1]); add=lambda x,y:(x[0]+y[0],x[1]+y[1]); two=(F(2),F(0))
A={1:(F(1),F(0)),2:(F(0),F(1)),3:(F(-1),F(-1))}; B={1:(F(-1,2),F(0)),2:(F(0),F(-1,2)),3:(F(1,2),F(1,2))}
chk['omega']=add(add(mul((0,1),(0,1)),(0,1)),(1,0))==(0,0)
chk['tetrahedral']=all(add(mul(sub(A[i],A[j]),sub(B[i],B[j])),mul(two,mul(sub(A[i],B[j]),sub(B[i],A[j]))))==(0,0) and mul(sub(A[i],B[j]),sub(B[i],A[j]))!=(0,0) for i,j in [(1,2),(1,3),(2,3)]) and all(add(A[j],mul(two,B[j]))==(0,0) for j in (1,2,3))
for a_,b_ in chk.items(): print(a_,b_)
print('all',all(chk.values()),len(chk))
