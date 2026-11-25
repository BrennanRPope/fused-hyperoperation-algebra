"""lean_mirror.py: every statement of lean/ShapeDirections.lean checked with sympy (exact) before compile."""
import sympy as sp
from fractions import Fraction as Fr
k,l,a,b,q,p,c,s,x,m,qq,A,w1,w2,X_,Y_=sp.symbols('k l a b q p c s x m qq A w1 w2 X Y')
W=lambda K,u,v:u+v+K*u*v; M=lambda u,v:u+v+u*v
ok=[]
def t(name,e): ok.append((name,sp.expand(e)==0))
t('weight_is_units', l*W(l*k,a,b)-W(k,l*a,l*b))
t('weight_to_star', k*W(k,a,b)-M(k*a,k*b))
t('qcount_is_weight', W(q-1,a,b)-(a+b+q*a*b-a*b))
z=sp.symbols('z'); path3=1+3*z+(2*k+1)*z*z+k*k*z**3
t('path3_sum', path3-(1+z+z+z+z*z*k+z*z*k+z*z+z**3*k*k))
t('ruler_num', a*(1+c*b)+b*(1+c*a)+p*a*b-W(p+2*c,a,b))
t('ruler_den', (1+c*a)*(1+c*b)-q*a*b-((1-(q+p*c+c*c)*a*b)+c*W(p+2*c,a,b)))
t('ruler_keeps_disc',(p+2*c)**2-4*(q+p*c+c*c)-(p*p-4*q))
t('zoom_disc',(l*p)**2-4*(l*l*q)-l*l*(p*p-4*q))
# jets
def comp4(f,g): return f[0]*(g[0]**2+2*g[1])+3*f[1]*g[0]+f[2]+g[2]
def comp5(f,g): return f[0]*(2*g[0]*g[1]+2*g[2])+f[1]*(3*g[0]**2+3*g[1])+4*f[2]*g[0]+f[3]+g[3]
E=[sp.Rational(1,2),sp.Rational(1,6),sp.Rational(1,24),sp.Rational(1,120)]
c5=sp.symbols('c5'); C=[s/2,s*(3*s-1)/12,s*(2*s-1)*(3*s-1)/48,c5]
t('ladder_jets_commute4',comp4(E,C)-comp4(C,E)); t('ladder_jets_commute5',comp5(E,C)-comp5(C,E))
c3,c4=sp.symbols('c3 c4'); G=[s/2,c3,c4,c5]
sol=sp.solve([comp4(E,G)-comp4(G,E),comp5(E,G)-comp5(G,E)],[c3,c4])
t('ladder_forces_c3',sol[c3]-s*(3*s-1)/12); t('ladder_forces_c4',sol[c4]-s*(2*s-1)*(3*s-1)/48)
l2,l3,l4=-s/2,s*(3*s+1)/12,-(s*(2*s+1)*(3*s+1))/48; r1,r2,r3=s,s*(s-1)/4,-(s*(s-1))/12
t('rung_rate_inverts1',r1+2*l2); t('rung_rate_inverts2',r2+2*l2*r1+3*l3); t('rung_rate_inverts3',r3+2*l2*r2+3*l3*r1+4*l4)
t('lens_jets_at_minus2',(-s)*(3*(-s)-1)/12-s*(3*s+1)/12); t('lens_jets_at_minus3',(-s)*(2*(-s)-1)*(3*(-s)-1)/48+(s*(2*s+1)*(3*s+1))/48)
t('rung_disc_is_height',s*s-4*(s*(s-1)/4)-s)
ok.append(('rung_leaves_plane', all(-(v*(v-1))/12!=0 for v in [Fr(2),Fr(-1),Fr(1,2)])))
ok.append(('rung_crossings', -(0*(0-1))/12==0 and -(2*0-1)/Fr(12)==Fr(1,12) and -(2*1-1)/Fr(12)==Fr(-1,12)))
det3=lambda a1,a2,a3,b1,b2,b3,c1,c2,c3_: a1*(b2*c3_-b3*c2)-a2*(b1*c3_-b3*c1)+a3*(b1*c2-b2*c1)
ok.append(('jet_rank', det3(12,0,0,0,12,0,12,3,-1)==-144 and det3(12,0,0,12,0,0,0,12,0)==0))
t('glue_commutes',M(M(a,q-1),b)-M(M(a,b),q-1)); t('glue_twice',M(M(a,q-1),qq-1)-M(a,q*qq-1))
y_=sp.symbols('y_'); t('glue_descends',M(M(a,x-1),M(b,y_-1))-M(M(a,b),x*y_-1))
psi3=lambda a1,a2,xx:3*xx**4+(a1*a1+4*a2)*xx**3
t('cusp_psi3',psi3(0,0,x)-3*x**4); t('node_psi3',psi3(1,0,x)-x**3*(3*x+1))
# general psi3 with a3=a4=a6=0: b4=b6=b8=0
ok.append(('node_three_points', Fr(-1,3)**2-4*(-(Fr(-1,3))**3)==Fr(-1,27)))
# node equation at x=-1/3: y^2 + x y - x^3 = y^2 - y/3 + 1/27 : disc = 1/9 - 4/27
ok.append(('node_disc_meaning', Fr(1,9)-Fr(4,27)==Fr(-1,27)))
ok.append(('periods_outnumber_curves', all(((6*g-6<g*(g+1))==(4<=g)) for g in range(2,40))))
def wpow(K,u,n):
    r=0
    for _ in range(n): r=W(K,r,u)
    return sp.expand(r)
t('pseries_two',wpow(k,x,2)-(2*x+k*x*x)); t('pseries_three',wpow(k,x,3)-(3*x+3*k*x*x+k*k*x**3))
t('pseries_five',wpow(k,x,5)-(5*x+10*k*x*x+10*k*k*x**3+5*k**3*x**4+k**4*x**5))
t('additive_at_two',wpow(2*m,x,2)-2*(x+m*x*x)); t('additive_at_three',wpow(3*m,x,3)-3*(x+3*m*x*x+3*m*m*x**3))
t('additive_at_five',wpow(5*m,x,5)-5*(x+10*m*x*x+50*m*m*x**3+125*m**3*x**4+125*m**4*x**5))
t('mult_at_three',wpow(k,x,3)-k*k*x**3-3*(x+k*x*x)); t('mult_at_two',wpow(k,x,2)-k*x*x-2*x)
t('fermat_three_a',(3*m+1)**2-1-3*(3*m*m+2*m)); t('fermat_three_b',(3*m+2)**2-1-3*(3*m*m+4*m+1))
# pair weight: tau2 identity sanity
for n,v in ok: print(("ok   " if v else "FAIL ")+n)
print(sum(v for _,v in ok),"of",len(ok))
