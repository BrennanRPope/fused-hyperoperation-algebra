"""Sympy check of every ring/field identity stated in lean/ThermoFused.lean."""
import sympy as sp
a,b,c,h,m,k,q,s,r,n,t,u,v,x,y,z,e1,e2,E1,E2,X1,X2,Qh,Qc,u1,u2,i,SA,SB=sp.symbols('a b c h m k q s r n t u v x y z e1 e2 E1 E2 X1 X2 Qh Qc u1 u2 i SA SB')
ma=lambda p,w:p+w+p*w; wa=lambda kk,p,w:p+w+kk*p*w; inv=lambda w:-w/(1+w)
vel=lambda p,w:(p+w)/(1+p*w); lens=lambda w:(1-w)/(1+w); st=lambda p,w:(p+w-2*p*w)/(1-p*w)
I=[("partition_square",a*b-1-ma(a-1,b-1)),
("boson_mode",(1+inv(-z))*(1-z)-1),("fermi num",(E1-1)*(E2+1)+(E1+1)*(E2-1)-2*(E1*E2-1)),
("fermi den",(E1+1)*(E2+1)+(E1-1)*(E2-1)-2*(E1*E2+1)),("vel walls",vel(1,v)-1),("vel -1",vel(-1,v)+1),
("engine_cascade",e1*Qh+e2*(Qh-e1*Qh)-wa(-1,e1,e2)*Qh),("carnot",1-c/h-wa(-1,1-m/h,1-c/m)),
("fridge",u1*Qc+u2*(Qc+u1*Qc)-ma(u1,u2)*Qc),("reverse",x/(1-x)-inv(-x)),("CA half",wa(-1,1-s,1-s)-(1-s*s)),
("ratio_series",1-u*v-wa(-1,1-u,1-v)),("decimation",(X1*X2+1)-(X1+X2)-(X1-1)*(X2-1)),
("ising_series",1-2*wa(-2,u,v)-(1-2*u)*(1-2*v)),("lens_vel",lens(vel(u,v))-lens(u)*lens(v)),
("involution",lens(lens(t))-t),("doppler",lens(x)*((1+x)/(1-x))-1),("stokes add",st(a,b)/(1-st(a,b))-a/(1-a)-b/(1-b)),
("stokes wall",st(1,b)-1),("shields",(n*r/(1+(n-1)*r))/(1-n*r/(1+(n-1)*r))-n*r/(1-r)),
("chain3",(1+m*z+m*z**2+z**3)-(1+z)*(1+(m-1)*z+z**2)),("wick num",(i*a+i*b+k*(i*a)*(i*b)-i*(a+b+i*k*a*b)).subs(i**2,-1)),
("jarzynski/hess",ma(a-1,b-1)-(a*b-1)),("qlog",(x*y-1)/k-wa(k,(x-1)/k,(y-1)/k))]
bad=[nm for nm,ex in I if sp.simplify(sp.expand(ex).subs(i**2,-1))!=0]
# lens fixed point: (1-t)/(1+t)-t vanishes on t^2+2t-1=0
fp=all(sp.simplify(lens(rt)-rt)==0 for rt in sp.solve(t**2+2*t-1,t))
print(f"{len(I)-len(bad)+fp} of {len(I)+1} Lean identities hold symbolically", bad)
