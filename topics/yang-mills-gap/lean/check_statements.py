"""Sympy/numpy check that every Lean statement in YangMillsGap.lean is true (the file is also compiled in the shared library since Nov 26)."""
import sympy as sp, numpy as np
k,d,x,n,t,a,b,c,e=sp.symbols('k d x n t a b c e')
w=lambda p,q: p+q+k*p*q
ok=[]
ok.append(sp.expand(4*d*(1+k*x+d*x*x)-((2*d*x+k)**2-(k*k-4*d)))==0)
ok.append(sp.expand(((n+sp.Rational(1,2))**2)-((-1-n+sp.Rational(1,2))**2))==0)
ok.append(sp.expand((1-t)**2-t**2-(1-2*t))==0)
ok.append(sp.expand(w(w(w(a,b),c),e)-w(w(w(e,c),b),a))==0)
rng=np.random.default_rng(0); A,B=rng.normal(size=(2,3,3)); kk=0.7; I=np.eye(3)
ok.append(np.allclose((I+kk*A)@(I+kk*B)-(I+kk*B)@(I+kk*A), kk*kk*(A@B-B@A)))
ok.append([w(p,q).subs(k,-2) for p,q in [(0,0),(0,1),(1,0),(1,1)]]==[0,1,1,0])
ok.append(sp.Rational(1)*2==sp.Rational(8,3)*sp.Rational(1,2)*sp.Rational(3,2))
for K,D in [(1,1),(0,2),(-3,3)]:
    xs=np.linspace(-50,50,20001); ok.append((1+K*xs+D*xs*xs>0).all())
print(ok, "ALL TRUE" if all(ok) else "SOMETHING FALSE")
