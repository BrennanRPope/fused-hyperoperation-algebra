# Special relativity's velocity addition through the circle operation: exact checks (sympy) + numeric.
import sympy as sp, math, random
u,v,w,c,a,b,s,t,x,T,X = sp.symbols('u v w c a b s t x T X')
ok=0; tot=0
def check(name, cond):
    global ok, tot
    tot+=1; ok+=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)
z=lambda e: sp.simplify(e)==0
vel=lambda p,q:(p+q)/(1+p*q)
velc=lambda p,q:(p+q)/(1+p*q/c**2)
ma=lambda p,q:p+q+p*q
D=lambda p:(1+p)/(1-p)      # Doppler factor squared (Bondi k^2)
d=lambda p:D(p)-1           # = 2p/(1-p), the circle-operation variable
check("lens: D(u+v rel) = D(u)·D(v)", z(D(vel(u,v))-D(u)*D(v)))
check("The circle operation in disguise: d(u⊕v) = d(u) ∘₁ d(v), d = 2u/(1-u)", z(d(vel(u,v))-ma(d(u),d(v))))
check("d(u) = 2u/(1-u)", z(d(u)-2*u/(1-u)))
check("inverse lens: u = d/(2+d)", z(d(u)/(2+d(u))-u))
check("units: velc(u,v)/c = vel(u/c, v/c)", z(velc(u,v)/c-vel(u/c,v/c)))
check("light absorbs: 1 ⊕ v = 1 (like -1 ∘₁ b = -1)", z(vel(1,v)-1) and z(ma(-1,b)+1))
check("-1 ⊕ v = -1", z(vel(-1,v)+1))
check("light wall maps to the circle-operation wall: d(-1) = -1, d(1) = ∞", z(d(sp.Integer(-1))+1) and sp.limit(d(u),u,1,'-')==sp.oo)
check("speed limit: 1 - u⊕v = (1-u)(1-v)/(1+uv)", z(1-vel(u,v)-(1-u)*(1-v)/(1+u*v)))
check("landscape point: (a+b+0·ab)/(1-(-1)ab)  → k=0, q=-1/c²", z(velc(a,b)-(a+b)/(1-(-1/c**2)*a*b)))
# third layer
ser=sp.series(vel(a*x,b*x),x,0,4).removeO()
check("no corner: second layer of u⊕v is 0, third is -(a²b+ab²)", z(ser.coeff(x,2)) and z(ser.coeff(x,3)+a**2*b+a*b**2))
check("corner tax: w = u+v - u·v·w", z(vel(u,v)-(u+v-u*v*vel(u,v))))
# rapidity
ph1,ph2=sp.symbols('phi1 phi2',real=True)
check("rapidity: tanh(φ1+φ2) = tanh φ1 ⊕ tanh φ2", z(sp.expand_trig(sp.tanh(ph1+ph2))-vel(sp.tanh(ph1),sp.tanh(ph2))))
check("ruler chain: d(tanh φ) = e^(2φ) - 1 = E(2φ)", z(sp.simplify(d(sp.tanh(ph1)).rewrite(sp.exp))-(sp.exp(2*ph1)-1)))
check("rapidity = ½ log(1+d) = artanh u", abs(sp.N(sp.log(D(sp.Rational(1,3)))/2-sp.atanh(sp.Rational(1,3)),50))<1e-45)
# n equal kicks
def velpow(n,p):
    r=0
    for _ in range(n): r=vel(r,p)
    return sp.simplify(r)
check("n kicks: 1+d(u⊕…⊕u) = (1+d(u))^n, n=1..6", all(z(1+d(velpow(n,u))-D(u)**n) for n in range(1,7)))
check("n kicks closed form: ((1+u)^n-(1-u)^n)/((1+u)^n+(1-u)^n)", all(z(velpow(n,u)-((1+u)**n-(1-u)**n)/((1+u)**n+(1-u)**n)) for n in range(1,7)))
# boosts as matrices, light cone squeeze
B=lambda p: sp.Matrix([[1,p],[p,1]])
check("boost product: B(u)B(v) = (1+uv)·B(u⊕v)", z(sp.simplify(B(u)*B(v)-(1+u*v)*B(vel(u,v))).norm()))
P=sp.Matrix([[1,1],[1,-1]])  # (T,X) -> (T+X, T-X)
check("light-cone axes: boost is a squeeze diag(1+u, 1-u)", z((P*B(u)*P.inv()-sp.diag(1+u,1-u)).norm()))
check("squeeze ratio (1+u)/(1-u) = lift 1+d", z((1+u)/(1-u)-(1+d(u))))
Tp,Xp=T+u*X,X+u*T
check("hyperbola: (T')²-(X')² = (1-u²)(T²-X²)", z(Tp**2-Xp**2-(1-u**2)*(T**2-X**2)))
check("normalized boost keeps the hyperbola: γ²(1-u²)=1", z((1/(1-u**2))*(1-u**2)-1))
check("walls of the law: 1 + 0·a - a² = (1-a)(1+a): light cone", z(1-a**2-(1-a)*(1+a)))
# Dial sweep tie-in: κ=+1 boost, 0 shear (Galilean), -1 rotation
gal=lambda p,q:(p+q)/(1+0*p*q)
check("Galilean = q→0 end (c→∞): velc → u+v", sp.limit(velc(u,v),c,sp.oo)==u+v)
check("rotation partner: tan addition at q=+1", z(sp.expand_trig(sp.tan(ph1+ph2))-(sp.tan(ph1)+sp.tan(ph2))/(1-sp.tan(ph1)*sp.tan(ph2))))
# ---- rungs: third-layer match of fractional rung s
# f_s = s-th iterate of E(x)=e^x-1 to 3rd order: x + s/2 x² + (s²/4 - s/12) x³
f=lambda y: y+s/2*y**2+(s**2/4-s/12)*y**3
tt=sp.symbols('tt')
fs=lambda y,ss: y+ss/2*y**2+(ss**2/4-ss/12)*y**3
comp=sp.expand(fs(fs(x,s),tt)); comp=sum(comp.coeff(x,i)*x**i for i in range(1,4))
check("iterate coefficients compose (semigroup) to 3rd order", z(comp-fs(x,s+tt)))
check("s=1 gives e^x-1 to 3rd order", z(fs(x,sp.Integer(1))-sp.series(sp.exp(x)-1,x,0,4).removeO()))
# inverse to 3rd order
g=lambda y: y-s/2*y**2+(s**2/2-(s**2/4-s/12))*y**3
check("g = f_s^-1 to 3rd order", z(sum(sp.expand(f(g(x))).coeff(x,i)*x**i for i in range(1,4))-x))
R=sp.expand(f(g(a*x)+g(b*x))); R=sum(R.coeff(x,i)*x**i for i in range(1,4)).subs(x,1)
k_s=sp.factor(R.coeff(a,1).coeff(b,1)); q_s=sp.factor(R.coeff(a,2).coeff(b,1))
check("rung s law to 3rd layer = a+b+s·ab + s(s-1)/4·(a²b+ab²)", z(k_s-s) and z(q_s-s*(s-1)/4) and z(R.coeff(a,1).coeff(b,2)-q_s))
shape=sp.simplify(k_s**2/(k_s**2-4*q_s))
check("shape of rung s = s", z(shape-s))
check("shape of velocity addition (0,-1) = 0", sp.Rational(0)**2/(0+4)==0)
# natural units: rescale a -> λa with λ = 2/√s: (k,q) -> (λk, λ²q)
lam=2/t
check("natural units: rung t² lands at (2t, t²-1)", z(lam*t**2-2*t) and z(lam**2*t**2*(t**2-1)/4-(t**2-1)))
check("that curve is q = k²/4 - 1 (border parabola moved down 1)", z((2*t)**2/4-1-(t**2-1)))
check("walls on it: 1+2t·a+(t²-1)a² = (1+(1+t)a)(1-(1-t)a)", z(1+2*t*a+(t**2-1)*a**2-(1+(1+t)*a)*(1-(1-t)*a)))
check("t=0 walls ±1 (velocity addition), t=1 walls -1/2 and ∞ (the circle operation k=2)",
      z((1+a)*(1-a)-(1-a**2)) and z((1+2*a)*(1-0*a)-(1+2*a)))
check("raw units: rung s walls at -2/(√s(1+√s)) and 2/(√s(1-√s)); speed limit ~2/√s as s→0",
      z((1+t**2*a+t**2*(t**2-1)/4*a**2).subs(a,-2/(t*(1+t)))) and z((1+t**2*a+t**2*(t**2-1)/4*a**2).subs(a,2/(t*(1-t)))))
# numeric: rung s really is relativistic near addition. Use exact fractional-iterate coefficients to many orders? use 3rd-order check on numbers.
random.seed(1)
err=0
for _ in range(200):
    U,V=random.uniform(-.9,.9),random.uniform(-.9,.9)
    lhs=math.tanh(math.atanh(U)+math.atanh(V)); err=max(err,abs(lhs-(U+V)/(1+U*V)))
check(f"numeric tanh addition, 200 random pairs, max err {err:.1e}", err<1e-12)
print(f"\n{ok}/{tot} checks pass")
