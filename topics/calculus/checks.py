"""Checks for 'Calculus from the circle operation' (continuous claims; the discrete ones are in Lean)."""
import sympy as sp, math
x,h,c=sp.symbols('x h c'); n=sp.symbols('n',positive=True,integer=True)
f=sp.Function('f'); g=sp.Function('g')
ma=lambda a,b:a+b+a*b
ms=lambda a,b:(a-b)/(1+b)
ok=[]
def chk(name,cond): ok.append((name,bool(cond))); print(('PASS ' if cond else 'FAIL ')+name)
# test functions
F=sp.sin(x)+x**2/3; G=sp.exp(x/2)-sp.Rational(1,3)
def lim(e): return sp.simplify(sp.limit(e,h,0))
# 1 circle-operation derivative = f'/(1+f) = (ln(1+f))'
chk("D+ f = f'/(1+f)", sp.simplify(lim(ms(F.subs(x,x+h),F)/h) - sp.diff(F,x)/(1+F))==0)
# 2 geometric (product) derivative, shifted by 1
e=((1+F.subs(x,x+h))/(1+F))**(1/h)
val=sp.N(sp.limit(e.subs(x,sp.Rational(1,2)),h,0)); tgt=sp.N(sp.exp((sp.diff(F,x)/(1+F)).subs(x,sp.Rational(1,2))))
chk("geometric derivative = exp(f'/(1+f))", abs(val-tgt)<1e-12)
# 3 both sides lifted (bigeometric shifted): (1+x) f'/(1+f)
chk("bigeometric = (1+x)f'/(1+f)", sp.simplify(lim(ms(F.subs(x,ma(x,h)),F)/h) - (1+x)*sp.diff(F,x)/(1+F))==0)
# 4 Jackson: input stepped by the circle operation, gap h(1+x); in X=1+x it is D_q with q=1+h
X,q=sp.symbols('X q',positive=True); Fx=sp.Function('Phi')
jack=(Fx(X*q)-Fx(X))/((q-1)*X); ours=(Fx((ma(X-1,h))+1)-Fx(X))/(ma(X-1,h)-(X-1))
chk("input circle-operation step = Jackson q-derivative, q=1+h", sp.simplify(jack.subs(q,1+h)-ours)==0)
# 5 product rule: D+(f (+) g) = D+f (+0) D+g i.e. rates add in the limit; exact step uses the circle operation with weight h
p,qq=sp.symbols('p qq')
chk("time-scale product: (1+hp)(1+hq) = 1+h(p+q+h p q)", sp.expand((1+h*p)*(1+h*qq)-(1+h*(p+qq+h*p*qq)))==0)
chk("D+(f(+)g) = D+f + D+g", sp.simplify(sp.diff(ma(F,G),x)/(1+ma(F,G)) - (sp.diff(F,x)/(1+F)+sp.diff(G,x)/(1+G)))==0)
# 6 chain rule: D+(f o g) = (D+f)(g) * g'
inner=x**2+x
chk("chain: D+(F o g) = D+F(g)*g'", sp.simplify(sp.diff(F.subs(x,inner),x)/(1+F.subs(x,inner)) - (sp.diff(F,x)/(1+F)).subs(x,inner)*sp.diff(inner,x))==0)
# 7 product integral: the circle operation n slices of r(x_i)dx -> exp(int r)-1
r=lambda t: math.cos(t)+0.3*t
for N in (10,100,1000,10000):
    acc=0.0; dx=1.0/N
    for i in range(N): acc=acc+r(i*dx)*dx+acc*r(i*dx)*dx
    last=acc
exact=math.exp(math.sin(1)+0.15)-1
chk("product integral of r over [0,1] -> exp(int r)-1 (N=10^4 err %.1e)"%abs(last-exact), abs(last-exact)<1e-3)
# 8 constant slope a/N combined by the circle operation N times -> e^a - 1
a=0.7; N=10**6; chk("(a/N)^(the circle operation N) -> e^a - 1", abs(((1+a/N)**N-1)-(math.exp(a)-1))<1e-5)
# 9 Newton forward with step x/N -> Taylor (sin at x=1.3)
def newton(fun,xv,N):
    hh=xv/N; vals=[fun(k*hh) for k in range(N+1)]; tot=0.0; d=vals[:]; Cnk=1
    for k in range(N+1):
        tot+=Cnk*d[0]; d=[d[i+1]-d[i] for i in range(len(d)-1)]; Cnk=Cnk*(N-k)//(k+1)
    return tot
chk("Newton forward rebuilds sin(1.3) with step 1.3/20", abs(newton(math.sin,1.3,20)-math.sin(1.3))<1e-9)
# layer weights C(N,k) h^k -> x^k/k!
Nn=10**5; xv=1.3
w=[math.comb(Nn,k)*(xv/Nn)**k for k in range(5)]
chk("cube layer C(N,k)(x/N)^k -> x^k/k!", all(abs(w[k]-xv**k/math.factorial(k))<1e-3 for k in range(5)))
# Newton nth difference over step h / h^n -> nth derivative
k=4; hh=1e-2; dvals=[math.exp(0.5+i*hh) for i in range(k+1)]
for _ in range(k): dvals=[dvals[i+1]-dvals[i] for i in range(len(dvals)-1)]
chk("4th difference / h^4 -> 4th derivative of e^x", abs(dvals[0]/hh**k-math.exp(0.5))<0.05)
print(sum(o for _,o in ok),"/",len(ok),"passed")
