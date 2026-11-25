"""Fused Hyperoperation Algebra: numeric checks of the rung identities.
Rung n operation = ordinary hyperoperation H_n conjugated by x -> 1+x:
  a (n) b = H_n(1+a, 1+b) - 1
  rung 1: a+b+1   rung 2: a+b+ab (mult-add)   rung 3: (1+a)^(1+b)-1   rung 4: (1+a)^^(1+b)-1
"""
import numpy as np, math, itertools
np.seterr(over="ignore")
rng = np.random.default_rng(1)
r1 = lambda a,b: a+b+1
r2 = lambda a,b: a+b+a*b
r3 = lambda a,b: (1+a)**(1+b)-1
L  = lambda a: np.log1p(a)          # log(1+a)
LL = lambda a: np.log(np.log1p(a))  # log log(1+a), needs a>0
ok = []
def check(name, cond): ok.append(cond); print(f"{'PASS' if cond else 'FAIL'}  {name}")
A = rng.uniform(0.1,3,1000); B = rng.uniform(0.1,3,1000); C = rng.uniform(0.1,3,1000)
check("rung2 = conj of multiply", np.allclose(r2(A,B), (1+A)*(1+B)-1))
check("rung3 right identity 0: a (3) 0 = a", np.allclose(r3(A,0*A), A))
check("rung3 left absorber 0: 0 (3) b = 0", np.allclose(r3(0*A,B), 0))
check("rung3 carries mult-add: a (3) b = a (+) ((1+a)^b - 1)", np.allclose(r3(A,B), r2(A,(1+A)**B-1)))
check("rung3 right-distributes over rung2: (a(+)b)(3)c = (a(3)c)(+)(b(3)c)", np.allclose(r3(r2(A,B),C), r2(r3(A,C),r3(B,C))))
check("rung3 tower law: (a(3)b)(3)c = a(3)(b(+)c)", np.allclose(r3(r3(A,B),C), r3(A,r2(B,C))))
check("rung3 sum law: a(3)(b(1)c) = (a(3)b)(+)(a(3)c)", np.allclose(r3(A,r1(B,C)), r2(r3(A,B),r3(A,C))))
check("L linearizes rung2: L(a(+)b) = L(a)+L(b)", np.allclose(L(r2(A,B)), L(A)+L(B)))
check("LL linearizes rung3: LL(a(3)b) = LL(a)+L(b)", np.allclose(LL(r3(A,B)), LL(A)+L(B)))
check("rung3 root inverse: c=(1+a)^(1/(1+b))-1 solves c(3)b=a", np.allclose(r3((1+A)**(1/(1+B))-1,B), A))
check("rung3 log inverse: b=L(a)/L(c)-1 solves c(3)b=a", np.allclose(r3(C, L(A)/L(C)-1), A))
check("rung3 not associative (sample 1,2,1)", r3(r3(1,2),1)!=r3(1,r3(2,1)))
print("   (1(3)2)(3)1 =", r3(r3(1,2),1), "  1(3)(2(3)1) =", r3(1,r3(2,1)), "  (note (1(3)1)(3)1 = 1(3)(1(3)1) = 15 by coincidence)")
check("rung3 not commutative (1,2)", r3(1,2)!=r3(2,1))
lid = [ (1+b)**(1/(1+b))-1 for b in (0.5,1,2)]
check("no left identity: needed e depends on b", len(set(np.round(lid,9)))==3)
print("   left 'identity' needed for b=0.5,1,2:", np.round(lid,4))

# Fixed-step flows of each rung (iterate the rung-n operation with a fixed step) and their ODE speeds
# rung2 step Y -> Y*c (Y=1+a): Y(t)=Y0 c^t, Y'=log(c) Y
# rung3 step Y -> Y^c: Y(t)=Y0^(c^t), Y'=log(c) Y log Y
lam=math.log(1.7); Y0=3.0
t=np.arange(0,6)
Yd = [Y0]
for _ in range(5): Yd.append(Yd[-1]**1.7)
Yc = np.exp(np.log(Y0)*np.exp(lam*t))
check("iterating Y->Y^c equals the flow of Y'=log(c)·Y·log Y", np.allclose(Yd, Yc))
# chain rule ladder F_{n+1}(Y) = Y F_n(log Y)
F = [lambda y: 1+0*y, lambda y: y]
for n in range(3): F.append((lambda Fn: (lambda y: y*Fn(np.log(y))))(F[-1]))
y = np.array([30.,1e3,1e6])
check("ladder F3 = Y log Y", np.allclose(F[2](y), y*np.log(y)))
check("ladder F4 = Y log Y loglog Y", np.allclose(F[3](y), y*np.log(y)*np.log(np.log(y))))
# Osgood: integral of dY/F up to big X grows without bound for ladder speeds, converges for Y^p
def osg(f, a, X):  # integrate in s=log Y to avoid huge numbers
    s=np.linspace(math.log(a), math.log(X), 200001); Y=np.exp(s)
    return np.trapezoid(1/(f(Y)/Y), s)
for X in (1e10,1e100,1e300):
    print(f"   Osgood integrals to {X:.0e}: Y logY: {osg(F[2],20,X):8.3f}   Y logY loglogY: {osg(F[3],20,X):7.3f}   Y^1.5: {osg(lambda y:np.exp(1.5*np.log(y)),20,X):.4f}   Y^2: {osg(lambda y:np.exp(2*np.log(y)),20,X):.4f}")
check("Osgood: Y logY integral grows like loglog X (diverges)", abs(osg(F[2],20,1e300) - (math.log(math.log(1e300))-math.log(math.log(20))))<1e-3)
check("Osgood: Y^2 integral converges to 1/20", abs(osg(lambda y:np.exp(2*np.log(y)),20,1e300)-0.05)<1e-6)
print(f"\n{sum(ok)}/{len(ok)} checks passed")
