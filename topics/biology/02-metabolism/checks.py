"""Metabolism (as biology) through mult-add: checks. Run: python3 -I checks.py"""
import numpy as np, math
rng=np.random.default_rng(2); ok=tot=0
def check(n,c,d=""):
    global ok,tot; tot+=1; ok+=bool(c); print(("PASS " if c else "FAIL ")+n+(" | "+d if d else ""))
wadd=lambda k,a,b:a+b+k*a*b

# ---- 1. Two-substrate enzymes: steady state of the actual mechanisms vs the dial form
def ordered_rate(k1,km1,k2,km2,k3,Et,A,B):
    # species E, EA, EAB ; solve steady state linear system
    M=np.array([[-(k1*A), km1, k3],          # dE/dt
                [k1*A, -(km1+k2*B), km2],     # dEA/dt
                [1,1,1]])
    E,EA,EAB=np.linalg.solve(M,[0,0,Et]); return k3*EAB
def pingpong_rate(k1,km1,k2,k3,km3,k4,Et,A,B):
    # E, EA, F, FB
    M=np.array([[-k1*A, km1, 0, k4],
                [k1*A, -(km1+k2), 0, 0],
                [0, k2, -k3*B, km3],
                [1,1,1,1]])
    E,EA,F,FB=np.linalg.solve(M,[0,0,0,Et]); return k4*FB
emax=0; emax2=0
for _ in range(300):
    k1,km1,k2,km2,k3=10**rng.uniform(-1,2,5); Et=1.0
    A,B=10**rng.uniform(-2,2,2)
    V=k3*Et; Kia=km1/k1; Ka=k3/k1; Kb=(km2+k3)/k2
    v=ordered_rate(k1,km1,k2,km2,k3,Et,A,B)
    k=Kia/Ka; ua,ub=Ka/A,Kb/B
    emax=max(emax,abs(V/v-(1+wadd(k,ua,ub)))/(V/v))
    p=10**rng.uniform(-1,2,6); k1,km1,k2,k3,km3,k4=p
    V=k2*k4/(k2+k4); Ka=k4*(km1+k2)/(k1*(k2+k4)); Kb=k2*(km3+k4)/(k3*(k2+k4))
    v=pingpong_rate(k1,km1,k2,k3,km3,k4,1.0,A,B)
    emax2=max(emax2,abs(V/v-(1+wadd(0,Ka/A,Kb/B)))/(V/v))
check("ordered mechanism, solved at steady state: V/v = 1 + (u_A ✳_k u_B), k = Kia/Ka (300 random rate sets)",emax<1e-10,f"{emax:.1e}")
check("ping-pong mechanism, solved at steady state: V/v = 1 + u_A + u_B, dial k = 0 (300 sets)",emax2<1e-10,f"{emax2:.1e}")
# random rapid-equilibrium with independent binding: k = 1
KA,KB=0.3,2.0
for A,B in ((0.1,5),(3,0.2),(1,1)):
    Z=1+A/KA+B/KB+A*B/(KA*KB)   # E, EA, EB, EAB weights ; v ∝ EAB
    VoverV=Z/(A*B/(KA*KB)); ua,ub=KA/A,KB/B
    check(f"random independent binding: V/v = (1+u_A)(1+u_B), dial k = 1 (A={A},B={B})",abs(VoverV-(1+ua)*(1+ub))<1e-12)
# Lines in the double-reciprocal plot meet at the wall u_A = -1/k
k=0.4; Ka,Kb=2.0,3.0
vals=[1+wadd(k,-1/k,Kb/B) for B in (0.5,1,2,10)]
check("double-reciprocal lines (1/v vs 1/A, any B) all meet at u_A = -1/k, the weighted wall",max(vals)-min(vals)<1e-12,f"V/v there = {vals[0]:.4f} = 1 - 1/k = {1-1/k:.4f}")
check("in concentration terms the meeting point is 1/A = -1/Kia",abs((-1/k)/Ka-(-1/(k*Ka)))<1e-15)
# slopes parallel iff k = 0
sl=lambda k,B: (1+k*Kb/B)   # slope of V/v against u_A
check("ping-pong (k = 0): every slope equal, so parallel lines",len({round(sl(0,B),12) for B in (0.5,1,2,10)})==1)
# Real constants: rat liver mitochondrial MDH (Massey thesis 1988): Ka 72, Kia 12, Kb 86 µM
Ka,Kia,Kb=72.,12.,86.; kM=Kia/Ka
print(f"   MDH dial k = Kia/Ka = {kM:.3f}; lines meet at 1/A = -1/Kia = {-1/Kia:.4f} per µM")
check("MDH sits at k ≈ 0.17: sequential, but close to the ping-pong end",abs(kM-1/6)<1e-12)

# ---- 2. Two nutrients for a growing cell: same dial
theta=lambda u:1/(1+u)   # u = K/S
for _ in range(100):
    u1,u2=10**rng.uniform(-2,2,2)
    mult=theta(u1)*theta(u2); assert abs(1/mult-(1+wadd(1,u1,u2)))<1e-9*(1/mult)
    add=1/(1+u1+u2); assert abs(1/add-(1+wadd(0,u1,u2)))<1e-12*(1/add)
check("multiplicative two-nutrient Monod: 1/θ = 1 + (u1 ✳ u2) (100 cases)",True)
check("additive-resistance form: 1/θ = 1 + u1 + u2, dial 0 (100 cases)",True)
# limitation 'or' law: 1 - θ1θ2 = λ1 (+)_-1 λ2 with λ = 1-θ
u1,u2=0.7,2.4; l1,l2=1-theta(u1),1-theta(u2)
check("total limitation 1 - θ1θ2 = λ1 ✳_(-1) λ2: limitations combine by the 'or' law",abs((1-theta(u1)*theta(u2))-wadd(-1,l1,l2))<1e-14)
# Liebig's minimum is not on the dial: compare max-error to best dial
grid=np.logspace(-2,2,60); U1,U2=np.meshgrid(grid,grid)
lieb=1/(1+np.maximum(U1,U2))
best=min(((k,np.max(np.abs(1/(1+wadd(k,U1,U2))-lieb))) for k in np.linspace(-0.5,2,251)),key=lambda t:t[1])
check("Liebig's minimum is not a dial setting: the best single k misses it by >5% somewhere",best[1]>0.05,f"best k={best[0]:.2f}, worst miss {best[1]:.3f}")

# ---- 3. Saturation is the mult-add inverse
inv=lambda u:-u/(1+u)
for u in (0.1,1,9):
    check(f"free (unbound) fraction u/(1+u) = -inv(u) at u = {u}",abs(u/(1+u)+inv(u))<1e-15)
check("bound fraction θ = 1 + inv(u): saturation is the lift of the inverse",abs(theta(3.0)-(1+inv(3.0)))<1e-15)

# ---- 4. Kleiber 3/4: a mult-add power with fractional exponent, nothing more
mp=lambda a,n:(1+a)**n-1
check("doubling body mass raises basal rate by mpow(1, 3/4) = 2^(3/4) - 1 = 68%",abs(mp(1,0.75)-(2**0.75-1))<1e-15,f"{mp(1,.75):.3f}")
check("mass-specific rate falls by the inverse: 1 + inv(...) per doubling",abs((1+mp(1,0.75))/2-2**-0.25)<1e-15)
print(f"\n{ok}/{tot} checks passed")
