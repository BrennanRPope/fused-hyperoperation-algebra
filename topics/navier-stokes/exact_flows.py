"""Navier-Stokes through the circle-operation lens: symbolic checks.
Run: python3 exact_flows.py
For each classic exact solution: verify NS (with a pressure), split the
nonlinearity into corners pot x pot / pot x vort / vort x vort, and test
whether the vort x vort Lamb vector is a gradient (curl-free)."""
import sympy as sp
x, y, z, t, r = sp.symbols('x y z t r', real=True)
nu, a, G, A, B, C, n, lam = sp.symbols('nu alpha Gamma A B C n lambda', positive=True)
X = (x, y, z)

def grad(f): return sp.Matrix([sp.diff(f, v) for v in X])
def div(u): return sum(sp.diff(u[i], X[i]) for i in range(3))
def curl(u): return sp.Matrix([sp.diff(u[2], y)-sp.diff(u[1], z), sp.diff(u[0], z)-sp.diff(u[2], x), sp.diff(u[1], x)-sp.diff(u[0], y)])
def lap(u): return sp.Matrix([sum(sp.diff(u[i], v, 2) for v in X) for i in range(3)])
def adv(u, w): return sp.Matrix([sum(u[j]*sp.diff(w[i], X[j]) for j in range(3)) for i in range(3)])
def simp(M): return sp.simplify(sp.Matrix(M))

print("== 1. Cole-Hopf: u = -2 nu grad log phi ==")
phi = sp.Function('phi')(x, y, z, t)
u = -2*nu*grad(sp.log(phi))
res = u.diff(t) + adv(u, u) - nu*lap(u)
target = -2*nu*grad((phi.diff(t) - nu*sum(phi.diff(v, 2) for v in X))/phi)
print("Burgers residual minus -2nu grad(heat/phi):", simp(res - target).T)
print("curl u:", simp(curl(u)).T)
# the cancelling cross term
f = sp.Function('f')(x)
print("d2 log f = f''/f - (f'/f)^2 :", sp.simplify(sp.diff(sp.log(f), x, 2) - (f.diff(x, 2)/f - (f.diff(x)/f)**2)) == 0)

def check(name, u, upot, p=None):
    print(f"\n== {name} ==")
    w = simp(curl(u)); uv = simp(u - upot)
    print("div u:", sp.simplify(div(u)))
    print("curl(upot) (should be 0):", simp(curl(upot)).T)
    lamb_vv = simp(uv.cross(w)); lamb_pv = simp(upot.cross(w))
    print("curl of vort x vort Lamb vector:", simp(curl(lamb_vv)).T)
    print("curl of pot x vort Lamb vector :", simp(curl(lamb_pv)).T)
    stretch = simp(adv(w, u))
    print("stretching (w.grad)u:", stretch.T)
    print("helicity density u.w:", sp.simplify(u.dot(w)))
    if p is not None:
        R = simp(u.diff(t) + adv(u, u) + grad(p) - nu*lap(u))
        print("NS residual with given pressure:", R.T)
    else:
        # NS holds iff curl of (u_t - nu lap u + (u.grad)u) vanishes
        R = simp(curl(u.diff(t) + adv(u, u) - nu*lap(u)))
        print("curl of NS momentum residual (0 => a pressure exists):", R.T)

rr2 = x**2 + y**2
# Lamb-Oseen
s = rr2/(4*nu*t)
uth = G/(2*sp.pi)*(1 - sp.exp(-s))/rr2
check("Lamb-Oseen vortex", sp.Matrix([-y*uth, x*uth, 0]), sp.zeros(3, 1))
# Burgers vortex (steady)
uthb = G/(2*sp.pi)*(1 - sp.exp(-a*rr2/(4*nu)))/rr2
upot = sp.Matrix([-a*x/2, -a*y/2, a*z])
check("Burgers vortex", upot + sp.Matrix([-y*uthb, x*uthb, 0]), upot)
# Taylor-Green 2D
E = sp.exp(-2*nu*t)
check("Taylor-Green (2D)", sp.Matrix([sp.sin(x)*sp.cos(y)*E, -sp.cos(x)*sp.sin(y)*E, 0]), sp.zeros(3, 1))
# ABC
E1 = sp.exp(-nu*t)
uabc = sp.Matrix([A*sp.sin(z)+C*sp.cos(y), B*sp.sin(x)+A*sp.cos(z), C*sp.sin(y)+B*sp.cos(x)])*E1
check("ABC (Beltrami)", uabc, sp.zeros(3, 1))
print("curl u - u (ABC eigen):", simp(curl(uabc) - uabc).T)
# Kolmogorov (forced; check unforced decaying version)
check("Kolmogorov shear", sp.Matrix([sp.sin(n*y)*sp.exp(-nu*n**2*t), 0, 0]), sp.zeros(3, 1))
# 3D Taylor-Green initial condition (not an exact solution)
utg = sp.Matrix([sp.sin(x)*sp.cos(y)*sp.cos(z), -sp.cos(x)*sp.sin(y)*sp.cos(z), 0])
w = simp(curl(utg))
print("\n== 3D Taylor-Green initial data (the classic blow-up test) ==")
print("curl of vort x vort Lamb:", simp(curl(utg.cross(w))).T)
print("helicity density:", sp.simplify(utg.dot(w)))
print("stretching (w.grad)u:", simp(adv(w, utg)).T)
