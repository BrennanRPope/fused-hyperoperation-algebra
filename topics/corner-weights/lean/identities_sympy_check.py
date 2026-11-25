"""Checks with sympy every identity stated in CornerWeights.lean (the Lean file could not be compiled here)."""
import sympy as sp
ok = []
def chk(name, e):
    r = sp.expand(e) == 0; ok.append(r); print(("PASS " if r else "FAIL ") + name)
def det3(m): return sp.Matrix(m).det()
t, a, b = sp.symbols('t alpha beta')
m = sp.symbols('m11 m12 m13 m21 m22 m23 m31 m32 m33'); M = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]]
I3 = sp.eye(3); Mm = sp.Matrix(M)
e1 = m[0] + m[4] + m[8]
e2 = (m[0]*m[4] - m[1]*m[3]) + (m[0]*m[8] - m[2]*m[6]) + (m[4]*m[8] - m[5]*m[7])
e3 = Mm.det()
chk("det3_I_add_tM", (I3 + t*Mm).det() - (1 + t*e1 + t**2*e2 + t**3*e3))
p1, p2, p3, q1, q2, q3 = sp.symbols('p1 p2 p3 q1 q2 q3')
chk("cauchy2_cleared", (p1+q2)*(p2+q1) - (p1+q1)*(p2+q2) - (p1-p2)*(q1-q2))
x = [p1, p2, p3]; y = [q1, q2, q3]
C = [[sp.prod([x[i] + y[l] for l in range(3) if l != j]) for j in range(3)] for i in range(3)]
target = (x[0]-x[1])*(x[0]-x[2])*(x[1]-x[2])*(y[0]-y[1])*(y[0]-y[2])*(y[1]-y[2])
print("   cauchy3 det / target =", sp.factor(det3(C) / target))
chk("cauchy3_cleared", det3(C) - target)
pi, pj, qi, qj = sp.symbols('pi pj qi qj')
chk("kp_pair_cross_ratio", (pi-pj)*(qi-qj) - ((pi+qj)*(pj+qi) - (pi+qi)*(pj+qj)))
chk("kdv_pair_velocity_num", (pi-1)*(pj+1) - (pj-1)*(pi+1) - 2*(pi-pj))
chk("kdv_pair_velocity_den", (pi+1)*(pj+1) - (pi-1)*(pj-1) - 2*(pi+pj))
chk("lucas_from_two_dets_num", (I3 + a*Mm).det() - (I3 + b*Mm).det() - (a-b)*(e1 + (a+b)*e2 + (a*a+a*b+b*b)*e3))
chk("lucas_from_two_dets_den", a*(I3 + b*Mm).det() - b*(I3 + a*Mm).det() - (a-b)*(1 - a*b*(e2 + (a+b)*e3)))
xx, yy = sp.symbols('x y')
N = xx + yy + (a+b)*xx*yy; Dd = 1 - a*b*xx*yy
chk("plane_resolvent_den", Dd + b*N - (1+b*xx)*(1+b*yy))
chk("plane_resolvent_num", N - (xx*(1+b*yy) + yy*(1+b*xx) + (a-b)*xx*yy))
k, d = sp.symbols('k d')
Nc = sp.Matrix([[k, 1], [-d, 0]])
chk("companion_sq", (Nc**2 - (k*Nc - d*sp.eye(2))).applyfunc(sp.expand).norm()**2 if False else sum((Nc**2 - (k*Nc - d*sp.eye(2))).applyfunc(sp.expand)))
U2, U3, U4 = k, k*k - d, k**3 - 2*k*d
chk("companion_cube", sum((Nc**3 - (U3*Nc - d*U2*sp.eye(2))).applyfunc(sp.expand)))
chk("companion_fourth", sum((Nc**4 - (U4*Nc - d*U3*sp.eye(2))).applyfunc(sp.expand)))
A_, c_ = sp.symbols('A c')
two = lambda u, v: u + v + k*A_*u*v
three = a + b + c_ + k*A_*(a*b + a*c_ + b*c_) + k**2*A_**3*a*b*c_
chk("size_pair_not_iterated", three - two(two(a, b), c_) - k**2*A_**2*(A_-1)*a*b*c_)
chk("hard_core_is_addition", three.subs(A_, 0) - (a + b + c_))
# offdiagonal blend: three-body term of diagonal + s*offdiag for a symmetric 3x3 with unit diagonal
s, c12, c13, c23 = sp.symbols('s c12 c13 c23')
B = sp.Matrix([[1, s*c12, s*c13], [s*c12, 1, s*c23], [s*c13, s*c23, 1]])
chk("blend_three_body", B.det() - (1 - s**2*c12**2)*(1 - s**2*c13**2)*(1 - s**2*c23**2)
    - (2*s**3*c12*c13*c23 - s**4*(c12**2*c13**2 + c12**2*c23**2 + c13**2*c23**2) + s**6*c12**2*c13**2*c23**2))
print(f"{sum(ok)} of {len(ok)} identities hold")
