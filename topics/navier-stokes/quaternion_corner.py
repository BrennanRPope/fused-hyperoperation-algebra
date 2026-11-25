"""The quaternion product of velocity and vorticity, u*w = -(u.w) + (u x w):
its number part is minus the helicity density, its arrow part is the Lamb
vector (the piece of the nonlinearity that the log move cannot absorb)."""
import sympy as sp
from sympy.algebras.quaternion import Quaternion as Q
x, y, z = sp.symbols('x y z', real=True); A, B, C = sp.symbols('A B C', positive=True)
def curl(u): return [sp.diff(u[2], y)-sp.diff(u[1], z), sp.diff(u[0], z)-sp.diff(u[2], x), sp.diff(u[1], x)-sp.diff(u[0], y)]
flows = {
 'ABC (Beltrami)': [A*sp.sin(z)+C*sp.cos(y), B*sp.sin(x)+A*sp.cos(z), C*sp.sin(y)+B*sp.cos(x)],
 'Taylor-Green 3D start': [sp.sin(x)*sp.cos(y)*sp.cos(z), -sp.cos(x)*sp.sin(y)*sp.cos(z), 0],
}
for name, u in flows.items():
    w = curl(u); q = Q(0, *u)*Q(0, *w)
    print(name, '| number part:', sp.simplify(q.a), '| arrow part:', [sp.simplify(c) for c in (q.b, q.c, q.d)])
# The circle operation in variables that square to zero and anticommute: the cube's layers
e = sp.symbols('e1:4', commutative=False)
print("Grassmann cube (1+e1)(1+e2)(1+e3) has 1,3,3,1 terms by degree: layers of Lambda(R^3) = pressure, velocity, vorticity, helicity")
