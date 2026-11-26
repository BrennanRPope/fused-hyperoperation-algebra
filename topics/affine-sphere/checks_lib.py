import sympy as sp
u, v = sp.symbols('u v', real=True)
def chord_surface(k, d):
    q = 1 + k * u + d * (u ** 2 + v ** 2)
    return [u / q, (u ** 2 + v ** 2) / q, v / q]
