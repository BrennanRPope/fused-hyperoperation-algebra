"""Exact orbit equation around a four-charge (STU) black hole for any visitor made of BPS / anti-BPS pieces:
L^2 (du/dphi)^2 = D^2 H + 2 D sum_j w_j s_j H/H_j + sum_{j != l} w_j w_l (s_j s_l - 1) H/(H_j H_l) - L^2 u^2,
D = E - sum_j w_j s_j,  H = prod (1 + k_i u).  Every term is a polynomial in u: its coefficients are the elementary
symmetric polynomials of the charges (e1 = sum, e2 = pairs, e3 = triples, e4 = product)."""
import sys; sys.path.insert(0, '..')
import sympy as sp
from misalign import W_stu
from precess import u, L, E
k = sp.symbols('k1:5', positive=True); w = sp.symbols('w1:5', positive=True)
def poly_F(ks, ws, ss):
    Hs = [1 + kk*u for kk in ks]; H = sp.prod(Hs); D = E - sum(a*b for a, b in zip(ws, ss))
    out = D**2*H + 2*D*sum(ws[j]*ss[j]*H/Hs[j] for j in range(4))
    out += sum(ws[j]*ws[l]*(ss[j]*ss[l] - 1)*H/(Hs[j]*Hs[l]) for j in range(4) for l in range(4) if j != l)
    return sp.expand(sp.cancel(out)) - L**2*u**2
if __name__ == "__main__":
    import itertools, random
    random.seed(1); ok = True
    for ss in itertools.product([1, -1], repeat=4):
        W = W_stu(list(k), list(w), list(ss)); P = poly_F(list(k), list(w), list(ss))
        vals = {**{kk: sp.Rational(random.randint(1, 9), random.randint(1, 5)) for kk in k}, **{ww: sp.Rational(random.randint(1, 9), 7) for ww in w},
                u: sp.Rational(1, 13), E: sp.Rational(9, 10), L: sp.Integer(3)}
        diff = sp.N((W - P).subs(vals), 30); ok &= abs(diff) < 1e-20
        print(ss, "degree in u:", sp.degree(P, u), " W - poly =", sp.N(diff, 5))
    print("exact polynomial form:", ok)
