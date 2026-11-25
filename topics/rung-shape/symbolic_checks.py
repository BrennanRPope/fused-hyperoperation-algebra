"""Exact checks (sympy) of every statement in lean/RungShape.lean."""
import sympy as sp
a, b, c, h, k, t, s, x, y, z, d = sp.symbols('a b c h k t s x y z d')
M = lambda p, q: p + q + p*q
W = lambda w, p, q: p + q + w*p*q
ok = []
def chk(name, cond):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
chk("twist_rung1", sp.expand(M(a+h, b+k) - M(a+h, b) - M(a, b+k) + M(a, b) - h*k) == 0)
chk("twist_dial", sp.expand(W(t, a+h, b+k) - W(t, a+h, b) - W(t, a, b+k) + W(t, a, b) - t*h*k) == 0)
chk("ruled_rung1", sp.expand(M(a, b) - (a + (1+a)*b)) == 0)
chk("ruled_rung1'", sp.expand(M(a, b) - (b + (1+b)*a)) == 0)
chk("axis_rung1", sp.expand(M(a, 0) - a) == 0)
# flow_heights_add: target is in the ideal of the two hypotheses
H1 = y*(2 - t*x) - 2*x; H2 = z*(2 - s*y) - 2*y
T = 2*(z*(2 - (s+t)*x)) - 2*(2*x)
chk("flow_heights_add (ideal membership, explicit combination)", sp.expand(T - ((2 - t*x)*H2 + 2*H1 + s*z*H1)) == 0)
G = sp.groebner([H1, H2], z, y, x, s, t, order='lex')
chk("flow_heights_add (Groebner reduction)", G.reduce(sp.expand(T))[1] == 0)
chk("flow_zero", sp.expand(x*(2 - 0*x) - 2*x) == 0)
chk("sheet_rung0", sp.simplify(a + (c - a) - c) == 0)
chk("sheet_rung1", sp.simplify(M(a, (c - a)/(1 + a)) - c) == 0)
chk("sheet_dial", sp.simplify(W(t, a, (c - a)/(1 + t*a)) - c) == 0)
chk("sheet_edges", M(c, 0) == c and M(0, c) == c)
chk("flow_closed", sp.simplify(2*x/(2 - t*x)*(2 - t*x) - 2*x) == 0)
# near-identity flow is the x^2/2 flow: d/dt phi_t(x) = phi^2/2, phi_0 = x
phi = 2*x/(2 - t*x)
chk("phi_t solves dphi/dt = phi^2/2", sp.simplify(sp.diff(phi, t) - phi**2/2) == 0)
chk("phi heights add: phi_s(phi_t(x)) = phi_{s+t}(x)", sp.simplify((2*phi/(2 - s*phi)) - 2*x/(2 - (s+t)*x)) == 0)
# and it matches E^t to second order: E^t(x) = x + t x^2/2 + ...
chk("phi matches the regular iterate to 2nd order (x + t x^2/2)", sp.expand(sp.series(phi, x, 0, 3).removeO() - (x + t*x**2/2)) == 0)
# counts
chk("addPow: d copies of a added = d*a (d = 0..12)", all(sum([a]*n) == n*a for n in range(13)))
mp = lambda p, n: 0 if n == 0 else M(mp(p, n-1), p)
chk("cube_count: 1 + mpow 1 d = 2^d (d = 0..20)", all(1 + mp(1, n) == 2**n for n in range(21)))
# saddle curvature at identity for the near-identity form a + b + s ab (+ higher terms)
Bs = a + b + s*a*b
fa, fb = sp.diff(Bs, a), sp.diff(Bs, b)
K = (sp.diff(Bs, a, 2)*sp.diff(Bs, b, 2) - sp.diff(Bs, a, b)**2)/(1 + fa**2 + fb**2)**2
chk("Gaussian curvature of the graph at the identity = -s^2/9", sp.simplify(K.subs({a: 0, b: 0}) + s**2/9) == 0)
chk("sheet_rung1_lens: a+b+ab = c  <=>  (1+a)(1+b) = 1+c", sp.expand((1+a)*(1+b) - (1+c) - (M(a, b) - c)) == 0)
print(f"\n{sum(ok)} of {len(ok)} exact checks pass")
