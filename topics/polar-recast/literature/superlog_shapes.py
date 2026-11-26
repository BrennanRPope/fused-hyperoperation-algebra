"""Super logarithm vs the shapes of the polar system (DRAFT, Nov 26 2026).
1. Each rung operation B_s(x,y) = E^s(E^-s x + E^-s y) is a one-dimensional formal group law;
   its logarithm is E^-s, i.e. the hypernatural log ln_s(X) = E^-s(X-1). (Lazard framework: known.)
2. Every plane-dial shape (x+y+kxy)/(1-dxy) has logarithm g with g'(x) = 1/(1+kx+dx^2): elementary.
3. Fit the plane dial to the rung at height s to third order and read off the discriminant k^2-4d.
4. The fourth-order gap measures how far the rung is from every conic."""
import sympy as sp
x, s, k, d = sp.symbols("x s k d")
N = 6
lam = x**2/2 - x**3/12 + x**4/48 - x**5/180 + sp.Rational(11, 8640)*x**6   # iterative log of e^x-1
# flow f_s = exp(s * lam * d/dx) x, truncated
term, f = x, x
for n in range(1, N):
    term = sp.expand(lam*sp.diff(term, x)) / n
    term = sp.series(term, x, 0, N).removeO()
    f += s**n * term
f = sp.expand(sp.series(f, x, 0, N).removeO())
# check f_1 = e^x - 1 and group law f_s(f_t) = f_{s+t}
t = sp.symbols("t")
chk1 = sp.simplify(sp.series(f.subs(s, 1) - (sp.exp(x)-1), x, 0, N).removeO())
fst = sp.series(f.subs(x, f.subs(s, t)), x, 0, N).removeO()
chk2 = sp.simplify(sp.expand(fst - f.subs(s, s+t)))
print("flow check f_1 = e^x-1:", chk1 == 0, "| f_s o f_t = f_(s+t):", chk2 == 0)
ruler = sp.expand(f.subs(s, -s))            # E^-s: the hypernatural log's series around X=1
coef = [sp.factor(ruler.coeff(x, n)) for n in range(1, N)]
print("ln_s series (u = X-1):", " + ".join(f"({c})u^{n}" for n, c in enumerate(coef, 1)))
g = sp.series(sp.integrate(sp.series(1/(1+k*x+d*x**2), x, 0, N).removeO(), x), x, 0, N).removeO()
gc = [sp.expand(g.coeff(x, n)) for n in range(1, N)]
sol = sp.solve([sp.Eq(gc[1], coef[1]), sp.Eq(gc[2], coef[2])], [k, d], dict=True)[0]
print("third-order fit: k =", sol[k], ", d =", sp.factor(sol[d]), ", k^2 - 4d =", sp.factor(sol[k]**2 - 4*sol[d]))
gap4 = sp.factor(coef[3] - gc[3].subs(sol))
print("fourth-order gap (rung minus best conic):", gap4)
for sv in [-1, -sp.Rational(1, 2), 0, sp.Rational(1, 2), 1, 2]:
    print(f"  s={sv}: k^2-4d={sp.factor(sol[k]**2-4*sol[d]).subs(s, sv)}, gap4={gap4.subs(s, sv)}")
