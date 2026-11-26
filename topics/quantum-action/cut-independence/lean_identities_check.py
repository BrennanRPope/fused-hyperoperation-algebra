"""Python mirror of every statement in lean/CutIndependence.lean (the file is also compiled in the shared library since Nov 26)."""
import sympy as sp
ok = []
b, c, d, e, k, p, q = sp.symbols("b c d e k p q")
ok.append(("swap_assoc", sp.expand(b*c + (b+d)*e - (d*e + b*(c+e))) == 0))
ok.append(("leftover_zero_iff", sp.solve(sp.Eq(p - q, 0), p) == [q]))
def u(kv, n):
    a, bb = 0, 1
    if n == 0: return 0
    for _ in range(n - 1): a, bb = bb, sp.expand(kv*bb - a)
    return bb
ok.append(("u_two", sp.expand(u(k, 2) - k) == 0))
ok.append(("cassini n<=12", all(sp.expand(u(k, n+1)**2 - u(k, n)*u(k, n+2) - 1) == 0 for n in range(13))))
ok.append(("dims_r3", (u(1, 1), u(1, 2), u(1, 3)) == (1, 1, 0)))
ok.append(("dims_r2", (u(0, 1), u(0, 2)) == (1, 0)))
ok.append(("w2_r3, w2_r2", u(1, 1)**2 + u(1, 2)**2 == 2 and u(0, 1)**2 == 1))
M = lambda kv: sp.Matrix([[kv, -1], [1, 0]])
ok.append(("step_closes_r3", M(1)**3 == -sp.eye(2)))
ok.append(("step_closes_r2", M(0)**2 == -sp.eye(2)))
ok.append(("u_at_two n<=40", all(u(2, n) == n for n in range(41))))
ok.append(("disc_factor", sp.expand(4 - k*k - (2 - k)*(2 + k)) == 0))
for name, v in ok: print(("PASS " if v else "FAIL ") + name)
print(f"{sum(v for _, v in ok)}/{len(ok)} pass")
