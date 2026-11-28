"""DRAFT twist checks. Units c = 1, linear order. Question (Brennan: "Test it"): given only the winding-rate pattern
(the four-leaf, E) and the wave's direction, does the chain itself force the tilt's twist B = i*E, without Einstein's
field equations?  Inputs allowed: E, the direction n = z, the chain's speed limit (c, from ../tilt S1-S3), and that
curvature comes from composing steps along paths (associative chain -> Bianchi identity)."""
import json, numpy as np, sympy as sp
R = {}
def ok(name, cond, info=""):
    R[name] = {"pass": bool(cond), "info": info}; print(("PASS " if cond else "FAIL ") + name, info)
t, x, y, z = sp.symbols('t x y z', real=True); X = [t, x, y, z]
def riem(h, a, b, c, d):
    D = lambda i, j, M: sp.diff(M, X[i], X[j])
    return sp.Rational(1, 2)*(D(b, c, h[a, d]) + D(a, d, h[b, c]) - D(a, c, h[b, d]) - D(b, d, h[a, c]))
eps = sp.LeviCivita
def EB(h):
    E = sp.Matrix(3, 3, lambda i, j: riem(h, 0, i+1, 0, j+1))
    B = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2)*sum(eps(i+1, p, q)*riem(h, p, q, 0, j+1) for p in (1, 2, 3) for q in (1, 2, 3)))
    return E, B
curl = lambda M: sp.Matrix(3, 3, lambda i, j: sum(eps(i+1, p+1, q+1)*sp.diff(M[q, j], X[p+1]) for p in range(3) for q in range(3)))

# B1: for ANY metric (all ten pieces arbitrary functions), dB/dt = curl E. This is the Bianchi identity: it only says the
#     curvature comes from composing steps along paths. No Einstein equation is used.
H = sp.Matrix(4, 4, lambda i, j: sp.Function(f'h{min(i,j)}{max(i,j)}')(t, x, y, z))
E, B = EB(H)
res_plus = sp.simplify(sp.diff(B, t) - curl(E)); res_minus = sp.simplify(sp.diff(B, t) + curl(E))
sign = 1 if res_plus == sp.zeros(3) else (-1 if res_minus == sp.zeros(3) else 0)
ok("B1 for every metric: dB/dt = curl E identically (Bianchi; no field equation)", sign == 1, f"sign found: {sign}")

# B2: the partner relation dE/dt = -curl B is NOT an identity: a rate-only field breaks it. That one is where Einstein's
#     vacuum equations come in, and we don't use it.
h2 = sp.zeros(4); h2[0, 0] = sp.Function('F')(t, x, y, z)
E2, B2 = EB(h2); partner = sp.simplify(sp.diff(E2, t) + curl(B2))
ok("B2 partner relation fails for a general rate field (needs the field equations)", partner != sp.zeros(3))

# B3: give only the rate pattern: a four-leaf E moving along z at speed s, nothing before it arrives. B1 then fixes B.
s_ = sp.symbols('s', positive=True); u = t - z/s_
e1, e2 = sp.Function('e1'), sp.Function('e2')            # plus and cross amplitudes of the four-leaf
Ein = sp.Matrix([[e1(u), e2(u), 0], [e2(u), -e1(u), 0], [0, 0, 0]])
rhs = curl(Ein)                                            # dB/dt must equal this
# B depends on u only (it arrives with the wave), so dB/dt = dB/du; integrate in u with B = 0 before arrival
U = sp.symbols('U')
Bder = sp.Matrix(3, 3, lambda i, j: sp.simplify(rhs[i, j]))
# each entry of curl E is (+-1/s) * d e/du; so B = (+-1/s) e
Bfix = sp.Matrix(3, 3, lambda i, j: sp.integrate(sp.simplify(Bder[i, j].subs(t, U + z/s_)), U) if Bder[i, j] != 0 else 0)
Bfix = Bfix.subs(U, u)
ok("B3 the four-leaf plus its direction forces a nonzero twist B (the tilt is not optional)", Bfix != sp.zeros(3),
   f"B_xx = {sp.simplify(Bfix[0,0])}, B_xy = {sp.simplify(Bfix[0,1])}")

# B4: at s = c (the chain's speed limit) the forced twist is exactly i*E: complex B = i * complex E, the 45 degree turn.
Ec = Ein[0, 0] - sp.I*Ein[0, 1]; Bc = (Bfix[0, 0] - sp.I*Bfix[0, 1]).subs(s_, 1)
ratio = sp.simplify(Bc/Ec.subs(s_, 1))
ok("B4 at the speed limit c the forced twist is B = i*E (same as Einstein's wave in ../tilt A8)", ratio == sp.I and sp.simplify(Bfix[0,0].subs(s_,1) - Ein[0,1].subs(s_,1)) == 0 and sp.simplify(Bfix[0,1].subs(s_,1) + Ein[0,0].subs(s_,1)) == 0, f"B/E = {ratio}; B_xx = E_xy, B_xy = -E_xx as in ../tilt A2")
ratio_s = sp.simplify((Bfix[0, 0] - sp.I*Bfix[0, 1])/Ec)

# B5: at any other speed the twist is i*E/s: the right size only at c.
ok("B5 moving at speed s gives B = i*E/s: twist = i*E only at s = c", sp.simplify(ratio_s - sp.I/s_) == 0, f"B/E = {ratio_s}")

# B6: cross-check, not an input: the partner relation (Einstein's vacuum equations) then holds only at s = c.
partner3 = sp.simplify(sp.diff(Ein, t) + curl(Bfix))
res_s = sp.simplify(partner3[0, 0])
ok("B6 cross-check: the forced twist also satisfies Einstein's partner relation, and only at s = c",
   sp.simplify(partner3.subs(s_, 1)) == sp.zeros(3) and sp.simplify(res_s.subs(s_, 2)) != 0, f"partner xx residual = {res_s}")

# B7: the four-leaf moving along z is not the Hessian of any rate field (curl E != 0), so it cannot live in the rate alone;
#     a rate pattern's curvature always has curl 0 (and so B stays 0, ../tilt A1).
ok("B7 a moving four-leaf has curl E != 0, while any pure rate field has curl E = 0", curl(Ein) != sp.zeros(3)
   and sp.simplify(curl(E2)) == sp.zeros(3))

sv = np.linspace(0.3, 3, 60)
json.dump({"checks": R, "speeds": list(sv), "twist_over_rate": list(1/sv),
           "partner_residual_factor": str(res_s)}, open(__file__.replace("checks.py", "results.json"), "w"), indent=1)
print(f"{sum(r['pass'] for r in R.values())} of {len(R)} pass")
