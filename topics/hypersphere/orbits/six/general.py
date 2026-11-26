import sympy as sp
r, E, L = sp.symbols('r E L', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r); V = sp.Function('V')(r)
# Static metric -A dt^2 + B dr^2 + r^2 dphi^2. Proper-time radial eq: B rdot^2 = G, G = E^2/A - 1 - L^2/r^2
G = E**2/A - 1 - L**2/r**2
G1 = sp.diff(G, r); G2 = sp.diff(G, r, 2)
# circular: G = 0, G1 = 0 -> solve for E^2, L^2
sol = sp.solve([sp.Eq(G, 0), sp.Eq(G1, 0)], [E**2, L**2], dict=True)[0]
E2s, L2s = sol[E**2], sol[L**2]
# radial proper freq^2 = -G2/(2B); phi proper freq = L/r^2
T = sp.simplify((-G2.subs({E**2: E2s, L**2: L2s}) / (2 * B)) / (L2s / r**4))
v2 = r * sp.diff(A, r) / (2 * A)      # local speed^2 measured by static observer
guess = (3 + r * sp.diff(A, r, 2) / sp.diff(A, r) - 4 * v2) / B
print('GR general: T - (3 + rA\'\'/A\' - 4v^2)/B =', sp.simplify(T - guess))
# SR particle in potential V: H = sqrt(1+p^2) + V(r); m=c=1
pr, Lm = sp.symbols('p_r L_m', positive=True)
H = sp.sqrt(1 + pr**2 + Lm**2 / r**2) + V
Hr = sp.diff(H, r).subs(pr, 0); Hrr = sp.diff(H, r, 2).subs(pr, 0); Hpp = sp.diff(H, pr, 2).subs(pr, 0); HL = sp.diff(H, Lm).subs(pr, 0)
Lsol = sp.solve(sp.Eq(Hr, 0), Lm)
for Ls in Lsol:
    Tsr = sp.simplify((Hrr * Hpp / HL**2).subs(Lm, Ls))
    vsr2 = sp.simplify((r * HL).subs(Lm, Ls) ** 2)
    g2 = (3 + r * sp.diff(V, r, 2) / sp.diff(V, r) - vsr2)
    print('SR general: T - (3 + rV\'\'/V\' - v^2) =', sp.simplify(Tsr - g2))
