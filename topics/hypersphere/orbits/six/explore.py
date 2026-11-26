import sympy as sp
r,pr,L,M,k = sp.symbols('r p_r L M k', positive=True)
def ratio(H, vlocal):
    # circular orbit: dH/dr = 0 at pr=0; epicyclic omega_r^2 = H_rr * H_prpr ; Omega = H_L
    Hr = sp.diff(H,r).subs(pr,0)
    Hrr = sp.diff(H,r,2).subs(pr,0)
    Hpp = sp.diff(H,pr,2).subs(pr,0)
    HL = sp.diff(H,L).subs(pr,0)
    return Hr, sp.simplify(Hrr*Hpp/HL**2), HL
cases = {}
# Coulomb SR, c=1, m=1
H = sp.sqrt(1+pr**2+L**2/r**2) - k/r
cases['coulomb'] = (H, lambda HL: r*HL, k)
# Nordstrom: H = sqrt(m(r)^2 + p^2), m = 1 - M/r
H = sp.sqrt((1-M/r)**2+pr**2+L**2/r**2)
cases['nordstrom'] = (H, lambda HL: r*HL, M)
# Schwarzschild coordinate-time Hamiltonian
f = 1-2*M/r
H = sp.sqrt(f*(1+f*pr**2+L**2/r**2))
cases['gr'] = (H, lambda HL: r*HL/sp.sqrt(f), M)
for name,(H,vf,par) in cases.items():
    Hr, R, HL = ratio(H, vf)
    # solve circular condition for L^2
    sol = sp.solve(sp.Eq(sp.numer(sp.together(Hr)),0), L)
    sol = [s for s in sol if s.is_positive is not False]
    print(name, 'L solutions', sol)
    for Ls in sol:
        Rs = sp.simplify(R.subs(L,Ls))
        v = sp.simplify(vf(HL).subs(L,Ls))
        print('  ratio^2 =', sp.factor(Rs), '  local v^2 =', sp.factor(sp.simplify(v**2)))

print('--- numeric checks')
Hc = cases['coulomb'][0]
Hr, R, HL = ratio(Hc, None)
for kv, rv in [(1,1),(1,3),(2,1.5)]:
    Ls = sp.sqrt(sp.Rational(1,2)*kv**2 + kv*sp.sqrt(kv**2+4*rv**2)/2)
    Rn = float(R.subs({k:kv, r:rv, L:Ls})); vn = float((r*HL).subs({k:kv,r:rv,L:Ls}))
    print('coulomb ratio^2', Rn, '1-v^2', 1-vn**2)
x, v = sp.symbols('x v', positive=True)
print('GR check', sp.simplify((1-2*x)*(1-4*x/(1-2*x)) - (1-6*x)))
