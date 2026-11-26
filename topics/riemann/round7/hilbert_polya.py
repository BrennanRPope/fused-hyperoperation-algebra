"""Round 7, item 2: Hilbert-Polya operators and the Euler-product dial.

Bender-Brody-Muller (PRL 2017): H = Delta^-1 (xp+px) Delta, Delta = 1 - e^{-ip}.
psi_s(x) = -zeta(s, x+1) satisfies Delta psi_s = x^-s, (xp+px) x^-s = i(2s-1) x^-s,
so H psi_s = i(2s-1) psi_s, and the boundary condition psi_s(0) = -zeta(s) = 0 picks s = rho.
"Eigenvalue" E = i(2 rho - 1) is real iff Re rho = 1/2.

The recipe never uses the Euler product.  For ANY Dirichlet series f = sum a_n n^-s:
psi_s(x) = -sum a_n (n+x)^-s,  Delta_f has Laplace symbol -1/A(t), A(t) = sum a_n e^{-nt}.
For periodic a_n (period m) this is a rational function of the unit shift:
   psi(x) - psi(x+m) = -sum_{r=1}^{m} a_r (x+r)^-s .
Boundary condition psi_s(0) = -f(s) = 0.  Off-line zeros of f give non-real "eigenvalues".
"""
import json
import mpmath as mp
mp.mp.dps = 30
OUT = {}

# --- 1. BBM for zeta: the intertwining identity at a test point ---
s = mp.mpc(0.5, 14.134725141734693790457251983562)
x = mp.mpf('0.7')
psi = lambda s, x: -mp.zeta(s, x+1)
lhs = psi(s, x) - psi(s, x-1)
print("BBM: Delta psi_s(x) - x^-s =", mp.nstr(abs(lhs - x**(-s)), 3))
print("     psi_s(0) at rho_1 =", mp.nstr(abs(psi(s, 0)), 3), "; E = i(2 rho-1) =", mp.nstr(1j*(2*s-1), 12))
OUT['bbm_zeta'] = {'intertwine_err': float(abs(lhs - x**(-s))), 'psi0': float(abs(psi(s, 0)))}

# --- 2. Davenport-Heilbronn ---
SQ5 = mp.sqrt(5)
XI = (mp.sqrt(10 - 2*SQ5) - 2)/(SQ5 - 1)
A = [1, XI, -XI, -1, 0]
def f_DH(s):
    return 5**(-s)*sum(A[r-1]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 6) if A[r-1] != 0)
def psi_DH(s, x):
    return -5**(-s)*sum(A[r-1]*mp.zeta(s, (x+r)/mp.mpf(5)) for r in range(1, 6) if A[r-1] != 0)
# functional equation check
def Lam_DH(s): return (5/mp.pi)**(s/2)*mp.gamma((s+1)/2)*f_DH(s)
s0 = mp.mpc(0.3, 7.1)
print("DH functional equation |Lam(s)-Lam(1-s)|/|Lam(s)| =", mp.nstr(abs(Lam_DH(s0)-Lam_DH(1-s0))/abs(Lam_DH(s0)), 3))
OUT['dh_fe'] = float(abs(Lam_DH(s0)-Lam_DH(1-s0))/abs(Lam_DH(s0)))
starts = [mp.mpc(0.808517, 85.699348), mp.mpc(0.650830, 114.163343), mp.mpc(0.574356, 166.479306), mp.mpc(0.724258, 176.702461)]
dh = []
for z0 in starts:
    r = mp.findroot(f_DH, z0)
    xx = mp.mpf('0.3')
    inter = (psi_DH(r, xx) - psi_DH(r, xx+5)) + sum(A[k-1]*(xx+k)**(-r) for k in range(1, 6))
    E = 1j*(2*r-1)
    print(f"DH zero {mp.nstr(r, 12)}: |f| {mp.nstr(abs(f_DH(r)), 3)}, psi(0)={mp.nstr(abs(psi_DH(r, 0)), 3)}, "
          f"intertwining err {mp.nstr(abs(inter), 3)}, E = {mp.nstr(E, 10)}")
    dh.append({'zero': [float(mp.re(r)), float(mp.im(r))], 'E': [float(mp.re(E)), float(mp.im(E))],
               'psi0': float(abs(psi_DH(r, 0))), 'inter': float(abs(inter))})
OUT['dh'] = dh

# --- 3. x^2+5y^2 (round 6 zero) via the genus formula ---
def chi_list(D):
    from sympy import jacobi_symbol, factorint
    q = abs(D)
    def kr(n):
        v = 1
        for p, e in factorint(n).items():
            if p == 2:
                if D % 2 == 0: return 0
                w = 1 if D % 8 in (1, 7) else -1
            else:
                w = int(jacobi_symbol(D % p, p)) if D % p else 0
            v *= w**e
        return v
    return [kr(n) for n in range(1, q+1)]
C4, C5, C20 = chi_list(-4), chi_list(5), chi_list(-20)
def Lser(s, c):
    q = len(c)
    return q**(-s)*sum(c[r-1]*mp.zeta(s, mp.mpf(r)/q) for r in range(1, q+1) if c[r-1] != 0)
def Z5(s): return mp.zeta(s)*Lser(s, C20) + Lser(s, C4)*Lser(s, C5)
direct = mp.nsum(lambda m, n: 0 if (m == 0 and n == 0) else (m*m+5*n*n)**(-2), [-mp.inf, mp.inf], [-mp.inf, mp.inf])
print("x^2+5y^2 genus formula at s=2:", mp.nstr(Z5(2), 15), " direct:", mp.nstr(direct, 15))
r5 = mp.findroot(Z5, mp.mpc(0.93297, 15.66825))
E5 = 1j*(2*r5-1)
print("x^2+5y^2 zero", mp.nstr(r5, 12), " E =", mp.nstr(E5, 10))
OUT['x2+5y2'] = {'genus_s2': float(mp.re(Z5(2))), 'direct_s2': float(direct), 'zero': [float(mp.re(r5)), float(mp.im(r5))],
                 'E': [float(mp.re(E5)), float(mp.im(E5))]}

# --- 4. semiclassics: the Berry-Keating cell sees only the conductor ---
# N(T) for degree 1, conductor q: (T/2pi) log(qT/(2 pi e)) + 7/8 + S(T); the xp cell with l_x l_p = 2pi/q gives the first term.
for q, name in [(1, 'zeta'), (5, 'DH')]:
    T = mp.mpf(200)
    print(f"Berry-Keating smooth count, conductor {q} ({name}) at T=200:", mp.nstr(T/(2*mp.pi)*mp.log(q*T/(2*mp.pi*mp.e)) + mp.mpf(7)/8, 8))
json.dump(OUT, open('hilbert_polya.json', 'w'), indent=1)
