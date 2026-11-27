"""DRAFT checks for the prior-art audit (Nov 27 2026). Run: python3 checks.py > checks_output.txt
Each check confirms a computed statement the audit's README leans on. numpy + sympy only."""
import numpy as np, sympy as sp

rng = np.random.default_rng(0)
results = []
def report(name, ok, detail):
    results.append(ok); print(f"{name}: {'PASS' if ok else 'FAIL'}  {detail}")

# K1. The repo's 3D law (round 12) is the kappa-stereographic addition with kappa = q
#     (Bachmann, Becigneul, Ganea 2020, eq. 2), and its quaternion form (u+v)(1 - q*conj(u)*v)^-1.
def qmul(a, b):
    a0, a1 = a[0], a[1:]; b0, b1 = b[0], b[1:]
    return np.concatenate([[a0*b0 - a1@b1], a0*b1 + b0*a1 + np.cross(a1, b1)])
def qinv(a): return np.array([a[0], *(-a[1:])]) / (a@a)
def conj(a): return np.array([a[0], *(-a[1:])])
def law_quat(u, v, q):
    U, V = np.concatenate([[0], u]), np.concatenate([[0], v])
    one = np.array([1., 0, 0, 0])
    return qmul(U + V, qinv(one - q*qmul(conj(U), V)))
def kappa_add(x, y, k):   # Bachmann et al. 2020, eq. (2)
    xy, x2, y2 = x@y, x@x, y@y
    return ((1 - 2*k*xy - k*y2)*x + (1 + k*x2)*y) / (1 - 2*k*xy + k*k*x2*y2)
def ungar_mobius(x, y):   # Ungar's Mobius addition in the unit ball
    xy, x2, y2 = x@y, x@x, y@y
    return ((1 + 2*xy + y2)*x + (1 - x2)*y) / (1 + 2*xy + x2*y2)
err1 = err2 = err3 = 0
for _ in range(500):
    q = rng.uniform(-2, 2); u, v = rng.normal(size=3)*0.4, rng.normal(size=3)*0.4
    w = law_quat(u, v, q)
    err1 = max(err1, abs(w[0]))                                   # stays a pure quaternion (3D)
    err2 = max(err2, np.max(np.abs(w[1:] - kappa_add(u, v, q))))  # equals kappa-addition with kappa = q
    err3 = max(err3, np.max(np.abs(kappa_add(u, v, -1.0) - ungar_mobius(u, v))))
report("K1 law = kappa-stereographic addition (kappa = q); q = -1 is Ungar's Mobius addition",
       max(err1, err2, err3) < 1e-12, f"real part {err1:.1e}, vs kappa-add {err2:.1e}, vs Ungar {err3:.1e} (500 random cases)")

# K2. Round 25's tide 3/2(qbar - q_loc), with q = q_vac + 8*pi*rho/3, is the textbook tidal split of a
#     spherical mass: radial minus transverse tidal eigenvalue = 4*pi*G*(rhobar - rho) (Newtonian; Sussman's
#     E ~ rho - <rho> in GR). Symbolic.
r, G = sp.symbols('r G', positive=True)
m = sp.Function('m')(r)
gr = -G*m/r**2                                    # radial gravity
radial = sp.diff(gr, r)                           # -Phi''
transverse = gr/r                                 # -Phi'/r
rho = sp.diff(m, r)/(4*sp.pi*r**2); rhobar = 3*m/(4*sp.pi*r**3)
gap = sp.simplify(radial - transverse - 4*sp.pi*G*(rhobar - rho))
qloc, qbar = 8*sp.pi*G*rho/3, 8*sp.pi*G*rhobar/3
gap2 = sp.simplify(sp.Rational(3, 2)*(qbar - qloc) - 4*sp.pi*G*(rhobar - rho))
report("K2 tide 3/2(qbar - q_loc) = 4 pi G (rhobar - rho) = textbook tidal split", gap == 0 and gap2 == 0,
       f"residuals {gap}, {gap2}")

# K3. The Schwarzschild split 1 - 6M/r = (1 - 2M/r)(1 - (2v)^2) holds, but the same split is NOT the
#     stability condition of every static round spacetime: in Schwarzschild-de Sitter the ISCO local speed
#     is not c/2 (the hypersphere page's general rule T = (N - 4v^2)/B is the correct general form; the
#     bridge survey's one-line summary of N2 is too broad).
M, rr, L = sp.symbols('M r Lambda', positive=True)
f = 1 - 2*M/rr
v2 = rr*sp.diff(f, rr)/(2*f)                      # local orbital speed^2 seen by static observers
split_ok = sp.simplify(f*(1 - 4*v2) - (1 - 6*M/rr)) == 0
def isco_speed(Lam):
    fS = 1 - 2/rr - Lam*rr**2/3                   # M = 1
    stab = rr*fS*sp.diff(fS, rr, 2) - 2*rr*sp.diff(fS, rr)**2 + 3*fS*sp.diff(fS, rr)  # >0 stable
    rI = sp.nsolve(stab, rr, 6.0)
    return float(sp.sqrt((rr*sp.diff(fS, rr)/(2*fS)).subs(rr, rI))), float(rI)
v0, r0 = isco_speed(0.0); v1, r1 = isco_speed(1e-4)
report("K3 Schwarzschild split exact; not general (SdS ISCO speed differs from 1/2)",
       split_ok and abs(v0 - 0.5) < 1e-10 and abs(v1 - 0.5) > 1e-4,
       f"Schwarzschild ISCO r={r0:.6f}, v={v0:.10f}; SdS (Lambda=1e-4, M=1) ISCO r={r1:.6f}, v={v1:.6f}")

# K4. Closed de Sitter (a = cosh t) spans exactly pi of conformal time (used by round 26's 'swings once,
#     then freezes' wave, which Abbassi et al. 2013 already solve for every n).
t = sp.symbols('t', real=True)
span = sp.integrate(1/sp.cosh(t), (t, -sp.oo, sp.oo))
report("K4 conformal time span of closed de Sitter = pi", sp.simplify(span - sp.pi) == 0, f"span = {span}")

# K5. The critical-branching residue 1 - (2/3)kappa (parabolic-fingerprints) is the classical iterative
#     residue 1 - b/a^2 of x - a x^2 + b x^3 with a = sigma^2/2, b = E[Z(Z-1)(Z-2)]/6.
s2, m3 = sp.symbols('sigma2 m3', positive=True)
a, b = s2/2, m3/6
kappa = m3/s2**2
report("K5 residue 1 - (2/3)kappa = classical 1 - b/a^2", sp.simplify((1 - b/a**2) - (1 - sp.Rational(2, 3)*kappa)) == 0,
       "symbolic identity")

# K6. The (k, d) plane law (a+b+k ab)/(1-d ab) is Hirzebruch's two-parameter Todd-genus formal group law
#     (x+y+mu1 xy)/(1+mu2 xy) with mu1 = k, mu2 = -d; chi_y genus at mu1 = y-1, mu2 = y.
x, y, k, d = sp.symbols('x y k d')
law = (x + y + k*x*y)/(1 - d*x*y)
todd = lambda m1, m2: (x + y + m1*x*y)/(1 + m2*x*y)
report("K6 plane law = two-parameter Todd FGL (mu1 = k, mu2 = -d)", sp.simplify(law - todd(k, -d)) == 0, "symbolic identity")

print(f"\n{sum(results)} / {len(results)} pass")
