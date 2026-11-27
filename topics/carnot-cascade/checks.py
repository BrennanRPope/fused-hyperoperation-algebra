"""
Carnot cascade and finite-time efficiencies: checks for the S1 reopen (closed-leads-revisit).

Coordinates: tau = T_c/T_h (the lift 1 - eta of Carnot), beta = ln(T_h/T_c), which adds across a cascade.
A law on the plane dial: d eta/d beta = c0 (1 + k eta + d eta^2) (a Riccati equation).
Engines in series always combine by the union rule eta = eta1 + eta2 - eta1 eta2 (k = -1), because heat passes on.
A "knob" of the union rule is tau -> tau^lam; Curzon-Ahlborn is lam = 1/2.

Sections
  A. Which max-power efficiencies are knob images of Carnot (cascade-stable)? Only 1 - tau^lam (Cauchy).
  B. Low-dissipation engines, derived: eta* = eta_C (1 + a)/(2(1 + a) - eta_C).
     Heat-type dissipation (a constant): eta_C/(2 - g eta_C), g = 1/(1 + a), the Schmiedl-Seifert form.
     Entropy-type dissipation (a = sqrt(rho tau)): symmetric case is exactly Curzon-Ahlborn (Esposito et al. 2010).
  C. The heat-type family lies on the plane dial at k = 2g - 2, d = -g(2 - g): the parabola d = k^2/4 - 1,
     from parity (-2, 0) (lower bound eta_C/2) to the velocity point (0, -1) (upper bound eta_C/(2 - eta_C)).
  D. Other finite-time efficiencies (entropy-type low-dissipation at rho != 1, ecological, Feynman ratchet) are off it.
  E. A fine cascade of max-power stages of any of these kinds tends to Curzon-Ahlborn
     (known for linear-response stages: Van den Broeck 2005). Splitting raises the lower edge, lowers the upper.
Needs sympy, numpy, mpmath.
"""
import json
import numpy as np
import sympy as sp
import mpmath as mp

results = {}
def ok(name, cond):
    results[name] = bool(cond); print(f"  [{'PASS' if cond else 'FAIL'}] {name}"); return cond

tau, g, a, rho, lam = sp.symbols('tau g a rho lam', positive=True)
eC = 1 - tau

print("== A. Knob images of Carnot ==")
union = lambda x, y: x + y - x*y
t1, t2 = sp.symbols('t1 t2', positive=True)
knob = lambda t: 1 - t**lam
ok("A1: eta = 1 - tau^lam composes by the union rule across a cascade (any lam)",
   sp.simplify(union(knob(t1), knob(t2)) - knob(t1*t2)) == 0)
ok("A2: Curzon-Ahlborn 1 - sqrt(tau) is lam = 1/2; two half-steps make Carnot",
   sp.simplify(union(1 - sp.sqrt(tau), 1 - sp.sqrt(tau)) - eC) == 0)
# cascade-stable means f(t1 t2) lifts multiply: 1 - eta is a continuous multiplicative map, so a power (Cauchy).
for name, f in [("eta_C/2", eC/2), ("eta_C/(2 - eta_C)", eC/(2 - eC)), ("eta_C/(2 - eta_C/2)", eC/(2 - eC/2))]:
    F = sp.Lambda(tau, f)
    gap = sp.simplify(union(F(t1), F(t2)) - F(t1*t2))
    results[f"_A3_{name}"] = gap != 0
ok("A3: eta_C/2, eta_C/(2 - eta_C), eta_C/(2 - eta_C/2) are NOT cascade-stable (not knob images)",
   all(results.pop(k) for k in list(results) if k.startswith("_A3_")))

print("\n== B. Low-dissipation engines at maximum power, derived ==")
Th, S, sh, x = sp.symbols('T_h DeltaS Sigma_h x', positive=True)
Tc = Th*tau
# optimum: t_c = a t_h with a = sqrt(T_c Sigma_c/(T_h Sigma_h)), t_h = 2(1 + a) T_h Sigma_h/((T_h - T_c) S)
sc = a**2*Th*sh/Tc
th = 2*(1 + a)*Th*sh/((Th - Tc)*S); tc = a*th
Qh = Th*(S - sh/th); Qc = Tc*(S + sc/tc)
P = lambda th_, tc_: (Th*(S - sh/th_) - Tc*(S + sc/tc_))/(th_ + tc_)
u, v = sp.symbols('u v', positive=True)
Puv = P(u, v)
stat = [sp.simplify(sp.diff(Puv, u).subs({u: th, v: tc})), sp.simplify(sp.diff(Puv, v).subs({u: th, v: tc}))]
ok("B1: the stated times are a stationary point of the power", all(s_ == 0 for s_ in stat))
eta_star = sp.simplify(1 - Qc/Qh)
ok("B2: eta* = eta_C (1 + a)/(2(1 + a) - eta_C)", sp.simplify(eta_star - eC*(1 + a)/(2*(1 + a) - eC)) == 0)
# numeric maximum (grid refinement) agrees, for entropy-type rho at three temperatures
def num_max(etaC, rho_):
    Tc_ = 1 - etaC; lo = np.array([-4., -4.]); hi = np.array([4., 4.])
    for _ in range(9):
        A_ = np.linspace(lo[0], hi[0], 121); B_ = np.linspace(lo[1], hi[1], 121)
        U, V = np.meshgrid(10**A_, 10**B_, indexing='ij')
        Qh_ = 1 - 1/U; Qc_ = Tc_*(1 + rho_/V); Pw = (Qh_ - Qc_)/(U + V)
        i, j = np.unravel_index(np.argmax(Pw), Pw.shape)
        c = np.array([A_[i], B_[j]]); w = (hi - lo)/10; lo, hi = c - w, c + w
    return 1 - Qc_[i, j]/Qh_[i, j]
errs = []
for etaC in (0.3, 0.6, 0.9):
    for r in (0.1, 1.0, 10.0):
        aa = np.sqrt((1 - etaC)*r)
        errs.append(abs(num_max(etaC, r) - etaC*(1 + aa)/(2*(1 + aa) - etaC)))
ok("B3: brute-force maximum of the power agrees to 1e-6 (9 cases)", max(errs) < 1e-6)
ok("B4: entropy-type, symmetric (rho = 1, a = sqrt(tau)): exactly Curzon-Ahlborn",
   sp.simplify(eta_star.subs(a, sp.sqrt(tau)) - (1 - sp.sqrt(tau))) == 0)
ok("B5: heat-type (a constant): eta* = eta_C/(2 - g eta_C), g = 1/(1 + a)",
   sp.simplify(eta_star - eC/(2 - eC/(1 + a))) == 0)
ok("B6: g = 0 gives eta_C/2, g = 1 gives eta_C/(2 - eta_C) (the two bounds)",
   sp.simplify((eC/(2 - g*eC)).subs(g, 0) - eC/2) == 0 and sp.simplify((eC/(2 - g*eC)).subs(g, 1) - eC/(2 - eC)) == 0)

print("\n== C. The heat-type family on the plane dial ==")
eta = sp.symbols('eta')
f = eC/(2 - g*eC)
deta_dbeta = sp.simplify(-tau*sp.diff(f, tau))
tau_of_eta = sp.solve(sp.Eq(f, eta), tau)[0]
rhs = sp.factor(sp.simplify(deta_dbeta.subs(tau, tau_of_eta)))
print("   d eta/d beta =", rhs)
target = sp.Rational(1, 2)*(1 + (2*g - 2)*eta - g*(2 - g)*eta**2)
ok("C1: d eta/d beta = (1/2)(1 + (2g - 2) eta - g(2 - g) eta^2): a plane law for every g", sp.simplify(rhs - target) == 0)
kk = 2*g - 2; dd = -g*(2 - g)
ok("C2: the family is the parabola d = k^2/4 - 1 (discriminant k^2 - 4d = 4 for every g)", sp.simplify(dd - (kk**2/4 - 1)) == 0)
ok("C3: ends: g = 0 at parity (k, d) = (-2, 0); g = 1 at the velocity point (0, -1); g = 1/2 at (-1, -3/4)",
   (kk.subs(g, 0), dd.subs(g, 0)) == (-2, 0) and (kk.subs(g, 1), dd.subs(g, 1)) == (0, -1)
   and (kk.subs(g, sp.Rational(1, 2)), dd.subs(g, sp.Rational(1, 2))) == (-1, -sp.Rational(3, 4)))
# Curzon-Ahlborn: d eta/d beta = (1/2)(1 - eta): k = -1, d = 0, same rate 1/2 as the family
ca = 1 - sp.sqrt(tau)
ok("C4: Curzon-Ahlborn: d eta/d beta = (1 - eta)/2, the union line at half speed (not on the parabola)",
   sp.simplify(-tau*sp.diff(ca, tau) - (1 - ca)/2) == 0)
# each family member is a Moebius map of Carnot's lift tau
ok("C5: every family member is a Moebius function of tau", sp.degree(sp.numer(sp.together(f)), tau) == 1
   and sp.degree(sp.denom(sp.together(f)), tau) == 1)
# composition law of the lower edge is parity, of the upper edge velocity addition
b1, b2 = sp.symbols('b1 b2', positive=True)
lo_ = lambda B: (1 - sp.exp(-B))/2; up_ = lambda B: sp.tanh(B/2)
ok("C6: lower edge composes by parity x + y - 2xy; upper edge by velocity addition (x + y)/(1 + xy)",
   sp.simplify(lo_(b1) + lo_(b2) - 2*lo_(b1)*lo_(b2) - lo_(b1 + b2)) == 0
   and sp.simplify(((up_(b1) + up_(b2))/(1 + up_(b1)*up_(b2)) - up_(b1 + b2)).rewrite(sp.exp)) == 0
   and sp.simplify((eC/(2 - eC)) - sp.tanh(-sp.log(tau)/2).rewrite(sp.exp)) == 0)

print("\n== D. Riccati test: which finite-time efficiencies are plane laws? ==")
mp.mp.dps = 40
def riccati_residual(expr):
    e = sp.lambdify(tau, expr, 'mpmath')
    ts = [mp.mpf(s_) for s_ in ('0.05', '0.2', '0.35', '0.5', '0.65', '0.8', '0.9')]
    xs = [e(t_) for t_ in ts]; ds = [-t_*mp.diff(e, t_) for t_ in ts]
    A_ = mp.matrix([[1, x_, x_**2] for x_ in xs[:3]]); c = mp.lu_solve(A_, mp.matrix(ds[:3]))
    return float(max(abs(ds[i] - (c[0] + c[1]*xs[i] + c[2]*xs[i]**2)) for i in range(3, 7)))
cases = {
    "Carnot": eC, "Curzon-Ahlborn": 1 - sp.sqrt(tau),
    "lower bound eta_C/2": eC/2, "upper bound eta_C/(2-eta_C)": eC/(2 - eC),
    "heat-type, symmetric eta_C/(2-eta_C/2)": eC/(2 - eC/2),
    "entropy-type rho = 0.1": eC*(1 + sp.sqrt(tau/10))/(2*(1 + sp.sqrt(tau/10)) - eC),
    "entropy-type rho = 10": eC*(1 + sp.sqrt(10*tau))/(2*(1 + sp.sqrt(10*tau)) - eC),
    "ecological (Angulo-Brown 1991)": 1 - sp.sqrt(tau*(1 + tau)/2),
    "Feynman ratchet (Tu 2008)": eC**2/(eC - tau*sp.log(tau)),
}
res = {k_: riccati_residual(v_) for k_, v_ in cases.items()}
for k_, v_ in res.items(): print(f"   {k_:42s} residual {v_:.1e}")
on = [k_ for k_, v_ in res.items() if v_ < 1e-30]; off = [k_ for k_, v_ in res.items() if v_ > 1e-4]
ok("D1: on the dial: Carnot, Curzon-Ahlborn, both bounds, the heat-type symmetric case", len(on) == 5)
ok("D2: off the dial (residual > 1e-4): entropy-type rho != 1, ecological, Feynman ratchet", len(off) == 4)
results["_D_residuals"] = res

print("\n== E. Fine cascades of max-power stages ==")
def cascade(fn, etaC, n):
    B = -np.log(1 - etaC); stage = fn(np.exp(-B/n)); return 1 - (1 - stage)**n
fam = {"lower eta_C/2": lambda t_: (1 - t_)/2, "upper eta_C/(2-eta_C)": lambda t_: (1 - t_)/(1 + t_),
       "heat-type symmetric": lambda t_: (1 - t_)/(2 - (1 - t_)/2),
       "entropy-type rho=10": lambda t_: (1 - t_)*(1 + np.sqrt(10*t_))/(2*(1 + np.sqrt(10*t_)) - (1 - t_)),
       "Feynman ratchet": lambda t_: (1 - t_)**2/((1 - t_) - t_*np.log(t_))}
etaC = 0.8; ca_v = 1 - np.sqrt(1 - etaC)
lim = {k_: cascade(fn, etaC, 100000) for k_, fn in fam.items()}
for k_, v_ in lim.items(): print(f"   {k_:24s} n=1: {cascade(fam[k_], etaC, 1):.5f}  n=10^5: {v_:.6f}   (Curzon-Ahlborn {ca_v:.6f})")
ok("E1: every family tends to Curzon-Ahlborn under fine cascading (to 1e-4 at n = 10^5)",
   all(abs(v_ - ca_v) < 1e-4 for v_ in lim.values()))
ok("E2: splitting raises the lower edge and lowers the upper edge (n = 2 vs n = 1)",
   cascade(fam["lower eta_C/2"], etaC, 2) > cascade(fam["lower eta_C/2"], etaC, 1)
   and cascade(fam["upper eta_C/(2-eta_C)"], etaC, 2) < cascade(fam["upper eta_C/(2-eta_C)"], etaC, 1))
ok("E3: Curzon-Ahlborn itself does not move under any split (n = 1, 2, 7, 50)",
   all(abs(cascade(lambda t_: 1 - np.sqrt(t_), etaC, n) - ca_v) < 1e-12 for n in (1, 2, 7, 50)))
# why: the limit is 1 - exp(-beta * eta'(beta = 0)), and every max-power efficiency here has slope 1/2 at eta_C -> 0
slopes = {k_: float(sp.limit(sp.diff(v_.subs(tau, sp.exp(-sp.Symbol('B'))), sp.Symbol('B')), sp.Symbol('B'), 0))
          for k_, v_ in cases.items() if k_ != "Carnot"}
print("   slopes d eta/d beta at beta = 0:", {k_: round(v_, 6) for k_, v_ in slopes.items()})
ok("E4: the reason: every max-power efficiency has slope 1/2 at small eta_C (linear response)",
   all(abs(s_ - 0.5) < 1e-12 for k_, s_ in slopes.items() if not k_.startswith("ecological")))
ok("E5: the ecological efficiency (not a max-power point) has slope 3/4, so its fine cascade tends to 1 - tau^(3/4)",
   abs(slopes["ecological (Angulo-Brown 1991)"] - 0.75) < 1e-12
   and abs(cascade(lambda t_: 1 - np.sqrt(t_*(1 + t_)/2), etaC, 100000) - (1 - (1 - etaC)**0.75)) < 1e-4)

n_pass = sum(v_ for k_, v_ in results.items() if not k_.startswith('_')); n_all = sum(1 for k_ in results if not k_.startswith('_'))
print(f"\n{n_pass}/{n_all} checks pass")
json.dump(results, open(__file__.replace('checks.py', 'checks.json'), 'w'), indent=1)
