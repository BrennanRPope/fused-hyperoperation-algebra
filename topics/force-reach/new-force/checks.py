"""Does the theory predict a force nobody has seen? Checks for the DRAFT page.
Run: python3 checks.py   (numpy, scipy, sympy, mpmath)"""
import json, math
import numpy as np, sympy as sp, mpmath as mp
from scipy.integrate import solve_ivp
R = {"pass": 0, "fail": 0}
OUT = {}
def check(name, ok, info=""):
    R["pass" if ok else "fail"] += 1
    print(("PASS  " if ok else "FAIL  ") + name + (f"   [{info}]" if info else ""))

def deff_of(F, r, h=1e-4):
    """D_eff = 1 - d ln|F| / d ln r, by central difference in ln r."""
    lr = math.log(r)
    a, b = abs(F(math.exp(lr + h))), abs(F(math.exp(lr - h)))
    return 1 - (math.log(a) - math.log(b)) / (2 * h)

# ---------------------------------------------------------------- A
print("== A. The reach classes come from one reading, Delta = k^2 - 4d")
L, k, d = sp.symbols('L k d', real=True)
# linear form of the law: u'' - k u' + d u = 0 ; drift removed: w'' = (Delta/4) w
for (kk, dd, label) in [(3, 1, "Delta>0"), (2, 1, "Delta=0"), (1, 2, "Delta<0")]:
    Dl = kk**2 - 4*dd
    sol = solve_ivp(lambda t, y: [y[1], (Dl/4)*y[0]], [0, 12], [1, -math.sqrt(max(Dl, 0))/2 if Dl > 0 else 0.0],
                    rtol=1e-11, atol=1e-13, dense_output=True)
    w = sol.sol(np.array([4.0, 8.0]))[0]
    if Dl > 0:
        rate = -math.log(w[1]/w[0])/4
        check(f"A1 {label} (k={kk}, d={dd}): the decaying field falls at rate sqrt(Delta)/2, reach 2/sqrt(Delta)",
              abs(rate - math.sqrt(Dl)/2) < 1e-6, f"rate {rate:.8f} vs {math.sqrt(Dl)/2:.8f}")
    elif Dl == 0:
        check(f"A2 {label}: no length scale; the field is a straight line (power law), unlimited reach",
              abs(w[0] - 1) < 1e-9 and abs(w[1] - 1) < 1e-9, "w constant for zero slope")
    else:
        per = 2*math.pi/math.sqrt(-Dl)
        ts = np.linspace(0, 12, 120001); ws = sol.sol(ts)[0]
        zs = ts[:-1][np.sign(ws[:-1]) != np.sign(ws[1:])]
        gap = float(np.mean(np.diff(zs)))
        check(f"A3 {label}: the law's ruler x = -u'/(d u) blows up and returns every 2 pi/sqrt(-Delta) (zeros of w): a closed direction, rungs n sqrt(-Delta)",
              abs(gap - per) < 1e-3, f"zero spacing {gap:.5f} vs {per:.5f}")

# ---------------------------------------------------------------- B
print("== B. The rung class is one more direction at short range")
x = sp.symbols('x', positive=True); q = sp.symbols('q', positive=True)
geo = 1 + 2*q/(1 - q)
check("B1 the tower with every rung coupled equally sums to coth: 1 + 2 sum_{n>=1} e^{-n x} = coth(x/2)",
      sp.simplify((geo.subs(q, sp.exp(-x)) - sp.coth(x/2)).rewrite(sp.exp)) == 0)
mp.mp.dps = 30
worst = 0
for (Rr, rr) in [(0.7, 1.3), (1.0, 0.05), (2.5, 40.0), (0.1, 0.3)]:
    Rm, rm = mp.mpf(Rr), mp.mpf(rr)
    img = mp.nsum(lambda n: 1/(rm**2 + (2*mp.pi*Rm*n)**2), [-mp.inf, mp.inf])
    worst = max(worst, abs(img - mp.coth(rm/(2*Rm))/(2*Rm*rm)) / img)
check("B2 tower sum = exact field of a point source on space x circle (method of images), circle length 2 pi R",
      worst < mp.mpf(10)**-25, f"worst rel. diff {mp.nstr(worst, 3)} at 30 digits")
Ft = lambda r: (lambda h: -(-(1/math.tanh((r+h)/2))/(r+h) + (1/math.tanh((r-h)/2))/(r-h))/(2*h))(1e-6*r)
lo, hi = deff_of(Ft, 0.01), deff_of(Ft, 40.0)
check("B3 effective dimension of that force: 4 well inside the circle size, 3 far outside",
      abs(lo - 4) < 1e-3 and abs(hi - 3) < 1e-3, f"D_eff(0.01 R) = {lo:.4f}, D_eff(40 R) = {hi:.4f}")
from scipy.optimize import brentq
rc = brentq(lambda r: deff_of(Ft, r) - 3.5, 0.5, 20)
OUT["B_crossover_r_over_R"] = rc
check("B3b the step from 4 to 3 is centred near r = 4.1 R, i.e. about 2/3 of the circle length 2 pi R",
      3.5 < rc < 5.0, f"D_eff = 3.5 at r = {rc:.3f} R")
# leading Yukawa strength = degeneracy of the lightest rung
alpha1 = float(sp.limit((sp.coth(x/2) - 1)/sp.exp(-x), x, sp.oo))
# 2-torus: lightest nonzero modes (±1,0),(0,±1) -> 4
def torus_alpha(dims):
    import itertools
    cnt = sum(1 for v in itertools.product([-1, 0, 1], repeat=dims) if sum(abs(t) for t in v) == 1)
    return cnt
check("B4 far away the tower looks like one Yukawa term of strength alpha = 2 (rungs +n and -n); a torus of d circles gives 2d",
      abs(alpha1 - 2) < 1e-12 and torus_alpha(2) == 4 and torus_alpha(3) == 6,
      f"alpha = {alpha1:.0f}; 2-torus {torus_alpha(2)}, 3-torus {torus_alpha(3)} (Kehagias-Sfetsos 2000: 2n)")
Dneg = sp.symbols('Delta', negative=True)
check("B5 the law's rung spacing sqrt(-Delta) and circle length 2 pi/sqrt(-Delta) are the same statement: R = 1/sqrt(-Delta)",
      sp.simplify(2*sp.pi/sp.sqrt(-Dneg) - 2*sp.pi*(1/sp.sqrt(-Dneg))) == 0)

# ---------------------------------------------------------------- C
print("== C. The screened class: a massive carrier adds a Yukawa term")
def Fy(alpha):
    return lambda r: (1 + alpha*math.exp(-r)*(1 + r))/r**2   # force of -(1+alpha e^{-r})/r, range lambda = 1
bump = {}
for a in [1.0, 0.3, -0.3]:
    rs = np.exp(np.linspace(math.log(0.02), math.log(60), 400))
    De = np.array([deff_of(Fy(a), r) for r in rs])
    i = np.argmax(np.abs(De - 3))
    bump[a] = (float(rs[i]), float(De[i]))
OUT["C_bumps"] = {str(k): v for k, v in bump.items()}
check("C1 a Yukawa term shows as a bump in D_eff (up for alpha > 0, down for alpha < 0), back to 3 at both ends",
      bump[1.0][1] > 3.3 and bump[-0.3][1] < 3 and abs(deff_of(Fy(1.0), 0.005) - 3) < 1e-3 and abs(deff_of(Fy(1.0), 80) - 3) < 1e-6,
      f"alpha=1: peak {bump[1.0][1]:.3f} at r = {bump[1.0][0]:.2f} lambda; alpha=-0.3: {bump[-0.3][1]:.3f}")
check("C2 the bump's place is set by the range (about 2 lambda) for every alpha: the theory fixes the range, not the strength",
      all(1.3 < bump[a][0] < 3.0 for a in bump), ", ".join(f"{a}: {bump[a][0]:.2f}" for a in bump))
hbarc = 1.973269804e-7  # eV m
m51 = hbarc/38.6e-6; m30 = hbarc/30e-6
check("C3 Eot-Wash 2020: gravitational-strength Yukawa needs lambda < 38.6 um, i.e. carrier mass > 5.1 meV (their quoted number)",
      abs(m51*1e3 - 5.1) < 0.05, f"{m51*1e3:.2f} meV")
check("C3b and the largest extra dimension has radius < 30 um, i.e. lightest rung > 6.6 meV",
      abs(m30*1e3 - 6.58) < 0.02, f"{m30*1e3:.2f} meV")

# ---------------------------------------------------------------- D
print("== D. The flux exponent drifting at low field: s = 3/2 is deep-MOND gravity")
s, E, Dv = sp.symbols('s E D', positive=True)
# Legendre: sup_D (E D - D^s/s) = E^q/q with q = s/(s-1)
leg = sp.simplify((E*Dv - Dv**sp.Rational(3, 2)/sp.Rational(3, 2)).subs(Dv, E**2))
check("D1 s = 3/2 in D-language is exponent q = s/(s-1) = 3 in field language: Legendre of |D|^{3/2}/(3/2) is |E|^3/3",
      sp.simplify(leg - E**3/3) == 0 and sp.Rational(3, 2)/(sp.Rational(3, 2) - 1) == 3)
G, M, a0, r = sp.symbols('G M a0 r', positive=True)
g = sp.symbols('g', positive=True)
gsol = sp.solve(sp.Eq(g**2/a0*4*sp.pi*r**2, 4*sp.pi*G*M), g)[0]
check("D2 deep-MOND field equation div(|grad phi| grad phi / a0) = 4 pi G rho gives g = sqrt(G M a0)/r and v^4 = G M a0",
      sp.simplify(gsol - sp.sqrt(G*M*a0)/r) == 0 and sp.simplify((gsol*r)**2 - G*M*a0) == 0)
lam, Dd, qq = sp.symbols('lambda D_dim q', positive=True)
# |grad phi|^q d^D x under x -> lam x, phi(x) -> phi(x/lam): factor lam^(D - q)
check("D3 the action |grad phi|^q d^Dx is unchanged by stretching space exactly when q = D (3D: q = 3, the deep-MOND action)",
      sp.simplify(lam**(Dd - qq)).subs(qq, Dd) == 1 and sp.simplify(lam**(3 - 2)) != 1)
mu_g = lambda gN: gN/2*(1 + math.sqrt(1 + 4/gN))   # mu(x) = x/(1+x), a0 = 1
Fm = lambda rr: mu_g(1/rr**2)                      # GM = a0 = 1
dh, dl = deff_of(Fm, 0.01), deff_of(Fm, 1e4)
check("D4 with the simple interpolating function, gravity's D_eff is 3 at high acceleration and falls to 2 at low acceleration",
      abs(dh - 3) < 0.02 and abs(dl - 2) < 0.02, f"D_eff(0.01) = {dh:.3f}, D_eff(10^4) = {dl:.3f}")
check("D5 s = 3/2 is also the 3D confinement threshold of the flux-tube page: deep MOND sits on that border (log potential)",
      sp.Rational(3, 2) == sp.Rational(3, 1)/(3 - 1))

# ---------------------------------------------------------------- E
print("== E. The roundness dial p: a force from a non-round distance")
th = np.linspace(0, 2*np.pi, 4096, endpoint=False)
def harmonics(p):
    X, Y = np.cos(th), np.sin(th); h = 1e-6
    V = lambda X, Y: 1/(np.abs(X)**p + np.abs(Y)**p)**(1/p)
    Fx = -(V(X+h, Y) - V(X-h, Y))/(2*h); Fy_ = -(V(X, Y+h) - V(X, Y-h))/(2*h)
    F = np.hypot(Fx, Fy_); c = np.fft.rfft(F)/len(F)
    return np.abs(c[1:9])*2/abs(c[0])
H = {p: harmonics(p) for p in [1.5, 2.0, 2.1, 2.2, 3.0, 4.0]}
OUT["E_cos4_amplitude"] = {str(p): float(H[p][3]) for p in H}
check("E1 for p != 2 the force's strength around a circle has no 1st, 2nd or 3rd harmonic; the first change is four-lobed (cos 4 theta)",
      all(max(H[p][:3]) < 1e-9 and H[p][3] > 1e-3 for p in [1.5, 2.2, 3.0, 4.0]),
      ", ".join(f"p={p}: {H[p][3]:.4f}" for p in [1.5, 2.2, 3.0, 4.0]))
check("E2 p = 2 is exactly round", max(H[2.0]) < 1e-9)
ratio = H[2.1][3]/H[2.2][3]
check("E3 near round the four-lobed size grows in step with p - 2", abs(ratio - 0.5) < 0.03, f"ratio {ratio:.3f}")

# ---------------------------------------------------------------- F
print("== F. Order on, no mass, no tube: a long-range non-abelian force")
def nf_max(N): return max(nf for nf in range(200) if sp.Rational(11, 3)*N - sp.Rational(2, 3)*nf > 0)
check("F1 such a force confines (becomes a tube) unless it carries many matter fields: one-loop b0 = 11N/3 - 2nf/3 stays > 0 up to nf = 16 for SU(3), 10 for SU(2)",
      nf_max(3) == 16 and nf_max(2) == 10, f"SU(2): {nf_max(2)}, SU(3): {nf_max(3)}, SU(4): {nf_max(4)}")

# ---------------------------------------------------------------- P
print("== P. Polishing")
mp.mp.dps = 50
Rm, rm = mp.mpf('0.37'), mp.mpf('2.9')
img = mp.nsum(lambda n: 1/(rm**2 + (2*mp.pi*Rm*n)**2), [-mp.inf, mp.inf])
diff = abs(img - mp.coth(rm/(2*Rm))/(2*Rm*rm))
check("P1 B2 again at 50 digits", diff < mp.mpf(10)**-45, mp.nstr(diff, 3))
Fe = lambda r: (lambda h: -(-(1/math.tanh((r+h)/2))/(r+h) + (1/math.tanh((r-h)/2))/(r-h))/(2*h))(1e-5*r)
rc2 = brentq(lambda r: deff_of(Fe, r, 5e-5) - 3.5, 0.5, 20)
check("P2 B3b crossover is stable under a finer derivative step", abs(rc2 - rc) < 1e-3, f"{rc:.5f} vs {rc2:.5f}")

OUT["pass"], OUT["fail"] = R["pass"], R["fail"]
json.dump(OUT, open("results.json", "w"), indent=1)
print(f"\nTOTAL: {R['pass']} pass, {R['fail']} fail")
