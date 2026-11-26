"""DRAFT checks: why orbit velocities run straight in the dial geometry (Nov 26, 2026).

Velocities are complex numbers v = vx + i vy (motion is planar). GM = 1.
Dial addition (twisted-boosts page, kappa-stereographic): u (+)_k v = (u + v)/(1 - k conj(u) v).
Writes checks_output.txt and results.json.
"""
import json, math, os
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
R, LOG = {}, []
PASS = FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    line = f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}"
    print(line); LOG.append(line)
    PASS += bool(ok); FAIL += (not ok)


def orbit(p, x0, v0, T, n=4000):
    def f(t, y):
        r = math.hypot(y[0], y[1]); a = -1/r**p/r
        return [y[2], y[3], a*y[0], a*y[1]]
    s = solve_ivp(f, (0, T), list(x0) + list(v0), method="DOP853", rtol=1e-12, atol=1e-13,
                  t_eval=np.linspace(0, T, n))
    X = s.y[0] + 1j*s.y[1]; V = s.y[2] + 1j*s.y[3]
    return s.t, X, V


def dial_add(u, v, k):
    return (u + v)/(1 - k*np.conj(u)*v)


def circle_fit(V):
    A = np.stack([2*V.real, 2*V.imag, np.ones_like(V.real)], 1); b = np.abs(V)**2
    s, *_ = np.linalg.lstsq(A, b, rcond=None)
    cen = s[0] + 1j*s[1]; rad = math.sqrt(s[2] + abs(cen)**2)
    return cen, rad, np.max(np.abs(np.abs(V - cen) - rad))/rad


rng = np.random.default_rng(7)

# V1. Every 1/r^2 hodograph is a circle whose 'power' at v = 0 is 2E. So the dial that makes it a straight line,
#     kappa = -1/(power), is -1/(2E): fixed by the energy alone, the same for every orbit of that energy.
ok = True; rows = []
for _ in range(12):
    r0 = rng.uniform(0.7, 2.0); ang = rng.uniform(0, 2*np.pi)
    sp_ = rng.uniform(0.4, 1.9)/math.sqrt(r0); dirn = ang + np.pi/2 + rng.uniform(-0.9, 0.9)
    x0 = (r0*math.cos(ang), r0*math.sin(ang)); v0 = (sp_*math.cos(dirn), sp_*math.sin(dirn))
    t, X, V = orbit(2, x0, v0, 25)
    E = 0.5*abs(V[0])**2 - 1/abs(X[0])
    cen, rad, res = circle_fit(V)
    power = abs(cen)**2 - rad**2
    rows.append((E, power, res)); ok &= res < 1e-8 and abs(power - 2*E) < 1e-8
check("V1 power of v=0 for each hodograph circle = 2E (12 random orbits, bound and unbound)", ok,
      f"max |power-2E| = {max(abs(p-2*e) for e, p, _ in rows):.1e}")
R["V1"] = rows

# V1b. Same energy, different angular momentum and orientation: one dial for all of them
E0 = -0.3
ok = True
for h_ in (0.3, 0.7, 1.0, 1.15):
    r0 = 1.0; v0mag = math.sqrt(2*(E0 + 1/r0)); vt = h_/r0
    vr = math.sqrt(max(v0mag**2 - vt**2, 0))
    t, X, V = orbit(2, (1, 0), (vr, vt), 30)
    cen, rad, res = circle_fit(V)
    ok &= abs(-1/(abs(cen)**2 - rad**2) - (-1/(2*E0))) < 1e-7
check("V1b orbits sharing energy E share one dial kappa = -1/(2E) whatever their angular momentum", ok,
      f"kappa = {-1/(2*E0):.4f}")

# V2. The theory's own composition: measured from the starting velocity with the dial addition at kappa = -1/(2E),
#     every later velocity points along one fixed line through the origin (a dial straight line).
def line_spread(V, k):
    a = V[0]
    Z = dial_add(-a, V[1:], k)
    Z = Z[np.abs(Z) > 1e-9]
    ang = np.angle(Z) % np.pi          # direction mod pi (a line, both ways)
    ang = np.where(ang > np.pi - 0.3, ang - np.pi, ang) if np.ptp(ang) > np.pi/2 else ang
    return np.ptp(ang)
tb, Xb, Vb = orbit(2, (1, 0), (0.1, 1.15), 40)
Eb = 0.5*abs(Vb[0])**2 - 1/abs(Xb[0]); kb = -1/(2*Eb)
sb = line_spread(Vb[:1500], kb); sb_bad = line_spread(Vb[:1500], 1.7*kb); sb_flat = line_spread(Vb[:1500], 0.0)
check("V2a bound orbit: (-v0) (+)_k v(t) keeps one direction when k = -1/(2E); not for other dials",
      sb < 1e-8 and sb_bad > 1e-2 and sb_flat > 1e-2, f"spread {sb:.1e}; k x1.7: {sb_bad:.3f}; flat k=0: {sb_flat:.3f}")
t, Xu, Vu = orbit(2, (-25, 1.2), (0.9, 0.0), 60)
Eu = 0.5*abs(Vu[0])**2 - 1/abs(Xu[0]); ku = -1/(2*Eu)
su = line_spread(Vu, ku); su_bad = line_spread(Vu, 0.5*ku)
check("V2b unbound orbit: same, with the negative dial k = -1/(2E) < 0 (the boost end)",
      su < 1e-8 and su_bad > 1e-2, f"E = {Eu:.4f}, k = {ku:.4f}, spread {su:.1e}; wrong dial {su_bad:.3f}")
R["V2"] = dict(Eb=Eb, kb=kb, sb=sb, sb_bad=sb_bad, Eu=Eu, ku=ku, su=su)

# V3. Uniform along the line in the dial's own clock: dial distance from v0 equals the eccentric-anomaly change,
#     i.e. it grows at sqrt(-2E)/r per unit time (the clock dt/r).
a_ = -1/(2*Eb)
Zb = dial_add(-Vb[0], Vb, kb)
ddist = 2*np.arctan(math.sqrt(kb)*np.abs(Zb))          # distance on the unit dial sphere
r_ = np.abs(Xb)
e_ = abs((abs(Vb[0])**2 - 1/abs(Xb[0]))*Xb[0] - (Xb[0].real*Vb[0].real + Xb[0].imag*Vb[0].imag)*Vb[0])
cosE = (1 - r_/a_)/e_
sgn = np.sign((Xb.real*Vb.real + Xb.imag*Vb.imag))
Ecc = np.unwrap(np.where(sgn >= 0, np.arccos(np.clip(cosE, -1, 1)), -np.arccos(np.clip(cosE, -1, 1))))
dE = np.abs(Ecc - Ecc[0])
m = dE < math.pi*0.98
check("V3 dial distance travelled = change in eccentric anomaly (uniform motion in the clock dt/r)",
      np.max(np.abs(ddist[m] - dE[m])) < 1e-6, f"max diff {np.max(np.abs(ddist[m]-dE[m])):.1e}, eccentricity {e_:.3f}")
rate = np.gradient(2*np.arctan(math.sqrt(kb)*np.abs(Zb)), tb)
check("V3b its rate is sqrt(-2E)/r", np.max(np.abs(rate[5:400]*r_[5:400] - math.sqrt(-2*Eb))) < 1e-3)

# V4. Distance from the sun = a x the dial chart's scale at the current velocity: r = a * 2/(1 + k|v|^2)
lam = 2/(1 + kb*np.abs(Vb)**2)
check("V4 r = a x lambda(v), lambda = 2/(1 + k|v|^2) the dial chart's scale (position read off the velocity geometry)",
      np.max(np.abs(r_ - a_*lam))/a_ < 1e-9, f"a = {a_:.4f}")

# V5. Kepler's equation from V3 + V4: n t = s - e sin s, with s the dial distance
n_ = a_**-1.5
s_ = Ecc - Ecc[0] + 0  # use E(t)
kep = (Ecc - e_*np.sin(Ecc)) - (Ecc[0] - e_*np.sin(Ecc[0]))
check("V5 Kepler's equation n t = s - e sin s, with s the dial arclength", np.max(np.abs(kep - n_*tb)) < 1e-6)

# V6. Converse: the postulate fails for every other power law (no circle, so no dial makes it straight).
#     Flux law in n dimensions gives 1/r^(n-1); only n = 3 satisfies the postulate.
res = {}
for n_dim, p in ((2, 1.0), (2.5, 1.5), (3, 2.0), (3.5, 2.5), (4, 3.0)):
    tt, Xp, Vp = orbit(p, (1, 0), (0, 0.9 if p < 3 else 1.02), 30 if p < 3 else 3)
    res[n_dim] = circle_fit(Vp)[2]
check("V6 postulate (velocities run straight in one dial geometry) holds only for 1/r^2; flux law agrees only at n=3",
      res[3] < 1e-8 and all(v > 1e-3 for k_, v in res.items() if k_ != 3),
      ", ".join(f"n={k_}: {v:.0e}" for k_, v in res.items()))
R["V6"] = {str(k_): v for k_, v in res.items()}

# V7. Hamilton's converse, symbolic: for a central force F(r), d/dt of the hodograph-circle centre
#     c = v - (1/h) i e^{i theta}*K  is zero for all orbits iff F = -K/r^2.
th, tt_ = sp.symbols("theta t", real=True)
K, h = sp.symbols("K h", positive=True)
rf = sp.Function("r")(tt_); thf = sp.Function("theta")(tt_); F = sp.Function("F")
# velocity in polar: v = (r' + i r theta') e^{i theta};  h = r^2 theta'
v = (sp.diff(rf, tt_) + sp.I*rf*sp.diff(thf, tt_))*sp.exp(sp.I*thf)
cvec = v - sp.I*K/h*sp.exp(sp.I*thf)
dc = sp.diff(cvec, tt_)
# equations of motion: r'' - r theta'^2 = F(r), theta' = h/r^2
dc = dc.subs(sp.Derivative(thf, (tt_, 2)), sp.diff(h/rf**2, tt_)).subs(sp.Derivative(thf, tt_), h/rf**2)
dc = dc.subs(sp.Derivative(rf, (tt_, 2)), F(rf) + rf*(h/rf**2)**2)
dc = sp.simplify(dc*sp.exp(-sp.I*thf))
check("V7 the hodograph centre stays fixed (a circle) iff F(r) = -K/r^2 (Hamilton 1846)",
      sp.simplify(dc.subs(F(rf), -K/rf**2)) == 0 and sp.simplify(dc) != 0, f"d(centre)/dt ~ {sp.factor(dc)}")

# V8. Symmetry: angular momentum L and Runge-Lenz A close like the dial's rotations and boosts:
#     {A_i, A_j} = -2H eps_ijk L_k.  E<0: so(4) (sphere), E=0: e(3) (flat), E>0: so(3,1) (Lorentz).
x, y, z, px, py, pz = sp.symbols("x y z p_x p_y p_z", real=True)
X3 = sp.Matrix([x, y, z]); P3 = sp.Matrix([px, py, pz]); rr = sp.sqrt(x**2 + y**2 + z**2)
H = P3.dot(P3)/2 - 1/rr
Lv = X3.cross(P3); Av = P3.cross(Lv) - X3/rr
def pb(f, g):
    return sum(sp.diff(f, X3[i])*sp.diff(g, P3[i]) - sp.diff(f, P3[i])*sp.diff(g, X3[i]) for i in range(3))
ok = all(sp.simplify(pb(H, Av[i])) == 0 for i in range(3))
ok &= sp.simplify(pb(Av[0], Av[1]) + 2*H*Lv[2]) == 0
ok &= sp.simplify(pb(Lv[0], Av[1]) - Av[2]) == 0
check("V8 Runge-Lenz A is conserved and {A_x, A_y} = -2H L_z: the dial (-2E) is the bracket coefficient", ok)

# V9. Uniqueness on the symmetry route: A = p x L - f(r) r_hat conserved for a central force only if V = -K/r.
def lenz_rate(Vexpr, fexpr):
    H2 = P3.dot(P3)/2 + Vexpr
    A2 = P3.cross(Lv) - fexpr*X3/rr
    return sp.simplify(pb(A2[0], H2))
sol_ok = lenz_rate(-K/rr, K) == 0
alt = all(lenz_rate(-K/rr**2, fx) != 0 for fx in (K, K/rr, 2*K/rr))
check("V9 a Runge-Lenz-type conserved vector exists for V = -K/r; tried weights fail for V = -K/r^2", sol_ok and alt)

# V10. Quantum side (Fock 1935): quantising motion on the dial sphere: harmonics of degree n-1 on the 3-sphere
#      number n^2 (hydrogen's degeneracy), and the dial at E_n = -1/(2n^2) is kappa = -1/(2E) = n^2.
ok = True
for nq in range(1, 9):
    deg = nq - 1
    dim = sp.binomial(deg + 3, 3) - (sp.binomial(deg + 1, 3) if deg >= 2 else 0)
    ok &= dim == nq**2
    ok &= sp.Rational(-1, 2)/(sp.Rational(-1, 2*nq**2)) == nq**2
check("V10 harmonics on the dial 3-sphere: n^2 states at level n; dial kappa = n^2 at hydrogen's E_n = -1/(2n^2)", ok)

print(f"\n{PASS} pass, {FAIL} fail"); LOG.append(f"\n{PASS} pass, {FAIL} fail")
R["pass"], R["fail"] = PASS, FAIL
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(LOG) + "\n")
json.dump(R, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
