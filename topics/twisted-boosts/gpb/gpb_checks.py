# Gravity Probe B through "twist = curvature x area".
# Geodetic precession = velocity-space twist (Thomas-type, q = -1/c^2) + space-curvature holonomy.
# Checks: exact Schwarzschild parallel transport (numeric ODE), Flamm-surface holonomy and curvature integral,
# the circle-operation combination of the two factors, GP-B numbers, frame dragging as a gravitomagnetic flux.
import numpy as np, sympy as sp, math
ok=0; tot=0
def check(name, cond):
    global ok, tot
    tot+=1; ok+=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)

# ---------- 1. exact Schwarzschild: spin carried on a circular geodesic ----------
# units G=c=M=1; equatorial plane; coordinates (t, r, phi); spin S^mu with S.u = 0
def precession_per_orbit(r, steps=200000):
    # Christoffels for ds^2 = -(1-2/r)dt^2 + dr^2/(1-2/r) + r^2 dphi^2 (theta = pi/2)
    f = 1-2/r
    Om = r**-1.5                       # dphi/dt on circular geodesic
    ut = 1/math.sqrt(1-3/r); uphi = Om*ut
    G = {}  # Gamma^a_bc for a,b,c in t,r,phi (0,1,2), only nonzero
    G[(0,0,1)] = G[(0,1,0)] = 1/(r*r*f)
    G[(1,0,0)] = f/(r*r)
    G[(1,2,2)] = -r*f
    G[(2,1,2)] = G[(2,2,1)] = 1/r
    u = np.array([ut,0.0,uphi])
    def rhs(S):
        dS = np.zeros(3)
        for (a,b,c),g in G.items(): dS[a] -= g*u[b]*S[c]
        return dS
    # initial spin: radial unit vector (orthogonal to u)
    S = np.array([0.0, math.sqrt(f), 0.0])
    tau_orbit = 2*math.pi/uphi; h = tau_orbit/steps
    # comoving orthonormal spatial frame: e_rhat = (0, sqrt f, 0); e_phi' = boosted phi direction
    v = r*Om/math.sqrt(f)*math.sqrt(1)   # local speed relative to static observers
    v = math.sqrt(1/(r-2))*1.0
    g = 1/math.sqrt(1-v*v)
    e_r = np.array([0, math.sqrt(f), 0])
    e_p = g*np.array([v/math.sqrt(f), 0, 1/r])      # boosted: gamma (v e_t_hat + e_phi_hat)
    gmet = np.diag([-f, 1/f, r*r])
    ang0 = math.atan2(e_p@gmet@S, e_r@gmet@S); prev = ang0; tot_ang = 0.0
    for _ in range(steps):
        k1=rhs(S); k2=rhs(S+h/2*k1); k3=rhs(S+h/2*k2); k4=rhs(S+h*k3)
        S = S + h/6*(k1+2*k2+2*k3+k4)
        a = math.atan2(e_p@gmet@S, e_r@gmet@S); d = (a-prev+math.pi)%(2*math.pi)-math.pi; tot_ang += d; prev = a
    return tot_ang   # rotation relative to the radial direction over one orbit
for r in (6.5, 10.0, 30.0, 200.0):
    rel = precession_per_orbit(r, 40000)
    exact = -2*math.pi*math.sqrt(1-3/r)
    check(f"Schwarzschild r={r}M: spin turns {rel:.6f} vs radial; known -2π√(1-3M/r) = {exact:.6f}", abs(rel-exact) < 1e-6)

# ---------- 2. the space part: holonomy of the curved space slice (Flamm's surface) ----------
# 2D metric dr^2/(1-2M/r) + r^2 dphi^2; Gaussian curvature K = -M/r^3
rr, M, ph = sp.symbols('r M phi', positive=True)
a2 = 1/(1-2*M/rr); b2 = rr**2
K = sp.simplify(-1/(2*sp.sqrt(a2*b2))*sp.diff(sp.diff(b2,rr)/sp.sqrt(a2*b2),rr))
check("Gaussian curvature of the space slice is -M/r³ (tidal, saddle-like)", sp.simplify(K + M/rr**3) == 0)
outside = sp.integrate((M/rr**3)*2*sp.pi*rr/sp.sqrt(1-2*M/rr), (rr, sp.Symbol('R',positive=True), sp.oo))
Rs = sp.Symbol('R', positive=True)
check("curvature outside the orbit: ∫ M/r³ dA = 2π(1 - √(1-2M/R)) exactly",
      abs(float((outside - 2*sp.pi*(1-sp.sqrt(1-2*M/Rs))).subs({M:1,Rs:sp.Rational(37,3)}))) < 1e-12)
def flamm_holonomy(r, n=200000):
    # parallel transport a tangent vector once round the circle r on the 2D slice
    f = 1-2/r; X = np.array([math.sqrt(f), 0.0]); h = 2*math.pi/n
    # Christoffels: G^r_phiphi = -r f, G^phi_rphi = 1/r ; velocity dphi = 1
    def rhs(X): return np.array([r*f*X[1], -X[0]/r])
    for _ in range(n):
        k1=rhs(X);k2=rhs(X+h/2*k1);k3=rhs(X+h/2*k2);k4=rhs(X+h*k3);X=X+h/6*(k1+2*k2+2*k3+k4)
    ang = math.atan2(r*X[1], X[0]/math.sqrt(f))
    return ang
for r in (4.0, 10.0, 50.0):
    got = flamm_holonomy(r, 20000)
    want = (2*math.pi*(1-math.sqrt(1-2/r)))   # turn relative to the fixed frame after a lap, as an angle mod 2π
    # vector vs radial ends at -2π√(1-2/r); report mod 2π
    check(f"space slice r={r}M: carried vector comes back turned by 2π(1-√(1-2M/r)) = {want:.6f}", abs(((got - (-2*math.pi*math.sqrt(1-2/r)))+math.pi)%(2*math.pi)-math.pi) < 1e-8)

# ---------- 3. the two factors fuse by the circle operation ----------
x = sp.symbols('x', positive=True)            # x = M/r
s = sp.sqrt(1-2*x)                           # space factor (Flamm)
v2 = x/(1-2*x)                               # local orbital speed² vs static observers
gam = 1/sp.sqrt(1-v2)
check("local orbital speed: v² = M/(r-2M), so γ = √((r-2M)/(r-3M))", sp.simplify(gam**2 - (1-2*x)/(1-3*x)) == 0)
check("exact factorization: √(1-3M/r) = √(1-2M/r) · (1/γ)  (squared exactly; and numerically on x in (0,1/3))", sp.simplify((1-3*x) - (1-2*x)*(1-v2)) == 0 and max(abs(float((sp.sqrt(1-3*x) - s/gam).subs(x,xx))) for xx in [0.001,0.05,0.1,0.2,0.3,0.33]) < 1e-14)
ps = 1 - s; pv = 1 - 1/gam; ptot = 1 - sp.sqrt(1-3*x)
check("precession fractions fuse: -p = (-p_space) ∘₁ (-p_velocity)  (since 1-p = (1-p_s)(1-p_v))", max(abs(float((-ptot - ((-ps)+(-pv)+(-ps)*(-pv))).subs(x,xx))) for xx in [0.001,0.05,0.1,0.2,0.3,0.33]) < 1e-14)
check("first order: p_space ≈ M/r, p_velocity ≈ M/(2r), total ≈ 3M/(2r)",
      sp.series(ps,x,0,2).removeO()==x and sp.series(pv,x,0,2).removeO()==x/2 and sp.series(ptot,x,0,2).removeO()==sp.Rational(3,2)*x)
check("velocity part to first order = Thomas: (γ-1) and (1-1/γ) both ≈ v²/2 = M/(2r)",
      sp.series(gam-1,x,0,2).removeO()==x/2)

# ---------- 4. Gravity Probe B numbers ----------
GM = 3.986004418e14; c = 299792458.0; Re = 6378.137e3
yr = 365.25*86400; mas = 180/math.pi*3600*1000
r = 7027.4e3   # GP-B semi-major axis (642 km mean altitude)
n_orb = math.sqrt(GM/r**3)                      # orbits: angular rate
space_rate = (GM/(c*c*r))*n_orb                 # 2πM/r per orbit -> times n/2π
vel_rate = 0.5*(GM/(c*c*r))*n_orb               # π v²/c² per orbit
geo = (space_rate+vel_rate)*yr*mas
exact = (1-math.sqrt(1-3*GM/(c*c*r)))*n_orb*yr*mas
print(f"   r = {r/1e3:.1f} km: space part {space_rate*yr*mas:.1f} mas/yr, velocity part {vel_rate*yr*mas:.1f} mas/yr, total {geo:.1f}, exact {exact:.1f}")
check("space part + velocity part = (3/2)(GM)^{3/2}/(c² r^{5/2})", abs(geo - 1.5*GM**1.5/(c*c*r**2.5)*yr*mas) < 1e-6)
check(f"GP-B geodetic: twist bookkeeping gives {geo:.0f} mas/yr; measured 6601.8 ± 18.3 (GR with Earth's shape: 6606.1)", abs(geo-6601.8) < 2*18.3)
check("without the space part (velocity twist alone) it would be 1/3 of that: ruled out by ~240 sigma", abs(vel_rate*yr*mas - 6601.8)/18.3 > 100)
check("without the velocity twist (space alone) it would be 2/3: ruled out by ~120 sigma", abs(space_rate*yr*mas - 6601.8)/18.3 > 100)

# ---------- 5. frame dragging: twist as a gravitomagnetic flux ----------
I = 8.034e37; w = 7.2921159e-5; G = 6.67430e-11
# polar orbit average of (G I / c² r³)[3 (ω·n) n - ω]: the sin² average gives ω/2
tt = np.linspace(0, 2*np.pi, 100001)[:-1]
avg = np.mean(3*np.sin(tt)**2 - 1)
check("polar-orbit average of 3 sin² - 1 is 1/2", abs(avg-0.5) < 1e-9)
fd = G*I*w/(c*c*r**3)*avg*yr*mas
print(f"   frame dragging {fd:.1f} mas/yr")
check(f"GP-B frame dragging: {fd:.1f} mas/yr vs measured 37.2 ± 7.2 (GR 39.2)", abs(fd-37.2) < 7.2)
# Stokes: local precession = half the circulation of the gravitomagnetic potential per unit area.
# Potential A = (G/c²) (J × r)/|r|³ (J = I ω). Small loop of radius ε centred at a point: circulation / area -> (∇×A)·n
J = I*w*np.array([0,0,1.0])
def A(p): return G/(c*c)*np.cross(J, p)/np.linalg.norm(p)**3
p0 = np.array([0.0, 0.0, r]); eps = 50.0; m = 4000
circ = 0.0
for k in range(m):
    th0, th1 = 2*np.pi*k/m, 2*np.pi*(k+1)/m
    q0 = p0 + eps*np.array([np.cos(th0), np.sin(th0), 0]); q1 = p0 + eps*np.array([np.cos(th1), np.sin(th1), 0])
    circ += A((q0+q1)/2) @ (q1-q0)
curl_z = circ/(np.pi*eps**2)
pred = G/(c*c*r**3)*(3*(J@p0/r)*(p0[2]/r) - J[2])   # (∇×A)_z = (G/c²r³)[3(J·n)n_z - J_z]
check("frame-dragging field = circulation per unit area of the gravitomagnetic potential (Stokes), over the pole", abs(curl_z/pred-1) < 1e-4)
print(f"\n{ok}/{tot} checks pass")
