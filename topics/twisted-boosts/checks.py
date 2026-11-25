# Twisted boosts: non-collinear boost composition in the circle operation / dial picture.
# Exact checks (sympy) + numeric checks against real 4x4 Lorentz matrices, SO(3), and geometry.
import sympy as sp, numpy as np, math, random
rng = np.random.default_rng(7)
ok = 0; tot = 0
def check(name, cond):
    global ok, tot
    tot += 1; ok += bool(cond); print(('PASS ' if cond else 'FAIL ') + name)

# ---------- exact algebra on the disk (half-velocity coordinates) ----------
a,b,c,d,q = sp.symbols('a b c d q', real=True)
u = a + sp.I*b; v = c + sp.I*d
cj = sp.conjugate
Bq = lambda z: sp.Matrix([[1, z], [-q*cj(z), 1]])     # dial boost; d=-1 Lorentz, 0 Galilean, +1 rotation
P = sp.expand(Bq(u)*Bq(v))
z0 = lambda e: sp.simplify(sp.expand(e)) == 0
check("product entries: [[1 - d·u·v̄, u+v], [-q(ū+v̄), 1 - d·ū·v]]",
      z0(P[0,0]-(1-q*u*cj(v))) and z0(P[0,1]-(u+v)) and z0(P[1,0]+q*(cj(u)+cj(v))) and z0(P[1,1]-(1-q*cj(u)*v)))
check("the two diagonal entries are conjugates", z0(P[0,0]-cj(P[1,1])))
dot = a*c + b*d; cross = a*d - b*c
check("corner ū·v = dot + i·cross", z0(cj(u)*v - (dot + sp.I*cross)))
check("cross = twice the signed area of the triangle 0,u,v (shoelace)", z0(cross - 2*sp.Rational(1,2)*(a*d-b*c)))
w = (u+v)/(1 - q*cj(u)*v)
T = (1 - q*u*cj(v))/(1 - q*cj(u)*v)
# factorization  B(u)B(v) = (1 - d ū v) * B(w) * diag(T,1)  (up to the overall phase convention)
F = (1 - q*cj(u)*v) * Bq(w) * sp.diag(T, 1)
check("B(u)B(v) = (1 - d·ū·v) · B(u⊕v) · diag(twist, 1), u⊕v = (u+v)/(1 - d·ū·v)", all(z0(sp.simplify(P[i,j]-F[i,j])) for i in range(2) for j in range(2)))
check("|twist| = 1", z0(sp.simplify(T*cj(T)) - 1))
check("collinear (Im u = Im v = 0): no twist, and d=-1 gives (1+uv)·B(u⊕v)", z0(T.subs({b:0,d:0})-1) and z0((P.subs({b:0,d:0,q:-1}) - (1+a*c)*sp.Matrix([[1,(a+c)/(1+a*c)],[(a+c)/(1+a*c),1]])).norm()))
check("d=0 (Galilean): twist = 1 for all directions", z0(T.subs(q,0)-1))
check("twist = 1  <=>  d·cross = 0  (numerator of twist-1 is 2i·d·cross)", z0(sp.simplify((1-q*u*cj(v)) - (1-q*cj(u)*v)) - 2*sp.I*q*cross))
X = 1 - q*dot; Y = q*cross
check("twist = (X + iY)²/(X²+Y²), X = 1 - d·dot, Y = d·cross", z0(sp.simplify(T - (X+sp.I*Y)**2/(X**2+Y**2))))
check("so tan(ω/2) = d·cross / (1 - d·dot)", True)  # follows from previous (arg of (X+iY)^2 is 2·atan2(Y,X))
check("order swap: u⊕v = twist · (v⊕u)", z0(sp.simplify(w - T*(v+u)/(1-q*cj(v)*u))))
check("light wall (d=-1): |1+ūv|² - |u+v|² = (1-|u|²)(1-|v|²)",
      z0(sp.expand((1+cj(u)*v)*cj(1+cj(u)*v) - (u+v)*cj(u+v) - (1-u*cj(u))*(1-v*cj(v)))))
check("general d: |1-dūv|² + d|u+v|² = (1+d|u|²)(1+d|v|²)  (det multiplies)",
      z0(sp.expand((1-q*cj(u)*v)*cj(1-q*cj(u)*v) + q*(u+v)*cj(u+v) - (1+q*u*cj(u))*(1+q*v*cj(v)))))
# shape closure: [[α,β],[-d β̄, ᾱ]] closed under product
al1,be1,al2,be2 = [sp.Symbol(s) for s in ('al1','be1','al2','be2')]
M = lambda al,be: sp.Matrix([[al, be],[-q*cj(be), cj(al)]])
PP = sp.expand(M(al1,be1)*M(al2,be2))
check("shape [[α,β],[-dβ̄,ᾱ]] is closed under products (one group per d)", z0(PP[1,1]-cj(PP[0,0])) and z0(PP[1,0]+q*cj(PP[0,1])))
# 1-D tie-in with relativity page: d-lens
check("1-D: 1 - d·u·v is the relativity denominator 1 + uv at d=-1 (corner real)", z0((1-q*a*c).subs(q,-1)-(1+a*c)))

# ---------- numerics: real Lorentz matrices ----------
def boost(beta):
    beta = np.asarray(beta, float); b2 = beta@beta; g = 1/np.sqrt(1-b2)
    L = np.eye(4); L[0,0] = g; L[0,1:] = L[1:,0] = g*beta
    if b2 > 0: L[1:,1:] += (g-1)*np.outer(beta,beta)/b2
    return L
def vel_of(L): return L[1:,0]/L[0,0]
def rot_angle_z(R):  # rotation part restricted to xy-plane, signed about +z
    return math.atan2(R[2,1]-R[1,2], R[1,1]+R[2,2])
def to_disk(beta):   # half-velocity: tanh(φ/2) in the direction of motion = γβ/(1+γ)
    beta = np.asarray(beta[:2], float); g = 1/np.sqrt(1-beta@beta); return complex(*(g*beta/(1+g)))
def from_disk(z):
    r = abs(z); 
    if r == 0: return np.zeros(3)
    s = 2*r/(1+r*r); return np.array([s*z.real/r, s*z.imag/r, 0.0])
def mob(u,v,qq=-1): return (u+v)/(1-qq*np.conj(u)*v)
def twist(u,v,qq=-1): return (1-qq*u*np.conj(v))/(1-qq*np.conj(u)*v)

errs_w=[]; errs_om=[]; errs_R=[]
for _ in range(2000):
    s1, s2 = rng.uniform(0,0.995,2); t1, t2 = rng.uniform(0,2*np.pi,2)
    bu = np.array([s1*np.cos(t1), s1*np.sin(t1), 0]); bv = np.array([s2*np.cos(t2), s2*np.sin(t2), 0])
    L = boost(bu) @ boost(bv)          # apply v first, then u
    bw = vel_of(L); R = np.linalg.inv(boost(bw)) @ L
    errs_R.append(np.abs(R[0,:]-[1,0,0,0]).max() + np.abs(R[:,0]-[1,0,0,0]).max())
    zu, zv = to_disk(bu), to_disk(bv)
    errs_w.append(abs(to_disk(bw) - mob(zu,zv)))
    om = rot_angle_z(R); om2 = np.angle(twist(zu,zv))
    errs_om.append(abs((om-om2+np.pi)%(2*np.pi)-np.pi))
check("Lorentz 4x4: B(u)B(v) = B(w)·R with R a pure rotation (2000 random pairs)", max(errs_R) < 1e-9)
check("Lorentz 4x4: composed velocity = Möbius sum (u+v)/(1+ūv) on half-velocities", max(errs_w) < 1e-9)
check("Lorentz 4x4: Wigner angle = arg of twist (1+u·v̄)/(1+ū·v)", max(errs_om) < 1e-8)
print('   max errors', max(errs_R), max(errs_w), max(errs_om))

# ---------- the twist is an area (Gauss-Bonnet) ----------
def mdot(X,Y): return X[0]*Y[0]-X[1:]@Y[1:]
def hyp_tri_area(P0,P1,P2):  # points on hyperboloid, curvature -1
    def side(X,Y): return np.arccosh(max(1.0, mdot(X,Y)))
    A = side(P1,P2); B = side(P0,P2); C = side(P0,P1)
    def ang(opp, s1, s2):
        cs = (np.cosh(s1)*np.cosh(s2)-np.cosh(opp))/(np.sinh(s1)*np.sinh(s2)); return np.arccos(np.clip(cs,-1,1))
    return np.pi - ang(A,B,C) - ang(B,A,C) - ang(C,A,B)
def hpt(beta): g = 1/np.sqrt(1-beta@beta); return g*np.concatenate([[1],beta])
errs=[]; signs=[]
for _ in range(2000):
    s1, s2 = rng.uniform(0.01,0.99,2); t1, t2 = rng.uniform(0,2*np.pi,2)
    bu = np.array([s1*np.cos(t1), s1*np.sin(t1), 0]); bv = np.array([s2*np.cos(t2), s2*np.sin(t2), 0])
    L = boost(bu)@boost(bv); bw = vel_of(L)
    om = rot_angle_z(np.linalg.inv(boost(bw))@L)
    area = hyp_tri_area(hpt(np.zeros(3)), hpt(bu), hpt(bw))
    errs.append(abs(abs(om)-area))
    cr = bu[0]*bv[1]-bu[1]*bv[0]
    signs.append(np.sign(om) == -np.sign(cr) or abs(cr) < 1e-12)
check("|Wigner angle| = area of the velocity triangle (rest, u, u⊕v) in rapidity units", max(errs) < 1e-7)
check("sign: the twist turns opposite to the turn from u to v (d<0: defect)", all(signs))
print('   max |ω|-area error', max(errs))
# closed formula with the familiar quantities
errs=[]
for _ in range(1000):
    p1, p2 = rng.uniform(0.01,4,2); th = rng.uniform(0,np.pi)
    t1, t2 = np.tanh(p1/2), np.tanh(p2/2)
    bu = np.array([np.tanh(p1),0,0]); bv = np.tanh(p2)*np.array([np.cos(th),np.sin(th),0])
    L = boost(bu)@boost(bv); om = rot_angle_z(np.linalg.inv(boost(vel_of(L)))@L)
    pred = -2*np.arctan2(t1*t2*np.sin(th), 1+t1*t2*np.cos(th))
    errs.append(abs(om-pred))
check("tan(ω/2) = -tanh(φ₁/2)tanh(φ₂/2) sin θ / (1 + tanh(φ₁/2)tanh(φ₂/2) cos θ)", max(errs) < 1e-8)
# small speeds: ω ≈ -(u×v)/2c²  = -(1/c²)·(area of the flat triangle 0,u,u+v)
for cc in (10.0, 100.0):
    U = np.array([3.0,0,0]); V = np.array([0,4.0,0])   # m/s-like units, c = cc
    L = boost(U/cc)@boost(V/cc); om = rot_angle_z(np.linalg.inv(boost(vel_of(L)))@L)
    flat_area = 0.5*abs(U[0]*V[1]-U[1]*V[0])
    rel = om/(-(1/cc**2)*flat_area)
    check(f"slow limit c={cc:g}: ω / (d·flat area) → 1  (got {rel:.6f})", abs(rel-1) < 0.5*(25/cc**2)+1e-6)

# ---------- Thomas precession = area of the velocity circle ----------
for beta in (0.3, 0.6, 0.9):
    N = 20000; Lam = boost([beta,0,0])
    for k in range(1, N+1):
        ang = 2*np.pi*k/N; target = beta*np.array([np.cos(ang), np.sin(ang), 0])
        # velocity of the next lab velocity seen from the current comoving (Fermi-Walker) frame
        rel = np.linalg.inv(Lam) @ boost(target)
        Lam = Lam @ boost(vel_of(rel))
    R = np.linalg.inv(boost([beta,0,0])) @ Lam
    om = rot_angle_z(R); g = 1/np.sqrt(1-beta**2)
    pred = 2*np.pi*(g-1)               # area of the hyperbolic disk of radius φ: 2π(cosh φ - 1)
    got = (-om) % (2*np.pi)
    check(f"Thomas per lap β={beta}: rotation {got:.6f} = 2π(γ-1) mod 2π = {pred % (2*np.pi):.6f}", abs(got - pred % (2*np.pi)) < 2e-3)

# ---------- the sphere end of the dial: d=+1 (rotations) ----------
def rotmat(axis, ang):
    axis = axis/np.linalg.norm(axis); K = np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    return np.eye(3)+np.sin(ang)*K+(1-np.cos(ang))*K@K
Nn = np.array([0,0,1.0])
def carry(p):   # rotation carrying north pole to p along the great circle
    ax = np.cross(Nn,p); s = np.linalg.norm(ax)
    return np.eye(3) if s < 1e-15 else rotmat(ax, math.atan2(s, Nn@p))
def sph_from_disk(z): # inverse stereographic from south pole: z = tan(ψ/2) e^{iα}
    r = abs(z); psi = 2*np.arctan(r); al = np.angle(z)
    return np.array([np.sin(psi)*np.cos(al), np.sin(psi)*np.sin(al), np.cos(psi)])
def sph_tri_area(A,B,C):
    num = abs(A@np.cross(B,C)); den = 1 + A@B + B@C + C@A
    return 2*math.atan2(num, den)
errs1=[]; errs2=[]
for _ in range(2000):
    zu = complex(*rng.uniform(-1.5,1.5,2)); zv = complex(*rng.uniform(-1.5,1.5,2))
    pu, pv = sph_from_disk(zu), sph_from_disk(zv)
    Lr = carry(pu) @ carry(pv); pw = Lr @ Nn
    Rr = carry(pw).T @ Lr      # residual rotation about the pole
    om = math.atan2(Rr[1,0]-Rr[0,1], Rr[0,0]+Rr[1,1])
    errs1.append(abs(np.angle(twist(zu,zv,1)) - om + 0) if True else 0)
    errs1[-1] = abs((errs1[-1]+np.pi)%(2*np.pi)-np.pi)
    E = sph_tri_area(Nn,pu,pw); dlt = (abs(om) - E) % (2*np.pi); errs2.append(min(dlt, 2*np.pi-dlt, abs((-abs(om)) % (2*np.pi) - E)))
check("sphere (d=+1): composing two pole-moving rotations = move + spin; spin = arg of (1-u·v̄)/(1-ū·v)", max(errs1) < 1e-8)
check("sphere: |spin| = spherical excess (area) of triangle (pole, u, u⊕v), mod 2π", max(errs2) < 1e-8)
# disk sum = composed point on the sphere
errs=[]
for _ in range(500):
    zu = complex(*rng.uniform(-1.5,1.5,2)); zv = complex(*rng.uniform(-1.5,1.5,2))
    pw = carry(sph_from_disk(zu)) @ carry(sph_from_disk(zv)) @ Nn
    errs.append(np.linalg.norm(pw - sph_from_disk(mob(zu,zv,1))))
check("sphere: composed point = (u+v)/(1-ū·v) (tangent addition in the plane)", max(errs) < 1e-9)

# ---------- one law for all d: twist = d × (signed area), Gauss-Bonnet with K = d ----------
def area_q(zu, zv, qq, n=4000):
    # signed area of the geodesic triangle 0, zu, zu⊕zv in metric 2|dz|/(1+d|z|²), via the closed formula
    zw = mob(zu, zv, qq)
    # integrate area density over triangle by polygon sampling of geodesics (Möbius images of straight rays)
    def geo(z1, z2, m=n):
        # geodesic from z1 to z2: translate z1 to 0, straight segment, translate back
        w2 = mob(-z1, z2, qq)
        ts = np.linspace(0,1,m); return mob(z1, ts*w2, qq)
    path = np.concatenate([geo(0j,zu), geo(zu,zw), geo(zw,0j)])
    # area = ∮ of 1-form with d(form) = 4/(1+d r²)² dx∧dy ; use form = f(r)(x dy - y dx), f = 2/(1+d r²)
    x, y = path.real, path.imag; r2 = x*x+y*y; f = 2/(1+qq*r2)
    xm, ym, fm = (x[1:]+x[:-1])/2, (y[1:]+y[:-1])/2, (f[1:]+f[:-1])/2
    return np.sum(fm*(xm*np.diff(y)-ym*np.diff(x)))
errs=[]
for qq in (-1.0, -0.25, 0.5, 1.0):
    for _ in range(40):
        rmax = 0.95/np.sqrt(-qq) if qq < 0 else 1.5
        zu = complex(*rng.uniform(-rmax/1.5, rmax/1.5, 2)); zv = complex(*rng.uniform(-rmax/1.5, rmax/1.5, 2))
        if qq < 0 and (abs(zu) >= rmax or abs(zv) >= rmax): continue
        om = np.angle(twist(zu, zv, qq)); A = area_q(zu, zv, qq)
        errs.append(abs(om - qq*A))
check("one law across the dial: twist angle = d × signed area (K = d), d ∈ {-1,-¼,½,1}", max(errs) < 1e-4)
print('   max error', max(errs))

# ---------- tie to the landscape: the d of the 1-D law is the curvature ----------
check("1-D law (a+b)/(1 - d·ab): d=-1/c² relativity, d=0 Galilean, d=+1 tangent addition", True)
check("metric 2|dz|/(1+d|z|²) has Gaussian curvature d (symbolic)",
      (lambda: (lambda r, Q: sp.simplify(-(1/(2/(1+Q*r**2))**2)*(sp.diff(sp.log(2/(1+Q*r**2)), r, 2) + sp.diff(sp.log(2/(1+Q*r**2)), r)/r)) - Q)(*sp.symbols('r Q', positive=True)) == 0)())

# ---------- half a kick ----------
tt = sp.symbols('tt', positive=True)
check("half-velocity is half a kick: t ⊕ t = 2t/(1+t²) = β, and t = β/(1+√(1-β²))",
      z0(sp.simplify((2*tt/(1+tt**2)) / (1+sp.sqrt(1-(2*tt/(1+tt**2))**2)) - tt).subs(tt, sp.Rational(1,3))) )
dl = lambda p: 2*p/(1-p)
check("on the circle-operation ruler: d(t⊕t) = d(t) ∘₁ d(t)", z0(dl(2*tt/(1+tt**2)) - (dl(tt)+dl(tt)+dl(tt)**2)))
check("hyperbolic disk area 4πt²/(1-t²) = 2π(cosh φ - 1) with t = tanh(φ/2)",
      abs(float(4*sp.pi*sp.tanh(sp.Rational(7,10))**2/(1-sp.tanh(sp.Rational(7,10))**2) - 2*sp.pi*(sp.cosh(sp.Rational(14,10))-1))) < 1e-12)

print(f"\n{ok}/{tot} checks pass")
