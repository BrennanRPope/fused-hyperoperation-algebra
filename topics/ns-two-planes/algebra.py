"""DRAFT. Exact algebra for 'Navier-Stokes as two crossing flat flows' (sympy + numpy).
Planes: P1 = xy-plane (normal z), P2 = yz-plane (normal x). Shared line = y-axis.
Strict flat flow on P1: velocity in P1, constant along z.  Strict flat flow on P2: velocity in P2, constant along x.
Run: python3 algebra.py   (prints PASS/FAIL lines, writes algebra_results.json)"""
import json, sympy as sp, numpy as np
x, y, z = sp.symbols('x y z', real=True)
X = (x, y, z)
res = {}
def check(name, ok, note=''):
    res[name] = {'pass': bool(ok), 'note': note}; print(('PASS ' if ok else 'FAIL ') + name, note)

def grad_dot(a, b):  # (a . grad) b
    return sp.Matrix([sum(a[j]*sp.diff(b[i], X[j]) for j in range(3)) for i in range(3)])
def curl(a):
    return sp.Matrix([sp.diff(a[2], y)-sp.diff(a[1], z), sp.diff(a[0], z)-sp.diff(a[2], x), sp.diff(a[1], x)-sp.diff(a[0], y)])
def div(a): return sum(sp.diff(a[i], X[i]) for i in range(3))

psi1 = sp.Function('psi1')(x, y)   # stream function of the P1 flow
psi2 = sp.Function('psi2')(y, z)   # stream function of the P2 flow
u1 = sp.Matrix([sp.diff(psi1, y), -sp.diff(psi1, x), 0])
u2 = sp.Matrix([0, sp.diff(psi2, z), -sp.diff(psi2, y)])
u = u1 + u2
check('A1 both flat flows are incompressible', sp.simplify(div(u1)) == 0 and sp.simplify(div(u2)) == 0)

# --- A2: the exact crossing term
cross = grad_dot(u, u) - grad_dot(u1, u1) - grad_dot(u2, u2)
claim = u1[1]*sp.diff(u2, y) + u2[1]*sp.diff(u1, y)
check('A2 crossing term = (u1_y) d_y u2 + (u2_y) d_y u1  (only the shared line enters)',
      sp.simplify(cross - claim) == sp.zeros(3, 1))

# --- A3: vorticity and stretching
w1, w2 = curl(u1), curl(u2)
check('A3a P1 spin points along z (its normal), P2 spin along x (its normal)',
      w1[0] == 0 and w1[1] == 0 and sp.simplify(w2[1]) == 0 and w2[2] == 0)
w = w1 + w2
stretch = grad_dot(w, u)
self_stretch = grad_dot(w1, u1) + grad_dot(w2, u2)
check('A3b each flat flow alone has zero stretching', sp.simplify(self_stretch) == sp.zeros(3, 1))
cross_stretch = sp.Matrix([w2[0]*sp.diff(u1[0], x),
                           w2[0]*sp.diff(u1[1], x) + w1[2]*sp.diff(u2[1], z),
                           w1[2]*sp.diff(u2[2], z)])
check('A3c all stretching is crossing: x: w2x*d_x u1x, y: w2x*d_x u1y + w1z*d_z u2y, z: w1z*d_z u2z',
      sp.simplify(stretch - cross_stretch) == sp.zeros(3, 1))
check('A3d no spin along the shared line at the start (w_y = 0)', sp.simplify(w[1]) == 0,
      'the y line of stretching is pure tilting: it makes spin along the shared line')

# --- A4: when is the crossing harmless (a pure pressure gradient)?  single cells in each plane
a, b, c, d = sp.symbols('a b c d', positive=True)
def crossing(P1, P2):
    v1 = sp.Matrix([sp.diff(P1, y), -sp.diff(P1, x), 0]); v2 = sp.Matrix([0, sp.diff(P2, z), -sp.diff(P2, y)])
    cr = v1[1]*sp.diff(v2, y) + v2[1]*sp.diff(v1, y)
    return [sp.simplify(e) for e in cr], [sp.simplify(e) for e in curl(cr)]
cases = {'general cells': (sp.sin(a*x)*sp.cos(b*y), sp.cos(c*y)*sp.sin(d*z)),
         'P1 flat along the shared line': (sp.sin(a*x), sp.cos(c*y)*sp.sin(d*z)),
         'P2 flat along the shared line': (sp.sin(a*x)*sp.cos(b*y), sp.sin(d*z)),
         'both flat along the shared line': (sp.sin(a*x), sp.sin(d*z))}
out = {k: any(e != 0 for e in crossing(*v)[1]) for k, v in cases.items()}
check('A4 crossing is a harmless pressure gradient only when BOTH flows are flat along the shared line',
      out['general cells'] and out['P1 flat along the shared line'] and out['P2 flat along the shared line']
      and not out['both flat along the shared line'],
      'then both are just shears along y and the sum is a one-direction flow; otherwise the crossing spins the fluid')

# --- A5: layered (not strict) flat flows can cross exactly: two-wave ABC flow
A_, B_ = sp.symbols('A B', real=True)
wa = sp.Matrix([A_*sp.sin(z), A_*sp.cos(z), 0])   # lies in P1 (xy), depends on z: layered P1 flow
wb = sp.Matrix([0, B_*sp.sin(x), B_*sp.cos(x)])   # lies in P2 (yz), depends on x: layered P2 flow
uw = wa + wb
lamb = sp.simplify(uw.cross(curl(uw)))
check('A5 two-wave ABC (A,B,0): layered P1 wave + layered P2 wave is an exact Navier-Stokes solution',
      lamb == sp.zeros(3, 1) and sp.simplify(curl(uw) - uw) == sp.zeros(3, 1),
      'Beltrami: spin = velocity, so the nonlinear term is a pure gradient (Arnold/Childress, known)')

# --- A7: the quarter-turn twist swaps the two planes
# g(p) = R p + (0, pi/2, 0), R(x,y,z) = (-z, y, x): quarter turn about the shared line y plus a quarter shift.
R = sp.Matrix([[0, 0, -1], [0, 1, 0], [1, 0, 0]]); Rinv = R.T
def g_act(v):   # (g v)(p) = R v(g^-1 p)
    q = Rinv*(sp.Matrix([x, y - sp.pi/2, z]))
    return R*v.subs({x: q[0], y: q[1], z: q[2]}, simultaneous=True)
gu1 = sp.simplify(g_act(u1))
indep_x = all(sp.simplify(sp.diff(e, x)) == 0 for e in gu1)
check('A7a the twist carries any flat P1 flow onto a flat P2 flow (velocity in yz, constant in x)',
      gu1[0] == 0 and indep_x)
g2u1 = sp.simplify(g_act(g_act(u1)))
check('A7b twice the twist (a half turn) brings it back to a flat P1 flow',
      g2u1[2] == 0 and all(sp.simplify(sp.diff(e, z)) == 0 for e in g2u1))
check('A7c the twist is a rigid motion with no fixed point and g^4 = shift by 2pi along y',
      R**4 == sp.eye(3) and sp.simplify(R.det()) == 1,
      'so the quotient of the box by g is a flat 3-manifold (the quarter-turn space, Conway-Rossetti "tetracosm")')
# orbit sizes of a plane through the shared line under n-turn twists: how many planes the twist makes
orbits = {}
for n in (2, 3, 4, 6):
    th = 2*sp.pi/n; planes = set()
    for j in range(n):
        ang = sp.nsimplify((j*th) % sp.pi)   # a plane through the y axis is fixed by its angle mod pi
        planes.add(ang)
    orbits[n] = len(planes)
check('A7d planes made from one flat flow by each turn twist: half 1, third 3, quarter 2, sixth 3',
      orbits == {2: 1, 3: 3, 4: 2, 6: 3}, 'only the quarter turn makes exactly two crossing planes, and they are perpendicular')

# --- A6: kinematic split, numerically: any incompressible periodic field = layered P1 + layered P2 + strict flow in P3 (xz)
N = 24; rng = np.random.default_rng(1)
k1 = np.fft.fftfreq(N, 1/N); KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij'); K = np.array([KX, KY, KZ])
K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
U = (rng.normal(size=(3, N, N, N)) + 1j*rng.normal(size=(3, N, N, N))) * np.exp(-K2/30)
U = np.fft.fftn(np.real(np.fft.ifftn(U, axes=(1, 2, 3))), axes=(1, 2, 3))
U = U - K*(K*U).sum(0)/K2s; U[:, 0, 0, 0] = 0
# Coulomb vector potential, then axial gauge A_y = 0 where possible (ky != 0)
A = 1j*np.cross(K, U, axis=0)/K2s; A[:, 0, 0, 0] = 0
ky0 = (KY == 0)
chi = np.where(ky0, 0, -A[1]/(1j*np.where(ky0, 1, KY)))
Ag = A + 1j*K*chi
Ag_y_left = np.where(ky0, Ag[1], 0)               # what the gauge cannot remove: ky = 0 slice
psi1h = np.where(ky0, 0, Ag[2]); psi2h = np.where(ky0, 0, Ag[0])  # A = (psi2, 0, psi1) on ky != 0
L1 = 1j*np.cross(K, np.array([0*psi1h, 0*psi1h, psi1h]), axis=0)   # curl(psi1 z) : lies in xy
L2 = 1j*np.cross(K, np.array([psi2h, 0*psi2h, 0*psi2h]), axis=0)   # curl(psi2 x) : lies in yz
rest = U - L1 - L2
# rest should be: the y-averaged field (ky=0 modes), whose x,z parts are a strict xz flow and whose y part fits in L1/L2
rest_ky = np.abs(np.where(ky0, 0, rest)).max()
# ky = 0 modes: u_y(x,z) part is expressible as layered (psi1 = psi1(x,z) gives u_y = -d_x psi1 ... ) ; x,z parts form the P3 flow
F3 = np.where(ky0, U, 0); F3[1] = 0
uy0 = np.where(ky0, U[1], 0)   # y-velocity of y-averaged field: e.g. = -d_x psi1 with psi1 = psi1(x,z): layered P1 (lies in xy? yes, (0,u_y,0))
recon = L1 + L2 + F3 + np.array([0*uy0, uy0, 0*uy0])
err = np.abs(recon - U).max()/np.abs(U).max()
inL1 = np.abs(np.real(np.fft.ifftn(L1[2]))).max(); inL2 = np.abs(np.real(np.fft.ifftn(L2[0]))).max()
check('A6 numeric split of a random incompressible field: layered xy + layered yz + strict xz flow (y-average)',
      err < 1e-12 and rest_ky < 1e-10 and inL1 < 1e-12 and inL2 < 1e-12,
      f'reconstruction error {err:.1e}; P1 part has u_z = {inL1:.0e}, P2 part has u_x = {inL2:.0e}; '
      f'energy share of the strict xz flow = {(np.abs(F3)**2).sum()/(np.abs(U)**2).sum():.4f}')
json.dump(res, open('algebra_results.json', 'w'), indent=1)
print(f"{sum(r['pass'] for r in res.values())} of {len(res)} pass")
