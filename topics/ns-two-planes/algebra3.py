"""DRAFT. Round 2: three mutually perpendicular planes (Brennan: "it might be three orthagonal planes").
Plane m = the coordinate plane with normal e_m (m = x, y, z): plane x = yz-plane, plane y = xz-plane, plane z = xy-plane.
Run: python3 algebra3.py -> algebra3_results.json"""
import json, itertools, numpy as np, sympy as sp
res = {}
def check(name, ok, note=''):
    res[name] = {'pass': bool(ok), 'note': note}; print(('PASS ' if ok else 'FAIL ') + name, note)

# ---------- T1: the symmetric three-plane split (numeric, random field)
N = 24; rng = np.random.default_rng(7)
k1 = np.fft.fftfreq(N, 1/N); KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij'); K = np.array([KX, KY, KZ])
K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
U = (rng.normal(size=(3, N, N, N)) + 1j*rng.normal(size=(3, N, N, N)))*np.exp(-K2/25)*(K2 > 0)
U = np.fft.fftn(np.real(np.fft.ifftn(U, axes=(1, 2, 3))), axes=(1, 2, 3)); U = U - K*(K*U).sum(0)/K2s
curl = lambda V: 1j*np.cross(K, V, axis=0)
W = curl(U)
E = np.eye(3)
parts = [curl(np.array([E[m][i]*W[m]/K2s for i in range(3)])) for m in range(3)]   # Biot-Savart of one spin component
err = np.abs(sum(parts) - U).max()/np.abs(U).max()
inplane = max(np.abs(parts[m][m]).max() for m in range(3))
# spin of plane m's flow = w_m e_m + grad(d_m psi_m)  (curl curl = -Laplacian + grad div); the gradients cancel in the sum
spins = [curl(parts[m]) for m in range(3)]
psi = [W[m]/K2s for m in range(3)]
gradpart = [np.array([-K[i]*K[m]*psi[m] for i in range(3)]) for m in range(3)]   # grad(d_m psi_m) in Fourier
own = max(np.abs(spins[m] - np.array([E[m][i]*W[m] for i in range(3)]) - gradpart[m]).max() for m in range(3))/np.abs(W).max()
sumerr = np.abs(sum(spins) - W).max()/np.abs(W).max()
check('T1 any incompressible flow = three layered flat flows, one per perpendicular plane (no plane singled out)',
      err < 1e-12 and inplane < 1e-12,
      f'reconstruction error {err:.1e}; each piece has zero velocity along its plane normal ({inplane:.0e}). '
      'It is the Biot-Savart law applied to one spin component at a time (Coulomb gauge vector potential)')
check('T2 the three planes\' spins add up to the true spin; plane m carries the m-th spin component w_m plus a gradient that cancels across planes',
      own < 1e-12 and sumerr < 1e-12,
      f'errors {own:.1e}, {sumerr:.1e}. For strict flat flows (constant along the normal) the gradient is zero, so plane m\'s spin is exactly w_m')

# ---------- T3: a plane never stretches its own spin; production splits over pairs of planes = crossing lines
u = [np.real(np.fft.ifftn(p, axes=(1, 2, 3))) for p in parts]
w = np.real(np.fft.ifftn(W, axes=(1, 2, 3)))
grad = [np.array([[np.real(np.fft.ifftn(1j*K[j]*parts[m][i])) for j in range(3)] for i in range(3)]) for m in range(3)]
T = np.zeros((3, 3, 3))
for i, j, m in itertools.product(range(3), repeat=3):
    T[i, j, m] = np.mean(w[i]*w[j]*grad[m][i, j])       # spin of plane i, along spin of plane j, gradient from plane m
Ptot = np.mean(np.einsum('ixyz,jxyz,ijxyz->xyz', w, w, sum(grad)))
check('T3a production = sum of the 27 plane terms', abs(T.sum() - Ptot) < 1e-10*max(1, abs(Ptot)))
check('T3b a plane\'s own flow never stretches its own spin (all terms with i = m vanish pointwise)',
      max(np.abs(grad[m][m]).max() for m in range(3)) < 1e-12,
      'plane m\'s flow has no velocity along e_m, so it cannot stretch w_m; every stretching term pairs two different planes')
lines = {2: T[0, :, 1].sum() + T[1, :, 0].sum(), 0: T[1, :, 2].sum() + T[2, :, 1].sum(), 1: T[0, :, 2].sum() + T[2, :, 0].sum()}
check('T3c so all stretching is carried by the three crossing lines (planes i and m meet on the third axis)',
      abs(sum(lines.values()) - Ptot) < 1e-10*max(1, abs(Ptot)),
      'identity checked on a random field (whose net stretching is near 0); the shares in real runs are in threeplanes.py')

# ---------- T4: no flat twist turns the three planes into each other on the cube
# R = cyclic permutation x->y->z->x (120 deg about the diagonal n). A screw (R, t) generates a free action with the
# lattice 2pi Z^3 only if  t.n not in (2pi/sqrt3) Z  (else (R, t+L) is a pure rotation: fixed point), and g^3 = (I, 3(t.n) n)
# must be a lattice vector, which forces t.n in (2pi/sqrt3) Z. Contradiction. Brute force confirms on a grid of t.
R = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
n = np.ones(3)/np.sqrt(3)
def free_and_closes(t):
    t = np.array(t, float)
    g3 = t + R@t + R@R@t
    closes = np.allclose(g3/(2*np.pi), np.round(g3/(2*np.pi)), atol=1e-9)
    # element (R^k, s): fixed point iff the axial part of s is 0 mod the lattice's axial steps 2pi/sqrt3
    s1 = t; s2 = t + R@t
    ok1 = abs(((s1@n)/(2*np.pi/np.sqrt(3)) + 0.5) % 1 - 0.5) > 1e-9
    ok2 = abs(((s2@n)/(2*np.pi/np.sqrt(3)) + 0.5) % 1 - 0.5) > 1e-9
    return closes and ok1 and ok2
grid = [2*np.pi*np.array(v)/12 for v in itertools.product(range(12), repeat=3)]
nfree = sum(free_and_closes(t) for t in grid)
check('T4a no screw that cycles the three planes (x->y->z) gives a flat space on the cube: 0 of 1728 shifts work',
      nfree == 0, 'the cube lattice has no room for it; three-fold twists need a hexagonal box (the tricosm)')
# quarter turn about y (round 1) ties planes x and z and keeps plane y
Rq = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
normals = {m: E[m] for m in range(3)}
img = {m: int(np.argmax(np.abs(Rq@E[m]))) for m in range(3)}
check('T4b the quarter turn swaps planes x and z and maps plane y to itself', img == {0: 2, 1: 1, 2: 0})
# Hantzsche-Wendt (didicosm): three screw half-turns about perpendicular axes; each keeps every coordinate plane
Rh = [np.diag(d) for d in ([1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
keeps = all(int(np.argmax(np.abs(r@E[m]))) == m for r in Rh for m in range(3))
check('T4c the half-turn space (Hantzsche-Wendt, three perpendicular screw half-turns) keeps each of the three planes',
      keeps, 'it fits the three-plane frame but never swaps planes; spin of each plane flips sign under the other two turns')
json.dump(res, open('algebra3_results.json', 'w'), indent=1)
print(f"{sum(r['pass'] for r in res.values())} of {len(res)} pass")
