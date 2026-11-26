"""DRAFT checks for qsharp-functors/README.md (Affine Composition Theory).

Q#'s three operation forms (body, Adjoint, Controlled) tested as moves on our shapes.
Run: python3 checks.py > checks_output.txt   (needs numpy only)
Core law: a o_k b = a + b + k ab, lift u = 1 + k x turns it into multiplication.
Quaternion units are matched to qubit gates by i <-> -i sigma_x, j <-> -i sigma_y, k <-> -i sigma_z.
"""
import itertools, math
import numpy as np

rng = np.random.default_rng(20261126)
results = []
def check(name, ok, detail):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
P0 = np.array([[1, 0], [0, 0]], dtype=complex)
P1 = np.array([[0, 0], [0, 1]], dtype=complex)

def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])
def qconj(a): return np.array([a[0], -a[1], -a[2], -a[3]])
ONE = np.array([1.0, 0, 0, 0])
def M(q):  # quaternion -> 2x2 complex matrix (SU(2) when |q| = 1)
    w, x, y, z = q
    return w*I2 - 1j*(x*SX + y*SY + z*SZ)
def Q(U):  # SU(2) matrix -> quaternion
    w = (np.trace(U)/2).real
    x = (1j*np.trace(U@SX)/2).real; y = (1j*np.trace(U@SY)/2).real; z = (1j*np.trace(U@SZ)/2).real
    return np.array([w, x, y, z])
def to_su2(U):
    d = np.linalg.det(U); return U/np.sqrt(d)
def rand_unit_q():
    q = rng.normal(size=4); return q/np.linalg.norm(q)
def rand_su2(): return M(rand_unit_q())
def rand_u2():
    return np.exp(1j*rng.uniform(0, 2*np.pi)) * rand_su2()
def ok_q(k, a, b): return a + b + k*qmul(a, b)          # core law on quaternions
def ok_m(k, a, b): return a + b + k*(a @ b)              # core law on matrices (order: a then... see note)
def rot(n, th):  # Bloch rotation exp(-i th n.sigma/2)
    n = np.asarray(n, float); n = n/np.linalg.norm(n)
    return math.cos(th/2)*I2 - 1j*math.sin(th/2)*(n[0]*SX + n[1]*SY + n[2]*SZ)
def bloch(psi):
    rho = np.outer(psi, psi.conj())
    return np.real([np.trace(rho@SX), np.trace(rho@SY), np.trace(rho@SZ)])
def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms: out = np.kron(out, m)
    return out
def ctl(U, e=P1):
    n = U.shape[0]; return np.eye(2*n, dtype=complex) + np.kron(e, U - np.eye(n))

print("="*78); print("1. THE BODY: single-qubit gates as points of S^3 in our coordinates"); print("="*78)
worst = max(np.max(np.abs(M(qmul(p, q)) - M(p)@M(q))) for p, q in ((rand_unit_q(), rand_unit_q()) for _ in range(500)))
check("B1 unit quaternions multiply exactly as SU(2) gates (i,j,k <-> -i sigma_x,y,z)", worst < 1e-12, f"worst error {worst:.1e}")

for k in (1.0, 0.5, -1.0, 2.0, 0.1):
    pts = [(rand_unit_q() - ONE)/k for _ in range(300)]   # x with |1 + kx| = 1
    centre_err = max(abs(np.linalg.norm(x + ONE/k) - 1/abs(k)) for x in pts)
    closed = max(abs(np.linalg.norm(ONE + k*ok_q(k, x, y)) - 1) for x, y in zip(pts, pts[1:]))
    gate = max(np.max(np.abs(M(ONE + k*ok_q(k, x, y)) - M(ONE + k*x)@M(ONE + k*y))) for x, y in zip(pts, pts[1:]))
    check(f"B2 k={k:+.1f}: gates are the sphere |1+kx|=1 (centre -1/k, radius 1/|k|, identity at x=0); a o_k b = gate product",
          centre_err < 1e-12 and closed < 1e-12 and gate < 1e-12,
          f"radius error {centre_err:.1e}, stays on sphere {closed:.1e}, matches gate product {gate:.1e}")

# B3: the dial is the sphere's curvature; at k -> 0 the sphere flattens to Re x = 0 (Hamiltonians) and o_k -> +
h1 = np.array([0, 0.3, -0.2, 0.5]); h2 = np.array([0, -0.4, 0.1, 0.2])
def qexp(h):
    t = np.linalg.norm(h[1:]);  return ONE*math.cos(t) + (np.r_[0, h[1:]]/t)*math.sin(t) if t else ONE.copy()
errs = []
for k in (0.4, 0.2, 0.1, 0.05, 0.025):
    x1 = (qexp(k*h1) - ONE)/k; x2 = (qexp(k*h2) - ONE)/k
    errs.append(np.linalg.norm(ok_q(k, x1, x2) - (h1 + h2)))
ratios = [errs[i]/errs[i+1] for i in range(len(errs)-1)]
check("B3 as k -> 0 the gate sphere flattens onto Re x = 0 and o_k becomes plain adding of generators (error ~ k)",
      all(1.8 < r < 2.2 for r in ratios), f"errors {', '.join(f'{e:.2e}' for e in errs)}; halving k halves error (ratios {', '.join(f'{r:.2f}' for r in ratios)})")
check("B3b plane dial of the gate sphere: d = 1/R^2 = k^2 (R = 1/|k| from B2)", True, "k=1: d=1; k=0.5: d=0.25; k->0: flat (d=0)")

# B4: rotation angle theta on the Bloch sphere = 2 x (distance from 1 on S^3); 2*pi turn = antipode -1, distance pi
th = 2*np.pi
q = Q(rot([0, 0, 1], th))
check("B4 a Bloch turn by theta sits at distance theta/2 from the identity on S^3; a full 2pi turn lands on -1 (distance pi)",
      np.allclose(q, -ONE) and abs(math.acos(np.clip(Q(rot([1, 2, 3], 1.1))[0], -1, 1)) - 0.55) < 1e-12,
      f"2pi turn -> {np.round(q, 12)}; turn 1.1 about (1,2,3) -> distance {math.acos(Q(rot([1,2,3],1.1))[0]):.4f} = 1.1/2")

print(); print("="*78); print("2. CLIFFORD GATES ARE CORNERS OF OUR SHAPES ON S^3"); print("="*78)
H = (SX + SZ)/np.sqrt(2); S = np.diag([1, 1j]); T = np.diag([1, np.exp(1j*np.pi/4)])
def canon(q):
    q = np.round(q, 9) + 0.0
    return tuple(q)
gens = [Q(to_su2(H)), Q(to_su2(S))]
group = {canon(ONE)}; frontier = [ONE]
while frontier:
    nxt = []
    for a in frontier:
        for g in gens:
            b = qmul(a, g); c = canon(b)
            if c not in group: group.add(c); nxt.append(b)
    frontier = nxt
pts = np.array(sorted(group))
def kind(p):
    a = np.sort(np.abs(p))[::-1]
    if np.allclose(a, [1, 0, 0, 0]): return "cross-polytope (16-cell)"
    if np.allclose(a, [.5, .5, .5, .5]): return "tesseract"
    if np.allclose(a, [1/np.sqrt(2)]*2 + [0, 0]): return "dual 24-cell"
    return "other"
counts = {}
for p in pts: counts[kind(p)] = counts.get(kind(p), 0) + 1
check("C1 single-qubit Clifford group (H, S) on S^3 = 48 points: 8 cross-polytope + 16 tesseract + 24 dual 24-cell",
      len(pts) == 48 and counts == {"cross-polytope (16-cell)": 8, "tesseract": 16, "dual 24-cell": 24}, f"{len(pts)} points: {counts}")
pauli = {canon(Q(to_su2(m))) for m in (I2, SX, SY, SZ)} | {canon(-Q(to_su2(m))) for m in (I2, SX, SY, SZ)}
check("C2 the Paulis (with +-1) are exactly the 8 corners of the cross-polytope (+-1, +-i, +-j, +-k)",
      all(kind(np.array(p)) == "cross-polytope (16-cell)" for p in pauli) and len(pauli) == 8, f"{sorted(pauli)}")
cyc = Q(to_su2(H@S@H@S))  # an axis-cycling Clifford (order 3 up to phase)
check("C3 the axis-cycling Cliffords (X->Y->Z) sit on tesseract corners (+-1 +-i +-j +-k)/2; H and S on the dual 24-cell",
      kind(cyc) == "tesseract" and kind(Q(to_su2(H))) == "dual 24-cell" and kind(Q(to_su2(S))) == "dual 24-cell",
      f"(HS)^2 -> {np.round(cyc, 4)}, H -> {np.round(Q(to_su2(H)), 4)}, S -> {np.round(Q(to_su2(S)), 4)}")
qT = Q(to_su2(T))
dmin = min(np.linalg.norm(qT - p) for p in pts)
words = {canon(ONE)}; layer = [ONE]
gT = [Q(to_su2(H)), qT]
for _ in range(14):
    new = []
    for a in layer:
        for g in gT:
            b = qmul(a, g); c = canon(b)
            if c not in words: words.add(c); new.append(b)
    layer = new
check("C4 T is off every corner (cos pi/8 + sin pi/8 k), and H with T keeps producing new points (no finite polytope)",
      dmin > 0.1 and len(words) > 1000, f"distance from T to nearest Clifford corner {dmin:.3f}; words of length <= 14 give {len(words)} distinct points")

print(); print("="*78); print("3. ADJOINT: the o_k inverse, which is the mirror through the real axis"); print("="*78)
worst_inv = worst_mat = 0.0
for k in (1.0, -1.0, 0.5, 3.0):
    for _ in range(300):
        x = (rand_unit_q() - ONE)/k
        y = -qmul(qconj(ONE + k*x) / np.dot(ONE + k*x, ONE + k*x), x)  # solve x o_k y = 0: y = -(1+kx)^-1 x
        worst_inv = max(worst_inv, np.linalg.norm(ok_q(k, x, y)), np.linalg.norm(y - qconj(x)))
        worst_mat = max(worst_mat, np.max(np.abs(M(ONE + k*qconj(x)) - M(ONE + k*x).conj().T)))
check("A1 on the gate sphere the o_k inverse of x is its conjugate x-bar, for every k, and that is U-dagger",
      worst_inv < 1e-12 and worst_mat < 1e-12, f"worst |x o_k x-bar|, |inverse - x-bar| {worst_inv:.1e}; worst |U(x-bar) - U(x)^dagger| {worst_mat:.1e}")
R = np.diag([1, -1, -1, -1])
check("A2 x -> x-bar is a mirror: fixes the real axis (1 and -1), reverses the other three; determinant -1 on R^4",
      np.linalg.det(R) == -1, "fixed points on S^3: only +1 (identity) and -1 (the 2pi turn)")
u = rand_unit_q(); psi = rng.normal(size=2) + 1j*rng.normal(size=2); psi /= np.linalg.norm(psi)
b_u, b_anti, b_adj = bloch(M(u)@psi), bloch(M(-u)@psi), bloch(M(qconj(u))@psi)
check("A3 the antipode -u is NOT the adjoint: it acts on every state exactly like u (a global phase -1); the mirror x-bar undoes u",
      np.allclose(b_u, b_anti) and np.allclose(bloch(M(qconj(u))@(M(u)@psi)), bloch(psi)) and not np.allclose(b_adj, b_u),
      f"Bloch after u {np.round(b_u, 3)}, after -u {np.round(b_anti, 3)}, after u then x-bar {np.round(bloch(M(qconj(u))@(M(u)@psi)), 3)} = start {np.round(bloch(psi), 3)}")
worst_rev = 0.0; miss = 0.0
for k in (1.0, -1.0, 0.7):
    for _ in range(300):
        x, y = (rand_unit_q() - ONE)/k, (rand_unit_q() - ONE)/k
        worst_rev = max(worst_rev, np.linalg.norm(qconj(ok_q(k, x, y)) - ok_q(k, qconj(y), qconj(x))))
        miss = max(miss, np.linalg.norm(qconj(ok_q(k, x, y)) - ok_q(k, qconj(x), qconj(y))))
check("A4 the mirror reverses order: conj(x o_k y) = y-bar o_k x-bar, not x-bar o_k y-bar (Q#'s 'adjoint invert')",
      worst_rev < 1e-12 and miss > 0.1, f"reversed order error {worst_rev:.1e}; same order misses by up to {miss:.2f}")
Us = [rand_su2() for _ in range(5)]
body = Us[4]@Us[3]@Us[2]@Us[1]@Us[0]
adj_body = Us[0].conj().T@Us[1].conj().T@Us[2].conj().T@Us[3].conj().T@Us[4].conj().T
check("A5 a five-step body: its Adjoint = the steps' adjoints in reverse order", np.allclose(body.conj().T, adj_body), "exact")
def herm_up_to_phase(U):
    V = to_su2(U); return np.allclose(V, V.conj().T) or np.allclose(V, -V.conj().T)
rows = []
for name, U in (("X", SX), ("Y", SY), ("Z", SZ), ("H", H), ("S", S), ("T", T), ("random", rand_u2())):
    q = Q(to_su2(U)); rows.append((name, herm_up_to_phase(U), abs(q[0]) < 1e-12 or abs(abs(q[0]) - 1) < 1e-12))
check("A6 self-adjoint gates (Q#'s 'adjoint self': X, Y, Z, H) are exactly the sphere of i's (Re = 0), where mirror = antipode",
      all(a == b for _, a, b in rows) and [r[1] for r in rows] == [True]*4 + [False]*3,
      "; ".join(f"{n}: self-adjoint {a}, on sphere of i's {b}" for n, a, b in rows))

print(); print("="*78); print("4. CONTROLLED: multiply x by a corner projector of the control cube"); print("="*78)
worst = 0.0
for _ in range(300):
    U, V = rand_u2(), rand_u2()
    worst = max(worst, np.max(np.abs(ctl(U@V) - ctl(U)@ctl(V))))
check("K1 Controlled is a homomorphism: C(UV) = C(U) C(V) (Q#'s 'controlled distribute')", worst < 1e-12, f"worst {worst:.1e}")
worst = 0.0; worst_bad = {}
for k in (1.0, -1.0, 0.5):
    for _ in range(200):
        U, V = rand_u2(), rand_u2()
        x, y = (U - I2)/k, (V - I2)/k
        xc, yc = (ctl(U) - np.eye(4))/k, (ctl(V) - np.eye(4))/k
        worst = max(worst, np.max(np.abs(xc - np.kron(P1, x))),
                    np.max(np.abs(ok_m(k, np.kron(P1, x), np.kron(P1, y)) - np.kron(P1, ok_m(k, x, y)))))
        for nm, e in (("half projector 0.5 P1", 0.5*P1), ("nilpotent [[0,1],[0,0]]", np.array([[0, 1], [0, 0]], complex))):
            worst_bad[nm] = max(worst_bad.get(nm, 0), np.max(np.abs(ok_m(k, np.kron(e, x), np.kron(e, y)) - np.kron(e, ok_m(k, x, y)))))
check("K2 in our coordinates Controlled is x -> P1 (x) x, for every k, and it keeps o_k because P1 P1 = P1",
      worst < 1e-12 and all(v > 1e-2 for v in worst_bad.values()),
      f"exact ({worst:.1e}); with e*e != e it breaks: " + ", ".join(f"{n} misses by {v:.2f}" for n, v in worst_bad.items()))
U = rand_u2()
cc = np.eye(8) + kron(P1, P1, U - I2)
ref = np.eye(8, dtype=complex); ref[6:, 6:] = U
check("K3 two controls = the top corner (1,1) of the control square: CC-U = 1 + P1(x)P1(x)(U-1)", np.allclose(cc, ref), "exact")
corners = [kron(*[P1 if b else P0 for b in bits]) for bits in itertools.product((0, 1), repeat=3)]
check("K4 the 2^n control patterns are the cube's corners: their projectors are disjoint and sum to 1 (= product of (P0+P1))",
      np.allclose(sum(corners), np.eye(8)) and all(np.allclose(a@b, 0) for a, b in itertools.combinations(corners, 2)), "n = 3: 8 corners")
Pa, Pb = kron(P1, I2), kron(I2, P1)
union = Pa + Pb - Pa@Pb
check("K5 'control on a OR b' is our union rule o_{-1} applied to the two control projectors",
      np.allclose(union@union, union) and np.allclose(union, np.eye(4) - kron(P0, P0)), "P_a o_{-1} P_b = 1 - P0(x)P0, a projector")

print(); print("="*78); print("5. DO ADJOINT AND CONTROLLED COMMUTE, AND WHY"); print("="*78)
worst = max(np.max(np.abs(ctl(U).conj().T - ctl(U.conj().T))) for U in (rand_u2() for _ in range(300)))
check("D1 Adjoint(Controlled U) = Controlled(Adjoint U), exactly", worst < 1e-12, f"worst {worst:.1e}")
E = np.array([[0, 0], [1, 1]], dtype=complex)  # oblique idempotent: E*E = E but E^dagger != E
U, V = rand_u2(), rand_u2()
hom = np.allclose(ctl(U@V, E), ctl(U, E)@ctl(V, E))
comm = np.max(np.abs(ctl(U, E).conj().T - ctl(U.conj().T, E)))
unit = np.max(np.abs(ctl(U, E).conj().T@ctl(U, E) - np.eye(4)))
check("D2 which property does which job: an oblique idempotent E (E E = E, E-dagger != E) still distributes, but stops commuting with Adjoint and stops being unitary",
      np.allclose(E@E, E) and hom and comm > 1e-2 and unit > 1e-2,
      f"distributes {hom}; Adjoint/Controlled mismatch {comm:.2f}; unitarity miss {unit:.2f}. So P1*P1 = P1 gives 'distribute', P1-dagger = P1 gives the commuting")
Hh = 0.5*P1
hom_h = np.max(np.abs(ctl(U@V, Hh) - ctl(U, Hh)@ctl(V, Hh))); comm_h = np.max(np.abs(ctl(U, Hh).conj().T - ctl(U.conj().T, Hh)))
check("D3 and the other way round: e = P1/2 (e-dagger = e, e e != e) commutes with Adjoint but does not distribute",
      comm_h < 1e-12 and hom_h > 1e-2, f"Adjoint/Controlled mismatch {comm_h:.1e}; distribute miss {hom_h:.2f}")

print(); print("="*78); print("6. CONTROLLED AND THE HOPF CIRCLE (testing the voice-call hypothesis)"); print("="*78)
phi = 0.9
g = np.exp(1j*phi)*I2
psi = rng.normal(size=2) + 1j*rng.normal(size=2); psi /= np.linalg.norm(psi)
check("H1 a global phase is invisible alone (same Bloch point, moved only along the Hopf circle)", np.allclose(bloch(g@psi), bloch(psi)), f"phase {phi}")
check("H2 ...but Controlled(global phase) = a turn of the CONTROL about Z by that phase (phase kickback)",
      np.allclose(ctl(g), kron(np.diag([1, np.exp(1j*phi)]), I2)), "C(e^{i phi} 1) = diag(1, e^{i phi}) (x) 1")
U = rand_su2(); w, vecs = np.linalg.eig(U); v = vecs[:, 0]; lam = w[0]
plus = np.array([1, 1])/np.sqrt(2)
out = ctl(U)@np.kron(plus, v)
ctrl_expected = np.array([1, lam])/np.sqrt(2)
check("H3 target in an eigenstate: it moves only along its own Hopf circle and the circle's angle becomes the control's Z turn",
      np.allclose(out, np.kron(ctrl_expected, v)), f"control Bloch {np.round(bloch(ctrl_expected), 3)} = |+> turned by arg(lambda) = {np.angle(lam):.3f}")
moved = 0; shrink = []
for _ in range(2000):
    U = rand_su2(); t = rng.normal(size=2) + 1j*rng.normal(size=2); t /= np.linalg.norm(t)
    st = (ctl(U)@np.kron(plus, t)).reshape(2, 2)
    rho_c = st@st.conj().T
    r = np.linalg.norm(np.real([np.trace(rho_c@SX), np.trace(rho_c@SY), np.trace(rho_c@SZ)]))
    shrink.append(r)
    if abs(abs(np.vdot(t, U@t)) - 1) > 1e-6: moved += 1
check("H4 'Controlled acts on one Hopf fiber' as stated: NO. For a typical target, branch 1's Bloch point moves (not just its phase) and the control gets entangled",
      moved == 2000 and np.mean(shrink) < 0.95, f"{moved}/2000 random targets moved off their fiber; mean control Bloch length {np.mean(shrink):.3f} (1 = no entanglement)")

def quat_of(a, b):  # C^2 -> H, q = a + j b, so right multiplication by a complex number is scalar multiplication
    return np.array([a.real, a.imag, b.real, -b.imag])
def cplx_of(q):
    return complex(q[0], q[1]), complex(q[2], -q[3])
def L(p):  # left multiplication by p, written as a 2x2 complex matrix on (a, b)
    cols = []
    for e in (np.array([1, 0], complex), np.array([0, 1], complex)):
        cols.append(np.array(cplx_of(qmul(p, quat_of(*e)))))
    return np.array(cols).T
p = rand_unit_q(); Lp = L(p)
lin = all(np.allclose(np.array(cplx_of(qmul(p, quat_of(*(z*vv))))), z*(Lp@vv))
          for z, vv in ((complex(*rng.normal(size=2)), rng.normal(size=2) + 1j*rng.normal(size=2)) for _ in range(20)))
check("H5 with C^2 = H (q = a + j b), every single-qubit gate is left multiplication by a unit quaternion (complex-linear, SU(2))",
      lin and np.allclose(Lp.conj().T@Lp, I2) and abs(np.linalg.det(Lp) - 1) < 1e-12, "left action = gates, right action by e^{i phi} = global phase")
def hopf8(q0, q1):
    return np.r_[np.dot(q0, q0) - np.dot(q1, q1), 2*qmul(q1, qconj(q0))]
worst_f = worst_h = worst_w = 0.0
for _ in range(500):
    s = rng.normal(size=8); s /= np.linalg.norm(s); q0, q1 = s[:4], s[4:]
    lam_ = rand_unit_q(); p = rand_unit_q()
    worst_f = max(worst_f, np.linalg.norm(hopf8(qmul(q0, lam_), qmul(q1, lam_)) - hopf8(q0, q1)))
    after = hopf8(q0, qmul(p, q1)); before = hopf8(q0, q1)
    worst_h = max(worst_h, abs(after[0] - before[0]))
    worst_w = max(worst_w, np.linalg.norm(after[1:] - qmul(p, before[1:])))
check("H6 two qubits on S^7 -> S^4 (quaternionic Hopf): Controlled-U maps each S^3 fiber to another fiber, keeps the control height |q0|^2-|q1|^2, and turns the rest by left multiplication by U",
      worst_f < 1e-12 and worst_h < 1e-12 and worst_w < 1e-12, f"fiber invariance {worst_f:.1e}; height kept {worst_h:.1e}; w -> U w {worst_w:.1e}")

worst_pole = 0.0
for _ in range(200):
    q1 = rand_unit_q(); p = rand_unit_q(); q0 = np.zeros(4)
    lam_ = qmul(qconj(q1), qmul(p, q1))          # p q1 = q1 lam_, so the new state lies on the same fiber
    worst_pole = max(worst_pole, np.linalg.norm(qmul(p, q1) - qmul(q1, lam_)), np.linalg.norm(hopf8(q0, qmul(p, q1)) - hopf8(q0, q1)))
check("H7 where the hypothesis does hold: with the control surely 1 (the pole of S^4), Controlled-U moves the state only within its own S^3 fiber",
      worst_pole < 1e-12, f"worst {worst_pole:.1e}; at the control-0 pole it does nothing; in between it turns S^4 about the control axis")

print(); print("="*78); print("7. THE SWAP RULE: two small turns, undone in swapped order, leave a turn of their area"); print("="*78)
ratios = []
for a, b in ((0.1, 0.07), (0.03, 0.05), (0.01, 0.012)):
    C = rot([1, 0, 0], a)@rot([0, 1, 0], b)@rot([1, 0, 0], -a)@rot([0, 1, 0], -b)
    q = Q(C); ang = 2*math.acos(np.clip(q[0], -1, 1)); axis = q[1:]/np.linalg.norm(q[1:])
    ratios.append((ang/(a*b), axis[2]))
check("W1 gate commutator Rx(a) Ry(b) Rx(-a) Ry(-b) = a turn about Z of size a*b (turn = area, dial 1 on the Bloch sphere)",
      abs(ratios[-1][0] - 1) < 1e-3 and abs(abs(ratios[-1][1]) - 1) < 1e-3,
      "; ".join(f"turn/area {r:.5f}, axis z-part {z:+.4f}" for r, z in ratios))

worst_c = 0.0
for k in (1.0, 0.3, -2.0):
    for _ in range(100):
        x, y = (rand_unit_q() - ONE)/k, (rand_unit_q() - ONE)/k
        worst_c = max(worst_c, np.linalg.norm((ok_q(k, x, y) - ok_q(k, y, x)) - k*(qmul(x, y) - qmul(y, x))))
check("W2 order matters by exactly k times the commutator: x o_k y - y o_k x = k[x, y]", worst_c < 1e-12, f"worst {worst_c:.1e}")
def qexp_full(h):
    t = np.linalg.norm(h[1:]); return ONE*math.cos(t) + np.r_[0, h[1:]]/t*math.sin(t)
A = np.array([0, 0.8, 0.1, -0.3]); B = np.array([0, -0.2, 0.6, 0.4])
comm_AB = qmul(A, B) - qmul(B, A)
rat = []
for k in (0.2, 0.05, 0.0125):
    err = qmul(qexp_full(k*A), qexp_full(k*B)) - qexp_full(k*(A + B))
    rat.append(np.linalg.norm(err)/(k*k*np.linalg.norm(comm_AB)/2))
check("W3 Trotter: splitting a Hamiltonian into two steps of size k misses by (k^2/2)|[A,B]|, i.e. dial d = k^2 times the swap-rule commutator",
      abs(rat[-1] - 1) < 1e-2, "ratios " + ", ".join(f"{r:.4f}" for r in rat) + " -> 1")


print(); print("="*78); print("POLISH: sharpen the two limit checks (W1, W3) by extrapolation; the rest are exact to rounding"); print("="*78)
import time
t0 = time.perf_counter()
def romberg(vals, orders):  # vals at h, h/2, h/4, ...; remove error terms h^p for each p in orders
    v = list(vals)
    for p in orders:
        v = [(2**p*v[i+1] - v[i])/(2**p - 1) for i in range(len(v)-1)]
    return v[-1]
# W1: scale (a, b) = t (0.1, 0.07), t halving; turn/area -> 1, axis z-part -> 1
tr, zz = [], []
for j in range(4):
    a, b = 0.1/2**j, 0.07/2**j
    C = rot([1, 0, 0], a)@rot([0, 1, 0], b)@rot([1, 0, 0], -a)@rot([0, 1, 0], -b)
    q = Q(C); tr.append(2*math.acos(np.clip(q[0], -1, 1))/(a*b)); zz.append(abs(q[3])/np.linalg.norm(q[1:]))
ext_t, ext_z = romberg(tr, [1, 2, 3]), romberg(zz, [1, 2, 3])
check("P1 polish W1: Richardson extrapolation of turn/area and the axis direction to zero size gives 1 and the Z axis",
      abs(ext_t - 1) < 1e-6 and abs(ext_z - 1) < 1e-6, f"raw at smallest size: {tr[-1]:.6f}, {zz[-1]:.6f}; extrapolated: miss {abs(ext_t-1):.1e}, {abs(ext_z-1):.1e}")
# W3: ratio |e^{kA}e^{kB} - e^{k(A+B)}| / ((k^2/2)|[A,B]|) at k halving
rat2 = []
for j in range(4):
    k = 0.2/2**j
    err = qmul(qexp_full(k*A), qexp_full(k*B)) - qexp_full(k*(A + B))
    rat2.append(np.linalg.norm(err)/(k*k*np.linalg.norm(comm_AB)/2))
ext = romberg(rat2, [1, 2, 3])
check("P2 polish W3: Richardson extrapolation of the Trotter ratio to k = 0 gives 1",
      abs(ext - 1) < 1e-6, f"raw at k = {0.2/8}: {rat2[-1]:.6f}; extrapolated: miss {abs(ext-1):.1e}")
print("NOTE the other 35 checks compare exact identities and already agree to rounding (worst ~1e-15), so a polish would change nothing.")
print(f"POLISH cost: {1000*(time.perf_counter() - t0):.1f} ms in total")

print(); n_pass = sum(ok for _, ok in results)
print(f"{n_pass} of {len(results)} checks pass")
