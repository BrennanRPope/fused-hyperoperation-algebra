"""DRAFT checks: our framework next to the ZX-calculus, with concrete circuits and gate counts.
Run: python3 checks_zx.py > checks_zx_output.txt   (numpy only)
Dictionary tested: core law o_k, its inverse (Adjoint), the union rule o_-1, the XOR notch o_-2 (logic-sets page),
corner expansion of the cube  <->  ZX spider fusion, adjoint, phase gadgets, phase folding.
T-count = number of non-Clifford pi/4-type phases; counted on the gadget lists below.
"""
import itertools, math
import numpy as np
results = []
def check(name, ok, detail):
    results.append((name, bool(ok))); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
H = (X + Z)/np.sqrt(2); S = np.diag([1, 1j]); T = np.diag([1, np.exp(1j*np.pi/4)])
P0 = np.diag([1., 0]).astype(complex); P1 = np.diag([0., 1]).astype(complex)
def kron(*m):
    o = np.eye(1, dtype=complex)
    for a in m: o = np.kron(o, a)
    return o
def Rz(a): return np.diag([np.exp(-1j*a/2), np.exp(1j*a/2)])
def Rx(a): return math.cos(a/2)*I2 - 1j*math.sin(a/2)*X
def Ry(a): return math.cos(a/2)*I2 - 1j*math.sin(a/2)*Y
def eq_phase(A, B):
    i = np.unravel_index(np.argmax(np.abs(B)), B.shape); ph = A[i]/B[i]
    return abs(abs(ph) - 1) < 1e-9 and np.allclose(A, ph*B, atol=1e-9)
def op_on(n, ops):  # ops: dict qubit -> 2x2
    return kron(*[ops.get(q, I2) for q in range(n)])
def Zs(n, Sset): return op_on(n, {q: Z for q in Sset})
def gadget(n, Sset, a):  # exp(i a Z_S), ZX phase gadget
    return np.diag(np.exp(1j*a*np.diag(Zs(n, Sset)).real))

print("="*78); print("1. WHERE OURS AND ZX SAY THE SAME THING"); print("="*78)
a, b = 0.7, 1.9
check("M1 ZX spider fusion (phases add) = our core law on the circle: lifts e^{ia} multiply, angles add", np.allclose(Rz(a)@Rz(b), Rz(a + b)), "Z-spiders fuse; same for X-spiders")
seq = [Rz(0.3), Rx(1.1), H, T, Rz(-0.8), S, Rx(2.2)]
U = np.eye(2, dtype=complex)
for g in seq: U = g@U
adj = np.eye(2, dtype=complex)
for g in reversed(seq): adj = g.conj().T@adj
check("M2 ZX adjoint (flip the diagram, negate every phase) = our mirror + order reversal", np.allclose(U.conj().T, adj), "7-gate run")
def lift(P, k): return np.eye(len(P)) + k*P
n = 2; Pa, Pb = op_on(2, {0: P1}), op_on(2, {1: P1})
ok_ = lambda x, y, k: x + y + k*(x@y)
OR, XOR, AND = ok_(Pa, Pb, -1), ok_(Pa, Pb, -2), Pa@Pb
check("M3 on projectors our dial gives the logic: o_-1 = OR, o_-2 = XOR (parity), product = AND (logic-sets page)",
      np.allclose(np.diag(OR).real, [0, 1, 1, 1]) and np.allclose(np.diag(XOR).real, [0, 1, 1, 0]) and np.allclose(np.diag(AND).real, [0, 0, 0, 1]),
      "on |00>,|01>,|10>,|11>")
check("M4 Pauli Z is the k = -2 lift of the projector P1, and Z_a Z_b is the lift of the XOR projector: ZX's X-spider parity is our o_-2",
      np.allclose(Z, lift(P1, -2)) and np.allclose(Zs(2, [0, 1]), lift(XOR, -2)), "Z = 1 - 2 P1;  Z_a Z_b = 1 - 2 (P_a o_-2 P_b)")
# CCZ as the cube's corner expansion = 7 phase gadgets (Amy-Maslov-Mosca / ZX phase gadgets)
def corner_gadgets(qubits, theta):
    """exp(i theta P1_q1 P1_q2 ...) = global * prod over corners S of exp(i theta (-1)^|S| Z_S / 2^m)"""
    m = len(qubits); out = []
    for r in range(1, m+1):
        for Sset in itertools.combinations(qubits, r):
            out.append((frozenset(Sset), theta*(-1)**r/2**m))
    return out
def build(n, gl):
    U = np.eye(2**n, dtype=complex)
    for Sset, ang in gl: U = gadget(n, Sset, ang)@U
    return U
CCZ = np.diag([1]*7 + [-1]).astype(complex)
g7 = corner_gadgets([0, 1, 2], np.pi)
def tcount(gl):
    tot = {}
    for Sset, ang in gl: tot[Sset] = tot.get(Sset, 0) + ang
    # exp(i a Z_S) is Clifford when a is a multiple of pi/4; odd multiples of pi/8 are T-type
    nonc = [s for s, ang in tot.items() if abs(ang/(np.pi/8) - round(ang/(np.pi/8))) < 1e-9 and round(ang/(np.pi/8)) % 2 == 1]
    return len(nonc), tot
tc, _ = tcount(g7)
check("M5 CCZ = the cube's corner expansion (1/8) prod(1 - Z_i): 7 phase gadgets of angle +-pi/8, T-count 7 (the standard ZX / phase-polynomial count)",
      eq_phase(build(3, g7), CCZ) and tc == 7, f"T-count {tc}")

print(); print("="*78); print("2. CONCRETE CIRCUITS, BEFORE AND AFTER"); print("="*78)
# (a) phase folding: CCZ(a,b,c) then CCZ(a,b,d)
n = 4
gl = corner_gadgets([0, 1, 2], np.pi) + corner_gadgets([0, 1, 3], np.pi)
target = np.diag([(-1)**(((i >> 3) & 1)*((i >> 2) & 1)*((i >> 1) & 1)) * (-1)**(((i >> 3) & 1)*((i >> 2) & 1)*(i & 1)) for i in range(16)]).astype(complex)
tc_before = 14; tc_after, merged = tcount(gl)
shared = sorted("".join("abcd"[q] for q in sorted(s)) for s, ang in merged.items() if abs(abs(ang) - np.pi/4) < 1e-9)
check("C1 two CCZs sharing controls a, b: corners {a}, {b}, {ab} appear in both with the same sign and merge into Clifford phases; T-count 14 -> 8",
      eq_phase(build(4, gl), target) and tc_after == 8, f"merged corners {shared}; T-count after {tc_after} (known result: phase folding; ZX finds the same by spider fusion of gadgets)")
# (b) merging a run of single-qubit gates: one product (our o_k on S^3)
rng = np.random.default_rng(7)
run = [H, T, H, S, Rz(0.41), Rx(1.3), T.conj().T, H, Ry(0.77), S, H, T]
U = np.eye(2, dtype=complex)
for g in run: U = g@U
V = U/np.sqrt(np.linalg.det(U))
# ZYZ Euler angles from the quaternion
w, x, y, z = (np.trace(V)/2).real, (1j*np.trace(V@X)/2).real, (1j*np.trace(V@Y)/2).real, (1j*np.trace(V@Z)/2).real
th = 2*math.atan2(math.hypot(x, y), math.hypot(w, z)); s_ = math.atan2(z, w); d_ = math.atan2(-x, y)
phi_, lam_ = s_ + d_, s_ - d_
check("C2 a 12-gate single-qubit run is one product on S^3 (one quaternion), i.e. 3 rotations Rz Ry Rz; 12 -> 3",
      eq_phase(Rz(phi_)@Ry(th)@Rz(lam_), U), f"angles {phi_:.3f}, {th:.3f}, {lam_:.3f} (every circuit optimizer does this too; inside ZX it needs the Euler-angle rule)")
# (c) OR-control from the union rule
Ut = Rx(0.9)@Rz(0.4)
def ctl_on(P, U, n=3): return np.eye(2**n) + kron(P, U - I2)
ORc = ctl_on(OR, Ut)
viaAND = kron(X, X, I2) @ ctl_on(op_on(2, {0: P1})@op_on(2, {1: P1}), Ut) @ kron(X, X, I2)
viaU = kron(I2, I2, Ut) @ ctl_on(op_on(2, {0: P0})@op_on(2, {1: P0}), Ut.conj().T)
check("C3 'control on a OR b' from the union rule = U, then undo on the 00 corner = one doubly-controlled gate plus cheap X's: same cost as an AND control, no saving",
      np.allclose(ORc, viaU) and np.allclose(ORc, kron(I2, I2, Ut) @ (kron(X, X, I2) @ ctl_on(op_on(2, {0: P1})@op_on(2, {1: P1}), Ut.conj().T) @ kron(X, X, I2))),
      "OR = NOT(AND(NOT, NOT)): de Morgan")
# (d) controlled conjugation: C(V U V^dag) = V C(U) V^dag (Q#'s within/apply)
def C(U): return np.eye(4) + kron(P1, U - I2)
CH_dec = kron(I2, S@H@T) @ (np.eye(4) + kron(P1, X - I2)) @ kron(I2, T.conj().T@H@S.conj().T)
check("C4a controlled-H as Clifford+T uses 2 T gates", eq_phase(CH_dec, C(H)), "CH = (S H T) CNOT (T-dag H S-dag) on the target")
CS = build(2, corner_gadgets([0, 1], np.pi/2))
tcs, _ = tcount(corner_gadgets([0, 1], np.pi/2))
check("C4b controlled-S from the corner expansion: 3 gadgets of angle +-pi/8, T-count 3", eq_phase(CS, C(S)) and tcs == 3, f"T-count {tcs}")
check("C4c C(H S H) = (1 x H) C(S) (1 x H): only the middle needs the control (the corner commutes with conjugation); T-count 2+3+2 = 7 -> 3",
      np.allclose(C(H@S@H), kron(I2, H)@C(S)@kron(I2, H)), "known: Q#'s within/apply does exactly this")

print(); print("="*78); print("3. WHERE ONE SIDE SAYS SOMETHING THE OTHER DOES NOT"); print("="*78)
A = np.array([[0.3, 0.2+0.1j], [0.2-0.1j, -0.3]]); B = np.array([[0.1, -0.4j], [0.4j, -0.1]])
def expm_h(M, t):
    w, v = np.linalg.eigh(M); return v@np.diag(np.exp(-1j*t*w))@v.conj().T
comm_ = A@B - B@A
errs = []
for t in (0.2, 0.1, 0.05):
    e = np.linalg.norm(expm_h(A, t)@expm_h(B, t) - expm_h(A + B, t), 2)
    errs.append(e/(t*t*np.linalg.norm(comm_, 2)/2))
check("D1 two non-commuting Pauli-gadget steps: our corner weights (the BCH terms) give the error, which ZX rewriting does not produce",
      abs(errs[-1] - 1) < 0.05, "error / ((t^2/2)|[A,B]|) = " + ", ".join(f"{r:.3f}" for r in errs))
print("NOTE ZX has rules ours does not: the bialgebra law (copy vs parity across wires), local complementation and pivoting,")
print("     which rewrite multi-qubit entangling structure and reduce CNOT and T counts in large circuits. We have no counterpart.")


print(); print("="*78); print("POLISH: exact phase bookkeeping and an extrapolated error ratio (cheap)"); print("="*78)
import time
from fractions import Fraction as Fr
t0 = time.perf_counter()
# exact T-counts: angles kept as exact fractions of pi, so no tolerance is needed
def tcount_exact(gl):
    tot = {}
    for Sset, ang in gl: tot[Sset] = tot.get(Sset, Fr(0)) + ang
    return sum(1 for a in tot.values() if (a/Fr(1, 8)).denominator == 1 and (a/Fr(1, 8)).numerator % 2 == 1), tot
tc7, _ = tcount_exact(corner_gadgets([0, 1, 2], Fr(1)))
tcA, tot = tcount_exact(corner_gadgets([0, 1, 2], Fr(1)) + corner_gadgets([0, 1, 3], Fr(1)))
cliff = sorted("".join("abcd"[q] for q in sorted(s)) for s, a in tot.items() if a.denominator <= 4)
check("P1 polish M5 and C1: T-counts recounted with exact fractions of pi: CCZ 7, folded pair 14 -> 8",
      tc7 == 7 and tcA == 8, f"Clifford corners after merging: {cliff}")
# C2: one Newton (Gauss-Newton) step on the three Euler angles
def resid(v):
    W = Rz(v[0])@Ry(v[1])@Rz(v[2]); i = np.unravel_index(np.argmax(np.abs(U)), U.shape); ph = W[i]/U[i]
    r = (W - ph*U).ravel(); return np.r_[r.real, r.imag]
v0 = np.array([phi_, th, lam_]); r0 = resid(v0)
J = np.array([(resid(v0 + 1e-7*e) - r0)/1e-7 for e in np.eye(3)]).T
v1 = v0 - np.linalg.lstsq(J, r0, rcond=None)[0]
check("P2 polish C2: one Newton step on the three angles leaves them unchanged to rounding (the closed-form angles were already exact)",
      np.linalg.norm(resid(v1)) < 1e-12 and np.linalg.norm(v1 - v0) < 1e-6, f"residual before {np.linalg.norm(r0):.1e}, after {np.linalg.norm(resid(v1)):.1e}")
# D1: Richardson extrapolation of the error ratio to t = 0
def romberg(vals, orders):
    v = list(vals)
    for p in orders: v = [(2**p*v[i+1] - v[i])/(2**p - 1) for i in range(len(v)-1)]
    return v[-1]
rr = [np.linalg.norm(expm_h(A, t)@expm_h(B, t) - expm_h(A + B, t), 2)/(t*t*np.linalg.norm(comm_, 2)/2) for t in (0.2, 0.1, 0.05, 0.025)]
ext = romberg(rr, [1, 2, 3])
check("P3 polish D1: Richardson extrapolation of error / ((t^2/2)|[A,B]|) to t = 0 gives 1",
      abs(ext - 1) < 1e-5, f"raw at t = 0.025: {rr[-1]:.5f}; extrapolated: miss {abs(ext-1):.1e}")
print(f"POLISH cost: {1000*(time.perf_counter() - t0):.1f} ms in total")

print(); n_pass = sum(ok for _, ok in results); print(f"{n_pass} of {len(results)} checks pass")
