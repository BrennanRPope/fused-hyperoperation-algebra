"""DRAFT physics bridges (Nov 28 2026): where Affine Composition Theory meets SO(10) grand unification, superstrings,
M-theory, Kaluza-Klein, Einstein-Cartan torsion, Regge / loop quantum gravity / spin foams, Page-Wootters and thermal time,
and Pin+/Pin- topology. Earlier rounds already checked most single pieces (round 36: the law's SO(10); torsion-dial: the
twist family; action/: Ponzano-Regge and Turaev-Viro; klein-universe/: Pin+/Pin-). These checks only compute what a
bridge needs that no earlier page computed.
B1 the 16 in two signatures: Weyl spinors of SO(10) (internal) and SO(9,1) (superstring spacetime), and which is real.
B2 one generation from the 16: hypercharges, anomalies, sin^2 theta_W = 3/8, the Z4 centre of Spin(10) = the X symmetry
   5(B-L) - 4Y mod 4, and the mod-16 count with and without the right-handed neutrino.
B3 the heterotic route: E8's 240 roots split under SO(10) x SU(4) as 45 + 15 + (10,6) + (16,4) + (16bar,4bar).
B4 dimension bookkeeping: 10, 11, 26, Witten's 4 + 7, the law's 4 + 10 (arithmetic only).
B5 the twist family on the 3-sphere = SU(2) is the Cartan-Schouten family: torsion (1 - f), curvature (2f - f^2)/4 (exact).
B6 the same family on the octonions' 7-sphere: the order gap is (1 - f) [a, b], so f = 0 has torsion 2 c_ijk (the
   octonion structure constants, the form of Englert's 1982 torsion); the law f = 1 has none.
B7 Pin+ and Pin-: a reflection's spinor map squares to +1 or -1 (exact, 4D gammas).
B8 a Page-Wootters clock and the law's fused shift a -> 1 + k a: a frozen global state, conditional states that evolve,
   and the fused shift with |k| = 1 as a clock hand.
B9 Wilson lines on the law's 10-sphere: no loops; only the antipodal Z2 acts freely; RP10 is Pin- only."""
import json, io, time, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ------------------------------------------------------------------ gamma matrices
s0 = np.eye(2, dtype=complex); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1]).astype(complex)
def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms: out = np.kron(out, m)
    return out
def euclid_gammas(n):                       # n = 2m Hermitian gammas, all squaring to +1
    m = n // 2; G = []
    for j in range(m):
        G.append(kron(*([sz] * j + [sx] + [s0] * (m - j - 1))))
        G.append(kron(*([sz] * j + [sy] + [s0] * (m - j - 1))))
    return G

# ------------------------------------------------------------------ B1
say("== B1. The 16: Weyl spinor of SO(10) and of SO(9,1)")
def spinor_report(G):
    n = len(G); d = G[0].shape[0]
    chi = np.eye(d, dtype=complex)
    for g in G: chi = chi @ g
    # normalise chirality so chi^2 = 1
    c2 = (chi @ chi)[0, 0]; chi = chi / np.sqrt(c2)
    weyl = int(round(np.trace((np.eye(d) + chi) / 2).real))
    # solve B G_a = eta G_a^* B for eta = +/-1 (B unique up to scale when it exists)
    out = {}
    for eta in (+1, -1):
        rows = []
        for g in G:
            rows.append(np.kron(g.T, np.eye(d)) - eta * np.kron(np.eye(d), g.conj()))   # vec(B g) - eta vec(g* B)
        M = np.vstack(rows)
        u, s, vh = np.linalg.svd(M)
        null = vh[np.abs(s) < 1e-9] if np.sum(np.abs(s) < 1e-9) else vh[-1:][:0]
        if np.sum(s < 1e-9) == 0: out[eta] = None; continue
        B = vh[-1].conj().reshape(d, d, order="F")
        B = B / np.abs(np.linalg.det(B)) ** (1 / d)
        c = (B @ B.conj())[0, 0]
        BBc = (B @ B.conj()) / c
        majorana = bool(np.allclose(BBc, np.eye(d), atol=1e-8) and c.real > 0)
        keeps_hand = bool(np.allclose(B @ chi @ np.linalg.inv(B), chi.conj(), atol=1e-8))
        out[eta] = dict(B_Bstar_sign=float(np.sign(c.real)), keeps_hand=keeps_hand, majorana=majorana)
    mw = any(v and v["majorana"] and v["keeps_hand"] for v in out.values())
    return dict(dim_dirac=d, dim_weyl=weyl, conj=out, majorana_weyl=mw)
GE = euclid_gammas(10)
GL = [1j * GE[0]] + GE[1:]                 # one direction squares to -1: signature (9,1)
r10 = spinor_report(GE); r91 = spinor_report(GL)
ok = r10["dim_weyl"] == 16 and r91["dim_weyl"] == 16 and r91["majorana_weyl"] and not r10["majorana_weyl"]
check("B1 (numbers; textbook, e.g. Polchinski vol. 2 App. B) both SO(10) and SO(9,1) have 32-component spinors that split into two 16s of opposite hand. In SO(9,1), the superstring's spacetime, a 16 can be real (Majorana-Weyl: 16 real numbers). In SO(10), the law's internal symmetry, conjugation swaps the hands, so the 16 is complex and its conjugate is the other 16bar. Same count, different reality: same arithmetic, not the same object",
      ok, f"SO(10): Weyl dim {r10['dim_weyl']}, Majorana-Weyl {r10['majorana_weyl']} ({r10['conj']}); SO(9,1): Weyl dim {r91['dim_weyl']}, Majorana-Weyl {r91['majorana_weyl']} ({r91['conj']})")
results["B1"] = dict(SO10=r10, SO91=r91)

# ------------------------------------------------------------------ B2
say("== B2. One generation from the 16; Z4 centre; mod 16")
H = Fr(1, 2)
W16 = [w for w in itertools.product((H, -H), repeat=5) if sum(1 for x in w if x < 0) % 2 == 0]
W10 = [tuple(Fr(s) if i == j else Fr(0) for i in range(5)) for j in range(5) for s in (1, -1)]
def charges(w):
    BL = -Fr(2, 3) * (w[0] + w[1] + w[2])
    T3L = (w[3] + w[4]) / 2; T3R = (w[3] - w[4]) / 2
    Y = T3R + BL / 2
    return dict(BL=BL, T3L=T3L, T3R=T3R, Y=Y, Q=T3L + Y)
rows = [charges(w) | dict(w=w) for w in W16]
def name(c, w):
    colour = len({x for x in w[:3]}) > 1       # quarks: colour weights not all equal
    if colour: return {Fr(1, 6): "Q", Fr(-2, 3): "u^c", Fr(1, 3): "d^c"}[c["Y"]]
    return {Fr(-1, 2): "L", Fr(1): "e^c", Fr(0): "nu^c"}[c["Y"]]
content = {}
for r in rows: content[name(r, r["w"])] = content.get(name(r, r["w"]), 0) + 1
Ys = [r["Y"] for r in rows]
quark = [r for r in rows if len({x for x in r["w"][:3]}) > 1]
anom = dict(TrY=sum(Ys), TrY3=sum(y ** 3 for y in Ys),
            TrSU2sqY=sum(r["Y"] for r in rows if r["T3L"] != 0) / 2,     # each doublet member counted; index 1/2 per doublet
            TrSU3sqY=sum(r["Y"] for r in quark) / 3)                     # crude: proportional, zero iff cancels
sin2 = sum(r["T3L"] ** 2 for r in rows) / sum(r["Q"] ** 2 for r in rows)
X = [5 * r["BL"] - 4 * r["Y"] for r in rows]
centre16 = {complex(np.round(np.exp(1j * np.pi * float(sum(w))), 12)) for w in W16}
centre10 = {complex(np.round(np.exp(1j * np.pi * float(sum(w))), 12)) for w in W10}
iX = {complex(np.round(1j ** int(x), 12)) for x in X}
mod16 = dict(with_nuR=(16 * 3) % 16, without_nuR=(15 * 3) % 16)
ok = (content == {"Q": 6, "u^c": 3, "d^c": 3, "L": 2, "e^c": 1, "nu^c": 1} and all(v == 0 for v in anom.values()) and sin2 == Fr(3, 8)
      and all(x % 4 == 1 for x in X) and centre16 == {1j} and centre10 == {-1} and iX == {1j} and mod16 == dict(with_nuR=0, without_nuR=13))
check("B2 (exact; Georgi 1975, Fritzsch-Minkowski 1975; Garcia-Etxebarria-Montero 2019) one hand of SO(10)'s spinor is one generation, Q, u^c, d^c, L, e^c, nu^c, with the measured hypercharges and every anomaly cancelled; sin^2 theta_W = 3/8 where SO(10) holds. Every state carries X = 5(B-L) - 4Y = 1 mod 4, and the element of the Z4 centre of Spin(10) acts on all 16 states as i and on the 10 as -1: the X symmetry's Z4 is Spin(10)'s centre. The Z4 anomaly is counted mod 16: 3 x 16 = 48 = 0 (safe), 3 x 15 = 45 = 13 mod 16 (anomalous without the right-handed neutrino)",
      ok, f"content {content}; anomalies {{{', '.join(f'{k}: {v}' for k, v in anom.items())}}}; sin^2 theta_W = {sin2}; X values {sorted(set(X))}; centre on 16 {centre16}, on 10 {centre10}; i^X {iX}; mod 16 {mod16}")
results["B2"] = dict(content=content, anomalies={k: str(v) for k, v in anom.items()}, sin2=str(sin2), X=sorted({str(x) for x in X}), mod16=mod16,
                     states=[dict(name=name(r, r["w"]), Y=str(r["Y"]), Q=str(r["Q"]), BL=str(r["BL"])) for r in rows])

# ------------------------------------------------------------------ B3
say("== B3. Heterotic route: E8 roots under SO(10) x SU(4)")
roots = []
for i, j in itertools.combinations(range(8), 2):
    for a, b in itertools.product((1, -1), repeat=2):
        v = [0.0] * 8; v[i] = a; v[j] = b; roots.append(tuple(v))
for s in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for x in s if x < 0) % 2 == 0: roots.append(s)
R = np.array(roots)
gram_ok = set(np.round(R @ R.T, 9).ravel()) <= {-2, -1, -0.5, 0, 0.5, 1, 2} and set(np.round((R ** 2).sum(1), 9)) == {2.0}
first, last = np.abs(R[:, :5]).sum(1) > 0, np.abs(R[:, 5:]).sum(1) > 0
so10 = int(np.sum(first & ~last)); su4 = int(np.sum(~first & last))
spin = np.abs(R).min(1) == 0.5
mixed_vec = int(np.sum(first & last & ~spin)); mixed_spin = int(np.sum(spin))
# spinor roots: first-5 part is a weight of 16 (even # minus) or 16bar (odd); last-3 part then fixes 4 or 4bar
par5 = (np.sum(R[spin][:, :5] < 0, 1) % 2)
n16 = int(np.sum(par5 == 0)); n16b = int(np.sum(par5 == 1))
dims = dict(adj_SO10=so10 + 5, adj_SU4=su4 + 3, rep_10_6=mixed_vec, rep_16_4=n16, rep_16b_4b=n16b)
ok = len(roots) == 240 and gram_ok and dims == dict(adj_SO10=45, adj_SU4=15, rep_10_6=60, rep_16_4=64, rep_16b_4b=64) and sum(dims.values()) == 248
check("B3 (exact from E8's root system; textbook branching, e.g. Slansky 1981) E8's 248 = 240 roots + 8 Cartan splits under SO(10) x SU(4) into (45,1) + (1,15) + (10,6) + (16,4) + (16bar,4bar). The heterotic string's E8 broken by an SU(4) bundle on the hidden 6 dimensions leaves SO(10), and matter comes in 16s labelled by SU(4)'s 4, which is how Braun-He-Ovrut-Pantev (2005) count three generations. The law's SO(10) is the same group with the same 16, reached another way (a point on a 10-sphere, not a bundle inside E8)",
      ok, f"{len(roots)} roots, lengths and products correct {gram_ok}; pieces {dims}; total {sum(dims.values())}")
results["B3"] = dims

# ------------------------------------------------------------------ B4
say("== B4. Dimension bookkeeping (arithmetic only)")
dimtab = dict(bosonic_string=26, superstring=10, M_theory=11, sugra_max_Nahm=11,
              Witten_KK_minimal=4 + (4 + 2 + 1), law_S10_internal=4 + 10, law_S3_internal=4 + 3,
              smallest_sphere_with_SO10_isometry=9, SO11_dim=11 * 10 // 2, SO10_dim=45)
ok = dimtab["Witten_KK_minimal"] == dimtab["M_theory"] and dimtab["law_S10_internal"] == 14 and dimtab["SO11_dim"] - dimtab["SO10_dim"] == 10
check("B4 (arithmetic; Witten 1981, Nahm 1978) the smallest hidden space whose symmetries hold SU(3) x SU(2) x U(1) is CP2 x S2 x S1, 4 + 2 + 1 = 7 dimensions, so 4 + 7 = 11, the largest dimension supergravity allows: the famous coincidence behind 11D. The law's 10-sphere as a Kaluza-Klein hidden space gives 4 + 10 = 14, which matches neither 10, 11 nor 26. Its full rotations are SO(11) (55); fixing the law's point 1 leaves SO(10) (45), the other 10 rotations moving 1. A Kaluza-Klein reading would make SO(11) the gauge group, broken to SO(10) by the choice of 1",
      ok, str(dimtab))
results["B4"] = dimtab

# ------------------------------------------------------------------ B5
say("== B5. Twist family on SU(2) = Cartan-Schouten family (exact)")
f = sp.symbols("f"); lam = 1 - f / 2
Xs = [sp.Matrix(m) * sp.I / 2 for m in ([[0, 1], [1, 0]], [[0, -sp.I], [sp.I, 0]], [[1, 0], [0, -1]])]
br = lambda A, B: A * B - B * A
def nabla(A, B): return lam * br(A, B)                     # Cartan-Schouten: nabla_X Y = lam [X, Y] on left-invariant fields
def T(A, B): return nabla(A, B) - nabla(B, A) - br(A, B)
def Rm(A, B, C): return nabla(A, nabla(B, C)) - nabla(B, nabla(A, C)) - nabla(br(A, B), C)
okT = all(sp.simplify(T(A, B) - (1 - f) * br(A, B)) == sp.zeros(2) for A in Xs for B in Xs)
okR = all(sp.simplify(Rm(A, B, C) + (2 * f - f ** 2) / 4 * br(br(A, B), C)) == sp.zeros(2) for A in Xs for B in Xs for C in Xs)
# the law's member f on the 3-sphere: order gap at second order equals (1 - f) [a, b] (numeric, quaternions as SU(2))
from scipy.linalg import expm, fractional_matrix_power
rng = np.random.default_rng(1)
def su2(v): return sum(c * np.array(X, dtype=complex) for c, X in zip(v, Xs))
def member(x, y, fv): return fractional_matrix_power(x, 1 - fv / 2) @ y @ fractional_matrix_power(x, fv / 2)
gaps = {}
for fv in (0.0, 0.5, 1.0, 1.5, 2.0):
    e = 1e-3; a, b = su2(rng.normal(size=3)), su2(rng.normal(size=3))
    x, y = expm(e * a), expm(e * b)
    gap = (member(x, y, fv) - member(y, x, fv)) / e ** 2
    gaps[fv] = float(np.linalg.norm(gap - (1 - fv) * (a @ b - b @ a)) / np.linalg.norm(a @ b - b @ a))
ok = okT and okR and max(gaps.values()) < 5e-3
check("B5 (exact algebra plus numbers; Cartan and Schouten 1926) on the law's 3-sphere, which is the group SU(2), the twist family x^(1-f/2) y x^(f/2) is the Cartan-Schouten family of connections, nabla_X Y = lam [X, Y] with lam = 1 - f/2: torsion (1 - f)[X, Y] and curvature -(2f - f^2)/4 [[X, Y], Z]. The law f = 1 is the torsion-free (Levi-Civita) member; f = 0 and 2 are the flat, torsionful (teleparallel, Weitzenboeck) members. Same thing, known since 1926",
      ok, f"torsion identity {okT}; curvature identity {okR}; member order gap vs (1-f)[a,b], relative error by f: {gaps}")
results["B5"] = dict(torsion=okT, curvature=okR, order_gap_rel_err=gaps)

# ------------------------------------------------------------------ B6
say("== B6. Same family on the octonions' 7-sphere: Englert-type torsion")
def cd_mul(a, b):
    n = len(a)
    if n == 1: return a * b
    h = n // 2; a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    cj = lambda x: np.concatenate([x[:1], -x[1:]])
    return np.concatenate([cd_mul(a1, b1) - cd_mul(cj(b2), a2), cd_mul(b2, a1) + cd_mul(a2, cj(b1))])
def oexp(v):                                 # v imaginary
    t = np.linalg.norm(v); out = np.zeros(8); out[0] = np.cos(t); out[1:] = np.sin(t) * v / t if t else 0; return out
def opow(x, s):                              # x unit octonion: x = exp(theta u)
    th = np.arccos(np.clip(x[0], -1, 1)); v = x[1:]; nv = np.linalg.norm(v)
    return oexp(s * th * v / nv) if nv > 1e-15 else np.eye(8)[0]
def omember(x, y, fv): return cd_mul(cd_mul(opow(x, 1 - fv / 2), y), opow(x, fv / 2))
E = np.eye(8)
c = np.zeros((7, 7, 7))
for i in range(7):
    for j in range(7):
        p = cd_mul(E[i + 1], E[j + 1]); c[i, j] = p[1:]
anti = np.allclose(c, -np.transpose(c, (1, 0, 2))) and np.allclose(c, np.transpose(c, (1, 2, 0)))
nonzero = int(np.sum(np.abs(c) > 0.5)) // 6
og = {}
for fv in (0.0, 0.5, 1.0, 2.0):
    e = 1e-3; a, b = rng.normal(size=7), rng.normal(size=7)
    x, y = oexp(e * a), oexp(e * b)
    gap = (omember(x, y, fv) - omember(y, x, fv)) / e ** 2
    comm = cd_mul(np.r_[0, a], np.r_[0, b]) - cd_mul(np.r_[0, b], np.r_[0, a])
    og[fv] = float(np.linalg.norm(gap - (1 - fv) * comm) / np.linalg.norm(comm))
comm_is_2c = all(np.allclose(cd_mul(E[i + 1], E[j + 1]) - cd_mul(E[j + 1], E[i + 1]), np.r_[0, 2 * c[i, j]]) for i in range(7) for j in range(7) if i != j)
ok = anti and nonzero == 7 and max(og.values()) < 5e-3 and comm_is_2c
check("B6 (numbers; Cartan-Schouten 1926 on S7; Englert 1982; octonionic form: Phys. Lett. B 127 (1983) 191, Loginov) on the octonions' 7-sphere the same family's order gap is (1 - f)[a, b], and for imaginary units [e_i, e_j] = 2 c_ijk e_k with c totally antisymmetric (7 lines of the Fano plane). So the flat member f = 0 (left multiplication, one of Cartan-Schouten's two absolute parallelisms of S7) has torsion 2 c_ijk, the octonion structure constants: the form of the torsion in Englert's 11D supergravity solution. The law f = 1 has none. Same family, different member",
      ok, f"c totally antisymmetric {anti}; Fano lines {nonzero}; commutator = 2c {comm_is_2c}; order gap vs (1-f)[a,b] by f: {og}")
results["B6"] = dict(antisym=anti, fano_lines=nonzero, order_gap_rel_err=og)

# ------------------------------------------------------------------ B7
say("== B7. Pin+ and Pin-")
g4 = euclid_gammas(4)
pin_plus = g4[1]                                 # spacelike gamma squaring to +1 (Cl with spacelike^2 = +1)
pin_minus = 1j * g4[1]                           # the other Clifford convention: spacelike^2 = -1
sq = dict(pin_plus=float(np.real((pin_plus @ pin_plus)[0, 0])), pin_minus=float(np.real((pin_minus @ pin_minus)[0, 0])))
evals = dict(pin_plus=sorted(set(np.round(np.linalg.eigvals(pin_plus), 9))), pin_minus=sorted(set(np.round(np.linalg.eigvals(pin_minus).imag, 9))))
ok = sq == dict(pin_plus=1.0, pin_minus=-1.0)
check("B7 (exact; textbook, Kirby-Taylor 1990; repo klein-universe U7) a mirror acts on spinors by a gamma matrix along its normal. In Pin+ that matrix squares to +1 (once round a Klein twist: +-1; twice: +1); in Pin- it squares to -1 (once round: +-i; twice: -1, forced). Pin+ in 4D has bordism group Z16 (Kirby-Taylor 1990), the same mod-16 that B2 counts, which is how Garcia-Etxebarria-Montero tie 16 states per generation to topological superconductors",
      ok, f"squares {sq}; eigenvalues pin+ {evals['pin_plus']}, pin- (imag parts) {evals['pin_minus']}")
results["B7"] = sq

# ------------------------------------------------------------------ B8
say("== B8. Page-Wootters clock and the fused shift")
N, dt = 48, 0.2
w0 = 2 * np.pi / (N * dt)
Hs = np.diag([0.0, 3 * w0])                      # periodic over the clock's cycle
psi0 = np.array([1, 1]) / np.sqrt(2)
U = lambda t: np.diag(np.exp(-1j * np.diag(Hs) * t))
Psi = np.zeros((N, 2), complex)
for t in range(N): Psi[t] = U(t * dt) @ psi0 / np.sqrt(N)
shift = np.roll(np.eye(N), 1, axis=0)
frozen = np.allclose(np.kron(shift, U(dt)) @ Psi.ravel(), Psi.ravel())
cond = [Psi[t] / np.linalg.norm(Psi[t]) for t in range(N)]
schro = max(np.linalg.norm(cond[t + 1] - U(dt) @ cond[t]) for t in range(N - 1))
# fused shift a -> 1 + k a with |k| = 1: a turn about a* = 1/(1 - k)
k = np.exp(2j * np.pi / N); a = 0.3 + 0.1j; astar = 1 / (1 - k); traj = [a]
for _ in range(N): traj.append(1 + k * traj[-1])
radius_spread = float(np.ptp([abs(z - astar) for z in traj])); returns = abs(traj[-1] - traj[0])
ok = frozen and schro < 1e-12 and radius_spread < 1e-12 and returns < 1e-12
check("B8 (numbers; Page-Wootters 1983; ours: the fused shift as a clock hand) the global state of clock + system is unchanged by one joint tick (frozen, the Wheeler-DeWitt picture), yet the system conditioned on each clock reading evolves by Schroedinger's law. The law's fused shift a -> 1 + k a with |k| = 1 turns every point about a* = 1/(1 - k) by arg k and returns after N ticks when k^N = 1, the same cycle as the clock's shift. A resemblance with shared arithmetic (rotation by a root of unity), not a derivation of quantum time",
      ok, f"frozen under joint tick {frozen}; conditional states vs Schroedinger {schro:.1e}; fused-shift radius spread {radius_spread:.1e}; returns after N={N} ticks {returns:.1e}")
results["B8"] = dict(frozen=frozen, schrodinger_err=schro, radius_spread=radius_spread, N=N)

# ------------------------------------------------------------------ B9
say("== B9. Can the heterotic trick (Wilson lines) live on the law's 10-sphere?")
from math import comb
def sw_structures(n):                         # RP^n: w1 = (n+1)a, w2 = C(n+1,2) a^2, with a^2 != 0 iff n >= 2
    w1 = (n + 1) % 2; w2 = comb(n + 1, 2) % 2
    orient = w1 == 0; spin = orient and w2 == 0
    pin_plus = w2 == 0; pin_minus = (w2 + w1 * w1) % 2 == 0
    return dict(orientable=orient, spin=spin, pin_plus=pin_plus, pin_minus=pin_minus)
tab = {n: sw_structures(n) for n in range(2, 13)}
antipode_moves_1 = bool(np.allclose(-np.eye(11)[0], np.eye(11)[0]) is False)
ok = tab[10] == dict(orientable=False, spin=False, pin_plus=False, pin_minus=True) and tab[2]["pin_minus"] and not tab[2]["pin_plus"] and antipode_moves_1
check("B9 (exact Stiefel-Whitney arithmetic; Lefschetz; ours: the reading) the heterotic models break SO(10) with Wilson lines, which need loops in the hidden space. The 10-sphere has none. The only group that can act freely on an even sphere is Z2, the antipodal map (Lefschetz fixed-point theorem), giving RP10: not orientable, no Spin or Pin+ structure, only Pin- (w1 = a, w2 = a^2, w2 + w1^2 = 0). And the antipodal map sends the law's point 1 to -1, so it is not a symmetry of the law as set up. So the heterotic road's way of reaching the middle stage has no copy on the law's sphere",
      ok, f"RP10 {tab[10]}; RP2 {tab[2]} (known: Pin-); antipode moves 1 {antipode_moves_1}; table {{{', '.join(f'{n}: ' + ''.join(c if v else '.' for c, v in zip('os+-', t.values())) for n, t in tab.items())}}} (letters o s + - = orientable, spin, Pin+, Pin-; dot = no)")
results["B9"] = tab

say(f"\n{passed} / {passed + failed} pass in {time.time() - T0:.1f} s")
results["summary"] = dict(passed=passed, failed=failed)
open(__file__.replace("checks.py", "checks_output.txt"), "w").write(_log.getvalue())
json.dump(results, open(__file__.replace("checks.py", "results.json"), "w"), indent=1, default=str)
