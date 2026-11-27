"""DRAFT round 5 checks: the strong force's group SU(3) built into spacetime from a hidden SU(3) space.

Same law as round 1-4 (a o_k b = a + b + k*a*b on matrices); a chain of short links round a small
loop gives curvature.  Here the chain is 12 x 12: 4 spacetime directions + 8 hidden directions
(the group space of SU(3) with its round, bi-invariant metric).  Textbook non-abelian Kaluza-Klein
(DeWitt 1964, Kerner 1968, Cho 1975); the loop reading is ours.

Run: python3 checks.py  (numpy, scipy).  Writes results.json and checks_output.txt.
"""
import json, math, itertools, sys, io
import numpy as np
from scipy.linalg import expm

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ------------------------------------------------------------- Gell-Mann basis, T_a = lambda_a / 2
L = np.zeros((8, 3, 3), complex)
L[0][0, 1] = L[0][1, 0] = 1
L[1][0, 1] = -1j; L[1][1, 0] = 1j
L[2][0, 0] = 1; L[2][1, 1] = -1
L[3][0, 2] = L[3][2, 0] = 1
L[4][0, 2] = -1j; L[4][2, 0] = 1j
L[5][1, 2] = L[5][2, 1] = 1
L[6][1, 2] = -1j; L[6][2, 1] = 1j
L[7] = np.diag([1, 1, -2]) / math.sqrt(3)
T = L / 2
f = np.zeros((8, 8, 8))                        # [T_a, T_b] = i f_abc T_c
for a, b, c in itertools.product(range(8), repeat=3):
    f[a, b, c] = (-2j * np.trace((T[a] @ T[b] - T[b] @ T[a]) @ T[c])).real
say("== J. The strong force's group from a hidden SU(3) space (12 x 12 chain)")
CA = np.einsum("abc,abd->cd", f, f)
check("J0 SU(3) bracket: f_123 = 1, f_458 = sqrt(3)/2, and sum_ab f_abc f_abd = 3 delta_cd (colour factor C_A = 3)",
      abs(f[0, 1, 2] - 1) < 1e-12 and abs(f[3, 4, 7] - math.sqrt(3) / 2) < 1e-12 and np.abs(CA - 3 * np.eye(8)).max() < 1e-12,
      f"C_A = {CA[0,0]:.12f}")

# ------------------------------------------------------------- the hidden space: g = exp(i theta.T)
def Ug(th):
    return expm(1j * np.einsum("a,aij->ij", th, T))
def rhoF(th):
    """right-invariant forms: dg g^-1 = i T_c rho^c_m dtheta^m, rho = (e^M - 1)/M, M_cb = -theta_a f_abc"""
    M = -np.einsum("a,abc->cb", th, f)
    big = np.zeros((16, 16)); big[:8, :8] = M; big[:8, 8:] = np.eye(8)
    return expm(big)[:8, 8:]
th0 = np.array([0.31, -0.22, 0.17, 0.40, -0.13, 0.25, 0.09, -0.36])
err = 0.0; hh = 1e-6
for m in range(8):
    d = np.zeros(8); d[m] = hh
    Mm = (Ug(th0 + d) - Ug(th0 - d)) / (2 * hh) @ np.linalg.inv(Ug(th0))
    for c in range(8):
        err = max(err, abs((-2j * np.trace(T[c] @ Mm)).real - rhoF(th0)[c, m]))
check("J1 right-invariant forms from the series (e^M - 1)/M match dg g^-1", err < 1e-8, f"{err:.1e}")
rhov = 1.3
def Rvec(th):
    return np.linalg.inv(rhoF(th))             # columns: right-invariant vector fields R_a
def gamH(th):
    r = rhoF(th); return (rhov**2 / 4) * r.T @ r
KK = Rvec(th0).T @ gamH(th0) @ Rvec(th0)
check("J2 the hidden space's 8 symmetry fields have constant products R_a.R_b = (rho^2/4) delta_ab",
      np.abs(KK - rhov**2 / 4 * np.eye(8)).max() < 1e-10, f"{np.abs(KK - rhov**2/4*np.eye(8)).max():.1e}")
def lie(a, b, th, h=1e-5):
    dA = np.zeros((8, 8)); dB = np.zeros((8, 8))
    for n in range(8):
        d = np.zeros(8); d[n] = h
        dA[:, n] = (Rvec(th + d)[:, a] - Rvec(th - d)[:, a]) / (2 * h)
        dB[:, n] = (Rvec(th + d)[:, b] - Rvec(th - d)[:, b]) / (2 * h)
    return dB @ Rvec(th)[:, a] - dA @ Rvec(th)[:, b]
signs, berr = [], 0.0
for a, b in ((0, 1), (3, 4), (5, 6)):
    br = lie(a, b, th0); target = Rvec(th0) @ f[a, b]
    s = float(br @ target / (target @ target)); signs.append(s)
    berr = max(berr, np.abs(br - s * target).max())
check("J3 the symmetry fields don't commute, with SU(3)'s own bracket: [R_a, R_b] = s f_abc R_c, one sign s for all pairs",
      berr < 1e-6 and max(abs(abs(s) - 1) for s in signs) < 1e-6 and len(set(np.sign(signs))) == 1,
      f"s = {signs[0]:+.6f}, {signs[1]:+.6f}, {signs[2]:+.6f}; worst {berr:.1e}")

# ------------------------------------------------------------- 12D metric, Christoffels, loops
etan = np.diag([-1.0, 1, 1, 1])
Acon = np.zeros((8, 4))                        # A^a_mu, constant: no derivative part anywhere
def metric12(th, hidden, A):
    if hidden == "SU3":
        Rv, gm = Rvec(th), gamH(th)
    else:
        Rv, gm = np.eye(8), np.eye(8) * rhov**2 / 4      # flat 8-torus: commuting symmetries
    Nm = A.T @ Rv.T
    G = np.zeros((12, 12))
    G[:4, :4] = etan + Nm @ gm @ Nm.T
    G[:4, 4:] = -Nm @ gm; G[4:, :4] = G[:4, 4:].T; G[4:, 4:] = gm
    return G
def make_gamma(hidden, A, h=1e-5):
    def Gam(*x):
        th = np.array(x[4:]); Gm = metric12(th, hidden, A); Gi = np.linalg.inv(Gm)
        dG = np.zeros((12, 12, 12))            # dG[i, j, l] = d_l G_ij (no dependence on spacetime coords)
        if hidden == "SU3":
            for n in range(8):
                d = np.zeros(8); d[n] = h
                dG[:, :, 4 + n] = (metric12(th + d, hidden, A) - metric12(th - d, hidden, A)) / (2 * h)
        lower = 0.5 * (dG + np.transpose(dG, (0, 2, 1)) - np.transpose(dG, (2, 0, 1)))
        return np.einsum("ad,dbc->abc", Gi, lower)
    return Gam
def holonomyN(x0, mu, nu, eps, Gf, N, sub):
    """chain of short links round a small square; each link is (1 + k a) with a = exp(-Gamma dx) - 1,
    joined with the law out = v + out + v out (later steps on the left)"""
    vals = []
    def edge(dirn, start, sign, size):
        dd = size / sub
        for s_ in range(sub):
            mid = start.copy(); mid[dirn] += sign * dd * (s_ + 0.5)
            vals.append(expm(-Gf(*mid)[:, dirn, :] * sign * dd) - np.eye(N))
    p0 = np.array(x0, float); p0[mu] -= eps / 2; p0[nu] -= eps / 2
    p1 = p0.copy(); p1[mu] += eps; p2 = p1.copy(); p2[nu] += eps; p3 = p2.copy(); p3[mu] -= eps
    edge(mu, p0, +1, eps); edge(nu, p1, +1, eps); edge(mu, p2, -1, eps); edge(nu, p3, -1, eps)
    out = np.zeros((N, N))
    for v in vals:
        out = v + out + v @ out
    return out
x12 = np.concatenate([[0.0, 0.2, -0.1, 0.3], th0])
def scalar_from_loops(hidden, A, eps=2e-3, sub=6):
    Gf = make_gamma(hidden, A); N = 12
    Rl = np.zeros((N, N, N, N)); planes = {}
    for mu, nu in itertools.combinations(range(N), 2):
        H = (8 * holonomyN(x12, mu, nu, eps / 2, Gf, N, sub) - holonomyN(x12, mu, nu, eps, Gf, N, sub)) / eps**2
        Rl[:, :, mu, nu] = -H; Rl[:, :, nu, mu] = H; planes[(mu, nu)] = float(np.abs(H).max())
    gin = np.linalg.inv(metric12(x12[4:], hidden, A))
    return float(np.einsum("bd,abad->", gin, Rl)), planes
def field_sq(A):
    """F^a_{mu nu} F^{a mu nu} for a constant field: only the commutator part f_abc A^b_mu A^c_nu"""
    F = np.einsum("abc,bm,cn->amn", f, A, A)
    return float(np.einsum("amn,mp,nq,apq->", F, etan, etan, F)), F

# J4: no field. The hidden space's own curvature, read from loops, is dim x C_A / rho^2: the colour factor is curvature.
R0, planes0 = scalar_from_loops("SU3", Acon)
CA_geo = R0 * rhov**2 / 8
check("J4 no field: curvature of the hidden SU(3) space from loops = 8 C_A / rho^2, so C_A read off the geometry = 3",
      abs(CA_geo - 3) < 1e-5, f"R = {R0:.6f} vs {24/rhov**2:.6f}; C_A = {CA_geo:.6f}")

# J5: a field with A^4 = c1 dx, A^5 = c2 dy (lambda_4, lambda_5: the colours SU(2) doesn't have)
c1v, c2v = 0.7, 0.5
A45 = Acon.copy(); A45[3, 1] = c1v; A45[4, 2] = c2v
FF, F45 = field_sq(A45)
R45, planes45 = scalar_from_loops("SU3", A45)
pred45 = 24 / rhov**2 - (rhov**2 / 16) * FF
check("J5 12D curvature from loops = R(SU3) - (rho^2/16) F^a F^a, field made ONLY of the [A,A] term (lambda_4 dx, lambda_5 dy)",
      abs(R45 - pred45) < 1e-5, f"{R45:.6f} vs {pred45:.6f}")
share8 = F45[7, 1, 2]**2 / np.sum(F45[:, 1, 2]**2)
check("J6 where that field lands: 1/4 in colour 3, 3/4 in colour 8 (the direction an SU(2) space doesn't have)",
      abs(share8 - 0.75) < 1e-12 and abs(F45[2, 1, 2]**2 / np.sum(F45[:, 1, 2]**2) - 0.25) < 1e-12,
      f"colour 8 share {share8:.4f}")

# J7: control. Same field on a flat 8-torus: every loop closes.
Rt, planest = scalar_from_loops("torus", A45)
check("J7 same constant field on a flat 8-torus (commuting symmetries): every loop closes, no field at all",
      max(planest.values()) < 1e-8 and abs(Rt) < 1e-8, f"largest leftover/area {max(planest.values()):.1e}")

# J8: sweep the field strength.
sweep = []
for c2s in (0.3, 0.6, 0.9):
    As = Acon.copy(); As[3, 1] = c1v; As[4, 2] = c2s
    rs, _ = scalar_from_loops("SU3", As)
    sweep.append((c1v * c2s, rs - 24 / rhov**2, -(rhov**2 / 16) * field_sq(As)[0]))
sweep = [(0.0, R0 - 24 / rhov**2, 0.0)] + sweep + [(c1v * c2v, R45 - 24 / rhov**2, pred45 - 24 / rhov**2)]
sweep.sort()
check("J8 sweep: the field term follows -(rho^2/8)(c1 c2)^2 at every strength", all(abs(a - b) < 1e-5 for _, a, b in sweep),
      ", ".join(f"{x:.2f}:{a:+.5f}" for x, a, _ in sweep))

# J9: LEP four-jet measurements of the colour factor (Stenzel, hep-ex/0410064)
lep = {"ALEPH": (2.93, 0.60), "OPAL": (3.02, 0.55)}
check("J9 colour factor from the hidden space's curvature vs LEP four-jet values", all(abs(CA_geo - v) < e for v, e in lep.values()),
      ", ".join(f"{k} {v}+-{e} ({abs(CA_geo-v)/e:.2f} sigma)" for k, (v, e) in lep.items()))

say(f"\n{passed} passed, {failed} failed")
results.update({"C_A_from_curvature": CA_geo, "R_SU3_loops": R0, "R_SU3_formula": 24 / rhov**2, "rho": rhov,
                "R12_field_loops": R45, "R12_field_pred": pred45, "F2": FF, "colour8_share": share8,
                "torus_R": Rt, "torus_c1c2": c1v * c2v, "torus_max_leftover": max(planest.values()), "sweep": sweep, "bracket_signs": signs,
                "planes_field": {f"{a}-{b}": v for (a, b), v in planes45.items()},
                "planes_torus": {f"{a}-{b}": v for (a, b), v in planest.items()},
                "lep": lep, "summary": {"passed": passed, "failed": failed}})
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
