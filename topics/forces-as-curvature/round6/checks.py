"""DRAFT round 6 checks: the field equations from loops.

Rounds 1-5 read each force's field off a chain of short links round a small loop (law a o_k b = a + b + k ab).
Here the same loops give the ACTION, and varying the action gives the FIELD EQUATIONS:
  L. strong force: Wilson's loop action -> Yang-Mills action and Yang-Mills equations (Wilson 1974, textbook);
  M. gravity: the curvature read from loops, used as an action for an expanding universe, varied ->
     Friedmann's equations (Einstein-Hilbert 1915, Friedmann 1922, textbook).
Run: python3 checks.py  (numpy, scipy).  Writes results.json and checks_output.txt.
"""
import json, math, itertools, io
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

# ======================================================================= L. Yang-Mills from loops
say("== L. Strong force: loop action -> Yang-Mills action -> Yang-Mills equations (4D, Euclidean, SU(3))")
Lm = np.zeros((8, 3, 3), complex)
Lm[0][0, 1] = Lm[0][1, 0] = 1; Lm[1][0, 1] = -1j; Lm[1][1, 0] = 1j; Lm[2][0, 0] = 1; Lm[2][1, 1] = -1
Lm[3][0, 2] = Lm[3][2, 0] = 1; Lm[4][0, 2] = -1j; Lm[4][2, 0] = 1j; Lm[5][1, 2] = Lm[5][2, 1] = 1
Lm[6][1, 2] = -1j; Lm[6][2, 1] = 1j; Lm[7] = np.diag([1, 1, -2]) / math.sqrt(3)
T = Lm / 2; Nc = 3; g = 0.9
rng = np.random.default_rng(20261127)
modes = [(rng.normal(size=4), rng.normal(size=(4, 8)) * 0.5, rng.uniform(0, 6.3)) for _ in range(3)]
def A_smooth(mu, x):
    """a smooth SU(3) field: A_mu = A_mu^a T_a, sum of a few waves"""
    c = sum(C[mu] * math.sin(kv @ x + ph) for kv, C, ph in modes)
    return np.einsum("a,aij->ij", c, T)
def comm(X, Y):
    return X @ Y - Y @ X
def dA(Af, mu, nu, x, h=1e-4):
    e = np.zeros(4); e[nu] = h
    return (Af(mu, x + e) - Af(mu, x - e)) / (2 * h)
def Fcont(Af, mu, nu, x):
    return dA(Af, nu, mu, x) - dA(Af, mu, nu, x) + 1j * g * comm(Af(mu, x), Af(nu, x))
def link(Af, x, mu, a, sub=6):
    """one lattice link: a chain of `sub` short pieces exp(i g A dx), joined by the law (here: matrix product of lifts)"""
    U = np.eye(3, dtype=complex); e = np.zeros(4); e[mu] = 1
    for s in range(sub):
        U = U @ expm(1j * g * (a / sub) * Af(mu, x + e * a * (s + 0.5) / sub))
    return U
def plaq(Af, x, mu, nu, a):
    em = np.zeros(4); em[mu] = a; en = np.zeros(4); en[nu] = a
    return link(Af, x, mu, a) @ link(Af, x + em, nu, a) @ link(Af, x + en, mu, a).conj().T @ link(Af, x, nu, a).conj().T
x0 = np.array([0.31, -0.22, 0.17, 0.40])
rowsL = []
for a in (0.08, 0.04, 0.02, 0.01):
    def corner(m, n):                      # put each loop's centre at x0
        c = x0.copy(); c[m] -= a / 2; c[n] -= a / 2; return c
    wil = sum(1 - np.trace(plaq(A_smooth, corner(m, n), m, n, a)).real / Nc for m, n in itertools.combinations(range(4), 2))
    cont = sum((g**2 * a**4 / (4 * Nc)) * 2 * np.trace(Fcont(A_smooth, m, n, x0) @ Fcont(A_smooth, m, n, x0)).real
               for m, n in itertools.combinations(range(4), 2))      # F^a F^a = 2 tr F^2
    rowsL.append((a, wil / a**4, cont / a**4, abs(wil - cont) / cont))
ordL = math.log(rowsL[0][3] / rowsL[-1][3]) / math.log(rowsL[0][0] / rowsL[-1][0])
check("L1 Wilson's loop action sum(1 - Re tr U/3) / a^4 -> (g^2/12) F^a F^a: the Yang-Mills action density",
      rowsL[-1][3] < 1e-3 and ordL > 1.8, f"rel. error {rowsL[-1][3]:.1e} at a = 0.01, order {ordL:.2f}")
# the commutator matters: drop it and the action is wrong
cont_nocomm = sum((g**2 / (4 * Nc)) * 2 * np.trace(X @ X).real for X in
                  [dA(A_smooth, n, m, x0) - dA(A_smooth, m, n, x0) for m, n in itertools.combinations(range(4), 2)])
check("L2 leaving out the k.ab ([A,A]) part gives the wrong action by an order-one amount",
      abs(cont_nocomm - rowsL[-1][1]) / rowsL[-1][1] > 0.05, f"{cont_nocomm:.4f} vs {rowsL[-1][1]:.4f}")

def link_force(Af, x, nu, a):
    """d/de_b of Wilson's action S = (2N/g^2) sum_P (1 - Re tr U_P / N) when the link U_nu(x) -> exp(i e_b T_b) U_nu(x):
    the 'force' on one link. It is a sum over the 6 loops that use this link."""
    beta = 2 * Nc / g**2
    U = link(Af, x, nu, a); en = np.zeros(4); en[nu] = a
    Stap = np.zeros((3, 3), complex)
    for mu in range(4):
        if mu == nu:
            continue
        em = np.zeros(4); em[mu] = a
        # loop going forward in mu, and backward in mu
        Stap += link(Af, x + en, mu, a) @ link(Af, x + em, nu, a).conj().T @ link(Af, x, mu, a).conj().T
        Stap += link(Af, x + en - em, mu, a).conj().T @ link(Af, x - em, nu, a).conj().T @ link(Af, x - em, mu, a)
    return np.array([-(beta / Nc) * np.trace(1j * T[b] @ U @ Stap).real for b in range(8)])
def DF(Af, nu, x, h=1e-3):
    """continuum Yang-Mills equation, Euclidean: (D_mu F_mu nu)^b, with D = d + i g [A, .]"""
    out = np.zeros((3, 3), complex)
    for mu in range(4):
        e = np.zeros(4); e[mu] = h
        out += (Fcont(Af, mu, nu, x + e) - Fcont(Af, mu, nu, x - e)) / (2 * h) + 1j * g * comm(Af(mu, x), Fcont(Af, mu, nu, x))
    return np.array([2 * np.trace(out @ T[b]).real for b in range(8)])
nu0 = 1
rowsF = []
for a in (0.08, 0.04, 0.02, 0.01):
    # the force is read in the colour frame of the link's start, so it carries an error ~ a (a frame turn);
    # two link sizes (Richardson, as in round 1) remove it and leave ~ a^2
    def f_at(aa):
        xs = x0.copy(); xs[nu0] -= aa / 2
        return link_force(A_smooth, xs, nu0, aa) * g / aa**3
    rowsF.append((a, 2 * f_at(a / 2) - f_at(a)))
DFc = DF(A_smooth, nu0, x0)
sgn = float(np.sign(rowsF[-1][1] @ DFc))
errsF = [float(np.abs(fl - sgn * DFc).max() / np.abs(DFc).max()) for _, fl in rowsF]
ordF = math.log(errsF[0] / errsF[-1]) / math.log(rowsF[0][0] / rowsF[-1][0])
check("L3 the force on one link (sum over its 6 loops) x g/a^3 -> D_mu F_mu nu: varying the loop action gives the Yang-Mills equations",
      errsF[-1] < 1e-3 and ordF > 1.8, f"rel. error {errsF[-1]:.1e} at a = 0.01 (two sizes), order {ordF:.2f}, sign {sgn:+.0f}")
# a solution: field in one colour direction with constant strength (F_12 = B T_3): no force at all
Bc = 0.8
def A_abelian(mu, x):
    return (Bc * x[0]) * T[2] if mu == 1 else np.zeros((3, 3), complex)
fsol = np.abs(link_force(A_abelian, x0, 1, 0.05)).max()
check("L4 an exact solution (constant field in one colour direction): every link feels zero force", fsol < 1e-10, f"{fsol:.1e}")
# rounds 4-5's pure-commutator field (A_1 = c1 T_4, A_2 = c2 T_5, constant) is NOT a solution: the k.ab term drives it
c1v, c2v = 0.7, 0.5
def A_comm(mu, x):
    return c1v * T[3] if mu == 0 else (c2v * T[4] if mu == 1 else np.zeros((3, 3), complex))
fc = link_force(A_comm, x0, 0, 0.02) * g / 0.02**3
DFcomm = DF(A_comm, 0, x0)
check("L5 rounds 4-5's constant field made only of k.ab is not at rest: its force is the [A,[A,A]] term, nonzero",
      np.abs(DFcomm).max() > 1e-3 and np.abs(fc - sgn * DFcomm).max() / np.abs(DFcomm).max() < 2e-3,
      f"|D F| = {np.abs(DFcomm).max():.4f}, loops agree to {np.abs(fc - sgn*DFcomm).max()/np.abs(DFcomm).max():.1e}")
results["yang_mills"] = {"action_rows": rowsL, "action_order": ordL, "force_rows": [(a, list(f)) for a, f in rowsF],
                         "DF_cont": list(DFc), "force_errs": errsF, "force_order": ordF, "sign": sgn,
                         "solution_force": fsol, "comm_field_DF": list(DFcomm), "comm_field_loops": list(fc)}

# ======================================================================= M. Gravity: Friedmann from loops
say("\n== M. Gravity: curvature from loops as an action for an expanding universe; varying it gives Friedmann's equations")
def holonomy4(t0, mu, nu, eps, Gam, sub=6):
    vals = []
    def edge(dirn, start, sign):
        dd = eps / sub
        for s_ in range(sub):
            mid = start.copy(); mid[dirn] += sign * dd * (s_ + 0.5)
            vals.append(expm(-Gam(mid[0])[:, dirn, :] * sign * dd) - np.eye(4))
    p0 = np.array([t0, 0, 0, 0.0]); p0[mu] -= eps / 2; p0[nu] -= eps / 2
    p1 = p0.copy(); p1[mu] += eps; p2 = p1.copy(); p2[nu] += eps; p3 = p2.copy(); p3[mu] -= eps
    edge(mu, p0, +1); edge(nu, p1, +1); edge(mu, p2, -1); edge(nu, p3, -1)
    out = np.zeros((4, 4))
    for v in vals:
        out = v + out + v @ out                 # the law, later steps on the left
    return out
def make_gam(a, ad, N, Nd):
    """ds^2 = -N(t)^2 dt^2 + a(t)^2 (dx^2 + dy^2 + dz^2): only a, a', N, N' go in; the loops supply a''"""
    def Gam(t):
        G = np.zeros((4, 4, 4)); A_, Ad, N_, Ndt = a(t), ad(t), N(t), Nd(t)
        G[0, 0, 0] = Ndt / N_
        for i in (1, 2, 3):
            G[0, i, i] = A_ * Ad / N_**2
            G[i, 0, i] = G[i, i, 0] = Ad / A_
        return G
    return Gam
def R_loops(t0, a, ad, N, Nd, eps=1e-3):
    Gam = make_gam(a, ad, N, Nd)
    Rl = np.zeros((4, 4, 4, 4))
    for mu, nu in itertools.combinations(range(4), 2):
        H = (8 * holonomy4(t0, mu, nu, eps / 2, Gam) - holonomy4(t0, mu, nu, eps, Gam)) / eps**2
        Rl[:, :, mu, nu] = -H; Rl[:, :, nu, mu] = H
    gin = np.diag([-1 / N(t0)**2] + [1 / a(t0)**2] * 3)
    return float(np.einsum("bd,abad->", gin, Rl))
Lam = 0.6
one = lambda t: 1.0; zero = lambda t: 0.0
# generic, non-solution expansion
ag = lambda t: 1 + 0.3 * t + 0.2 * t**2; agd = lambda t: 0.3 + 0.4 * t; agdd = lambda t: 0.4
tt = 0.7
Rg = R_loops(tt, ag, agd, one, zero)
Rg_formula = 6 * (agdd(tt) / ag(tt) + agd(tt)**2 / ag(tt)**2)
check("M1 curvature of an expanding universe from loops = 6(a''/a + a'^2/a^2)", abs(Rg - Rg_formula) < 1e-5,
      f"{Rg:.7f} vs {Rg_formula:.7f}")
tgrid = np.linspace(0, 2, 81)
wq = np.ones(81); wq[1:-1:2] = 4; wq[2:-1:2] = 2; wq *= (tgrid[1] - tgrid[0]) / 3   # Simpson
def action(a, ad, N, Nd):
    """S = integral of sqrt(-g) (R - 2 Lambda) per unit comoving volume, R read from loops"""
    return float(sum(w * N(t) * a(t)**3 * (R_loops(t, a, ad, N, Nd) - 2 * Lam) for t, w in zip(tgrid, wq)))
t0b, wb = 1.0, 0.25
bump = lambda t: math.exp(-((t - t0b) / wb)**2); bumpd = lambda t: -2 * (t - t0b) / wb**2 * bump(t)
def dS_dN(a, ad, d=2e-3):
    Sp = action(a, ad, lambda t: 1 + d * bump(t), lambda t: d * bumpd(t))
    Sm = action(a, ad, lambda t: 1 - d * bump(t), lambda t: -d * bumpd(t))
    return (Sp - Sm) / (2 * d)
def dS_da(a, ad, d=2e-3):
    Sp = action(lambda t: a(t) + d * bump(t), lambda t: ad(t) + d * bumpd(t), one, zero)
    Sm = action(lambda t: a(t) - d * bump(t), lambda t: ad(t) - d * bumpd(t), one, zero)
    return (Sp - Sm) / (2 * d)
EN = lambda a, ad: float(sum(w * (6 * a(t) * ad(t)**2 - 2 * Lam * a(t)**3) * bump(t) for t, w in zip(tgrid, wq)))
dNg = dS_dN(ag, agd); ENg = EN(ag, agd)
check("M2 varying the loop action in the clock rate N gives 6 a a'^2 - 2 Lambda a^3: Friedmann's first equation, H^2 = Lambda/3",
      abs(dNg - ENg) / abs(ENg) < 1e-3, f"{dNg:.6f} vs {ENg:.6f}")
Ea = lambda a, ad, add: float(sum(w * (6 * ad(t)**2 + 12 * a(t) * add(t) - 6 * Lam * a(t)**2) * bump(t) for t, w in zip(tgrid, wq)))
dag = dS_da(ag, agd); Eag = Ea(ag, agd, agdd)
check("M3 varying it in the size a gives 6a'^2 + 12 a a'' - 6 Lambda a^2: Friedmann's second equation",
      abs(dag - Eag) / abs(Eag) < 1e-3, f"{dag:.6f} vs {Eag:.6f}")
H = math.sqrt(Lam / 3)
ads = lambda t: math.exp(H * t); adsd = lambda t: H * math.exp(H * t)
dNs, das = dS_dN(ads, adsd), dS_da(ads, adsd)
scale = abs(EN(lambda t: math.exp(1.2 * H * t), lambda t: 1.2 * H * math.exp(1.2 * H * t)))
check("M4 de Sitter expansion a = e^{Ht}, H = sqrt(Lambda/3), makes the loop action stationary both ways",
      abs(dNs) / scale < 1e-3 and abs(das) / scale < 1e-3, f"dS/dN {dNs:.1e}, dS/da {das:.1e} (vs {scale:.3f} for H x 1.2)")
sweepM = []
for fac in (0.6, 0.8, 1.0, 1.2, 1.4):
    Hf = fac * H
    af = lambda t, Hf=Hf: math.exp(Hf * t); afd = lambda t, Hf=Hf: Hf * math.exp(Hf * t)
    sweepM.append((fac, dS_dN(af, afd), EN(af, afd)))
check("M5 sweep of expansion rates: only H = sqrt(Lambda/3) is stationary; the loop answer follows Friedmann at every rate",
      all(abs(s - e) < 1e-3 * max(1, abs(e)) for _, s, e in sweepM) and abs(sweepM[2][1]) < 1e-3,
      ", ".join(f"{f:.1f}H:{s:+.3f}" for f, s, _ in sweepM))
results["friedmann"] = {"Lambda": Lam, "R_loops": Rg, "R_formula": Rg_formula, "dS_dN": dNg, "E_N": ENg, "dS_da": dag,
                        "E_a": Eag, "deSitter_dN": dNs, "deSitter_da": das, "scale": scale, "sweep": sweepM}

say(f"\n{passed} passed, {failed} failed")
results["summary"] = {"passed": passed, "failed": failed}
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
