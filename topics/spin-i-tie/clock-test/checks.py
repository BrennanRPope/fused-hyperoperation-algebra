"""Clock test checks (DRAFT, Nov 27 2026).

Run:  python3 checks.py > checks_output.txt     (writes results.json)

Part A, a clock in place of the energy zero
  K1  Page-Wootters: a timeless state H_total Psi = 0 gives Schrodinger evolution for the system
  K2  moving a constant between clock and system changes nothing (why-i's A4, now derived)
  K3  with real numbers and separate quarter turns the same construction fails; a shared turn fixes it
  K4  gravity as the clock (Born-Oppenheimer / WKB): the i comes from the gravitational phase,
      error falls like 1/M; a real state holds both i and -i
Part B, is the law's fused shift the clock?
  K5  de Sitter minisuperspace: the real (Hartle-Hawking-type) state does not oscillate below the
      law's slice size 1/sqrt(q) and oscillates above it, at q = 1/4, 1, 4
  K6  the clock that the oscillation defines grows as a = cosh(sqrt(q) tau)/sqrt(q) (the fused shift
      with real k, round 16), and matter carried along it turns by the fused shift with imaginary k
Part C, can a small clock nobody can read make the probabilities?
  P1  an unknown global clock phase changes no probability
  P2  an unknown relative phase turns interference into a mixture with weights |a|^2 (weights put in, not made)
  P3  Zurek's envariance: local phases are invisible; equal amplitudes -> equal odds; fine-graining gives 2/3, 1/3
  P4  Bell: a local clock phase shared at the source gives CHSH <= 2; quantum mechanics 2 sqrt 2
  P5  the rest-energy (Compton) clock: its absolute phase cancels between paths; gravity leaves m g h t / hbar
"""
import json
import numpy as np
import sympy as sp
from scipy.linalg import expm, null_space
from scipy.integrate import solve_ivp

rng = np.random.default_rng(271127)
results, passed, total = {}, 0, 0


def check(name, ok, detail, **extra):
    global passed, total
    total += 1
    passed += bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    results[name] = {"pass": bool(ok), "detail": detail, **extra}


def realify(M):
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


# ---------------------------------------------------------------- K1, K2
d, n = 8, 3                       # clock levels, system levels
w = 2 * np.pi / d                 # clock tick so that d ticks close the loop
HC = -np.diag(np.arange(d)) * w   # ideal clock (Page-Wootters / Salecker-Wigner); its energy balances the system's
ES = np.array([1, 3, 6]) * w      # system energies on the clock's ladder
U = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]
HS = U @ np.diag(ES) @ U.conj().T
c0 = np.ones(d) / np.sqrt(d)
psi0 = rng.normal(size=n) + 1j * rng.normal(size=n)
psi0 /= np.linalg.norm(psi0)


def clock_state(t, HCx):
    return expm(-1j * HCx * t) @ c0


def history(HCx, HSx, ps):
    # Psi = sum over a full loop of clock states (x) system states, with H_total = HC (x) 1 + 1 (x) HS
    Psi = np.zeros(d * n, complex)
    for t in range(d):
        Psi += np.kron(clock_state(t, HCx), expm(-1j * HSx * t) @ ps)
    return Psi / np.linalg.norm(Psi)


Htot = np.kron(HC, np.eye(n)) + np.kron(np.eye(d), HS)  # constraint: clock energy + system energy = 0
Psi = history(HC, HS, psi0)
res_constraint = np.linalg.norm(Htot @ Psi)


def conditional(Psi, t, HCx):
    ct = clock_state(t, HCx)
    v = np.kron(ct.conj(), np.eye(n)) @ Psi
    return v / np.linalg.norm(v)


fid = [abs(np.vdot(conditional(Psi, t, HC), expm(-1j * HS * t) @ psi0)) for t in np.linspace(0, d - 1, 15)]
check("K1", res_constraint < 1e-12 and min(fid) > 1 - 1e-12,
      f"timeless state: |H_total Psi| = {res_constraint:.1e}; conditioned on the clock reading t, the system "
      f"follows e^(-i H_S t) psi0 with fidelity >= {min(fid):.15f} at 15 readings (fractional t included)")

# K2: move c between clock and system: HC + c, HS - c leaves the balance HC + HS unchanged
probs0 = [np.abs(conditional(Psi, t, HC)) ** 2 for t in np.linspace(0, 7, 8)]
dmax = 0
for c in [0.37, 2.0, -5.5]:
    Psi_c = history(HC + c * np.eye(d), HS - c * np.eye(n), psi0)
    Htot_c = np.kron(HC + c * np.eye(d), np.eye(n)) + np.kron(np.eye(d), HS - c * np.eye(n))
    assert np.allclose(Htot_c, Htot)
    for t, p0 in zip(np.linspace(0, 7, 8), probs0):
        # same physical state, read with the shifted clock
        dmax = max(dmax, np.abs(np.abs(conditional(Psi, t, HC + c * np.eye(d))) ** 2 - p0).max())
check("K2", dmax < 1e-12,
      f"booking a constant c (0.37, 2, -5.5) to both clock and system leaves the constraint identical and moves "
      f"every conditional probability by at most {dmax:.1e}: the energy zero is a convention because only the "
      f"balance is constrained (why-i's A4 and A6, here derived from the constraint)")

# ---------------------------------------------------------------- K3
JC = realify(1j * np.eye(d)); JS = realify(1j * np.eye(n))
HCr, HSr = realify(HC.astype(complex)), realify(HS)          # real symmetric energies
AC, AS = JC @ HCr * -1, JS @ HSr * -1                          # rotation generators, e^{A t} real rotations
assert np.allclose(AC, -AC.T) and np.allclose(AS, -AS.T)
c0r = np.concatenate([c0, np.zeros(d)])
p0r = np.concatenate([psi0.real, psi0.imag])
Psi_r = sum(np.kron(expm(AC * t) @ c0r, expm(AS * t) @ p0r) for t in range(d))
Psi_r /= np.linalg.norm(Psi_r)
Htot_r = np.kron(HCr, np.eye(2 * n)) + np.kron(np.eye(2 * d), HSr)
viol_plain = np.linalg.norm(Htot_r @ Psi_r)
# balanced product: impose JC (x) 1 = 1 (x) JS (Bang-Cho-Baek; Stueckelberg's shared J)
Jl, Jr = np.kron(JC, np.eye(2 * n)), np.kron(np.eye(2 * d), JS)
Pbal = (np.eye(4 * d * n) - Jl @ Jr) / 2                     # projector onto Jl = Jr  (Jl Jr = -1 there)
assert np.allclose(Pbal @ Pbal, Pbal)
frac_bal = np.linalg.norm(Pbal @ Psi_r) ** 2
Psi_b = Pbal @ Psi_r / np.linalg.norm(Pbal @ Psi_r)
viol_bal = np.linalg.norm(Htot_r @ Psi_b)
Panti = np.eye(4 * d * n) - Pbal
Psi_a = Panti @ Psi_r / np.linalg.norm(Panti @ Psi_r)
H_diff = np.kron(HCr, np.eye(2 * n)) - np.kron(np.eye(2 * d), HSr)
ev_d, V_d = np.linalg.eigh(H_diff)
cw = (V_d.T @ Psi_a) ** 2
off = np.abs(ev_d - 2 * np.pi * np.round(ev_d / (2 * np.pi)))
viol_anti_sum = float(np.sqrt(cw[off > 1e-8].sum()))          # weight off the full-turn ladder
mirror_turns = sorted(set(np.round(ev_d[cw > 1e-10] / (2 * np.pi)).astype(int).tolist()))
anti_breaks = np.linalg.norm(Htot_r @ Psi_a)
# conversely: a generic real solution of the constraint does not give one Schrodinger flow
K = null_space(Htot_r)
Pr = K @ rng.normal(size=K.shape[1]); Pr /= np.linalg.norm(Pr)


def cond_real(P, t):
    v = np.kron((expm(AC * t) @ c0r)[None, :], np.eye(2 * n)) @ P
    return v / np.linalg.norm(v)


v0 = cond_real(Pr, 0)
flow_err_plain = max(np.linalg.norm(cond_real(Pr, t) - expm(AS * t) @ v0) for t in [1, 2, 3])
Pb2 = Pbal @ Pr; Pb2 /= np.linalg.norm(Pb2)
v0b = cond_real(Pb2, 0)
flow_err_bal = max(np.linalg.norm(cond_real(Pb2, t) - expm(AS * t) @ v0b) for t in [1, 2, 3])
check("K3", viol_plain > 0.3 and viol_bal < 1e-10 and viol_anti_sum < 1e-10 and anti_breaks > 0.3 and flow_err_plain > 0.1 and flow_err_bal < 1e-10,
      f"real numbers, plain product, separate quarter turns: the clocked history breaks the constraint "
      f"(|H_total Psi| = {viol_plain:.2f}); its balanced part ({frac_bal:.3f} of the weight) obeys it to {viol_bal:.1e}, "
      f"and the rest breaks it ({anti_breaks:.2f}) and obeys H_C - H_S = 0 up to whole turns instead (off-ladder weight {viol_anti_sum:.1e}, turns {mirror_turns}): a mirror sector where the system runs with "
      f"-i against the clock. A generic real solution of the constraint read by the clock leaves one rotation flow "
      f"by {flow_err_plain:.2f}; inside the shared-turn (balanced) part by {flow_err_bal:.1e}")

# ---------------------------------------------------------------- K4
Hm = np.array([[0.3, 0.8], [0.8, -0.5]])        # matter, real symmetric
Wg = 0.5
chi0 = np.array([0.6, 0.8])
taus = np.linspace(0.5, 4, 8)
errs = {}
for M in [10.0, 100.0, 1000.0]:
    e_m, V = np.linalg.eigh(Hm)
    pm = np.sqrt(2 * M * (M * Wg - e_m))        # exact branch momenta of  (-1/2M d^2/da^2 - M W + Hm) Psi = 0
    worst = 0
    for tau in taus:
        a = tau * np.sqrt(2 * Wg)               # WKB time tau = a / sqrt(2W)
        branch = V @ (np.exp(1j * pm * a) * (V.T @ chi0))
        branch *= np.exp(-1j * M * np.sqrt(2 * Wg) * a)   # strip the common gravitational phase e^{i M S}
        target = expm(-1j * Hm * tau) @ chi0
        worst = max(worst, np.linalg.norm(branch - target))
    errs[M] = worst
slope = np.polyfit(np.log10(list(errs)), np.log10(list(errs.values())), 1)[0]
# real state: cos branch = (forward + backward)/2 ; backward branch follows e^{+i Hm tau}
M = 1000.0
e_m, V = np.linalg.eigh(Hm); pm = np.sqrt(2 * M * (M * Wg - e_m))
tau = 2.0; a = tau * np.sqrt(2 * Wg)
back = V @ (np.exp(-1j * pm * a) * (V.T @ chi0)) * np.exp(1j * M * np.sqrt(2 * Wg) * a)
err_back = np.linalg.norm(back - expm(+1j * Hm * tau) @ chi0)
real_state = V @ (np.cos(pm * a) * (V.T @ chi0))
check("K4", slope < -0.95 and errs[1000.0] < 5e-3 and err_back < 5e-3 and np.isrealobj(real_state),
      f"gravity as the clock: along the e^(+iMS) branch the matter state follows e^(-i H_m tau) with error "
      f"{', '.join(f'{v:.1e}' for v in errs.values())} for M = 10, 100, 1000 (slope {slope:.2f}: error ~ 1/M); "
      f"the e^(-iMS) branch follows e^(+i H_m tau) (error {err_back:.1e}); the real solution cos(pa) is their sum: "
      f"the i in Schrodinger's equation is the orientation of the gravitational phase (Banks 1985; Kiefer-Singh 1991)",
      errs={str(k): float(v) for k, v in errs.items()})

# ---------------------------------------------------------------- K5
eps = 0.02
k5 = {}
for q in [0.25, 1.0, 4.0]:
    a_t = 1 / np.sqrt(q)
    f = lambda a, y: [y[1], (a ** 2 - q * a ** 4) / eps ** 2 * y[0]]
    # regular, real start at a = 0 (Hartle-Hawking-type: real, smooth at a = 0); scaled to avoid overflow
    sol = solve_ivp(f, [0, 2.2 * a_t], [1.0, 0.0], max_step=a_t / 4000, rtol=1e-10, atol=1e-14, dense_output=True)
    A = sol.t; Y = sol.y[0]
    below = Y[A < 0.98 * a_t]; above = Y[A > 1.02 * a_t]
    zc_below = int(np.sum(np.diff(np.sign(below)) != 0))
    zc_above = int(np.sum(np.diff(np.sign(above)) != 0))
    # WKB count of zero crossings above the turning point: (1/pi) * integral sqrt(q a^4 - a^2)/eps
    from scipy.integrate import quad
    wkb = quad(lambda a: np.sqrt(max(q * a ** 4 - a ** 2, 0)) / eps, 1.02 * a_t, 2.2 * a_t)[0] / np.pi
    k5[q] = (zc_below, zc_above, wkb)
ok5 = all(v[0] == 0 and v[1] >= 8 and abs(v[1] - v[2]) <= 2 for v in k5.values())
check("K5", ok5,
      f"real state of the de Sitter constraint (eps = {eps}): zero crossings below / above the turning point and the "
      f"WKB count above: " + "; ".join(f"q = {q}: {v[0]} / {v[1]} (WKB {v[2]:.1f})" for q, v in k5.items())
      + "; the turning point is a = 1/sqrt(q), the law's slice size (round 16): no oscillation, so no clock and "
        "no phase, below it; a clock and a phase above it",
      k5={str(k): list(map(float, v)) for k, v in k5.items()})

# ---------------------------------------------------------------- K6
a_s, q_s, t_s = sp.symbols("a q tau", positive=True)
p = sp.sqrt(q_s * a_s ** 4 - a_s ** 2)
dtau_da = sp.simplify(a_s / p)                       # from H = -(p^2 + a^2 - q a^4)/(2a): adot = -p/a
tau_of_a = sp.acosh(sp.sqrt(q_s) * a_s) / sp.sqrt(q_s)
ok_tau = sp.simplify(sp.diff(tau_of_a, a_s) - dtau_da) == 0
a_of_tau = sp.cosh(sp.sqrt(q_s) * t_s) / sp.sqrt(q_s)
ok_frw = sp.simplify(sp.diff(a_of_tau, t_s) ** 2 - (q_s * a_of_tau ** 2 - 1)) == 0
# numerically: fused shift chains reproduce both the growth (real k) and the matter turn (imaginary k)
qv, T, N = 1.0, 1.5, 400000
# round 16: forward run plus inverse run of the shift with k = sqrt(q), balanced at size 1/sqrt(q)
grow = ((1 + np.sqrt(qv) * T / N) ** N + (1 - np.sqrt(qv) * T / N) ** N) / (2 * np.sqrt(qv))
turn = np.linalg.matrix_power(np.eye(2) + (-1j) * Hm * (T / N), N)
err_grow = abs(grow - np.cosh(np.sqrt(qv) * T) / np.sqrt(qv))
err_turn = np.abs(turn - expm(-1j * Hm * T)).max()
check("K6", ok_tau and ok_frw and err_grow < 1e-5 and err_turn < 1e-5,
      f"the clock read off the oscillating state, tau = arccosh(sqrt(q) a)/sqrt(q), solves dtau/da = a/p exactly: "
      f"{ok_tau}; inverted, a = cosh(sqrt(q) tau)/sqrt(q) obeys the closed de Sitter law adot^2 = q a^2 - 1: {ok_frw}. "
      f"The fused shift with real k = sqrt(q) (forward + inverse run) gives that growth to {err_grow:.1e}; with "
      f"imaginary k it gives the matter turn e^(-i H_m tau) to {err_turn:.1e}")

# ---------------------------------------------------------------- P1, P2, P3
psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi)
B = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[0]
P_known = np.abs(B.conj().T @ psi) ** 2
thetas = np.linspace(0, 2 * np.pi, 360, endpoint=False)
P_avg = np.mean([np.abs(B.conj().T @ (np.exp(1j * th) * psi)) ** 2 for th in thetas], axis=0)
dP1 = np.abs(P_avg - P_known).max()
check("P1", dP1 < 1e-14,
      f"averaging over a completely unknown global clock phase changes every outcome probability by {dP1:.1e}: "
      f"a hidden global phase is exactly what A4 says drops out, so it cannot be the dice")

aa, bb = 0.8 * np.exp(0.3j), 0.6
I_det = np.abs(aa + bb) ** 2
I_avg = np.mean([np.abs(aa * np.exp(1j * th) + bb) ** 2 for th in thetas])
check("P2", abs(I_avg - (abs(aa) ** 2 + abs(bb) ** 2)) < 1e-12 and abs(I_det - I_avg) > 0.1,
      f"an unknown relative phase between two paths washes the interference out: {I_det:.4f} with the phase known, "
      f"{I_avg:.4f} averaged = |a|^2 + |b|^2. The weights |a|^2 are the norm that rotations keep; the averaging "
      f"supplies the mixture, not the squares (decoherence, not a derivation of the Born rule)")

# envariance: reduced state independent of the Schmidt phases; equal amplitudes -> swap symmetry; fine-graining
amps = np.array([np.sqrt(2 / 3), np.sqrt(1 / 3)])
red = []
for phis in [(0, 0), (0.7, -2.1), (3.0, 1.0)]:
    st = np.zeros((2, 2), complex)
    for kk in range(2):
        st[kk, kk] = amps[kk] * np.exp(1j * phis[kk])     # sum_k a_k e^{i phi_k} |k>|e_k>
    red.append(st @ st.conj().T)
spread = max(np.abs(r - red[0]).max() for r in red)
# fine-grain: |0>|e0> -> (|0>|e0a> + |0>|e0b>)/sqrt2 ; three equal-amplitude branches, two of them carry outcome 0
fine = np.zeros((2, 3)); fine[0, 0] = fine[0, 1] = np.sqrt(2 / 3) / np.sqrt(2); fine[1, 2] = np.sqrt(1 / 3)
branch_amps = fine[np.nonzero(fine)]
odds0 = 2 / 3 if np.allclose(branch_amps, branch_amps[0]) else None
check("P3", spread < 1e-15 and odds0 == 2 / 3,
      f"Zurek's envariance: the system's local state does not depend on the branch phases (spread {spread:.1e}); "
      f"after fine-graining sqrt(2/3)|0> + sqrt(1/3)|1> into three equal-amplitude branches "
      f"({', '.join(f'{x:.4f}' for x in branch_amps)}), symmetry gives outcome 0 the odds 2/3. "
      f"Unknowable phases do lead to the Born rule here, inside quantum mechanics (with Zurek's assumptions)")

# ---------------------------------------------------------------- P4
lam = np.linspace(0, 2 * np.pi, 200000, endpoint=False)


def E_local(al, be):
    A_ = np.sign(np.cos(2 * (al - lam)))       # each particle reads the shared clock phase lambda it left with
    B_ = -np.sign(np.cos(2 * (be - lam)))
    return np.mean(A_ * B_)


E_q = lambda al, be: -np.cos(2 * (al - be))
sets = (0, np.pi / 4, np.pi / 8, 3 * np.pi / 8)
def chsh(E):
    a1, a2, b1, b2 = sets
    return abs(E(a1, b1) - E(a1, b2) + E(a2, b1) + E(a2, b2))
S_loc, S_q = chsh(E_local), chsh(E_q)
# best a local clock can do over a grid of settings
best = 0
grid = np.linspace(0, np.pi / 2, 9)
for a1 in grid:
    for a2 in grid:
        for b1 in grid:
            for b2 in grid:
                v = abs(E_local(a1, b1) - E_local(a1, b2) + E_local(a2, b1) + E_local(a2, b2))
                best = max(best, v)
check("P4", S_loc <= 2 + 1e-9 and best <= 2 + 1e-3 and abs(S_q - 2 * np.sqrt(2)) < 1e-12,
      f"a clock phase carried by each particle from the source (local hidden variable): CHSH = {S_loc:.4f} at the "
      f"quantum-optimal settings, best over a 9^4 grid of settings {best:.4f} <= 2; quantum mechanics "
      f"{S_q:.4f} = 2 sqrt 2; loophole-free experiments (2015) measured above 2 (Hensen et al.: 2.42 +/- 0.20). "
      f"A local small clock is ruled out", S_local=float(S_loc), S_best=float(best), S_quantum=float(S_q))

# ---------------------------------------------------------------- P5
m, c, g, h, t, hb = sp.symbols("m c g h t hbar", positive=True)
phase_low = m * c ** 2 * t / hb
phase_high = m * c ** 2 * (1 + g * h / c ** 2) * t / hb   # weak-field redshift of the Compton clock
diff = sp.simplify(phase_high - phase_low)
check("P5", sp.simplify(diff - m * g * h * t / hb) == 0,
      f"two paths at heights 0 and h: the Compton clock phases mc^2 t/hbar each run at ~1e20 turns per second for an "
      f"atom, but they cancel; the difference is {diff}, the gravity-induced quantum phase measured for neutrons by "
      f"Colella-Overhauser-Werner (1975). Whether this counts as reading the Compton clock is disputed "
      f"(Muller-Peters-Chu 2010 vs Wolf et al. 2011)")

# ================================================================ Part D: a clock shared between the particles
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1.0])
paulis = [sx, sy, sz]


def max_chsh(rho):
    # Horodecki 1995: the best CHSH over all settings is 2 sqrt(m1 + m2), m1, m2 the two largest eigenvalues of T^T T
    T = np.array([[np.real(np.trace(rho @ np.kron(p, q))) for q in paulis] for p in paulis])
    ev = np.sort(np.linalg.eigvalsh(T.T @ T))[::-1]
    return 2 * np.sqrt(ev[0] + ev[1])


# D1: any shared classical clock = a mixture of deterministic local rules; the best of all 16 is 2
best_det = 0
for a0 in (1, -1):
    for a1 in (1, -1):
        for b0 in (1, -1):
            for b1 in (1, -1):
                best_det = max(best_det, abs(a0 * b0 - a0 * b1 + a1 * b0 + a1 * b1))
check("D1", best_det == 2,
      f"a clock shared by both particles but read by each with its own setting, however the clock runs (one phase, "
      f"many hands, global, relational): every such model is a mixture of the 16 deterministic answer rules, and the "
      f"best of them scores CHSH = {best_det}. Sharing a classical clock is still a local hidden variable (Bell 1964); "
      f"it adds nothing to the earlier wall (probability-entanglement/, checks F3, G1, I6)")

# D2: Toner-Bacon: two shared random clock hands on the sphere + 1 bit sent from one side reproduce the singlet exactly
N = 2_000_000
l1 = rng.normal(size=(N, 3)); l1 /= np.linalg.norm(l1, axis=1)[:, None]
l2 = rng.normal(size=(N, 3)); l2 /= np.linalg.norm(l2, axis=1)[:, None]


def unit(th):
    return np.array([np.sin(th), 0, np.cos(th)])


def tb(a, b, bit=True):
    sa1, sa2 = np.sign(l1 @ a), np.sign(l2 @ a)
    A_ = -sa1
    cbit = sa1 * sa2 if bit else np.ones(N)          # the one bit Alice sends
    B_ = np.sign((l1 + cbit[:, None] * l2) @ b)
    return np.mean(A_ * B_), np.mean(A_), np.mean(B_)


angs = [0.0, 0.4, 1.0, 1.7, 2.6, np.pi]
tb_err = max(abs(tb(unit(0.0), unit(x))[0] + np.cos(x)) for x in angs)
marg = max(max(abs(tb(unit(0.3), unit(x))[1]), abs(tb(unit(0.3), unit(x))[2])) for x in angs)
th4 = [0, np.pi / 2, np.pi / 4, -np.pi / 4]
Ets = lambda i, j, bit=True: tb(unit(th4[i]), unit(th4[2 + j]), bit)[0]
S_tb = abs(Ets(0, 0) + Ets(0, 1) + Ets(1, 0) - Ets(1, 1))
S_nobit = abs(Ets(0, 0, False) + Ets(0, 1, False) + Ets(1, 0, False) - Ets(1, 1, False))
check("D2", tb_err < 3e-3 and marg < 3e-3 and abs(S_tb - 2 * np.sqrt(2)) < 6e-3 and S_nobit <= 2 + 6e-3,
      f"Toner-Bacon (2003): two shared random clock hands plus ONE bit sent from one particle's side to the other's "
      f"give the singlet's correlation -cos(angle) to {tb_err:.1e} ({N:,} runs), fair coins on each side "
      f"(bias {marg:.1e}), and CHSH = {S_tb:.4f} (2 sqrt 2 = 2.8284). With the bit switched off the same clock "
      f"scores {S_nobit:.4f} <= 2. The price of a shared clock that reaches 2 sqrt 2 is a signal between the sides, "
      f"faster than light when they are far apart", S_TB=float(S_tb), S_nobit=float(S_nobit))

# D3: Pironio's floor: 2 sqrt 2 on CHSH needs sqrt 2 - 1 = 0.414 bits on average; mixing gives exactly that
p_bit = np.sqrt(2) - 1
S_mix = (1 - p_bit) * 2 + p_bit * 4                   # local optimum 2, one-bit PR-box strategy 4
check("D3", abs(S_mix - 2 * np.sqrt(2)) < 1e-12,
      f"the cheapest classical route: use the one-bit strategy that scores 4 in a fraction {p_bit:.4f} of runs and the "
      f"best local rule otherwise: CHSH = {S_mix:.4f} = 2 sqrt 2 at {p_bit:.4f} bits per pair on average. Pironio "
      f"(PRA 68, 062102, 2003) proved this is also the minimum. So a shared clock needs at least 0.41 bits of signal "
      f"per pair, and nothing in a clock picture says why it should stop at 2 sqrt 2 rather than go to 4")

# D4: a shared QUANTUM clock (Page-Wootters with the clock entangled to both particles)
#     is ordinary quantum mechanics: the reading-t state reaches 2 sqrt 2 only if the generator couples the particles
dC = 16
wC = 2 * np.pi / dC
HCq = -np.diag(np.arange(dC)) * wC
c0q = np.ones(dC) / np.sqrt(dC)
gq = 1                                                 # coupling on the clock's ladder
H_coupled = gq * wC * (np.kron(sx, sx) + np.kron(sy, sy)) / 2 + 2 * wC * np.eye(4)   # shifted: spectrum on the ladder
H_free = wC * (np.kron(sz, np.eye(2)) + np.kron(np.eye(2), sz)) / 2 + 2 * wC * np.eye(4)
ab0 = np.zeros(4, complex); ab0[1] = 1                 # |0 1>, a product state


def pw(Hsys):
    Psi = sum(np.kron(expm(-1j * HCq * t) @ c0q, expm(-1j * Hsys * t) @ ab0) for t in range(dC))
    return Psi / np.linalg.norm(Psi)


def reading(Psi, t):
    ct = expm(-1j * HCq * t) @ c0q
    v = np.kron(ct.conj(), np.eye(4)) @ Psi
    v /= np.linalg.norm(v)
    return np.outer(v, v.conj())


Hc_tot = np.kron(HCq, np.eye(4)) + np.kron(np.eye(dC), H_coupled)
Psi_c, Psi_f = pw(H_coupled), pw(H_free)
ts = np.linspace(0, dC - 1, 61)
S_c = [max_chsh(reading(Psi_c, t)) for t in ts]
S_f = [max_chsh(reading(Psi_f, t)) for t in ts]
check("D4", np.linalg.norm(Hc_tot @ Psi_c) < 1e-10 and abs(max(S_c) - 2 * np.sqrt(2)) < 1e-6 and max(S_f) <= 2 + 1e-9,
      f"one quantum clock shared by both particles (Page-Wootters, timeless state |H Psi| = "
      f"{np.linalg.norm(Hc_tot @ Psi_c):.1e}): read at time t, the pair's best CHSH rises from 2 to {max(S_c):.4f} "
      f"= 2 sqrt 2 when the clock's partner generator couples the two particles, and never passes {max(S_f):.4f} "
      f"when it turns each one separately. So the shared clock is not what crosses Bell's line; the joint turn is, and "
      f"the odds read off it are the Born rule, assumed. This is quantum mechanics, not a new mechanism. The law's "
      f"fused shift acting on the joint state behaves the same way: product turn, CHSH <= 2; coupled turn, 2 sqrt 2",
      S_coupled_max=float(max(S_c)), S_free_max=float(max(S_f)), S_curve=[float(x) for x in S_c],
      ts=[float(x) for x in ts])

# ================================================================ Part T: helical (twisted) time
def Jz_spin(j):
    ms = np.arange(j, -j - 1, -1)
    return np.diag(ms)


def Jx_spin(j):
    ms = np.arange(j, -j - 1, -1); dim = len(ms)
    Jp = np.zeros((dim, dim))
    for k in range(1, dim):
        m = ms[k]
        Jp[k - 1, k] = np.sqrt(j * (j + 1) - m * (m + 1))
    return (Jp + Jp.T) / 2


# T1: space turning steadily as time runs (helical time, pitch Omega): the rotating-frame Hamiltonian H - Omega J_z;
#     after a stretch of time with Omega T = theta, a spin-m component has turned by e^{-i m theta}
th_scan = np.linspace(0, 4 * np.pi, 401)
phase_half = np.angle(np.exp(-1j * 0.5 * th_scan))
f_half = {name: complex(expm(-1j * th * Jz_spin(0.5))[0, 0]) for name, th in
          [("quarter", np.pi / 2), ("half", np.pi), ("full", 2 * np.pi), ("double", 4 * np.pi)]}
f_one = {name: complex(expm(-1j * th * Jz_spin(1.0))[0, 0]) for name, th in
         [("quarter", np.pi / 2), ("half", np.pi), ("full", 2 * np.pi)]}
ok_t1 = (abs(f_half["half"] - (-1j)) < 1e-12 and abs(f_half["full"] + 1) < 1e-12 and abs(f_half["double"] - 1) < 1e-12
         and abs(f_one["quarter"] - (-1j)) < 1e-12 and abs(f_one["half"] + 1) < 1e-12)
check("T1", ok_t1,
      f"helical time, space turning by theta per loop: a spin-1/2 component picks up e^(-i theta/2): quarter turn "
      f"{f_half['quarter']:.3f}, half turn {f_half['half']:.0f}, full turn {f_half['full'].real:+.0f}, two turns "
      f"{f_half['double'].real:+.0f}; a spin-1 (m = 1) component picks up -i already at a quarter turn and -1 at a half "
      f"turn. The 'i once round' of the half-turn space (klein-universe/, U7) is the spin-1/2 case of this; the sign "
      f"at a full turn was measured with neutrons (Rauch et al. 1975; Werner et al. 1975). Which twist gives i "
      f"depends on the particle's spin",
      f_half={k: [v.real, v.imag] for k, v in f_half.items()}, f_one={k: [v.real, v.imag] for k, v in f_one.items()})

# T2: can the twist BE quantum mechanics' i?  The half-turn J = R(pi) squares to (-1)^(2j)
sq = {j: float(np.real((expm(-1j * np.pi * Jz_spin(j)) @ expm(-1j * np.pi * Jz_spin(j)))[0, 0])) for j in (0.5, 1.0, 1.5, 2.0)}
Rz = expm(-1j * np.pi * Jz_spin(0.5))
comm_axis = np.linalg.norm(Rz @ sz - sz @ Rz)
comm_cross = np.linalg.norm(Rz @ sx - sx @ Rz)
check("T2", sq[0.5] == -1 and sq[1.5] == -1 and sq[1.0] == 1 and sq[2.0] == 1 and comm_axis < 1e-12 and comm_cross > 1,
      f"a half-turn of space squares to (-1)^(2j): {', '.join(f'spin {j}: {v:+.0f}' for j, v in sq.items())}. So it "
      f"is a quarter turn J with J^2 = -1 only for fermions; for photons and other integer spins it squares to +1. "
      f"It also singles out the twist's axis: it commutes with turns about that axis ({comm_axis:.0e}) but not across "
      f"it ({comm_cross:.1f}). Quantum mechanics' i must square to -1 for every particle and commute with every "
      f"observable (Stueckelberg 1960). So a twist of space cannot be the i; at most it adds a spin-dependent phase "
      f"on top of it (Hestenes reads the Dirac i as the spin plane for this reason, and meets the same axis problem)")

# T3: closed helical time (mapping torus): period T, twist theta. Spin-m modes have energies (2 pi n + m theta)/T
Tper = 1.0
ns = np.arange(-3, 4)
lev = lambda m, th: np.sort((2 * np.pi * ns + m * th) / Tper)
flow_full_int = np.allclose(np.sort(lev(1, 2 * np.pi))[:-1], np.sort(lev(1, 0))[1:])
anti = np.sort((2 * np.pi * (ns + 0.5)) / Tper)
flow_full_half = np.allclose(lev(0.5, 2 * np.pi), anti)
flow_4pi_half = np.allclose(np.sort(lev(0.5, 4 * np.pi))[:-1], np.sort(lev(0.5, 0))[1:])
check("T3", flow_full_int and flow_full_half and flow_4pi_half,
      f"closed helical time (time a circle of period T, space turned by theta once round: a mapping torus): spin-m "
      f"levels sit at (2 pi n + m theta)/T. A full turn maps integer-spin levels onto themselves ({flow_full_int}) but "
      f"moves spin-1/2 levels by half a spacing, onto the antiperiodic ladder fermions use on a thermal circle "
      f"({flow_full_half}); spin-1/2 needs two full turns to come back ({flow_4pi_half}). The twist shifts energies "
      f"(a rotation chemical potential, e^(-beta(H - Omega J)) after Wick rotation); it does not create the phase")

# T4: twists act locally (or rotate both particles together): Tsirelson's 2 sqrt 2 does not move at any pitch
singlet = np.array([0, 1, -1, 0], complex) / np.sqrt(2)
rho_s = np.outer(singlet, singlet.conj())
S_tw, fid_joint = [], []
for th in np.linspace(0, 4 * np.pi, 25):
    R = expm(-1j * th * Jx_spin(0.5) * 0.7 - 1j * th * Jz_spin(0.5) * 0.3)
    U = np.kron(R, expm(-1j * 1.3 * th * Jz_spin(0.5)))
    S_tw.append(max_chsh(U @ rho_s @ U.conj().T))
    fid_joint.append(abs(np.vdot(singlet, np.kron(R, R) @ singlet)))
check("T4", max(abs(np.array(S_tw) - 2 * np.sqrt(2))) < 1e-12 and min(fid_joint) > 1 - 1e-12,
      f"twisting each particle by its own amount, at 25 pitches over two full turns: the singlet's best CHSH stays "
      f"{min(S_tw):.6f} to {max(S_tw):.6f}; turning both together leaves the singlet unchanged (overlap >= "
      f"{min(fid_joint):.12f}). A twist is a turn, turns keep lengths, so no pitch moves the Born odds or the "
      f"2 sqrt 2 ceiling; a twist angle shared as a hidden variable is a classical clock (D1: <= 2)")

# T5: the law's own helix: a complex dial k = s + i w makes the fused shift a log-spiral (grow while turning).
s_, w_ = sp.symbols("s omega", real=True)
k2 = sp.expand((s_ + sp.I * w_) ** 2)
im_k2 = sp.im(k2)
check("T5", sp.simplify(im_k2 - 2 * s_ * w_) == 0,
      f"the law's fused shift with a complex dial k = s + i omega runs as e^(k tau): it grows at rate s while it turns "
      f"at rate omega, a log-spiral (a screw in scale and phase) of pitch omega/s. Its k^2 = {k2} has imaginary part "
      f"{im_k2}, so the sign table k^2 = -q_v q_s with real q's allows only s = 0 (pure turn: the quantum phase) or "
      f"omega = 0 (pure growth: de Sitter). In the law as it stands a single dial cannot be a helix; growth and turn "
      f"must live in different sectors (geometry vs matter). Ours, conditional on round 16's sign table")

# T6: how much helix the sky allows for space itself, and the internal (Compton) helix
H0 = 67.4e3 / 3.0857e22                     # s^-1
t0 = 13.8e9 * 3.156e7                       # s
vort = 4.7e-11                              # Saadeh et al. PRL 117, 131302 (2016), vector mode, 95%
Omega_max = vort * H0
turn_age = Omega_max * t0
hbar_MeVs = 6.582119569e-22
pitch_e = 0.51099895 / hbar_MeVs
pitch_p = 938.27209 / hbar_MeVs
check("T6", Omega_max < 1e-27 and turn_age < 1e-10,
      f"a helix of space itself (the cosmos turning steadily as time runs) is limited by the CMB's vorticity bound "
      f"(sigma_V/H)_0 < 4.7e-11 (Saadeh et al. 2016): Omega < {Omega_max:.1e} rad/s, under {turn_age:.1e} rad over the "
      f"universe's whole age. A pitch of one radian per Hubble time, the natural size if the law's sqrt(q) set it, is "
      f"excluded by a factor {1/vort:.0e}. The helix that does run fast is internal, the Compton clock: "
      f"{pitch_e:.2e} rad/s for an electron, {pitch_p:.2e} for a proton, {pitch_e/H0:.1e} and {pitch_p/H0:.1e} radians per "
      f"Hubble time: set by each particle's mass, so it singles out no universal pitch",
      Omega_max=float(Omega_max), pitch_e=float(pitch_e), pitch_p=float(pitch_p), H0=float(H0))

# T7: torsion as the law's order gap A o_k B - B o_k A = k [A, B] (torsion-dial/), for the generators of a helix
def affine_gen(kind):
    # 5x5 affine generators acting on (t, x, y, z, 1)
    G = np.zeros((5, 5))
    if kind == "Pt": G[0, 4] = 1
    if kind == "Px": G[1, 4] = 1
    if kind == "Pz": G[3, 4] = 1
    if kind == "Jz": G[1, 2], G[2, 1] = -1, 1
    return G


kk = 0.7
gap = lambda A, B: np.linalg.norm((A + B + kk * A @ B) - (B + A + kk * B @ A))
g_time = gap(affine_gen("Pt"), affine_gen("Jz"))       # helical time: step in time, turn in space
g_screw = gap(affine_gen("Pz"), affine_gen("Jz"))      # screw along its own axis
g_edge = gap(affine_gen("Px"), affine_gen("Jz"))       # step across the turning axis
check("T7", g_time < 1e-15 and g_screw < 1e-15 and g_edge > 0.5,
      f"torsion read as the law's order gap k[A, B]: a step in time and a turn of space commute, gap {g_time:.0e}; a "
      f"screw along its own axis, gap {g_screw:.0e}; a step across the turning axis, gap {g_edge:.2f}. So steady "
      f"helical time carries no torsion in the law's sense (as klein-universe/ found for shape twists); torsion "
      f"needs the twist concentrated on a line, as on a spinning cosmic string read as a timelike screw dislocation "
      f"(Tod 1994; Letelier 1995), or spin density in Einstein-Cartan (torsion-dial/, T8)")

print(f"\n{passed}/{total} pass")
results["summary"] = f"{passed}/{total} pass"
json.dump(results, open("results.json", "w"), indent=1, default=float)
