"""DRAFT checks: is probability entanglement? time from entanglement; small loops; whole-path model.
Run: python3 checks.py > checks_output.txt   (writes results.json for figures.py)"""
import json, math, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from lib import *

rng = np.random.default_rng(20261125)
R = {}
passed = failed = 0
def check(name, ok, detail=""):
    global passed, failed
    ok = bool(ok); passed += ok; failed += (not ok)
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""))

# ------------------------------------------------------------------ A. probability from entanglement
print("== A. probability from entanglement (known: purification, Zurek envariance; union by branch counting)")
p = rng.dirichlet(np.ones(5))
psi = np.zeros(25, complex)
for i in range(5): psi[i * 5 + i] = math.sqrt(p[i])
rA = ptrace(rho_of(psi), [5, 5], [0])
check("A1 purification: any distribution p is the reduced state of a pure entangled state; Shannon H(p) = entanglement entropy",
      np.allclose(np.diag(rA).real, p) and abs(vn(rA) + (p * np.log(p)).sum()) < 1e-12)

def env_state(c):
    n = len(c); v = np.zeros(n * n, complex)
    for i in range(n): v[i * n + i] = c[i]
    return v
X = np.array([[0, 1], [1, 0]])
for c, expect in [([1 / math.sqrt(2)] * 2, True), ([math.sqrt(.7), math.sqrt(.3)], False)]:
    v = env_state(np.array(c)); swapped = np.kron(X, np.eye(2)) @ v; undone = np.kron(np.eye(2), X) @ swapped
    ok = np.allclose(undone, v)
    check(f"A2 envariance: swap on system undone by swap on environment iff Schmidt weights equal (weights {np.round(np.abs(c)**2,2)} -> {ok})",
          ok == expect)
# fine graining: amplitudes sqrt(m_i/M) -> M equal branches
m = [3, 5]; M = sum(m)
v = np.zeros((2, M), complex)
col = 0
for i, mi in enumerate(m):
    for j in range(mi): v[i, col] = 1 / math.sqrt(M); col += 1
check("A3 fine-graining: amplitudes sqrt(3/8), sqrt(5/8) become 8 equal branches; counting gives 3/8, 5/8",
      np.allclose(np.abs(v) ** 2, np.where(np.abs(v) > 0, 1 / M, 0)) and np.allclose((np.abs(v) ** 2).sum(1), np.array(m) / M))
ok = True
for a, Ma, b, Mb in [(3, 8, 2, 5), (1, 3, 1, 3), (7, 10, 4, 9), (5, 6, 1, 7)]:
    count = sum(1 for i in range(Ma) for j in range(Mb) if i < a or j < b)
    ok &= Fr(count, Ma * Mb) == Fr(a, Ma) + Fr(b, Mb) - Fr(a, Ma) * Fr(b, Mb)
check("A4 union rule by branch counting: independent systems, equal branches on an Ma x Mb grid, 'A or B' = p o_{-1} q exactly", ok)

# ------------------------------------------------------------------ B. union rule <-> no correlation
print("== B. the union rule is the signature of NO correlation (ours, stated in our terms; elementary)")
P = random_projector(rng, 2, 1); Q = random_projector(rng, 3, 1)
PA = np.kron(P, np.eye(3)); QB = np.kron(np.eye(2), Q)
J = PA + QB - PA @ QB
check("B1 commuting projectors: P v Q = P o_{-1} Q at operator level (P + Q - PQ is a projector)", np.allclose(J @ J, J) and np.allclose(J, J.conj().T))
ra, rb = rho_of(haar_state(rng, 2)), rho_of(haar_state(rng, 3))
rho = np.kron(ra, rb)
pa, qb, pj = [np.trace(rho @ O).real for O in (PA, QB, J)]
check("B2 product state: Prob(a or b) = p o_{-1} q", abs(pj - circ(pa, qb, -1)) < 1e-12)
pts = []; bound_ok = True
for trial in range(4000):
    kind = trial % 3
    if kind == 0: psi = haar_state(rng, 4); rho = rho_of(psi)
    elif kind == 1:
        w = rng.dirichlet(np.ones(4)); rho = sum(w[i] * np.kron(rho_of(haar_state(rng, 2)), rho_of(haar_state(rng, 2))) for i in range(4))
    else:
        w = rng.dirichlet(np.ones(2) * .5); rho = w[0] * rho_of(haar_state(rng, 4)) + w[1] * np.eye(4) / 4
    I = mutual_info(rho, 2, 2)
    best = 0
    for _ in range(12):
        P1 = np.kron(random_projector(rng, 2, 1), np.eye(2)); Q1 = np.kron(np.eye(2), random_projector(rng, 2, 1))
        pa, qb, pab = [np.trace(rho @ O).real for O in (P1, Q1, P1 @ Q1)]
        pj = np.trace(rho @ (P1 + Q1 - P1 @ Q1)).real
        d = pj - circ(pa, qb, -1)
        bound_ok &= abs(d + (pab - pa * qb)) < 1e-12 and abs(d) <= math.sqrt(2 * I) + 1e-12
        best = max(best, abs(d))
    pts.append((I, best, negativity(rho, 2, 2) > 1e-9, kind))
check("B3 defect = Prob(a or b) - p o_{-1} q = -covariance, and |defect| <= sqrt(2 I(A:B)) (Wolf-Verstraete-Hastings-Cirac 2008) in 48,000 tests", bound_ok)
R["B3"] = pts
zero_I = [b for I, b, e, k in pts if I < 1e-10]
check("B4 union rule holds for all tested local yes/no pairs only when I(A:B)=0 (all states with I>1e-3 show a defect)",
      all(b > 1e-6 for I, b, e, k in pts if I > 1e-3))
err = []
for _ in range(300):
    P1 = random_projector(rng, 3, 1); Q1 = random_projector(rng, 3, 1)
    Jn = P1 + Q1 - P1 @ Q1
    u, s, _ = np.linalg.svd(np.hstack([P1, Q1])); r = (s > 1e-9).sum(); Pj = u[:, :r] @ u[:, :r].conj().T
    err.append((np.linalg.norm(P1 @ Q1 - Q1 @ P1), np.linalg.norm(Jn - Pj)))
err = np.array(err)
Pc = np.diag([1, 0, 0]).astype(complex); Qc = np.diag([0, 1, 0]).astype(complex)
u, s_, _ = np.linalg.svd(np.hstack([Pc, Qc])); Pjc = u[:, :2] @ u[:, :2].conj().T
check("B5 one system, non-commuting tests: P + Q - PQ is not the 'or' projector whenever P, Q don't commute; exact when they do",
      np.all(err[err[:, 0] > 1e-3, 1] > 1e-4) and np.allclose(Pc + Qc - Pc @ Qc, Pjc), f"corr(commutator, error) {np.corrcoef(err[:,0],err[:,1])[0,1]:.2f}")
R["B5"] = err.tolist()

# ------------------------------------------------------------------ C. typicality
print("== C. typicality: small parts of a big entangled whole look random AND independent (known: Popescu-Short-Winter 2006, Page 1993)")
typ = []
for dR in [2, 4, 8, 16, 32, 64, 128, 256]:
    Is, Ds, Ts = [], [], []
    for _ in range(60):
        psi = haar_state(rng, 4 * dR); rho = ptrace(rho_of(psi), [2, 2, dR], [0, 1])
        Is.append(mutual_info(rho, 2, 2))
        ra = ptrace(rho, [2, 2], [0]); Ts.append(np.abs(np.linalg.eigvalsh(rho - np.eye(4) / 4)).sum())
        best = 0
        for _ in range(6):
            P1 = np.kron(random_projector(rng, 2, 1), np.eye(2)); Q1 = np.kron(np.eye(2), random_projector(rng, 2, 1))
            best = max(best, abs(np.trace(rho @ (P1 @ Q1)).real - np.trace(rho @ P1).real * np.trace(rho @ Q1).real))
        Ds.append(best)
    typ.append((dR, float(np.mean(Is)), float(np.mean(Ds)), float(np.mean(Ts))))
for dR, I, D, T in typ: print(f"   d_env={dR:4d}  <I(A:B)>={I:.5f} nats  <max union defect>={D:.4f}  <||rho_AB - I/4||_1>={T:.4f}  PSW bound sqrt(4/d)={math.sqrt(4/dR):.4f}")
check("C1 as the environment grows, I(A:B) and the union-rule defect fall toward 0 and rho_AB -> I/4 within the Popescu-Short-Winter bound",
      typ[-1][1] < typ[0][1] / 50 and typ[-1][2] < typ[0][2] / 5 and all(T <= math.sqrt(4 / d) for d, I, D, T in typ if d >= 8))
slope = np.polyfit(np.log([t[0] for t in typ[2:]]), np.log([t[1] for t in typ[2:]]), 1)[0]
check("C2 mutual information of two small parts falls like 1/d_env (Page-type scaling; monogamy: entanglement with the environment uses it up)", -1.2 < slope < -0.8, f"slope {slope:.2f}")
R["C"] = typ

# ------------------------------------------------------------------ D. time from entanglement
print("== D. time from entanglement (known: Page-Wootters 1983; 't Hooft cogwheel)")
N = 24; Tclk = 2 * math.pi
theta = sp.symbols('theta')
def history(psi0, H, N=N):
    """Page-Wootters history state: sum_t |t>_C (x) e^{-iHt}|psi0>, t = 0..N-1 in units of T/N (here H has period 2pi)"""
    w, V = np.linalg.eigh(H); out = []
    for t in range(N):
        U = V @ np.diag(np.exp(-1j * w * Tclk * t / N)) @ V.conj().T; out.append(U @ psi0)
    return np.array(out).reshape(-1) / math.sqrt(N)
Hs = np.diag([0.0, 1.0])     # qubit, level spacing commensurate with the clock
ok = True
for ang in [0.3, 1.1, 2.5]:
    psi0 = np.array([math.cos(ang / 2), math.sin(ang / 2)], complex); Psi = history(psi0, Hs)
    for t in [0, 5, 17]:
        cond = Psi.reshape(N, 2)[t]; cond = cond / np.linalg.norm(cond)
        target = np.diag(np.exp(-1j * np.diag(Hs) * Tclk * t / N)) @ psi0
        ok &= abs(abs(np.vdot(cond, target)) - 1) < 1e-12
# constraint: (H_C + H_S)|Psi> = 0 with H_C the clock's generator
Hc_shift = np.roll(np.eye(N), 1, axis=0)  # |t> -> |t+1>
w_c = -np.angle(np.linalg.eigvals(Hc_shift))
check("D1 Page-Wootters: conditioning the static history state on clock reading t gives e^{-iHt}|psi0> exactly", ok)
ents = []
for ang in np.linspace(0, math.pi, 25):
    psi0 = np.array([math.cos(ang / 2), math.sin(ang / 2)], complex); Psi = history(psi0, Hs)
    rS = ptrace(rho_of(Psi), [N, 2], [1]); rC = ptrace(rho_of(Psi), [N, 2], [0])
    ents.append((ang, vn(rC, 2), float(np.var(np.diag(Hs) @ (np.abs(psi0) ** 2) * 0 + np.abs(psi0) ** 2 @ np.diag(Hs) ** 2 - (np.abs(psi0) ** 2 @ np.diag(Hs)) ** 2))))
check("D2 no change, no time: clock-system entanglement is 0 exactly for stationary states and >0 whenever the system changes",
      ents[0][1] < 1e-9 and ents[-1][1] < 1e-9 and all(e[1] > 1e-3 for e in ents[1:-1]))
check("D3 the entanglement is largest (1 bit) when the qubit changes most (equal superposition)", abs(max(e[1] for e in ents) - 1) < 1e-9 and abs(ents[12][1] - 1) < 1e-9)
R["D"] = ents
n = 7; S = np.roll(np.eye(n), 1, axis=0)
ev = np.linalg.eigvals(S); Hcog = np.sort(np.mod(-np.angle(ev), 2 * math.pi))
w, V = np.linalg.eig(S); Hm = V @ np.diag(-1j * np.log(w) * 1j * 1j) @ np.linalg.inv(V)
Ucheck = V @ np.diag(np.exp(1j * np.angle(w))) @ np.linalg.inv(V)
check("D4 't Hooft cogwheel: a deterministic n-step loop (clock mod n) IS a quantum system with equally spaced energies 2*pi*k/n",
      np.allclose(Hcog, 2 * math.pi * np.arange(n) / n) and np.allclose(Ucheck, S))

# ------------------------------------------------------------------ E. small loops -> probability
print("== E. a fast small loop read coarsely gives probability (known: Poincare, Hopf, Keller: method of arbitrary functions)")
mafs = []
for om in [0.5, 1, 2, 4, 8, 16, 32]:
    T = rng.normal(5.0, 0.6, 400000)          # smooth spread of reading times
    ph = np.mod(om * T, 2 * math.pi)
    mafs.append((om, float(abs((ph < 2 * math.pi / 3).mean() - 1 / 3))))
check("E1 fast loop: Prob(phase in a 1/3 sector) -> 1/3 as the loop speeds up, whatever the (smooth) spread of reading times",
      mafs[-1][1] < 0.004 and mafs[0][1] > 0.05, " ".join(f"{o}:{e:.3f}" for o, e in mafs))
R["E1"] = mafs
lastdig = []
for jit in [0.3, 3, 30, 300]:
    ticks = np.floor(rng.normal(1e6, jit, 200000)).astype(np.int64)
    h = np.bincount(ticks % 10, minlength=10) / len(ticks); lastdig.append((jit, h.tolist()))
check("E2 computer clock: last digit of the tick counter is uniform when the reading time is spread over many ticks, not when it is sharp",
      max(abs(np.array(lastdig[-1][1]) - .1)) < .005 and max(abs(np.array(lastdig[0][1]) - .1)) > .2)
R["E2"] = lastdig
def union_defect(r, times):
    a = np.mod(times, 1) < 0.3; b = np.mod(r * times, 1) < 0.4
    return float((a | b).mean() - circ(a.mean(), b.mean(), -1))
gold = (1 + 5 ** .5) / 2
Ts = rng.uniform(0, 1, 400000) * 2000
dg, d2 = union_defect(gold, Ts), union_defect(2.0, Ts)
check("E3 two independent loops (incommensurate speeds) combine by the union rule; locked loops (speed ratio 2) do not",
      abs(dg) < 0.003 and abs(d2) > 0.05, f"golden {dg:.4f}, ratio-2 {d2:.4f}")
R["E3"] = {"golden": dg, "two": d2}
# finite o_1 clock
pr = 1019
def order2(p):
    x, k = 2, 1
    while x != 1: x = x * 2 % p; k += 1
    return k
orb = [(pow(2, k, pr) - 1) % pr for k in range(1, pr)]
acc = 0; ok = True
for k in range(1, 40):
    acc = (acc + 1 + acc * 1) % pr; ok &= acc == (pow(2, k, pr) - 1) % pr
check("E4 finite o_1 clock: 1 o 1 o ... o 1 (k times) = 2^k - 1 mod p", ok)
check(f"E5 for p = {pr} (2 a primitive root) the o_1 clock visits every element except the wall -1 once per loop", order2(pr) == pr - 1 and len(set(orb)) == pr - 1 and (pr - 1) not in orb)
hist = np.bincount(np.array(orb) % 10, minlength=10) / len(orb)
check("E6 its last digit is close to uniform", max(abs(hist - .1)) < .02, np.round(hist, 3).tolist())
R["E6"] = {"p": pr, "orbit": orb[:200], "lastdigit": hist.tolist()}
# Poisson law for hitting small windows: the o_1 clock with step a (x -> (1+a)(1+x) - 1 mod p) vs rotation
PM = 2 ** 31 - 1; GM = 16807          # 16807 is a primitive root mod 2^31 - 1 (Lehmer / MINSTD), so step a = 16806
eps = 0.002; W = int(eps * PM)
def hit_clock(x0, lo=0, hi=W, maxn=20000):
    y = (x0 + 1) % PM                 # work with y = 1 + x, where the o_1 step is multiplication by 1 + a
    for n in range(1, maxn):
        y = y * GM % PM; x = (y - 1) % PM
        if lo <= x < hi: return n
    return maxn
starts = rng.integers(1, PM - 2, 3000)
hd = np.array([hit_clock(int(s)) for s in starts])
alpha = math.sqrt(2) - 1
def hit_rot(x0, maxn=20000):
    ns = np.arange(1, maxn); xs = np.mod(x0 + ns * alpha, 1); i = np.argmax(xs < eps); return ns[i]
hr = np.array([hit_rot(x) for x in rng.uniform(0, 1, 3000)])
def ks_exp(h, mean):
    h = np.sort(h); F = np.arange(1, len(h) + 1) / len(h); return float(np.max(np.abs(F - (1 - np.exp(-h / mean)))))
kd, kr = ks_exp(hd, 1 / eps), ks_exp(hr, 1 / eps)
check("E7 Poisson law: waiting time for a chaotic o_1 clock (step a = 16806 mod 2^31 - 1) to enter a small window is exponential",
      kd < 0.04 and abs(hd.mean() * eps - 1) < 0.1, f"KS {kd:.3f}, mean*eps {hd.mean()*eps:.2f}")
check("E8 an irrational rotation (non-mixing loop) does NOT give exponential waiting (three-gap structure)", kr > 0.15, f"KS {kr:.3f}")
R["E7"] = {"clock": hd.tolist(), "rotation": hr.tolist(), "eps": eps}
tgrid = [100, 300, 600]
h1 = hd[:2000]
h2 = np.array([hit_clock(int(s), PM // 2, PM // 2 + W // 2) for s in starts[:2000]])
h12 = np.minimum(h1, h2)
devs = [abs((h12 <= t).mean() - circ((h1 <= t).mean(), (h2 <= t).mean(), -1)) for t in tgrid]
check("E9 two separate windows of the chaotic clock: 'hit either by time t' = (hit 1) o_{-1} (hit 2) within sampling error", max(devs) < 0.03, np.round(devs, 4).tolist())
R["E9"] = {"t": tgrid, "h1": [float((h1 <= t).mean()) for t in tgrid], "h2": [float((h2 <= t).mean()) for t in tgrid], "h12": [float((h12 <= t).mean()) for t in tgrid]}

# ------------------------------------------------------------------ F. Born weights from a loop, and the Bell wall
print("== F. can a loop give Born weights |a|^2? (known: Kochen-Specker 1967 qubit model; Bell 1964; Gleason 1957)")
phis = rng.uniform(-math.pi / 2, math.pi / 2, 400000)
wts = np.cos(phis)            # loop whose dwell time is proportional to cos(phase)
ok = True; born = []
for th in np.linspace(0, math.pi, 13):
    pplus = (wts * (np.cos(phis - th) > 0)).sum() / wts.sum(); born.append((th, float(pplus)))
    ok &= abs(pplus - math.cos(th / 2) ** 2) < 0.004
check("F1 one qubit: a loop with dwell density cos(phase) read by 'which half-circle' reproduces the Born rule cos^2(theta/2)", ok)
R["F1"] = born
a_, b_ = 0.6, 0.8 * np.exp(0.4j)
ph = rng.uniform(0, 2 * math.pi, 200000)
avg = np.mean(np.abs(a_ + b_ * np.exp(1j * ph)) ** 2)
check("F2 averaging a fast relative-phase loop removes interference: <|a + b e^{i phi}|^2> = |a|^2 + |b|^2 (the square itself is not explained)",
      abs(avg - (abs(a_) ** 2 + abs(b_) ** 2)) < 0.005)
best = 0
for sA in itertools.product([1, -1], repeat=2):
    for sB in itertools.product([1, -1], repeat=2):
        best = max(best, abs(sA[0] * sB[0] + sA[0] * sB[1] + sA[1] * sB[0] - sA[1] * sB[1]))
check("F3 Bell: every local loop/hidden-variable model has CHSH <= 2 (max over all 16 deterministic strategies = 2)", best == 2)

# ------------------------------------------------------------------ G. whole-path model
print("== G. a whole-path model: outcomes from a loop phase whose weight depends on both ends of the path")
def E_shared(D):            # measurement-independent shared phase, A=sgn cos(phi-a), B=-sgn cos(phi-b)
    D = abs((D + math.pi) % (2 * math.pi) - math.pi); return -(1 - 2 * D / math.pi)
def chsh(E, a=0, a2=math.pi / 2, b=math.pi / 4, b2=-math.pi / 4):
    return abs(E(a - b) + E(a - b2) + E(a2 - b) - E(a2 - b2))
Eq = lambda D: -math.cos(D)
grid = np.linspace(-math.pi, math.pi, 49)
bs = max(chsh(E_shared, a, a2, b, b2) for a in [0] for a2 in grid for b in grid for b2 in grid)
check("G1 shared loop phase, settings independent of the phase: best CHSH = 2 (straight-line correlations)", abs(bs - 2) < 1e-9, f"{bs:.4f}")
def wp_model(a, b, n=200000, target=lambda D: math.sin(D / 2) ** 2):
    """phase phi drawn with weight depending on BOTH settings: mass target(D) on the arcs where the two readings agree"""
    phi = rng.uniform(0, 2 * math.pi, 4 * n)
    A = np.sign(np.cos(phi - a)); B = -np.sign(np.cos(phi - b)); agree = A == B
    D = abs((a - b + math.pi) % (2 * math.pi) - math.pi)
    pa = target(D)
    fa = agree.mean()
    wt = np.where(agree, pa / max(fa, 1e-12), (1 - pa) / max(1 - fa, 1e-12))
    keep = rng.uniform(0, wt.max(), 4 * n) < wt
    return A[keep][:n], B[keep][:n], phi[keep][:n]
ok = True; Ecurve = []
for D in np.linspace(0.1, math.pi - 0.1, 9):
    A, B, _ = wp_model(0.3, 0.3 - D); Ecurve.append((D, float((A * B).mean())))
    ok &= abs((A * B).mean() + math.cos(D)) < 0.01 and abs(A.mean()) < 0.01 and abs(B.mean()) < 0.01
check("G2 whole-path weight w(phi | a, b): local readings A(phi,a), B(phi,b) reproduce the singlet E = -cos(a-b) and fair marginals", ok)
R["G2"] = Ecurve
Ewp = lambda D: float((lambda r: (r[0] * r[1]).mean())(wp_model(0, -D, 100000)))
cw = chsh(Ewp)
check("G3 ... and reach CHSH = 2*sqrt(2) with locality kept; what is dropped is measurement independence", abs(cw - 2 * math.sqrt(2)) < 0.03, f"{cw:.3f}")
pr_target = lambda D: 0.0 if D < math.pi / 2 - 1e-9 else 1.0
Epr = lambda D: float((lambda r: (r[0] * r[1]).mean())(wp_model(0, -D, 100000, pr_target)))
cpr = abs(Epr(-math.pi / 4) + Epr(math.pi / 4) + Epr(math.pi / 4) - Epr(3 * math.pi / 4))
check("G4 the same freedom gives a PR box (CHSH 4): whole-path weighting alone does NOT explain Tsirelson's 2*sqrt(2)", cpr > 3.9, f"{cpr:.3f}")
def pathsum_E(a, b):
    """two branches |01>-|10>, each a loop phase composed along the path; weight = |sum of phasors|^2"""
    P = {}
    for A in (1, -1):
        for B in (1, -1):
            sa = np.array([math.cos(a / 2), math.sin(a / 2)]) if A == 1 else np.array([-math.sin(a / 2), math.cos(a / 2)])
            sb = np.array([math.cos(b / 2), math.sin(b / 2)]) if B == 1 else np.array([-math.sin(b / 2), math.cos(b / 2)])
            amp = (sa[0] * sb[1] - sa[1] * sb[0]) / math.sqrt(2); P[(A, B)] = amp ** 2
    return sum(A * B * v for (A, B), v in P.items()), sum(P.values())
okn = all(abs(pathsum_E(a, b)[1] - 1) < 1e-12 and abs(pathsum_E(a, b)[0] + math.cos(a - b)) < 1e-12 for a, b in rng.uniform(0, 6, (50, 2)))
mx = max(chsh(lambda D: pathsum_E(0, -D)[0], *rng.uniform(-math.pi, math.pi, 4)) for _ in range(20000))
check("G5 path-sum version (phases composed along two branches, weight = |sum|^2) reproduces the singlet; random search never exceeds 2*sqrt(2)",
      okn and mx <= 2 * math.sqrt(2) + 1e-9 and mx > 2.8, f"max {mx:.4f}")
# measurement dependence: mutual information between phi (binned) and settings
settings = [(0, math.pi / 4), (0, -math.pi / 4), (math.pi / 2, math.pi / 4), (math.pi / 2, -math.pi / 4)]
bins = 72; Hcond = 0; marg = np.zeros(bins)
for a, b in settings:
    _, _, ph = wp_model(a, b, 200000); h = np.histogram(np.mod(ph, 2 * math.pi), bins, (0, 2 * math.pi))[0].astype(float); h /= h.sum()
    marg += h / 4; Hcond += -(h[h > 0] * np.log2(h[h > 0])).sum() / 4
MI = -(marg[marg > 0] * np.log2(marg[marg > 0])).sum() - Hcond
check("G6 how much the loop phase 'knows' the settings at the CHSH angles (bits, 72 bins; includes small binning noise)", 0 < MI < 0.5, f"{MI:.3f} bits")
R["G6"] = MI
# same statistics from a Bohm-like nonlocal model: phi uniform, B also reads a
def bohm_like(a, b, n=200000):
    phi = rng.uniform(0, 2 * math.pi, n); A = np.sign(np.cos(phi - a))
    D = a - b; pagree = math.sin(D / 2) ** 2
    B = np.where(rng.uniform(size=n) < pagree, A, -A)  # B's reading uses a (nonlocal)
    return A, B
diff = 0
for a, b in settings:
    A1, B1, _ = wp_model(a, b, 200000); A2, B2 = bohm_like(a, b)
    for x, y in itertools.product((1, -1), repeat=2):
        diff = max(diff, abs(((A1 == x) & (B1 == y)).mean() - ((A2 == x) & (B2 == y)).mean()))
check("G7 a nonlocal (Bohm-like) model with settings-independent phase gives the same outcome table: same data, two stories", diff < 0.01, f"max diff {diff:.4f}")
# entropy-proof winding vs diffusing local phase
steps = 4000; walkers = 400
noise = rng.normal(0, 0.05, (walkers, steps))
base = np.linspace(0, 2 * math.pi * 3, steps)          # loop wound 3 times along the path
path = base + np.cumsum(noise, 1) - np.cumsum(noise, 1)[:, -1:] * np.linspace(0, 1, steps)  # noise, ends pinned (closed loop)
wind = np.round((path[:, -1] - path[:, 0]) / (2 * math.pi)).astype(int)
local_var = np.var(np.cumsum(noise, 1)[:, -1])
check("G8 'not susceptible to entropy': a closed loop's winding number (3) survives noise in every one of 400 runs while the open phase spreads",
      np.all(wind == 3) and local_var > 5, f"open-phase variance {local_var:.1f} rad^2")
psi = haar_state(rng, 4); Pq = random_projector(rng, 4, 2)
check("G9 a single global phase is unobservable (e^{i alpha}|psi> gives identical probabilities); only phase differences between paths count",
      abs(np.vdot(psi, Pq @ psi) - np.vdot(np.exp(.7j) * psi, Pq @ (np.exp(.7j) * psi))) < 1e-12)
t1, t2 = sp.symbols('t1 t2', real=True)
Ei = lambda t: sp.exp(sp.I * t) - 1
check("G10 composing path segments: loop phasors E(i t) = e^{it} - 1 compose by o_1 (height flow at imaginary height)",
      sp.simplify(sp.expand(Ei(t1) + Ei(t2) + Ei(t1) * Ei(t2) - Ei(t1 + t2))) == 0)

# ------------------------------------------------------------------ H. AdS3 / CFT: cross-ratio as a modular-time Poisson clock
print("== H. the cross-ratio is a probability: 'at least one tick' of a unit-rate clock running on A's modular time (ours, reading; ingredients known)")
ok = True; ok2 = True
for _ in range(2000):
    pts_ = np.sort(rng.uniform(0, 2 * math.pi, 4)); a1, a2, b1, b2 = pts_
    x = xr_circle(a1, a2, b1, b2); t = -math.log(1 - x)
    tz = abs(zeta(b2, a1, a2) - zeta(b1, a1, a2))
    ok &= abs(t - tz) < 1e-9 * max(1, t)
check("H1 t = -log(1 - x) equals the width of B measured in A's modular (Rindler) time coordinate zeta (2000 random arcs)", ok)
def kin_measure(a1, a2, b1, b2, n=600):
    u = np.linspace(a1, a2, n + 1); v = np.linspace(b1, b2, n + 1)
    uc = (u[1:] + u[:-1]) / 2; vc = (v[1:] + v[:-1]) / 2
    U, V = np.meshgrid(uc, vc)
    return float((1 / (4 * np.sin((V - U) / 2) ** 2)).sum() * (u[1] - u[0]) * (v[1] - v[0]))
errs = []
for _ in range(20):
    a1, a2, b1, b2 = np.sort(rng.uniform(0, 2 * math.pi, 4))
    if min(b1 - a2, 2 * math.pi - b2 + a1) < 0.3: continue
    x = xr_circle(a1, a2, b1, b2); errs.append(abs(kin_measure(a1, a2, b1, b2) + math.log(1 - x)))
check("H2 t is also the Crofton (kinematic-space) measure of geodesics joining A to B (Czech-Lamprou-McCandlish-Sully 2015)", max(errs) < 1e-4, f"max err {max(errs):.1e}")
def poisson_geodesics(delta=0.02):
    mass = math.pi / math.tan(delta / 2)
    n = rng.poisson(mass)
    u = rng.uniform(0, 2 * math.pi, n)
    r = rng.uniform(0, 1, n); c0 = 1 / math.tan(delta / 2)
    s = 2 * np.arctan(1 / (c0 - r * c0))            # shorter-arc length, density ~ 1/(4 sin^2(s/2)) on (delta, pi)
    return np.mod(u, 2 * math.pi), np.mod(u + s, 2 * math.pi)
def inarc(x, lo, hi): return (x >= lo) & (x <= hi)
A_ = (0.3, 1.4); B1 = (2.4, 3.3); B2 = (3.3, 4.6)
xs = {k: xr_circle(A_[0], A_[1], *B) for k, B in [("B1", B1), ("B2", B2), ("B12", (B1[0], B2[1]))]}
trials = 3000; links = {"B1": [], "B2": [], "B12": []}
for _ in range(trials):
    u, v = poisson_geodesics()
    for k, B in [("B1", B1), ("B2", B2), ("B12", (B1[0], B2[1]))]:
        L = (inarc(u, *A_) & inarc(v, *B)) | (inarc(v, *A_) & inarc(u, *B)); links[k].append(L.sum())
fr = {k: float((np.array(v) > 0).mean()) for k, v in links.items()}
se = max(math.sqrt(xs[k] * (1 - xs[k]) / trials) for k in xs)
check("H3 Monte Carlo: a random (Poisson) spray of geodesics links A and B with probability x (3000 universes)",
      all(abs(fr[k] - xs[k]) < 4 * se for k in xs), " ".join(f"{k}: sim {fr[k]:.3f} vs x {xs[k]:.3f}" for k in xs))
corr = np.corrcoef(links["B1"], links["B2"])[0, 1]
check("H4 gluing: links to B1 and to B2 are independent, so x(A,B1B2) = x(A,B1) o_{-1} x(A,B2) is the ordinary union rule",
      abs(corr) < 0.06 and abs(xs["B12"] - circ(xs["B1"], xs["B2"], -1)) < 1e-12, f"count correlation {corr:.3f}")
R["H3"] = {"A": A_, "B1": B1, "B2": B2, "x": xs, "sim": fr}
# free-fermion lattice: CMI of adjacent intervals = (c/3) t, c = 1
ff = []
for la, lb, lc in [(20, 10, 20), (30, 30, 30), (40, 15, 60), (25, 50, 25), (60, 20, 40), (50, 50, 80), (80, 30, 30)]:
    A = list(range(0, la)); B = list(range(la, la + lb)); C = list(range(la + lb, la + lb + lc))
    cmi = ff_entropy(A + B) + ff_entropy(B + C) - ff_entropy(B) - ff_entropy(A + B + C)
    t = math.log((la + lb) * (lb + lc) / (lb * (la + lb + lc)))
    ff.append((la, lb, lc, cmi, t / 3))
for r in ff: print(f"   lA,lB,lC={r[:3]}  CMI={r[3]:.4f}  t/3={r[4]:.4f}")
check("H5 real quantum state (free fermions, c = 1): conditional mutual information I(A:C|B) = t/3, so x = 1 - exp(-3 I/c)",
      all(abs(r[3] - r[4]) < 0.02 for r in ff), f"max |diff| {max(abs(r[3]-r[4]) for r in ff):.4f}")
R["H5"] = ff
eta1, eta2 = sp.symbols('eta1 eta2', real=True)
E = lambda t: sp.exp(t) - 1
check("H6 odds X = x/(1-x) = E(t) = e^t - 1 has the form of a Doppler shift z = e^eta - 1; adding rapidities composes shifts by o_1",
      sp.simplify(E(eta1) + E(eta2) + E(eta1) * E(eta2) - E(eta1 + eta2)) == 0)
ok = True
for _ in range(200):
    a1, a2, b1, b2 = np.sort(rng.uniform(0, 2 * math.pi, 4)); s = rng.uniform(-2, 2)
    def flow(u):  # move u by s in A's modular time (zeta -> zeta + s), staying in the complement
        z = zeta(u, a1, a2) + s
        # invert: sin((u-a1)/2)/sin((u-a2)/2) = +/- e^z on complement; solve numerically by bisection on (a2, a1+2pi)
        lo, hi = a2 + 1e-12, a1 + 2 * math.pi - 1e-12
        for _ in range(200):
            mid = (lo + hi) / 2
            if zeta(mid, a1, a2) < z: hi = mid   # zeta decreases from +inf to -inf across the complement
            else: lo = mid
        return (lo + hi) / 2
    nb1, nb2 = flow(b1), flow(b2)
    ok &= abs(xr_circle(a1, a2, nb1, nb2) - xr_circle(a1, a2, b1, b2)) < 1e-7
check("H7 A's modular flow (its own time) slides B along without changing x: x is a clock interval, not a clock reading", ok)


# ------------------------------------------------------------------ I. the helix (universal cover of the loop)
print("== I. the helix: the loop's universal cover (known: covering spaces, Aharonov-Bohm 1959, Berry 1984); our lift/flow reading")
th = np.cumsum(rng.normal(0.9, 0.6, 400))           # a phase path that keeps turning
z = np.exp(1j * th)
principal = np.angle(z)
lifted = np.unwrap(principal)
check("I1 rolling up: the helix height theta maps to the circle point e^{i theta}; unrolling (lifting) recovers theta up to one constant turn",
      np.allclose(lifted - th, (lifted - th)[0]) and abs(((lifted - th)[0] / (2 * math.pi)) - round((lifted - th)[0] / (2 * math.pi))) < 1e-9)
# o_1 on the circle and its lift: principal log is not additive, the lift along the path is
a_s, b_s = 2.5, 1.9
Ea, Eb = np.exp(1j * a_s) - 1, np.exp(1j * b_s) - 1
comp = Ea + Eb + Ea * Eb
l1 = lambda w: np.log(1 + w)          # lift l_1 (principal branch)
check("I2 the lift l_1 = log(1 + x) turns o_1 into + on the line, but on the circle its principal branch drops a turn (2.5 + 1.9 > pi)",
      abs(l1(comp) - 1j * (a_s + b_s)) > 6 and abs(l1(comp) + 2j * math.pi - 1j * (a_s + b_s)) < 1e-12)
check("I3 ... and the dropped piece is exactly 2*pi*i times the winding number: the helix is what keeps the count", abs((1j * (a_s + b_s) - l1(comp)) / (2j * math.pi) - 1) < 1e-12)
# finite circle (tan-law) clock: p+1 elements, cyclic; its index is the helix mod p+1
p_ = 19
def tadd(x, y):
    INF = None
    if x is INF: return INF if y == 0 else (-pow(y, -1, p_)) % p_
    if y is INF: return tadd(y, x)
    den = (1 - x * y) % p_
    if den == 0: return INF
    return (x + y) * pow(den, -1, p_) % p_
gens = []
for g in range(1, p_):
    x, seen = g, [0]
    while x != 0:
        seen.append(x); x = tadd(x, g)
        if len(seen) > p_ + 2: break
    if len(seen) == p_ + 1: gens.append(g)
check(f"I4 the finite tan-law clock over F_{p_} (p = 3 mod 4) is a loop of p+1 = {p_+1} steps; counting steps on the integers is its helix",
      len(gens) > 0, f"generators {gens[:6]}")
# Aharonov-Bohm: phase depends only on winding number times flux
flux = 0.7
def ab_phase(path_xy):
    ang = np.unwrap(np.arctan2(path_xy[:, 1], path_xy[:, 0])); return flux * (ang[-1] - ang[0]) / (2 * math.pi)
ok = True
for n in [1, 2, 3]:
    for _ in range(5):
        s = np.linspace(0, 2 * math.pi * n, 3000); r = 1 + 0.4 * np.sin(rng.integers(1, 7) * s + rng.uniform(0, 6))
        ok &= abs(ab_phase(np.c_[r * np.cos(s), r * np.sin(s)]) - n * flux) < 1e-9
check("I5 Aharonov-Bohm: a path winding n times around a hidden flux picks up n x flux, whatever its shape (a winding count noise can't change)", ok)
# Bell applies to any hidden variable space, helix included
best = 0
for _ in range(20000):
    h = rng.normal(0, 10)                       # helix height: phase AND winding
    f = lambda s, h=h: np.sign(np.sin(h * (1 + s) + math.floor(h / (2 * math.pi)) * s + 0.3))
    A0, A1, B0, B1 = f(0.0), f(1.3), f(2.1), f(-0.7)
    best = max(best, abs(A0 * B0 + A0 * B1 + A1 * B0 - A1 * B1))
check("I6 a shared helix (phase + winding) read locally, settings independent of it: still CHSH <= 2 (Bell holds for any hidden variable)", best <= 2)

json.dump(R, open("results.json", "w"))
print(f"\n{passed} passed, {failed} failed")
