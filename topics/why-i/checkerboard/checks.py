"""Checks for topics/why-i/checkerboard (Nov 25, 2026). Run: python3 checks.py  (numpy, scipy, sympy)
Each check tests an identity or a limit that the write-up claims."""
import numpy as np, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cb import circ, move_matrix, swap_matrix, evolve, kernel
from scipy.special import iv, comb
rng = np.random.default_rng(1125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
DIALS = [(-1.0, "k = -1 (telegraph)"), (-1j, "k = -i (checkerboard)")]

# ---------- A. the chain identity ----------
A1 = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)); B1 = rng.normal(size=(4, 4)); C1 = rng.normal(size=(4, 4))
I4 = np.eye(4)
ok = all(np.allclose(I4 + k*circ(A1, B1, k), (I4 + k*A1) @ (I4 + k*B1)) for k in [-1, -1j, 0.3, 2])
check("A1 1 + k(a o_k b) = (1 + k a)(1 + k b) for non-commuting matrices, at every dial", ok)
ok = True
for k in [-1, -1j]:
    lhs = circ(circ(A1, B1, k), C1, k)
    rhs = A1 + B1 + C1 + k*(A1 @ B1 + A1 @ C1 + B1 @ C1) + k**2*(A1 @ B1 @ C1)
    ok &= np.allclose(lhs, rhs)
check("A2 a o b o c = sum over the 7 non-empty ordered sub-chains, weight k^(size - 1): each once, order kept", ok)

# ---------- B. the checkerboard and the telegraph are that chain ----------
n, N, eps, mu = 19, 8, 0.1, 1.3
M = move_matrix(n); X = swap_matrix(n); Id = np.eye(2*n)
ok = True
for k, _ in DIALS:
    t = -eps*mu*X
    S = M @ (Id + k*t)
    tj = [np.linalg.matrix_power(np.linalg.inv(M), j) @ t @ np.linalg.matrix_power(M, j) for j in range(N)]
    chain = tj[N - 1]
    for j in range(N - 2, -1, -1): chain = circ(chain, tj[j], k)
    ok &= np.allclose(np.linalg.matrix_power(S, N), np.linalg.matrix_power(M, N) @ (Id + k*chain))
check("B1 [ours, exact at finite eps] N steps of (turn, move) = the free moves times 1 + k(t_N o_k ... o_k t_1), "
      "where t_j is the turn element seen from step j", ok)
# B2: brute force sum over turn-time subsets
N = 12
for k, lab in DIALS:
    w = -k*eps*mu
    amp = np.zeros((2*N + 1, 2), dtype=complex)
    for S_ in itertools.product([0, 1], repeat=N):
        d = 0; x = 0
        for turn in S_:
            if turn: d = 1 - d
            x += 1 if d == 0 else -1
        amp[x + N, d] += w**sum(S_)
    ev = evolve(N, eps, mu, k)[-1]; c0 = ev.shape[0]//2
    ok = np.allclose(amp, ev[c0 - N:c0 + N + 1])
    check(f"B2 {lab}: summing (turn weight)^R over all 2^12 choices of turn times gives the chain's amplitude at every site "
          "(every checkerboard path once)", ok)
# B3: Feynman-Hibbs / Jacobson-Schulman path counts
def count_formula(N, x, R, d_end):
    # the starting direction R counts as a virtual first right step, so a turn before the first move is counted
    r = (N + x)//2 + 1; l = N - r + 1
    if (N + x) % 2 or r < 1 or l < 0: return 0
    if d_end == 0:                      # start R, end R: even number of turns R = 2s
        if R % 2: return 0
        s = R//2
        if s == 0: return 1 if l == 0 else 0
        return int(comb(r - 1, s, exact=True)*comb(l - 1, s - 1, exact=True)) if r >= 1 and l >= 1 else 0
    else:                               # start R, end L: odd R = 2s + 1
        if R % 2 == 0: return 0
        s = (R - 1)//2
        return int(comb(r - 1, s, exact=True)*comb(l - 1, s, exact=True)) if r >= 1 and l >= 1 else 0
counts = {}
for S_ in itertools.product([0, 1], repeat=N):
    d = 0; x = 0
    for turn in S_:
        if turn: d = 1 - d
        x += 1 if d == 0 else -1
    counts[(x, sum(S_), d)] = counts.get((x, sum(S_), d), 0) + 1
ok = all(counts.get((x, R, d), 0) == count_formula(N, x, R, d) for x in range(-N, N + 1) for R in range(N + 1) for d in (0, 1))
check("B3 the number of paths with R turns ending at x is C(r-1, s) C(l-1, s-1) (end moving right, R = 2s) or "
      "C(r-1, s) C(l-1, s) (end moving left, R = 2s+1); r counts the right steps plus the starting direction, l the left steps (Feynman-Hibbs problem 2-6)", ok)
# B4: normalized telegraph = pure chain with turn weight a eps/(1 - a eps), times (1 - a eps)^N
a = 1.3
p = evolve(N, eps, a, -1, normalized=True)[-1]
q = evolve(N, eps/(1 - a*eps), a, -1)[-1]*(1 - a*eps)**N
check("B4 Kac's probability chain (continue 1 - a eps, turn a eps) is the k = -1 chain with turn weight a eps/(1 - a eps), times (1 - a eps)^N",
      np.allclose(p, q))
P = evolve(60, 0.05, a, -1, normalized=True)
check("B5 at k = -1 (normalized) every entry is a probability and they sum to 1 at every step: a true measure on paths",
      np.all(P.real >= -1e-15) and np.allclose(P.real.sum(axis=(1, 2)), 1) and np.allclose(P.imag, 0))
U = evolve(60, 0.05, 1.3, -1j, normalized=True)
Npure = evolve(60, 0.05, 1.3, -1j)
check("B6 at k = -i the step (cos m eps, i sin m eps) keeps the total |amplitude|^2 = 1 exactly (Meyer's quantum walk); "
      "the pure chain (1, i m eps) drifts by (1 + m^2 eps^2)^(N/2), which goes to 1 as eps -> 0",
      np.allclose((np.abs(U)**2).sum(axis=(1, 2)), 1) and np.isclose((np.abs(Npure[-1])**2).sum(), (1 + (1.3*.05)**2)**60))

# ---------- C. the continuum limit ----------
T, m_ = 2.0, 1.3
xs_test = np.array([-1.2, -0.5, 0.0, 0.4, 1.1])
for k, lab in DIALS:
    w = -k*m_
    errs = []
    for N in [40, 80, 160, 320]:
        e_ = T/N; ev = evolve(N, e_, m_, k)[-1]; c0 = ev.shape[0]//2
        GR, GL = kernel(xs_test, T, w)
        idx = np.round(xs_test/e_).astype(int)
        idx = idx + ((N - idx) % 2)                  # sites reachable after N steps have the parity of N
        xx = idx*e_; GR, GL = kernel(xx, T, w)
        lat = ev[c0 + idx]/(2*e_)                     # amplitude per unit length (every other site is reachable)
        errs.append(max(np.abs(lat[:, 0] - GR).max(), np.abs(lat[:, 1] - GL).max()))
    check(f"C1 {lab}: amplitude/(2 eps) converges to the continuum kernel inside the light cone; error halves with eps",
          all(errs[i + 1] < 0.6*errs[i] for i in range(3)), f"errors {[f'{e:.1e}' for e in errs]}")
# C2: the atom on the light cone
for k, lab in DIALS:
    w = -k*m_
    ev = evolve(400, T/400, m_, k)[-1]; c0 = ev.shape[0]//2
    check(f"C2 {lab}: the never-turn path carries amplitude 1 to x = t (pure chain); with Kac's normalization it carries (1 - a eps)^N -> e^(-a t)",
          np.isclose(ev[c0 + 400, 0], 1) and np.isclose(evolve(400, T/400, m_, -1, normalized=True)[-1][c0 + 400, 0].real, (1 - m_*T/400)**400) and abs((1 - m_*T/400)**400/np.exp(-m_*T) - 1) < 0.01)
# C3: the kernels solve the continuum equations (finite differences inside the cone)
h = 1e-4
for k, lab in DIALS:
    w = -k*m_
    x0, t0 = 0.3, 1.7
    GR, GL = kernel(x0, t0, w)
    dRt = (kernel(x0, t0 + h, w)[0] - kernel(x0, t0 - h, w)[0])/(2*h); dRx = (kernel(x0 + h, t0, w)[0] - kernel(x0 - h, t0, w)[0])/(2*h)
    dLt = (kernel(x0, t0 + h, w)[1] - kernel(x0, t0 - h, w)[1])/(2*h); dLx = (kernel(x0 + h, t0, w)[1] - kernel(x0 - h, t0, w)[1])/(2*h)
    check(f"C3 {lab}: the kernel solves d_t G_R = -d_x G_R + w G_L and d_t G_L = d_x G_L + w G_R (w = -k m)",
          abs(dRt + dRx - w*GL) < 1e-6 and abs(dLt - dLx - w*GR) < 1e-6)
# C4: second-order forms: telegraph equation (k=-1, with the damping put back) and Klein-Gordon (k=-i)
def u_tel(x, t, a):
    GR, GL = kernel(x, t, a); return (np.exp(-a*t)*(GR + GL)).real
x0, t0, a = 0.3, 1.7, 1.3
ut = (u_tel(x0, t0 + h, a) - u_tel(x0, t0 - h, a))/(2*h)
utt = (u_tel(x0, t0 + h, a) - 2*u_tel(x0, t0, a) + u_tel(x0, t0 - h, a))/h**2
uxx = (u_tel(x0 + h, t0, a) - 2*u_tel(x0, t0, a) + u_tel(x0 - h, t0, a))/h**2
check("C4 k = -1: Kac's density u = p_R + p_L solves the telegraph equation u_tt + 2a u_t = u_xx", abs(utt + 2*a*ut - uxx) < 1e-3, f"residual {utt + 2*a*ut - uxx:.1e}")
GL = lambda x, t: kernel(x, t, 1j*m_)[1]
gtt = (GL(x0, t0 + h) - 2*GL(x0, t0) + GL(x0, t0 - h))/h**2; gxx = (GL(x0 + h, t0) - 2*GL(x0, t0) + GL(x0 - h, t0))/h**2
check("C4b k = -i: each component solves the Klein-Gordon equation psi_tt - psi_xx + m^2 psi = 0 (the 1+1 Dirac equation squared)",
      abs(gtt - gxx + m_**2*GL(x0, t0)) < 1e-3)
# C5: momentum-space proof step: S(p)^(t/eps) -> exp(t(-i p sigma_z - k mu sigma_x)), uniformly bounded
sz = np.diag([1., -1.]); sx = np.array([[0, 1.], [1, 0]])
from scipy.linalg import expm
ok = True; bounds = []
for k, lab in DIALS:
    for pm in [0.0, 0.7, 3.0]:
        A = -1j*pm*sz - k*m_*sx
        es = []
        for N in [100, 1000, 10000]:
            e_ = T/N; Sp = np.diag([np.exp(-1j*pm*e_), np.exp(1j*pm*e_)]) @ (np.eye(2) - k*e_*m_*sx)
            es.append(np.abs(np.linalg.matrix_power(Sp, N) - expm(T*A)).max())
        ok &= es[1] < es[0]/5 and es[2] < es[1]/5
check("C5 proof step: in momentum space one step is 1 + eps A + O(eps^2), so N = t/eps steps converge to exp(t A) for every momentum (Lie-Trotter)", ok)
ps = np.linspace(-np.pi/0.01, np.pi/0.01, 2001); e_ = 0.01
nrm = max(np.linalg.norm(np.diag([np.exp(-1j*pm*e_), np.exp(1j*pm*e_)]) @ (np.eye(2) + 1j*e_*m_*sx), 2) for pm in ps)
check("C6 proof step: the step's size is at most sqrt(1 + m^2 eps^2) at every momentum, so N = t/eps steps stay below e^(m^2 eps t/2) -> 1 "
      "(dominated convergence then gives the kernel as a limit of distributions)", nrm <= np.sqrt(1 + (m_*e_)**2) + 1e-12, f"{nrm:.8f}")

# ---------- D. slow limits: diffusion and Schrodinger ----------
a = 2.0; e_ = 0.002
pp = np.array([0.05, 0.1, 0.2])
rate = []
for pm in pp:
    Sp = np.diag([np.exp(-1j*pm*e_), np.exp(1j*pm*e_)]) @ np.array([[1 - a*e_, a*e_], [a*e_, 1 - a*e_]])
    lam = np.max(np.abs(np.linalg.eigvals(Sp)))
    rate.append(-np.log(lam)/e_)
check("D1 k = -1, long times: the slow mode decays at rate p^2/(2a), i.e. diffusion with D = 1/(2a) (Kac 1974; c = 1)",
      np.allclose(rate, pp**2/(2*a), rtol=0.02), f"{np.round(rate, 5)} vs {pp**2/(2*a)}")
m_ = 2.0
omg = []
for pm in pp:
    Sp = np.diag([np.exp(-1j*pm*e_), np.exp(1j*pm*e_)]) @ np.array([[np.cos(m_*e_), 1j*np.sin(m_*e_)], [1j*np.sin(m_*e_), np.cos(m_*e_)]])
    ph = np.angle(np.linalg.eigvals(Sp)); omg.append(np.min(np.abs(ph))/e_)
omg = np.array(omg)
check("D2 k = -i: the step's frequencies are sqrt(m^2 + p^2) (relativistic), which is m + p^2/(2m) for small p: Schrodinger's free particle",
      np.allclose(omg, np.sqrt(m_**2 + pp**2), rtol=1e-4) and np.allclose(omg - m_, pp**2/(2*m_), rtol=0.01), f"{np.round(omg - m_, 6)} vs {pp**2/(2*m_)}")
check("D3 the turn rate a continued to i m sends the diffusion constant 1/(2a) to -i/(2m): the heat equation becomes Schrodinger's "
      "(Nelson's continuation in mass, here at the level of the chain)", np.isclose(1/(2*(1j*m_)), -1j/(2*m_)))
# D4 telegraph variance exactly
a = 1.3; e_ = 0.002; Nst = 1000
Pm = 0.5*evolve(Nst, e_, a, -1, normalized=True, start_dir=0) + 0.5*evolve(Nst, e_, a, -1, normalized=True, start_dir=1)
c0 = Pm.shape[1]//2; xg = (np.arange(Pm.shape[1]) - c0)*e_
tt = Nst*e_; var = (Pm[-1].real.sum(axis=1)*xg**2).sum()
check("D4 k = -1: the mean square distance is t/a - (1 - e^(-2at))/(2a^2), which becomes the diffusive t/a at long times",
      np.isclose(var, tt/a - (1 - np.exp(-2*a*tt))/(2*a**2), rtol=3e-3), f"{var:.5f} vs {tt/a - (1 - np.exp(-2*a*tt))/(2*a**2):.5f}")

# ---------- E. why the quantum dial has a real measure here (Ichinose 1982; Ichinose-Tamura) ----------
ok = True; tv = []
for N in [50, 100, 200, 400]:
    e_ = 2.0/N; mq = 1.3
    Aq = evolve(N, e_, mq, -1j)[-1]; Ac = evolve(N, e_, mq, -1)[-1]
    ok &= np.all(np.abs(Aq) <= Ac.real + 1e-12)
    tv.append(Ac.real.sum())
check("E1 path by path, the k = -i weight (i m eps)^R has the same size as the k = -1 weight (m eps)^R; summed, the sizes are "
      "(1 + m eps)^N -> e^(m t), finite as eps -> 0: a complex path measure of bounded total variation (Cameron's no-go does not bite here)",
      ok and np.allclose(tv, [(1 + 1.3*2.0/N)**N for N in [50, 100, 200, 400]]) and abs(tv[-1]/np.exp(1.3*2.0) - 1) < 0.01,
      f"total size {np.round(tv, 3)} vs e^(mt) = {np.exp(2.6):.3f}")
print(f"\n{sum(o for _, o in results)}/{len(results)} pass")
