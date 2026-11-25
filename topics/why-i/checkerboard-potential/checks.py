"""Checks for topics/why-i/checkerboard-potential (Nov 25, 2026). Run: python3 checks.py"""
import numpy as np, itertools, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "checkerboard"))
from cbv import evolve_v, reference, schrodinger, heat_killing
from cb import circ, move_matrix, swap_matrix
rng = np.random.default_rng(1125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
DIALS = [(-1.0, "k = -1"), (-1j, "k = -i")]

# ---------- A. the potential is one more element in the chain ----------
n, eps, mu = 21, 0.1, 1.3
xs = (np.arange(n) - n//2)*eps
Vx = 1.5*np.cos(2*xs) + 0.7*xs
M = move_matrix(n); X = swap_matrix(n); Id = np.eye(2*n)
Vop = np.kron(np.diag(Vx), np.eye(2))           # state index 2*site + direction
ok = True
for k, _ in DIALS:
    t = -eps*mu*X; v = eps*Vop
    ok &= np.allclose(t @ v, v @ t) and np.allclose((Id + k*v) @ (Id + k*t), Id + k*circ(v, t, k))
check("A1 at one site the potential and the turn commute (no twist), so one tick's factor is 1 + k(v o_k t)", ok)
ok = True
for k, _ in DIALS:
    N = 7; t = -eps*mu*X; v = eps*Vop
    S = M @ (Id + k*v) @ (Id + k*t)
    elems = []
    Minv = np.linalg.inv(M)
    for j in range(N):
        P = np.linalg.matrix_power(Minv, j); Q = np.linalg.matrix_power(M, j)
        elems += [P @ t @ Q, P @ v @ Q]                # tick j: turn first, then potential
    chain = elems[-1]
    for e in reversed(elems[:-1]): chain = circ(chain, e, k)
    ok &= np.allclose(np.linalg.matrix_power(S, N), np.linalg.matrix_power(M, N) @ (Id + k*chain))
check("A2 [ours, exact at finite eps] N ticks with a potential = free moves times 1 + k(chain of all 2N turn and potential elements)", ok)

# ---------- B. every path still counted once ----------
N = 12; n = 2*N + 3; xs = (np.arange(n) - n//2)*eps; Vx = 1.5*np.cos(2*xs) + 0.7*xs
psi0 = np.zeros((n, 2)); psi0[n//2, 0] = 1
for k, lab in DIALS:
    amp = np.zeros((n, 2), dtype=complex)
    for S_ in itertools.product([0, 1], repeat=N):
        d = 0; s = n//2; wgt = 1.0 + 0j
        for turn in S_:
            if turn: d = 1 - d; wgt *= -k*eps*mu
            wgt *= 1 + k*eps*Vx[s]
            s += 1 if d == 0 else -1
        amp[s, d] += wgt
    check(f"B1 {lab}: summing over all 2^12 checkerboard paths, weight (-k eps mu)^R times the product of (1 + k eps V) along the path, "
          "gives the chain exactly: each path once", np.allclose(amp, evolve_v(psi0, N, eps, mu, k, Vx)))
# B2 grouping: for one fixed path, the sum over potential-insertion subsets is the product along the path
path_V = Vx[n//2 + np.array([0, 1, 2, 1, 0, -1, 0, 1, 2, 3, 2, 3])]
for k, lab in DIALS:
    tot = sum(np.prod([k*eps*path_V[j] for j in range(N) if S_[j]]) for S_ in itertools.product([0, 1], repeat=N))
    check(f"B2 {lab}: for one path, the 2^12 ways to insert potential elements sum to the product of (1 + k eps V) along it (the Dyson series in V, path by path)",
          np.isclose(tot, np.prod(1 + k*eps*path_V)))

# ---------- C. smooth limits with the potential ----------
L = 8.0
def grid(e_):
    nn = int(round(L/e_)); return (np.arange(nn) - nn//2)*e_, nn
Vfun = lambda x: 3.0*np.exp(-x**2) + 0.5*x
sig, p0, T = 0.4, 2.0, 1.5
for k, lab, mu_, norm_tel in [(-1.0, "k = -1 (telegraph with killing)", 1.3, True), (-1j, "k = -i (Dirac with V)", 1.3, False)]:
    errs = []
    for e_ in [0.04, 0.02, 0.01]:
        x, nn = grid(e_)
        phi = np.exp(-x**2/(2*sig**2))*(np.exp(1j*p0*x) if k == -1j else 1)
        psi0 = np.stack([phi, 0.4*phi], axis=1)
        Nt = int(round(T/e_))
        ch = evolve_v(psi0, Nt, e_, mu_, k, Vfun(x), linear=False, normalized=True)
        ref = reference(psi0, Nt*e_, e_, mu_, k, Vfun(x), normalized_telegraph=norm_tel, steps=4*Nt)
        errs.append(np.abs(ch - ref).max()/np.abs(ref).max())
    check(f"C1 {lab}: the chain converges to an independent continuum solver (Strang/FFT); error halves with eps",
          errs[1] < 0.6*errs[0] and errs[2] < 0.6*errs[1], f"relative errors {[f'{e:.1e}' for e in errs]}")
# C2 plain (linear) chain converges to the same limit as the exponential one: the gap closes like eps
gaps = []
for e_ in [0.01, 0.005, 0.0025]:
    x, nn = grid(e_); phi = np.exp(-x**2/(2*sig**2))*np.exp(1j*p0*x); psi0 = np.stack([phi, 0.4*phi], axis=1)
    Nt = int(round(1.5/e_))
    a_ = evolve_v(psi0, Nt, e_, 1.3, -1j, Vfun(x), linear=True, normalized=True)
    b_ = evolve_v(psi0, Nt, e_, 1.3, -1j, Vfun(x), linear=False, normalized=True)
    gaps.append(np.abs(a_ - b_).max()/np.abs(b_).max())
check("C2 the plain factor 1 - i eps V and the phase e^(-i eps V) differ by O(eps^2) per tick, so by O(eps) over a fixed time: same limit",
      gaps[1] < 0.6*gaps[0] and gaps[2] < 0.6*gaps[1], f"gaps {[f'{g:.1e}' for g in gaps]}")
x, nn = grid(0.005); phi = np.exp(-x**2/(2*sig**2))*np.exp(1j*p0*x); psi0 = np.stack([phi, 0.4*phi], axis=1)
# C3 unitarity / sub-stochasticity
U = evolve_v(psi0, 300, 0.005, 1.3, -1j, Vfun(x), linear=False, normalized=True, keep=True)
nrm = (np.abs(U)**2).sum(axis=(1, 2))
Vpos = 3*np.exp(-x**2) + 0.2
P0 = np.stack([np.exp(-x**2), np.exp(-x**2)], axis=1); P0 /= P0.sum()
P = evolve_v(P0, 300, 0.005, 1.3, -1.0, Vpos, linear=False, normalized=True, keep=True).real
check("C3 k = -i with the phase e^(-i eps V): total |amplitude|^2 stays exactly constant (unitary, an electric quantum walk); "
      "k = -1 with e^(-eps V), V >= 0: probabilities stay >= 0 and their total only falls (killing)",
      np.allclose(nrm, nrm[0]) and np.all(P >= -1e-15) and np.all(np.diff(P.sum(axis=(1, 2))) <= 1e-15))

# ---------- D. long times: Schrodinger with V and diffusion with killing (Feynman-Kac) ----------
def nonrel(m, e_, x0, w=1.0, L2=6.0):
    nn = int(round(L2/e_)); x = (np.arange(nn) - nn//2)*e_
    Vh = 0.5*m*w**2*x**2; s0 = 1/np.sqrt(m*w)
    chi0 = np.exp(-(x - x0)**2/(2*s0**2)); chi0 /= np.sqrt((np.abs(chi0)**2).sum()*e_)
    psi0 = np.stack([chi0/np.sqrt(2), -chi0/np.sqrt(2)], axis=1)     # positive-energy spinor at rest (H = p sz - m sx)
    Nq = int(round(np.pi/w/e_))
    psi = evolve_v(psi0, Nq, e_, m, -1j, Vh, linear=False, normalized=True)
    chi_c = (psi[:, 0] - psi[:, 1])/np.sqrt(2)
    chi_s = schrodinger(chi0, Nq*e_, e_, m, Vh, steps=4000)
    ov = abs(np.vdot(chi_s, chi_c))/np.sqrt((np.abs(chi_s)**2).sum()*(np.abs(chi_c)**2).sum())
    return ov, (x*np.abs(chi_c)**2).sum()/(np.abs(chi_c)**2).sum()
# slow swing: amplitude 0.1, top speed 0.1 (light speed = 1)
ov10, xc10 = nonrel(10.0, 0.005, 0.1); ov40, xc40 = nonrel(40.0, 0.00125, 0.1)
# fast swing: amplitude 0.5, top speed 0.5, i.e. half light speed: the chain follows Dirac, not Schrodinger
ovf, xcf = nonrel(40.0, 0.00125, 0.5)
check("D1 k = -i, slow packet in a harmonic V (omega = 1, top speed 0.1): after half a period the chain's packet has swung from x = 0.1 to -0.1 "
      "and matches Schrodinger's with V, better as m grows; a fast swing (top speed 0.5) does not, because it is relativistic",
      ov40 > ov10 > 0.99 and abs(xc40 + 0.1) < 0.005 and ovf < 0.95,
      f"slow: overlap {ov10:.4f} at m = 10, {ov40:.4f} at m = 40, centre {xc40:.3f}; fast: overlap {ovf:.4f}")
a = 20.0; e_ = 0.002; nn = int(round(6.0/e_)); x = (np.arange(nn) - nn//2)*e_
c = 4.0; Vk = 0.5*c*x**2
u0 = np.exp(-x**2/(2*0.3**2)); u0 /= u0.sum()*e_
P0 = np.stack([u0/2, u0/2], axis=1)*e_
Tk = 1.0; Nk = int(round(Tk/e_))
Pt = evolve_v(P0, Nk, e_, a, -1.0, Vk, linear=False, normalized=True).real
u_chain = Pt.sum(axis=1)/e_; u_ref = heat_killing(u0, Nk*e_, e_, 1/(2*a), Vk, steps=4000)
rel = np.abs(u_chain - u_ref).sum()/np.abs(u_ref).sum()
check("D2 k = -1, fast turning (a = 20) with killing V = 2 x^2: the chain's density matches the Feynman-Kac equation u_t = (1/2a) u'' - V u",
      rel < 0.05, f"relative L1 difference {rel:.3f}")
# D3 slow decay rates = oscillator levels; on the lattice D = (1 - a eps)/(2a) (the persistent walk), which -> 1/(2a)
def slow_rates(a, e_, c, half_width=2.0):
    nn = int(round(2*half_width/e_)); x = (np.arange(nn) - nn//2)*e_; Vk = 0.5*c*x**2
    Mm = np.zeros((2*nn, 2*nn))
    for j in range(2*nn):
        b = np.zeros(2*nn); b[j] = 1
        Mm[:, j] = evolve_v(b.reshape(nn, 2), 1, e_, a, -1.0, Vk, linear=False, normalized=True).real.reshape(-1)
    lam = np.sort(np.abs(np.linalg.eigvals(Mm)))[::-1]
    return np.unique(np.round(-np.log(lam[:12])/e_, 4))[:3]
a = 20.0
r1 = slow_rates(a, 0.01, c); r2 = slow_rates(a, 0.004, c)
Om_lat = lambda e_: np.sqrt(2*c*(1 - a*e_)/(2*a))
ok = np.allclose(r1, Om_lat(0.01)*(np.arange(3) + .5), rtol=0.03) and np.allclose(r2, Om_lat(0.004)*(np.arange(3) + .5), rtol=0.03)
check("D3 k = -1: the slowest decay rates with V = (c/2) x^2 are Omega (n + 1/2), Feynman-Kac's oscillator in imaginary time. On the lattice "
      "Omega = sqrt(2 c D) with D = (1 - a eps)/(2a), which tends to sqrt(c/a)", ok,
      f"eps 0.01: {r1} vs {np.round(Om_lat(0.01)*(np.arange(3)+.5), 4)}; eps 0.004: {r2} vs {np.round(Om_lat(0.004)*(np.arange(3)+.5), 4)}; limit {np.round(np.sqrt(c/a)*(np.arange(3)+.5), 4)}")

# ---------- E. the measure at the quantum dial survives ----------
N = 200; e_ = 0.01; mq = 1.3; n = 2*N + 3; x = (np.arange(n) - n//2)*e_; Vb = 3*np.cos(3*x) + 2*x
psi0 = np.zeros((n, 2)); psi0[n//2, 0] = 1
Aq = evolve_v(psi0, N, e_, mq, -1j, Vb, linear=False)            # turn weight i m eps, phase e^(-i eps V)
Ac = evolve_v(psi0, N, e_, mq, -1.0, np.zeros(n))                 # sizes: telegraph weights with no potential
check("E1 with the phase e^(-i eps V) each path's weight keeps its size (m eps)^R, so the total size is still (1 + m eps)^N -> e^(m t): "
      "the bounded complex path measure survives any real potential", np.all(np.abs(Aq) <= Ac.real + 1e-12) and np.isclose(Ac.real.sum(), (1 + mq*e_)**N))
Al = evolve_v(psi0, N, e_, mq, -1j, Vb, linear=True)
bound = (1 + mq*e_)**N*np.prod(np.sqrt(1 + (e_*np.abs(Vb).max())**2)**N)
check("E2 with the plain factor 1 - i eps V the sizes grow at most by (1 + eps^2 max V^2)^(N/2) -> 1 on the light cone (V is bounded there)",
      np.abs(Al).sum() <= bound, f"{np.abs(Al).sum():.3f} <= {bound:.3f}")
# E3 Feynman-Kac at k = -1 by Monte Carlo over random paths (turn with prob a eps, weight e^(-eps V))
a = 1.3; N = 100; e_ = 0.02; n = 2*N + 3; x = (np.arange(n) - n//2)*e_; Vp = 1 + np.sin(2*x)**2
psi0 = np.zeros((n, 2)); psi0[n//2, 0] = 1
P = evolve_v(psi0, N, e_, a, -1.0, Vp, linear=False, normalized=True).real
f = np.cos(x)
exact = (P.sum(axis=1)*f).sum()
Np = 20000; d = np.zeros(Np, int); s = np.full(Np, n//2); W = np.ones(Np)
for _ in range(N):
    flip = rng.random(Np) < a*e_; d = np.where(flip, 1 - d, d)
    W *= np.exp(-e_*Vp[s]); s += np.where(d == 0, 1, -1)
est = (W*f[s]).mean(); se = (W*f[s]).std()/np.sqrt(Np)
check("E3 k = -1: averaging random Kac paths with weight exp(-sum eps V) reproduces the chain (Feynman-Kac for the telegraph; random evolutions)",
      abs(est - exact) < 4*se, f"paths {est:.4f} +/- {se:.4f}, chain {exact:.4f}")

# ---------- G. a vector potential is a phase on the moves, and it is a gauge choice in 1+1 ----------
N = 40; e_ = 0.05; n = 203; x = (np.arange(n) - n//2)*e_
chi = np.sin(1.7*x) + 0.3*x**2
theta = np.roll(chi, -1) - chi                     # link phase from site s to s+1 = chi(s+1) - chi(s)
theta[-1] = 0; chi_ring = chi                      # the packet never reaches the seam
phi = np.exp(-x**2/0.2); psi0 = np.stack([phi, 0.5*phi], axis=1).astype(complex)
withA = evolve_v(psi0*np.exp(1j*chi)[:, None], N//2, e_, 1.3, -1j, Vfun(x), linear=False, normalized=True, theta=theta)
gauge = evolve_v(psi0, N//2, e_, 1.3, -1j, Vfun(x), linear=False, normalized=True)*np.exp(1j*chi)[:, None]
check("G1 a vector potential enters as link phases on the moves (Feynman-Hibbs's e^(i A dx)); a pure-gauge one only multiplies the "
      "state by e^(i chi(x)), as gauge invariance requires", np.allclose(withA, gauge))
print(f"\n{sum(o for _, o in results)}/{len(results)} pass")
