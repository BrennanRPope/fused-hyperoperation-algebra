"""Checks for topics/why-i/e-i-pi/g-dials (Nov 25, 2026): candidate dials for Newton's G.
Run: python3 checks.py  (numpy only). Writes checks_output.txt next to this file."""
import numpy as np, itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
results, lines = [], []
def check(name, ok, detail=""):
    results.append(bool(ok)); s = ("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""); lines.append(s); print(s)

# CODATA 2018 / Planck 2018
G, hbar, c, kB = 6.67430e-11, 1.054571817e-34, 2.99792458e8, 1.380649e-23
Lam = 1.1056e-52                      # m^-2 (Planck 2018)
mp, Msun = 1.67262192e-27, 1.98847e30
lP = np.sqrt(hbar*G/c**3); EP = np.sqrt(hbar*c**5/G)

# ---------------- D: G needs a partner ----------------
# dimensions (M, L, T) of G, hbar, c, Lambda, mass
dim = {"G": (-1, 3, -2), "hbar": (1, 2, -1), "c": (0, 1, -1), "Lam": (0, -2, 0), "m": (1, 0, 0)}
A = np.array([dim[k] for k in ["G", "hbar", "c", "Lam"]], float).T
_, s, vt = np.linalg.svd(A); null = vt[-1]/vt[-1][0]
check("D1 G alone has units, so no power of G is a pure number; with hbar, c and Lambda there is exactly ONE pure number, G hbar Lambda / c^3",
      np.linalg.matrix_rank(A) == 3 and np.allclose(null, [1, 1, -3, 1]), f"exponents (G, hbar, c, Lambda) = {np.round(null, 6)}")
A2 = np.array([dim[k] for k in ["G", "hbar", "c", "m"]], float).T; null2 = np.linalg.svd(A2)[2][-1]; null2 = null2/null2[0]
check("D2 with a mass m instead of Lambda the one pure number is G m^2/(hbar c) (gravitational fine-structure constant)",
      np.allclose(null2, [1, -1, -1, 2]), f"exponents (G, hbar, c, m) = {np.round(null2, 6)}")
aL = hbar*G*Lam/c**3; aG = G*mp**2/(hbar*c)
check("D3 sizes: G hbar Lambda/c^3 = (l_P)^2 Lambda ~ 2.9e-122; for the proton G m^2/(hbar c) ~ 5.9e-39",
      2.8e-122 < aL < 3.0e-122 and 5.8e-39 < aG < 6.0e-39, f"{aL:.3e}, {aG:.3e}; l_P = {lP:.4e} m, E_P = {EP/1.602176634e-10:.4e} GeV")

# ---------------- S: the stabilized Poincare-Heisenberg algebra (Mendes 1994; Chryssomalakos-Okon 2004) ----------------
def so_G(gdiag):
    """orthogonal algebra of a diagonal (possibly degenerate) metric: [M_AB, M_CD] = G_BC M_AD - G_AC M_BD - G_BD M_AC + G_AD M_BC"""
    n = len(gdiag); idx = list(itertools.combinations(range(n), 2)); pos = {p: a for a, p in enumerate(idx)}; d = len(idx)
    def M(a, b):  # returns (index, sign)
        if a == b: return None
        return (pos[(a, b)], 1.0) if a < b else (pos[(b, a)], -1.0)
    Gm = np.diag(gdiag); C = np.zeros((d, d, d))
    for i, (A_, B_) in enumerate(idx):
        for j, (C_, D_) in enumerate(idx):
            for coef, (a, b) in [(Gm[B_, C_], (A_, D_)), (-Gm[A_, C_], (B_, D_)), (-Gm[B_, D_], (A_, C_)), (Gm[A_, D_], (B_, C_))]:
                if coef and a != b:
                    k, sg = M(a, b); C[i, j, k] += coef*sg
    return C, idx
def jacobi_ok(C):
    n = C.shape[0]; E = np.eye(n); br = lambda x, y: np.einsum('i,j,ijk->k', x, y, C)
    return all(np.allclose(br(E[i], br(E[j], E[k])) + br(E[j], br(E[k], E[i])) + br(E[k], br(E[i], E[j])), 0)
               for i, j, k in itertools.combinations(range(n), 3))
def H2_dim(C):
    """dim H^2(g, g) (Chevalley-Eilenberg, adjoint coefficients), vectorized"""
    n = C.shape[0]; pairs = np.array(list(itertools.combinations(range(n), 2))); trip = np.array(list(itertools.combinations(range(n), 3)))
    P_, Q_ = pairs[:, 0], pairs[:, 1]; d = np.eye(n)
    # d1: phi_ab(e_j) = delta_ja e_b ;  (d phi)(e_i, e_j) = [e_i, phi e_j] - [e_j, phi e_i] - phi([e_i, e_j])
    D1 = (np.einsum('ja,ibk->ijkab', d, C) - np.einsum('ia,jbk->ijkab', d, C) - np.einsum('ija,bk->ijkab', C, d))
    d1 = D1[P_, Q_].reshape(len(pairs)*n, n*n)
    # C2 basis: Psi[c] antisymmetric in its first two slots
    B = np.zeros((len(pairs)*n, n, n, n))
    for r, (p, q) in enumerate(pairs):
        for b in range(n): B[r*n + b, p, q, b] = 1; B[r*n + b, q, p, b] = -1
    I, J, L = trip[:, 0], trip[:, 1], trip[:, 2]
    out = (np.einsum('tmk,ctm->ctk', C[I], B[:, J, L, :]) - np.einsum('tmk,ctm->ctk', C[J], B[:, I, L, :])
           + np.einsum('tmk,ctm->ctk', C[L], B[:, I, J, :])
           - np.einsum('tm,cmtk->ctk', C[I, J], B[:, :, L, :]) + np.einsum('tm,cmtk->ctk', C[I, L], B[:, :, J, :])
           - np.einsum('tm,cmtk->ctk', C[J, L], B[:, :, I, :]))
    d2 = out.reshape(len(pairs)*n, len(trip)*n).T
    return len(pairs)*n - np.linalg.matrix_rank(d2, tol=1e-9) - np.linalg.matrix_rank(d1, tol=1e-9)

# 1+1: axes (0 time, 1 space, 2 'x-axis', 3 'p-axis'); J = M01, X_mu = M_mu2, P_mu = M_mu3, I = M23
def stab11(g2, g3): return so_G([-1, 1, g2, g3])
corners = {(0, 0): None, (1, 0): None, (0, 1): None, (1, 1): None, (1, -1): None}
ok_j = True
for (a, b) in list(corners):
    C, idx = stab11(float(a), float(b)); ok_j &= jacobi_ok(C); corners[(a, b)] = H2_dim(C)
C00, idx = stab11(0., 0.); pos = {p: k for k, p in enumerate(idx)}
X0, X1, P0, P1, I_, J_ = pos[(0, 2)], pos[(1, 2)], pos[(0, 3)], pos[(1, 3)], pos[(2, 3)], pos[(0, 1)]
heis = np.isclose(C00[X0, P0, I_], 1) and np.isclose(C00[X1, P1, I_], -1) and np.allclose(C00[:, I_, :], 0) and np.allclose(C00[X0, X1], 0)
check("S1 the two-scale algebra (Lorentz + positions X + momenta P + I) obeys Jacobi for every (g2, g3); at g2 = g3 = 0 it is Poincare-Heisenberg: "
      "[X_mu, P_nu] ~ eta_mu_nu I, I central, positions commute", ok_j and heis)
check("S2 stability (same test as for c and hbar, H^2(g,g)): Poincare-Heisenberg (g2 = g3 = 0) CAN be deformed; turning on both scales gives a "
      "rigid algebra (H^2 = 0); one scale alone is not enough", corners[(0, 0)] > 0 and corners[(1, 1)] == 0 and corners[(1, -1)] == 0
      and corners[(1, 0)] > 0 and corners[(0, 1)] > 0, f"dim H^2 at (g2, g3) = (0,0) {corners[(0,0)]}, (1,0) {corners[(1,0)]}, (0,1) {corners[(0,1)]}, (1,1) {corners[(1,1)]}, (1,-1) {corners[(1,-1)]}")
C11, _ = stab11(0.3, 0.7)
check("S3 the scales show up where expected: [X0, X1] = -g2 J (positions stop commuting: a minimal length l^2), [P0, P1] = -g3 J (curved momenta / de Sitter: 1/L^2), "
      "[X, I] = g2 P, [P, I] = -g3 X", np.isclose(C11[X0, X1, J_], -0.3) and np.isclose(C11[P0, P1, J_], -0.7)
      and np.isclose(C11[X0, I_, P0], 0.3) and np.isclose(C11[P0, I_, X0], -0.7))
C31, idx31 = so_G([-1, 1, 1, 1, 0, 0]); C31s, _ = so_G([-1, 1, 1, 1, 1, 1])
h31 = (H2_dim(C31), H2_dim(C31s))
check("S4 the same in 3+1 (15 generators): Poincare-Heisenberg can be deformed, the two-scale algebra is rigid (Mendes 1994; Chryssomalakos-Okon 2004)",
      jacobi_ok(C31s) and h31[0] > 0 and h31[1] == 0, f"dim H^2 = {h31}")
# invariant: rescale X -> aX, P -> P/a keeps [X,P] = I: g2 -> a^2 g2, g3 -> g3/a^2, product fixed
a = 2.7; Cr, _ = stab11(0.3*a**2, 0.7/a**2)
T = np.eye(6); T[X0, X0] = T[X1, X1] = a; T[P0, P0] = T[P1, P1] = 1/a     # new basis X' = aX, P' = P/a
Ct = np.einsum('ia,jb,abk,kl->ijl', T, T, C11, np.linalg.inv(T))
check("S5 with hbar fixed by [x, p] = i hbar, the only pure number the rigid algebra carries is g2*g3 = l^2/L^2 (rescaling positions and momenta oppositely leaves it, and the algebra, unchanged)",
      np.allclose(Ct, Cr))
LdS2 = 3/Lam
check("S6 if l is the Planck length and L the de Sitter radius (L^2 = 3/Lambda), l^2/L^2 = G hbar Lambda/(3 c^3): the SAME pure number as D1, "
      "the coupling MacDowell-Mansouri and Freidel-Starodubtsev expand in (link ours)", np.isclose(lP**2/LdS2, aL/3), f"l_P^2/L^2 = {lP**2/LdS2:.3e}")
check("S7 so G = c^3 l^2 / hbar: if stability's length l is the Planck length, G is not a dial of its own but l's square in other units",
      np.isclose(c**3*lP**2/hbar, G))

# ---------------- I: Barbero-Immirzi gamma: a dial on the Hodge-star direction ----------------
def hodge(sig):
    eta = np.diag(sig); inv = np.linalg.inv(eta); pairs = list(itertools.combinations(range(4), 2))
    eps = np.zeros((4,)*4)
    for p in itertools.permutations(range(4)):
        eps[p] = np.linalg.det(np.eye(4)[list(p)])
    S = np.zeros((6, 6))
    for col, (a_, b_) in enumerate(pairs):
        F = np.zeros((4, 4)); F[a_, b_] = 1; F[b_, a_] = -1
        Fup = inv @ F @ inv; star = 0.5*np.einsum('mnrs,rs->mn', eps, Fup)
        S[:, col] = [star[i, j] for (i, j) in pairs]
    return S
SL, SE = hodge([-1, 1, 1, 1]), hodge([1, 1, 1, 1])
check("I1 on curvature 2-forms the Lorentzian Hodge star squares to -1 (eigenvalues +-i, three each): it plays i; the Euclidean one squares to +1 (known; Riemann-Silberstein in other clothes)",
      np.allclose(SL @ SL, -np.eye(6)) and np.allclose(SE @ SE, np.eye(6)) and np.allclose(sorted(np.linalg.eigvals(SL).imag), [-1]*3 + [1]*3))
ks = np.linspace(-5, 5, 2001)
detL = np.array([np.linalg.det(np.eye(6) + k*SL) for k in ks])
check("I2 Holst's operator 1 + (1/gamma)*star: det = (1 + 1/gamma^2)^3 in Lorentzian signature, never 0 for real gamma (Barbero's real dial, min 1 at gamma = infinity, the Palatini end) ...",
      np.allclose(detL, (1 + ks**2)**3) and detL.min() >= 1 - 1e-9)
zL = [np.linalg.det(np.eye(6) + (1/g)*SL.astype(complex)) for g in [1j, -1j]]
P = 0.5*(np.eye(6) - 1j*SL)
check("I3 ... and 0 exactly at gamma = +-i: the self-dual (Ashtekar) value is the WALL of this lift, where it becomes a rank-3 projector (reading ours)",
      np.allclose(zL, 0) and np.allclose(P @ P, P) and np.linalg.matrix_rank(P) == 3)
detE = np.array([np.linalg.det(np.eye(6) + k*SE) for k in ks])
check("I4 Euclidean signature: det = (1 - 1/gamma^2)^3, wall at gamma = +-1, which is where Benedetti-Speziale (2011) find a UV fixed point gamma^2 = 1",
      np.allclose(detE, (1 - ks**2)**3) and np.isclose(np.linalg.det(np.eye(6) + SE), 0))
js = np.arange(1, 4001)/2; root_j = np.sqrt(js*(js + 1))
def solve(weights):
    f = lambda g: np.sum(weights*np.exp(-2*np.pi*g*root_j)) - 1
    lo, hi = 0.01, 2.0
    for _ in range(200):
        mid = 0.5*(lo + hi); lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5*(lo + hi)
gDL = solve(2*np.ones_like(js)); gGM = solve(2*js + 1); gABCK = np.log(2)/(np.pi*np.sqrt(3))
check("I5 black-hole state counting fixes gamma: sum_j 2 exp(-2 pi gamma sqrt(j(j+1))) = 1 gives 0.2375 (Domagala-Lewandowski, Meissner 2004)",
      abs(gDL - 0.2375) < 5e-4, f"gamma = {gDL:.5f}")
check("I6 counting each puncture with its 2j+1 states instead gives 0.274 (Ghosh-Mitra); spin 1/2 only gives ln2/(pi sqrt3) = 0.127 (ABCK 1998)",
      abs(gGM - 0.274) < 1e-3 and abs(gABCK - 0.1274) < 1e-4, f"{gGM:.5f}, {gABCK:.5f}")
Amin = 4*np.sqrt(3)*np.pi*gDL*lP**2
check("I7 gamma multiplies G in the area gap, A_min = 4 sqrt3 pi gamma hbar G/c^3: gamma rescales G's quantum of area but is not a dial on G itself",
      np.isclose(Amin, 8*np.pi*gDL*hbar*G/c**3*np.sqrt(0.75)), f"A_min = {Amin:.3e} m^2")

# ---------------- J: Jacobson / Sakharov: G as inverse entropy per area ----------------
eta = c**3/(4*G*hbar)
rs = 2*G*Msun/c**2; S_sun = eta*4*np.pi*rs**2
check("J1 Jacobson: entropy per area eta = c^3/(4 G hbar), so G = c^3/(4 hbar eta): G -> 0 means infinitely many states per area (stiff spacetime). Solar-mass hole: S ~ 1.05e77 k_B",
      np.isclose(c**3/(4*hbar*eta), G) and 1.0e77 < S_sun < 1.1e77, f"S = {S_sun:.3e} k_B")
NA, NB = 37.0, 112.0
check("J2 state counts of two horizon patches multiply, so (N - 1) composes by mult-add at dial 1 and entropy log N = log(1 + (N-1)) is its flat coordinate (reading ours; Tsallis's rule is the general dial, gr-bridges J5)",
      np.isclose(NA*NB - 1, (NA - 1) + (NB - 1) + (NA - 1)*(NB - 1)) and np.isclose(np.log(NA*NB), np.log(NA) + np.log(NB)))

# ---------------- W: the horizon as a wall ----------------
k_w = -2*G/c**2; x = np.linspace(0, 2*c**2/(2*G), 5)
gtt = 1 + k_w*(Msun/np.array([1e9, rs, 0.5*rs]))
check("W1 Schwarzschild g_tt = 1 + k x with x = M/r and k = -2G/c^2: the lift hits 0 at the horizon r = 2GM/c^2 (the wall -1/k) and is negative inside, where time and radius trade roles (reading ours)",
      gtt[0] > 0 and np.isclose(gtt[1], 0) and gtt[2] < 0)
kap = 0.5*c**2*(2*G*Msun/(c**2*rs**2)); TH = hbar*kap/(2*np.pi*c*kB)
check("W2 the slope of g_tt at the wall is the surface gravity c^4/(4GM); with the 2 pi period it gives Hawking's 6.17e-8 K (same as e-i-pi C23)",
      np.isclose(kap, c**4/(4*G*Msun)) and abs(TH - 6.170e-8) < 2e-10, f"T = {TH:.4e} K")
N = lambda r: np.sqrt(1 - rs/r); r1, r2, r3 = 3*rs, 7*rs, 40*rs
n12, n23 = N(r1)/N(r2), N(r2)/N(r3)
check("W3 redshift factors between shells multiply, so (factor - 1) composes by mult-add at dial 1, but that dial is set by c (1/c^2 in the potential), not by G: masses themselves do not compose this way (dial-map G1-G2)",
      np.isclose(N(r1)/N(r3) - 1, (n12 - 1) + (n23 - 1) + (n12 - 1)*(n23 - 1)))

# ---------------- other candidates ----------------
r, M = 1e7, 6e24
pot = [cc**2*(np.sqrt(1 - 2*G*M/(r*cc**2)) - 1) for cc in [1e5, 1e6, 1e8]]
check("E1 Ehlers' frame theory dials lambda = 1/c^2 with G held fixed: c^2 (sqrt(g_tt) - 1) -> -GM/r as c -> infinity, so G survives the Newton-Cartan end; it is a c dial, not a G dial",
      abs(pot[-1] + G*M/r) < 1e-6*G*M/r and abs(pot[0] + G*M/r) > abs(pot[-1] + G*M/r))
scal = np.array([1.0, 1e-3, 1e-6]); hb_, G_ = hbar*scal, G*scal
EPs, lPs = np.sqrt(hb_*c**5/G_), np.sqrt(hb_*G_/c**3)
check("R1 relative locality (Amelino-Camelia, Freidel, Kowalski-Glikman, Smolin 2011): hbar, G -> 0 with hbar/G fixed keeps the Planck energy and sends the Planck length to 0, isolating E_P as a dial",
      np.allclose(EPs, EPs[0]) and np.allclose(lPs/lPs[0], scal))
kap_ = 1.0; pa, pb = 0.3, 0.5
ma = kap_*(1 - np.exp(-pa/kap_)); mb = kap_*(1 - np.exp(-pb/kap_)); mab = kap_*(1 - np.exp(-(pa + pb)/kap_))
check("K1 kappa-Poincare: massless momenta moving together compose as a + b + k ab with k = -1/E_P (why-i/kappa-poincare K1, re-run)",
      np.isclose(mab, ma + mb - ma*mb/kap_))
beta = 0.04; dps = np.linspace(0.1, 30, 200001); dx = 0.5*hbar*(1/dps + beta*dps)
check("U1 GUP [x, p] = i hbar (1 + beta p^2) gives a minimal length dx_min = hbar sqrt(beta) (Kempf-Mangano-Mann 1995); beta ~ l_P^2/hbar^2 inherits G through l_P",
      np.isclose(dx.min(), hbar*np.sqrt(beta), rtol=1e-6))
gs = 0.7
def flow_inv(w, t): return np.exp(-2*t)*(w - 1/gs) + 1/gs       # 1/g under beta(g) = 2g - 2g^2/g*
g0 = 1e-4; ts = np.linspace(0, 12, 4001); gt = 1/flow_inv(1/g0, ts)
dg = np.gradient(gt, ts); beta_g = 2*gt - 2*gt**2/gs
check("A1 asymptotic safety, toy flow beta(g) = 2g - 2g^2/g* for g = G k^2: g -> 0 in the infrared (Gaussian point), g -> g* in the ultraviolet (Weinberg 1979, Reuter 1998)",
      np.allclose(dg[5:-5], beta_g[5:-5], rtol=1e-3, atol=1e-8) and gt[0] < 1e-3 and abs(gt[-1] - gs) < 1e-6)
w = 3.3; t1, t2 = 0.4, 1.1
lam1 = np.exp(-2*t1); s1 = (1/gs)*(1 - lam1)
check("A2 in the variable 1/g those flow maps are stretch-and-shift maps x -> lam x + s with lam = e^(-2t), composing like the order dial's affine group with fixed point 1/g* (weak fit: any constant-coefficient Riccati flow does this)",
      np.isclose(flow_inv(flow_inv(w, t1), t2), flow_inv(w, t1 + t2)) and np.isclose(flow_inv(w, t1), lam1*w + s1))

n_ok = sum(results); lines.append(""); lines.append(f"{n_ok}/{len(results)} pass"); print(lines[-1])
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(lines) + "\n")
