"""Quantum noise through the mult-add cube: numerical checks.
Density-matrix simulations in plain numpy (no quantum library)."""
import numpy as np, itertools
rng = np.random.default_rng(7)
I = np.eye(2); X = np.array([[0,1],[1,0]],complex); Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1,-1]).astype(complex)
P = [I, X, Y, Z]

def apply(kraus, rho): return sum(K @ rho @ K.conj().T for K in kraus)
def bitflip(p): return [np.sqrt(1-p)*I, np.sqrt(p)*X]
def phaseflip(p): return [np.sqrt(1-p)*I, np.sqrt(p)*Z]
def depol(p): return [np.sqrt(1-p)*I] + [np.sqrt(p/3)*M for M in (X,Y,Z)]
def ampdamp(g): return [np.array([[1,0],[0,np.sqrt(1-g)]]), np.array([[0,np.sqrt(g)],[0,0]])]
def compose(k1, k2): return [B @ A for A in k1 for B in k2]   # k1 first, then k2

def ptm(kraus):
    """Pauli transfer matrix R_ij = tr(P_i E(P_j))/2"""
    return np.real(np.array([[np.trace(Pi @ apply(kraus, Pj))/2 for Pj in P] for Pi in P]))

def fit(name, chan, comb, pars):
    p, q = pars
    a = ptm(compose(chan(p), chan(q))); b = ptm(chan(comb(p, q)))
    return np.abs(a-b).max()

out = []
def say(s): print(s); out.append(s)

say("== 1. How two noise channels combine ==")
p, q = 0.13, 0.27
say(f"bit-flip  p+q-2pq : max PTM error {fit('bf', bitflip, lambda p,q: p+q-2*p*q, (p,q)):.2e}")
say(f"bit-flip  p+q-pq  : max PTM error {fit('bf', bitflip, lambda p,q: p+q-p*q, (p,q)):.2e}   <- the call's rule, wrong for bit-flips")
say(f"phase-flip p+q-2pq: max PTM error {fit('pf', phaseflip, lambda p,q: p+q-2*p*q, (p,q)):.2e}")
say(f"depolarizing p+q-(4/3)pq: max PTM error {fit('dp', depol, lambda p,q: p+q-4/3*p*q, (p,q)):.2e}")
say(f"amplitude damping g+h-gh: max PTM error {fit('ad', ampdamp, lambda p,q: p+q-p*q, (p,q)):.2e}   <- the call's rule fits here (T1 decay)")

say("\n== 2. Every rule is mult-add with a weight k: a + b + k*a*b, and (1 + k*a) multiplies ==")
for name, k, chan in [("amplitude damping", -1, ampdamp), ("bit-flip", -2, bitflip), ("depolarizing", -4/3, depol)]:
    n = 6; pp = 0.05
    K = chan(pp)
    for _ in range(n-1): K = compose(K, chan(pp))
    R = ptm(K)
    pred = (1 - (1 + k*pp)**n) / (-k)   # n-fold weighted mult-add
    if name == "amplitude damping": got = 1 - R[1,1]**2   # R_xx = sqrt(1-g)
    elif name == "bit-flip": got = (1 - R[3,3]) / 2
    else: got = (1 - R[1,1]) * 3/4
    say(f"{name:18s} k={k:+.3f}: 6 stacked at 0.05 -> predicted {pred:.6f}, density matrix {got:.6f}, wall at p*=-1/k={-1/k:.3f}")

say("\n== 3. Pauli channels live in a tetrahedron inside the cube [-1,1]^3; composing multiplies coordinates ==")
def pauli_chan(pr): return [np.sqrt(pr[i])*P[i] for i in range(4)]
ok = True
for _ in range(2000):
    a = rng.dirichlet([1,1,1,1]); b = rng.dirichlet([1,1,1,1])
    Ra, Rb = ptm(pauli_chan(a)), ptm(pauli_chan(b)); Rc = ptm(compose(pauli_chan(a), pauli_chan(b)))
    la, lb, lc = np.diag(Ra)[1:], np.diag(Rb)[1:], np.diag(Rc)[1:]
    ok &= np.allclose(la*lb, lc)
    # tetrahedron test: 1 +- ... >= 0
    s = lc
    ok &= min(1+s[0]+s[1]+s[2], 1+s[0]-s[1]-s[2], 1-s[0]+s[1]-s[2], 1-s[0]-s[1]+s[2]) >= -1e-12
say(f"2000 random pairs: coordinates multiply and stay in the tetrahedron: {ok}")

say("\n== 4. 3-qubit bit-flip code, full density-matrix simulation ==")
def kron(*ms):
    r = np.array([[1]], complex)
    for m in ms: r = np.kron(r, m)
    return r
def op_on(M, i, n=3): return kron(*[M if j == i else I for j in range(n)])
# encoder isometry V: |0> -> |000>, |1> -> |111>
V = np.zeros((8,2), complex); V[0,0] = 1; V[7,1] = 1
proj = {}
for s in itertools.product([0,1], repeat=2):  # syndrome of Z1Z2, Z2Z3
    Pm = np.zeros((8,8))
    for idx in range(8):
        b = [(idx >> (2-j)) & 1 for j in range(3)]
        if (b[0]^b[1], b[1]^b[2]) == s: Pm[idx,idx] = 1
    proj[s] = Pm
fix = {(0,0): None, (1,0): 0, (1,1): 1, (0,1): 2}
def logical_channel(px, pz):
    single = compose(bitflip(px), phaseflip(pz))
    R = np.zeros((4,4))
    for j, Pj in enumerate(P):
        rho = V @ Pj @ V.conj().T
        for i in range(3):
            K = [kron(*[k if m == i else I for m in range(3)]) for k in single]
            rho = apply(K, rho)
        out_rho = np.zeros((8,8), complex)
        for s, Pm in proj.items():
            r = Pm @ rho @ Pm
            if fix[s] is not None: Xi = op_on(X, fix[s]); r = Xi @ r @ Xi
            out_rho += r
        dec = V.conj().T @ out_rho @ V
        for i2, Pi in enumerate(P): R[i2, j] = np.real(np.trace(Pi @ dec))/2
    return R
f = lambda e: (3*e - e**3)/2
for px, pz in [(0.1, 0.02), (0.2, 0.1), (0.3, 0.05)]:
    R = logical_channel(px, pz); ex, ez = 1-2*px, 1-2*pz
    # lambda_Z = polarization against X errors, lambda_X = against Z errors
    say(f"px={px}, pz={pz}: logical lambda_Z {R[3,3]:.6f} vs f(ex)={f(ex):.6f};  logical lambda_X {R[1,1]:.6f} vs ez^3={ez**3:.6f}")

say("\n== 5. Shor code and concatenation (code capacity, independent X and Z flips, perfect decoding) ==")
# Track polarizations (x, z) = (1-2*px, 1-2*pz) against X and Z errors.
# Physical X errors: inner majority, then outer parity -> they come out as a logical Z (phase) error.
# Physical Z errors: inner parity, then outer majority -> they come out as a logical X (bit) error.
def shor(x, z):
    return f(z**3), f(x)**3          # (new x, new z): note the swap of roles
def shor_mc(px, pz, N=400000):
    ex = rng.random((N, 9)) < px; ez = rng.random((N, 9)) < pz
    zbar = (ex.reshape(N,3,3).sum(2) >= 2).sum(1) % 2 == 1      # X errors -> logical Z
    xbar = ((ez.reshape(N,3,3).sum(2) % 2) == 1).sum(1) >= 2    # Z errors -> logical X
    return xbar, zbar
for pp in [0.05, 0.1, 0.2]:
    x = 1-2*pp; X1, Z1 = shor(x, x)
    xb, zb = shor_mc(pp, pp)
    say(f"p={pp}: logical X-error {xb.mean():.5f} (MC) vs {(1-X1)/2:.5f} (map);  logical Z-error {zb.mean():.5f} (MC) vs {(1-Z1)/2:.5f} (map)")
# two levels by Monte Carlo (81 qubits): level-1 logical errors feed level 2
def shor_level2(pp, N=60000):
    ex = rng.random((N, 9, 9)) < pp; ez = rng.random((N, 9, 9)) < pp
    zb1 = (ex.reshape(N,9,3,3).sum(3) >= 2).sum(2) % 2 == 1
    xb1 = ((ez.reshape(N,9,3,3).sum(3) % 2) == 1).sum(2) >= 2
    zb2 = (xb1.reshape(N,3,3).sum(2) >= 2).sum(1) % 2 == 1     # level-1 X errors act as X at level 2
    xb2 = ((zb1.reshape(N,3,3).sum(2) % 2) == 1).sum(1) >= 2
    return xb2.mean(), zb2.mean()
for pp in [0.03, 0.05, 0.08]:
    x = z = 1-2*pp
    x, z = shor(x, z); x, z = shor(x, z)
    a, b = shor_level2(pp)
    say(f"two Shor levels, p={pp}: X {a:.5f} (MC) vs {(1-x)/2:.5f} (map); Z {b:.5f} (MC) vs {(1-z)/2:.5f} (map)")
def survives(pp, levels=80):
    x = z = 1-2*pp
    for _ in range(levels): x, z = shor(x, z)
    return x > 0.99 and z > 0.99
lo, hi = 1e-6, 0.5
for _ in range(60):
    mid = (lo+hi)/2
    lo, hi = (mid, hi) if survives(mid) else (lo, mid)
say(f"concatenated Shor threshold (equal independent X and Z flips, hierarchical decoding): p_th = {lo:.6f}")
def th(step):
    lo, hi = 0, 1
    for _ in range(60):
        m = (lo+hi)/2; x = m
        for _ in range(60): x = step(x)
        lo, hi = (lo, m) if x > 0.99 else (m, hi)
    return hi
tx = th(lambda x: f(f(x)**9)); tz = th(lambda z: f(f(z**3))**3)
say(f"separately: bit-flips survive below p_x = {(1-tx)/2:.6f}, phase-flips below p_z = {(1-tz)/2:.6f} (the safe region is a box)")
say("concatenated 3-bit repetition, bit flips only: x -> f(x) pushes every x>0 to 1, threshold p=1/2")

say("\n== 6. Repetition code: failing corners of the n-cube by layer ==")
from math import comb
for n in [3,5,7,9]:
    t = (n-1)//2
    layers = [comb(n,j) if j > t else 0 for j in range(n+1)]
    say(f"n={n}: failing corners per layer {layers}; lowest failing layer = {t+1} = distance-based")

say("\n== 7. MacWilliams = flipping the cube's signs (repetition code n=3 and its dual) ==")
import sympy as sp
x_, y_ = sp.symbols('x y')
WC = x_**3 + y_**3; Wd = x_**3 + 3*x_*y_**2
say(f"W_dual(x,y) == W_C(x+y, x-y)/|C| : {sp.simplify(WC.subs({x_: x_+y_, y_: x_-y_}, simultaneous=True)/2 - Wd) == 0}")
pp = sp.Rational(1,10)
und = WC.subs({x_:1-pp, y_:pp}) - (1-pp)**3
via_dual = sp.Rational(2,8)*Wd.subs({x_:1, y_:1-2*pp}) - (1-pp)**3
say(f"undetected error at p=1/10: direct {und} , via dual code in polarization {via_dual}")


say("\n== 8. Surface (toric) code, independent X flips, minimum-weight matching (pymatching) ==")
import scipy.sparse as sps, pymatching
def rep(n):
    return sps.csc_matrix((np.ones(2*n), (np.repeat(np.arange(n),2), np.array([[i,(i+1)%n] for i in range(n)]).ravel())), shape=(n,n))
def toric(L):
    Hr = rep(L)
    H = sps.hstack([sps.kron(Hr, sps.eye(L)), sps.kron(sps.eye(L), Hr.T)], format="csc")
    H1 = sps.csc_matrix(([1], ([0],[0])), shape=(1,L)); H0 = sps.csc_matrix(np.ones((1,L)))
    Lg = sps.block_diag([sps.kron(H1, H0), sps.kron(H0, H1)], format="csc")
    return H, Lg
surf = {}
for L in [6, 10, 14, 18]:
    H, Lg = toric(L); m = pymatching.Matching(H); n = H.shape[1]
    row = []
    for pp in np.round(np.arange(0.07, 0.131, 0.01), 3):
        N = 6000 if L < 14 else 3000
        err = (rng.random((N, n)) < pp).astype(np.uint8)
        syn = (err @ H.T.toarray()) % 2
        corr = m.decode_batch(syn)
        fail = (((err + corr) % 2) @ Lg.T.toarray() % 2).any(1).mean()
        row.append((float(pp), float(fail)))
    surf[L] = row
    say(f"L={L}: " + ", ".join(f"p={a:.2f}:{b:.3f}" for a, b in row))
np.save("surface_results.npy", surf, allow_pickle=True)
say("Curves cross near p = 0.10 (the known matching threshold is about 10.3%; the best possible decoder reaches about 10.9%).")

open("checks_output.txt", "w").write("\n".join(out) + "\n")
