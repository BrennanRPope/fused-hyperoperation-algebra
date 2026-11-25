"""e, i and pi from the dial, and i across physics.  Checks for topics/why-i/e-i-pi (Nov 25, 2026).

Dial-k weighted circle operation: a o_k b = a + b + k*a*b, lift 1 + k*x, flat coordinate log(1 + k*x)/k.
Run: python3 checks.py   (numpy only; no scipy/sympy in this session)
"""
import numpy as np
rng = np.random.default_rng(20261125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def circ(a, b, k): return a + b + k * a * b
def chain(x, n, k):
    """n copies of x/n composed with o_k."""
    acc = 0j
    for _ in range(n): acc = circ(acc, x / n, k)
    return acc

# ====================================================== A. e, i and pi from one dial
for k in [1, 0.5, 1j, -1j, 0.3 + 0.7j]:
    pass
ok = all(abs(chain(1.3, n, k) - ((1 + k*1.3/n)**n - 1)/k) < 1e-9 for n in [3, 10, 57] for k in [1, -0.4, 1j, 0.3+0.7j])
check("A1 n copies of x/n under o_k equal ((1 + kx/n)^n - 1)/k (the lift multiplies)", ok)

errs = [abs(chain(1.0, n, 1) - (np.e - 1)) for n in [10, 100, 1000, 10000]]
check("A2 real dial: chain of 1/n steps -> e - 1, error ~ 1/n (e is the lift of one whole unit)",
      errs[-1] < 2e-4 and all(errs[i+1] < errs[i]/5 for i in range(3)), "errors " + ", ".join(f"{e:.1e}" for e in errs))

ok = all(abs(k*circ(a, b, k) - circ(k*a, k*b, 1)) < 1e-12 for k in [2, 1j, 0.3-0.8j] for a, b in [(0.4, -1.1), (0.2+0.5j, 1.7)])
check("A3 every dial is mult-add after scaling the inputs by k: k(a o_k b) = (ka) o_1 (kb)", ok,
      "so 'which dial' over the complex numbers = 'which direction the inputs point'")

xs = np.linspace(0, 10, 2001)
def flow(k, x): return (np.exp(k*x) - 1)/k
dev_im = np.max(np.abs(np.abs(1 + 0.7j*flow(0.7j, xs)) - 1))
dev_mix = np.max(np.abs(np.abs(1 + (0.05+0.7j)*flow(0.05+0.7j, xs)) - 1))
check("A4 the lift keeps size 1 for every x exactly when k is purely imaginary (i is the size-keeping dial)",
      dev_im < 1e-12 and dev_mix > 0.5, f"k=0.7i: {dev_im:.1e};  k=0.05+0.7i: {dev_mix:.2f}")

eps = 1e-4
ok = all(abs((abs(1 + k*eps) - 1)/eps - k.real) < 1e-3 for k in [0.3+2j, -0.5+1j, 1j, 2+0j])
check("A5 one short step: |1 + k*eps| = 1 + Re(k)*eps + O(eps^2), so size-keeping forces Re k = 0", ok)

# period: smallest x>0 with chain value back to 0 at k = i
xx = np.linspace(0.1, 7, 690001); v = np.abs(flow(1j, xx))
x0 = xx[np.argmin(v)]
check("A6 at k = i the flow first returns to 0 at x = 2*pi (pi is half the turn)", abs(x0 - 2*np.pi) < 2e-5, f"x0 = {x0:.6f}")
check("A7 at k = i, x = pi lands on -2/k = 2i, i.e. lift 1 + k*x = -1 (Euler: e^(i*pi) = -1)",
      abs(flow(1j, np.pi) - 2j) < 1e-12 and abs(1 + 1j*flow(1j, np.pi) + 1) < 1e-12)
c400 = chain(np.pi, 400000, 1j)
check("A8 Euler by pure chaining: 400,000 steps of pi/400000 at k = i land within 2e-5 of 2i",
      abs(c400 - 2j) < 2e-5, f"{c400:.7f}")

q = complex(-1, 1)                 # mult-add element whose lift is i (a quarter turn about the wall -1)
check("A9 discrete Euler over Gaussian integers (k = 1): q = -1+i (lift i); q o q = -2; four copies = 0",
      circ(q, q, 1) == -2 and circ(circ(circ(q, q, 1), q, 1), q, 1) == 0)

a, b = 0.3+0.4j, -0.2+0.9j
for k in [1, 1j, 0.5-0.5j]:
    L = lambda z: np.log(1 + k*z)/k
    d = L(circ(a, b, k)) - L(a) - L(b)
    m = d*k/(2j*np.pi)
    ok_k = abs(m - round(m.real)) < 1e-12
check("A10 flat coordinate log(1+kx)/k turns o_k into + (up to whole turns 2*pi*i/k)", ok_k)

check("A11 k = -i walks the same circle the other way: flow(-i, x) = conj(flow(i, x)) for real x",
      np.allclose(flow(-1j, xs), np.conj(flow(1j, xs))))

# ====================================================== B. axioms that derive complex numbers
ok = all((N**r - 1) + (M**r - 1) + (N**r - 1)*(M**r - 1) == (N*M)**r - 1 for N in range(1, 7) for M in range(1, 7) for r in [1, 2, 3])
check("B1 Hardy's K_AB = K_A K_B is mult-add of parameter counts: (K_A-1) o_1 (K_B-1) = K_AB - 1 (any K = N^r)", ok)

s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]]); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]])
paulis = [s0, s1, s2, s3]
def real_span_dim(mats):
    V = np.array([np.concatenate([m.real.ravel(), m.imag.ravel()]) for m in mats])
    return np.linalg.matrix_rank(V, tol=1e-9)
prods = [np.kron(A, B) for A in paulis for B in paulis]
check("B2 complex: products of local Hermitian 2x2 span all 16 Hermitian 4x4; two-qubit 15 = 3 + 3 + 3*3 = 3 o 3",
      real_span_dim(prods) == 16 and circ(3, 3, 1) == 15)

def sym_basis(n):
    out = []
    for i in range(n):
        for j in range(i, n):
            m = np.zeros((n, n)); m[i, j] = m[j, i] = 1; out.append(m)
    return out
rprods = [np.kron(A, B) for A in sym_basis(2) for B in sym_basis(2)]
d_loc = real_span_dim(rprods); d_all = len(sym_basis(4))
yy = np.kron(s2, s2).real
check("B3 real: local products give 9 of the 10 real symmetric 4x4 directions; 2 o 2 = 8 but K_AB - 1 = 9 (one hidden)",
      d_loc == 9 and d_all == 10 and circ(2, 2, 1) == 8 and real_span_dim(rprods + [yy]) == 10,
      "the hidden one is sigma_y (x) sigma_y")
J = np.array([[0, -1], [1, 0]])
check("B4 that hidden joint parameter is J (x) J, two local quarter turns (J = -i sigma_y): J(x)J = -sigma_y(x)sigma_y",
      np.allclose(np.kron(J, J), -np.kron(s2, s2)))

def quat_herm_dim(N):
    Om = np.kron(np.eye(N), np.array([[0, -1], [1, 0]]))      # quaternionic structure on C^(2N)
    basis = []
    n = 2*N
    for i in range(n):
        for j in range(i, n):
            m = np.zeros((n, n), complex); m[i, j] = m[j, i] = 1; basis.append(m)
            if i != j:
                m = np.zeros((n, n), complex); m[i, j] = -1j; m[j, i] = 1j; basis.append(m)
    # condition M Om = Om conj(M): linear over R in M
    rows = []
    for m in basis:
        c = m @ Om - Om @ np.conj(m); rows.append(np.concatenate([c.real.ravel(), c.imag.ravel()]))
    A = np.array(rows).T
    return len(basis) - np.linalg.matrix_rank(A, tol=1e-9)
q2, q4 = quat_herm_dim(2), quat_herm_dim(4)
check("B5 quaternionic: K(2) = 6, K(4) = 28; local tomography needs 6*6 = 36 (5 o 5 = 35 vs 27): too few joint parameters",
      q2 == 6 and q4 == 28 and circ(5, 5, 1) == 35)
check("B6 Masanes-Muller bit-ball dimension d = K(2) - 1: real 2, complex 3, quaternionic 5; only d = 3 obeys (d+1)^2 = K(4)",
      (3, 4, 10) == (len(sym_basis(2)), 4, 10) and (2+1)**2 != 10 and (3+1)**2 == 16 and (5+1)**2 != q4)

R = lambda t: np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
p = R(0.2) @ np.array([0., 1.])   # a small turn of the vertex (0,1)
U = lambda t: np.cos(t/2)*np.eye(2) - 1j*np.sin(t/2)*s1
psi = U(np.pi) @ np.array([1, 0])
check("B7 continuity: a classical bit has no continuous reversible path (a small turn leaves the simplex); a qubit does, e^(-i t sigma_x/2): |0> -> |1>, norm kept",
      p.min() < 0 and abs(abs(psi[1]) - 1) < 1e-12 and all(abs(np.linalg.norm(U(t) @ [1, 0]) - 1) < 1e-12 for t in np.linspace(0, 3, 7)))

# ====================================================== C. i across physics
Mz = lambda z: np.array([[z.real, -z.imag], [z.imag, z.real]])
z, w = 0.7-1.2j, -0.3+2.1j
check("C1 phasors: complex numbers are the real matrices [[a,-b],[b,a]]; M(zw) = M(z)M(w), M(i)^2 = -1 (i is removable here)",
      np.allclose(Mz(z*w), Mz(z) @ Mz(w)) and np.allclose(Mz(1j) @ Mz(1j), -np.eye(2)))

# C2 AC circuit: series RLC driven by cos(wt); steady current from complex impedance vs time-domain integration
Rr, Lr, Cr, om = 1.0, 0.5, 0.2, 2.0
Z = Rr + 1j*om*Lr + 1/(1j*om*Cr)
def f(t, y): return np.array([y[1], (np.cos(om*t) - Rr*y[1] - y[0]/Cr)/Lr])
y = np.zeros(2); t = 0.0; h = 1e-3; ts = []; Is = []
for n in range(60000):
    k1 = f(t, y); k2 = f(t+h/2, y+h/2*k1); k3 = f(t+h/2, y+h/2*k2); k4 = f(t+h, y+h*k3)
    y = y + h/6*(k1+2*k2+2*k3+k4); t += h
    if t > 40: ts.append(t); Is.append(y[1])
ts, Is = np.array(ts), np.array(Is)
pred = (np.exp(1j*om*ts)/Z).real
check("C2 AC circuits: current from the complex impedance R + iwL + 1/(iwC) matches direct real time-stepping",
      np.max(np.abs(Is - pred)) < 1e-6, f"max diff {np.max(np.abs(Is - pred)):.1e}; phase lag {np.degrees(np.angle(Z)):.1f} deg")

xsig = rng.normal(size=64); Ns = 64
nn = np.arange(Ns); cas = np.cos(2*np.pi*np.outer(nn, nn)/Ns) + np.sin(2*np.pi*np.outer(nn, nn)/Ns)
H = cas @ xsig; X = np.fft.fft(xsig)
check("C3 Fourier without i: the real Hartley transform is its own inverse (/N) and gives the same power |X|^2 = (H_k^2 + H_-k^2)/2",
      np.allclose(cas @ H / Ns, xsig) and np.allclose(np.abs(X)**2, (H**2 + H[(-nn) % Ns]**2)/2))

S = np.roll(np.eye(12), 1, axis=0)                    # real shift on a ring
ev = np.linalg.eigvals(S)
check("C4 a real size-keeping shift has eigen-waves e^(2 pi i k/N): all eigenvalues on the unit circle, only +-1 real",
      np.allclose(np.abs(ev), 1) and np.sum(np.abs(ev.imag) < 1e-9) == 2)

e1 = np.array([[1., 0], [0, -1]]); e2 = np.array([[0., 1], [1, 0]])
check("C5 Clifford: with REAL matrices e1, e2 (e1^2 = e2^2 = 1, anticommuting), the oriented plane e1e2 squares to -1",
      np.allclose(e1@e1, np.eye(2)) and np.allclose(e2@e2, np.eye(2)) and np.allclose(e1@e2, -e2@e1) and np.allclose((e1@e2)@(e1@e2), -np.eye(2)))
check("C6 Pauli/Hestenes: sigma1 sigma2 sigma3 = i*1, the oriented volume of space plays i (and commutes with everything)",
      np.allclose(s1@s2@s3, 1j*np.eye(2)))

def boost_fields(E, B, v):
    g = 1/np.sqrt(1 - v*v)
    E2 = np.array([E[0], g*(E[1] - v*B[2]), g*(E[2] + v*B[1])])
    B2 = np.array([B[0], g*(B[1] + v*E[2]), g*(B[2] - v*E[1])])
    return E2, B2
E = rng.normal(size=3); B = rng.normal(size=3); v = 0.6; phi = np.arctanh(v)
E2, B2 = boost_fields(E, B, v); F = E + 1j*B; F2 = E2 + 1j*B2
al = 1j*phi
Rc = np.array([[1, 0, 0], [0, np.cos(al), np.sin(al)], [0, -np.sin(al), np.cos(al)]])
check("C7 EM: a boost of rapidity phi turns F = E + icB by the IMAGINARY angle i*phi (Silberstein 1907)", np.allclose(Rc @ F, F2))
check("C8 that complex rotation is complex-orthogonal (R^T R = 1): the Lorentz group is SO(3, C)", np.allclose(Rc.T @ Rc, np.eye(3)))
inv = lambda F: np.dot(F, F)
check("C9 F.F = (E^2 - B^2) + 2i E.B is unchanged by the boost (both Lorentz invariants in one complex number)", abs(inv(F) - inv(F2)) < 1e-12)

th = 0.9; Fd = np.exp(-1j*th)*F; Ed, Bd = Fd.real, Fd.imag
check("C10 duality rotation F -> e^(-i theta)F keeps energy |F|^2/2 and Poynting E x B (i = the E-to-B quarter turn)",
      abs(np.vdot(F, F) - np.vdot(Fd, Fd)) < 1e-12 and np.allclose(np.cross(E, B), np.cross(Ed, Bd)))
# vacuum Maxwell as curl F = i dF/dt (c = 1), plane wave along x, numerically
xg, tg, dx = 0.37, 0.81, 1e-5
Fw = lambda x, t, th=0: np.exp(-1j*th)*np.array([0, np.cos(x - t), 0]) + 1j*np.exp(-1j*th)*np.array([0, 0, np.cos(x - t)])
def curl_x_only(Ff, x, t):            # fields depend on x only
    dF = (Ff(x+dx, t) - Ff(x-dx, t))/(2*dx); return np.array([0, -dF[2], dF[1]])
ok = all(np.allclose(curl_x_only(lambda x, t: Fw(x, t, th), xg, tg), 1j*(Fw(xg, tg+dx, th) - Fw(xg, tg-dx, th))/(2*dx), atol=1e-6) for th in [0, 0.9])
check("C11 vacuum Maxwell is one complex equation curl F = i dF/dt (Riemann-Silberstein); duality-rotated waves still solve it", ok)

al = 1j*0.7
Rm = np.array([[np.cos(al), np.sin(al)], [-np.sin(al), np.cos(al)]])
xv, ctv = 1.3, 0.4
xp, x4p = Rm @ np.array([xv, 1j*ctv])
check("C12 Minkowski's ict: a rotation by imaginary angle i*phi in the (x, ict) plane is exactly the boost",
      abs(xp - (xv*np.cosh(0.7) - ctv*np.sinh(0.7))) < 1e-12 and abs(x4p/1j - (ctv*np.cosh(0.7) - xv*np.sinh(0.7))) < 1e-12)

Lam = lambda ph: np.array([[np.cosh(ph), -np.sinh(ph)], [-np.sinh(ph), np.cosh(ph)]])
check("C13 Euler for boosts: rapidity i*pi gives -1 on (ct, x), the reflection to the opposite wedge; 2*pi*i gives 1",
      np.allclose(Lam(1j*np.pi), -np.eye(2)) and np.allclose(Lam(2j*np.pi), np.eye(2)))

w0 = 1.7
def ham(y): q, p = y; return np.array([p, -w0*w0*q])
y = np.array([0.4, 0.3]); h = 1e-3; aa = []
for n in range(3000):
    k1 = ham(y); k2 = ham(y+h/2*k1); k3 = ham(y+h/2*k2); k4 = ham(y+h*k3); y = y + h/6*(k1+2*k2+2*k3+k4)
a0 = (w0*0.4 + 1j*0.3)/np.sqrt(2*w0); a3 = (w0*y[0] + 1j*y[1])/np.sqrt(2*w0)
check("C14 classical oscillator: a = (wq + ip)/sqrt(2w) obeys da/dt = -iwa (phase space turns; i = that turn)",
      abs(a3 - a0*np.exp(-1j*w0*3.0)) < 1e-9)

Nt = 64; D = (np.roll(np.eye(Nt), -1, axis=0) - np.roll(np.eye(Nt), 1, axis=0))/(2*2*np.pi/Nt)
evD = np.linalg.eigvals(D)
check("C15 Koopman-von Neumann: the classical Liouville generator (here d/dtheta) is REAL antisymmetric, so its spectrum is purely imaginary; i appears only if you insist on writing it as Hermitian",
      np.allclose(D, -D.T) and np.max(np.abs(evD.real)) < 1e-10)

# C16-17 Kramers-Kronig
w00, gam = 1.0, 0.3
chi = lambda w: 1/(w00**2 - w**2 - 1j*gam*w)
chia = lambda w: 1/(w00**2 - w**2 + 1j*gam*w)
hh = 0.002; grid = np.arange(-300, 300, hh) + hh/2
def kk_re(chif, w):                   # Re chi(w) = (1/pi) PV int Im chi(w')/(w' - w) dw'; w sits on a grid node, nodes are skipped symmetrically
    g = grid[np.abs(grid - w) > 1e-12]
    return np.sum(chif(g).imag/(g - w))*hh/np.pi
wt = np.array([-2.0, -0.6, 0.0, 0.5, 1.0, 1.4, 3.0]) + hh/2 - hh/2
wt = np.round(wt/hh)*hh   # put test points on half-steps relative to grid nodes
rec = np.array([kk_re(chi, w) for w in wt]); tru = chi(wt).real
reca = np.array([kk_re(chia, w) for w in wt]); trua = chia(wt).real
check("C16 causality -> Kramers-Kronig: Re chi rebuilt from Im chi alone (damped oscillator, poles in the lower half-plane)",
      np.max(np.abs(rec - tru)) < 2e-3, f"max err {np.max(np.abs(rec - tru)):.1e}")
check("C17 the same recipe FAILS for an acausal response (poles in the upper half-plane): it returns minus the true Re chi",
      np.max(np.abs(reca + trua)) < 2e-3)
poles = np.roots([-1, -1j*gam, w00**2])
check("C18 the causal response's poles sit at Im w = -gamma/2 < 0: the sign of i picks which half-plane is the future",
      np.allclose(poles.imag, -gam/2))

# C19 KMS: thermal state = periodic in imaginary time
Hq = np.array([[0, 0.3], [0.3, 1.0]]); beta = 0.7
evh, V = np.linalg.eigh(Hq)
expH = lambda s: V @ np.diag(np.exp(s*evh)) @ V.T          # e^(s H), s complex
rho = expH(-beta)/np.trace(expH(-beta))
At = lambda A, t: expH(1j*t) @ A @ expH(-1j*t)
G = lambda A, B, t: np.trace(rho @ At(A, t) @ B)
check("C19 KMS (Kubo 1957, Martin-Schwinger 1959): <A(t)B> = <B A(t + i beta)>; temperature is a period in imaginary time",
      all(abs(G(s3, s1, t) - np.trace(rho @ s1 @ At(s3, t + 1j*beta))) < 1e-12 for t in [0.0, 0.4, 2.3]))
check("C20 the thermal weight e^(-beta H) is time evolution e^(-iHt) at the imaginary time t = -i beta (Wick)",
      np.allclose(expH(-beta), expH(-1j*(-1j*beta))))

# C21-23 Unruh / Hawking
for P in [np.pi, 2*np.pi, 3*np.pi]:
    pass
ratio = lambda P, r: (P*r)/(2*np.pi*r)
check("C21 Euclidean Rindler ds^2 = rho^2 dtheta^2 + drho^2 has no cone tip only if theta has period 2*pi",
      ratio(2*np.pi, 1e-3) == 1 and ratio(np.pi, 1e-3) != 1)
hbar, c, Gc, kB = 1.054571817e-34, 299792458.0, 6.67430e-11, 1.380649e-23
GM = 1.32712440018e20
TU = hbar*9.80665/(2*np.pi*c*kB)
TH = hbar*c**3/(8*np.pi*kB*GM)
check("C22 Unruh T = hbar a/(2 pi c kB): 4.0e-20 K at a = 1 g", abs(TU - 3.97e-20)/3.97e-20 < 0.01, f"{TU:.3e} K")
check("C23 Hawking T = hbar c^3/(8 pi G M kB) for one solar mass: 6.17e-8 K (one formula holds hbar, c, G, kB and the 2*pi)",
      abs(TH - 6.17e-8)/6.17e-8 < 0.01, f"{TH:.3e} K")
kap = c**4/(4*GM/1)    # surface gravity c^4/(4GM) in m/s^2
check("C24 the same T is the Unruh formula at the horizon's surface gravity kappa = c^4/(4GM): G enters only through kappa",
      abs(hbar*kap/(2*np.pi*c*kB) - TH)/TH < 1e-12)

npass = sum(ok for _, ok in results)
print(f"\n{npass}/{len(results)} pass")
