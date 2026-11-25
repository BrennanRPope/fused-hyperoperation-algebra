"""Does the classical limit force the imaginary dial, or only allow it?  topics/why-i/e-i-pi/chain-resolves (Nov 25, 2026).

Harmonic oscillator, m = omega = 1. Real-time action S = int (x'^2/2 - x^2/2) dt.
Path weight exp(-kappa*S/hbar): kappa = 1 is the real dial (Laplace), kappa = -i the imaginary dial (Feynman's e^{iS/hbar}).
Checkerboard steps as on topics/why-i/checkerboard (k = -1 Kac, k = -i Feynman).
Run: python3 checks.py   (numpy only)
"""
import numpy as np
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def Scl(x, y, t):            # classical action x -> y in time t (m = omega = 1)
    return ((x*x + y*y)*np.cos(t) - 2*x*y)/(2*np.sin(t))
def SE(x, y, t):             # Euclidean classical action
    return ((x*x + y*y)*np.cosh(t) - 2*x*y)/(2*np.sinh(t))
def xcl(x0, x1, T, t):       # classical path
    return (x0*np.sin(T - t) + x1*np.sin(t))/np.sin(T)

x0, x1 = 0.0, 1.0
def Smid(xm, T): return Scl(x0, xm, T/2) + Scl(xm, x1, T/2)

# R1 reduced action in the middle point
ok = True
for T in [0.6*np.pi, 1.4*np.pi]:
    xs = xcl(x0, x1, T, T/2); h = 1e-4
    d2 = (Smid(xs + h, T) - 2*Smid(xs, T) + Smid(xs - h, T))/h**2
    d1 = (Smid(xs + h, T) - Smid(xs - h, T))/(2*h)
    ok &= abs(d2 - 2/np.tan(T/2)) < 1e-4 and abs(d1) < 1e-7
check("R1 the action, as a function of the middle point, is stationary at the classical path; its curvature is 2 cot(T/2), which turns negative past T = pi", ok)

# R2 real dial, short time: Gaussian peaked at the classical point, width ~ sqrt(hbar)
T = 0.6*np.pi; xs = xcl(x0, x1, T, T/2); g = np.linspace(xs - 3, xs + 3, 600001)
def width(hb):
    w = np.exp(-(Smid(g, T) - Smid(xs, T))/hb); w /= w.sum()
    mu = (g*w).sum(); return mu, np.sqrt(((g - mu)**2*w).sum())
m1, w1 = width(0.04); m2, w2 = width(0.01)
check("R2 real dial, T < pi: e^(-S/hbar) piles up on the classical middle point; width halves when hbar quarters (Laplace)",
      abs(m1 - xs) < 1e-6 and abs(w1/w2 - 2) < 1e-3, f"peak {m1:.6f} vs classical {xs:.6f}; widths {w1:.4f}, {w2:.4f}")

def window_integral(T, hb, lo, hi, kappa=-1j, n=2_000_001, taper=0.15):
    x = np.linspace(lo, hi, n); dx = x[1] - x[0]
    w = np.ones_like(x); L = taper*(hi - lo)
    a = x - lo < L; w[a] = np.sin(np.pi/2*(x[a] - lo)/L)**2
    b = hi - x < L; w[b] = np.sin(np.pi/2*(hi - x[b])/L)**2
    return np.sum(w*np.exp(-kappa*Smid(x, T)/hb))*dx
def fresnel(T, hb):
    A = 2/np.tan(T/2); xs = xcl(x0, x1, T, T/2)
    return np.sqrt(2*np.pi*hb/abs(A))*np.exp(1j*np.pi/4*np.sign(A))*np.exp(1j*Smid(xs, T)/hb)

hb = 0.01
for T, tag in [(0.6*np.pi, "R3"), (1.4*np.pi, "R4")]:
    xs = xcl(x0, x1, T, T/2)
    near = window_integral(T, hb, xs - 2.5, xs + 2.5); far = window_integral(T, hb, xs + 0.8, xs + 3.0)
    ex = fresnel(T, hb)
    extra = "" if tag == "R3" else " - still works past T = pi, and the phase has turned by a quarter (e^(-i pi/2))"
    check(f"{tag} imaginary dial, T = {T/np.pi:.1f} pi: a window round the classical point gives the whole Fresnel value; a window away from it gives ~0 (stationary phase)" + extra,
          abs(near - ex)/abs(ex) < 2e-3 and abs(far)/abs(ex) < 2e-3, f"|near-exact|/|exact| = {abs(near-ex)/abs(ex):.1e}, |far|/|exact| = {abs(far)/abs(ex):.1e}")

T = 1.4*np.pi; xs = xcl(x0, x1, T, T/2)
vals = [np.exp(-(Smid(xs + d, T) - Smid(xs, T))/0.01) for d in [0.1, 0.5, 1.0]]
check("R5 real dial, T > pi: e^(-S/hbar) GROWS away from the classical point (the classical path is a maximum along this direction), so the sum diverges",
      vals[0] < vals[1] < vals[2] and vals[2] > 1e10, f"weights at 0.1, 0.5, 1.0 away: {vals[0]:.1e}, {vals[1]:.1e}, {vals[2]:.1e}")

# R6-R8 discretized whole-path action: Hessian eigenvalues (Morse index)
def hessian(T, N=200):
    e = T/N; n = N - 1
    A = (np.diag(np.full(n, 2/e)) - np.diag(np.full(n - 1, 1/e), 1) - np.diag(np.full(n - 1, 1/e), -1)) - np.eye(n)*e
    return A
Ts = [0.5*np.pi, 0.9*np.pi, 1.1*np.pi, 1.9*np.pi, 2.1*np.pi, 2.9*np.pi, 3.1*np.pi]
negs = [int(np.sum(np.linalg.eigvalsh(hessian(t)) < 0)) for t in Ts]
check("R6 Morse index: the number of directions in which the classical path is a MAXIMUM of the action = number of half-periods passed (floor(T/pi))",
      negs == [int(np.floor(t/np.pi)) for t in Ts], f"T/pi = {[round(t/np.pi,1) for t in Ts]} -> {negs}")
tops = [(np.linalg.eigvalsh(hessian(t, 200)).max(), np.linalg.eigvalsh(hessian(t, 400)).max()) for t in Ts]
check("R7 there are always large positive directions (wiggly paths cost kinetic action; the top one grows like 4N/T as slices refine), so a dial with Re(kappa) < 0 diverges at every T",
      all(a > 0 and abs(b/a - 2) < 0.02 for a, b in tops))
def allowed_angles(T):
    lam = np.linalg.eigvalsh(hessian(T)); ang = np.linspace(-np.pi, np.pi, 3601)
    ok = [a for a in ang if np.min((np.exp(1j*a)*lam).real) >= -1e-12]
    return min(ok), max(ok), len(ok)
lo1, hi1, n1 = allowed_angles(0.6*np.pi); lo2, hi2, n2 = allowed_angles(1.4*np.pi)
check("R8 dials kappa = e^(i alpha) that keep the Gaussian path sum from blowing up: the whole half-plane |alpha| <= 90 deg for T < pi; ONLY alpha = +-90 deg (purely imaginary) for T > pi",
      abs(lo1 + np.pi/2) < 2e-3 and abs(hi1 - np.pi/2) < 2e-3 and n2 == 2 and abs(abs(lo2) - np.pi/2) < 1e-9,
      f"T=0.6pi: [{np.degrees(lo1):.1f}, {np.degrees(hi1):.1f}] deg; T=1.4pi: {n2} angles at {np.degrees(lo2):.1f}, {np.degrees(hi2):.1f} deg")
ph = [np.sum(np.sign(np.linalg.eigvalsh(hessian(t))))*np.pi/4 for t in Ts]
rel = [(p - ph[0])/(-np.pi/2) for p in ph]
check("R9 Maslov: at the imaginary dial each passed conjugate point turns the amplitude by a quarter turn, e^(-i pi/2) (i again, one per caustic)",
      np.allclose(rel, negs))

# R10-R11 Hamilton-Jacobi from either dial
h = 1e-5; x, y, t = 0.37, -0.81, 0.9
St = (Scl(y, x, t + h) - Scl(y, x, t - h))/(2*h); Sx = (Scl(y, x + h, t) - Scl(y, x - h, t))/(2*h)
check("R10 imaginary dial -> classical mechanics: the classical action obeys Hamilton-Jacobi dS/dt + (dS/dx)^2/2 + x^2/2 = 0 (Madelung/WKB as hbar -> 0)",
      abs(St + Sx**2/2 + x*x/2) < 1e-8)
St = (SE(y, x, t + h) - SE(y, x, t - h))/(2*h); Sx = (SE(y, x + h, t) - SE(y, x - h, t))/(2*h)
check("R11 real dial (Wick-rotated, where it converges for all times) -> dS/dtau + (dS/dx)^2/2 - x^2/2 = 0: Hamilton-Jacobi of the UPSIDE-DOWN potential (Hopf-Cole, Freidlin-Wentzell)",
      abs(St + Sx**2/2 - x*x/2) < 1e-8)
lamE = [np.linalg.eigvalsh(hessian(t) + 2*np.eye(199)*t/200).min() for t in [0.5, 3.0, 10.0]]
check("R12 the Euclidean action has no conjugate points (its Hessian stays positive at every duration): the real dial converges for all times, but on the inverted mechanics",
      all(l > 0 for l in lamE))

# R13 caustic / two branches: the real dial keeps the minimum, the imaginary dial keeps both
s1, s2, tt = -1.0, 1.0, 1.0
xg = np.linspace(-3, 3, 6001)
S1 = (xg - s1)**2/(2*tt); S2 = (xg - s2)**2/(2*tt); Smin = np.minimum(S1, S2)
errs = []
for hb in [0.1, 0.01, 0.001]:
    Sreal = -hb*np.log(np.exp(-(S1 - Smin)/hb) + np.exp(-(S2 - Smin)/hb)) + Smin
    errs.append(np.max(np.abs(Sreal - Smin)))
check("R13 two classical routes: the real dial's -hbar log(sum) tends to the MINIMUM of the two actions (one branch kept, the other forgotten; viscosity solution)",
      errs[-1] <= 0.001*np.log(2) + 1e-12 and errs[0] > errs[1] > errs[2])
hb = 0.01
I = np.abs(np.exp(1j*S1/hb) + np.exp(1j*S2/hb))**2
mid = (np.abs(xg) < 0.5)
zc = np.sum(np.diff(np.sign(np.cos((S1 - S2)[mid]/hb/2))) != 0)
spacing = 2*np.pi*hb*tt/(s2 - s1)
check("R14 the imaginary dial keeps both routes: |e^(iS1/h) + e^(iS2/h)|^2 swings between 0 and 4 with fringe spacing 2 pi hbar t/d (interference)",
      I[mid].max() > 3.99 and I[mid].min() < 0.01 and abs(1.0/zc*1.0 - spacing) < 0.2*spacing, f"spacing {spacing:.4f}; fringes in |x|<0.5: {zc}")

# R15-R17 checkerboard chain: reversibility
def step(psi, k, eps, mu):
    if np.isclose(k, -1): c_, t_ = 1 - eps*mu, eps*mu
    else: c_, t_ = np.cos(eps*mu), 1j*np.sin(eps*mu)
    turned = np.stack([c_*psi[:, 0] + t_*psi[:, 1], c_*psi[:, 1] + t_*psi[:, 0]], axis=1)
    return np.stack([np.roll(turned[:, 0], 1), np.roll(turned[:, 1], -1)], axis=1)
n, eps, mu, N = 401, 0.05, 1.0, 300
def run(k, start):
    psi = np.zeros((n, 2), complex); psi[start, 0] = 1; out = [psi]
    for _ in range(N): psi = step(psi, k, eps, mu); out.append(psi)
    return np.array(out)
U1, U2 = run(-1j, 200), run(-1j, 202)
K1, K2 = run(-1, 200), run(-1, 202)
dU = np.linalg.norm((U1 - U2).reshape(N + 1, -1), axis=1)
dK = np.abs(K1 - K2).reshape(N + 1, -1).sum(axis=1)
check("R15 checkerboard at k = -i (Feynman): two different starts stay exactly as distinguishable after 300 steps (nothing lost)",
      np.max(np.abs(dU - dU[0])) < 1e-12)
totK = np.abs(K1).reshape(N + 1, -1).sum(axis=1)
check("R16 checkerboard at k = -1 (Kac): total probability is kept, but the two starts blur together (distance falls)",
      np.max(np.abs(totK - 1)) < 1e-12 and dK[-1] < 0.5*dK[0], f"L1 distance {dK[0]:.2f} -> {dK[-1]:.3f}")
c_, t_ = 1 - eps*mu, eps*mu
Mk = np.array([[c_, t_], [t_, c_]]); Mu = np.array([[np.cos(eps*mu), 1j*np.sin(eps*mu)], [1j*np.sin(eps*mu), np.cos(eps*mu)]])
check("R17 Kac's turn step can't be run backward as a process (its inverse has negative entries); Feynman's can (inverse = conjugate transpose)",
      np.linalg.inv(Mk).min() < 0 and np.allclose(np.linalg.inv(Mu), Mu.conj().T))
p = K1[-1]; ent = -np.sum(np.abs(p)*np.log(np.abs(p) + 1e-300))
check("R18 the Kac chain's entropy rises step by step (information lost); the Feynman chain's step is invertible at every resolution eps",
      ent > 2.0 and all(abs(np.linalg.det(np.array([[np.cos(e), 1j*np.sin(e)], [1j*np.sin(e), np.cos(e)]]))) - 1 < 1e-12 for e in [0.5, 0.05, 0.005]))

npass = sum(ok for _, ok in results)
print(f"\n{npass}/{len(results)} pass")
