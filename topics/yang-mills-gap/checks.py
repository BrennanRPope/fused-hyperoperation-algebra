"""Yang–Mills mass gap through Affine Composition Theory: Pass 1 checks.

Run: python3 -I checks.py  ->  checks_output.txt, results.json
Each check prints PASS/FAIL. A polishing step re-runs the key numbers at higher
precision or finer grids and reports how far they moved.

Parts (fundamental to derived):
 A  the law: only the circle class (Delta < 0) closes up; gap = -Delta = 4 x curvature reading
 B  loops: the plaquette holonomy is the composition round a square; order gap = field strength's [A,A]
 C  one plaquette: the rung ladder (strong) turns into an oscillator (weak); never closes
 D  2D exact: string tension = log ratio of rung weights; zeros never reach real beta > 0
 E  twist: ladder shifted by a fractional rung; gap closes at half a rung (theta = pi)
 F  dial gas strips: a gap that closes at a continuous phase change (k_c = 3 + 2 sqrt 2)
 G  height flow: the gap at the end of the line, exp(-2/x) x^(1/3), zero to every order
 H  strong coupling: lowest glueball = tube of 4 plaquettes, m = -4 ln u
"""
import json, math, sys
import numpy as np
from scipy import integrate, special, linalg, optimize
import mpmath as mp

OUT = []
RES = {}
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    if ok:
        npass += 1
    else:
        nfail += 1
    line = f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else "")
    print(line)
    OUT.append(line)


def say(s):
    print(s)
    OUT.append(s)

# ---------------------------------------------------------------- A: the law
say("== A. Which dial laws close up into a circle (and so have a gap)?")
# law a o b = (a + b + k ab) / (1 - d ab); generator field v(x) = 1 + k x + d x^2; Delta = k^2 - 4d


def ruler_length(k, d, lo=-np.inf, hi=np.inf):
    f = lambda x: 1.0 / (1 + k * x + d * x * x)
    val, err = integrate.quad(f, lo, hi, limit=400)
    return val


def law(a, b, k, d):
    return (a + b + k * a * b) / (1 - d * a * b)


circle_settings = [(0.0, 1.0), (1.0, 1.0), (-2.0, 3.0), (3.0, 4.0), (0.5, 0.1)]
for k, d in circle_settings:
    D = k * k - 4 * d
    P = ruler_length(k, d)
    pred = 2 * math.pi / math.sqrt(-D)
    check(f"A1 circumference k={k}, d={d}: {P:.12f} vs 2pi/sqrt(-Delta)", abs(P - pred) < 1e-9 * pred,
          f"Delta={D}")
RES["A_circumferences"] = {f"{k},{d}": ruler_length(k, d) for k, d in circle_settings}


def L(x, k, d):
    """ruler (invariant length) from 0 to x, closed form for Delta < 0, on the branch through 0."""
    D = k * k - 4 * d
    s = math.sqrt(-D)
    return (2 / s) * (math.atan((2 * d * x + k) / s) - math.atan(k / s))

# the law becomes addition modulo the circumference
rng = np.random.default_rng(1)
worst = 0
for k, d in circle_settings:
    D = k * k - 4 * d
    P = 2 * math.pi / math.sqrt(-D)
    for _ in range(200):
        a, b = rng.normal(0, 2, 2)
        c = law(a, b, k, d)
        r = (L(c, k, d) - L(a, k, d) - L(b, k, d)) / P
        worst = max(worst, abs(r - round(r)))
check("A2 the law is addition of rulers modulo the circumference (1000 random pairs)", worst < 1e-9, f"worst {worst:.1e}")

# the characters chi_n = exp(2 pi i n L / P) are eigenfunctions of the invariant Laplacian (v d/dx)^2
worst = 0
for k, d in circle_settings:
    D = k * k - 4 * d
    P = 2 * math.pi / math.sqrt(-D)
    for n in [1, 2, 3]:
        for x0 in [-1.3, 0.2, 2.7]:
            h = 1e-3
            mp.mp.dps = 40
            v = lambda x: 1 + k * x + d * x * x
            chi = lambda x: mp.e ** (2j * mp.pi * n * mp.mpf(L(float(x), k, d)) / P) if False else None
            # use mpmath for the derivatives via closed form ruler
            s = mp.sqrt(-D)
            Lm = lambda x: (2 / s) * (mp.atan((2 * d * x + k) / s) - mp.atan(k / s))
            ch = lambda x: mp.exp(2j * mp.pi * n * Lm(x) / P)
            Df = lambda f: (lambda x: (1 + k * x + d * x * x) * mp.diff(f, x))
            lap = Df(Df(ch))(mp.mpf(x0))
            ratio = -lap / ch(mp.mpf(x0))
            worst = max(worst, abs(complex(ratio) - n * n * (-D)) / (n * n * (-D)))
check("A3 invariant Laplacian eigenvalues are n^2 (-Delta): gap = -Delta = 4 x curvature reading", worst < 1e-12,
      f"worst rel {worst:.1e}")

# line class and border: the ruler between walls is infinite (no closing up, continuous spectrum)
# border (k = 2, d = 1): v = (1 + x)^2, ruler from the wall is 1/eps - 1/11 exactly
lens = [1 / eps - 1 / 11.0 for eps in [1e-2, 1e-4, 1e-6]]
num = ruler_length(2.0, 1.0, -1 + 1e-2, 10.0)
check("A4 border (Delta=0): ruler up to the double wall is 1/eps - 1/11, grows without bound", abs(num - lens[0]) < 1e-6,
      f"quad {num:.4f} vs {lens[0]:.4f}")
for k, d, name in [(3.0, 1.0, "line"), (1.0, 0.0, "circle operation k=1")]:
    D = k * k - 4 * d
    walls = sorted(np.roots([d, k, 1]).real) if d != 0 else [-1 / k]
    w = walls[-1]
    lens = [ruler_length(k, d, w + eps, 10.0) for eps in [1e-2, 1e-4, 1e-6]]
    grows = lens[2] > lens[1] > lens[0] and (lens[2] - lens[1]) > 1.0
    check(f"A4 {name} (Delta={D}): ruler up to the wall grows without bound", grows,
          ", ".join(f"{x:.2f}" for x in lens))

# compactness forbids fractional rungs: exp(i nu L) single valued on the circle only for whole nu
k, d = 1.0, 1.0
P = 2 * math.pi / math.sqrt(3)
okw = abs(np.exp(2j * math.pi * 1.0) - 1) < 1e-12
okf = abs(np.exp(2j * math.pi * 0.5) - 1) > 1.9
check("A5 going once round the circle returns a whole rung to itself and sends a half rung to minus itself", okw and okf)

# ---------------------------------------------------------------- B: loops
say("== B. The plaquette: composition round a square")
# abelian circle-class links: holonomy = law-composition of signed link values; only total ruler matters
k, d = 1.0, 1.0


def inv(a, k, d):
    return -a / (1 + k * a)  # inverse under the law: a o inv(a) = 0 when d = 0; general below


def inv_general(a, k, d):
    # solve law(a, b) = 0 -> a + b + k a b = 0 -> b = -a / (1 + k a)
    return -a / (1 + k * a)

worst = 0
for _ in range(300):
    a1, a2, a3, a4 = rng.normal(0, 0.7, 4)
    h1 = law(law(law(a1, a2, k, d), inv_general(a3, k, d), k, d), inv_general(a4, k, d), k, d)
    h2 = law(law(law(inv_general(a4, k, d), a2, k, d), a1, k, d), inv_general(a3, k, d), k, d)
    worst = max(worst, abs(h1 - h2) / (1 + abs(h1)))
check("B1 abelian plaquette: order of the four links does not matter", worst < 1e-9, f"{worst:.1e}")

# non-abelian: su(2) links U = exp(a A_mu(x)); log of plaquette / a^2 -> F = dA - dA + [A, A]
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
T = [1j * s / 2 for s in (sx, sy, sz)]


def Afield(x, y, mu, coupled=True):
    # a smooth su(2) potential; components along different generators so [A_x, A_y] != 0
    if mu == 0:
        return (0.7 + 0.3 * y) * T[0] + 0.2 * x * T[2]
    return (0.5 - 0.4 * x) * T[1] + (0.1 * y * T[2] if coupled else 0)


def plaquette(x, y, a):
    U1 = linalg.expm(a * Afield(x + a / 2, y, 0))
    U2 = linalg.expm(a * Afield(x + a, y + a / 2, 1))
    U3 = linalg.expm(a * Afield(x + a / 2, y + a, 0))
    U4 = linalg.expm(a * Afield(x, y + a / 2, 1))
    return U1 @ U2 @ np.linalg.inv(U3) @ np.linalg.inv(U4)


def Fexact(x, y):
    h = 1e-6
    dAy = (Afield(x + h, y, 1) - Afield(x - h, y, 1)) / (2 * h)
    dAx = (Afield(x, y + h, 0) - Afield(x, y - h, 0)) / (2 * h)
    Ax, Ay = Afield(x, y, 0), Afield(x, y, 1)
    return dAy - dAx + (Ax @ Ay - Ay @ Ax), Ax @ Ay - Ay @ Ax

errs = []
for a in [0.1, 0.05, 0.025]:
    x0, y0 = 0.3, -0.2
    lg = linalg.logm(plaquette(x0, y0, a)) / a ** 2
    F, comm = Fexact(x0 + a / 2, y0 + a / 2)
    errs.append(np.abs(lg - F).max())
check("B2 log(plaquette)/a^2 -> dA_y - dA_x + [A_x, A_y] at the centre, error shrinks in step with a", 1.8 < errs[0] / errs[1] < 2.2 and 1.8 < errs[1] / errs[2] < 2.2,
      ", ".join(f"{e:.1e}" for e in errs))
F, comm = Fexact(0.3, -0.2)
check("B3 the [A_x, A_y] piece (the order gap / holonomy reading) is nonzero here and is a fixed part of F",
      np.abs(comm).max() > 0.05, f"|[A,A]| = {np.abs(comm).max():.3f} of |F| = {np.abs(F).max():.3f}")
RES["B_errors"] = errs

# ---------------------------------------------------------------- C: one plaquette
say("== C. One plaquette, quantum: rung ladder at strong coupling, oscillator at weak")
# U(1): H = 2 g^2 p^2 + (1/g^2)(1 - cos theta) in the rung basis n
NMAX = 60


def u1_plaquette_levels(g, nmax=NMAX):
    n = np.arange(-nmax, nmax + 1)
    H = np.diag(2 * g * g * n ** 2 + 1 / g ** 2) - (np.eye(len(n), k=1) + np.eye(len(n), k=-1)) / (2 * g * g)
    return np.linalg.eigvalsh(H)

# SU(2): characters chi_j, j = 0, 1/2, 1, ...; H = (g^2/2) 4 j(j+1) + (1/g^2)(1 - chi_{1/2}/2)


def su2_plaquette_levels(g, jmax2=120):
    js = np.arange(0, jmax2 + 1) / 2
    H = np.diag(2 * g * g * js * (js + 1) + 1 / g ** 2)
    off = -1 / (2 * g * g)
    for i in range(len(js) - 1):
        H[i, i + 1] = H[i + 1, i] = off
    return np.linalg.eigvalsh(H)

gs = np.geomspace(0.15, 4, 60)
gapU = [np.diff(u1_plaquette_levels(g)[:2])[0] for g in gs]
gapS = [np.diff(su2_plaquette_levels(g)[:2])[0] for g in gs]
RES["C"] = {"g": gs.tolist(), "gap_u1": gapU, "gap_su2": gapS}
g = 4.0
check("C1 U(1) strong coupling: gap -> 2 g^2 x (rung 1 - rung 0) = 2 g^2", abs(np.diff(u1_plaquette_levels(g)[:2])[0] / (2 * g * g) - 1) < 0.01,
      f"{np.diff(u1_plaquette_levels(g)[:2])[0]:.4f} vs {2*g*g}")
g = 0.15
check("C2 U(1) weak coupling: gap -> 2 (oscillator, sqrt(2 x 2g^2 x 1/g^2))", abs(np.diff(u1_plaquette_levels(g)[:2])[0] - 2) < 0.05,
      f"{np.diff(u1_plaquette_levels(g)[:2])[0]:.4f}")
g = 4.0
check("C3 SU(2) strong coupling: gap -> 2 g^2 x Casimir(1/2) = 1.5 g^2", abs(np.diff(su2_plaquette_levels(g)[:2])[0] / (1.5 * g * g) - 1) < 0.01)
check("C4 the gap never closes along the whole coupling line (one plaquette); it dips a little in the crossover", min(gapU) > 1.5 and min(gapS) > 1.0,
      f"min U(1) {min(gapU):.3f}, min SU(2) {min(gapS):.3f}")
# polish: truncation
g = 0.15
check("C5 polish: rung cut 60 -> 120 moves the weak gap by < 1e-10",
      abs(np.diff(u1_plaquette_levels(g, 60)[:2])[0] - np.diff(u1_plaquette_levels(g, 120)[:2])[0]) < 1e-10)

# ---------------------------------------------------------------- D: 2D exact
say("== D. Two dimensions, exact: rung weights, string tension, zeros")


def c_u1(n, beta):
    return special.iv(n, beta)


def c_su2_num(j, beta):
    # Haar on class angle: (2/pi) sin^2 t dt, chi_j = sin((2j+1)t)/sin t, weight exp(beta cos t)
    f = lambda t: (2 / math.pi) * math.sin(t) ** 2 * (math.sin((2 * j + 1) * t) / math.sin(t)) * math.exp(beta * math.cos(t))
    return integrate.quad(f, 0, math.pi, limit=200)[0]


def c_su2(j, beta):
    return 2 * (2 * j + 1) * special.iv(2 * j + 1, beta) / beta

worst = max(abs(c_su2_num(j, b) - c_su2(j, b)) / c_su2(j, b) for j in [0, .5, 1, 1.5] for b in [0.5, 2, 6])
check("D1 SU(2) rung weights c_j = 2(2j+1) I_{2j+1}(beta)/beta (vs Haar integral)", worst < 1e-8, f"{worst:.1e}")
xs = np.linspace(0.01, 60, 3000)
check("D2 on real beta > 0 rung 0 outweighs every other rung (U(1): I_0 > I_n)",
      all((special.iv(0, xs) > special.iv(n, xs)).all() for n in range(1, 8)))
check("D3 same for SU(2) per dimension: c_0 > c_j/(2j+1)",
      all((c_su2(0, xs) > c_su2(j, xs) / (2 * j + 1)).all() for j in [.5, 1, 1.5, 2]))
# string tension sigma_R = -ln(c_R/(d_R c_0)); ratio sigma_n/sigma_1 from n (strong) to n^2 (weak)
betas = np.geomspace(0.05, 200, 80)
ratio2 = [math.log(special.iv(0, b) / special.iv(2, b)) / math.log(special.iv(0, b) / special.iv(1, b)) for b in betas]
RES["D_ratio2"] = {"beta": betas.tolist(), "ratio": ratio2}
b = 1e-6
r_strong = math.log(special.iv(0, b) / special.iv(2, b)) / math.log(special.iv(0, b) / special.iv(1, b))
pred = 2 + math.log(2) / math.log(2 / b)
check("D4 U(1): sigma_2/sigma_1 -> 2 at strong coupling (counting), slowly: 2 + ln2/ln(2/beta)", abs(r_strong - pred) < 1e-3 and abs(ratio2[0] - (2 + math.log(2) / math.log(2 / 0.05))) < 0.01,
      f"beta=1e-6: {r_strong:.4f} vs {pred:.4f}; beta=0.05: {ratio2[0]:.4f}")
check("D5 ... and -> 4 = 2^2 at weak coupling (Casimir rung ladder)", abs(ratio2[-1] - 4) < 0.05, f"{ratio2[-1]:.3f}")
sig_su2 = lambda j, b: -math.log(c_su2(j, b) / ((2 * j + 1) * c_su2(0, b)))
r = sig_su2(1, 200) / sig_su2(0.5, 200)
check("D6 SU(2) weak coupling: sigma_1/sigma_1/2 -> C2(1)/C2(1/2) = 8/3", abs(r - 8 / 3) < 0.03, f"{r:.4f}")
b = 400.0
check("D7 U(1) sigma_1 ~ 1/(2 beta) at weak coupling: closes only as beta -> infinity",
      abs(math.log(special.iv(0, b) / special.iv(1, b)) * 2 * b - 1) < 0.01)

# Fisher zeros of the torus partition function Z_N = sum_n I_n(beta)^N (U(1), N plaquettes)


def Z_u1(beta, N, nmax=25):
    return sum((1 if n == 0 else 2) * special.iv(n, beta) ** N for n in range(nmax))


def find_zeros(N, re=(0.2, 14), im=(0.2, 14), grid=220):
    X, Y = np.meshgrid(np.linspace(*re, grid), np.linspace(*im, grid))
    B = X + 1j * Y
    # work with Z / I_0^N to avoid overflow
    W = 1 + sum(2 * (special.iv(n, B) / special.iv(0, B)) ** N for n in range(1, 25))
    A = np.abs(W)
    zs = []
    for i in range(1, grid - 1):
        for j in range(1, grid - 1):
            if A[i, j] == A[i - 1:i + 2, j - 1:j + 2].min() and A[i, j] < 0.3:
                f = lambda z: 1 + sum(2 * (mp.besseli(n, z) / mp.besseli(0, z)) ** N for n in range(1, 25))
                try:
                    z = complex(mp.findroot(f, mp.mpc(B[i, j])))
                    if re[0] < z.real < re[1] and im[0] * 0.5 < z.imag < im[1] and all(abs(z - w) > 1e-6 for w in zs):
                        zs.append(z)
                except Exception:
                    pass
    return zs

mp.mp.dps = 30
zeros = {}
for N in [8, 16, 32]:
    zs = find_zeros(N)
    zeros[N] = zs
    say(f"   N={N}: {len(zs)} zeros found in the box; smallest |Im beta| = {min(abs(z.imag) for z in zs):.4f}")
RES["D_zeros"] = {N: [[z.real, z.imag] for z in zs] for N, zs in zeros.items()}
# BKW: zeros accumulate where |I_1/I_0| = 1 ; on that curve Im beta stays away from 0 for finite Re beta
check("D8 zeros found for N = 8, 16, 32 and none on the real axis",
      all(len(zeros[N]) > 3 and min(abs(z.imag) for z in zeros[N]) > 0.1 for N in zeros))


def equimod_gap(z):
    return abs(special.iv(1, z) / special.iv(0, z)) - 1

dev = []
for N in [8, 16, 32]:
    dev.append(np.median([abs(equimod_gap(z)) for z in zeros[N]]))
check("D9 zeros crowd onto the curve |c_1| = |c_0| (top two rungs tie) as N grows", dev[2] < dev[1] < dev[0],
      ", ".join(f"{x:.3f}" for x in dev))

# where can the top two rungs tie? Only near the zeros of I_0 on the imaginary axis (beta = i j_{0,m})
Xg, Yg = np.meshgrid(np.linspace(0.01, 40, 800), np.linspace(-1, 1, 81))
ratio = np.abs(special.iv(1, Xg + 1j * Yg) / special.iv(0, Xg + 1j * Yg))
check("D10 in a band |Im beta| <= 1 round the whole positive real axis rung 0 is strictly the largest (max |I_1/I_0| < 1)",
      ratio.max() < 1, f"max {ratio.max():.4f}")
edge = []
for Y in [5, 10, 20, 40, 80, 160]:
    xs_ = np.linspace(0.01, 12, 4000); ys_ = np.linspace(Y - 3.2, Y + 3.2, 400)
    Xe, Ye = np.meshgrid(xs_, ys_)
    re_ = np.abs(special.iv(1, Xe + 1j * Ye) / special.iv(0, Xe + 1j * Ye))
    m = (re_ >= 1).any(axis=0)
    edge.append([Y, float(xs_[m].max())])
slope = np.polyfit(np.log([e[0] for e in edge]), [e[1] for e in edge], 1)[0]
RES["D_tie_edge"] = edge
check("D10b the tie region hugs the imaginary axis: its right edge creeps out only like ln|Im beta| (tie needs e^(-2 Re beta) ~ Re beta / |beta|^2)",
      0.6 < slope < 1.05, f"edge at Im=5..160: {', '.join(f'{e[1]:.2f}' for e in edge)}; slope vs ln Im = {slope:.3f}")
j0 = special.jn_zeros(0, 6)
allz = [complex(*z) for z in RES["D_zeros"][32]]
dist = [min(abs(z - 1j * j) for j in j0) for z in allz]
check("D11 every zero found (N=32) sits within 2.5 of a zero of rung 0's weight, i j_{0,m} = 2.40i, 5.52i, 8.65i, ...",
      max(dist) < 2.5 and all(z.real < 2.5 for z in allz), f"max distance {max(dist):.2f}")
RES["D_j0"] = j0.tolist()

# ---------------------------------------------------------------- E: twist
say("== E. A twist shifts the ladder by a fractional rung; the gap closes at half a rung")


def twisted_levels(theta, nmax=20):
    n = np.arange(-nmax, nmax + 1)
    return np.sort((n + theta / (2 * math.pi)) ** 2)

ths = np.linspace(0, 2 * math.pi, 401)
gaps = [twisted_levels(t)[1] - twisted_levels(t)[0] for t in ths]
RES["E"] = {"theta": ths.tolist(), "gap": gaps}
check("E1 gap = 1 - theta/pi for theta in [0, pi] (units g^2 L / 2)",
      max(abs(gaps[i] - (1 - ths[i] / math.pi)) for i in range(201)) < 1e-12)
check("E2 gap closes exactly at theta = pi (half a rung) and reopens after", gaps[200] < 1e-12 and gaps[150] > 0.2 and gaps[250] > 0.2)
# each of the two lowest levels still has its own excitation gap at theta = pi
lv = twisted_levels(math.pi)
check("E3 at theta = pi the two lowest levels tie, but the next level sits 2 rungs-worth above them (each vacuum stays gapped)",
      abs(lv[0] - lv[1]) < 1e-12 and abs(lv[2] - lv[0] - 2.0) < 1e-12, f"levels {lv[:4]}")
# zeros in complex theta pinch onto theta = pi: Z_N = sum_n exp(-N (n + theta/2pi)^2 / (2 beta))


def Ztw(theta, N, beta=1.0, nmax=8):
    return sum(mp.e ** (-N * (n + theta / (2 * mp.pi)) ** 2 / (2 * beta)) for n in range(-nmax, nmax + 1))

pinch = []
for N in [10, 20, 40, 80]:
    # nearest zero to theta = pi on the line Re theta = pi
    z = mp.findroot(lambda t: Ztw(t, N), mp.mpc(math.pi, 2 * math.pi * 2 * math.pi / N))
    pinch.append(complex(z))
RES["E_pinch"] = [[z.real, z.imag] for z in pinch]
check("E4 nearest zero sits on Re theta = pi and its height falls like 1/N (pinch onto the real axis)",
      all(abs(z.real - math.pi) < 1e-8 for z in pinch) and abs(pinch[0].imag / pinch[-1].imag - 8) < 0.1,
      ", ".join(f"{z.imag:.4f}" for z in pinch))
say("   (exact: Im theta = 2 pi^2 beta / N x (2/1)... fitted ratio " + f"{pinch[0].imag/pinch[-1].imag:.4f})")

# ---------------------------------------------------------------- F: dial gas strips
say("== F. Our dial gas on strips: a gap that closes at a continuous phase change")
from scipy.sparse.linalg import eigsh


def dialgas_gap(k, W, z=None):
    """row transfer matrix of the dial gas Z = sum z^|S| k^(edges inside S) on a cylinder of width W.
    On the symmetric line z = k^-2 it is the Ising model with K = ln(k)/4."""
    if z is None:
        z = k ** -2
    S = np.array([[(s >> i) & 1 for i in range(W)] for s in range(2 ** W)], float)
    nocc = S.sum(1)
    horiz = (S * np.roll(S, -1, axis=1)).sum(1)
    Dvec = z ** (nocc / 2) * k ** (horiz / 2)        # split the row weight between neighbouring rows
    V = k ** (S @ S.T)                                # vertical bonds between rows
    Tm = Dvec[:, None] * V * Dvec[None, :]
    ev = np.sort(np.linalg.eigvalsh(Tm))[::-1]
    return math.log(ev[0] / ev[1]), ev

kc = 3 + 2 * math.sqrt(2)
Ws = [4, 6, 8, 10]
fss = {}
for k in [kc * 0.8, kc, kc * 1.25]:
    fss[k] = [dialgas_gap(k, W)[0] for W in Ws]
RES["F"] = {"W": Ws, "kc": kc, "gaps": {f"{k:.4f}": v for k, v in fss.items()}}
wg = [W * g for W, g in zip(Ws, fss[kc])]
check("F1 at k_c = 3 + 2 sqrt 2 the gap x width is nearly constant (gap closes like 1/W: no gap in the infinite plane)",
      max(wg) - min(wg) < 0.06 * np.mean(wg), ", ".join(f"{x:.4f}" for x in wg))
check("F2 at k_c, W x gap approaches pi/4 (the conformal value for this model's first excitation)", abs(wg[-1] - math.pi / 4) < 0.02,
      f"{wg[-1]:.4f} vs {math.pi/4:.4f}")
lo = fss[kc * 0.8]
check("F3 below k_c (disordered side) the gap settles at a positive value as W grows", abs(lo[-1] - lo[-2]) < 0.02 and lo[-1] > 0.2,
      ", ".join(f"{x:.4f}" for x in lo))
hi = [dialgas_gap(kc * 1.6, W)[0] for W in Ws + [12]]
check("F4 well above k_c (k = 1.6 k_c, ordered side) the top two levels merge exponentially: two ordered states, not a particle",
      all(hi[i + 1] < 0.5 * hi[i] for i in range(len(hi) - 1)), ", ".join(f"{x:.1e}" for x in hi))
g2 = {}
for W in [8, 10, 12]:
    _, ev = dialgas_gap(kc * 1.6, W)
    g2[W] = math.log(ev[0] / ev[2])
ext1 = (100 * g2[10] - 64 * g2[8]) / 36
ext2 = (144 * g2[12] - 100 * g2[10]) / 44
check("F5 ... while the gap above that pair stays open: a + b/W^2 fits agree, limit about 0.9",
      abs(ext1 - ext2) < 0.03 and ext2 > 0.5, f"from W=8,10: {ext1:.3f}; from W=10,12: {ext2:.3f}")
RES["F"]["ordered"] = {"two_lowest": hi, "excitation_extrapolated": [ext1, ext2]}
# polish: W = 12 at k_c
g12c = dialgas_gap(kc, 12)[0]
check("F6 polish: W = 12 at k_c keeps W x gap within 1% of W = 10", abs(12 * g12c - wg[-1]) < 0.01 * wg[-1], f"{12*g12c:.4f}")
RES["F"]["W12_kc"] = 12 * g12c

# ---------------------------------------------------------------- G: height flow
say("== G. The end of the line: the height flow's gap, exp(-2/x) x^(1/3), zero to every order")
mp.mp.dps = 50


ABEL_C = [mp.mpf(-1)/36, mp.mpf(1)/540, mp.mpf(1)/7776, mp.mpf(-71)/435456, mp.mpf(8759)/163296000, mp.mpf(31)/20995200]


def abel(x, n=3000):
    """Abel (Fatou) coordinate of E(x) = e^x - 1 on x > 0: alpha(E(x)) = alpha(x) + 1."""
    y = mp.mpf(x)
    for _ in range(n):
        y = mp.log(1 + y)
    # asymptotic series (coefficients solved with sympy): -2/y + (1/3) log y + sum c_i y^i
    return -2 / y + mp.log(y) / 3 + sum(c * y ** (i + 1) for i, c in enumerate(ABEL_C)) + n


def heights_to_one(x0):
    return abel(1) - abel(x0)

xs0 = [mp.mpf(1) / m for m in (4, 8, 16, 32, 64, 128, 256)]
const = [mp.e ** (-heights_to_one(x)) * mp.e ** (2 / x) * x ** (-mp.mpf(1) / 3) for x in xs0]
check("G1 the number of heights to reach strong coupling is 2/x0 - (1/3) ln x0 + C: e^{-s*} e^{2/x0} x0^{-1/3} settles (drift ~ x0)",
      abs(const[-1] - const[-2]) < 1e-3 * abs(const[-1]) and abs(const[-1]-const[-2]) < 0.6 * abs(const[-2]-const[-3]), ", ".join(mp.nstr(c, 8) for c in const))
RES["G_const"] = [float(c) for c in const]
# abel equation holds
check("G2 Abel coordinate obeys alpha(e^x - 1) = alpha(x) + 1", abs(abel(mp.e ** mp.mpf(0.3) - 1) - abel(0.3) - 1) < 1e-20)
# zero to every order: m(x)/x^p -> 0 for every p
mvals = [(x, mp.e ** (-heights_to_one(x))) for x in xs0]
p = 20
rr = [m / x ** p for x, m in mvals]
check("G3 the gap m(x0) = e^{-s*} beats every power: m/x0^20 falls toward 0 once x0 < 1/32",
      rr[-1] < rr[-2] < rr[-3] and rr[-1] < 1e-30, ", ".join(mp.nstr(v, 3) for v in rr[-3:]))
# continuum line: zoom e^s while letting x run keeps m x zoom fixed (invariance along the flow)
x = mp.mpf("0.05")
m0 = mp.e ** (-heights_to_one(x))
x1 = mp.log(1 + x)  # one height back toward the UV
m1 = mp.e ** (-heights_to_one(x1)) * mp.e  # zoom in by e^1
check("G4 stepping one height toward weak coupling and zooming by e leaves the gap fixed (the continuum line)", abs(m1 / m0 - 1) < 1e-18, mp.nstr(abs(m1 / m0 - 1), 3))
# sign: with x < 0 the flow goes to zero at large heights: no growth, no gap from this mechanism
y = mp.mpf(-0.5)
for _ in range(50):
    y = mp.e ** y - 1
check("G5 the other sign (x < 0) flows to zero coupling: weak at long distance, no gap is generated", -0.05 < y < 0)

# ---------------------------------------------------------------- H: strong coupling glueball
say("== H. Strong coupling in 3 or 4 dimensions: lowest glueball = a tube of 4 plaquettes")
# smallest surface joining two plaquettes at time distance t: a 1x1 tube of 4t plaquettes -> correlation ~ u^{4t}
for name, u in [("U(1)", lambda b: special.iv(1, b) / special.iv(0, b)), ("SU(2)", lambda b: c_su2(.5, b) / (2 * c_su2(0, b)))]:
    b = 0.3
    m = -4 * math.log(u(b))
    say(f"   {name}: beta = {b}: u = {u(b):.5f}, m = -4 ln u = {m:.4f} (rung-1 weight cost x 4 plaquettes)")
# check the counting with a small transfer-matrix: the tube count of minimal surfaces between two plaquettes at distance t in Z^3
# (enumerate closed surfaces of area 4t + 2 containing two given parallel plaquettes: only the straight tube)
say("   (counting argument, not a numerical check; the Hamiltonian version is C1/C3: one plaquette = 4 links in rung 1)")

say(f"\nTOTAL: {npass} pass, {nfail} fail")
with open("checks_output.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
with open("results.json", "w") as f:
    json.dump(RES, f, default=lambda o: str(o), indent=1)
