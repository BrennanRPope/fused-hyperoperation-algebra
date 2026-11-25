"""DRAFT stress test of the clock reading x = 1 - e^{-t} and of gluing by o_{-1}:
heated, rotating ('pressure'), chemical-potential, heavy-particle and local-quench states.
Run: python3 checks2.py > checks2_output.txt  (writes results2.json)"""
import json, math
import numpy as np
from lib import circ, haar_state, rho_of, ptrace, vn

rng = np.random.default_rng(1125)
R = {}; passed = failed = 0
def check(name, ok, detail=""):
    global passed, failed
    ok = bool(ok); passed += ok; failed += (not ok)
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""))

def S_from_C(C):
    w = np.clip(np.linalg.eigvalsh(C), 1e-15, 1 - 1e-15)
    return float(-(w * np.log(w) + (1 - w) * np.log(1 - w)).sum())
def toeplitz_C(nk_func, Lmax, nk=400000):
    k = (np.arange(nk) + 0.5) / nk * 2 * np.pi - np.pi
    n = nk_func(k); d = np.arange(Lmax + 1)
    return np.array([np.mean(n * np.exp(1j * k * dd)) for dd in d])
def block_S(c, sites):
    s = np.asarray(sites); D = s[:, None] - s[None, :]
    C = np.where(D >= 0, c[np.abs(D)], np.conj(c[np.abs(D)]))
    return S_from_C(C)
def cmi(Sf, la, lb, lc):
    A = list(range(0, la)); B = list(range(la, la + lb)); C = list(range(la + lb, la + lb + lc))
    return Sf(A + B) + Sf(B + C) - Sf(B) - Sf(A + B + C)
def t_line(la, lb, lc, beta=None):
    f = (lambda l: l) if beta is None else (lambda l: math.sinh(math.pi * l / beta))
    return math.log(f(la + lb) * f(lb + lc) / (f(lb) * f(la + lb + lc)))
CONFIGS = [(40, 20, 40), (60, 30, 50), (50, 50, 50), (30, 15, 80), (80, 40, 30)]
fermi = lambda e: 1 / (1 + np.exp(np.clip(e, -700, 700)))

# ---------------------------------------------------------------- 1. heated state (thermal, infinite line)
print("== 1. heated state: free fermions at temperature 1/beta (c = 1, velocity 2 for hopping 1)")
heat = {}
for beta in [15.0, 40.0]:
    c = toeplitz_C(lambda k: fermi(-2 * beta * np.cos(k)), 260)
    Sf = lambda s: block_S(c, s)
    rows = []
    for la, lb, lc in CONFIGS:
        I = cmi(Sf, la, lb, lc); tb = t_line(la, lb, lc, 2 * beta); tv = t_line(la, lb, lc)
        rows.append((la, lb, lc, I, tb / 3, tv / 3))
    heat[beta] = rows
    for r in rows: print(f"   beta={beta}  {r[:3]}  CMI={r[3]:.4f}  t_thermal/3={r[4]:.4f}  t_vacuum/3={r[5]:.4f}")
    check(f"H1 beta={beta}: CMI = t_beta/3 with t from the thermal (sinh) cross-ratio; the clock survives heating",
          all(abs(r[3] - r[4]) < 0.01 for r in rows), f"max diff {max(abs(r[3]-r[4]) for r in rows):.4f}")
    check(f"H2 beta={beta}: the vacuum clock would be wrong (heating changes which clock, not whether there is one)",
          max(abs(r[3] - r[5]) for r in rows) > 0.03, f"max vacuum miss {max(abs(r[3]-r[5]) for r in rows):.3f}")
R["heat"] = {str(b): v for b, v in heat.items()}
ok = True
for _ in range(500):
    beta = rng.uniform(.5, 5); w = np.sort(rng.uniform(-3, 3, 4)); a1, a2, b1, b2 = w
    xs = math.sinh(math.pi*(a2-a1)/beta)*math.sinh(math.pi*(b2-b1)/beta)/(math.sinh(math.pi*(b1-a1)/beta)*math.sinh(math.pi*(b2-a2)/beta))
    z = np.exp(2 * math.pi * w / beta); xp = (z[1]-z[0])*(z[3]-z[2])/((z[2]-z[0])*(z[3]-z[1]))
    ok &= abs(xs - xp) < 1e-9 * max(1, xs)
    zeta = lambda u: math.log(abs(math.sinh(math.pi*(u-a1)/beta)/math.sinh(math.pi*(u-a2)/beta)))
    ok &= abs(-math.log(1 - xs) - abs(zeta(b2) - zeta(b1))) < 1e-8
check("H3 thermal cross-ratio = vacuum cross-ratio after w -> e^{2 pi w/beta}; t = width in the interval's thermal modular time (Cardy-Tonni weight)", ok)

# ---------------------------------------------------------------- 2. 'pressure' states
print("== 2. 'pressure' states: rotating / boosted (beta_L != beta_R) and chemical potential")
bL, bR = 15.0, 60.0
c = toeplitz_C(lambda k: np.where(k > 0, fermi(-2 * bR * np.cos(k)), fermi(-2 * bL * np.cos(k))), 260)
Sf = lambda s: block_S(c, s); rot = []
for la, lb, lc in CONFIGS:
    I = cmi(Sf, la, lb, lc); tL = t_line(la, lb, lc, 2 * bL); tR = t_line(la, lb, lc, 2 * bR)
    rot.append((la, lb, lc, I, (tL + tR) / 6, tL / 3, tR / 3))
for r in rot: print(f"   {r[:3]}  CMI={r[3]:.4f}  (t_L+t_R)/6={r[4]:.4f}  t_L/3={r[5]:.4f}  t_R/3={r[6]:.4f}")
check("P1 SPINNING state: rotating/boosted (left movers 4x hotter than right movers): CMI = (t_L + t_R)/6, two clocks at half rate each",
      all(abs(r[3] - r[4]) < 0.01 for r in rot), f"max diff {max(abs(r[3]-r[4]) for r in rot):.4f}")
ok = True
for la, lb, lc in CONFIGS:
    tL = t_line(la, lb, lc, 2 * bL); tR = t_line(la, lb, lc, 2 * bR)
    x = 1 - math.exp(-(tL + tR) / 2); yL, yR = 1 - math.exp(-tL / 2), 1 - math.exp(-tR / 2)
    ok &= abs(x - circ(yL, yR, -1)) < 1e-12 and abs((1 - x) - math.sqrt((1 - (1 - math.exp(-tL))) * (1 - (1 - math.exp(-tR))))) < 1e-12
check("P2 so x = y_L o_{-1} y_R: the link probability is the union of a left-mover clock and a right-mover clock, each at half rate", ok)
R["rot"] = rot
mu_rows = []
for kF in [0.3 * math.pi, 0.4 * math.pi]:
    mu = -2 * math.cos(kF); v = 2 * math.sin(kF)
    for beta in [None, 40.0]:
        f = (lambda k: (np.abs(k) < kF).astype(float)) if beta is None else (lambda k, b=beta: fermi(b * (-2 * np.cos(k) - mu)))
        c = toeplitz_C(f, 330); Sf = lambda s: block_S(c, s)
        for la, lb, lc in [(60, 40, 60), (80, 60, 70), (100, 50, 100)]:
            I = cmi(Sf, la, lb, lc); t = t_line(la, lb, lc, None if beta is None else v * beta)
            mu_rows.append((kF / math.pi, beta or 0, la, lb, lc, I, t / 3))
for r in mu_rows: print(f"   kF={r[0]:.1f}pi beta={r[1]}  {r[2:5]}  CMI={r[5]:.4f}  t/3={r[6]:.4f}")
check("P3 chemical potential (filling 30%, 40%): CMI = t/3 with the thermal length set by the Fermi velocity; mu itself drops out (spectral flow)",
      all(abs(r[5] - r[6]) < 0.02 for r in mu_rows), f"max diff {max(abs(r[5]-r[6]) for r in mu_rows):.4f} (Friedel ripples)")
R["mu"] = mu_rows

# ---------------------------------------------------------------- 3. black hole on a circle (holographic BTZ): the plateau
print("== 3. black hole on a circle (holographic BTZ, units c/3): the entanglement plateau")
rp = 3.0                                          # horizon radius; beta = 2 pi / rp
f_conn = lambda th: math.log(math.sinh(rp * th / 2))
S_btz = lambda th: min(f_conn(th), f_conn(2 * math.pi - th) + math.pi * rp)
from math import isclose
lo, hi = math.pi, 2 * math.pi - 1e-9
for _ in range(200):
    mid = (lo + hi) / 2
    if f_conn(mid) < f_conn(2 * math.pi - mid) + math.pi * rp: lo = mid
    else: hi = mid
th_star = lo
def t_circ(a, b, c_, f):     # adjacent arcs of sizes a,b,c with single-interval function f (in units c/3)
    return f(a + b) + f(b + c_) - f(b) - f(a + b + c_)
plat = []
for tot in np.linspace(0.6, 2 * math.pi - 0.05, 120):
    a = b = c_ = tot / 3
    plat.append((tot, t_circ(a, b, c_, S_btz), t_circ(a, b, c_, f_conn)))
below = [p for p in plat if p[0] < th_star]; above = [p for p in plat if p[0] > th_star + 1e-6]
check(f"B1 BTZ on a circle: CMI = (c/3) t while A+B+C is smaller than the plateau size {th_star:.3f} (of 2 pi)",
      all(abs(p[1] - p[2]) < 1e-12 for p in below))
check("B2 beyond it the whole region's entropy switches to 'complement + horizon' and CMI rises above the clock value: the geometric reading breaks",
      all(p[1] > p[2] + 1e-9 for p in above), f"largest excess {max(p[1]-p[2] for p in above):.3f} (units c/3)")
fp = lambda th: rp / 2 / math.tanh(rp * th / 2)
lump = fp(th_star) + fp(2 * math.pi - th_star)
gl_ok = True
for _ in range(300):
    a, b, c1, c2 = rng.uniform(0.2, 2.0, 4)
    I = lambda X, Y, Z: S_btz(X + Z) + S_btz(Y + Z) - S_btz(Z) - S_btz(X + Y + Z) if X + Y + Z < 2 * math.pi else None
    if a + b + c1 + c2 >= 2 * math.pi - 1e-9: continue
    I12, I1, I2 = I(a, c1 + c2, b), I(a, c1, b), I(a, c2, b + c1)
    x = lambda i: 1 - math.exp(-i)
    gl_ok &= min(I12, I1, I2) >= -1e-12 and abs(x(I12) - circ(x(I1), x(I2), -1)) < 1e-12
check("B3 what replaces it: the link-weight built from the true entropies gains a positive lump on geodesics of length theta*, plus the complement's weight beyond; the entropic clock still glues by o_(-1)",
      lump > 0 and gl_ok, f"lump weight {lump:.3f} (units c/3)")
R["btz"] = {"plat": plat, "th_star": th_star, "rp": rp}

# ---------------------------------------------------------------- 4. heavy particle (conical defect)
print("== 4. heavy primary state (conical defect, alpha = sqrt(1 - 24h/c)) on a circle")
cone = {}
for alpha in [0.5, 0.8]:
    f_one = lambda th, al=alpha: math.log(math.sin(al * th / 2))
    S_cone = lambda th, al=alpha: min(math.log(math.sin(al * th / 2)), math.log(math.sin(al * (2 * math.pi - th) / 2)))
    rows = [(tot, t_circ(tot/3, tot/3, tot/3, S_cone), t_circ(tot/3, tot/3, tot/3, f_one)) for tot in np.linspace(.6, 2*math.pi-.05, 120)]
    cone[alpha] = rows
    check(f"C1 alpha={alpha}: up to half the circle, CMI = (c/3) t_alpha with t from the cross-ratio of the points e^(i alpha u)",
          all(abs(r[1] - r[2]) < 1e-12 for r in rows if r[0] <= math.pi))
    check(f"C2 alpha={alpha}: beyond half the circle the pure-state flip S(A) = S(complement) takes over and CMI departs from the clock",
          any(abs(r[1] - r[2]) > 1e-3 for r in rows if r[0] > math.pi + .1))
    ok = True
    for _ in range(200):
        u = np.sort(rng.uniform(0, math.pi, 5)); z = np.exp(1j * alpha * u)
        cr = lambda p, q, r_, s: abs((z[q]-z[p])*(z[s]-z[r_])/((z[r_]-z[p])*(z[s]-z[q])))
        x1, x2, x12 = cr(0, 1, 2, 3), cr(0, 1, 3, 4), cr(0, 1, 2, 4)
        ok &= abs(x12 - circ(x1, x2, -1)) < 1e-9
    check(f"C3 alpha={alpha}: gluing by o_(-1) still holds for the defect's cross-ratios (projective identity)", ok)
R["cone"] = {str(a): v for a, v in cone.items()}

# ---------------------------------------------------------------- 5. local quench (free fermions on a ring)
print("== 5. local excitation: add one fermion wave packet to the ring's ground state and let it split and fly apart")
N = 602; n = np.arange(N)
ks = 2 * math.pi * np.arange(-(N // 2), N // 2 + 1)[: N] / N
ks = 2 * math.pi * (np.arange(N) - N // 2) / N
Uk = np.exp(1j * np.outer(n, ks)) / math.sqrt(N)
occ = (np.abs(ks) < math.pi / 2).astype(float)
C0 = (Uk * occ) @ Uk.conj().T
x0, sig, k0 = N // 2, 6.0, 0.6 * math.pi
phi = np.exp(-(n - x0) ** 2 / (2 * sig ** 2)) * np.cos(k0 * (n - x0)); phi = phi / np.linalg.norm(phi)
pt = phi - C0 @ phi; pt = pt / np.linalg.norm(pt)
eps_k = -2 * np.cos(ks); vg = 2 * math.sin(k0)
def C_at(tau):
    psi = Uk @ (np.exp(-1j * eps_k * tau) * (Uk.conj().T @ pt))
    return C0 + np.outer(psi, psi.conj()), psi
def S_ring(C, sites): s = np.asarray(sites) % N; return S_from_C(C[np.ix_(s, s)])
def cmi_ring(C, a, b, c_, start):
    A = list(range(start, start + a)); B = list(range(start + a, start + a + b)); Cc = list(range(start + a + b, start + a + b + c_))
    return S_ring(C, A + B) + S_ring(C, B + Cc) - S_ring(C, B) - S_ring(C, A + B + Cc)
t_ring = lambda a, b, c_: math.log(math.sin(math.pi*(a+b)/N)*math.sin(math.pi*(b+c_)/N)/(math.sin(math.pi*b/N)*math.sin(math.pi*(a+b+c_)/N)))
la, lb, lc = 50, 100, 50; start = x0 - lb // 2 - la      # A left of centre, B straddles the centre, C right
vac = cmi_ring(C0, la, lb, lc, start)
check("Q0 ground state of the ring: CMI = t/3 (chord form)", abs(vac - t_ring(la, lb, lc) / 3) < 0.01, f"{vac:.4f} vs {t_ring(la,lb,lc)/3:.4f}")
quench = []
for tau in np.linspace(0, 80, 41):
    C, psi = C_at(tau); quench.append((tau, cmi_ring(C, la, lb, lc, start) - vac))
inAC = [q for q in quench if lb / 2 + 3 * sig < vg * q[0] < lb / 2 + la - 3 * sig]
early = [q for q in quench if vg * q[0] < lb / 2 - 3 * sig]
late = [q for q in quench if vg * q[0] > lb / 2 + la + 3 * sig]
check("Q1 while the two halves of the particle sit in B (not yet out): no change, the clock reading holds",
      all(abs(q[1]) < 0.02 for q in early), f"max {max(abs(q[1]) for q in early):.3f}")
check("Q2 once one half is in A and the other in C: CMI jumps by 2 ln 2 = 1.386 (one shared bit counted from both sides): the reading breaks",
      all(abs(q[1] - 2 * math.log(2)) < 0.06 for q in inAC), f"values {[round(q[1],3) for q in inAC][:4]}")
check("Q3 after both halves leave A and C: back to the clock value", all(abs(q[1]) < 0.05 for q in late), f"max {max(abs(q[1]) for q in late):.3f}" if late else "none")
R["quench"] = {"rows": quench, "vg": vg, "la": la, "lb": lb, "lc": lc, "sig": sig}
C, psi = C_at(25.0); R["quench"]["density"] = (np.abs(psi) ** 2).tolist()
# the atom: extra CMI only when (left half, right half) lies in the A x C rectangle
Cq, _ = C_at(25.0); dpos = vg * 25.0
cases = {"A x C": (x0 - int(dpos) - 20, 40, 2 * int(dpos) - 40, 40),
         "A x B": (x0 - int(dpos) - 20, 40, 2 * int(dpos) + 10, 40),
         "both outside": (x0 + int(dpos) + 40, 40, 40, 40)}
extra = {}
for k_, (st, a, b, c_) in cases.items():
    extra[k_] = cmi_ring(Cq, a, b, c_, st) - cmi_ring(C0, a, b, c_, st)
check("Q4 the extra weight is a point lump: present only when the pair (left half, right half) falls in the A x C rectangle",
      abs(extra["A x C"] - 2 * math.log(2)) < 0.08 and abs(extra["A x B"]) < 0.08 and abs(extra["both outside"]) < 0.03,
      " ".join(f"{k}: {v:.3f}" for k, v in extra.items()))
R["quench"]["cases"] = extra

# ---------------------------------------------------------------- 6. what never breaks: chain rule and strong subadditivity
print("== 6. the general law: x := 1 - exp(-3 I(A:C|B)/c) glues by o_(-1) in ANY state (chain rule) and is a probability (strong subadditivity)")
glue = []
def glue_test(Sf, a, b, c1, c2, start=0, name=""):
    A = list(range(start, start + a)); B = list(range(start + a, start + a + b))
    C1 = list(range(start + a + b, start + a + b + c1)); C2 = list(range(start + a + b + c1, start + a + b + c1 + c2))
    I = lambda X, Y, Z: Sf(X + Z) + Sf(Y + Z) - Sf(Z) - Sf(X + Y + Z)
    I12 = I(A, C1 + C2, B); I1 = I(A, C1, B); I2 = I(A, C2, B + C1)
    x = lambda i: 1 - math.exp(-3 * i)
    glue.append((name, x(I12), circ(x(I1), x(I2), -1)))
c_heat = toeplitz_C(lambda k: fermi(-2 * 8 * np.cos(k)), 260)
for _ in range(6):
    a, b, c1, c2 = rng.integers(10, 60, 4); glue_test(lambda s: block_S(c_heat, s), a, b, c1, c2, name="heated")
Cq, _ = C_at(25.0)
for _ in range(6):
    a, b, c1, c2 = rng.integers(10, 60, 4); st = rng.integers(0, N)
    glue_test(lambda s: S_ring(Cq, s), a, b, c1, c2, st, name="quench")
okg = all(abs(g[1] - g[2]) < 1e-9 for g in glue)
# random 4-qubit mixed states
ssa_ok = True
for _ in range(300):
    psi = haar_state(rng, 2 ** 6); rho = ptrace(rho_of(psi), [2] * 6, [0, 1, 2, 3])
    Sx = lambda keep: vn(ptrace(rho, [2] * 4, keep)) if len(keep) < 4 else vn(rho)
    I12 = Sx([0, 1]) + Sx([1, 2, 3]) - Sx([1]) - Sx([0, 1, 2, 3])
    I1 = Sx([0, 1]) + Sx([1, 2]) - Sx([1]) - Sx([0, 1, 2]); I2 = Sx([0, 1, 2]) + Sx([1, 2, 3]) - Sx([1, 2]) - Sx([0, 1, 2, 3])
    x = lambda i: 1 - math.exp(-i)
    ssa_ok &= I1 >= -1e-10 and I2 >= -1e-10 and abs(x(I12) - circ(x(I1), x(I2), -1)) < 1e-10
    glue.append(("random qubits", x(I12), circ(x(I1), x(I2), -1)))
check("G1 gluing x(A:C1C2|B) = x(A:C1|B) o_(-1) x(A:C2|BC1) holds exactly in heated, quenched and random 4-qubit states (it IS the chain rule)", okg and ssa_ok)
check("G2 strong subadditivity makes every such x a genuine probability in [0, 1) (300 random mixed states)", ssa_ok)
R["glue"] = glue

# ---------------------------------------------------------------- 7. compression (pressure at low heat) vs heat
print("== 7. compression: heavy matter squeezed toward the black-hole threshold (alpha -> 0), and dense fermions at zero temperature")
def window_pure(alpha):     # largest A+B+C (equal thirds) for which CMI = clock value, pure heavy state
    f1 = lambda th: math.log(math.sin(alpha * th / 2)); Sc = lambda th: min(f1(th), f1(2 * math.pi - th))
    ok = [tot for tot in np.linspace(.3, 2 * math.pi - .01, 2000) if abs(t_circ(tot/3, tot/3, tot/3, Sc) - t_circ(tot/3, tot/3, tot/3, f1)) < 1e-9]
    return max(ok)
def window_btz(r):
    fc = lambda th: math.log(math.sinh(r * th / 2)); Sb = lambda th: min(fc(th), fc(2 * math.pi - th) + math.pi * r)
    ok = [tot for tot in np.linspace(.3, 2 * math.pi - .01, 2000) if abs(t_circ(tot/3, tot/3, tot/3, Sb) - t_circ(tot/3, tot/3, tot/3, fc)) < 1e-9]
    return max(ok)
alphas = [0.9, 0.6, 0.3, 0.1, 0.02]; rps = [0.002, 0.02, 0.2, 0.5, 1.0, 2.0, 4.0]
wp = [(a, window_pure(a)) for a in alphas]; wb = [(r, window_btz(r)) for r in rps]
for a, w in wp: print(f"   compressed, alpha={a}: clock holds for A+B+C up to {w/math.pi:.3f} pi")
for r, w in wb: print(f"   heated BTZ, r+={r} (T = r+/2pi): clock holds up to {w/math.pi:.3f} pi")
check("K1 squeezing heavy matter all the way to the black-hole threshold never moves the window: the clock holds for exactly half the circle",
      all(abs(w - math.pi) < 0.01 for a, w in wp))
check("K2 heat widens the window toward the whole circle (the plateau); right at the threshold heat and squeezing meet at half",
      abs(wb[0][1] - math.pi) < 0.02 and all(wb[i][1] < wb[i + 1][1] + 1e-9 for i in range(len(wb) - 1)) and wb[-1][1] > 1.85 * math.pi)
ok = True
for a in alphas:
    for tot in [0.5, 1.5, 2.8]:
        ta = t_circ(tot/3, tot/3, tot/3, lambda th: math.log(math.sin(a * th / 2)))
        ok &= ta > 0
lim = t_circ(1, 1, 1, lambda th: math.log(math.sin(0.001 * th / 2))); flat = math.log(2 * 2 / (1 * 3))
check("K3 as alpha -> 0 the squeezed clock becomes the flat-space (line) clock: t -> log[(a+b)(b+c)/(b(a+b+c))]", ok and abs(lim - flat) < 1e-5)
R["window"] = {"pure": wp, "btz": wb}
dense = []
for fill in [0.5, 0.8, 0.9, 0.95, 0.98, 0.99]:
    kF = fill * math.pi
    c = toeplitz_C(lambda k, kF=kF: (np.abs(k) < kF).astype(float), 170); Sf = lambda s_: block_S(c, s_)
    I = cmi(Sf, 60, 40, 60); dense.append((fill, I, t_line(60, 40, 60) / 3, 1 / (1 - fill)))
for d in dense: print(f"   filling {d[0]:.2f} (hole spacing {d[3]:.0f} sites): CMI={d[1]:.4f}  clock t/3={d[2]:.4f}")
check("K4 dense fermions at zero temperature: the clock holds while blocks are much bigger than the spacing between holes",
      all(abs(d[1] - d[2]) < 0.01 for d in dense if d[3] <= 10))
check("K5 squeezed until the spacing approaches the block size, it breaks (the continuum description itself fails)",
      abs(dense[-1][1] - dense[-1][2]) > 0.02 and all(abs(dense[i][1]-dense[i][2]) <= abs(dense[i+1][1]-dense[i+1][2]) + 1e-4 for i in range(len(dense)-1)), f"miss at 99% filling {abs(dense[-1][1]-dense[-1][2]):.3f}")
R["dense"] = dense

# ---------------------------------------------------------------- 8. tunnelling (optional): the Gamow factor as a 'no-tick' probability
print("== 8. tunnelling through a barrier: transmission T plays the role of 1 - x (no tick), so blocking glues by o_(-1) in the WKB limit")
def transmit(E, layers):    # exact 1D transfer matrix, hbar = 2m = 1; layers: list of (V, width)
    regions = [(0.0, None)] + layers + [(0.0, None)]
    ks = [np.sqrt(complex(E - V)) for V, _ in regions]
    M = np.eye(2, dtype=complex); x = 0.0
    for j in range(len(regions) - 1):
        k1, k2 = ks[j], ks[j + 1]
        D = lambda k, x: np.array([[np.exp(1j * k * x), np.exp(-1j * k * x)], [1j * k * np.exp(1j * k * x), -1j * k * np.exp(-1j * k * x)]])
        M = np.linalg.solve(D(k2, x), D(k1, x)) @ M
        if regions[j + 1][1] is not None: x += regions[j + 1][1]
    return float(abs(1 / M[1, 1]) ** 2) if abs(M[1, 1]) > 0 else 0.0
tun = []
for L in [1, 2, 3, 4, 6]:
    T1 = transmit(1.0, [(5.0, L)]); T2 = transmit(1.0, [(5.0, L)]); T12 = transmit(1.0, [(5.0, 2 * L)])
    tun.append((L, -math.log(T12), -math.log(T1) - math.log(T2)))
for r in tun: print(f"   slice width {r[0]}: -ln T(both) = {r[1]:.3f}   -ln T1 - ln T2 = {r[2]:.3f}")
rel = [abs(r[1] - r[2]) / r[1] for r in tun]
check("T1 the exponent of tunnelling adds over adjacent slices, so 'blocked' = blocked1 o_(-1) blocked2 up to edge factors that fade for thick barriers",
      rel[-1] < 0.1 and rel[-1] < rel[0], f"relative miss {rel[0]:.2f} -> {rel[-1]:.3f}")
R["tunnel"] = tun

# ---------------------------------------------------------------- 9. the environment knobs: does k move?
print("== 9. the hypercube in a bath with knobs (heat beta, pressure/density mu, spin Omega): does the gluing dial k move?")
import sympy as sp
def k_eff(x1, x2, x12): return (x12 - x1 - x2) / (x1 * x2)
ks_ = []
def triples(Sf, cfgs, scale=3.0, name=""):
    for a, b, c1, c2 in cfgs:
        A = list(range(0, a)); B = list(range(a, a + b)); C1 = list(range(a + b, a + b + c1)); C2 = list(range(a + b + c1, a + b + c1 + c2))
        I = lambda X, Y, Z: Sf(X + Z) + Sf(Y + Z) - Sf(Z) - Sf(X + Y + Z)
        x = lambda i: 1 - math.exp(-scale * i)
        x1, x2, x12 = x(I(A, C1, B)), x(I(A, C2, B + C1)), x(I(A, C1 + C2, B))
        if x1 > 1e-4 and x2 > 1e-4: ks_.append((name, k_eff(x1, x2, x12)))
cfgs = [(30, 10, 20, 25), (40, 15, 15, 30), (25, 20, 30, 20), (50, 5, 10, 10)]
envs = {"vacuum": lambda k: (np.abs(k) < math.pi / 2).astype(float),
        "hot (beta=15)": lambda k: fermi(-2 * 15 * np.cos(k)),
        "spinning (15 / 60)": lambda k: np.where(k > 0, fermi(-2 * 60 * np.cos(k)), fermi(-2 * 15 * np.cos(k))),
        "dense (90% full)": lambda k: (np.abs(k) < 0.9 * math.pi).astype(float),
        "dense + warm (mu, beta=40)": lambda k: fermi(40 * (-2 * np.cos(k) + 2 * math.cos(0.7 * math.pi)))}
for nm, f in envs.items():
    c = toeplitz_C(f, 200); triples(lambda s_: block_S(c, s_), cfgs, name=nm)
Cq, _ = C_at(25.0)
triples(lambda s_: S_ring(Cq, [q + x0 - 60 for q in s_]), cfgs, name="quench")
for al in [0.5, 0.05]:
    f1 = lambda th, al=al: math.log(math.sin(al * th / 2))
    for a, b, c1, c2 in [(0.4, 0.3, 0.5, 0.6), (0.2, 0.5, 0.3, 0.3)]:
        I = lambda X, Y, Z: f1(X + Z) + f1(Y + Z) - f1(Z) - f1(X + Y + Z)
        x = lambda i: 1 - math.exp(-i)
        ks_.append((f"squeezed alpha={al}", k_eff(x(I(a, c1, b)), x(I(a, c2, b + c1)), x(I(a, c1 + c2, b)))))
for a, b, c1, c2 in [(1.0, 0.5, 1.5, 2.5), (2.0, 0.3, 1.0, 2.5)]:
    I = lambda X, Y, Z: S_btz(X + Z) + S_btz(Y + Z) - S_btz(Z) - S_btz(X + Y + Z)
    x = lambda i: 1 - math.exp(-i)
    ks_.append(("black hole past plateau", k_eff(x(I(a, c1, b)), x(I(a, c2, b + c1)), x(I(a, c1 + c2, b)))))
for nm in dict.fromkeys(n_ for n_, _ in ks_):
    vals = [k for n_, k in ks_ if n_ == nm]; print(f"   {nm:28s} k_eff = {np.mean(vals):+.6f}  (spread {np.ptp(vals):.1e}, {len(vals)} triples)")
check("D1 k does not move: every environment (hot, spinning, dense, warm+dense, squeezed, quenched, past the plateau) glues at k = -1 exactly",
      all(abs(k + 1) < 1e-6 for _, k in ks_), f"max |k+1| {max(abs(k+1) for _,k in ks_):.1e}")
R["keff"] = ks_
x_, y_, lam, sc = sp.symbols('x y lambda s', positive=True)
f = lambda z: 1 - (1 - z) ** lam
auto_ok = True
for _ in range(2000):
    xa, ya, la_ = rng.uniform(0, .999), rng.uniform(0, .999), rng.uniform(.05, 20)
    g = lambda z: 1 - (1 - z) ** la_
    auto_ok &= abs(g(circ(xa, ya, -1)) - circ(g(xa), g(ya), -1)) < 1e-12
check("D2 turning a knob changes the clock's RATE: x -> 1 - (1-x)^lambda, and every such change commutes with the union rule (automorphism)", auto_ok)
check("D3 the only way to move k is to change the ceiling: y = s*x glues by o_k with k = -1/s (x stays a probability, so s = 1)",
      sp.simplify(sc * (x_ + y_ - x_ * y_) - (sc * x_ + sc * y_ + (-1 / sc) * (sc * x_) * (sc * y_))) == 0)
run = []
for size in np.linspace(0.05, 4, 40):        # equal thirds, size in units of beta (line)
    a = size / 3
    tb = math.log(math.sinh(math.pi*2*a)**2 / (math.sinh(math.pi*a) * math.sinh(math.pi*3*a)))
    tv = math.log((2*a)**2 / (a * 3*a)); run.append((size, tb / tv))
check("D4 the rate runs with scale: thermal clock / vacuum clock -> 1 for regions much smaller than beta and -> 0 for much larger (the clock runs, k doesn't)",
      run[0][1] > 0.99 and run[-1][1] < 0.01 and all(run[i][1] >= run[i + 1][1] for i in range(len(run) - 1)))
R["run"] = run
ok = True
for _ in range(300):
    beta = rng.uniform(.5, 5); w = np.sort(rng.uniform(-3, 3, 4)); a1, a2, b1, b2 = w
    xs = math.sinh(math.pi*(a2-a1)/beta)*math.sinh(math.pi*(b2-b1)/beta)/(math.sinh(math.pi*(b1-a1)/beta)*math.sinh(math.pi*(b2-a2)/beta))
    E = np.expm1(2 * math.pi * w / beta); xp = (E[1]-E[0])*(E[3]-E[2])/((E[2]-E[0])*(E[3]-E[1]))
    ok &= abs(xs - xp) < 1e-6 * max(1, xs)
check("D5 heat acts on positions by our height flow: the hot cross-ratio is the cold cross-ratio of E(2 pi w / beta) = e^(2 pi w/beta) - 1", ok)

# ---------------------------------------------------------------- 10. the bath as a bigger cube: tracing it out
print("== 10. the bath as a bigger cube: system corners inside a system+bath cube; forget (sum over) the bath and see what dial is left")
def fold(xs, k):
    out = 0.0
    for v in xs: out = circ(out, v, k)
    return out
def corner_sum(xs, k):
    import itertools as it
    n = len(xs); tot = 0.0
    for r in range(1, n + 1):
        for T in it.combinations(range(n), r): tot += k ** (r - 1) * np.prod([xs[i] for i in T])
    return tot
def bath_union(a, levels):          # levels: list of (prob, intensity lambda); events fire independently with prob lambda*a_i GIVEN the bath level
    return sum(pb * (1 - np.prod([1 - lam * ai for ai in a])) for pb, lam in levels)
ok = True; rows = []
for _ in range(200):
    pi_on = rng.uniform(.05, 1); Lam = rng.uniform(.2, 1); n_ = rng.integers(2, 6); a = rng.uniform(0, 1, n_)
    levels = [(pi_on, Lam), (1 - pi_on, 0.0)]
    x = [pi_on * Lam * ai for ai in a]                 # each event's marginal probability
    U = bath_union(a, levels)
    ok &= abs(U - fold(x, -1 / pi_on)) < 1e-12 and abs(U - corner_sum(x, -1 / pi_on)) < 1e-12
check("E1 two-state bath (on with probability pi, off otherwise): forgetting it turns o_(-1) into o_k' with k' = -1/pi, exactly, for 2 to 5 events (corner sum too)", ok)
mc_ok = True
pi_on, Lam, a = 0.3, 0.8, np.array([.5, .4, .7])
on = rng.uniform(size=400000) < pi_on
fires = (rng.uniform(size=(400000, 3)) < Lam * a) & on[:, None]
mc = fires.any(1).mean(); pred = fold(list(pi_on * Lam * a), -1 / pi_on)
check("E2 Monte Carlo of the same bath agrees", abs(mc - pred) < 0.003, f"sim {mc:.4f} vs o_(-1/pi) {pred:.4f}; plain union would give {fold(list(pi_on*Lam*a), -1):.4f}")
def kprime(beta, Delta, mu=0.0, Omega=0.0, J=0.0):
    return -(1 + np.exp(beta * (Delta - mu - Omega * J)))
kn = []
for beta in [0.0, 0.5, 1, 2, 4]:
    for mu in [0.0, 1.0, 2.0]:
        pon = 1 / (1 + math.exp(beta * (1.0 - mu)))
        kn.append((beta, mu, -1 / pon, kprime(beta, 1.0, mu)))
check("E3 with the bath's 'on' level at energy Delta in a bath of temperature 1/beta, chemical potential mu, spin Omega: k' = -(1 + e^{beta(Delta - mu - Omega J)}) = -1/(Fermi occupation)",
      all(abs(r[2] - r[3]) < 1e-9 for r in kn))
check("E4 the knobs move k': infinitely hot -> k' = -2; cold with Delta > mu -> k' -> -infinity (bath almost never on); cold with Delta < mu -> k' -> -1 (bath always on: plain union rule)",
      abs(kprime(0, 1) + 2) < 1e-12 and kprime(50, 1, 0) < -1e10 and abs(kprime(50, 1, 2) + 1) < 1e-12)
R["kprime"] = kn
# three-level bath: no single k' works
pairs_k, trip_k = [], []
for _ in range(100):
    probs = rng.dirichlet(np.ones(3)); lams = np.sort(rng.uniform(0.05, 1, 3)); a = rng.uniform(.1, 1, 3)
    levels = list(zip(probs, lams)); x = [sum(p * l for p, l in levels) * ai for ai in a]
    U2 = bath_union(a[:2], levels); k2 = k_eff(x[0], x[1], U2)
    U3 = bath_union(a, levels)
    # solve fold(x, k) = U3 for k by bisection
    lo, hi = -60.0, 0.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if fold(x, mid) > U3: hi = mid
        else: lo = mid
    pairs_k.append(k2); trip_k.append((lo + hi) / 2)
gap = np.max(np.abs(np.array(pairs_k) - np.array(trip_k)))
check("E5 a bath with three or more levels leaves no single dial: the k' that fits two events misfits three (the reduced rule leaves the o_k family)",
      gap > 0.05, f"largest pair-vs-triple mismatch {gap:.2f}")
ok = True
for _ in range(100):
    probs = rng.dirichlet(np.ones(3)); lams = rng.uniform(0, 1, 3); a = rng.uniform(0, 1, 4)
    for p, l in zip(probs, lams):
        ok &= abs(bath_union(a, [(1.0, l)]) - fold([l * ai for ai in a], -1)) < 1e-12
check("E6 remembering the bath (conditioning on its state) restores k = -1 in every bath state: k moves only when the bath is forgotten", ok)

json.dump(R, open("results2.json", "w"))
print(f"\n{passed} passed, {failed} failed")
