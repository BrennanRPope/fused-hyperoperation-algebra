"""Lee-Yang refactor across roundness: the far-zone test at every p (DRAFT, numerical exploration, not a certificate).

The p = 3 certificate's far zone rests on an exact contour formula. It holds at every roundness p > 1, not only whole p:
    G_c(s) = C_c(s) + |J_c(s)| sin(s - pi c/2 + arg J_c(s)),
    C_c(s) = -Im int_0^inf (1 - e^(i pi p/2) y^p)^c e^(-s y) dy                (from turning the path at 0 upward)
    J_c(s) = int_0^inf y^c w(y)^c e^(-s y) dy,  w(y) = (1 - (1 + i y)^p)/(-i y)  (from turning the path at 1 upward)
with the continuous branch of log w. This script
  R1 checks the formula against the defining integral at non-whole roundness;
  R2 evaluates the three far-zone inequalities with the exact C and J (no remainder bounds) on a grid of
     roundness p, steps b in (0, p - 1] and s in [2, 300], and records where they take over;
  R3 checks that far out the tight inequality's margin is the kink margin 1 - |sin(pi p/2)|/p^(p-1);
  R4 counts the zeros of G at b below the takeover point (the near zone a certificate still has to box).
"""
import mpmath as mp, math, time, json
mp.mp.dps = 15
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
t0 = time.perf_counter()
def logw(y, p):
    q = mp.exp(-p*mp.log(1 + 1j*y))
    return -1j*mp.pi + p*mp.log(1 + 1j*y) + mp.log(1 - q) - (mp.log(y) - 1j*mp.pi/2)
def J(c, s, p, d=0):
    return mp.quad(lambda y: (-y)**d*y**c*mp.exp(c*logw(y, p))*mp.exp(-s*y), [0, 1/s, 5/s, 20/s, 60/s])
def C(c, s, p, d=0):
    th = mp.pi*p/2
    f = lambda y: (-y)**d*mp.im((1 - mp.exp(1j*th)*y**p)**c)*mp.exp(-s*y)
    try:
        return -mp.quad(f, [0, 1/s, 5/s, 20/s, 60/s])
    except ZeroDivisionError:   # mpmath's error estimate divides by zero when two refinements agree exactly
        return -mp.quad(f, [0, 1/s, 5/s, 20/s, 60/s], method="gauss-legendre")
def G(c, s, p): return mp.quad(lambda x: (1 - x**p)**c*mp.cos(s*x), [0, 0.5, 0.9, 1])
# R1
gaps = []
for p, c, s in [(1.5, 0.4, 5.0), (2.5, 1.3, 7.0), (3.0, 2.0, 21.0), (3.7, 2.2, 11.0), (5.3, 3.1, 15.0)]:
    j = J(c, s, p); gaps.append(float(abs(G(c, s, p) - C(c, s, p) - abs(j)*mp.sin(s - mp.pi*c/2 + mp.arg(j)))))
OUT["R1_gaps"] = gaps
check("R1 the contour formula holds at non-whole roundness (p = 1.5, 2.5, 3, 3.7, 5.3)", max(gaps) < 1e-8, f"largest gap {max(gaps):.0e}")
# R2, R3
def criteria(p, b, s):
    a = b + 1
    Jb, Jb1, Ja = J(b, s, p), J(b, s, p, 1), J(a, s, p)
    Cb, Cb1, Ca = C(b, s, p), C(b, s, p, 1), C(a, s, p)
    eps = Cb/abs(Jb); r = Jb1/Jb
    epsp = Cb1/abs(Jb) - Cb*mp.re(r)/abs(Jb)
    A = 1 - abs(eps)
    if A <= 0: return float(A), -1.0, -1.0
    B = (1 + mp.im(r))*mp.sqrt(1 - eps**2) - abs(epsp)
    Cc = mp.cos(mp.asin(abs(eps)) + abs(mp.arg(Ja/Jb))) - abs(Ca)/abs(Ja)
    return float(A), float(B), float(Cc)
S = [2*1.09**k for k in range(59)]   # 2 .. ~320
PS = [1.25, 1.5, 2.5, 3.0, 3.5, 5.0, 7.0]; FR = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
def scan_one(p):
    T = 0; worstC = 9; rows = []
    for f in FR:
        b = f*(p - 1); last = 0
        for s in S:
            A, B, Cc = criteria(p, b, s)
            if min(A, B, Cc) <= 0: last = s
            if s > 200 and f == 1.0: worstC = min(worstC, Cc)
        rows.append(dict(b=b, takeover=last)); T = max(T, last)
    kink = 1 - abs(math.sin(math.pi*p/2))/p**(p - 1)
    return p, dict(takeover=T, rows=rows, far_margin=worstC, kink_margin=kink)
from multiprocessing import Pool
with Pool(len(PS)) as pool:
    scan = dict(pool.map(scan_one, PS))
for p, v in scan.items():
    print(f"   p = {p}: takeover s ~ {v['takeover']:.1f}, far margin {v['far_margin']:.3f} vs kink margin {v['kink_margin']:.3f}", flush=True)
OUT["scan"] = scan
check("R2 with the exact contour pieces, all three far-zone inequalities hold at every sampled s past a takeover point, "
      "for 7 roundnesses and 6 steps each", all(v["takeover"] < 300 for v in scan.values()),
      "takeover s: " + ", ".join(f"p={p}: {v['takeover']:.1f}" for p, v in scan.items()))
check("R3 far out (s > 200, top step b = p - 1) the tight margin is the kink margin 1 - |sin(pi p/2)|/p^(p-1), for p <= 5",
      all(abs(v["far_margin"] - v["kink_margin"]) < 0.02 for p, v in scan.items() if p <= 5),
      ", ".join(f"p={p}: {v['far_margin']:.3f} vs {v['kink_margin']:.3f}" for p, v in scan.items()))
# R4: zeros of G at b below the takeover
def zeros_below(p, b, T):
    zs = []; prev = G(b, 0.05, p); s = 0.05; h = 0.05
    while s < T:
        s2 = s + h; v = G(b, s2, p)
        if prev*v < 0: zs.append(float(mp.findroot(lambda x: G(b, x, p), (s, s2), solver="illinois", verify=False)))
        prev, s = v, s2
    return zs
def nz_one(p):
    T = scan[p]["takeover"]
    return p, max(len(zeros_below(p, r["b"], T)) for r in scan[p]["rows"])
with Pool(len(PS)) as pool:
    nz = dict(pool.map(nz_one, PS))
OUT["near_zone_zeros"] = nz
check("R4 near zone: zeros of G at b below the takeover point, the part a certificate still boxes (most over b)",
      True, ", ".join(f"p={p}: {n}" for p, n in nz.items()))
OUT["seconds"] = time.perf_counter() - t0
print(f"{sum(RES)}/{len(RES)} pass, {OUT['seconds']:.0f} s")
json.dump(OUT, open("general_p_far.json", "w"), indent=1, default=float)
