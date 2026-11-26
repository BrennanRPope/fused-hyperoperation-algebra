"""Lee-Yang near zone for every roundness 1.1 <= p <= 5 (DRAFT; rigorous interval computation in Python, not in Lean).

The far zone (far_bounds_p.py) proves the climb b -> b + 1 for s >= T*(p). This script boxes everything below:
for every p in [1.1, 5] and every step b in [0, p - 1], every zero of G_b in [0, TN] is isolated, and G at a = b + 1
has a definite sign there that alternates from zero to zero. TN sits past T*(p), so near and far zones share a zero.

How a whole range of p is covered at once
  * Exact values at nodes. At rational (p, c), G_c(s) = sum_n Re(i^n) s^n m_n / n! with m_n = B((n+1)/p, c+1)/p.
    Moments are interval Gamma values; the sum is an exact big-integer Horner sum, so each value and s-derivative
    comes with an error under 2^-200 (e^s + 1).
  * A surrogate in (p, c). On each big box (a strip of p times a unit range of c) the values at Chebyshev-like nodes
    are interpolated in both variables. The interpolation error is bounded with no fitting:
      in c:  |d^n G/dc^n| <= int_0^1 |log(1 - x^p)|^n <= int_0^1 |log(1 - x)|^n = n!            (1 - x <= 1 - x^p)
      in p:  (1 - x^p)^c is analytic for Re p > 0 with |(1 - x^p)^c| <= 2^c, so Cauchy gives
             |d^n G/dp^n| / n! <= 2^c / rho^n for any rho < (left end of the strip).
    So |G - surrogate| <= max|omega_c| + Lambda_c max|omega_p| 2^c / rho^n_p, about 1e-9 or less everywhere.
  * Small boxes. Each small box in (p, b) gets the surrogate's exact Taylor expansion at its centre, so its range over
    the box is a tight enclosure (no cancellation), plus the error terms above, plus float rounding.
  * The zero test on each small box is the p = 3 certificate's: cells of width 1/8 in s, Taylor models in s with
    7 derivatives and |d^7 G/ds^7| <= 1/8; a run of cells holds exactly one zero if the ends have opposite signs and
    the slope keeps one sign; G at b + 1 must keep one sign on the run (or, failing that, on every sub-cell of width
    1/128 where G_b may vanish), and the signs must alternate. A box that fails is split in four, up to 5 times;
    sub-boxes wholly past the line b = p - 1 are dropped.

  Z1 every zero isolated on every small box;  Z2 alternation;  Z3 overlap with the far zone;
  Z4 surrogate against direct quadrature at random points (polish);  Z5 the boxes containing p = 3 are certified,
  as the separate p = 3 certificate says they should be.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, math, time, sys, random
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
from multiprocessing import Pool
iv = mp.iv
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
t0 = time.perf_counter()
# Two ranges share this code. NZ_RANGE=1.1-5 (default) is the run in section "N3, near zone"; NZ_RANGE=5-8 extends it.
# Past p = 5 the values near the grid's end fall to about 1e-9, so that run uses more nodes, more digits and more
# s-derivatives (|d^DER G/ds^DER| <= 1/(DER + 1) either way).
CONFIGS = {
    "1.1-5": dict(P=200, DER=7, NC=16, tag="near_zone_p", lo=Fr(11, 10), hi=Fr(5),
                  STRIPS=[(Fr(11, 10), Fr(2), 16, 42), (Fr(2), Fr(3), 12, 31), (Fr(3), Fr(4), 10, 50), (Fr(4), Fr(5), 9, 76)]),
    "5-8": dict(P=420, DER=10, NC=24, tag="near_zone_p58", lo=Fr(5), hi=Fr(8),
                STRIPS=[(Fr(5), Fr(6), 14, 112), (Fr(6), Fr(7), 14, 152), (Fr(7), Fr(8), 14, 205)]),
}
CFG_NAME = os.environ.get("NZ_RANGE", "1.1-5"); CFG = CONFIGS[CFG_NAME]
H = Fr(1, 8); DER = CFG["DER"]; P = CFG["P"]
U = 2.0**-53

# ---------------- exact values at a node ----------------
def node_eval(args):
    """Enclosures of d^j G/ds^j, j < DER, at s = kH, k = 0..TN/H, for rational p, c: arrays mid, rad (k, j)."""
    p, c, TN = args
    T = float(TN); NT = max(40, int(T) + 10)
    while True:
        tl = mp.mpf(T)**(NT + 1)/mp.factorial(NT + 1)/(1 - T/(NT + 2))
        if NT + 2 > T and tl < mp.mpf(2)**(-P): break
        NT += 10
    iv.prec = P + 80; mp.mp.prec = P + 80
    pp = iv.mpf(p.numerator)/p.denominator; cc = iv.mpf(c.numerator)/c.denominator
    gc = iv.gamma(cc + 1); ks = []
    for n in range(NT + DER + 1):
        al = iv.mpf(n + 1)/pp
        x = iv.gamma(al)*gc/(pp*iv.gamma(al + cc + 1))*mp.mpf(2)**P
        k = int(mp.nint(x.mid)); assert (x - k).a > -1 and (x - k).b < 1   # |m_n - k 2^-P| < 2^-P
        ks.append(k)
    W = [1]*(NT + 1)
    for n in range(NT - 1, -1, -1): W[n] = W[n + 1]*(n + 1)              # NT!/n!
    sh = 3                                                                # H = 2^-3
    coefs = []
    for j in range(DER):
        row = []
        for n in range(NT + 1):
            q = (n + j) % 4; sg = 1 if q == 0 else (-1 if q == 2 else 0)
            row.append((sg*ks[n + j]*W[n]) << (sh*(NT - n)) if sg else 0)
        coefs.append(row)
    nk = int(TN/H); denom = W[0] << (P + sh*NT)
    mid = np.zeros((nk + 1, DER)); rad = np.zeros((nk + 1, DER))
    for k in range(nk + 1):
        e = (math.exp(k/8) + 1)*2.0**(-P)*1.01          # moment rounding (sum s^n/n! = e^s) + series tail
        for j in range(DER):
            acc = 0
            for v in reversed(coefs[j]): acc = acc*k + v
            q = acc/denom                                  # correctly rounded
            mid[k, j] = q; rad[k, j] = abs(q)*1.2e-16 + e
    return mid, rad

# ---------------- interpolation machinery (exact rationals) ----------------
def nodes(n):   # Chebyshev-like nodes on [-1, 1], rounded to multiples of 2^-24
    return [Fr(round(math.cos((2*i + 1)*math.pi/(2*n))*2**24), 2**24) for i in range(n)]
def lagrange_polys(us):   # exact coefficients A[i][m] of l_i(u) = sum_m A[i][m] u^m
    n = len(us); A = []
    for i in range(n):
        poly = [Fr(1)]; den = Fr(1)
        for m in range(n):
            if m == i: continue
            new = [Fr(0)]*(len(poly) + 1)
            for d, a in enumerate(poly): new[d + 1] += a; new[d] -= a*us[m]
            poly = new; den *= us[i] - us[m]
        A.append([a/den for a in poly])
    return A
def taylor_rows(A, u0):   # M[alpha][i] = l_i^(alpha)(u0)/alpha!, exact, then as floats
    n = len(A); M = np.zeros((n, n))
    for i in range(n):
        for al in range(n):
            s = Fr(0); pw = Fr(1)
            for m in range(al, n):
                s += A[i][m]*math.comb(m, al)*pw; pw *= u0
            M[al, i] = float(s)
    return M
def omega_lebesgue(us, pieces=4096):   # rigorous max |prod (u - u_i)| and max sum |l_i(u)| over [-1, 1]
    iv.prec = 60; om = 0.0; lb = 0.0; n = len(us); ui = [iv.mpf(u.numerator)/u.denominator for u in us]
    dens = []
    for i in range(n):
        d = iv.mpf(1)
        for m in range(n):
            if m != i: d *= ui[i] - ui[m]
        dens.append(d)
    for k in range(pieces):
        x = iv.mpf([-1 + 2*k/pieces, -1 + 2*(k + 1)/pieces])
        fac = [x - u for u in ui]
        w = iv.mpf(1)
        for f in fac: w *= f
        om = max(om, float(abs(w).b))
        L = iv.mpf(0)
        for i in range(n):
            pr = iv.mpf(1)
            for m in range(n):
                if m != i: pr *= fac[m]
            L += abs(pr/dens[i])
        lb = max(lb, float(L.b))
    return om, lb

# ---------------- big boxes ----------------
STRIPS = CFG["STRIPS"]
NC = CFG["NC"]
def big_boxes():
    out = []
    for si, (pl, ph, npn, TN) in enumerate(STRIPS):
        cmax = math.ceil(ph)            # b <= p - 1 lies in [0, ceil(ph) - 1], a = b + 1 in [1, ceil(ph)]
        for k in range(cmax):
            out.append(dict(strip=si, pl=pl, ph=ph, cl=Fr(k), ch=Fr(k + 1), np=npn, nc=NC, TN=TN))
    return out

def far_T(pl, ph):
    D = json.load(open("far_bounds_p.json"))
    Ts = [s["T"] for s in D["F2_strips"] if s["p"][1] > float(pl) + 1e-12 and s["p"][0] < float(ph) - 1e-12]
    return max(Ts)

# ---------------- enclosures on a small box ----------------
class Surrogate:
    def __init__(self, bb, V, R, up, uc):
        self.bb = bb; self.V = V; self.R = R; self.absV = np.abs(V)
        self.up, self.uc = up, uc
        self.Ap, self.Ac = lagrange_polys(up), lagrange_polys(uc)
        pl, ph, cl, ch = bb["pl"], bb["ph"], bb["cl"], bb["ch"]
        self.pm, self.hp = (pl + ph)/2, (ph - pl)/2
        self.cm, self.hc = (cl + ch)/2, (ch - cl)/2
        omp, lbp = omega_lebesgue(up); omc, lbc = omega_lebesgue(uc)
        rho = pl - Fr(1, 16)
        eps_c = omc*float(self.hc)**len(uc)                                     # max|omega_c| * n!/n!
        eps_p = lbc*omp*float(self.hp)**len(up)*2.0**float(ch)/float(rho)**len(up)
        self.eps = (eps_c + eps_p)*1.001
        self.cache_p, self.cache_c = {}, {}
    def M(self, which, x0):
        cache, A, m, h = (self.cache_p, self.Ap, self.pm, self.hp) if which == "p" else (self.cache_c, self.Ac, self.cm, self.hc)
        if x0 not in cache: cache[x0] = taylor_rows(A, (x0 - m)/h)
        return cache[x0]
    def enclose(self, p0, hpw, c0, hcw):
        """(lo, hi) arrays (k, j) enclosing d^j G/ds^j over p in p0 +- hpw, c in c0 +- hcw, at every grid s."""
        Mp, Mc = self.M("p", p0), self.M("c", c0)
        n_p, n_c = Mp.shape[0], Mc.shape[0]
        C = np.einsum("ai,ikTj,bk->abTj", Mp, self.V, Mc, optimize=True)
        aMp, aMc = np.abs(Mp), np.abs(Mc)
        S = np.einsum("ai,ikTj,bk->abTj", aMp, self.absV, aMc, optimize=True)
        Rr = np.einsum("ai,ikTj,bk->abTj", aMp, self.R, aMc, optimize=True)
        err = (n_p*n_c + 6)*U*1.01*S + Rr
        h = float(hpw/self.hp); k = float(hcw/self.hc)
        wa = h**np.arange(n_p); wb = k**np.arange(n_c)
        Wt = np.outer(wa, wb)
        dev = np.einsum("ab,abTj->Tj", Wt, np.abs(C)) - np.abs(C[0, 0])
        tot_err = np.einsum("ab,abTj->Tj", Wt, err)
        spread = (dev + tot_err)*(1 + 1e-12) + self.eps + 1e-300
        return C[0, 0] - spread, C[0, 0] + spread

def cell_range(lo, hi, k, j):   # enclosure of d^j G/ds^j over [kH, (k+1)H] from the data at kH
    hh = 1/8
    var = 0.0
    for m in range(1, DER - j):
        var += max(abs(lo[k, j + m]), abs(hi[k, j + m]))*hh**m/math.factorial(m)
    var += hh**(DER - j)/(math.factorial(DER - j)*(DER + 1))
    var *= 1 + 1e-12
    return lo[k, j] - var, hi[k, j] + var
NS = 16                             # sub-cells per cell, used to pin down where G_b's zero sits inside a run
def sub_range(lo, hi, k, t0, t1):   # enclosure of G over [kH + t0, kH + t1], 0 <= t0 < t1 <= H, mean-value form
    tc = (t0 + t1)/2; w = (t1 - t0)/2
    plo = sum(lo[k, m]*tc**m/math.factorial(m) for m in range(DER))
    phi = sum(hi[k, m]*tc**m/math.factorial(m) for m in range(DER))
    dp = sum(max(abs(lo[k, m]), abs(hi[k, m]))*t1**(m - 1)/math.factorial(m - 1) for m in range(1, DER))
    rem = t1**DER/(math.factorial(DER)*(DER + 1))
    var = (dp*w + rem)*(1 + 1e-12) + 1e-300
    return plo - var, phi + var
def zero_sign(Eb, Ea, kl, kr):      # sign of G_a on every sub-cell of the run where G_b might vanish (0 if mixed)
    lob, hib = Eb; loa, hia = Ea
    sg = None; hs = 1/(8*NS)
    for q in range(kl, kr):
        for r in range(NS):
            l, h = sub_range(lob, hib, q, r*hs, (r + 1)*hs)
            if l > 0 or h < 0: continue
            al, ah = sub_range(loa, hia, q, r*hs, (r + 1)*hs)
            s = 1 if al > 0 else (-1 if ah < 0 else 0)
            if s == 0 or (sg is not None and s != sg): return 0
            sg = s
    return sg or 0
def certify(Eb, Ea, nk, Tstar):
    lob, hib = Eb; loa, hia = Ea
    cells = []
    for k in range(nk):
        l, h = cell_range(lob, hib, k, 0); cells.append(l > 0 or h < 0)
    zeros = []; ok = True; i = 0
    while i < nk:
        if cells[i]: i += 1; continue
        i2 = i
        while i2 < nk and not cells[i2]: i2 += 1
        kl, kr = i, i2
        if kr >= nk:                # a run cut off by the grid's end: fine only if it lies wholly in the far zone
            if kl/8 < Tstar: ok = False
            break
        sl = 1 if lob[kl, 0] > 0 else (-1 if hib[kl, 0] < 0 else 0)
        sr = 1 if lob[kr, 0] > 0 else (-1 if hib[kr, 0] < 0 else 0)
        dfn = sl*sr < 0
        mono = all((lambda r: r[0] > 0 or r[1] < 0)(cell_range(lob, hib, q, 1)) for q in range(kl, kr))
        ga = [cell_range(loa, hia, q, 0) for q in range(kl, kr)]
        sg = 1 if all(r[0] > 0 for r in ga) else (-1 if all(r[1] < 0 for r in ga) else 0)
        if not sg and dfn and mono: sg = zero_sign(Eb, Ea, kl, kr)
        if not (dfn and mono and sg): ok = False
        zeros.append((kl/8, kr/8, sg)); i = i2
    alt = ok and all(zeros[q][2]*zeros[q + 1][2] < 0 for q in range(len(zeros) - 1))
    reach = any(z[0] >= Tstar for z in zeros)
    return ok, alt, reach, zeros

SURR = {}
def small_box_job(args):
    si, p0, hpw, b0, hbw, depth = args
    pl, ph, npn, TN = STRIPS[si]
    kb = math.floor(b0 - hbw); ka = kb + 1
    Sb, Sa = SURR[(si, kb)], SURR[(si, ka)]
    Eb = Sb.enclose(p0, hpw, b0, hbw); Ea = Sa.enclose(p0, hpw, b0 + 1, hbw)
    nk = int(TN/H)
    Ts = far_T_cache[(si, p0 - hpw, p0 + hpw)]
    ok, alt, reach, zeros = certify(Eb, Ea, nk, Ts)
    if ok and alt and reach:
        return [dict(p=[float(p0 - hpw), float(p0 + hpw)], b=[float(b0 - hbw), float(b0 + hbw)], zeros=len(zeros),
                     first=zeros[0][:2] if zeros else None, Tstar=Ts, depth=depth, ok=True)]
    if depth >= 5:
        return [dict(p=[float(p0 - hpw), float(p0 + hpw)], b=[float(b0 - hbw), float(b0 + hbw)], zeros=len(zeros),
                     isolated=ok, alternate=alt, overlap=reach, depth=depth, ok=False)]
    out = []
    for dp in (-1, 1):
        for db in (-1, 1):
            nb0 = b0 + db*hbw/2
            if nb0 - hbw/2 >= p0 + dp*hpw/2 + hpw/2 - 1: continue   # wholly past the line b = p - 1
            out += small_box_job((si, p0 + dp*hpw/2, hpw/2, nb0, hbw/2, depth + 1))
    return out
far_T_cache = {}
class FT(dict):
    def __missing__(self, key):
        si, a, b = key; v = far_T(a, b); self[key] = v; return v
far_T_cache = FT()

def init_worker(surr):
    global SURR; SURR = surr

if __name__ == "__main__":
    BB = big_boxes()
    jobs = []; meta = []
    for bi, bb in enumerate(BB):
        up, uc = nodes(bb["np"]), nodes(bb["nc"])
        for i, u in enumerate(up):
            for k, v in enumerate(uc):
                p = (bb["pl"] + bb["ph"])/2 + u*(bb["ph"] - bb["pl"])/2
                c = (bb["cl"] + bb["ch"])/2 + v*(bb["ch"] - bb["cl"])/2
                jobs.append((p, c, bb["TN"])); meta.append((bi, i, k))
    print(f"{len(BB)} big boxes, {len(jobs)} nodes", flush=True)
    import pickle
    vals = [None]*len(jobs)
    for bi in reversed(range(len(BB))):   # one big box at a time (hardest first), cached, so a long run can resume
        idx = [q for q, m in enumerate(meta) if m[0] == bi]
        fn = f"/tmp/nz/nodes_{CFG['tag']}_{bi}.pkl"
        if os.path.exists(fn): got = pickle.load(open(fn, "rb"))
        else:
            with Pool(4) as pool: got = pool.map(node_eval, [jobs[q] for q in idx], chunksize=4)
            os.makedirs("/tmp/nz", exist_ok=True); pickle.dump(got, open(fn, "wb"))
        for q, g in zip(idx, got): vals[q] = g
        print(f"  big box {bi + 1}/{len(BB)} {time.perf_counter() - t0:.0f}s", flush=True)
    print(f"nodes done {time.perf_counter() - t0:.0f}s", flush=True)
    SURR = {}
    for bi, bb in enumerate(BB):
        nk = int(bb["TN"]/H) + 1
        V = np.zeros((bb["np"], bb["nc"], nk, DER)); R = np.zeros_like(V)
        for (b2, i, k), (m, r) in zip(meta, vals):
            if b2 == bi: V[i, k] = m; R[i, k] = r
        SURR[(bb["strip"], int(bb["cl"]))] = Surrogate(bb, V, R, nodes(bb["np"]), nodes(bb["nc"]))
    pickle.dump(SURR, open(f"/tmp/nz/surr_{CFG['tag']}.pkl", "wb"))
    OUT["interp_eps"] = {f"strip {s} c in [{k},{k + 1}]": S.eps for (s, k), S in SURR.items()}
    print("interpolation error bounds:", {k: f"{v:.1e}" for k, v in OUT["interp_eps"].items()}, flush=True)
    # small boxes
    W = Fr(1, 20)
    sjobs = []
    for si, (pl, ph, npn, TN) in enumerate(STRIPS):
        npieces = int((ph - pl)/W)
        for q in range(npieces):
            a, b = pl + q*W, pl + (q + 1)*W
            kb = 0
            while kb*W < b - 1:
                sjobs.append((si, (a + b)/2, W/2, kb*W + W/2, W/2, 0)); kb += 1
    print(f"{len(sjobs)} small boxes", flush=True)
    with Pool(4) as pool:              # forked workers share SURR copy-on-write
        res = pool.map(small_box_job, sjobs, chunksize=4)
    boxes = [b for r in res for b in r]
    bad = [b for b in boxes if not b["ok"]]
    OUT["boxes"] = boxes
    print(f"small boxes done {time.perf_counter() - t0:.0f}s", flush=True)
    check(f"Z1 every zero of G at b in [0, TN] is isolated, on boxes covering {float(CFG['lo'])} <= p <= {float(CFG['hi'])} and 0 <= b <= p - 1",
          not [b for b in bad if not b.get("isolated", True)],
          f"{len(boxes)} boxes after splitting, zeros per box {min(b['zeros'] for b in boxes)} to {max(b['zeros'] for b in boxes)}")
    check("Z2 G at a = b + 1 has a definite sign at every such zero, alternating", not [b for b in bad if not b.get("alternate", True)])
    check("Z3 every box reaches a zero past the far zone's T*, so near and far zones share a zero",
          not [b for b in bad if not b.get("overlap", True)])
    OUT["failed"] = bad
    # Z4: surrogate against direct quadrature
    random.seed(1); mp.mp.dps = 30; gaps = []
    for _ in range(12):
        si = random.randrange(len(STRIPS)); pl, ph, npn, TN = STRIPS[si]
        p = float(pl) + random.random()*float(ph - pl); b = random.random()*(p - 1); s = random.random()*float(TN)
        kb = math.floor(b); S = SURR[(si, kb)]
        pf, bf = Fr(p).limit_denominator(10**6), Fr(b).limit_denominator(10**6)
        lo, hi = S.enclose(pf, Fr(0), bf, Fr(0))
        k = int(s*8); sv = k/8
        pp, bq = mp.mpf(pf.numerator)/pf.denominator, mp.mpf(bf.numerator)/bf.denominator
        g = mp.quad(lambda x: (1 - x**pp)**bq*mp.cos(sv*x), mp.linspace(0, 1, 30))
        inside = lo[k, 0] - 1e-12 <= float(g) <= hi[k, 0] + 1e-12
        gaps.append(dict(p=p, b=b, s=sv, quad=float(g), lo=float(lo[k, 0]), hi=float(hi[k, 0]), inside=bool(inside)))
    OUT["Z4"] = gaps
    check("Z4 polish: the surrogate's enclosure contains the directly integrated G at 12 random points", all(g["inside"] for g in gaps),
          f"largest enclosure width {max(g['hi'] - g['lo'] for g in gaps):.1e}")
    # Z5: p = 3
    if CFG["lo"] <= 3 <= CFG["hi"]:
        at3 = [b for b in boxes if b["p"][0] <= 3 <= b["p"][1] and b["ok"]]
        OUT["Z5_p3_zero_counts"] = sorted(set(b["zeros"] for b in at3))
        check("Z5 boxes containing p = 3 are certified, as the separate p = 3 certificate says",
              bool(at3), f"{len(at3)} boxes, zero counts {OUT['Z5_p3_zero_counts']}")
    else:   # Z5 instead: the top step b = p - 1 is covered all the way along the range
        top = [b for b in boxes if b["ok"] and b["b"][1] >= b["p"][0] - 1]
        span = sorted((b["p"][0], b["p"][1]) for b in top); reach = float(CFG["lo"]); gaps = 0
        for a, c in span:
            if a > reach + 1e-12: gaps += 1
            reach = max(reach, c)
        check("Z5 the top step b = p - 1 is covered with no gap in p", gaps == 0 and reach >= float(CFG["hi"]) - 1e-12,
              f"{len(top)} boxes touch the top step")
    OUT["elapsed_s"] = time.perf_counter() - t0
    json.dump(OUT, open(CFG["tag"] + ".json", "w"), indent=1, default=str)
    print(f"{sum(RES)}/{len(RES)} pass, {OUT['elapsed_s']:.0f} s")
