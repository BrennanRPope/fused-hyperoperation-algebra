"""Lee-Yang far zone as formulas in the roundness p (DRAFT; interval arithmetic, hand-derived bounds, not checked in Lean).

The climb from step b to step a = b + 1 (a <= p) needs three inequalities at the zeros of G_b. The p = 3
certificate proved them for s >= 21 with constants worked out for p = 3. Here every constant is a formula in
(p, c, s), with two free split parameters kappa (inner window y1 = kappa/(p+2)) and lam (tail split).

  J_c(s)  = p^c Gamma(c+1) (s - i beta)^-(c+1) (1 + E),   beta = c (p-1)/2,            |E| <= Ebound
  Ebound  = Kt (k+1)(k+2) (R/s)^(k+1) / s^2                                            inner window [0, y1]
          + e^(-lam s y1) (R/sigma)^(k+1) (1 + min(S1, S2))                            outer tail
      R = sqrt(s^2 + beta^2), sigma = (1 - lam) s, k = c (+1 for J'), N = ceil(c (p-1)),
      S1 = sum_j C(N,j) (k+1)...(k+j)/sigma^j,  S2 = sum_j C(N,j) (k+1)...(k+2j)/(2 sigma^2)^j
      Kt = c H exp(c max(0, eta2 + K3 y1) y1^2),  H = |eta2| + K3 y1,  eta2 = (p-1)(5-p)/24
      K3 = K2 + a1 V + V^2 y1/2 + (a1 + V y1)^3/(3(1 - Z)),  a1 = (p-1)/2,  a2 = |(p-1)(p-2)|/6,
      K2 = |(p-1)(p-2)(p-3)|/24 max(1, (1+y1^2)^((p-4)/2)),  V = a2 + K2 y1,  Z = (a1 + V y1) y1 < 1/2
  |C_c(s)| <= c|sin(pi p/2)| Gamma(p+1)/s^(p+1) + |c(c-1)|/2 mu(c) Gamma(2p+1)/s^(2p+1) + e^(-lam s 2^(-1/p)) (...)
The three far-zone margins (A, B, C) are then lower-bounded by formulas in these. Every bound decreases in s once
s >= (p+1)/(lam 2^(-1/p)), so positive margins at one s = T certify every s >= T.

  F1 checks every bound against the exact J, J', C, C' (quadrature) across roundness, steps and s;
  F2 certifies the takeover T*(p) on boxes covering 1.1 <= p <= 8 and every step 0 <= b <= p-1, with
     interval arithmetic (each box: all p and all b in it), splitting boxes that need it;
  F3 certifies T* at sample roundness outside that range (p near 1, and p = 10, 12, 16), on step boxes;
  F4 compares T*(p) with p^2, with the numerical takeover from general_p_far.json and with the p = 3 certificate's 21.
"""
import math, time, json, sys
import mpmath as mp
from multiprocessing import Pool
iv = mp.iv
iv.dps = 20
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
t0 = time.perf_counter()

# ---------------- interval helpers ----------------
def I(x): return x if isinstance(x, type(iv.mpf(1))) else iv.mpf(x)
def imax(x, y): x, y = I(x), I(y); return iv.mpf([max(x.a, y.a), max(x.b, y.b)])
def iabs(x):
    x = I(x); lo, hi = x.a, x.b
    return iv.mpf([0 if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi))])
X0 = mp.mpf("1.46163214496836234126")           # where Gamma has its minimum on (0, inf)
GMIN = iv.mpf("0.885603194410888")               # just below Gamma(X0) = 0.88560319441088870...
def igamma(x):
    """Tight enclosure of Gamma on an interval x > 0: Gamma decreases before X0 and increases after."""
    x = I(x); lo, hi = iv.gamma(iv.mpf(x.a)), iv.gamma(iv.mpf(x.b))
    if x.b <= X0: return iv.mpf([hi.a, lo.b])
    if x.a >= X0: return iv.mpf([lo.a, hi.b])
    return iv.mpf([GMIN.a, max(lo.b, hi.b)])
def mu(c): return imax(iv.mpf(2)**(2 - c), iv.mpf(1.5)**(c - 2))
def rising(k, j):
    out = iv.mpf(1)
    for i in range(1, j + 1): out = out*(k + i)
    return out

# ---------------- the bounds ----------------
def J_params(p, c, kappa):
    y1 = kappa/(p + 2)
    a1 = (p - 1)/2; a2 = iabs((p - 1)*(p - 2))/6
    K2 = iabs((p - 1)*(p - 2)*(p - 3))/24*imax(1, (1 + y1**2)**((p - 4)/2))   # |M - 1 - i a1 y + a2 y^2| <= K2 y^3
    V = a2 + K2*y1                                                             # |M - 1 - i a1 y| <= V y^2
    Z = (a1 + V*y1)*y1                                                         # |M - 1| <= Z
    K3 = K2 + a1*V + V**2*y1/2 + (a1 + V*y1)**3/(3*(1 - Z))                    # |log M - i a1 y - eta2 y^2| <= K3 y^3
    eta2 = (p - 1)*(5 - p)/24
    H = iabs(eta2) + K3*y1
    Kt = c*H*iv.exp(c*imax(0, eta2 + K3*y1)*y1**2)                             # |M^c - e^(i beta y)| <= Kt y^2
    return y1, Kt, Z
def E_bound(p, c, s, e, kappa, lam):
    y1, Kt, Z = J_params(p, c, kappa)
    if not Z.b < 0.5: return None
    k = c + e; beta = c*(p - 1)/2; R = iv.sqrt(s**2 + beta**2)
    inner = Kt*(k + 1)*(k + 2)*(R/s)**(k + 1)/s**2
    sig = (1 - lam)*s
    N = math.ceil(float((c*(p - 1)).b))
    S1 = sum(math.comb(N, j)*rising(k, j)/sig**j for j in range(N + 1))           # |M| <= (1 + y)^(p-1)
    S2 = sum(math.comb(N, j)*rising(k, 2*j)/(2*sig**2)**j for j in range(N + 1))  # |M| <= (1 + y^2/2)^(p-1)
    S = S2 if S2.b < S1.b else S1
    tail = iv.exp(-lam*s*y1)*(R/sig)**(k + 1)*(S + 1)
    return inner + tail
def C_bound(p, c, s, e, lam):
    th = iv.pi*p/2; y2 = iv.mpf(2)**(-1/p); sig = (1 - lam)*s
    Lc = imax(1, iv.mpf(2)**(c - 1))
    lead = c*iabs(iv.sin(th))*igamma(p + 1 + e)/s**(p + 1 + e)
    second = iabs(c*(c - 1))/2*mu(c)*igamma(2*p + 1 + e)/s**(2*p + 1 + e)
    tail = iv.exp(-lam*s*y2)*((Lc + 1)*igamma(1 + e)/sig**(1 + e) + Lc*igamma(p*c + 1 + e)/sig**(p*c + 1 + e)
                              + c*igamma(p + 1 + e)/sig**(p + 1 + e))
    return lead + second + tail
def asb(x): return x/iv.sqrt(1 - x**2)
def env_s_over(s, beta):        # sup over s' >= s of s'/(s'^2 + beta^2)
    if s.a >= beta.b: return s/(s**2 + beta**2)
    return imax(s/(s**2 + beta**2), 1/(2*imax(beta, 1e-30)))
def margins(p, b, s, kappa, lam):
    a = b + 1
    Eb, Eb1, Ea = E_bound(p, b, s, 0, kappa, lam), E_bound(p, b, s, 1, kappa, lam), E_bound(p, a, s, 0, kappa, lam)
    if Eb is None or Ea is None or Eb1 is None: return None
    if not (Eb.b < 0.5 and Ea.b < 0.5): return None
    bb, ba = b*(p - 1)/2, a*(p - 1)/2
    Jb = p**b*igamma(b + 1)*(s**2 + bb**2)**(-(b + 1)/2)*(1 - Eb)
    Ja = p**a*igamma(a + 1)*(s**2 + ba**2)**(-(a + 1)/2)*(1 - Ea)
    Cb, Cb1, Ca = C_bound(p, b, s, 0, lam), C_bound(p, b, s, 1, lam), C_bound(p, a, s, 0, lam)
    q = (Eb + Eb1)/(1 - Eb)
    Rb = iv.sqrt(s**2 + bb**2)
    psip = 1 - (b + 1)*bb/(s**2 + bb**2) - (b + 1)*q/Rb          # lower bound on psi'
    ReJ = (b + 1)*(env_s_over(s, bb) + q/Rb)                       # upper bound on |Re J'/J|
    eps = Cb/Jb
    if not eps.b < 1: return None
    epsp = (Cb1 + Cb*ReJ)/Jb
    k = (p - 1)/(2*s)
    hk = (a + 1)*k/(1 + (k*b)**2)
    if not (k*b).b <= 1: hk = imax(hk, (a + 1)/(2*imax(b, 1e-30)))   # envelope, so the bound decreases in s
    dl = k*a + hk + asb(Ea) + asb(Eb)                               # atan x <= x
    if not dl.b < 1.55: return None
    A = 1 - eps
    B = psip*iv.sqrt(1 - eps**2) - epsp
    C = iv.sqrt(1 - eps**2)*iv.cos(dl) - eps*iv.sin(dl) - Ca/Ja
    return float(A.a), float(B.a), float(C.a)
def monotone_from(p, lam):
    return float(((p + 1)/(lam*iv.mpf(2)**(-1/p))).b)
PARAMS = [(k, l) for k in (0.5, 0.75, 1, 1.25, 1.5, 2) for l in (0.4, 0.5, 0.6, 0.7)]
TS = [2*1.04**i for i in range(700)]
def ok_at(P, B, T, params=PARAMS):
    for kap, lam in params:
        if T < monotone_from(P, lam): continue
        m = margins(P, B, iv.mpf(T), kap, lam)
        if m is not None and min(m) > 0: return (kap, lam, m)
    return None
def minT(P, B, params=PARAMS, cap=None):
    """Smallest grid T with all three margins positive (bisection over the grid; margins grow with T)."""
    hi = len(TS) - 1 if cap is None else min(len(TS) - 1, next(i for i, t in enumerate(TS) if t >= cap))
    if ok_at(P, B, TS[hi], params) is None: return None
    lo = 0
    while lo < hi:
        mid = (lo + hi)//2
        if ok_at(P, B, TS[mid], params) is not None: hi = mid
        else: lo = mid + 1
    return TS[lo]

# ---------------- F1: bounds against exact values ----------------
def logw(y, p):
    q = mp.exp(-p*mp.log(1 + 1j*y))
    return -1j*mp.pi + p*mp.log(1 + 1j*y) + mp.log(1 - q) - (mp.log(y) - 1j*mp.pi/2)
def Jx(c, s, p, d=0):
    return mp.quad(lambda y: (-y)**d*y**c*mp.exp(c*logw(y, p))*mp.exp(-s*y), [0, 1/s, 5/s, 20/s, 60/s, 200/s])
def Cx(c, s, p, d=0):
    th = mp.pi*p/2
    f = lambda y: (-y)**d*mp.im((1 - mp.exp(1j*th)*y**p)**c)*mp.exp(-s*y)
    try: return -mp.quad(f, [0, 1/s, 5/s, 20/s, 60/s, 200/s])
    except ZeroDivisionError: return -mp.quad(f, [0, 1/s, 5/s, 20/s, 60/s, 200/s], method="gauss-legendre")
def f1_one(p):
    mp.mp.dps = 30
    worst = {"E": 0.0, "E1": 0.0, "C": 0.0, "C1": 0.0}; bad = []
    for c in [0.2, 0.7, 1.0, 1.6, 2.5, 4.0, 6.0, 8.0]:
        if c > p: continue
        beta = c*(p - 1)/2
        for s in [10.0, 30.0, 100.0, 300.0]:
            Et = abs(Jx(c, s, p)/(p**c*mp.gamma(c + 1)*(s - 1j*beta)**(-c - 1)) - 1)
            E1t = abs(-Jx(c, s, p, 1)/(p**c*mp.gamma(c + 2)*(s - 1j*beta)**(-c - 2)) - 1)
            P, Cc, S = iv.mpf(p), iv.mpf(c), iv.mpf(s)
            for kap, lam in [(1, 0.5), (0.5, 0.6)]:
                Eb, E1b = E_bound(P, Cc, S, 0, kap, lam), E_bound(P, Cc, S, 1, kap, lam)
                if Eb is None: continue
                Cb, C1b = C_bound(P, Cc, S, 0, lam), C_bound(P, Cc, S, 1, lam)
                for key, t, bd in [("E", Et, Eb), ("E1", E1t, E1b), ("C", abs(Cx(c, s, p)), Cb), ("C1", abs(Cx(c, s, p, 1)), C1b)]:
                    r = float(t/bd.a); worst[key] = max(worst[key], r)
                    if r > 1: bad.append((key, p, c, s, kap, lam, float(t), float(bd.a)))
    return p, worst, bad

# ---------------- F2: boxes ----------------
def box_T(pl, ph, fl, fh, depth=0):
    """Certified takeover on the box p in [pl, ph], b = f (p-1) with f in [fl, fh]; splits when needed."""
    P = iv.mpf([pl, ph]); B = iv.mpf([fl, fh])*(P - 1)
    pm = (pl + ph)/2; Pm, Bm = iv.mpf(pm), iv.mpf((fl + fh)/2)*(pm - 1)
    T0 = minT(Pm, Bm) or TS[-1]
    hit = ok_at(Pm, Bm, T0)
    order = PARAMS if hit is None else [hit[:2]] + [q for q in PARAMS if q != hit[:2]]
    i0 = next(i for i, t in enumerate(TS) if t >= T0)
    for i in range(i0, min(i0 + 12, len(TS))):
        r = ok_at(P, B, TS[i], order)
        if r is not None:
            return [dict(p=[pl, ph], f=[fl, fh], T=TS[i], kappa=r[0], lam=r[1], margins=r[2], depth=depth)]
    if depth >= 6: return [dict(p=[pl, ph], f=[fl, fh], T=None, depth=depth)]
    if (ph - pl)/max(pl - 1, 1e-9) > (fh - fl)/2:
        m = (pl + ph)/2; return box_T(pl, m, fl, fh, depth + 1) + box_T(m, ph, fl, fh, depth + 1)
    m = (fl + fh)/2; return box_T(pl, ph, fl, m, depth + 1) + box_T(pl, ph, m, fh, depth + 1)
def f2_strip(args):
    pl, ph = args
    out = []
    for j in range(8): out += box_T(pl, ph, j/8, (j + 1)/8)
    return out

if __name__ == "__main__":
    PS1 = [1.25, 1.5, 2.5, 3.0, 3.7, 5.0, 7.0]
    with Pool(4) as pool: r1 = pool.map(f1_one, PS1)
    worst = {k: max(w[k] for _, w, _ in r1) for k in ["E", "E1", "C", "C1"]}
    bad = [b for _, _, bs in r1 for b in bs]
    OUT["F1_worst_ratio"] = worst; OUT["F1_violations"] = bad
    check("F1 every bound holds against exact J, J', C, C' (p = 1.25 to 7, c up to 8, s = 10 to 300)", not bad,
          "largest actual/bound " + ", ".join(f"{k} {v:.3f}" for k, v in worst.items()))
    edges = [round(1.1 + 0.05*i, 4) for i in range(18)] + [round(2.0 + 0.1*i, 4) for i in range(61)]
    strips = list(zip(edges[:-1], edges[1:]))
    with Pool(4) as pool: r2 = pool.map(f2_strip, strips)
    boxes = [b for st in r2 for b in st]
    failed = [b for b in boxes if b["T"] is None]
    per_strip = [dict(p=list(st), T=max(b["T"] for b in r) if all(b["T"] for b in r) else None, boxes=len(r))
                 for st, r in zip(strips, r2)]
    OUT["F2_boxes"] = boxes; OUT["F2_strips"] = per_strip
    check("F2 the takeover T*(p) is certified on boxes covering 1.1 <= p <= 8 and every step 0 <= b <= p - 1",
          not failed, f"{len(boxes)} boxes, {len(failed)} uncertified, every T at or past the monotone point")
    # F3
    samples = [1.01, 1.03, 1.05, 10.0, 12.0, 16.0]
    def f3(p):
        bx = []
        for j in range(8): bx += box_T(p, p, j/8, (j + 1)/8)
        return p, [b["T"] for b in bx]
    with Pool(4) as pool: r3 = pool.map(f3, samples)
    OUT["F3"] = {str(p): Ts for p, Ts in r3}
    check("F3 T* certified at p = 1.01, 1.03, 1.05, 10, 12, 16 on step boxes (split where needed)", all(all(Ts) for _, Ts in r3),
          ", ".join(f"p {p}: {max(Ts):.0f}" for p, Ts in r3 if all(Ts)))
    # F4
    num = json.load(open("general_p_far.json"))["scan"]
    rows = []
    for st in per_strip:
        p_hi = st["p"][1]; rows.append((st["p"], st["T"], st["T"]/p_hi**2 if st["T"] else None))
    big = [(st["p"], st["T"]) for st in per_strip if st["p"][0] >= 3]
    ratio = max(T/(pr[1])**2 for pr, T in big) if all(T for _, T in big) else None
    cmp = {}
    for key, v in num.items():
        p = float(key)
        st = next((s for s in per_strip if s["p"][0] <= p <= s["p"][1]), None)
        Tc = max(s["T"] for s in per_strip if s["p"][0] <= p <= s["p"][1]) if st else OUT["F3"].get(key) and max(OUT["F3"][key])
        cmp[key] = dict(numeric_takeover=v["takeover"], certified=Tc)
    OUT["F4_compare"] = cmp; OUT["F4_ratio_p2"] = ratio
    T3 = max(s["T"] for s in per_strip if s["p"][0] <= 3.0 <= s["p"][1])
    check("F4 for 3 <= p <= 8 the certified takeover stays under 4 p^2; at p = 3 it is near the p = 3 certificate's 21",
          ratio is not None and ratio < 4, f"largest T*/p^2 {ratio:.2f}; p = 3: {T3:.1f}")
    OUT["elapsed_s"] = time.perf_counter() - t0
    json.dump(OUT, open("far_bounds_p.json", "w"), indent=1, default=str)
    print(f"{sum(RES)}/{len(RES)} pass, {OUT['elapsed_s']:.0f} s")
