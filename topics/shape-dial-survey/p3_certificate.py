"""Lee-Yang at roundness p = 3: a certificate that G_a(s) = int_0^1 (1-x^3)^a cos(sx) dx has only real zeros for
every 0 < a <= 3 (ball dimension n = 3a + 1 <= 10 = 1 + p^2).

Proof structure (README, section N3, p = 3):
  * a <= 1: Polya (1918).
  * Climb in whole steps b -> a = b + 1, for every b in [0, 2]: if G_b has only real simple zeros and G_a has
    alternating signs at the consecutive zeros of G_b, then G_a has only real (interlacing) zeros
    (Krein's partial-fraction theorem for Cartwright-class functions; G_a/G_b -> 0 up the imaginary axis).
  * The alternation is certified here, for the whole continuum b in [0, 2]:
      FAR  (s >= T2 = 21): explicit inequalities from the exact contour formula
             G_c(s) = C_c(s) + |J_c(s)| sin(s - pi c/2 + arg J_c(s)),
             C_c(s) = -int_0^inf Im(1+iy^3)^c e^(-sy) dy,  J_c(s) = int_0^inf y^c (3+3iy-y^2)^c e^(-sy) dy,
           with Watson-type remainder bounds; every bound is monotone in s, so checking s = T2 covers [T2, inf).
      NEAR (0 <= s <= 25): every zero of G_b isolated and the sign of G_a there fixed, on boxes of width 0.01 in b,
           using Taylor models in s (|d^j G/ds^j| <= 1/(j+1)) and in b (|d^2 G/db^2| <= L2, |d^3 G/db^3| <= L3).
Arithmetic: exact big-integer Horner at dyadic s, coefficients to 2^-900 by exact recurrences; interval
arithmetic (mpmath.iv) for the far-zone inequalities; float bookkeeping with explicit slack.
"""
import json, time, math
from fractions import Fraction as Fr
import mpmath as mp
iv = mp.iv
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
t0 = time.perf_counter()
L2, L3 = 0.75, 2.11            # int_0^1 log^2(1-x^3) = 0.74752..., int |log(1-x^3)|^3 = 2.10389...
T2, TN, H, NT = 21, 25, Fr(1, 8), 150   # far zone from T2; near zone to TN; cell width; series length
P = 900; DER = 7; SLACK = 1e-25

# ---------- FAR ZONE: interval check of the explicit inequalities at s = T2 over b-boxes ----------
def far_margins(bl, bh, t):
    t = iv.mpf(t); b = iv.mpf([bl, bh]); a = b + 1
    def S(c, extra):  # sum_j C(4,j) Gamma(c+3+j+extra) / (Gamma(c+1+extra) t^(2+j))
        return sum(math.comb(4, j)*iv.gamma(c + 3 + j + extra)/(iv.gamma(c + 1 + extra)*t**(2 + j)) for j in range(5))
    def absc(c):  # |c(c-1)| as an interval
        v = c*(c - 1); lo, hi = v.a, v.b
        return iv.mpf([0 if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi))])
    alpha = lambda c: (absc(c) + 2*c/3)/2
    Eb, Eb1, Ea = alpha(b)*S(b, 0), alpha(b)*S(b, 1), alpha(a)*S(a, 0)
    RC = lambda c: absc(c)/2*(720/t**7 + 362880/t**10)
    RC1 = lambda c: absc(c)/2*(5040/t**8 + 3628800/t**11)
    Jb = 3**b*iv.gamma(b + 1)*t**(-b - 1)*(1 - Eb); Ja = 3**a*iv.gamma(a + 1)*t**(-a - 1)*(1 - Ea)
    Cb, Cb1, Ca = 6*b/t**4 + RC(b), 24*b/t**5 + RC1(b), 6*a/t**4 + RC(a)
    q = (b/t + Eb + Eb1)/(1 - Eb); betap = (b + 1)/t*q; ReJ = (b + 1)/t*(1 + q)
    eps = Cb/Jb; epsp = (Cb1 + Cb*ReJ)/Jb
    asb = lambda x: x/iv.sqrt(1 - x**2)          # >= arcsin x
    dl = 2*a/t + asb(Ea) + asb(Eb)
    A = 1 - eps
    B = (1 - betap)*iv.sqrt(1 - eps**2) - epsp
    C = iv.sqrt(1 - eps**2)*iv.cos(dl) - eps*iv.sin(dl) - Ca/Ja   # cos(arcsin eps + dl) - |C_a|/|J_a|
    return float(min(A.a, 1)), float(B.a), float(C.a), float(Eb.b), float(Ea.b)
iv.dps = 30
far = []
nb = 40
for i in range(nb):
    bl, bh = 2*i/nb, 2*(i + 1)/nb
    far.append(dict(b=[bl, bh], margins=far_margins(bl, bh, T2)[:3]))
OUT["far"] = dict(T2=T2, boxes=far, worst=[min(f["margins"][k] for f in far) for k in range(3)])
check(f"F1 far zone: for every b in [0, 2] and every s >= {T2}, the three inequalities hold (|eps| < 1, one simple zero "
      "per half-turn, and G at a = b + 1 keeps the sign the phase predicts)",
      all(min(f["margins"]) > 0 for f in far),
      f"worst margins {', '.join(f'{w:.3f}' for w in OUT['far']['worst'])} over {nb} boxes")

# ---------- NEAR ZONE: exact-arithmetic Taylor models ----------
mp.mp.prec = 1000
def moments(c):  # m_n(c) = B((n+1)/3, c+1)/3 for n = 0..NT+DER, c rational; m_(n+3) = m_n (n+1)/(n+3c+4)
    iv.prec = 1000; cc = iv.mpf(c.numerator)/c.denominator
    base = [iv.gamma(iv.mpf(n + 1)/3)*iv.gamma(cc + 1)/(3*iv.gamma(iv.mpf(n + 1)/3 + cc + 1)) for n in range(3)]
    m = [mp.mpf(x.mid) for x in base]; cm = mp.mpf(c.numerator)/c.denominator
    for n in range(3, NT + DER + 1): m.append(m[n - 3]*(n - 2)/(n - 2 + 3*cm + 3))
    return m
FACT = [mp.factorial(n) for n in range(NT + 1)]
def deriv_ints(m):  # integer coefficients of d^j G/ds^j = sum_n Re(i^(n+j)) s^n m_(n+j)/n!, scaled by 2^P
    out = []
    for j in range(DER):
        row = []
        for n in range(NT + 1):
            k = (n + j) % 4; sgn = 1 if k == 0 else (-1 if k == 2 else 0)
            row.append(int(mp.nint(sgn*m[n + j]/FACT[n]*mp.mpf(2)**P)) if sgn else 0)
        out.append(row)
    return out
def horner(row, num, den):  # exact: sum row_n num^n den^(NT-n), then divide
    acc = 0; dpow = 1
    for n in range(NT, -1, -1):
        acc = acc*num + row[n]*dpow; dpow *= den
    return acc
def err_bound(t):  # truncation + rounding, per value
    t = max(float(t), 1.0)
    tail = math.exp((NT + 1)*math.log(t) - math.lgamma(NT + 2))/(1 - t/(NT + 2))
    rnd = 2.0**(-P + 1)*(NT + 1)*math.exp(NT*math.log(t))*1.5
    return tail + rnd + SLACK
def values(rows, t):  # list over j of float approximations at rational t (dyadic)
    num, den = t.numerator, t.denominator; scale = mp.mpf(2)**P*mp.mpf(den)**NT
    return [float(mp.mpf(horner(r, num, den))/scale) for r in rows]
ETA = Fr(1, 10**12)
def box_enclosures(c0, delta, ts):
    """For G_c with c in [c0 - delta, c0 + delta]: at each grid s, enclosure (mid, rad) of d^j G/ds^j, j < DER."""
    R0, Rp, Rm = (deriv_ints(moments(c)) for c in (c0, c0 + ETA, c0 - ETA))
    d, eta = float(delta), float(ETA); res = {}
    for t in ts:
        e = err_bound(t); v0, vp, vm = values(R0, t), values(Rp, t), values(Rm, t)
        enc = []
        for j in range(DER):
            D = (vp[j] - vm[j])/(2*eta); Drad = eta**2*L3/6 + e/eta
            enc.append((v0[j], d*(abs(D) + Drad) + d*d*L2/2 + e + 1e-13*abs(v0[j])))
        res[t] = enc
    return res
def cell_range(enc, j, h):  # enclosure of d^j G/ds^j over [t, t + h], from the Taylor data at t
    h = float(h); mid, rad = enc[j]
    var = sum((abs(enc[j + m][0]) + enc[j + m][1])*h**m/math.factorial(m) for m in range(1, DER - j))
    var += h**(DER - j)/(math.factorial(DER - j)*(DER + 1))
    return mid - rad - var, mid + rad + var
def certify_box(b0, delta):
    ts = [H*k for k in range(int(TN/H) + 1)]
    Eb = box_enclosures(b0, delta, ts); Ea = box_enclosures(b0 + 1, delta, ts)
    cells = []
    for t in ts[:-1]:
        lo, hi = cell_range(Eb[t], 0, H); cells.append((t, lo > 0 or hi < 0))
    zeros, k, ok = [], 0, True
    while k < len(cells):
        if cells[k][1]: k += 1; continue
        k2 = k
        while k2 < len(cells) and not cells[k2][1]: k2 += 1
        tl, tr = cells[k][0], cells[k2 - 1][0] + H
        sl = Eb[tl][0]; sr = Eb[tr][0]
        definite = (abs(sl[0]) > sl[1]) and (abs(sr[0]) > sr[1]) and (sl[0]*sr[0] < 0)
        mono = all((lambda r: r[0] > 0 or r[1] < 0)(cell_range(Eb[cells[i][0]], 1, H)) for i in range(k, k2))
        ga = [cell_range(Ea[cells[i][0]], 0, H) for i in range(k, k2)]
        sgn = 1 if all(r[0] > 0 for r in ga) else (-1 if all(r[1] < 0 for r in ga) else 0)
        if not (definite and mono and sgn): ok = False
        zeros.append(dict(cell=[float(tl), float(tr)], Ga_sign=sgn)); k = k2
    alt = ok and all(zeros[i]["Ga_sign"]*zeros[i + 1]["Ga_sign"] < 0 for i in range(len(zeros) - 1))
    reach = any(z["cell"][0] >= T2 for z in zeros)
    return dict(b=[float(b0 - delta), float(b0 + delta)], n_zeros=len(zeros), isolated=ok, alternate=alt,
                past_T2=reach, first_zero=zeros[0]["cell"] if zeros else None, zeros=[z["cell"] for z in zeros])
DELTA = Fr(1, 200); NBOX = int(2/(2*DELTA))
boxes = []
for i in range(NBOX):
    b0 = DELTA*(2*i + 1); r = certify_box(b0, DELTA); boxes.append(r)
    if i % 25 == 0 or not (r["isolated"] and r["alternate"]): print("   box", r["b"], r["n_zeros"], r["isolated"], r["alternate"], f"{time.perf_counter()-t0:.0f}s", flush=True)
OUT["near"] = dict(TN=TN, cell=float(H), delta=float(DELTA), boxes=boxes)
check(f"N1 near zone: on all {NBOX} boxes of width 0.01 covering b in [0, 2], every zero of G_b in [0, {TN}] is isolated "
      "(sign change plus nonzero slope on its cells)", all(r["isolated"] for r in boxes),
      f"{min(r['n_zeros'] for r in boxes)} to {max(r['n_zeros'] for r in boxes)} zeros per box")
check("N2 near zone: G at a = b + 1 has a definite sign at each of those zeros, and the signs alternate",
      all(r["alternate"] for r in boxes))
check(f"N3 overlap: in every box the near zone reaches a zero past s = {T2}, so near and far zones share a zero",
      all(r["past_T2"] for r in boxes))
# ---------- P1: sanity, exact contour formula vs the defining integral ----------
mp.mp.dps = 30
def Gq(c, s): return mp.quad(lambda x: (1 - x**3)**c*mp.cos(s*x), [0, 1])
def Jq(c, s): return mp.quad(lambda y: y**c*(3 + 3j*y - y*y)**c*mp.exp(-s*y), [0, 1/s, 5/s, 30/s, mp.inf])
def Cq(c, s): return -mp.quad(lambda y: mp.im((1 + 1j*y**3)**c)*mp.exp(-s*y), [0, 1/s, 5/s, 30/s, mp.inf])
gaps = []
for c, s in [(0.5, 4.0), (1.5, 9.3), (2.0, 21.0), (3.0, 21.0)]:
    Jv = Jq(c, s); gaps.append(float(abs(Gq(c, s) - Cq(c, s) - abs(Jv)*mp.sin(s - mp.pi*c/2 + mp.arg(Jv)))))
OUT["contour_gap"] = max(gaps)
check("P1 polish: the contour formula reproduces the defining integral (4 spot checks)", max(gaps) < 1e-25, f"gap {max(gaps):.0e}")
OUT["seconds"] = time.perf_counter() - t0
print(f"{sum(RES)}/{len(RES)} pass, {OUT['seconds']:.0f} s")
json.dump(OUT, open("p3_certificate.json", "w"), indent=1, default=float)
