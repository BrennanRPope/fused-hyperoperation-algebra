"""Round 10 (DRAFT): the exact size of a break, proved without counting (README section 14).

Setting (as in Round 9): centres 0, 1, 2 on the line; class 1 = {lam, mu} about 0, class 2 = {nu} about 1 (letter N),
class 3 = {kappa} about 2 (letter K); all multipliers in (0, 1). The collision surface is
    lam = mu (2 - kappa),  nu = kappa / (2 - kappa),
where (a lam)(nu w) = (a mu)(kappa w) for every class-1 run a and every N/K word w.

Theorem P12 (proved in README section 14; this script checks its hypotheses and computes the numbers):
  on the collision surface, if
     (I)  the three hull images are disjoint: 2 lam < 1 - nu and 1 + nu < 2 - 2 kappa,
     (II) kappa < (3 - sqrt 5)/2, i.e. (2 - kappa)(1 - kappa) > 1,
     (P)  some odd prime p divides the numerators of both mu and kappa (mu, kappa rational),
  then the reduced pieces S' (no mu-containing class-1 run followed by K) give pairwise different maps, word by word,
  and dim K = min(1, s_F) EXACTLY, where s_F is the root of
     (P1 - 1)(Q - 1) - P1 Q (mu kappa)^s = 1,   P1 = 1/((1 - lam^s)(1 - mu^s)),  Q = 1/(1 - nu^s - kappa^s).

This script
  E1  the 3-adic break mu = 3/80, kappa = 3/40 (lam = 231/3200, nu = 3/77): exact collision, valuations, (I), (II),
      the Lemma A' disc version too; certified D*, certified exact dimension s_F, the drop.
  E2  sanity only (the proof does not use it): all reduced words up to length 8 give different maps (exact fractions),
      and the distinct-map counts of the 4-map system match the family model.
  E3  where Theorem P12 applies on the surface: grid of (mu, kappa), conditions (I), (II), s_F and D* (floats, figure).
  E4  points with the exact answer converging to the original break: mu_p = p/(20p + 1), kappa_p = p/(10p + 1).
  E5  the original break (19/200, 1/20, 1/19, 1/10): (P) fails (no shared prime). The same argument still works for
      the reduced pieces whose lam-power is at most k, as long as no ratio lam^a / mu^b (1 <= a <= k) lands in the
      two bad windows; that gives a certified lower bound within about 1e-6 of s_F. Not exact.
Writes exact.json."""
import json, math, itertools
from fractions import Fraction as F
import numpy as np
from mpmath import iv
import breaks as Bk
import general as G
import plane as P
iv.dps = 30
ivq = Bk.ivq

def surface(mu, kap):
    mu, kap = F(mu), F(kap)
    return mu * (2 - kap), mu, kap / (2 - kap), kap   # lam, mu, nu, kap

def vp(q, p):
    q = F(q); v = 0; n, d = q.numerator, q.denominator
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v

def odd_shared_primes(mu, kap):
    a, b = F(mu).numerator, F(kap).numerator
    g = math.gcd(a, b); out = []; p = 3
    while g > 1 and p <= g:
        if g % p == 0:
            out.append(p)
            while g % p == 0: g //= p
        p += 2
    if g > 1 and g % 2 == 1: out.append(g)
    return out

def conditions(lam, mu, nu, kap):
    lam, mu, nu, kap = map(F, (lam, mu, nu, kap))
    return dict(in_unit=all(0 < x < 1 for x in (lam, mu, nu, kap)),
                I_class1_vs_N=2 * lam < 1 - nu, I_N_vs_K=1 + nu < 2 - 2 * kap,
                II=(2 - kap) * (1 - kap) > 1,
                on_surface=(lam * nu == mu * kap and lam * (1 - nu) == 2 * mu * (1 - kap)))

# interval-arithmetic weights
def Dfun(s, lam, mu, nu, kap):        # D* is the root of (1 - lam^s)(1 - mu^s) - nu^s - kap^s  (increasing in s)
    return (1 - ivq(lam)**s) * (1 - ivq(mu)**s) - ivq(nu)**s - ivq(kap)**s

def reduced_weight(s, lam, mu, nu, kap, k=None):
    """sum over the reduced pieces S' of |piece|^s; with k, only class-1 runs whose lam-power is at most k.
       pieces u.v: u = lam^al mu^be nonempty; be = 0 -> any nonempty v; be >= 1 -> v starts with N."""
    L, M_, N_, K_ = ivq(lam)**s, ivq(mu)**s, ivq(nu)**s, ivq(kap)**s
    Lam = 1 / (1 - L) if k is None else sum((L**a for a in range(k + 1)), iv.mpf(0))
    Mu = 1 / (1 - M_); Q = 1 / (1 - N_ - K_)
    return (Lam - 1) * (Q - 1) + Lam * (Mu - 1) * N_ * Q

def family_weight(s, lam, mu, nu, kap):   # f - P1 Q (mu kap)^s, the Round 9 form (must agree with reduced_weight)
    P1 = 1 / ((1 - ivq(lam)**s) * (1 - ivq(mu)**s)); Q = 1 / (1 - ivq(nu)**s - ivq(kap)**s)
    return (P1 - 1) * (Q - 1) - P1 * Q * (ivq(mu) * ivq(kap))**s

def sQ_float(nu, kap):
    from scipy.optimize import brentq
    return brentq(lambda s: 1 - float(nu)**s - float(kap)**s, 1e-9, 5)

def certified_roots(lam, mu, nu, kap, k=None):
    sq = sQ_float(nu, kap)
    lo0 = sq + 1e-7
    Dlo, Dhi = Bk.certify_decreasing(lambda s: -Dfun(s, lam, mu, nu, kap), 0.02, 1.5)
    Flo, Fhi = Bk.certify_decreasing(lambda s: reduced_weight(s, lam, mu, nu, kap, k) - 1, lo0, 1.5)
    return (Dlo, Dhi), (Flo, Fhi)

def hom(c, l): return (F(l), F(c) * (1 - F(l)))
comp = Bk.comp

def reduced_words_distinct(lam, mu, nu, kap, n):
    """sanity: every word over the reduced pieces of total length <= n gives a different map (exact fractions)"""
    H = {"l": hom(0, lam), "m": hom(0, mu), "N": hom(1, nu), "K": hom(2, kap)}
    pieces = []
    for ln in range(2, n + 1):
        for al in range(0, ln):
            for be in range(0, ln - al):
                if al + be == 0: continue
                vlen = ln - al - be
                if vlen < 1: continue
                for v in itertools.product("NK", repeat=vlen):
                    if be >= 1 and v[0] == "K": continue
                    pieces.append("l" * al + "m" * be + "".join(v))
    maps = {}
    words_by_len = {0: [("", (F(1), F(0)))]}
    total = 0
    for ln in range(1, n + 1):
        cur = []
        for pc in pieces:
            if len(pc) > ln: continue
            pm = (F(1), F(0))
            for a in pc: pm = comp(pm, H[a])
            for w, wm in words_by_len.get(ln - len(pc), []):
                cur.append((pc + "|" + w, comp(pm, wm)))   # piece first (outer), then the rest
        words_by_len[ln] = cur
        for w, m in cur:
            total += 1
            if m in maps: return False, total, (maps[m], w)
            maps[m] = w
    return True, total, None

if __name__ == "__main__":
    OUT = {}
    # E1: the 3-adic break
    lam, mu, nu, kap = surface(F(3, 80), F(3, 40))
    assert (lam, nu) == (F(231, 3200), F(3, 77))
    cond = conditions(lam, mu, nu, kap)
    vals = {name: vp(x, 3) for name, x in zip(["lam", "mu", "nu", "kap"], (lam, mu, nu, kap))}
    x = comp(hom(0, lam), hom(1, nu)); y = comp(hom(0, mu), hom(2, kap))
    eps = max(lam, mu, nu, kap); R = 2 * (1 + eps) / (1 - eps); r2 = eps * (R + 1)
    disc = dict(eps=str(eps), R=str(R), r2=str(r2), threshold=str((1 - kap) / (3 - kap)),
                lemmaA=eps * 9 < 1, r2_below=(r2 < (1 - kap) / (3 - kap)))
    (Dlo, Dhi), (Flo, Fhi) = certified_roots(lam, mu, nu, kap)
    # agreement of the two weight formulas at a few points
    agree = all(abs(float((reduced_weight(iv.mpf(s), lam, mu, nu, kap) - family_weight(iv.mpf(s), lam, mu, nu, kap)).mid)) < 1e-20
                for s in (0.4, 0.5, 0.7))
    # lower bound for comparison: the collision-free subsystem with one class-1 map (Theorem P8 type)
    lb = max(Bk.certify_decreasing(lambda s, c=c: -(1 - ivq(c)**s - ivq(nu)**s - ivq(kap)**s), 0.02, 1.5)[0] for c in (lam, mu))
    OUT["E1"] = dict(params=[str(v) for v in (lam, mu, nu, kap)], conditions=cond, v3=vals, shared_primes=odd_shared_primes(mu, kap),
                     collision=[str(v) for v in x], collision_equal=(x == y), disc_version=disc,
                     Dstar=[Dlo, Dhi], exact_dim=[Flo, Fhi], drop_at_least=Dlo - Fhi, formulas_agree=agree,
                     old_style_lower_bound=lb, rho_min=str(2 - kap), rho_threshold=str(1 / (1 - kap)))
    print("E1 params", OUT["E1"]["params"], cond, vals, "collision", x, x == y)
    print("   disc version", disc)
    print("   D* in", (Dlo, Dhi), " EXACT dim in", (Flo, Fhi), " drop >=", Dlo - Fhi, " old lower bound", lb)

    # E2: sanity only
    ok, total, bad = reduced_words_distinct(lam, mu, nu, kap, 8)
    EX3 = [((F(0), F(0)), [(lam, F(0)), (mu, F(0))]), ((F(1), F(0)), [(nu, F(0))]), ((F(2), F(0)), [(kap, F(0))])]
    dc = G.distinct_counts(EX3, 9)
    model = Bk.__dict__.get("_model_cache")
    import sympy as sp
    t = sp.symbols("t"); P1t = 1 / (1 - t)**2; Qt = 1 / (1 - 2 * t); St = (P1t - 1) * (Qt - 1)
    model = [int(v) for v in sp.Poly(sp.series(Qt * P1t / (1 - (St - t**2 * P1t * Qt)), t, 0, 10).removeO(), t).all_coeffs()[::-1][1:]]
    OUT["E2_sanity"] = dict(reduced_words_upto8=total, all_distinct=ok, counterexample=bad, distinct_counts=dc,
                            family_model=model, match=(dc == model))
    print("E2 sanity: reduced words <= 8:", total, "all distinct:", ok, "| counts", dc, "model match", dc == model)

    # E3: where Theorem P12 applies (floats; the figure)
    from scipy.optimize import brentq
    grid = []
    for m_ in np.linspace(0.004, 0.30, 75):
        for k_ in np.linspace(0.004, 0.50, 75):
            l_ = m_ * (2 - k_); n_ = k_ / (2 - k_)
            if not (l_ < 1 and n_ < 1): continue
            I1 = 2 * l_ < 1 - n_; I2 = 1 + n_ < 2 - 2 * k_; II = (2 - k_) * (1 - k_) > 1
            ok_ = I1 and I2 and II
            sF = D = None
            if ok_:
                sq = brentq(lambda s: 1 - n_**s - k_**s, 1e-9, 5)
                def wF(s):
                    Lm = 1 / (1 - l_**s); Mu = 1 / (1 - m_**s); Q = 1 / (1 - n_**s - k_**s)
                    return (Lm - 1) * (Q - 1) + Lm * (Mu - 1) * n_**s * Q - 1
                sF = brentq(wF, sq + 1e-10, 5)
                D = brentq(lambda s: (1 - l_**s) * (1 - m_**s) - n_**s - k_**s, 1e-6, 5)
            grid.append([float(m_), float(k_), bool(I1), bool(I2), bool(II), sF, D])
    OUT["E3_region"] = grid
    okpts = [g for g in grid if g[5] is not None]
    print("E3 grid", len(grid), "points where P12's (I),(II) hold:", len(okpts),
          "drop range", min(g[6] - g[5] for g in okpts), max(g[6] - g[5] for g in okpts),
          "max s_F", max(g[5] for g in okpts))

    # E4: exact answers converging to the original break
    seq = []
    for p in [3, 5, 7, 11, 13, 101, 1009]:
        m_, k_ = F(p, 20 * p + 1), F(p, 10 * p + 1)
        L4 = surface(m_, k_); c4 = conditions(*L4)
        assert c4["on_surface"] and c4["I_class1_vs_N"] and c4["I_N_vs_K"] and c4["II"] and p in odd_shared_primes(m_, k_)
        _, (a_, b_) = certified_roots(*L4)
        seq.append(dict(p=p, mu=str(m_), kappa=str(k_), exact_dim=[a_, b_]))
    OUT["E4_sequence"] = seq
    print("E4", [(s_["p"], round(s_["exact_dim"][0], 9)) for s_ in seq])

    # E5: the original break
    lam0, mu0, nu0, kap0 = Bk.LAM, Bk.MU, Bk.NU, Bk.KAP
    c0 = conditions(lam0, mu0, nu0, kap0)
    badN = ((1 - nu0) / (1 + nu0), (1 + nu0) / (1 - nu0))          # rho J_N meets J_N
    badK = ((1 - nu0) / 2, (1 + nu0) / (2 - 2 * kap0))              # rho J_K meets J_N
    def bad(r): return badN[0] <= r <= badN[1] or badK[0] <= r <= badK[1]
    first_bad = None
    for a in range(1, 200):
        for b in range(1, 200):
            r = lam0**a / mu0**b
            if bad(r):
                first_bad = (a, b, float(r)); break
        if first_bad: break
    kmax = first_bad[0] - 1
    # bad pairs with a <= kmax?  none by construction; double-check with exact fractions
    assert not any(bad(lam0**a / mu0**b) for a in range(1, kmax + 1) for b in range(1, 200))
    (D0lo, D0hi), (F0lo, F0hi) = certified_roots(lam0, mu0, nu0, kap0)
    _, (L0lo, L0hi) = certified_roots(lam0, mu0, nu0, kap0, k=kmax)
    OUT["E5_original"] = dict(conditions=c0, shared_primes=odd_shared_primes(mu0, kap0), bad_N=[str(v) for v in badN],
                              bad_K=[str(v) for v in badK], first_bad=first_bad, kmax=kmax,
                              upper=[F0lo, F0hi], new_lower=[L0lo, L0hi], gap_at_most=F0hi - L0lo,
                              round9_lower=0.4367041, Dstar=[D0lo, D0hi])
    print("E5 original: shared odd primes", odd_shared_primes(mu0, kap0), "bad windows", [float(v) for v in badN], [float(v) for v in badK])
    print("   first bad (a, b, rho)", first_bad, "-> k =", kmax, " new lower bound", (L0lo, L0hi), " upper", (F0lo, F0hi),
          " gap <=", F0hi - L0lo)
    json.dump(OUT, open("exact.json", "w"), indent=1, default=str)
