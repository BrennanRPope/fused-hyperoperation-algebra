"""Round 6 (DRAFT): no accidental collisions for the original examples (README section 10).

The original examples: law 1 has walls (attracting 1, repelling Q) and multipliers 1/5, 3/20 (letters a, b);
law 2 has walls (attracting -1, repelling -Q) and multiplier 1/4 (letter c); Q = 2, 3, 12.
A law-i map acts as z_i -> lam * z_i with z_1 = (y-1)/(y-Q), z_2 = (y+1)/(y+Q).

What this script checks (every certified claim uses exact fractions; samples are only illustrations):
  R1  real ping-pong: law-1 runs map [-1,1] into J1 near 1, law-2 runs into J2 near -1; J1, J2 disjoint, 0 in neither.
  R2  the real ratio window for the last step, and that (3/4)^n, 5 (3/4)^n, (3/4)^n / 5 all miss it (walls +-3).
  T1  2-adic constants for even Q (Theorem P5): v2(kappa) >= 3 and |1 - kappa|_2 = 1, kappa = 4Q/(Q+1)^2.
  T2  illustration: the 2-adic size of z_1 at the end of random words that end in law 2, started at -Q.
      One value for Q = 2 and 12 (the proof says so); several for Q = 3 (why the 2-adic proof fails there).
  T3  5-adic constants for walls +-3 (Theorem P6): the four walls are distinct mod 5, and the residues used.
  T4  illustration of P6: 5-adic size of z_1 for random words in the allowed class is 1 or 5; outside it, more.
  L1  certified lower bound for the dimension at walls +-3 from the provably free sub-system of P6
      (python-flint ball arithmetic + Collatz-Wielandt vector, as in exact_bounds.py)."""
import random, json, math
from fractions import Fraction as F
import exact_bounds as EB

def v(x, p):
    x = F(x)
    if x == 0: return math.inf
    n, d, k = x.numerator, x.denominator, 0
    while n % p == 0: n //= p; k += 1
    while d % p == 0: d //= p; k -= 1
    return k

def law_apply(y, p, q, lam):
    """apply z -> lam z in z = (y-p)/(y-q); points are Fractions or None (= infinity)"""
    if y == q: return q
    z = F(1) if y is None else (y - p) / (y - q)
    z2 = lam * z
    return None if z2 == 1 else (p - q * z2) / (1 - z2)

def word_apply(word, y, Q):
    """word: list of syllables, innermost first; ('A', alpha, beta) or ('C', k)"""
    for s in word:
        if s[0] == 'A': y = law_apply(y, F(1), F(Q), F(1, 5) ** s[1] * F(3, 20) ** s[2])
        else:           y = law_apply(y, F(-1), F(-Q), F(1, 4) ** s[1])
    return y

def z1(y, Q): return F(1) if y is None else (y - 1) / (y - Q)

# ---------- R1, R2: real ping-pong ----------
def real_part(Q):
    Q = F(Q)
    zmax1 = F(2, 5) / (Q + 1)          # z_1(-1) = 2/(Q+1), largest law-1 multiplier 1/5
    zmax2 = F(2, 4) / (Q + 1)          # z_2(1) = 2/(Q+1), law-2 multiplier 1/4
    J1 = (EB.yc(zmax1, F(1), Q), F(1)); J2 = (F(-1), EB.yc(zmax2, F(-1), -Q))
    ok = J1[0] > 0 and J2[1] < 0
    zs = sorted([z1(J2[0], Q), z1(J2[1], Q)])            # z_1 is monotone on J2
    window = (zs[0] / zs[1], zs[1] / zs[0])
    return dict(J1=[str(x) for x in J1], J2=[str(x) for x in J2], disjoint_and_miss_0=ok,
                window=[str(x) for x in window], window_float=[float(x) for x in window])

def window_misses(window, n_max=80):
    lo, hi = F(window[0]), F(window[1]); bad = []
    for n in range(1, n_max):
        for name, r in (("(3/4)^n", F(3, 4) ** n), ("5(3/4)^n", 5 * F(3, 4) ** n), ("(3/4)^n/5", F(3, 4) ** n / 5)):
            if lo <= r <= hi: bad.append((name, n))
    # beyond n_max all three are < (3/4)^80 * 5 < lo, and they decrease, so the finite check covers every n
    tail_ok = 5 * F(3, 4) ** n_max < lo
    return bad, tail_ok

# ---------- T1: 2-adic constants for even Q ----------
def two_adic_constants(Q):
    Q = F(Q); kappa = 4 * Q / (Q + 1) ** 2
    return dict(kappa=str(kappa), v2_kappa=v(kappa, 2), v2_one_minus_kappa=v(1 - kappa, 2),
                M_exp=v((Q + 1) / (2 * Q), 2), s_star_exp=v((Q + 1) / (2 * Q), 2),
                criterion=v(kappa, 2) >= 1 and v(1 - kappa, 2) == 0)

# ---------- random words ----------
def rand_word(rng, nsyl, law2_rule=None, max_e=4, max_k=6):
    """random word ending (outermost) in law 2; law2_rule(k, prev) decides whether a law-2 run of length k
    may sit directly on the previous law-1 run prev = (alpha, beta) or None"""
    word = []; law = rng.choice('AC'); prev = None
    while len(word) < nsyl or law != 'C' or not word:
        if law == 'A':
            while True:
                a, b = rng.randint(0, max_e), rng.randint(0, max_e)
                if a + b: break
            word.append(('A', a, b)); prev = (a, b); law = 'C'
        else:
            while True:
                k = rng.randint(1, max_k)
                if law2_rule is None or law2_rule(k, prev): break
            word.append(('C', k)); prev = None; law = 'A'
            if len(word) >= nsyl: break
    if word[-1][0] != 'C': word.append(('C', 1))
    return word

def p6_rule(k, prev):
    if k % 2 == 1: return True
    if k % 10 == 0: return False
    return prev is None or prev[0] + prev[1] >= 2

def bad_rule(k, prev):  # forces the excluded pattern: even law-2 run on a one-letter law-1 run
    return prev is None or (k % 2 == 0 and prev[0] + prev[1] == 1) or k % 2 == 1

def sample(Q, p, base, rule, n=400, seed=1):
    rng = random.Random(seed); vals = {}
    for _ in range(n):
        w = rand_word(rng, rng.randint(1, 7), rule)
        val = v(z1(word_apply(w, base, Q), Q), p)
        vals[val] = vals.get(val, 0) + 1
    return dict(sorted(vals.items()))

import itertools
def all_words(nsyl, max_e=2, max_k=4, rule=None):
    """every word that ends (outermost) in law 2, with nsyl runs, run exponents up to max_e / max_k"""
    A_runs = [('A', a, b) for a in range(max_e + 1) for b in range(max_e + 1) if a + b]
    C_runs = [('C', k) for k in range(1, max_k + 1)]
    for start in 'AC':
        laws = [('A' if (i % 2 == 0) == (start == 'A') else 'C') for i in range(nsyl)]
        if laws[-1] != 'C': continue
        for combo in itertools.product(*[A_runs if l == 'A' else C_runs for l in laws]):
            if rule:
                ok = True
                for i, s in enumerate(combo):
                    if s[0] == 'C':
                        prev = combo[i - 1] if i > 0 else None
                        if not rule(s[1], None if prev is None else (prev[1], prev[2])): ok = False
                if not ok: continue
            yield list(combo)

# ---------- T3: 5-adic constants at walls +-3 ----------
def five_adic_constants(Q=3):
    walls = [1, Q, -1, -Q]; res = sorted(x % 5 for x in walls)
    z2Q = F(Q + 1, 2 * Q); z21 = F(2, Q + 1)
    r = lambda x: (x.numerator * pow(x.denominator, -1, 5)) % 5
    return dict(walls_mod5=res, distinct=len(set(res)) == 4,
                z2_Q_mod5=r(z2Q), z2_1_mod5=r(z21), minus_z2_Q_mod5=(-r(z2Q)) % 5,
                odd_run_lands_on=('infinity' if (-r(z2Q)) % 5 == 1 else 'finite'),
                odd_run_avoids_1_and_Q=((-r(z2Q)) % 5 != r(z21)) and ((-r(z2Q)) % 5 != r(z2Q)))

# ---------- L1: certified lower bound for the free sub-system at walls +-3 ----------
class SubSetup(EB.Setup):
    """Same switching matrix as exact_bounds.Setup, but law-2 runs of length k are dropped when 10 | k, or when
    k is even and the previous law-1 run has one letter (bucket 0). Only the 'low' matrix is used."""
    def T(self, t, which):
        assert which == 'low'
        L, nb = self.L, self.nb; n = L * nb; M = [[EB.arb(0)] * n for _ in range(n)]
        tt = EB.A(t) if isinstance(t, F) else EB.arb(t)
        for (j, k, i), (za, zb) in self.zr.items():
            z = zb; accb = [EB.arb(0)] * nb
            for deg, lam in self.blk[i]:
                if i == 1 and (deg % 10 == 0 or (deg % 2 == 0 and j == 0 and k == 0)): continue
                base = lam * ((1 - z) / (1 - lam * z)) ** 2
                accb[self.bucket(deg)] += EB.A(base) ** tt
            for b2 in range(nb): M[i * nb + b2][j * nb + k] = accb[b2]
        return M

def sub_lower(Q=3, nb=6, digits=4):
    from scipy.optimize import brentq
    cl = [(F(1), F(Q), [F(1, 5), F(3, 20)]), (F(-1), F(-Q), [F(1, 4)])]
    S = SubSetup(cl, nb=nb)
    lo = brentq(lambda t: S.sr(t, 'low') - 1, 1e-3, 3)
    tl = F(math.floor(lo * 10 ** digits), 10 ** digits) - F(1, 10 ** digits)
    ok, margin = S.certify(tl, 'low')
    return dict(t_low=str(tl), float=lo, certified=ok, margin=margin)

if __name__ == "__main__":
    out = {}
    for Q in (2, 3, 12):
        rp = real_part(Q); out[f"real_{Q}"] = rp
        print(f"R1 walls ±{Q}: J1={rp['J1']} J2={rp['J2']} disjoint & miss 0: {rp['disjoint_and_miss_0']}; "
              f"window {rp['window_float']}")
    bad, tail = window_misses(out["real_3"]["window"]); out["R2"] = dict(bad=bad, tail_ok=tail)
    print("R2 walls ±3: forms inside the window:", bad, "| tail covered:", tail)
    for Q in (2, 12, 3):
        c = two_adic_constants(Q); out[f"T1_{Q}"] = c; print(f"T1 walls ±{Q}:", c)
    for Q in (2, 12, 3):
        s = sample(Q, 2, F(-Q), None); out[f"T2_{Q}"] = {str(k): x for k, x in s.items()}
        print(f"T2 walls ±{Q}: 2-adic valuations of z1 at the end of 400 random law-2-ending words: {s}")
    c5 = five_adic_constants(3); out["T3"] = c5; print("T3 walls ±3:", c5)
    s_in = sample(3, 5, F(-3), p6_rule); s_out = sample(3, 5, F(-3), bad_rule, seed=2)
    out["T4"] = dict(allowed={str(k): x for k, x in s_in.items()}, excluded_pattern={str(k): x for k, x in s_out.items()})
    print("T4 walls ±3: 5-adic valuations of z1, allowed class:", s_in, "| with the excluded pattern:", s_out)
    lb = sub_lower(3); out["L1"] = lb; print("L1 walls ±3: certified lower bound for the free sub-system:", lb)
    json.dump(out, open("nocollide.json", "w"), indent=1, default=str)
