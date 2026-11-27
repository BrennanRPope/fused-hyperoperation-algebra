"""checks20.py: round 20 (DRAFT) of topics/new-ties: the 19-adic step for five letters, and the referee's recheck.

Q1  the five-letter casework: 12 of 3,838 tie patterns cancel modulo 19, all of shape N x1 x2 x3 m^M N with tied
    terms 0, 1, 4, a bare last N, and front L's adding to A3 (Proposition 20.1).
Q2  the identity T = 19^A3 (1 + 20^M u), u = T3 / 19^(A3 - 1) = -1 mod 19, and g = A3 + 1 + v19(M - mu(u)), with T and
    mu computed independently of adic20.py, on random words of every pattern.
Q3  the digit argument (max over 0 <= M <= B of v19(M - mu) is the largest k with mu mod 19^k <= B) against brute force,
    for every pattern with free front gaps <= 2 and M <= 3000.
Q4  the family N K l K K m^M N: u = 220,209 = -1 mod 361; the first M with g = 4, ..., 8 is mu mod 19^k (38, 5,092,
    25,669, 1,459,200, 43,552,883), and the exact maxima of g over M <= 10^100, 10^200, 10^300 are 80, 159, 236.
Q5  Proposition 20.3: the exact maxima of g over M <= 10^2, 10^100, 10^200, 10^300 are 9, 87, 164, 242 (adic20.json, over
    all 13,517 fronts with free front gaps <= 30): recomputed for every front with free front gaps <= 12 and for the
    fronts where each maximum is reached (the full recomputation is `python3 adic20.py all`, about 4 minutes).
Q6  the two-logarithm 19-adic bound (Bugeaud and Laurent 1996, as restated; Proposition 20.2): constant 0.337, and for the
    family g <= 10,771 for M <= 10^13 and <= 652,697 for M <= 10^100; every u is a 19-unit with a prime factor other than
    2 and 5 (so 20 and -u are multiplicatively independent), and 20, -u are both 1 mod 19 (so the bound's g is 1).
Q7  the size step survives: with g <= A3 + 1 + log_19 B + 4, the factor 19^g is at most 19^9 * B, so it changes the
    other word's gap-sum bound by at most log(19^9 B)/log 20 (for B = 10^100, at most 86), a fixed amount, not a factor.
Writes checks20_output.txt and checks20.json.
"""
import json, math, random, sys, time, tracemalloc
import adic20 as A

tracemalloc.start()
t_start = time.time()
LOG19, LOG20 = math.log(19), math.log(20)
results = []


def record(name, ok, detail):
    results.append(dict(name=name, ok=bool(ok), detail=detail))
    print(("PASS " if ok else "FAIL ") + name + ": " + detail, flush=True)


# --- independent arithmetic (not from adic20) ---
def T_of(items):
    T = 1; s = 0
    for (r, L, M) in items:
        s += L; T = 19 ** s + r * 200 ** L * 20 ** M * T
    return T


def v19(n):
    n = abs(n); k = 0
    while n % 19 == 0: n //= 19; k += 1
    return k


def mu_mod(u, k):
    """mu mod 19^k with 20^mu = -1/u in Z_19, by digit search (slow, independent of adic20.mu_digits)"""
    mod = 19 ** (k + 1); target = (-pow(u, -1, mod)) % mod
    mu = 0
    for j in range(k):
        q = 19 ** (j + 2)
        for d in range(19):
            c = mu + d * 19 ** j
            if (pow(20, c, q) - target) % q == 0: mu = c; break
        else:
            raise AssertionError("no digit")
    return mu


ties, bad = A.casework(4, 5)

# Q1
try:
    ok = ties == 3838 and len(bad) == 12
    for p in bad:
        ok &= p["tied"] == [0, 1, 4] and p["rs"][3] == 19 and p["Ls"][3] == 0 and sum(p["Ls"]) == p["tmin"]
    # four letters (three items) has no cancelling pattern (F3, round 19)
    t3, b3 = A.casework(3, 4)
    ok &= len(b3) == 0
    record("Q1 five-letter casework: cancelling patterns all of one shape", ok,
           "%d tie patterns, %d cancel, all with tied terms 0, 1, 4 and a bare last N; four letters: %d ties, %d cancel"
           % (ties, len(bad), t3, len(b3)))
except Exception as e:
    record("Q1 five-letter casework", False, repr(e))

# Q2
try:
    rnd = random.Random(2020); n = 0
    for p in bad:
        for _ in range(25):
            Ms = tuple(rnd.randint(0, 9) for _ in range(A.n_free(p)))
            front = A.pattern_front(p, Ms)
            A3 = sum(L for (_, L, _) in front)
            T3 = T_of(front)
            assert v19(T3) == A3 - 1
            u = T3 // 19 ** (A3 - 1)
            assert u % 19 == 18
            mu = mu_mod(u, 8)
            for M in (rnd.randint(0, 400), rnd.randint(0, 5000)):
                W = front + ((19, 0, M),)
                T = T_of(W)
                assert T == 19 ** A3 * (1 + 20 ** M * u)
                vm = v19(M - mu)
                if vm < 8:
                    assert v19(T) == A3 + 1 + vm
                n += 1
    record("Q2 T = 19^A3 (1 + 20^M u) and g = A3 + 1 + v19(M - mu(u))", True,
           "%d random words of the 12 patterns, T and mu computed independently, 0 exceptions" % n)
except Exception as e:
    record("Q2 identity", False, repr(e))

# Q3
try:
    nf = 0; worst = 0
    for p in bad:
        for Ms, front in A.fronts(p, 2):
            A3, u = A.u_of(front)
            mod = 19 ** (A3 + 14); base = 19 ** A3; pw = 1; best = -1
            for M in range(3001):
                T = (base + base * pw * u) % mod
                g = v19(T) if T else 10 ** 9
                assert g < A3 + 14
                best = max(best, g); pw = pw * 20 % mod
            mu = A.mu_digits(u, 14)
            k, Mk = A.max_v_upto(mu, 14, 3000)
            assert best == A3 + 1 + k, (front, best, A3 + 1 + k)
            nf += 1; worst = max(worst, best)
    record("Q3 the digit argument equals brute force", True,
           "%d fronts (free front gaps <= 2), every M <= 3000: exact max g agrees; largest %d" % (nf, worst))
except Exception as e:
    record("Q3 digit argument", False, repr(e))

# Q4
fam_out = {}
try:
    fam = next(p for p in bad if p["rs"] == (10, 10, 10, 19))
    front = A.pattern_front(fam, ())
    T3 = T_of(front)
    assert T3 == 220209 and T3 % 361 == 360
    mu = A.mu_digits(T3, A.K_for(10 ** 300))
    assert mu % 19 ** 8 == mu_mod(T3, 8)
    firsts = [mu % 19 ** k for k in range(2, 7)]
    assert firsts == [38, 5092, 25669, 1459200, 43552883], firsts
    for k, M in zip(range(2, 7), firsts):
        assert v19(T_of(front + ((19, 0, M),))) == 2 + k
    maxima = {e: 2 + A.max_v_upto(mu, A.K_for(10 ** e), 10 ** e)[0] for e in (100, 200, 300)}
    assert maxima == {100: 80, 200: 159, 300: 236}, maxima
    fam_out = dict(u=T3, first_M=firsts, maxima=maxima)
    record("Q4 the family N K l K K m^M N", True,
           "u = 220209 = -1 mod 361; first M with g = 4..8: %s (checked on T itself); max g for M <= 10^100, 10^200, "
           "10^300: %s" % (firsts, maxima))
except Exception as e:
    record("Q4 family", False, repr(e))

# Q5
q5 = {}
try:
    Bs = {2: 10 ** 2, 100: 10 ** 100, 200: 10 ** 200, 300: 10 ** 300}
    Kmax = A.K_for(10 ** 300)
    mx = {e: 0 for e in Bs}; nf = 0; excess = 0
    jr = json.load(open("adic20.json"))["reduction"]["by_B"]
    argmax = {jr[str(e)]["at"]["word"].split(" ")[0] for e in Bs}
    for p in bad:
        for Ms, front in A.fronts(p, 30):
            if sum(Ms) > 12 and A.name(front) not in argmax:
                continue
            A3, u = A.u_of(front)
            mu = A.mu_digits(u, Kmax)
            nf += 1
            for e, B in Bs.items():
                k = A.max_v_upto(mu, A.K_for(B), B)[0]
                mx[e] = max(mx[e], A3 + 1 + k)
                excess = max(excess, k - math.ceil(math.log(B) / LOG19))
    want = {e: json.load(open("adic20.json"))["reduction"]["by_B"][str(e)]["max_g"] for e in Bs}
    ok = mx == want == {2: 9, 100: 87, 200: 164, 300: 242}
    q5 = dict(fronts=nf, max_g=mx, largest_excess_over_ceil_log19B=excess)
    record("Q5 Proposition 20.3: the maxima are reached, and no front with free front gaps <= 12 beats them", ok,
           "%d fronts (all with free front gaps <= 12, plus the fronts where adic20.json says each maximum is reached; "
           "the full run over all 13,517 fronts with gaps <= 30 is adic20.py all, about 4 min); max g for M <= 10^2, 10^100, 10^200, 10^300: %s (adic20.json: %s); largest k - ceil(log_19 B) = %d"
           % (nf, mx, want, excess))
except Exception as e:
    record("Q5 Proposition 20.3", False, repr(e))

# Q6
try:
    C = 24 * 19 * 1 / ((19 - 1) * LOG19 ** 4)
    assert abs(C - 0.337) < 0.0005
    b13 = A.bl_bound_g(1, 220209, 10 ** 13)
    b100 = A.bl_bound_g(1, 220209, 10 ** 100)
    assert 10769 < b13 < 10772 and 652695 < b100 < 652698, (b13, b100)
    # independent recomputation of the same formula
    h2 = math.log(220209); Ep = 10 ** 100 / h2 + 1 / LOG20
    E = max(math.log(Ep) + math.log(LOG19) + 0.4, 10, 10 * LOG19)
    assert abs(1 + C * E * E * LOG20 * h2 - b100) < 1
    nu = 0
    for p in bad:
        for Ms, front in A.fronts(p, 6):
            A3, u = A.u_of(front)          # asserts u = -1 mod 19 and a prime factor other than 2, 5
            assert u % 19 == 18 and 20 % 19 == 1
            nu += 1
    record("Q6 the two-logarithm 19-adic bound, applied", True,
           "constant 24*19/(18 (log 19)^4) = %.4f; family: g <= %.0f for M <= 10^13, <= %.0f for M <= 10^100 "
           "(x 18 under the conservative reading); %d fronts: u is a 19-unit, -1 mod 19, with a prime factor other "
           "than 2 and 5" % (C, b13, b100, nu))
except Exception as e:
    record("Q6 Bugeaud-Laurent", False, repr(e))

# Q7
try:
    loss = {}
    for e in (100, 200, 300):
        B = 10 ** e
        gmax = 4 + 1 + math.log(B) / LOG19 + 4          # A3 <= 4, plus the largest excess seen (Q5), rounded up
        loss[e] = math.ceil(gmax * LOG19 / LOG20)
    ok = loss[100] <= 86 and q5.get("largest_excess_over_ceil_log19B", 99) <= 4
    record("Q7 the size step survives a logarithmic g", ok,
           "with g <= A3 + 1 + log_19 B + 4 (A3 <= 4), 19^g changes the other word's gap-sum bound by at most %s for "
           "B = 10^100, 10^200, 10^300: an additive amount, not a factor" % loss)
except Exception as e:
    record("Q7 size step", False, repr(e))

npass = sum(r["ok"] for r in results)
cur, peak = tracemalloc.get_traced_memory()
summary = "%d/%d checks pass" % (npass, len(results))
lines = [("PASS " if r["ok"] else "FAIL ") + r["name"] + ": " + r["detail"] for r in results]
lines += [summary]
open("checks20_output.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(results=results, family=fam_out, q5=q5), open("checks20.json", "w"), indent=1, default=str)
print(summary)
print("TIME %.0f s, PEAK %.0f MB" % (time.time() - t_start, peak / 2 ** 20))
sys.exit(0 if npass == len(results) else 1)
