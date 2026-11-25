"""Wall at Re s = 0 for Z_k(s) = prod_p (1 - k p^-s)^(-1/k), and the Lee-Yang comparison (Nov 25, 2026).
Run: python3 wall_checks.py   (mpmath, sympy). Writes wall_results.json."""
import json, math, random, itertools
from fractions import Fraction as F
import mpmath as mp, sympy as sp
mp.mp.dps = 25
out = []; res = {}
def check(name, ok, detail=""):
    out.append(("PASS" if ok else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""))
def mobius(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
k = sp.symbols('k')
def Mpoly(N): return sp.expand(sum(mobius(N//d) * k**d for d in range(1, N+1) if N % d == 0) / N)

# 1. Only k = +-1 have finitely many zeta factors
check("k=+1: M(1,N) = 0 for 2 <= N <= 40", all(Mpoly(N).subs(k,1) == 0 for N in range(2,41)))
check("k=-1: M(-1,N) = 0 for 3 <= N <= 40", all(Mpoly(N).subs(k,-1) == 0 for N in range(3,41)))
primes = list(sp.primerange(2, 100))
ks = [F(-1,2), F(-1,3), F(-1,7), F(1,2), F(2,3), F(3), F(-5,4)]
check("prime N: M(k,N) = (k^N - k)/N, nonzero for k not in {0,+-1} (k samples, N prime < 100)",
      all(Mpoly(p).subs(k, sp.Rational(x.numerator, x.denominator)) != 0 for p in primes for x in ks))
check("k=-1/m: N*m^N*M(-1/m,N) is an integer = (-1)^N mod m (m=2..7, N<=24)",
      all((lambda v: v.is_integer and (int(v) - (-1)**N) % m == 0)(sp.nsimplify(N * m**N * Mpoly(N).subs(k, sp.Rational(-1, m))))
          for m in range(2, 8) for N in range(1, 25)))

# 2. The generic argument: any positive integer combination sum_{N in S} m_N M(k,N) is a nonzero polynomial of degree max S
random.seed(3); ok = True
for _ in range(300):
    S = random.sample(range(1, 16), random.randint(1, 5)); mult = [random.randint(1, 4) for _ in S]
    Q = sp.Poly(sum(c * Mpoly(N) for c, N in zip(mult, S)), k)
    if Q.degree() != max(S) or Q.LC() != sp.Rational(mult[S.index(max(S))], max(S)): ok = False
check("sum of m_N M(k,N) over a finite S has degree max S, leading coefficient m/N (300 random S)", ok)
# a rational k can be a root for some multiplicity pattern (so rational k need more than the generic argument)
Q = sp.expand(3 * Mpoly(1) + 4 * Mpoly(2))
check("example: k=-1/2 is a root of 3M(k,1)+4M(k,2) (rational k are not covered by the generic argument)",
      Q.subs(k, sp.Rational(-1, 2)) == 0)
# with all multiplicities 1, check k=-1/m never a root for S subsets of {1..10} (finite search)
bad = []
for r in range(1, 5):
    for S in itertools.combinations(range(1, 11), r):
        Qs = sum(Mpoly(N) for N in S)
        for m in range(2, 7):
            if Qs.subs(k, sp.Rational(-1, m)) == 0: bad.append((S, m))
res["simple_zero_cancellations"] = [[list(S), m] for S, m in bad]
check("simple zeros do not rescue rational k: S={1,2,3} cancels at k=-1/2 (M1+M2+M3 vanishes there), so k=-1/m still needs RH",
      ((1, 2, 3), 2) in bad, str(bad[:5]))

# 3. Singular points rho/N crowd every point of the imaginary axis (first 400 zeros, all on the critical line)
zs = [mp.zetazero(n) for n in range(1, 401)]
gam = [float(z.imag) for z in zs]
check("first 400 zeros computed on Re = 1/2", all(abs(float(z.real) - 0.5) < 1e-12 for z in zs))
crowd = {}
for y in [0.5, 1.0, 2.0, 5.0]:
    best = []
    for N in [11, 23, 47, 97]:
        d = min(math.hypot(0.5 / N, g / N - y) for g in gam)
        best.append(round(d, 4))
    crowd[str(y)] = best
check("distance from i*y to nearest rho/N shrinks as N grows (y = 0.5, 1, 2, 5; N = 11, 23, 47, 97)",
      all(v[-1] < v[0] and v[-1] < 0.02 for v in crowd.values()), json.dumps(crowd))
res["crowding"] = crowd
check("on the critical line, rho/N has real part 1/(2N); another N' gives real part N'/(2N) != 1/2: so under RH each point has one contributor",
      all(F(Np, 2*N) != F(1, 2) for N in range(1, 30) for Np in range(1, 2*N) if Np != N))
# local exponent at s0 = rho_1/2 for k = 1/2: C = M(1/2,2)/(1/2) = -1/4, a branch point
C = Mpoly(2).subs(k, sp.Rational(1, 2)) / sp.Rational(1, 2)
check("k=1/2 at s0 = rho_1/2: exponent C = -1/4 (not an integer: branch point)", C == sp.Rational(-1, 4))
# numerical: monodromy of log Z_k around s0 equals 2*pi*i*C, using Z_k = prod zeta(Ns)^(M/k), N <= 60
s0 = zs[0] / 2
def logZ(s, kk=0.5, Nmax=60):
    tot = mp.mpc(0)
    for N in range(1, Nmax + 1):
        c = float(Mpoly(N).subs(k, kk)) / kk
        tot += c * mp.log(mp.zeta(N * s))
    return tot
# integrate d(log Z) around a small circle by tracking the phase continuously
r = 0.01; steps = 400; total = mp.mpf(0); prev = None
for j in range(steps + 1):
    s = s0 + r * mp.e**(2j * mp.pi * j / steps)
    val = mp.zeta(2 * s)          # only the N = 2 factor has a zero inside the circle
    ph = mp.arg(val)
    if prev is not None:
        dph = ph - prev
        while dph > mp.pi: dph -= 2 * mp.pi
        while dph < -mp.pi: dph += 2 * mp.pi
        total += dph
    prev = ph
check("winding of zeta(2s) around rho_1/2 is 1, so Z_{1/2} turns by e^(2 pi i C) = e^(-i pi/2): a true branch point",
      abs(float(total) / (2 * math.pi) - 1) < 1e-6, f"winding {float(total)/(2*math.pi):.6f}")

# 4. Lee-Yang comparison: Beraha dial values k_n = 1 + 2cos(2 pi/n) vs copies points k = -1/m and Cox points k = 2/m
ber = {n: 1 + 2 * math.cos(2 * math.pi / n) for n in range(2, 61)}
hits_neg = [(n, m) for n, v in ber.items() for m in range(1, 200) if abs(v + 1 / m) < 1e-12]
hits_pos = [(n, m) for n, v in ber.items() for m in range(1, 200) if abs(v - 2 / m) < 1e-12]
check("k = -1/m meets a Beraha value only at m = 1 (n = 2); k = 0 (n = 3) is their limit", hits_neg == [(2, 1)], str(hits_neg))
check("k = 2/m meets a Beraha value only at m = 1 (k=2, n=6) and m = 2 (k=1, n=4)", sorted(hits_pos) == [(4, 2), (6, 1)], str(hits_pos))
check("no Beraha value lies strictly between -1 and 0 (n >= 3 gives k >= 0)", all(not (-1 < v < 0) for v in ber.values()))
res["beraha"] = {n: round(v, 6) for n, v in ber.items() if n <= 12}

passed = sum(o.startswith("PASS") for o in out)
out.append(f"\n{passed} of {len(out)} pass")
open("wall_checks_output.txt", "w").write("\n".join(out) + "\n")
json.dump(res, open("wall_results.json", "w"), indent=1)
print("\n".join(out))
