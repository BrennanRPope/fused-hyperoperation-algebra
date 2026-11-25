"""Checks for the power set with one dial per element (Nov 25, 2026).
Pure Python (fractions, math, random). Run: python3 checks.py"""
from fractions import Fraction as F
from itertools import product
import math, random

results = []
def check(name, ok, detail=""):
    results.append(ok); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))

def subsets(n): return list(product([0,1], repeat=n))
def W(ts, a):
    w = 1
    for t, b in zip(ts, a): w *= t if b else 1
    return w
def multadd(a, b): return a + b + a*b
def chain(ts):
    x = 0
    for t in ts: x = multadd(x, t)
    return x

# 1. count identity, exact
random.seed(1)
for n in range(0, 7):
    ts = [F(random.randint(0, 9), random.randint(1, 5)) for _ in range(n)]
    tot = sum(W(ts, a) for a in subsets(n))
    prod = 1
    for t in ts: prod *= 1 + t
    check(f"count n={n}: sum of weights = prod(1+t) = 1 + chain", tot == prod == 1 + chain(ts))
check("all dials 1 -> 2^n, all 0 -> 1, all 1/2 -> 1.5^n",
      sum(W([1]*5, a) for a in subsets(5)) == 32 and sum(W([0]*5, a) for a in subsets(5)) == 1
      and sum(W([F(1,2)]*5, a) for a in subsets(5)) == F(3,2)**5)

# 2. merging = union convolution, brute force
def conv(F1, F2, n):
    out = {c: 0 for c in subsets(n)}
    for a in subsets(n):
        for b in subsets(n):
            c = tuple(x | y for x, y in zip(a, b)); out[c] += F1(a) * F2(b)
    return out
for n in range(1, 5):
    ts = [F(random.randint(0, 6), random.randint(1, 4)) for _ in range(n)]
    ss = [F(random.randint(0, 6), random.randint(1, 4)) for _ in range(n)]
    c = conv(lambda a: W(ts, a), lambda b: W(ss, b), n)
    ms = [multadd(t, s) for t, s in zip(ts, ss)]
    check(f"merge n={n}: W_t * W_s = W_(t ∘₁ s)", all(c[x] == W(ms, x) for x in c))
# coarse-graining: block of 2 elements present if either present
t, s = F(2, 3), F(5, 7)
inblock = t + s + t*s
check("block of two elements carries dial t ∘₁ s", inblock == multadd(t, s))

# 3. the half
h = math.sqrt(2) - 1
c = conv(lambda a: W([h]*4, a), lambda b: W([h]*4, b), 4)
check("half: W_h * W_h = P(S) (every subset weight 1), n=4", all(abs(v - 1) < 1e-12 for v in c.values()),
      f"h = {h:.6f}")
# uniqueness via zeta transform: F*F = 1 and F >= 0 forces zeta F(C) = sqrt2^|C|, F = h^|C|
def zeta(Fn, n): return {C: sum(Fn[A] for A in subsets(n) if all(a <= c for a, c in zip(A, C))) for C in subsets(n)}
def mobius(G, n): return {C: sum((-1)**(sum(C)-sum(A)) * G[A] for A in subsets(n) if all(a <= c for a, c in zip(A, C))) for C in subsets(n)}
n = 4
G = {C: math.sqrt(2)**sum(C) for C in subsets(n)}
Fm = mobius(G, n)
check("unique nonneg square root (zeta/Mobius) is the product dial sqrt2-1",
      all(abs(Fm[C] - h**sum(C)) < 1e-12 for C in subsets(n)))
# zeta is multiplicative over union convolution
ts = [F(1,2), F(3,1), F(2,5)]; ss = [F(1,1), F(1,3), F(4,1)]
cc = conv(lambda a: W(ts, a), lambda b: W(ss, b), 3)
z1 = zeta({a: W(ts, a) for a in subsets(3)}, 3); z2 = zeta({a: W(ss, a) for a in subsets(3)}, 3); z3 = zeta(cc, 3)
check("zeta(F*G) = zeta(F) zeta(G)", all(z3[C] == z1[C]*z2[C] for C in subsets(3)))
# real powers: dial 2^a - 1 and m-th roots
for m in (2, 3, 5):
    t = 2**(1/m) - 1; x = 0
    for _ in range(m): x = multadd(x, t)
    check(f"{m}-th root dial 2^(1/{m})-1 combined by the circle operation {m} times = 1", abs(x - 1) < 1e-12, f"t = {t:.6f}")
# halving commutes with merging (u = log(1+t) line)
t1, t2 = 0.7, 2.3
half = lambda t: math.sqrt(1 + t) - 1
check("half(t ∘₁ s) = half(t) ∘₁ half(s)", abs(half(multadd(t1, t2)) - multadd(half(t1), half(t2))) < 1e-12)

# partner involution and the angle picture
partner = lambda t: (1 - t) / (1 + t)
check("partner(0)=1, partner(1)=0, partner(h)=h, involution",
      partner(0) == 1 and partner(1) == 0 and abs(partner(h) - h) < 1e-12 and abs(partner(partner(0.37)) - 0.37) < 1e-12)
th = 0.3
check("angle: t=tan(th) has partner tan(45deg - th)", abs(partner(math.tan(th)) - math.tan(math.pi/4 - th)) < 1e-12)
check("half dial = tan(22.5 deg)", abs(math.tan(math.pi/8) - h) < 1e-12)
third = 2**(1/3) - 1
ang = math.degrees(math.atan(third))
check("angle picture is for halves only: the third dial sits at 14.57 deg, not 15", abs(ang - 15) > 0.3, f"{ang:.3f} deg")

# 4. complement
ts = [F(2), F(1,3), F(5,2)]
ok = all(W(ts, tuple(1-b for b in a)) == (ts[0]*ts[1]*ts[2]) * W([1/t for t in ts], a) for a in subsets(3))
check("complement turns dials t into 1/t", ok)
p = lambda t: t / (1 + t)
hp_union = 1 - 1/math.sqrt(2); hp_inter = 1/math.sqrt(2)
check("union half in probability: 1 - 1/sqrt2, joined twice gives 1/2", abs(p(h) - hp_union) < 1e-12 and abs(1 - (1-hp_union)**2 - 0.5) < 1e-12)
check("intersection half: 1/sqrt2, intersected twice gives 1/2; dial sqrt2+1 = 1/h",
      abs(hp_inter**2 - 0.5) < 1e-12 and abs(p(math.sqrt(2) + 1) - hp_inter) < 1e-12)
check("halves mirror about 1/2: 0.2929 + 0.7071 = 1", abs(hp_union + hp_inter - 1) < 1e-12)
# odds-against the circle operation for intersection
pa, pb = F(1,3), F(2,5)
against = lambda q: (1 - q) / q
check("intersection: odds-against the circle operation", against(pa*pb) == multadd(against(pa), against(pb)))
check("complement-fixed dial is 1 only (t = 1/t, t > 0)", all(abs(1/t - t) > 1e-9 for t in [0.2, 0.5, h, 2, 3]))

# 5. infinite sets: Bernoulli product measure, density, Borel-Cantelli
random.seed(7)
N = 200000
A = [random.random() < hp_union for _ in range(N)]
B = [random.random() < hp_union for _ in range(N)]
U = [a or b for a, b in zip(A, B)]
dA = sum(A)/N; dU = sum(U)/N
check("infinite set (sampled 200k): half subsets have density 0.293, union of two has 0.5",
      abs(dA - hp_union) < 0.005 and abs(dU - 0.5) < 0.005, f"{dA:.4f}, {dU:.4f}")
# summable dials: chain converges, random subset finite
ts = [1/(k*k) for k in range(1, 100001)]
logZ = sum(math.log1p(t) for t in ts)
check("summable dials 1/k^2: infinite chain converges (prod(1+1/k^2) = sinh(pi)/pi)", abs(math.exp(logZ) - math.sinh(math.pi)/math.pi) < 1e-3, f"{math.exp(logZ):.5f}")
ps = [t/(1+t) for t in ts]
expected = sum(ps)
check("summable dials: expected size of the random subset is finite", expected < 2, f"{expected:.4f}")

# 6. no half step by iteration: brute force every f on {0..4} with f(0)=1
from itertools import product as P2
found = [f for f in P2(range(5), repeat=5) if f[0] == 1 and all(f[f[k]] == 2**k for k in range(2))]
check("no f with f(0)=1 has f(f(0))=1 and f(f(1))=2 (searched all f on 0..4)", not found)

print(f"\n{sum(results)} of {len(results)} pass")
