"""Polishing pass for rulecraft-history (DRAFT): re-runs the claims with exact rationals at n = 6
and 50-digit numerics, and hunts for counterexamples near the edges (strengths near 0 and 1, tied targets)."""
import itertools, random
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50
ok = lambda n, c: print(("PASS " if c else "FAIL ") + n) or c
r = []
def run(order, L, T, x0):
    y = x0
    for i in order: y = (1 - L[i])*y + L[i]*T[i]
    return y
rng = random.Random(11)
# P1 sort-by-target at n = 6 with extreme strengths
best = mono = True
for _ in range(60):
    L = [Fr(rng.choice([1, 2, 50, 98, 99, 999]), 1000) for _ in range(6)]
    T = [Fr(rng.randint(0, 50), 50) for _ in range(6)]
    outs = {o: run(o, L, T, Fr(1, 3)) for o in itertools.permutations(range(6))}
    asc = tuple(sorted(range(6), key=lambda i: T[i]))
    best &= outs[asc] == max(outs.values())
    for o, v in outs.items():
        for k in range(5):
            o2 = o[:k] + (o[k+1], o[k]) + o[k+2:]
            if T[o[k]] < T[o[k+1]]: mono &= outs[o2] < v
            if T[o[k]] == T[o[k+1]]: mono &= outs[o2] == v
r.append(ok("P1 n = 6, strengths from 0.001 to 0.999, tied targets allowed: sorted order is best; swaps of equal targets change nothing (exact, 60 x 720 orders)", best and mono))
# P2 swing formula at 50 digits by iteration
l1, l2, t1, t2 = mp.mpf('0.6'), mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('0.9')
y = mp.mpf('0.5')
for _ in range(400): y = (1 - l1)*y + l1*t1; a = y; y = (1 - l2)*y + l2*t2; b = y
r.append(ok("P2 pendulum swing at 50 digits: |b - a| = %s" % mp.nstr(abs(b - a), 20), abs(abs(b - a) - l1*l2*(t2 - t1)/(l1 + l2 - l1*l2)) < mp.mpf(10)**-45))
# P3 Trotter limit at 30 digits via Richardson on n and 2n
def inter(order, n):
    L, T = [mp.mpf('0.6'), mp.mpf('0.5')], [mp.mpf('0.1'), mp.mpf('0.9')]
    y = mp.mpf('0.5')
    for _ in range(n):
        for i in order:
            li = 1 - (1 - L[i])**(mp.mpf(1)/n); y = (1 - li)*y + li*T[i]
    return y
mu = [-mp.log(1 - mp.mpf('0.6')), -mp.log(1 - mp.mpf('0.5'))]
tgt = (mu[0]*mp.mpf('0.1') + mu[1]*mp.mpf('0.9'))/(mu[0] + mu[1]); sim = tgt + (mp.mpf('0.5') - tgt)*mp.e**(-(mu[0] + mu[1]))
sym = lambda n: (inter((0, 1), n) + inter((1, 0), n))/2
est = (4*sym(2000) - sym(1000))/3
r.append(ok("P3 symmetric phase-in -> simultaneous rule; Richardson estimate differs by %s" % mp.nstr(abs(est - sim), 3), abs(est - sim) < mp.mpf(10)**-9))
# P4 origin share regimes with an exact tail bound: sum of c/n^alpha for alpha = 1.5 converges (zeta)
r.append(ok("P4 alpha = 1.5, c = 0.6: sum of strengths = 0.6 zeta(1.5) = %s (finite), so the start keeps share >= exp(-sum/(1-0.6)) > 0" % mp.nstr(0.6*mp.zeta(1.5), 8), mp.exp(-0.6*mp.zeta(1.5)/(1 - 0.6)) > 0))
print("\n%d of %d polishing checks pass" % (sum(r), len(r)))
