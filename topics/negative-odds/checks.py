"""Negative odds and the wall at t = -1. Exact checks with Fractions, numeric with numpy.
Odds t, probability p = t/(1+t). Union of independent items: t (+) s = t + s + ts (k = 1)."""
from fractions import Fraction as F
import numpy as np, cmath, math
ok = 0; fails = []
def check(name, cond):
    global ok
    if cond: ok += 1
    else: fails.append(name)
    print(("PASS " if cond else "FAIL ") + name)

ma = lambda t, s: t + s + t*s              # odds union law, k = 1
prob = lambda t: t/(1+t)
odds = lambda p: p/(1-p)
rewind = lambda t: -t/(1+t)                # quasi-inverse ("way back")
mirror = lambda t: 1/t                     # complement
third = lambda t: -1 - t                   # p -> 1/p
inter = lambda t, s: 1/ma(1/t, 1/s)        # intersection via complement

# 1. union of independent events has odds t (+) s
import random
random.seed(1)
for _ in range(200):
    p, q = F(random.randint(1,99),100), F(random.randint(1,99),100)
    u = 1-(1-p)*(1-q)
    assert odds(u) == ma(odds(p), odds(q))
check("union of independent events: odds compose by t+s+ts (200 exact cases)", True)

# 2. the wall -1 swallows every finite odds
check("wall absorbs union: (-1)(+)s = -1 for s in a sample", all(ma(F(-1), F(s,7)) == -1 for s in range(-30,30)))
# 3. wall is fixed by the complement mirror and absorbs intersection too
check("mirror fixes -1 and +1 only (on reals)", mirror(F(-1)) == -1 and mirror(F(1)) == 1)
check("wall absorbs intersection: inter(-1, s) = -1 for s != 0", all(inter(F(-1), F(s,5)) == -1 for s in range(-20,21) if s not in (0,-5)))
# 4. regions under the mirror: (-1,0) <-> (-inf,-1)
check("mirror swaps p<0 (t in (-1,0)) with p>1 (t<-1)", all(prob(F(-1, n)) < 0 and prob(mirror(F(-1,n))) > 1 for n in range(2,40)))
check("mirror is complement: prob(1/t) = 1 - prob(t)", all(prob(mirror(F(n,7))) == 1 - prob(F(n,7)) for n in range(-30,30) if n not in (0,-7)))
# 5. rewinds of real chances land in (-1,0); rewind probability = minus original odds
pts = [F(n,9) for n in range(1,60)]
check("rewind of any positive odds lies in (-1,0)", all(-1 < rewind(t) < 0 for t in pts))
check("rewind cancels: t (+) rewind(t) = 0", all(ma(t, rewind(t)) == 0 for t in pts))
check("prob(rewind(t)) = -t  (the undo of chance p has 'chance' minus its odds)", all(prob(rewind(t)) == -t for t in pts))
check("example: a -1/3 'chance' cancels a 1/4 chance (p+q-pq = 0)", F(1,4) + F(-1,3) - F(1,4)*F(-1,3) == 0)
# 6. the sign sector: t = -2 (p = 2) is its own undo
check("t = -2 is self-inverse: (-2)(+)(-2) = 0, i.e. 200% twice = never", ma(F(-2), F(-2)) == 0 and prob(F(-2)) == 2)
check("t < -1 has 1+t < 0: sign flips with each such item", all((1+ma(F(-n-1,1), F(-m-1,1))) > 0 for n in range(1,9) for m in range(1,9)))
# 7. the six-map symmetry permutes never (0), certain (inf), wall (-1)
def onP(f, x):
    try: return f(x)
    except ZeroDivisionError: return 'inf'
# compose as Mobius matrices
M = {'id':np.array([[1,0],[0,1]]), 'mirror':np.array([[0,1],[1,0]]), 'rewind':np.array([[-1,0],[1,1]]), 'third':np.array([[-1,-1],[0,1]])}
def norm(A):
    A = np.array(A, dtype=float); i = np.flatnonzero(np.abs(A.ravel())>1e-12)[0]; return np.round(A/A.ravel()[i], 9)
gen = [M['mirror'], M['rewind']]
group = [np.eye(2)]
changed = True
while changed:
    changed = False
    for g in list(group):
        for h in gen:
            c = norm(g @ h)
            if not any(np.allclose(c, norm(x)) for x in group): group.append(c); changed = True
check("mirror and rewind generate a group of exactly 6 maps (cross-ratio group)", len(group) == 6)
r3 = lambda t: rewind(mirror(t))
check("rewind after mirror = -1/(1+t), order 3", all(r3(r3(r3(F(n,7)))) == F(n,7) for n in range(1,20)))
check("order-3 map fixes the complex cube roots of unity (t^2+t+1=0)", abs((lambda w: -1/(1+w) - w)(cmath.exp(2j*math.pi/3))) < 1e-12)
check("third involution t -> -1-t sends p to 1/p", all(prob(third(F(n,7))) == 1/prob(F(n,7)) for n in range(1,30)))
check("rewind swaps certain(inf) and wall: rewind(t) -> -1 as t -> inf", abs(rewind(1e12) + 1) < 1e-9)
check("third swaps never and wall: third(0) = -1, third(-1) = 0", third(F(0)) == -1 and third(F(-1)) == 0)
check("p = 2 point (t=-2) is fixed by rewind; p = 1/2 (t=1) fixed by mirror; p = inf (t=-1) fixed by mirror",
      rewind(F(-2)) == -2 and mirror(F(1)) == 1 and mirror(F(-1)) == -1)

# 8. the clock: t(T) = e^{rT} - 1, union = adding time; negative odds = negative time
r = 0.7
clk = lambda T: math.expm1(r*T)
check("clock: t(T1)(+)t(T2) = t(T1+T2)", all(abs(ma(clk(a), clk(b)) - clk(a+b)) < 1e-9 for a in np.linspace(-3,3,13) for b in np.linspace(-3,3,13)))
check("clock: rewind(t(T)) = t(-T)", all(abs(rewind(clk(a)) - clk(-a)) < 1e-12 for a in np.linspace(-3,3,13)))
check("clock: negative time gives odds in (-1,0), wall at T -> -inf", all(-1 < clk(-a) < 0 for a in np.linspace(0.1,30,50)) and abs(clk(-60)+1) < 1e-12)
# simulation: event = at least one tick of rate-r Poisson clock in window T
rng = np.random.default_rng(0)
T = 1.3; ticks = rng.poisson(r*T, 200000)
p_sim = np.mean(ticks > 0)
check("Poisson simulation: P(tick in T) has odds e^{rT}-1 (to 0.01)", abs(odds(p_sim) - clk(T)) < 0.01)
check("no real time reaches t < -1: needs T with imaginary part pi/r", abs(cmath.exp(r*(1.0 + 1j*math.pi/r)) - 1 - (-math.exp(r)-1)) < 1e-9)

# 9. sign of k = time direction (Tsallis q <-> 2-q duality): x_{-k}(H) = -x_k(-H)
xk = lambda k, H: math.expm1(k*H)/k
check("x_{-k}(H) = -x_k(-H) for a grid of k, H", all(abs(xk(-k,H) + xk(k,-H)) < 1e-12 for k in (0.3,1,2.5) for H in np.linspace(-2,2,9)))
flip = lambda k, x: x/(1+k*x)
check("flip x/(1+kx) = -(rewind in law k) and carries law k to law -k",
      all(abs(flip(1, x) + rewind(x)) < 1e-12 for x in (0.2, 1.5, -0.4)) and
      all(abs(flip(1, ma(a,b)) - (flip(1,a)+flip(1,b)-flip(1,a)*flip(1,b))) < 1e-12 for a in (0.2,1.1) for b in (0.5,3.0)))

# 10. thermodynamic arrow (Crooks/detailed balance) acts by the mirror, not by negation
# two-state system, energy gap E, inverse temperature b: forward/back rate ratio = odds = e^{-bE}
b, E = 1.0, 0.8
t_fwd = math.exp(-b*E)
check("detailed balance: reversing the process (or E -> -E, Feynman-Stueckelberg) sends odds t -> 1/t",
      abs(math.exp(-b*(-E)) - 1/t_fwd) < 1e-12)
check("so time reversal keeps the sign of t: t>0 -> 1/t>0", 1/t_fwd > 0)
# Crooks: P(+s)/P(-s) = e^{s} for a Gaussian work distribution with var = 2 mean (detailed FT)
mu = 1.2; xs = np.array([0.5,1.0,2.0])
g = lambda x: np.exp(-(x-mu)**2/(4*mu))
check("Crooks ratio P(+s)/P(-s) = e^s (Gaussian with variance 2*mean)", np.allclose(g(xs)/g(-xs), np.exp(xs)))
# fixed points of log-odds reversal lam -> -lam (mod 2 pi i): lam = 0 (t=1) and lam = i pi (t=-1)
check("time-reversal fixed points on the complex log-odds: t = 1 and t = -1 (log(-1) = i pi = -i pi mod 2 pi i)",
      abs(cmath.exp(-1j*math.pi) - cmath.exp(1j*math.pi)) < 1e-12)

# 11. Fermi-Dirac: occupation n = t/(1+t) = our prob; pole at t=-1 = Lee-Yang zero of Z = 1+t; Matsubara
check("fermion mode: Z = 1+t vanishes at the wall t = -1", 1 + (-1) == 0)
beta = 2.0
for n in range(-3,4):
    w = (2*n+1)*math.pi/beta
    assert abs(cmath.exp(-beta*(1j*w)) + 1) < 1e-12
check("t = e^{-beta*E} hits -1 exactly at E = i(2n+1)pi/beta (fermionic Matsubara points)", True)
check("boson mode: n = t/(1-t) has its pole at t = +1 (bosonic Matsubara 2n pi/beta)", all(abs(cmath.exp(-beta*1j*2*n*math.pi/beta) - 1) < 1e-12 for n in range(-3,4)))

# 12. game mechanic: odds multipliers m = 1+t, final p = 1 - 1/prod(m)
mult = [F(4,3), F(2), F(3,4)]   # 1/4 chance, coin flip, a -1/3 curse
m = math.prod(mult); p_final = 1 - 1/m
direct = F(0)
for x in [F(1,4), F(1,2), F(-1,3)]: direct = direct + x - direct*x
check("multiplier rule equals stacking p+q-pq: 1/4, 1/2, -1/3 -> 1/2", p_final == direct == F(1,2))
# doom/destiny with saved information: m = c * eps^n
def resolve(tokens):
    c, n = F(1), 0
    for tok in tokens:
        if tok == 'doom': n += 1
        elif tok == 'destiny': n -= 1
        else: c *= tok
    if n > 0: return 'impossible'
    if n < 0: return 'certain'
    return max(F(0), min(F(1), 1 - 1/c))
check("doom+destiny cancel and ordinary odds decide (x2, x4/3 -> p = 5/8)", resolve(['doom','destiny',F(2),F(4,3)]) == F(5,8))
check("two doom vs one destiny -> impossible", resolve(['doom','doom','destiny',F(5)]) == 'impossible')
check("plain 0*inf would be undefined; counting copies always gives an answer (2000 random hands)",
      all(resolve(random.choices(['doom','destiny',F(2),F(3,4),F(4,3)], k=6)) is not None for _ in range(2000)))

print(f"\n{ok} passed, {len(fails)} failed {fails}")
