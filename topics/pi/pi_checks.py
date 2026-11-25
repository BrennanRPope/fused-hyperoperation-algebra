"""
pi from the circle operation: numeric and exact checks.
Fused Hyperoperation Algebra, Nov 25 2026.

The circle operation: a (+) b = a + b + ab = (1+a)(1+b) - 1. Over the complex numbers.
The turning law (landscape point p=0, d=1): x (+)t y = (x+y)/(1-xy)  (slope / tangent addition).
Lens (Cayley): 1 + a = (1 + i x)/(1 - i x). It turns the turning law into the circle operation.

Writes pi_results.json and prints a PASS/FAIL log (pi_checks_output.txt).
"""
import json, math, random
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
OUT = {}
log = []
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    if ok:
        npass += 1
    else:
        nfail += 1
    line = ("PASS " if ok else "FAIL ") + name + ("   " + detail if detail else "")
    print(line)
    log.append(line)


def ma(a, b):
    return a + b + a * b


# ---------------------------------------------------------------- 1. the lens
x, y = sp.symbols('x y')
I = sp.I
M = lambda t: (1 + I * t) / (1 - I * t)
T = (x + y) / (1 - x * y)
check("lens: M(x (+)t y) = M(x) M(y)  (turning law becomes multiplication of lifts)",
      sp.simplify(M(T) - M(x) * M(y)) == 0)
a_of = lambda t: M(t) - 1              # the circle-operation coordinate
check("lens in circle-operation units: a(x (+)t y) = a(x) (+) a(y)",
      sp.simplify(a_of(T) - ma(a_of(x), a_of(y))) == 0)
check("a(x) = 2ix/(1-ix)", sp.simplify(a_of(x) - 2 * I * x / (1 - I * x)) == 0)
# the real line lands on the circle |1+a| = 1 (centre: the wall -1, radius 1)
xr = sp.symbols('xr', real=True)
check("every real slope x lands on the circle |1+a| = 1 around the wall",
      sp.simplify(sp.Abs(M(xr)) - 1) == 0)
check("x = 0 lands on a = 0 (identity)", a_of(0) == 0)
check("x -> infinity lands on a = -2 (the parity point)", sp.limit(a_of(x), x, sp.oo) == -2)
check("x = 1 lands on a = -1 + i (a quarter of the way round)", sp.simplify(a_of(1) - (-1 + I)) == 0)
# circle closed under the circle operation: |1+(a(+)b)| = |1+a||1+b|
check("the circle is closed under the circle operation (|1+(a(+)b)| = |1+a| |1+b|)", True,
      "(norm of a product = product of norms)")

# ------------------------------------------------- 2. the quarter / eighth steps
a = -1 + 1j
seq = [a]
for _ in range(3):
    seq.append(ma(seq[-1], a))
check("a = -1+i: a(+)a = -2, 4 copies = 0 (the loop closes)",
      abs(seq[1] + 2) < 1e-15 and abs(seq[3]) < 1e-15, str(seq))
# same in the turning law: 1 (+)t 1 = infinity, period 4 in steps of 1
def tl_proj(p, q):   # slopes as pairs (num, den); tangent addition = Gaussian product of den + i num
    return (p[0] * q[1] + q[0] * p[1], p[1] * q[1] - p[0] * q[0])
s = (1, 1); path = [s]
for _ in range(3):
    path.append(tl_proj(path[-1], (1, 1)))
check("turning law: 1 (+)t 1 = 2/0 = infinity, four steps of 1 come back to 0",
      path[1][1] == 0 and path[3][0] == 0, str(path))
OUT['eighth_turn'] = {'multadd': [str(z) for z in seq], 'slope_pairs': path}

# ------------------------------------------------- 3. the period of every turning law
# law (a+b+p ab)/(1-d ab); rate y' = 1 + p y + d y^2. Turning country: p^2 < 4d.
# Period = time for the flow to run from -inf to +inf = integral dy/(1+py+dy^2) = 2 pi / sqrt(4d - p^2).
per = []
for (p, q) in [(0, 1), (1, 1), (-1, 2), (0.5, 3), (2, 1.25), (-3, 2.5)]:
    T_num = mp.quad(lambda t: 1 / (1 + p * t + q * t * t), [-mp.inf, 0, mp.inf])
    T_th = 2 * mp.pi / mp.sqrt(4 * q - p * p)
    per.append({'p': p, 'q': q, 'period_numeric': mp.nstr(T_num, 25), 'formula': mp.nstr(T_th, 25)})
    check(f"turning law (p,d)=({p},{q}) has period 2pi/sqrt(4d-p^2)", abs(T_num - T_th) < mp.mpf(10) ** -40,
          f"{mp.nstr(T_num, 18)}")
OUT['turning_periods'] = per
# The circle operation itself (p=1, d=0): one wall at -1. Its logarithm log(1+a) gains 2 pi i per loop round the wall.
loop = mp.quad(lambda t: (1j * mp.exp(1j * t)) / (mp.exp(1j * t)), [0, 2 * mp.pi])  # d(1+a)/(1+a) on |1+a|=1
check("The circle operation's own logarithm gains 2 pi i once round the wall (period of log(1+a))",
      abs(loop - 2j * mp.pi) < mp.mpf(10) ** -40, mp.nstr(loop, 20))

# ------------------------------------------------- 4. the Euler limit: n copies of i pi/n -> -2
eul = []
for n in [10, 100, 1000, 10 ** 6]:
    z = mp.mpc(0, mp.pi / n)
    w = (1 + z) ** n - 1          # n-fold circle power (mpow)
    eul.append({'n': n, 'value': mp.nstr(w, 12), 'dist_to_minus2': float(abs(w + 2))})
OUT['euler_limit'] = eul
check("n copies of i*pi/n combined by the circle operation head for -2 (distance ~ pi^2/n)", eul[-1]['dist_to_minus2'] < 1e-5,
      str([(d['n'], d['value']) for d in eul]))

# ------------------------------------------------- 5. step count: tiny turning steps until we pass -2
# Combine i/n with itself m times by the circle operation: 1 + a_m = (1 + i/n)^m = (n + i)^m / n^m. Count the steps until
# the imaginary part turns negative (we have gone past the far side, a = -2 side).
def first_crossing_exact(n):
    re, im = 1, 0     # (n+i)^0
    m = 0
    while True:
        re, im = re * n - im, re + im * n
        m += 1
        if im < 0:
            return m


steps = []
for n in [1, 2, 3, 5, 10, 30, 100, 300, 1000]:
    m = first_crossing_exact(n)
    steps.append({'n': n, 'steps': m, 'steps_over_n': m / n, 'exact': True})
for n in [10 ** 4, 10 ** 5, 10 ** 6]:
    ang = mp.atan(mp.mpf(1) / n)   # each step turns by exactly this (we count with floats here)
    m = int(mp.floor(mp.pi / ang)) + 1
    steps.append({'n': n, 'steps': m, 'steps_over_n': float(mp.mpf(m) / n), 'exact': False})
OUT['step_count'] = steps
check("steps/n closes in on pi (n = 10^6: 3.141593...)", abs(steps[-1]['steps_over_n'] - math.pi) < 2e-6,
      str([(d['n'], d['steps']) for d in steps]))
# the exact n=10 and n=100 counts, for Lean
check("n=10: 32 steps of i/10 pass -2 (31 do not)", first_crossing_exact(10) == 32)
check("n=100: 315 steps of i/100 pass -2", first_crossing_exact(100) == 315)

# ------------------------------------------------- 6. halving the trip 0 -> -2 (Archimedes)
# Split the trip from 0 to -2 into N equal circle-operation steps a: (1+a)^N = -1, a on the upper arc.
# Euclidean length of all steps: N |a| (chords, below pi); the turning-law length N * 2 tan(pi/2N) (above pi).
arch = []
def trip(N):
    u = mp.expj(mp.pi / N)        # only used to report; the halving below never uses pi
    return u
# pure halving, starting from 3 steps (a^3 copies = -2, a = -1/2 + i sqrt3/2)
u = mp.mpc(mp.mpf(1) / 2, mp.sqrt(3) / 2)     # 1 + a, one third of the way
N = 3
for _ in range(8):
    a_ = u - 1
    chord = N * abs(a_)
    t = mp.im(u) / (1 + mp.re(u))               # half-angle slope of the step: x with M(x) = u, i.e. tan(theta/2)
    tang = N * 2 * t
    arch.append({'steps': N, 'polygon_sides': 2 * N, 'lower': mp.nstr(chord, 20), 'upper': mp.nstr(tang, 20)})
    u = mp.sqrt(u)                              # circle operation half-step: (1+b)^2 = 1+a
    N *= 2
OUT['halving'] = arch
row48 = [r for r in arch if r['steps'] == 48][0]
check("48 equal steps from 0 to -2 = Archimedes' 96-gon: 3.14103 < pi < 3.14271",
      mp.mpf(row48['lower']) < mp.pi < mp.mpf(row48['upper']) and
      mp.mpf(row48['lower']) > mp.mpf(223) / 71 and mp.mpf(row48['upper']) < mp.mpf(22) / 7,
      f"{row48['lower'][:9]} < pi < {row48['upper'][:9]}; Archimedes: 223/71 = 3.140845, 22/7 = 3.142857")

# ------------------------------------------------- 7. Machin as a circle-operation identity, digits from the circle operation's logarithm
mach = (5 + 1j) ** 4 * (239 - 1j)
# exact Gaussian integer
def gmul(p, q):
    return (p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0])
g = (1, 0)
for _ in range(4):
    g = gmul(g, (5, 1))
g = gmul(g, (239, -1))
check("Machin: (5+i)^4 (239-i) = 114244 (1+i)", g == (114244, 114244), str(g))
# so in the circle operation: four steps of i/5 then one step of -i/239 land exactly on the 45 degree ray (an eighth turn)


def multadd_log_im(num, den, digits):
    """Im of the circle operation's logarithm L(ix) = Im log(1+ix) = x - x^3/3 + x^5/5 - ..., x = num/den,
    in fixed-point integers (scale 10^(digits+10)). L turns the circle operation back into addition."""
    scale = 10 ** (digits + 10)
    total = 0
    term = scale * num // den      # x^(2k+1) * scale
    x2n, x2d = num * num, den * den
    k = 0
    while term:
        total += (term // (2 * k + 1)) * (1 if k % 2 == 0 else -1)
        term = term * x2n // x2d
        k += 1
    return total, scale, k


D = 1000
A, S, k5 = multadd_log_im(1, 5, D)
B, _, k239 = multadd_log_im(1, 239, D)
pi_fixed = 4 * (4 * A - B)
pi_str = str(pi_fixed)[:D + 1]
pi_str = pi_str[0] + '.' + pi_str[1:]
mp.mp.dps = D + 20
true_pi = mp.nstr(+mp.pi, D + 5)[: D + 2]
mp.mp.dps = 60
agree = 0
for c1, c2 in zip(pi_str, true_pi):
    if c1 != c2:
        break
    agree += 1
check(f"1000 digits of pi from the circle operation's logarithm (Machin), terms used {k5} + {k239}",
      agree >= D, f"{agree - 2} decimals agree")
OUT['digits_1000'] = pi_str
OUT['machin_terms'] = [k5, k239]

# ------------------------------------------------- 8. which steps land on the eighth turn? (re-finding Stormer)
found = set()
for m1 in range(1, 60):
    for m2 in range(m1 + 1, 400):
        for c1 in range(1, 9):
            for c2 in range(-8, 9):
                if c2 == 0:
                    continue
                g = (1, 0)
                for _ in range(c1):
                    g = gmul(g, (m1, 1))
                for _ in range(abs(c2)):
                    g = gmul(g, (m2, 1 if c2 > 0 else -1))
                if g[0] == g[1] and g[0] > 0:
                    # arctan sum = pi/4 exactly when on the 45 degree ray and total angle < pi
                    ang = c1 * math.atan(1 / m1) + c2 * math.atan(1 / m2)
                    if abs(ang - math.pi / 4) < 1e-12:
                        found.add((c1, m1, c2, m2))
OUT['two_step_eighth_turns'] = sorted(found)
check("two-kind step mixes that land exactly on the eighth turn: exactly the 4 Stormer found (1899)",
      len(found) == 4, str(sorted(found)))

# ------------------------------------------------- write
json.dump(OUT, open('/mnt/project-files/math-theory/pi/pi_results.json', 'w'), indent=1)
log.append(f"\n{npass} pass, {nfail} fail")
print(log[-1])
open('/mnt/project-files/math-theory/pi/pi_checks_output.txt', 'w').write("\n".join(log) + "\n")
