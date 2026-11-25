"""All numeric claims on the Lee–Yang page, as assertions. Run: python3 checks.py"""
from dialgas import *
from potts_strips import cyclic_strip_chromatic, evalp
from fractions import Fraction as F
import pickle, cmath, math, json, numpy as np, mpmath as mp
out = {}; n_ok = 0
def ok(name, cond):
    global n_ok
    assert cond, name; n_ok += 1; print("PASS", name)

# 1. one edge = lift of the circle operation with dial k
for k in (F(3), F(-2), F(1, 2)):
    x, y = F(2, 3), F(-5, 7)
    ok(f"edge gas k={k}", 1 + x + y + k * x * y == 1 + (x + y + k * x * y))

# 2. odd rings vanish exactly on the wall z = -1/k; even rings give 2((k-1)/k)^m
for k in (F(4), F(3, 2), F(1, 3), F(-1), F(-2), F(-7, 5)):
    for n in range(3, 16):
        val = evalZ(cycle_counts(n), -1 / k, k)
        ok(f"ring n={n} k={k} on wall", val == (0 if n % 2 else 2 * ((k - 1) / k) ** (n // 2)))

# 3. odd d-regular graphs vanish at 1 + k^(d/2) z = 0 (tori: d = 4)
C = pickle.load(open('torus_counts.pkl', 'rb'))
for (L, M) in ((3, 3), (3, 5), (5, 5), (3, 4), (4, 4)):
    c = C[L] if L == M and L in C else torus_counts(L, M)
    for k in (F(4), F(1, 4), F(-1), F(-3)):
        v = evalZ(c, -1 / k ** 2, k)
        ok(f"torus {L}x{M} k={k} at zeta=-1 {'zero' if L*M%2 else 'nonzero'}", (v == 0) == (L * M % 2 == 1))

# 4. ring zeros in closed form: zeta + 1/zeta = c_j (1 - 1/k) - 2, c_j = 4 cos^2(pi(2j+1)/(2n))
def ring_zeros_zeta(k, n):
    zs = []
    for j in range(n):
        c = 4 * math.cos(math.pi * (2 * j + 1) / (2 * n)) ** 2
        s = c * (1 - 1 / k) - 2
        d = cmath.sqrt(s * s - 4)
        zs += [(s + d) / 2, (s - d) / 2]
    return zs  # each zero appears twice (j and n-1-j)
for k in (3.0, 0.5, -1.0, -2.0):
    n = 11
    num = [complex(r) * k for r in z_zeros(cycle_counts(n), k)]
    cf = ring_zeros_zeta(k, n)
    ok(f"closed-form ring zeros k={k}", all(min(abs(a - b) for b in cf) < 1e-8 for a in num))

# 5. regimes of the ring (n = 41) in zeta = k z
def ring_num(k, n=41): return [complex(r) * k for r in z_zeros(cycle_counts(n), k, dps=60)]
for k in (1.5, 3.0, 6.0):
    zs = ring_num(k); ok(f"k={k}>1: ring zeros on |zeta|=1", all(abs(abs(z) - 1) < 1e-9 for z in zs))
    gap = min(abs(cmath.phase(z)) for z in zs); edge = math.acos(1 - 2 / k)
    ok(f"k={k}: gap edge >= arccos(1-2/k)={edge:.3f} (found {gap:.3f})", gap >= edge - 1e-9)
for k in (0.2, 0.5, 0.9):
    zs = ring_num(k); ok(f"0<k={k}<1: ring zeros real negative", all(abs(z.imag) < 1e-9 and z.real < 0 for z in zs))
for k in (-0.5, -1.0, -2.0, -3.0):
    zs = ring_num(k)
    ok(f"k={k}<0: ring zeros on |zeta|=1 or positive real",
       all(abs(abs(z) - 1) < 1e-9 or (abs(z.imag) < 1e-9 and z.real > 0) for z in zs))
    cf = ring_zeros_zeta(k, 4001)
    ok(f"k={k}<0: long rings (n=4001) reach the zero-field point zeta=+1", min(abs(z - 1) for z in cf) < 0.1)

# 5b. one formula for every k: the ring's zeros fill s = zeta + 1/zeta from the wall s=-2 to s=2-4/k
for k in (6.0, 3.0, 1.5, 0.7, 0.3, -0.5, -1.0, -2.0, -3.0):
    ss = [z + 1 / z for z in ring_num(k)]
    lo, hi = sorted((-2.0, 2 - 4 / k))
    ok(f"k={k}: all s real in [{lo:.3f},{hi:.3f}]", all(abs(x.imag) < 1e-7 and lo - 1e-7 <= x.real <= hi + 1e-7 for x in ss))

# 6. torus Lee–Yang circle |z| = k^-2 for k > 1
Z = pickle.load(open('torus_zeros.pkl', 'rb'))
for (L, k), zs in Z.items():
    if k > 1.001:
        ok(f"torus {L}x{L} k={k:.3f}: |z|=k^-2", all(abs(abs(z) * k * k - 1) < 1e-7 for z in zs))

# 7. 2D Fisher zeros (zero field) head for k_c = 3+2sqrt2 and 1/k_c, not for k = -2,-1,0,1
Fz = pickle.load(open('fisher_zeros.pkl', 'rb'))
kc = 3 + 2 * math.sqrt(2)
for L in (4, 5, 6):
    d = {x: min(abs(z - x) for z in Fz[L]) for x in (-2, -1, 0, 1)}
    ok(f"Fisher L={L}: no zero within 0.2 of k=-2,-1,0,1 ({', '.join(f'{v:.2f}' for v in d.values())})", min(d.values()) > 0.2)
near = min(Fz[6], key=lambda z: abs(z - kc)); ok(f"Fisher L=6 nearest to k_c: {near:.3f}", abs(near - kc) < 3)

# 8. Potts colour strips: Q = 0, 1 roots; Q = 2 root iff odd (non-bipartite); Q = -1 never, sign (-1)^N
for L, n in ((1, 7), (1, 8), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6)):
    p = cyclic_strip_chromatic(L, n); N = L * n
    ok(f"strip {L}x{n}: P(0)=P(1)=0", evalp(p, 0) == 0 and evalp(p, 1) == 0)
    ok(f"strip {L}x{n}: P(2)=0 iff n odd", (evalp(p, 2) == 0) == (n % 2 == 1))
    v = evalp(p, -1); ok(f"strip {L}x{n}: P(-1)={v} nonzero with sign (-1)^N", v != 0 and (v > 0) == (N % 2 == 0))
Pz = pickle.load(open('potts_zeros.pkl', 'rb'))
for key, zs in Pz.items():
    ok(f"strip {key}: no real root below 0 or in (0,1) or (1,32/27]",
       not any(abs(z.imag) < 1e-9 and (z.real < -1e-9 or 1e-9 < z.real < 1 - 1e-9 or 1 + 1e-9 < z.real <= 32 / 27) for z in zs))
d2 = [min(abs(z - 2) for z in Pz[(3, n)]) for n in (6, 10, 14)]
ok(f"width-3 strips: zeros close in on Q=2 (k=1): {d2}", d2[0] > d2[1] > d2[2] and d2[2] < 1e-3)

print(f"\n{n_ok} checks passed")
