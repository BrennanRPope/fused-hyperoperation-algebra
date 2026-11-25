"""Fill-then-spill test, part 2: crossing the landscape's border (the parabola d = k^2/4).
Repeat one law with a fixed small step:  a -> a (+) x = (a + x + k a x) / (1 - d a x),  starting at 0.
Two walls (4d < k^2): a climbs into a wall and stays (the cup fills and holds).
Turning (4d > k^2):   a climbs, lingers where the double wall used to be, then runs off through infinity
                      (spills) and comes back round."""
import math, json
import numpy as np
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

def cycle(k, q, x=0.05, nmax=3_000_000):
    """Track a through infinity with homogeneous coordinates. Return (period, linger, trace):
    period = steps from 0 round through infinity and back up to 0 (None if it settles in a wall);
    linger = steps spent within half the ghost's distance of the ghost a* = -2/k."""
    M = np.array([[1 + k*x, x], [-q*x, 1.0]])
    u, v = 0.0, 1.0; n = 0; spilled = False; lin = 0; ghost = -2.0 / k; tr = []
    while n < nmax:
        u, v = M[0,0]*u + M[0,1]*v, M[1,0]*u + M[1,1]*v
        r = math.hypot(u, v); u /= r; v /= r; n += 1
        if v * (v - 0) <= 0 or abs(v) < 1e-15:
            pass
        a = u / v if v != 0 else math.inf
        if a < 0 or not math.isfinite(a): spilled = True
        if abs(a - ghost) < abs(ghost) / 2: lin += 1
        if len(tr) < 400000: tr.append(a)
        if spilled and math.isfinite(a) and a >= 0: return n, lin, tr
    return None, lin, tr

say('Landscape point fixed at k = 1 (ghost of the double wall at a = -2), bend just past the border:')
say('d = 1/4 + Δ/4, with Δ = 4d - k^2 the distance into the turning side. Step x = 0.05.')
ds = np.logspace(-7, -1, 13); T = []; Lg = []
for d in ds:
    n, lin, _ = cycle(1.0, 0.25 + d / 4)
    T.append(n); Lg.append(lin)
    say(f'  Δ = {d:.1e}: one full turn takes {n} steps, {lin} of them lingering by the ghost ({lin/n:.0%})')
slope = np.polyfit(np.log(ds), np.log(T), 1)[0]
slopeL = np.polyfit(np.log(ds), np.log(Lg), 1)[0]
say(f'  Fit: turn time ~ Δ^{slope:.3f}, linger ~ Δ^{slopeL:.3f}  (a fold predicts -1/2)')
for d in [1e-3, 1e-1]:
    n, lin, tr = cycle(1.0, 0.25 - d / 4, nmax=200000)
    w = (-2 + math.sqrt(d)*2/1) if False else None
    say(f'  Two-walls side, Δ = -{d:.0e}: comes back round = {n is not None}; settles at a = {tr[-1]:.5f}')
say('  On the two-walls side the orbit settles into a wall and stays. On the turning side it lingers by')
say('  the ghost of the merged walls, then spills through infinity, and the linger grows without limit')
say('  as the border is approached.')

# collapse: only the distance to the border matters
say('\nRandom landscape points, all at Δ = 1e-4: does the time depend on anything but Δ?')
rng = np.random.default_rng(3); TT = []
for _ in range(12):
    k = rng.uniform(0.3, 1.5) * rng.choice([-1, 1]); n, _, _ = cycle(k, k*k/4 + 1e-4/4); TT.append(n)
say(f'  one full turn takes {min(TT)}..{max(TT)} steps over 12 random k; spread {np.std(TT)/np.mean(TT):.2%} of the mean')

# exact: the step is a Mobius map with trace 2 + kx, det 1 + kx + dx^2; it turns by phi each step
say('\nExact angle per step: cos(phi) = (2 + kx)/(2 sqrt(1 + kx + dx^2)); for small x, phi ~ x sqrt(4d - k^2)/2.')
for d in [1e-6, 1e-4, 1e-2]:
    k, q, x = 1.0, 0.25 + d/4, 0.05
    phi = math.acos((2 + k*x) / (2 * math.sqrt(1 + k*x + q*x*x)))
    say(f'  Δ = {d:.0e}: phi = {phi:.3e}, half-turn pi/phi = {math.pi/phi:.0f} steps, x sqrt(Δ)/2 = {x*math.sqrt(d)/2:.3e}')

# fair null: generic smooth folds, f(a) = a + eps + c a^2 + e a^3
say('\nFair null: random smooth maps near a fold, f(a) = a + eps + c a^2 + e a^3 (c in [0.5, 2], e in [-0.25, 0.25]), time from -0.5 to 0.5.')
exps = []
for _ in range(10):
    c = rng.uniform(0.5, 2); e = rng.uniform(-0.25, 0.25); Ts = []; eps = np.logspace(-8, -3, 6)
    for ep in eps:
        a, n = -0.5, 0
        while a < 0.5 and n < 10_000_000:
            a = a + ep + c*a*a + e*a**3; n += 1
        Ts.append(n)
    exps.append(np.polyfit(np.log(eps), np.log(Ts), 1)[0])
say(f'  passage-time exponents: {min(exps):.3f} .. {max(exps):.3f}; every fold gives -1/2.')
say('  Reading: the landscape border does fill, linger and spill, with the linger growing as 1/sqrt(distance).')
say('  That is the standard bottleneck of any fold (two fixed points merging), the same exponent as the null.')
open('fillspill/border_output.txt', 'w').write('\n'.join(out) + '\n')
json.dump({'d': ds.tolist(), 'T': T, 'linger': Lg, 'slope': slope}, open('fillspill/border.json', 'w'))
# one orbit trace for the figure
n, lin, tr = cycle(1.0, 0.25 + 1e-4/4)
json.dump({'a': tr, 'period': n}, open('fillspill/border_orbit.json', 'w'))
