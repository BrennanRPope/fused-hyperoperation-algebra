"""Fill-then-spill test, part 1: the rungs.
Does the passage from one whole rung to the next fill (slow down near the rung) and then spill?"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fractional-rungs'))
from regular_iteration import Es, coeffs_at
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

# (a) landscape coordinates of rung s: k = s, d = s(s-1)/4 ; check d from the third coefficient
say('(a) Rung s on the landscape')
for s in [-1, -0.5, 0, 0.25, 0.5, 1, 1.5, 2, 2.5]:
    c = coeffs_at(s)
    q = -s * s / 2 + 3 * c[3]
    say(f'  s = {s:5}: corner k = {2*c[2]:.6f}, bend d = {q:.6f}, s(s-1)/4 = {s*(s-1)/4:.6f}, shape k^2/(k^2-4d) = '
        + (f'{(2*c[2])**2/((2*c[2])**2-4*q):.6f}' if abs((2*c[2])**2-4*q) > 1e-12 else 'border'))
say('  The rung curve is the parabola (s, s(s-1)/4); its shape number is s itself, so it moves at')
say('  constant speed in the unit-free coordinate and crosses the border once, at s = 0, linearly.')

# (b) flow speed: v(s) = d/ds E^s(x)
say('\n(b) Speed of the flow s -> E^s(x) between whole rungs')
S = np.linspace(-2.0, 2.0, 801); h = 1e-3
def speed(x, theta=lambda s: s, dtheta=lambda s: 1.0):
    vals = []
    for s in S:
        t = theta(s)
        f1 = Es(t + h, np.array([x]))[0]; f0 = Es(t - h, np.array([x]))[0]
        vals.append((f1 - f0) / (2 * h) * dtheta(s))
    return np.array(vals)

def stair_index(v, Sg=None):
    """Fit log(speed) = smooth trend (degree-6 polynomial) + A cos(2 pi s) + B sin(2 pi s).
    A < 0 means the flow slows at every whole rung (a staircase). Returns A."""
    Sg = S if Sg is None else Sg
    lv = np.log(np.maximum(v, 1e-300))
    X = np.column_stack([Sg**j for j in range(7)] + [np.cos(2*np.pi*Sg), np.sin(2*np.pi*Sg)])
    coef, *_ = np.linalg.lstsq(X, lv, rcond=None)
    return float(coef[7])

res = {}
for x in [0.5, 1.0, 2.0]:
    if x == 2.0: S = np.linspace(-2.0, 1.5, 701)
    v = speed(x)
    mono = bool(np.all(np.diff(v) > 0))
    res[x] = stair_index(v)
    say(f'  x = {x}: speed rises the whole way (no slow-down at any rung): {mono}; rung-locked wobble A = {res[x]:+.5f}')
    S = np.linspace(-2.0, 2.0, 801)

# Control: a valid but staircased fractional iterate E^{theta(s)}, theta(s) = s - 0.9 sin(2 pi s)/(2 pi)
th = lambda s: s - 0.9 * math.sin(2 * math.pi * s) / (2 * math.pi)
dth = lambda s: 1 - 0.9 * math.cos(2 * math.pi * s)
vc = speed(1.0, th, dth)
say(f'  Control (staircased interpolation through the same rungs): rung-locked wobble A = {stair_index(vc):+.3f}; '
    f'speed has {int(np.sum((vc[1:-1] < vc[:-2]) & (vc[1:-1] < vc[2:])))} local minima')

# Fair null: random valid interpolations theta(s) = s + periodic wiggle (all satisfy E^{s+1} = E o E^s)
rng = np.random.default_rng(7); idx = []
for _ in range(40):
    a1, b1, a2, b2 = rng.normal(0, 0.012, 4)
    t = lambda s: s + a1 * math.sin(2*math.pi*s) + b1 * (math.cos(2*math.pi*s) - 1) + a2 * math.sin(4*math.pi*s) + b2 * (math.cos(4*math.pi*s) - 1)
    dt = lambda s: 1 + 2*math.pi*(a1*math.cos(2*math.pi*s) - b1*math.sin(2*math.pi*s)) + 4*math.pi*(a2*math.cos(4*math.pi*s) - b2*math.sin(4*math.pi*s))
    idx.append(stair_index(speed(1.0, t, dt)))
idx = np.array(idx)
say(f'  Null: 40 random interpolations through the same rungs give rung-locked wobbles A from {idx.min():+.3f} to {idx.max():+.3f} '
    f'(spread {idx.std():.3f}); the regular one sits at {res[1.0]:+.5f}.')
say('  So a rung-locked pause is a free choice of interpolation, not something the rungs carry; the standard choice has none.')
say('  Reading: the standard (regular) flow between rungs never pauses at a rung. A pause would have to be put in by hand.')
open('fillspill/rungs_output.txt', 'w').write('\n'.join(out) + '\n')
json.dump({'S': S.tolist(), 'v1': speed(1.0).tolist(), 'vc': vc.tolist()}, open('fillspill/rungs_speed.json', 'w'))
