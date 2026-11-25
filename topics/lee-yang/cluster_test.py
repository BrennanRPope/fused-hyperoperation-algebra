"""Cluster test: do landscape laws sit at the Beraha dial values k_n = 1 + 2cos(2pi/n)?

Beraha numbers B_n = 4cos^2(pi/n) are the colour counts Q where chromatic zeros pile up.
On the dial k = Q - 1 they are k_n = 1 + 2cos(2pi/n):
  n=2,3,4 -> -1, 0, 1 (these are the inputs the landscape was built on, so not a test)
  n=5..   -> 1.618, 2, 2.247, 2.414, 2.532, 2.618, ... -> 3   (the predictions)
Census = every law placed at a fixed (k, q) on the project's pages. Parameter families
(a whole line or curve of points) can land anywhere, so they are listed but not scored.
"""
import csv, math, random, itertools
from fractions import Fraction as F

# ---------------- census ----------------
# (law, k, q, page)  -- fixed placements only; duplicates across pages folded by name
C = []
def add(page, k, q, *names):
    for n in names: C.append((n, float(k), float(q), page))

add('landscape', 1, 0, 'Compound growth, returns, interest', 'Redshifts stacking', 'K-theory, Todd genus',
    'Claim-count models (PGF)', 'PUE overheads', 'Tsallis q=0 / Campbell cost')
add('landscape', 0, 0, 'Shannon and Renyi entropy, hazards, log-returns', 'Burgers shocks')
add('landscape', -1, 0, 'Union of failures, survival, decrements', 'Amplitude damping (T1)',
    'Gini-Simpson, Kraft slack, erasures', 'Slashing, validator churn')
add('landscape', -2, 0, 'Parity, bit flips, binary channel cascade', 'Bit-flip and phase-flip noise')
add('landscape', -1.5, 0, '3-letter channel')
add('landscape', -4/3, 0, '4-letter channel, depolarizing qubit')
add('landscape', -1.125, 0, 'Depolarizing qutrit')
add('landscape', 2, 0, 'Leveraged compounding, 2x')
add('landscape', -0.5, 0, 'Partial damage, half per hit')
add('landscape', 0, -1, 'Velocity addition, stacked reflections')
add('landscape', 0, 1, 'Slope addition (tangent)')
add('landscape', 0, -0.5, 'Sine addition (third order)')
add('landscape', 2, 1, 'Euler characteristic genus (cusp)')
add('landscape', 2, 0.5, 'Rung 2')
add('landscape', -1, 0.5, 'Rung -1')
add('landscape', 0.5, -0.0625, 'Half rung')
add('thermodynamics', 1, 0, 'Partition functions multiply', 'Third law', 'Unattainability', "Pauling's ice",
    "Landauer's bound", 'Equipartition', 'Fermion grand partition', 'Fridges in cascade',
    'Jarzynski equality', "Hess's law")
add('thermodynamics', 0, 0, 'Free energy, heat capacity add', 'Heat capacity from Var E')
add('thermodynamics', -1, 0, 'Engines in series, Carnot', 'Curzon-Ahlborn', 'Otto cycle, heat pumps',
    'Second law composes')
add('thermodynamics', -2, 0, 'Ising bonds in series')
add('thermodynamics', 0, -1, 'Fermi function', 'Negative temperature', 'Ising bonds in parallel',
    'Doppler temperature')
add('thermodynamics', -2, 1, "Stokes' pile, gray plates")
add('thermodynamics', 0, 0.5, 'Kaniadakis kappa-sum (kappa=1)')
add('thermodynamics', 0, 1, 'Fermi function after Wick rotation')
add('simple-laws', 0, 0, 'Adding counts', 'Relative errors add (first order)', 'Torques add', 'Buoyancy',
    'Galilean velocity addition', 'Equal areas', 'Lenses in contact', 'Mirror and thin-lens equation',
    'Springs side by side', 'Springs end to end', 'Forces and momenta add', 'Series resistors',
    'Parallel resistors')
add('simple-laws', 1, 0, 'Multiplying through +1', 'Non-empty choices', 'Relative error of a product',
    'Rectangle and gnomon', 'Box and cube volume', 'Pressure gains stacking')
add('simple-laws', -1, 0, 'Volume losses stacking', 'Law of cooling')
add('simple-laws', 0, 1, 'Pythagoras', "Euclid's triples", 'Tangent addition', 'Sine and cosine addition')
add('simple-laws', 0, -0.5, 'Sine addition on its own')
add('simple-laws', 0, -1, 'Telegraph line', 'Relativistic velocity addition')
add('statistics', 0, 0, 'Disjoint events add', 'Mean, variance, cumulants add', 'CLT shape numbers',
    'Log-odds add')
add('statistics', -1, 0, 'Independent union', 'Inclusion-exclusion')
add('statistics', -2, 0, 'Parity, XOR, piling-up')
add('statistics', 1, 0, 'Generating functions multiply', 'Moments of a sum', 'Bayes odds')
add('statistics', 0, -1, 'Pooling opinions')

# fold obvious cross-page duplicates of the same law
FOLD = {'Tangent addition': 'Slope addition (tangent)', 'Relativistic velocity addition':
        'Velocity addition, stacked reflections', 'Sine addition on its own': 'Sine addition (third order)',
        'Parity, XOR, piling-up': 'Parity, bit flips, binary channel cascade'}
seen, laws = set(), []
for n, k, q, p in C:
    key = FOLD.get(n, n)
    if key in seen: continue
    seen.add(key); laws.append((key, k, q, p))

FAMILIES = [('Tsallis entropy', 'k = 1 - q, the whole line', 'thermodynamics'),
            ('Lee-Yang three-spin chain', '(2a, 1)', 'thermodynamics'),
            ('Law of cosines', '(-2cos C, 1)', 'simple-laws'),
            ('Rotation in tilted axes', '(2cos theta, 1)', 'simple-laws'),
            ('Lenses a distance d apart', '(-d, 0)', 'simple-laws'),
            ('Changing units', '(lam k, lam^2 q)', 'simple-laws'),
            ('Dependent union', '(-lift, 0)', 'statistics'),
            ('AMH dependence', 'k + q = -1', 'statistics'),
            ('Elliptic curve over F_p', '(-a_p, p)', 'elliptic'),
            ('Noisy channels, d letters', '(-d/(d-1), 0) and (-d^2/(d^2-1), 0)', 'information'),
            ('Kaniadakis kappa-sum', '(0, kappa^2/2)', 'thermodynamics')]

with open('cluster/landscape_census.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['law', 'k', 'q', 'page'])
    for r in laws: w.writerow(r)
    for n, pt, p in FAMILIES: w.writerow([n, 'family', pt, p])

out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

say(f'Census: {len(laws)} fixed placements, {len(FAMILIES)} parameter families (not scored).')
vals = sorted({round(k, 6) for _, k, _, _ in laws})
say('Distinct dial values k:', vals)
from collections import Counter
cnt = Counter(round(k, 6) for _, k, _, _ in laws)
say('Laws per k:', dict(sorted(cnt.items())))

# ---------------- predictions ----------------
def beraha_k(n): return 1 + 2 * math.cos(2 * math.pi / n)
P = [beraha_k(n) for n in range(5, 13)]           # 8 predicted values, 1.618 .. 2.732
say('\nPredicted dial values (n = 5..12):', [round(x, 4) for x in P])
DELTA = 0.01                                      # placements are exact numbers
def hits(pred, values):
    return [v for v in values if any(abs(v - p) <= DELTA for p in pred)]

lawk = [k for _, k, _, _ in laws]
H_laws = hits(P, lawk); H_vals = hits(P, vals)
say(f'Test A (window above 1, tolerance {DELTA}): {len(H_laws)} laws hit, {len(H_vals)} distinct value(s):',
    sorted(set(H_vals)))
inwin = [(n, k, q) for n, k, q, _ in laws if 1.5 < k <= 3]
say('  every law with 1.5 < k <= 3:', inwin)
nonsimple = [p for n, p in zip(range(5, 13), P) if n != 6]
say('  hits at the seven non-integer predictions (n = 5, 7..12):', len(hits(nonsimple, lawk)))

# Null A1: random predicted sets of the same size in the same span
random.seed(1)
N = 200_000
ge_laws = ge_vals = 0
for _ in range(N):
    R = [random.uniform(1.5, 3.0) for _ in range(len(P))]
    if len(hits(R, lawk)) >= len(H_laws): ge_laws += 1
    if len(hits(R, vals)) >= len(H_vals): ge_vals += 1
say(f'  Null A1, random 8-point sets in [1.5, 3]: P(>= {len(H_laws)} law hits) = {ge_laws/N:.3f}, '
    f'P(>= {len(H_vals)} value hits) = {ge_vals/N:.3f}')

# Null A2: "simple numbers" -- fractions with denominator <= 4 in [1.5, 3]
S = sorted({F(a, b) for b in range(1, 5) for a in range(1, 13) if F(3, 2) <= F(a, b) <= 3})
Sf = [float(x) for x in S]
say(f'  Null A2, simple numbers (denominator <= 4) in [1.5, 3], {len(S)} values:', [str(x) for x in S])
best = max(len(hits(list(c), lawk)) for c in itertools.combinations(Sf, len(P)))
frac = sum(len(hits(list(c), lawk)) >= len(H_laws) for c in itertools.combinations(Sf, len(P))) / \
       math.comb(len(Sf), len(P))
say(f'    8-point subsets of them hit up to {best} laws; share with >= {len(H_laws)} hits: {frac:.3f}')
say('    Every hit is at k = 2, which the simple-number rule predicts as surely as Beraha does.')

# Test B: whole dial, Beraha incl. n = 2,3,4 vs. small integers
PB = [beraha_k(n) for n in range(2, 13)]
Z = [-2, -1, 0, 1, 2]
say('\nTest B (whole dial, not independent of how the landscape was built):')
say(f'  Beraha n = 2..12 ({len(PB)} values) hit {len(hits(PB, lawk))} of {len(lawk)} laws, '
    f'{len(hits(PB, vals))} of {len(vals)} distinct values')
say(f'  Integers -2..2 ({len(Z)} values) hit {len(hits(Z, lawk))} laws, {len(hits(Z, vals))} distinct values')
miss = [(n, k, q) for n, k, q, _ in laws if abs(k + 2) < 1e-9]
say(f'  k = -2 (colour count Q = -1, where rings never vanish) holds {len(miss)} laws, more than k = 2 '
    f'({cnt[2.0]}):', [m[0] for m in miss])

# Likelihood comparison: each law's k drawn from a spike mixture on a predicted set plus a uniform background
def loglik(pred, eps=0.2, lo=-3.0, hi=3.0):
    L = 0.0
    for v in lawk:
        on = any(abs(v - p) <= DELTA for p in pred)
        dens = (1 - eps) * (on / len(pred)) / (2 * DELTA) + eps / (hi - lo)
        L += math.log(dens)
    return L
say(f'  Log-likelihood, spike-plus-background model: Beraha n=2..12 {loglik(PB):.1f}, integers -2..2 {loglik(Z):.1f}'
    f'  (difference {loglik(Z) - loglik(PB):.1f} in favour of integers)')

# Same-size comparison: Beraha n = 2..6 (5 values) against the integers -2..2 (5 values)
PB5 = [beraha_k(n) for n in range(2, 7)]
say(f'  Same size, 5 values each: Beraha n=2..6 {[round(x,3) for x in PB5]} hits {len(hits(PB5, lawk))} laws; '
    f'integers -2..2 hit {len(hits(Z, lawk))}. Log-likelihoods {loglik(PB5):.1f} vs {loglik(Z):.1f}.')

# Test C: invariance. A change of units a -> lam*a slides (k, q) to (lam k, lam^2 q).
say('\nTest C (units): on the q = 0 axis every k > 0 is k = 1 in other units, so only sign(k) and k = 0')
say('  are properties of the law. Off the axis the unit-free number is k^2/q:')
for n, k, q, _ in laws:
    if abs(k - 2) < 1e-9:
        inv = 'sign +' if q == 0 else f'k^2/q = {k*k/q:g}'
        say(f'    {n}: (k, q) = ({k:g}, {q:g}) -> {inv}; same law at ({k/2:g}, {q/4:g}) in units of 2')
say('  The three hits at k = 2 all move to k = 1 under a = 2a\'. None is pinned to 2.')
open('cluster/cluster_output.txt', 'w').write('\n'.join(out) + '\n')
