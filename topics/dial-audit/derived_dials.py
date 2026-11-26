"""The derived-dial catalogue ("secondary colours"): each entry names a derived dial, the primaries it is built from,
how, why it is derived, and where it is checked. Read by figures.py (the wheel) and checks_derived.py.
The seven primaries are k, d, s, p, t, r, q (r and q are partly derived; see README section 6 and 7)."""

PRIMARIES = ['k', 'd', 's', 'p', 't', 'r', 'q']
PRIMARY_NAMES = {'k': 'homothety', 'd': 'curvature dial', 's': 'height', 'p': 'roundness', 't': 'odds',
                 'r': 'similarity ratio', 'q': 'order'}

# name, symbol, parents, recipe, why derived, where checked
DERIVED = [
    ('curvature', '−Δ/4', ['k', 'd'], 'd − k²/4', 'a function of k and d; adds no law of its own', 'affine-names DL1; dial-audit A7'),
    ('discriminant', 'Δ', ['k', 'd'], 'k² − 4d; its sign sorts the laws into three kinds', 'a function of k and d', 'landscape; dial-audit A7'),
    ('shape number', 'd/k²', ['k', 'd'], 'the one number of the (k, d) plane a unit change keeps', 'k and d with the unit divided out', 'dial-audit A6, Lean ratio_unit_free'),
    ('cubic form', 'C', ['k', 'd'], 'zero on every (k, d) law (their unit curves are conics)', 'fixed at 0 by k and d (Brennan, Nov 26); nonzero only as a reading of shapes p reaches', 'dial-map note; affine-names §4'),
    ('Tsallis q', 'q_T', ['k'], '1 − k (Tsallis 1988; Borges 2004)', 'a relabelling of k', 'dial-survey finding 2'),
    ('copies allowed', 'm', ['k', 't'], 'k = −1/m is m one-copy slots, each at odds t/m', 'k read through the odds', 'dial-audit D4; dial-survey finding 4'),
    ('rung reading', 's·k', ['k', 's'], "the homothety reading of rung s of ∘ₖ", 'the first generator coefficient of a rung law', 'dial-audit B1'),
    ('face power', 'p on faces', ['p'], "p applied to a polytope's face normals (the family is a discrete choice)", 'the ℓp sum taken over face normals', 'shape-dials (settled Nov 26)'),
    ('Wythoff depth', 'depth', ['t'], "odds on a family's corner shapes, joined by a Minkowski sum", 'a weighted Minkowski sum of shapes p and the family already give', 'shape-dials §6, depth_checks 20/20'),
    ('chance', 't/(1+t)', ['t'], 'the probability reading of the odds', 'a function of t', 'affine-names §6'),
    ('prime dial', 't_p = p^−σ', ['t'], 'odds on the primes, set by one exponent', 'per-element odds on primes', 'dial-survey finding 1'),
    ('dimension', 'D', ['r'], 'solves Moran’s Σ rᵢᴰ = 1', 'a reading of r', 'ratio-dial-survey; dial-audit E1'),
    ('piece weights', 'wᵢ = rᵢᴰ', ['r', 't'], 'relative odds between pieces (which piece the natural measure picks)', 'r with D divided out; the categorical cousin of t', 'dial-audit E1–E4'),
    ('trial weight', 'rˢ at fixed s', ['r', 'k'], 'acts as the homothety k at a fixed trial exponent (not at D)', 'k in disguise', 'ratio-dial-survey C3'),
    ('growth degree', 'β', ['r'], 'doubling weight 2^−β is ratio ½ read at D = β', 'a setting read through r', 'affine-names §6 R4'),
    ('holonomy', 'ω', ['d', 'q'], 'd × the area two moves span (Gauss–Bonnet)', 'a reading of d and two directions', 'twisted-boosts; dial-audit F3, F4'),
    ('counting q', '[n]_q', ['q', 'k'], 'the n-fold ∘ₖ power of 1 at k = q − 1 (Jackson 1910)', 'k in disguise', 'dial-audit F1, F2, Lean qint_is_circ_power'),
    ('face weight', 'β/α', ['q', 'd'], 'the order survey’s face weight q = β/α', 'a reading of the curvature dial', 'order-dial-survey'),
]
