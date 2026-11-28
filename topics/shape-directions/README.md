# The shape's directions

> **Names:** some words or symbols on this page (mult-add, rung, twist) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Page: https://claude.ai/artifact/XJxHp4hRkWYBJ45ZWA7pqe

Question (Brennan): is there another whole dimension underneath, like d was?

## Findings
- One-input laws have exactly three global shapes, by number of periods of the log: 0 (addition), 1 (the circle operation ∘₁ (called mult-add before 1.0)), 2 (torus). Torsion counts 1, n, n² (checked n = 3, 5).
- Weight k is units: F_{λk}(a,b) = F_k(λa,λb)/λ. The colour count k + 1 is the same dial.
- Landscape bend (p,d): rulers x/(1+cx) keep D = p²−4d, zoom scales D by λ², so only sign(D) survives (line / circle / addition). Same d as the twisted boosts' curvature.
- Rung height s: rate jets meet the plane on D = s to third order; fourth-order gap −s(s−1)/12, zero only at s = 0, 1. Every rung is addition in its own ruler: a route, not a shape.
- Nome: the one continuous direction (rank-2 lattice mod scale = 1 complex number, told apart by j). Invisible near 0 (gluing commutes).
- Doppler rapidity is a coordinate on one law; the boost twist is curvature × area (Helmholtz–Lie), no new dial.
- New directions: g waves give g(g+1)/2 period numbers (g nomes + pair weights), more than the 3g−3 curves from g = 4. On graphs the weight survives rescaling. Over whole numbers the weight's primes split additive (p | k) from multiplicative.

## Files
- jets.py → jets.json (rung jets; rate valid through a⁶)
- checks.py → checks_output.txt, results.json (29/29 pass)
- lean_mirror.py (sympy mirror of every Lean statement, 39/39)
- lean/ShapeDirections.lean (36 theorems, core Lean 4, no sorry, ∘₁ notation; not compiled here; intended path FusedAlgebra/Topics/ShapeDirections.lean, imports only Base)
- page/index.html (published page)

## Several waves, followed (second pass)
- KP fixes each wave by two points (p, q) on a line; Hirota's two-wave condition forces the pair weight to be the cross ratio (p₁−p₂)(q₁−q₂)/((p₁−q₂)(q₁−p₂)). Verified: tau with these weights solves KP for g = 1..4 (sympy exact), and perturbing a weight breaks it.
- A = 1 exactly when a wave is trivial (p = q): two real waves always interact. A = 0 when they share a point (resonance). Spreading waves apart sends every A to 1: g separate circle operations.
- Count: g(g−1)/2 weights vs 2g−3 numbers from points (exact Jacobian rank, g = 2..8). Excess (g−2)(g−3)/2, the same as Riemann period numbers over curves. KdV (q = −p): g−1 numbers, excess from three waves.
- Four waves: one explicit rule (waves/relation4.txt), irreducible, S4-symmetric, degree 2 in each weight, 393 terms; vanishes identically on cross ratios (relcheck2.py). Five weights fix the sixth up to two values. All-2s array breaks it.
- Lean: lean/SeveralWaves.lean (23 theorems; mirrored in waves/mirror_check.py). Intended path FusedAlgebra/Topics/SeveralWaves.lean.
- Files: waves/rank.py (rank.json), kp.py (kp.json), sixth.py, interp4.py → relation4.txt, relcheck2.py (relcheck.json), mirror_check.py, tolean.py.

## Rule breakers (third pass)
- Graph-shaped arrays (weight on a graph's pairs, 1 elsewhere): every graph without a 4-cycle obeys the four-wave rule for all weights (path, star, triangle, paw and their subgraphs). The 4-cycle gives rule = ∏(w−1)²; 4-cycle + diagonal gives k(k−1)⁸; K₄ gives −k³(k−1)⁸(k+2). Breaking the (closed) rule means not even a limit of point-waves.
- One weight k on every pair: two waves any k; three waves any k ≠ 0, 1 (third wave at a root of t² − (k+1)t + 1; real iff k > 1 or k ≤ −3); four waves only k = −2 (roots forced to be cube roots of unity: waves (0,∞),(1,−1/2),(ω,−ω/2),(ω²,−ω²/2), the tetrahedron's vertex–antipode axes, complex only); five or more waves never.
- Lean: lean/RuleBreakers.lean (14 theorems, imports SeveralWaves); mirror waves/breakers_mirror.py (14/14). Checks: waves/breakers.py, breakers2.py, breakers_graphs.json.
