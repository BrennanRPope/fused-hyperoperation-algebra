# The shape's directions

Page: https://claude.ai/artifact/XJxHp4hRkWYBJ45ZWA7pqe

Question (Brennan): is there another whole dimension underneath, like q was?

## Findings
- One-input laws have exactly three global shapes, by number of periods of the log: 0 (addition), 1 (the circle operation ∘₁ (called mult-add before 1.0)), 2 (torus). Torsion counts 1, n, n² (checked n = 3, 5).
- Weight k is units: F_{λk}(a,b) = F_k(λa,λb)/λ. q_count = k+1 is the same dial.
- Landscape bend (p,q): rulers x/(1+cx) keep D = p²−4q, zoom scales D by λ², so only sign(D) survives (line / circle / addition). Same q as the twisted boosts' curvature.
- Rung height s: rate jets meet the plane on D = s to third order; fourth-order gap −s(s−1)/12, zero only at s = 0, 1. Every rung is addition in its own ruler: a route, not a shape.
- Nome: the one continuous direction (rank-2 lattice mod scale = 1 complex number, told apart by j). Invisible near 0 (gluing commutes).
- Doppler rapidity is a coordinate on one law; the boost twist is curvature × area (Helmholtz–Lie), no new dial.
- New directions: g waves give g(g+1)/2 period numbers (g nomes + pair weights), more than the 3g−3 curves from g = 4. On graphs the weight survives rescaling. Over whole numbers the weight's primes split additive (p | k) from multiplicative.

## Files
- jets.py → jets.json (rung jets; rate valid through a⁶)
- checks.py → checks_output.txt, results.json (29/29 pass)
- lean_mirror.py (sympy mirror of every Lean statement, 39/39)
- lean/ShapeDirections.lean (36 theorems, core Lean 4, no sorry, ✳ notation; not compiled here; intended path FusedAlgebra/Topics/ShapeDirections.lean, imports only Base)
- page/index.html (published page)
