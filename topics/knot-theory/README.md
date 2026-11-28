# Knot theory through the circle operation (Nov 25, 2026)

> **Names:** some words or symbols on this page (wall, twist) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Files
- tangles.py: tangle diagrams (crossing tiles + joins), Kauffman bracket by the state sum over the
  cube, writhe by walking each component, and the (f, g) tangle algebra.
- checks.py -> checks_output.txt: 63/63 pass. Also writes jones.json (Jones polynomials of the
  seven rational knots, used by the page).
- lean/KnotTheory.lean: 45 theorems, compiles clean (see lean/README.md).
- index.html: the page (inline SVG drawn by script).

Findings
- Tangle sum is the dial: in x = g/f, x(S+T) = x(S) ∘_d x(T) with weight d = −A² − A⁻² (the loop
  value). The corner term d·g₁g₂ is the loop trapped when both sides are "down". Lift 1 + d·x =
  ⟨D(T)⟩/f is closing the tangle, and the lift law is ⟨D(S+T)⟩ = ⟨D(S)⟩⟨D(T)⟩ (connected sum).
  Ingredients classical (Kauffman 1987, Goldman–Kauffman 1997); the dial reading is ours.
- In y = A²g/f a crossing is y ↦ qy + 1 with q = −A⁻⁴ = −t, so [n] has y = [n]_q. The step turns
  the plane about the wall 1/(1 − q); twists sit on a regular polygon; walls merge at A⁴ = −1
  where d = 0, q = 1 and y is Conway's fraction (additive). Stacking is 1/y ↦ q⁻¹(1/y) + 1.
  These two steps are Morier-Genoud–Ovsienko's q-deformed SL(2,Z) (known, 2020).
- Rational knots 3_1 … 6_3 built from words: Jones matches the Knot Atlas; |V(−1)| = numerator.
- Landscape bend does not appear: 2-tangles only reach the dial line.
- Torus knots: Δ_T(p,q) = [q]_(t^p)/[q]_t; T(2,n): [n]_(−t). "Circle-operation value" P − 1 under
  connected sum is an honest relabelling.
- Correction to the q page: the twist pins q = −t, so the counting point is t = −1 (determinant).

Verdict: not new knot theory; a clean lens (trapped loop = corner, closing = lift, Conway's
fraction = dial at weight 0). Open question only: do 6-ended tangles (5-dim bracket space) reach
a landscape bend?
