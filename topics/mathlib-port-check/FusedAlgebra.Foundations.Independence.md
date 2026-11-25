# FusedAlgebra.Foundations.Independence

Theorems (7): EQUIV 4, MATHLIB 3. No WEAKER / DIFFERENT / MISSING / VACUOUS.

The four independence theorems went from 8-part conjunctions to 4-part ones. A1–A4 are now the AddCommGroup instance of each model (ℤ, Matrix (Fin 2) (Fin 2) ℤ, ℤ×ℤ×ℤ), so nothing is lost. Same models, same walls (0, every w, -I, (-1,0,0)) and same counterexamples.
- M4 changed from `(m a b + m a c) + -a` to `m a b + m a c - a`, which are equal by sub_eq_add_neg. This matches the port's MARing.madd_add field, so the independence is relative to the same axiom list.
- Tri.mul is identical componentwise (p, q, r map to .1, .2.1, .2.2). Mat.mul is the ordinary 2×2 product.
- wAdd 2 a b = a + b + 2*a*b exactly.

Route: in no_M1 and no_M3 the holding axioms now come from new general lemmas (m2/m3/m4_multAdd for any ring, wAdd_axioms), not from componentwise computation. This is minor.

New port theorems with no original: 5 (axioms_of_maRing, m2_multAdd, m3_multAdd, m4_multAdd, wAdd_axioms), plus the abbrev maRingOfAxioms.
