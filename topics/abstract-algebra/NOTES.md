# Abstract algebra through the circle operation (Nov 25, 2026)

Page: "Circle-Operation Algebra Atlas" (claude.ai artifact DTBbDjwUepJAUmU9AQEAhV).

## For the library thread
New file to fold in: `lean/AbstractAlgebra.lean` -> `lean-library/FusedAlgebra/Topics/AbstractAlgebra.lean`.
Add `FusedAlgebra.Topics.AbstractAlgebra` to build.sh (before `FusedAlgebra Check`) and an import line to
FusedAlgebra.lean. Imports only FusedAlgebra.Base. Tested: full library + this file builds with core Lean
4.37.0-rc1, Check passes (653 theorems then; 53 written here plus 9 auto record-field lemmas).
Namespace `AbstractAlgebra`. Defines records CommMonoidLaws, CommGroupLaws, OpIso; defs fadd, geom, tateRel,
tateSetoid (a def, not an instance), TateQ, shiftIso, scaleIso.

## Proved (Lean) vs checked (checks.py, 33/33 pass)
Proved: monoid laws; shift iso (R,∘₁,0)≅(R,·,1) with wall->0; units = quasi-regular; Int units {0,-2};
field units off wall; scale iso for k≠0; Int scale embedding; weight −2 = XOR; nilpotent weight ⇒ group;
Fin 8/9/12 exact group criterion; Fin 12 units; Klein (Fin4,∘₂) not iso (Fin4,+); Fin 9 k=3 and Fin 8 k=4
cyclic, Fin 8 k=2 not; fused ring (fadd, ∘₁); mpow endomorphisms; Tate congruence + quotient monoid +
inverses; slope-addition walls mod 3,5,7,11,13.
Checked only: criterion for all n≤40; principal-units identification n≤30; ∘₁-irreducibles = p−1; finite
fields cyclic; Bennett commutative rung; ftet not distributive over fexp; landscape counts p−(Δ/p).
