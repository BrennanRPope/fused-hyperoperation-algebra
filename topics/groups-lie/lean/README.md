# GroupsLie.lean

- 52 theorems (50 about the mathematics plus 2 record-extension helpers, `M2.ext'` and `Fld.ext'`).
- Compiles clean with core Lean 4.37.0-rc1 against the shared library build (`/tmp/lib/build`):
  no errors, no warnings, no `sorry`. `#print axioms` on mob_mul, circle7_closed, f13_cyclic,
  Bl_invariant, mersenne and T_dilation showed only propext, Classical.choice, Quot.sound
  (the `#print` lines were removed afterwards).
- Imports only `FusedAlgebra.Base`. Namespace `GroupsLie`. Defines records `M2` (2×2 Int matrices)
  and `Fld` (fields α + βx + γx²), both deriving DecidableEq.
- Intended module path: `FusedAlgebra/Topics/GroupsLie.lean`
  (add `FusedAlgebra.Topics.GroupsLie` to build.sh and an import line to FusedAlgebra.lean).

Compile:

    export PATH=/opt/lean-4.37.0-rc1-linux/bin:$PATH
    LEAN_PATH=/tmp/lib/build lean GroupsLie.lean

Sections: (1) landscape law as a Möbius product, mob a = I + aX, X² = pX − qI;
(2) the θ-lift (1 + aθ)(1 + bθ), multiplicative norm, rate law, invariant ruler;
(3) Lie's fields α + βx + γx²: bracket, Jacobi, ad-invariant form β² − 4αγ;
(4) chart changes keep p² − 4d; (5) the dial in the ax+b group; (6) the matrix circle-operation commutator
= bracket; (7) finite clocks mod 5, 7, 11, 13 by `decide`, Mersenne powers, circle counts.
