# CategoryTheory.lean

- 36 theorems (plus the explicit isomorphisms `optionThese`, `optionTheseK`, `unitL`, `assoc` as definitions
  carrying their round-trip proofs), namespace `CategoryTheory`, core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Imports only `FusedAlgebra.Base`.
- Compiles clean: `LEAN_PATH=/tmp/lib/build lean lean/CategoryTheory.lean` (Nov 25, 2026). A `#print axioms`
  run on pair_unique, allThese_complete, align_fst, align_zip, reduced_product, circle_wall, sphere_rchi,
  torus_cells, s2_s2, assoc, unitL showed only propext, Classical.choice, Quot.sound (several use none).
- Intended module path: `FusedAlgebra/Topics/CategoryTheory.lean` (add to build.sh and FusedAlgebra.lean).

Contents: `These α β` (A ✳ B); Option(These α β) ≅ Option α × Option β with both round trips; Option Empty
is a point; functor laws; naturality; unit (These Empty β ≅ β); symmetry; associativity (7 corners each way);
tagged dial `TheseK` with 1 + K×(A ✳ₖ B) ≅ (1 + K×A) × (1 + K×B); Par: Kleisli composition, partial
projections, pairing, the product's universal property (fst_pair, snd_pair, pair_unique), Empty as zero object;
counting (all values listed, length = a ✳ b, completeness); wedge + smash; align: length = max, both inputs
recovered, zip = the both-part; cell polynomials: evaluation is multiplicative under products, reduced count
combines by the circle operation ∘₁ (called mult-add before 1.0) at every t, spheres at (−1)ⁿ, circle is the wall, torus cells 1-2-1, χ̃(S²×S²) = 3, suspension split count.
