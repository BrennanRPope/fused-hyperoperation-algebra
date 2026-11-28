# Category theory through the circle operation (Nov 25, 2026)

> **Names:** some words or symbols on this page (mult-add, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Page title: "Give each side a 'nothing' and the circle operation ∘₁ (called mult-add before 1.0) becomes a plain product" (source index.html, not yet published).

## Files
- index.html: the page (Option grid with wedge/smash, These/align/outer join in code, Par product diagram,
  the cube as associativity, reduced Euler characteristic line with the circle as wall, known/ours, verdict).
- checks.py -> checks_output.txt: 29/29 pass (enumeration of bijections, pentagon/hexagon, Par universal property
  by brute force, align vs zip_longest, cell polynomials, simplicial homology of a torus and its suspension).
- lean/CategoryTheory.lean: 36 theorems, compiles clean (see lean/README.md).

## Findings (all textbook unless marked)
1. A ∘₁ B = A + B + A×B is Haskell's `These` / Rust's `EitherOrBoth`; Option(A ∘₁ B) ≅ Option A × Option B
   (the lift as sets). A₊ × B₊ minus the base point = wedge (a+b) + smash (ab).
2. A ∘₁ B is the categorical product in Par (sets and partial functions); ∅ is a zero object and the unit;
   Par ≃ pointed sets carries ∘₁ to ×. (FinSet, ∘₁, ∅) is symmetric monoidal (pentagon, hexagon checked).
3. align/zip_longest returns These; both inputs come back, zip is the ab corner. Full outer joins and diffs are These-valued.
4. Σ(X×Y) ≃ ΣX ∨ ΣY ∨ Σ(X∧Y) (standard; topology page did not cover it; homology of ΣT² checked).
5. Ours as a reading: χ̃ = χ − 1 is the circle-operation coordinate (χ̃(X×Y) = χ̃X ∘₁ χ̃Y, at every t for cell
   polynomials); the circle is the wall, spheres sit at ±1, a point at 0.
6. Ours, small: tagged Option 1 + K×A lifts the dial A + B + k·AB for k ≥ 0; negative k has no set model.

## Verdict
No new mathematics; a firm underpinning (the circle operation is a real categorical product, the lift a functor, the wall the
circle). A practical lead for software: treat merge/sync/diff as functions out of These and hold them to the
monoidal laws.
