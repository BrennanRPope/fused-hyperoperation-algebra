# FusedAlgebra.Topics.Topology.HopfInvariantOne — port remarks

- Part 1 (the power of 2 in `3 ^ n - 1`) no longer uses the hand-built predicate `OneMod`; it is
  stated with mathlib's `multiplicity 2` and proved from the lifting-the-exponent lemma
  `Nat.two_pow_sub_pow`. The general helper lemmas `oneMod_oddPow` and `oneMod_sq` (for an
  arbitrary base `A`) are not restated: they are instances of mathlib's
  `emultiplicity_pow_sub_pow_of_prime` and `Int.two_pow_two_pow_sub_pow_two_pow`, and only the
  base `3` is needed.
- Part 3: `Trunc` (`α x + β x²` in `ℤ[x]/(x³)`) is now a `NonUnitalCommRing`, so the file no
  longer needs its own `madd` and `mpow`: mult-add is the library's `✳` and the mult-add power is
  the library's `MultAdd.mpow Trunc.x k` (Basic's `mpow` is defined for non-unital semirings, so
  this fits without change). `Trunc.multAdd_eq` records that `✳` is the old `madd` formula.
- `Elt` and `Trunc` are kept as structures over `ℤ` (they model K-theory classes); nothing was
  generalised here.
- The topology itself (existence and properties of Adams operations) is still only present as
  the hypotheses of `hopf_invariant_one`, exactly as in the original.
