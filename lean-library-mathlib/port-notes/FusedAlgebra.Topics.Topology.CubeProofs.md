# FusedAlgebra.Topics.Topology.CubeProofs — port remarks

- Almost everything in this module is in mathlib, so the port is short: `IsPrime` is `Nat.Prime`,
  `MEq p x y` is `x ≡ y [ZMOD p]` (`Int.ModEq`), `bsum f m` is `∑ k ∈ Finset.range m, f k`,
  `layer` is `Nat.choose`. Modules that used these helpers should use the mathlib names directly
  (see the table).
- `Int.ModEq` unfolds the other way round: `Int.modEq_iff_dvd : a ≡ b [ZMOD n] ↔ n ∣ b - a`,
  whereas `MEq p x y` was `p ∣ x - y`. `Int.ModEq.dvd` / `Int.emod_emod_of_dvd` and
  `Int.modEq_zero_iff_dvd` are usually what you want.
- The theorems take `(hp : p.Prime)` explicitly, as the originals did, rather than mathlib's
  `[Fact p.Prime]`.
- `blockCount`, `digitProd` and `cubeSum` are gone: the expressions are written inline (a product
  of two `Nat.choose`, a `Finset.prod` over the digits, a `Finset.sum`), matching mathlib's
  statements of Lucas's theorem and the binomial theorem.
- `fermat` is now stated for every integer `a` (the original was for `a : ℕ` cast to `ℤ`);
  `binomial` holds in any commutative semiring.
- New: `mpow_prime_modEq`, the mult-add reading of `frobenius` (`mpow x p ≡ x ^ p [ZMOD p]`).
- Only the converse half of `aks_middle` needed a real proof.
