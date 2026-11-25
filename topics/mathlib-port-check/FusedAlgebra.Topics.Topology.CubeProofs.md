# FusedAlgebra.Topics.Topology.CubeProofs: review

Theorems (23): SAME 4, GENERAL 2, EQUIV 1, MATHLIB 15, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked: `IsPrime` = `Nat.Prime` (`Nat.prime_def_lt''`, exact match). `MEq p x y` (`p ∣ x - y`) = `x ≡ y [ZMOD p]`: the relation is symmetric, so this is the same; only `Int.modEq_iff_dvd` unfolds as `p ∣ y - x`. `bsum` = `Finset.sum` over `range`. `layer` = `Nat.choose` (same Pascal recursion). `digitProd` unfolds to the port's `∏ i ∈ range f, (n / p^i % p).choose (k / p^i % p)`. The small checks use the same numbers (C(9,3)=84; C(10,4)=210 with digit product C(3,4)·C(1,0)=0).

Table claims that are loose (no finding on statements):
- `blockCount_eq` is listed as covered by `lucas_step`. `lucas_step` doesn't state it. It was a helper of the hand-made Lucas proof, and its content is `Nat.mul_add_div` / `Nat.mul_add_mod` / `Nat.div_eq_of_lt`. Marked DROPPED-OK, medium.
- `divmod_of_eq` has no single mathlib lemma. It is covered by a combination of four lemmas (medium).

Proof-route changes (your own proofs replaced by mathlib):
- `prime_dvd_middle`: was absorption + coprimality; now `Nat.Prime.dvd_choose_self`.
- `lucas_step`, `lucas`: the hand induction (carry cases on the cube) is gone; now `Choose.choose_modEq_choose_mod_mul_choose_div` and `Choose.lucas_theorem`.
- `binomial`: was a hand induction; now `add_pow`.
- `frobenius`: was binomial + middle layers vanish mod p; now `add_pow_char` in `ZMod p`.
- `fermat`: was induction a → a+1 through `frobenius` (the "one mult-add unit per step" reading); now `ZMod.pow_card`.
- `aks_middle`: the forward half is still proved here, along the same lines (absorption + one Lucas step).

New port theorems with no original: 1 (`mpow_prime_modEq`).
