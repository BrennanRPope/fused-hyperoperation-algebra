# FusedAlgebra.Topics.Topology.KummerWilson: review

Theorems (30): SAME 9, GENERAL 7, EQUIV 2, MATHLIB 12. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked:
- `val` = mathlib `multiplicity`. They agree for all p and x: both are 0 at x = 0 and p ≤ 1.
- `dsum p x` = `(p.digits x).sum`. They agree for p ≥ 2, and every use has p ≥ 2 or p prime.
- `carries`: identical recursive definition.
- `Val p x e` is still written out in `kummer`.
- `S` became `signedSum`: the same sum, over any commutative ring.
- Small checks have the same numbers: 6! = 720; carries 3 5 4 0 = 2; C(9,4) = 126; v₃(126) = 2.

Statements that match mathlib, as the user asked:
- `wilson`: `((p-1)! : ℤ) ≡ -1 [ZMOD p]`, the same as the original `MEq p (fact (p-1)) (-1)`. The proof is now `ZMod.wilsons_lemma`.
- `wilson_iff`: the same iff, for 2 ≤ n. The proof is now `Nat.prime_iff_fac_equiv_neg_one`.
- `kummer`: the same two-part statement (p^c divides C(j+k,k), p^(c+1) does not).
- `layer_val_dsum` became `sub_one_mul_multiplicity_choose_add_sum_digits`: the same equation, with no truncated subtraction.
- `fermat_unit`: `Int.ModEq.pow_card_sub_one_eq_one` (IsCoprime n p, n : ℤ) implies the original (a : ℕ, ¬p ∣ a).

Low/medium confidence:
- `Val.mul` (low): the table names `Nat.multiplicity_mul`, and I am not sure a declaration has that name. The content is mathlib's `multiplicity_mul (hp : Prime p) (hfin : FiniteMultiplicity p (a*b))`. Settle by `#check Nat.multiplicity_mul` / `#check multiplicity_mul`. The port does not use it.
- `bsum_meq` (medium): `Int.ModEq.sum` is believed to exist. Settle by `#check Int.ModEq.sum`. The port does not use it.

Proof-route changes (important to you):
- Wilson: your cube route (inclusion–exclusion `S_pow` on the (p−1)-cube + `alt_layer` + `fermat_unit`) is no longer the proof of `wilson`. The pieces survive: `signedSum_pow` (now derived from mathlib's forward-difference lemmas, not your induction), `neg_one_pow_mul_choose_pred_prime` (your induction kept), and `fermat_unit` (now mathlib). Nothing in the port assembles them. The docstring describes the route but does not prove it.
- `wilson_iff`: your converse via `dvd_fact` is replaced by mathlib.
- Kummer: your proof of `layer_val_dsum` (absorption + `val_spec` + `Val.mul` + `dsum_succ`) is replaced by mathlib's `sub_one_mul_multiplicity_choose_eq_sub_sum_digits'`. Your carries/digit-sum identity (`carries_dsum` → `sub_one_mul_carries_add_sum_digits`, any base p ≥ 2) is still proved by your own induction, and it is what links Kummer to `carries`. `dsum_succ` (`sum_digits_succ`) is ported but no longer feeds Kummer.
- `S_pow`: your induction is replaced by `fwdDiff_iter_eq_factorial` / `fwdDiff_iter_pow_eq_zero_of_lt`. The generalisation is to any commutative ring.

New port theorems with no original: 2 (`multiplicity_choose_eq_carries`, `carries_eq_card`).
