# FusedAlgebra.Topics.Topology.KummerWilson — port remarks

- `val p x` is mathlib's ℕ-valued `multiplicity p x`, `dsum p x` is `(p.digits x).sum`, `fact` is
  `Nat.factorial`, `Val p x e` is the statement `multiplicity p x = e`.
- `carries` (recursive, with a carry-in) is kept as a definition: mathlib states Kummer's theorem
  with the count `#{i ∈ Ico 1 b | p ^ i ≤ k % p ^ i + n % p ^ i}` instead. The new lemma
  `carries_eq_card` shows the two agree for a prime base, and `multiplicity_choose_eq_carries`
  gives Kummer's theorem as an equation.
- The digit-sum identities `sub_one_mul_carries_add_sum_digits` (was `carries_dsum`) and
  `sum_digits_succ` (was `dsum_succ`) hold in every base `p ≥ 2`, not only primes, as in the
  original; mathlib only has the prime case, so their proofs are ported rather than replaced.
- `multiplicity_mul_self_left` (`multiplicity p (p * m) = multiplicity p m + 1` for `2 ≤ p`,
  `0 < m`, `p` not necessarily prime) is a general lemma that would belong next to mathlib's
  `multiplicity` API; it lives in this file.
- `wilson` and `wilson_iff` are now corollaries of `ZMod.wilsons_lemma` and
  `Nat.prime_iff_fac_equiv_neg_one`. The original derived Wilson from inclusion–exclusion on the
  cube; the ingredients of that route are still here (`signedSum_pow`,
  `neg_one_pow_mul_choose_pred_prime`, and `Int.ModEq.pow_card_sub_one_eq_one` for `fermat_unit`)
  but the route itself is not re-proved.
- `signedSum_pow` (was `S_pow`) is derived from mathlib's forward-difference lemmas
  (`fwdDiff_iter_eq_factorial`, `fwdDiff_iter_pow_eq_zero_of_lt`) and holds in any commutative
  ring.
