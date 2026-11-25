module

public import FusedAlgebra.Topics.Topology.CubeProofs

/-!
# Kummer and Wilson on the mult-add cube

Two more classical facts read off the layer counts `n.choose k` of the cube.

## Kummer: carries

* `MultAdd.Cube.carries p a b c`: the number of carries when `a` and `b` are added in base `p`,
  starting with carry `c`.
* `MultAdd.Cube.sub_one_mul_carries_add_sum_digits`: every carry turns `p` units of digit sum
  into `1`, `(p - 1) * carries + s(a + b + c) = s(a) + s(b) + c` (any base `p ≥ 2`).
* `MultAdd.Cube.sum_digits_succ`: adding `1` trades each trailing digit `p - 1` for a carry.
* `MultAdd.Cube.sub_one_mul_multiplicity_choose_add_sum_digits`:
  `(p - 1) * v_p ((j + k).choose k) + s(j + k) = s(j) + s(k)`.
* `MultAdd.Cube.kummer`: the exact power of a prime `p` dividing `(j + k).choose k` is the number
  of carries when `j` and `k` are added in base `p` (Kummer 1852).

## Wilson: inclusion–exclusion on the cube

* `MultAdd.Cube.signedSum n f`: the signed corner sum `∑ k, (-1) ^ k * n.choose k * f k`.
* `MultAdd.Cube.signedSum_pow`: `∑ k, (-1) ^ k * n.choose k * (c - k) ^ m` is `n !` for `m = n`
  and `0` for `m < n`.
* `MultAdd.Cube.wilson`, `MultAdd.Cube.wilson_iff`: `(p - 1)! ≡ -1` mod `p` for a prime `p`
  (Lagrange 1771), and for `2 ≤ n` the congruence holds exactly when `n` is prime.

This is the mathlib form of `FusedAlgebra.Topics.Topology.KummerWilson`. Everything here is known
mathematics, and the headline results are in mathlib:

* Kummer's theorem is `multiplicity_choose'` (carries counted as a `Finset` cardinality) and
  `sub_one_mul_multiplicity_choose_eq_sub_sum_digits'` (the digit-sum form);
* the inclusion–exclusion identity is `fwdDiff_iter_eq_factorial` and
  `fwdDiff_iter_pow_eq_zero_of_lt` (iterated forward differences of `x ^ m`);
* Wilson's theorem is `ZMod.wilsons_lemma`, with converse `Nat.prime_iff_fac_equiv_neg_one`.

The core-Lean file built the `p`-adic valuation `val` and the digit sum `dsum` by hand; they are
now `multiplicity p x` and `(p.digits x).sum`. The recursive carry count `carries` (with its
carry-in argument) is kept, since mathlib has no such definition, and the digit-sum identity for
it is proved here in every base `p ≥ 2`.
-/

@[expose] public section

namespace MultAdd.Cube

open Finset Nat fwdDiff

/-! ### Kummer: carries -/

/-- `carries p a b c`: how many carries happen when `a` and `b` are added in base `p`, starting
with carry `c`. -/
def carries (p a b c : ℕ) : ℕ :=
  if _h : 2 ≤ p ∧ 0 < a + b then
    (a % p + b % p + c) / p + carries p (a / p) (b / p) ((a % p + b % p + c) / p)
  else 0
termination_by a + b
decreasing_by
  have h1 := Nat.div_le_self a p
  have h2 := Nat.div_le_self b p
  rcases Nat.eq_zero_or_pos a with ha | ha
  · have := Nat.div_lt_self (n := b) (by omega) (by omega : 1 < p); omega
  · have := Nat.div_lt_self ha (by omega : 1 < p); omega

theorem carries_eq {p : ℕ} (hp : 2 ≤ p) (a b c : ℕ) (hab : 0 < a + b) :
    carries p a b c = (a % p + b % p + c) / p +
      carries p (a / p) (b / p) ((a % p + b % p + c) / p) := by
  rw [carries]; simp [hp, hab]

@[simp]
theorem carries_zero (p c : ℕ) : carries p 0 0 c = 0 := by
  rw [carries]; simp

/-- The digit sum, one digit at a time. -/
theorem sum_digits_eq {p : ℕ} (hp : 2 ≤ p) (x : ℕ) :
    (p.digits x).sum = x % p + (p.digits (x / p)).sum := by
  rcases Nat.eq_zero_or_pos x with rfl | hx
  · simp
  · rw [Nat.digits_def' hp hx, List.sum_cons]

/-- Adding digit by digit: each carry turns `p` units of digit sum into `1`. -/
theorem sub_one_mul_carries_add_sum_digits {p : ℕ} (hp : 2 ≤ p) (a b c : ℕ) (hc : c ≤ 1) :
    (p - 1) * carries p a b c + (p.digits (a + b + c)).sum =
      (p.digits a).sum + (p.digits b).sum + c := by
  induction hs : a + b using Nat.strong_induction_on generalizing a b c with
  | _ s ih =>
  subst hs
  rcases Nat.eq_zero_or_pos (a + b) with h0 | hpos
  · obtain rfl : a = 0 := by omega
    obtain rfl : b = 0 := by omega
    rcases Nat.eq_zero_or_pos c with rfl | hc0
    · simp
    · rw [carries_zero, Nat.zero_add, Nat.digits_of_lt p c (by omega) (by omega)]
      simp
  · rw [carries_eq hp a b c hpos]
    have hp0 : 0 < p := by omega
    have ra := Nat.mod_add_div a p
    have rb := Nat.mod_add_div b p
    have ma := Nat.mod_lt a hp0
    have mb := Nat.mod_lt b hp0
    generalize hr : a % p + b % p + c = r
    have rr := Nat.mod_add_div r p
    have hc1 : r / p ≤ 1 := by
      have : r / p < 2 := Nat.div_lt_of_lt_mul (by omega)
      omega
    have hsum : a + b + c = p * (a / p + b / p) + r := by
      rw [Nat.mul_add]; omega
    have hdiv : (a + b + c) / p = a / p + b / p + r / p := by
      rw [hsum, Nat.mul_add_div hp0]
    have hmod : (a + b + c) % p = r % p := by
      rw [hsum, Nat.mul_add_mod]
    have hlt : a / p + b / p < a + b := by
      have h1 := Nat.div_le_self a p
      have h2 := Nat.div_le_self b p
      rcases Nat.eq_zero_or_pos a with ha | ha
      · have := Nat.div_lt_self (n := b) (by omega) (by omega : 1 < p); omega
      · have := Nat.div_lt_self ha (by omega : 1 < p); omega
    have IH := ih _ hlt (a / p) (b / p) (r / p) hc1 rfl
    rw [sum_digits_eq hp (a + b + c), hdiv, hmod, sum_digits_eq hp a, sum_digits_eq hp b]
    have e1 : (p - 1) * (r / p + carries p (a / p) (b / p) (r / p)) =
        (p - 1) * (r / p) + (p - 1) * carries p (a / p) (b / p) (r / p) := Nat.mul_add _ _ _
    have e2 : (p - 1) * (r / p) + r / p = p * (r / p) := by
      rw [Nat.sub_one_mul]
      have : r / p ≤ p * (r / p) := Nat.le_mul_of_pos_left _ hp0
      omega
    omega

/-- Taking out one factor `p` lowers the multiplicity by one. -/
theorem multiplicity_mul_self_left {p m : ℕ} (hp : 2 ≤ p) (hm : 0 < m) :
    multiplicity p (p * m) = multiplicity p m + 1 := by
  have hf : FiniteMultiplicity p m := Nat.finiteMultiplicity_iff.2 ⟨by omega, hm⟩
  refine multiplicity_eq_of_dvd_of_not_dvd ?_ ?_
  · rw [pow_succ']
    exact Nat.mul_dvd_mul_left p (pow_multiplicity_dvd p m)
  · rw [pow_succ', Nat.mul_dvd_mul_iff_left (by omega)]
    exact hf.not_pow_dvd_of_multiplicity_lt (Nat.lt_succ_self _)

/-- Adding `1`: the digit sum goes up by `1`, minus `p - 1` for each trailing digit `p - 1`. -/
theorem sum_digits_succ {p : ℕ} (hp : 2 ≤ p) (x : ℕ) :
    (p.digits (x + 1)).sum + (p - 1) * multiplicity p (x + 1) = (p.digits x).sum + 1 := by
  induction x using Nat.strong_induction_on with
  | _ x ih =>
  have hp0 : 0 < p := by omega
  have rx := Nat.mod_add_div x p
  have mx := Nat.mod_lt x hp0
  rcases Nat.lt_or_ge (x % p + 1) p with h | h
  · have e : x + 1 = p * (x / p) + (x % p + 1) := by omega
    have hdiv : (x + 1) / p = x / p := by
      rw [e, Nat.mul_add_div hp0, Nat.div_eq_of_lt h, Nat.add_zero]
    have hmod : (x + 1) % p = x % p + 1 := by
      rw [e, Nat.mul_add_mod, Nat.mod_eq_of_lt h]
    have hnd : ¬p ∣ x + 1 := fun hd ↦ by
      have := Nat.mod_eq_zero_of_dvd hd
      omega
    rw [multiplicity_eq_zero_of_not_dvd hnd, sum_digits_eq hp (x + 1), hdiv, hmod,
      sum_digits_eq hp x]
    omega
  · have e : x + 1 = p * (x / p + 1) := by rw [Nat.mul_add]; omega
    have hx : 0 < x := by
      rcases Nat.eq_zero_or_pos x with rfl | h0
      · simp at h; omega
      · exact h0
    have IH := ih (x / p) (Nat.div_lt_self hx (by omega))
    have hdiv : (x + 1) / p = x / p + 1 := by
      rw [e, Nat.mul_div_cancel_left _ hp0]
    have hmod : (x + 1) % p = 0 := by
      rw [e, Nat.mul_mod_right]
    rw [sum_digits_eq hp (x + 1), hdiv, hmod, sum_digits_eq hp x, e,
      multiplicity_mul_self_left hp (show 0 < x / p + 1 from Nat.succ_pos _), Nat.mul_add,
      Nat.mul_one]
    omega

/-- The cube count and the digit sums:
`(p - 1) * v_p ((j + k).choose k) + s(j + k) = s(j) + s(k)`. This is mathlib's
`sub_one_mul_multiplicity_choose_eq_sub_sum_digits'`, without the truncated subtraction. -/
theorem sub_one_mul_multiplicity_choose_add_sum_digits {p : ℕ} (hp : p.Prime) (k j : ℕ) :
    (p - 1) * multiplicity p ((j + k).choose k) + (p.digits (j + k)).sum =
      (p.digits j).sum + (p.digits k).sum := by
  have := Fact.mk hp
  have h1 := sub_one_mul_multiplicity_choose_eq_sub_sum_digits' (p := p) (k := k) (n := j)
  have h2 := sub_one_mul_carries_add_sum_digits hp.two_le j k 0 (Nat.zero_le _)
  rw [Nat.add_zero] at h2
  omega

/-- The multiplicity of a prime `p` in `(j + k).choose k` is the number of carries when `j` and
`k` are added in base `p`. -/
theorem multiplicity_choose_eq_carries {p : ℕ} (hp : p.Prime) (j k : ℕ) :
    multiplicity p ((j + k).choose k) = carries p j k 0 := by
  have h1 := sub_one_mul_multiplicity_choose_add_sum_digits hp k j
  have h2 := sub_one_mul_carries_add_sum_digits hp.two_le j k 0 (Nat.zero_le _)
  rw [Nat.add_zero] at h2
  have p2 := hp.two_le
  exact Nat.eq_of_mul_eq_mul_left (show 0 < p - 1 by omega) (by omega)

/-- **Kummer's theorem (1852).** The exact power of a prime `p` that divides `(j + k).choose k`,
the number of corners of the `(j + k)`-cube with `k` switches on, is the number of carries when
`j` and `k` are added in base `p`. Mathlib's form of it is `multiplicity_choose'`. -/
theorem kummer {p : ℕ} (hp : p.Prime) (j k : ℕ) :
    p ^ carries p j k 0 ∣ (j + k).choose k ∧ ¬p ^ (carries p j k 0 + 1) ∣ (j + k).choose k := by
  have hf : FiniteMultiplicity p ((j + k).choose k) :=
    Nat.finiteMultiplicity_iff.2 ⟨hp.ne_one, Nat.choose_pos (Nat.le_add_left k j)⟩
  rw [← multiplicity_choose_eq_carries hp]
  exact ⟨pow_multiplicity_dvd _ _, hf.not_pow_dvd_of_multiplicity_lt (Nat.lt_succ_self _)⟩

/-- The recursive carry count agrees with mathlib's count of the positions `i` at which a carry
happens. -/
theorem carries_eq_card {p : ℕ} (hp : p.Prime) (j k b : ℕ) (hb : Nat.log p (j + k) < b) :
    carries p j k 0 = #{i ∈ Ico 1 b | p ^ i ≤ k % p ^ i + j % p ^ i} := by
  have := Fact.mk hp
  rw [← multiplicity_choose_eq_carries hp, multiplicity_choose' hb]

/-! ### Wilson: inclusion–exclusion on the cube -/

section SignedSum

variable {R : Type*} [CommRing R]

/-- Signed corner sum over the `n`-cube: a corner with `k` switches on counts with sign
`(-1) ^ k` and weight `f k`. -/
def signedSum (n : ℕ) (f : ℕ → R) : R :=
  ∑ k ∈ range (n + 1), (-1) ^ k * (n.choose k : R) * f k

theorem signedSum_add (n : ℕ) (f g : ℕ → R) :
    signedSum n (fun k ↦ f k + g k) = signedSum n f + signedSum n g := by
  simp only [signedSum, mul_add, sum_add_distrib]

theorem signedSum_mul (n : ℕ) (c : R) (f : ℕ → R) :
    signedSum n (fun k ↦ c * f k) = c * signedSum n f := by
  simp only [signedSum, mul_sum]
  exact sum_congr rfl fun k _ ↦ by ring

@[simp]
theorem signedSum_zero (n : ℕ) : signedSum n (fun _ ↦ (0 : R)) = 0 := by
  simp [signedSum]

theorem signedSum_sub (n : ℕ) (f g : ℕ → R) :
    signedSum n (fun k ↦ f k - g k) = signedSum n f - signedSum n g := by
  simp only [signedSum, mul_sub, sum_sub_distrib]

theorem signedSum_sum (n : ℕ) (g : ℕ → ℕ → R) (m : ℕ) :
    signedSum n (fun k ↦ ∑ i ∈ range m, g i k) = ∑ i ∈ range m, signedSum n (g i) := by
  simp only [signedSum, mul_sum]
  exact sum_comm

/-- Pascal's rule for the signed sum: adding an axis takes a difference. -/
theorem signedSum_succ (n : ℕ) (f : ℕ → R) :
    signedSum (n + 1) f = signedSum n f - signedSum n (fun k ↦ f (k + 1)) := by
  have h : signedSum n f = ∑ k ∈ range (n + 1), (-1) ^ (k + 1) * (n.choose (k + 1) : R) * f (k + 1)
      + f 0 := by
    rw [signedSum, sum_range_succ', sum_range_succ, Nat.choose_succ_self]
    simp
  rw [h, signedSum, signedSum, sum_range_succ', sub_eq_add_neg, ← sum_neg_distrib, add_right_comm,
    ← sum_add_distrib]
  simp only [Nat.choose_zero_right, pow_zero, Nat.cast_one, mul_one, one_mul, add_left_inj]
  refine sum_congr rfl fun k _ ↦ ?_
  rw [Nat.choose_succ_succ, Nat.cast_add]
  ring

/-- **Inclusion–exclusion on the cube.** `∑ k, (-1) ^ k * n.choose k * (c - k) ^ m` is `n !` when
`m = n` and `0` when `m < n`. With `c = n` it counts the orderings of `n` things: corner `K` of
the `n`-cube is the set of values a map is made to miss, and the signed count of maps that miss
them leaves exactly the maps that miss nothing. In mathlib this is the `n`-th forward difference
of `x ^ m`: `fwdDiff_iter_eq_factorial` and `fwdDiff_iter_pow_eq_zero_of_lt`. -/
theorem signedSum_pow (n m : ℕ) (hm : m ≤ n) (c : R) :
    signedSum n (fun k ↦ (c - k) ^ m) = if m = n then (n ! : R) else 0 := by
  have e : signedSum n (fun k ↦ (c - k) ^ m) = (Δ_[1])^[n] (fun r : R ↦ r ^ m) (c - n) := by
    rw [fwdDiff_iter_eq_sum_shift, ← sum_range_reflect, signedSum]
    refine sum_congr rfl fun k hk ↦ ?_
    have hk' : k ≤ n := Nat.lt_succ_iff.1 (mem_range.1 hk)
    simp only [Nat.add_sub_cancel, Nat.sub_sub_self hk', Nat.choose_symm hk', Nat.cast_sub hk',
      zsmul_eq_mul, nsmul_eq_mul, mul_one]
    push_cast
    ring
  rw [e]
  split_ifs with h
  · subst h
    rw [fwdDiff_iter_eq_factorial]
    simp
  · rw [fwdDiff_iter_pow_eq_zero_of_lt (lt_of_le_of_ne hm h)]
    simp

end SignedSum

/-- Along the `(p - 1)`-cube, the layer counts alternate `±1` mod `p`. -/
theorem neg_one_pow_mul_choose_pred_prime {p : ℕ} (hp : p.Prime) (k : ℕ) (hk : k ≤ p - 1) :
    (-1 : ℤ) ^ k * ((p - 1).choose k : ℤ) ≡ 1 [ZMOD p] := by
  have := Fact.mk hp
  rw [← ZMod.intCast_eq_intCast_iff]
  push_cast
  obtain ⟨q, rfl⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.pos; omega⟩
  rw [Nat.add_sub_cancel] at hk ⊢
  induction k with
  | zero => simp
  | succ k ih =>
    have hd : (((q + 1).choose (k + 1) : ℕ) : ZMod (q + 1)) = 0 :=
      (ZMod.natCast_eq_zero_iff _ _).2 (hp.dvd_choose_self k.succ_ne_zero (by omega))
    rw [Nat.choose_succ_succ, Nat.cast_add] at hd
    have IH := ih (by omega)
    linear_combination IH - (-1) ^ k * hd

/-- **Wilson's theorem (Lagrange 1771).** `(p - 1)! ≡ -1` mod `p` for every prime `p`. This is
mathlib's `ZMod.wilsons_lemma`. On the cube: inclusion–exclusion on the `(p - 1)`-cube
(`signedSum_pow`) counts the `(p - 1)!` orderings; mod `p` every layer but the last adds exactly
`1` (`neg_one_pow_mul_choose_pred_prime` and Fermat), so the count is `p - 1 ≡ -1`. -/
theorem wilson {p : ℕ} (hp : p.Prime) : ((p - 1)! : ℤ) ≡ -1 [ZMOD p] := by
  have := Fact.mk hp
  rw [← ZMod.intCast_eq_intCast_iff]
  push_cast
  exact ZMod.wilsons_lemma p

/-- **Wilson's test.** For `2 ≤ n`, `(n - 1)! ≡ -1` mod `n` exactly when `n` is prime. This is
mathlib's `Nat.prime_iff_fac_equiv_neg_one`. -/
theorem wilson_iff (n : ℕ) (hn : 2 ≤ n) : ((n - 1)! : ℤ) ≡ -1 [ZMOD n] ↔ n.Prime := by
  rw [Nat.prime_iff_fac_equiv_neg_one (by omega), ← ZMod.intCast_eq_intCast_iff]
  push_cast
  rfl

/-! ### Small checks -/

example : 6 ! = 720 ∧ 720 % 7 = 6 := by decide

example : carries 3 5 4 0 = 2 ∧ Nat.choose 9 4 = 126 ∧ multiplicity 3 126 = 2 := by
  refine ⟨?_, by decide, multiplicity_eq_of_dvd_of_not_dvd (by decide) (by decide)⟩
  simp [carries]

end MultAdd.Cube
