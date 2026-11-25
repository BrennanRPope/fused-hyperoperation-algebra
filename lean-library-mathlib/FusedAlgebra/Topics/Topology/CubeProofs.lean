module

public import FusedAlgebra.Cube

/-!
# Classical proofs on the mult-add cube: Lucas, the AKS congruence, Fermat

The `n`-cube behind `1 + mpow x n = (1 + x) ^ n` has `n.choose k` corners with `k` switches on.
This file reads four classical facts off those layer counts.

* `MultAdd.Cube.binomial`: `(1 + x) ^ n = ∑ k, n.choose k * x ^ k`, the cube's corner pieces add
  up to the box.
* `MultAdd.Cube.prime_dvd_middle`: a prime `p` divides every middle layer of the `p`-cube.
* `MultAdd.Cube.aks_middle`: for `2 ≤ n`, `n` divides every middle layer of the `n`-cube exactly
  when `n` is prime (the identity behind the AKS primality test).
* `MultAdd.Cube.lucas_step`, `MultAdd.Cube.lucas`: Lucas's theorem (1878).
* `MultAdd.Cube.frobenius`, `MultAdd.Cube.mpow_prime_modEq`: `(1 + x) ^ p ≡ 1 + x ^ p` mod `p`,
  that is, the mult-add `p`-th power is the plain `p`-th power mod `p`.
* `MultAdd.Cube.fermat`: Fermat's little theorem, `a ^ p ≡ a` mod `p`.

This is the mathlib form of `FusedAlgebra.Topics.Topology.CubeProofs`. All of it is known
mathematics and almost all of it is already in mathlib, so the statements here are short
corollaries, kept so the library's names still exist:

* the binomial theorem is `add_pow` (through `MultAdd.Cube.one_add_mpow_eq_sum_choose`);
* the prime middle is `Nat.Prime.dvd_choose_self`;
* Lucas's theorem is `Choose.choose_modEq_choose_mod_mul_choose_div` and `Choose.lucas_theorem`;
* the freshman's dream is `add_pow_char`, and Fermat's little theorem is `ZMod.pow_card`.

Only the converse half of `aks_middle` is proved here. The core-Lean file built primes
(`IsPrime`), congruences (`MEq`) and finite sums (`bsum`) by hand; they are now `Nat.Prime`,
`Int.ModEq` (`a ≡ b [ZMOD p]`) and `Finset.sum` over `Finset.range`.
-/

@[expose] public section

namespace MultAdd.Cube

open Finset

/-! ### The binomial theorem as the cube's corner count -/

/-- **The cube's corner count.** `(1 + x) ^ n = ∑ k, n.choose k * x ^ k`: a box of side `1 + x`
cut at `1` on every axis falls into the `n`-cube's corner pieces. This is mathlib's `add_pow`. -/
theorem binomial {R : Type*} [CommSemiring R] (n : ℕ) (x : R) :
    (1 + x) ^ n = ∑ k ∈ range (n + 1), (n.choose k : R) * x ^ k := by
  rw [← one_add_mpow, one_add_mpow_eq_sum_choose]

/-! ### The middle of the prime cube -/

/-- **The middle of the prime cube.** For a prime `p`, every middle layer of the `p`-cube has a
multiple of `p` corners. This is mathlib's `Nat.Prime.dvd_choose_self`. -/
theorem prime_dvd_middle {p k : ℕ} (hp : p.Prime) (hk0 : 0 < k) (hkp : k < p) :
    p ∣ p.choose k :=
  hp.dvd_choose_self hk0.ne' hkp

/-! ### Lucas's theorem -/

/-- **Lucas, one step.** `n.choose k ≡ (n / p).choose (k / p) * (n % p).choose (k % p)` mod `p`:
split the axes into blocks of `p` plus a remainder; mod `p` only the corners whose blocks are
each all-off or all-on count. This is mathlib's `Choose.choose_modEq_choose_mod_mul_choose_div`. -/
theorem lucas_step {p : ℕ} (hp : p.Prime) (n k : ℕ) :
    (n.choose k : ℤ) ≡ (n / p).choose (k / p) * (n % p).choose (k % p) [ZMOD p] := by
  have := Fact.mk hp
  rw [mul_comm]
  exact Choose.choose_modEq_choose_mod_mul_choose_div

/-- **Lucas's theorem.** For `n`, `k` with at most `f` base-`p` digits, `n.choose k` is congruent
mod `p` to the product over the digits of `nᵢ.choose kᵢ`. This is mathlib's
`Choose.lucas_theorem`. -/
theorem lucas {p : ℕ} (hp : p.Prime) (f n k : ℕ) (hn : n < p ^ f) (hk : k < p ^ f) :
    (n.choose k : ℤ) ≡ ∏ i ∈ range f, ((n / p ^ i % p).choose (k / p ^ i % p) : ℤ) [ZMOD p] := by
  have := Fact.mk hp
  simpa using Choose.lucas_theorem hn hk

/-! ### Primes are exactly where the middle vanishes (the AKS identity) -/

/-- **The AKS identity on the cube.** For `2 ≤ n`: every middle layer of the `n`-cube has a
multiple of `n` corners exactly when `n` is prime. Equivalently `(1 + x) ^ n ≡ 1 + x ^ n` mod `n`:
the mult-add `n`-th power is the plain `n`-th power, the identity behind the AKS test. -/
theorem aks_middle (n : ℕ) (hn : 2 ≤ n) :
    (∀ k, 0 < k → k < n → n ∣ n.choose k) ↔ n.Prime := by
  refine ⟨fun hall ↦ ?_, fun hp k hk0 hkn ↦ prime_dvd_middle hp hk0 hkn⟩
  by_contra hnp
  obtain ⟨p, hp, hpn⟩ := Nat.exists_prime_and_dvd (show n ≠ 1 by lia)
  have hpne : p ≠ n := fun h ↦ hnp (h ▸ hp)
  have hp2 := hp.two_le
  have hplt : p < n := lt_of_le_of_ne (Nat.le_of_dvd (by lia) hpn) hpne
  obtain ⟨t, ht⟩ := hall p hp.pos hplt
  obtain ⟨q, rfl⟩ : ∃ q, n = q + 1 := ⟨n - 1, by lia⟩
  obtain ⟨j, rfl⟩ : ∃ j, p = j + 1 := ⟨p - 1, by lia⟩
  -- absorption: `(q + 1) * q.choose j = (q + 1).choose (j + 1) * (j + 1)`
  have hab := Nat.add_one_mul_choose_eq q j
  rw [ht, mul_assoc] at hab
  have hcancel : q.choose j = t * (j + 1) := Nat.eq_of_mul_eq_mul_left q.succ_pos hab
  -- Lucas: `q.choose j ≡ 1` mod `p`, because `q % p = j` and `j / p = 0`
  have hl := lucas_step hp q j
  have hqmod : q % (j + 1) = j := by
    obtain ⟨m, hm⟩ := hpn
    have hq : q = j + (j + 1) * (m - 1) := by
      cases m with
      | zero => lia
      | succ m => rw [Nat.add_sub_cancel]; lia
    rw [hq, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt j.lt_succ_self]
  rw [hqmod, Nat.div_eq_of_lt j.lt_succ_self, Nat.mod_eq_of_lt j.lt_succ_self,
    Nat.choose_zero_right, Nat.choose_self, hcancel] at hl
  have h0 : ((t * (j + 1) : ℕ) : ℤ) ≡ 0 [ZMOD (j + 1 : ℕ)] :=
    (Int.modEq_zero_iff_dvd).2 ⟨t, by push_cast; ring⟩
  have h1 : (1 : ℤ) ≡ 0 [ZMOD (j + 1 : ℕ)] := by simpa using hl.symm.trans h0
  have h2 : ((j + 1 : ℕ) : ℤ) ∣ 1 := Int.modEq_zero_iff_dvd.1 h1
  have := Int.le_of_dvd one_pos h2
  lia

/-! ### Fermat's little theorem -/

/-- **Mult-add powers are plain powers mod `p`.** `(1 + x) ^ p ≡ 1 + x ^ p`: in the `p`-cube only
the two end corners survive mod `p`. This is mathlib's `add_pow_char` in `ZMod p`. -/
theorem frobenius {p : ℕ} (hp : p.Prime) (x : ℤ) : (1 + x) ^ p ≡ 1 + x ^ p [ZMOD p] := by
  have := Fact.mk hp
  rw [← ZMod.intCast_eq_intCast_iff]
  push_cast
  rw [add_pow_char, one_pow]

/-- The same in mult-add terms: mod a prime `p`, the mult-add `p`-th power is the plain `p`-th
power. -/
theorem mpow_prime_modEq {p : ℕ} (hp : p.Prime) (x : ℤ) : mpow x p ≡ x ^ p [ZMOD p] := by
  have h := frobenius hp x
  rw [← one_add_mpow] at h
  exact Int.ModEq.add_left_cancel' 1 h

/-- **Fermat's little theorem.** `a ^ p ≡ a` mod `p`. This is mathlib's `ZMod.pow_card`. -/
theorem fermat {p : ℕ} (hp : p.Prime) (a : ℤ) : a ^ p ≡ a [ZMOD p] := by
  have := Fact.mk hp
  rw [← ZMod.intCast_eq_intCast_iff]
  push_cast
  exact ZMod.pow_card _

/-! ### Small checks -/

/-- `9` is not prime, and the `9`-cube shows it: layer `3` has `84` corners, not a multiple
of `9`. -/
example : Nat.choose 9 3 = 84 ∧ ¬ 9 ∣ 84 := by decide

/-- Lucas in base `7`: `10.choose 4 = 210`, and `10 = 13₇`, `4 = 04₇` gives
`3.choose 4 * 1.choose 0 = 0 ≡ 210`. -/
example : Nat.choose 10 4 = 210 ∧
    ∏ i ∈ range 2, (10 / 7 ^ i % 7).choose (4 / 7 ^ i % 7) = 0 ∧ 210 % 7 = 0 := by decide

end MultAdd.Cube
