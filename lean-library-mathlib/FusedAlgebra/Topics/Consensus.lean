module

public import FusedAlgebra.Basic

/-!
# Consensus protocols as mult-add

Alternatives to proof of work, read through weighted mult-add `wAdd k a b = a + b + k * a * b`:
`k = 1` is plain mult-add, `k = -1` is "either one happens".

## Main results

* Proof of work. An attacker with share `q`, `z` blocks behind, catches up with probability
  `r ^ z` where `r = q / (1 - q)`; with `a = r - 1` this is `1 + mpow a z`. `step_formula`,
  `half_is_identity`, `depth_useless`: at `q = 1 / 2` the step is the mult-add identity, and
  waiting more blocks buys no security.
* BFT voting. `bft_safe`, `bft_unsafe`: two quorums of size `n - f` overlap in more than `f`
  validators exactly when `3 * f < n`. `three_quorums`, `one_quorum`: the `1 / 4` and `1 / 2`
  thresholds.
* Committees. `maj3`: a majority of three goes bad on the top two layers of the cube;
  `maj3_drift`, `maj3_shrinks`: the pull toward `0` below `1 / 2`.
* Sybil resistance. `sybil_corner`, `sybil_gain`: in the multiplicative coordinate, splitting a
  stake gains the corner square `s₁ * s₂`; `additive_split`: an additive weight is split-proof.
* Slashing. `caught`: `m` independent watchers catch a double vote with chance
  `1 - (1 - d) ^ m`.
* Verifiable delay functions. `dbl`, `one_add_dbl`, `pietrzak_halves`: squaring is `a ✳ a`, and
  `n` squarings lift to the power `2 ^ n`.

The probability models themselves (that a race is a gambler's ruin, that a committee is a
binomial draw) are not part of this file. The numerical examples (`quarter_six`,
`coin_age_flip`, `fused_age_flip`, `linear_keeps`, `merge_cut`) stay over `ℚ`; everything else
is stated in the generality its proof allows.
-/

@[expose] public section

namespace ConsensusFused

open MultAdd

variable {R : Type*} {K : Type*}

/-! ### The algebra -/

/-- `k = -1` is "at least one of two independent events". This is `MultAdd.wAdd_neg_one`. -/
theorem either [CommRing R] (a b : R) : wAdd (-1) a b = 1 - (1 - a) * (1 - b) :=
  wAdd_neg_one a b

/-! ### Proof of work: depth is a mult-add power, `1 / 2` is the identity -/

/-- The step `a = r - 1` written without division. -/
theorem step_formula [Field K] (q : K) (h : q ≠ 1) :
    (q / (1 - q) - 1) * (1 - q) = 2 * q - 1 := by
  have h1 : 1 - q ≠ 0 := sub_ne_zero.2 h.symm
  field_simp; ring

/-- At the 51% line the step is exactly the mult-add identity. -/
theorem half_is_identity [Field K] [NeZero (2 : K)] : (1 / 2 : K) / (1 - 1 / 2) - 1 = 0 := by
  have h : (1 : K) - 1 / 2 = 1 / 2 := by field_simp; ring
  rw [h, div_self (one_div_ne_zero two_ne_zero), sub_self]

/-- Powers of the identity do nothing: at `q = 1 / 2` waiting more blocks buys no security. -/
@[simp]
theorem depth_useless [CommSemiring R] (k : R) (n : ℕ) : wpow k 0 n = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [wpow_succ, ih, wAdd_zero]

/-- Catch-up probability as a mult-add power, checked exactly for `q = 1 / 4`, `z = 6`:
`r = 1 / 3`, `a = -2 / 3`, and `1 + mpow a 6 = (1 / 3) ^ 6`. -/
theorem quarter_six : 1 + wpow 1 (-2 / 3 : ℚ) 6 = 1 / 729 := by
  rw [wpow_one, one_add_mpow]; norm_num

/-! ### BFT voting: where `1 / 3` comes from

Two quorums of size `Q` among `n` validators overlap in at least `2 * Q - n`. Safety: the overlap
must hold an honest validator, `2 * Q - n > f`. Liveness: a quorum must form without the `f`
faulty ones, `Q = n - f`. -/

theorem bft_safe (n f : ℕ) (h : 3 * f < n) : 2 * (n - f) > n + f := by lia

theorem bft_unsafe (n f : ℕ) (h : n ≤ 3 * f) : 2 * (n - f) ≤ n + f := by lia

/-- Three quorums sharing an honest validator need `f < n / 4`. -/
theorem three_quorums (n f : ℕ) (h : 4 * f < n) : 3 * (n - f) > 2 * n + f := by lia

/-- The longest-chain rule needs only a majority: one "quorum" against the rest. -/
theorem one_quorum (n f : ℕ) (h : 2 * f < n) : n - f > f := by lia

/-! ### Committees and delegation: the layers of the cube

Three independent members, each faulty with chance `p`: the cube's layers hold `(1 - p) ^ 3`,
`3 * p * (1 - p) ^ 2`, `3 * p ^ 2 * (1 - p)`, `p ^ 3` (Pascal's `1, 3, 3, 1`). A majority of
three goes bad on the top two layers. -/

section Committees

variable [CommRing R]

/-- The chance that a majority of three independent members is faulty. -/
def maj3 (p : R) : R := p ^ 3 + 3 * p ^ 2 * (1 - p)

theorem layers_sum (p : R) :
    (1 - p) ^ 3 + 3 * p * (1 - p) ^ 2 + 3 * p ^ 2 * (1 - p) + p ^ 3 = 1 := by ring

theorem maj3_closed (p : R) : maj3 p = 3 * p ^ 2 - 2 * p ^ 3 := by unfold maj3; ring

theorem maj3_mirror (p : R) : maj3 (1 - p) = 1 - maj3 p := by unfold maj3; ring

/-- The pull toward `0` below `1 / 2` and toward `1` above it. -/
theorem maj3_drift (p : R) : maj3 p - p = -(p * (1 - p) * (1 - 2 * p)) := by unfold maj3; ring

end Committees

theorem maj3_half [Field K] [NeZero (2 : K)] : maj3 (1 / 2 : K) = 1 / 2 := by
  unfold maj3; field_simp; ring

theorem maj3_shrinks [Field K] [LinearOrder K] [IsStrictOrderedRing K] (p : K) (h0 : 0 < p)
    (h1 : p < 1 / 2) : maj3 p < p := by
  have a : 0 < 1 - p := by linarith
  have b : 0 < 1 - 2 * p := by linarith
  rw [← sub_neg, maj3_drift, neg_lt_zero]
  positivity

/-! ### The Sybil corner

Weighting stake by `log (1 + s)` makes weight additive under mult-add. Real stake adds with `+`.
In the multiplicative coordinate `W = 1 + s`, splitting `s₁ + s₂` into two accounts turns
`1 + s₁ + s₂` into `(1 + s₁) * (1 + s₂)`: the gain is the corner square `s₁ * s₂`, never
negative. -/

theorem sybil_corner [CommRing R] (s₁ s₂ : R) :
    (1 + s₁) * (1 + s₂) - (1 + (s₁ + s₂)) = s₁ * s₂ := by ring

theorem sybil_gain [CommRing R] [PartialOrder R] [IsOrderedRing R] (s₁ s₂ : R) (h₁ : 0 ≤ s₁)
    (h₂ : 0 ≤ s₂) : 1 + (s₁ + s₂) ≤ (1 + s₁) * (1 + s₂) := by
  rw [← sub_nonneg, sybil_corner]
  exact mul_nonneg h₁ h₂

/-- An additive weight is split-proof: `n` equal accounts weigh what one does. This is
`map_nsmul` for the additive map `f`. -/
theorem additive_split {S : Type*} [NonAssocSemiring R] [NonAssocRing S] (f : R → S)
    (hf : ∀ a b, f (a + b) = f a + f b) (s : R) (n : ℕ) : f (n * s) = n * f s := by
  simpa using map_nsmul (AddMonoidHom.mk' f hf) n s

/-! ### Fusing stake with age fails

Two parties. The attacker holds `1 / 5` of the stake and has been silent for `10` epochs; the
honest side holds `4 / 5` and last produced `1` epoch ago. -/

/-- Coin age `s * (1 + t)`: the 20% attacker outweighs the 80% honest side. -/
theorem coin_age_flip : (4 / 5 : ℚ) * (1 + 1) < 1 / 5 * (1 + 10) := by norm_num

/-- Fused exponent `(1 + s) ^ (1 + t) - 1`: the same flip, more strongly. -/
theorem fused_age_flip : (1 + 4 / 5 : ℚ) ^ (1 + 1) - 1 < (1 + 1 / 5) ^ (1 + 10) - 1 := by
  norm_num

/-- Plain stake keeps the honest majority. -/
theorem linear_keeps : (1 / 5 : ℚ) < 4 / 5 := by norm_num

/-! ### Slashing and long-range attacks: `k = -1` powers -/

/-- `m` independent watchers each catching a double vote with chance `d`: caught with chance
`1 - (1 - d) ^ m`. The same holds for the share of an old validator set that has left after `t`
epochs at churn `c`. This is `MultAdd.wpow_neg_one`. -/
theorem caught [CommRing R] (d : R) (m : ℕ) : wpow (-1) d m = 1 - (1 - d) ^ m :=
  wpow_neg_one d m

/-! ### Verifiable delay functions: doubling up the cube

A VDF computes `x ^ 2 ^ T` by `T` squarings. In mult-add coordinates squaring is `a ✳ a`.
Pietrzak's proof checks the two halves: `(x ^ 2 ^ m) ^ 2 ^ n = x ^ 2 ^ (m + n)`. -/

/-- `n` squarings in mult-add coordinates: `a ↦ a ✳ a`, iterated `n` times. -/
def dbl [Add R] [Mul R] (n : ℕ) (a : R) : R := (fun x ↦ x ✳ x)^[n] a

@[simp]
theorem dbl_zero [Add R] [Mul R] (a : R) : dbl 0 a = a := rfl

theorem dbl_succ [Add R] [Mul R] (n : ℕ) (a : R) : dbl (n + 1) a = dbl n a ✳ dbl n a :=
  Function.iterate_succ_apply' _ n a

/-- `n` squarings lift to the power `2 ^ n`. -/
theorem one_add_dbl [Semiring R] (a : R) (n : ℕ) : 1 + dbl n a = (1 + a) ^ 2 ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [dbl_succ, one_add_multAdd, ih, pow_succ 2, pow_mul, sq]

/-- `n` squarings are the mult-add power `2 ^ n`. -/
theorem dbl_eq_mpow [Ring R] (a : R) (n : ℕ) : dbl n a = mpow a (2 ^ n) := by
  rw [mpow_eq, ← one_add_dbl, add_sub_cancel_left]

/-- Pietrzak's two halves. This is `pow_mul` and `pow_add`. -/
theorem pietrzak_halves [Monoid R] (x : R) (m n : ℕ) :
    (x ^ 2 ^ m) ^ 2 ^ n = x ^ 2 ^ (m + n) := by
  rw [← pow_mul, ← pow_add]

/-! ### Energy, from published figures (TWh per year) -/

/-- Ethereum before the Merge about 23 TWh, after about 0.0026 TWh (CCRI 2022): over 99.98%
less. -/
theorem merge_cut : (23 - 26 / 10000 : ℚ) / 23 > 9998 / 10000 := by norm_num

end ConsensusFused
