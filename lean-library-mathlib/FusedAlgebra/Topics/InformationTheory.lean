module

public import FusedAlgebra.Basic

/-!
# Information theory as weighted mult-add

A small kit of mult-add lemmas proves a family of known information-theory laws. The weighted
mult-add `wAdd k a b = a + b + k * a * b` says how the corner cell `a * b` of the mult-add square
is counted, and the weight is read off the law:

* `k = -1`: failure chances of independent parts (`union_rule`, `series_system`), Gini–Simpson
  entropy (`gini_simpson`), Kraft slack (`kraft_slack`), erasure channels (`bec_capacity`);
* `k = -2`: binary symmetric channels (`bsc_cascade`, `bsc_polarization`), Matsui's piling-up
  lemma (`piling_up`), parity checks (`n_flips`, `parity3`, `check_node`);
* `k = -(1 + c)`: the `q`-ary symmetric channel (`qary_cascade`);
* `k = 1 - q`: Tsallis pseudo-additivity (`tsallis`) and Rényi additivity (`renyi_additive`);
* `k = 0`: Shannon additivity (`shannon_additive`), for any `L` with `L (x * y) = L x + L y`.

## Main definitions

* `InfoFused.outer`: all products `xᵢ * yⱼ` of two lists, the joint distribution of two
  independent sources.
* `InfoFused.shan`: Shannon's sum `∑ x * L x` with `L` in place of the logarithm.
* `InfoFused.vadd`: velocity addition on (numerator, denominator) pairs.

All of this is known mathematics; what the file records is that each law is a statement about
weighted mult-add. Nothing here needs real numbers: where a law needs `log`, it is stated for any
`L` that turns products into sums. The statements are purely algebraic and are not connected to
mathlib's measure-theoretic entropy. The core-Lean original worked over `ℚ`; here the identities
hold in any commutative ring and the inequalities in any linearly ordered commutative ring (or
field of characteristic zero, where a `1 / 2` or `1 / 3` appears).
-/

@[expose] public section

namespace InfoFused

open MultAdd

variable {R : Type*}

/-! ### The common kit -/

section Kit

variable [CommSemiring R]

/-- Any map that turns products into sums turns `wAdd k` into `+` in the coordinate `1 + k * a`.
Taking `L = log` gives `log (1 + k * a)`. -/
theorem linearize {G : Type*} [Add G] (L : R → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (k a b : R) : L (1 + k * wAdd k a b) = L (1 + k * a) + L (1 + k * b) := by
  rw [one_add_mul_wAdd, hL]

/-- All products `xᵢ * yⱼ`, row by row: the cells of the outer-product grid. -/
def outer (xs ys : List R) : List R := (xs ×ˢ ys).map fun p ↦ p.1 * p.2

@[simp]
theorem outer_nil (ys : List R) : outer [] ys = [] := rfl

theorem outer_cons (x : R) (xs ys : List R) :
    outer (x :: xs) ys = ys.map (x * ·) ++ outer xs ys := by
  simp [outer, List.product_cons, Function.comp_def]

/-- The grid's total is the product of the two edges. -/
theorem outer_sum (xs ys : List R) : (outer xs ys).sum = xs.sum * ys.sum := by
  induction xs with
  | nil => simp
  | cons x xs ih => rw [outer_cons, List.sum_append, List.sum_map_mul_left, ih, List.sum_cons,
      List.map_id', add_mul]

/-- Shannon's sum `∑ x * L x`, with `L` in place of `log` (it is `-H` when `L = log`). -/
def shan (L : R → R) (xs : List R) : R := (xs.map fun x ↦ x * L x).sum

@[simp]
theorem shan_nil (L : R → R) : shan L [] = 0 := rfl

@[simp]
theorem shan_cons (L : R → R) (x : R) (xs : List R) :
    shan L (x :: xs) = x * L x + shan L xs := by
  simp [shan]

theorem shan_append (L : R → R) (xs ys : List R) :
    shan L (xs ++ ys) = shan L xs + shan L ys := by
  simp [shan]

theorem shan_map (L : R → R) (hL : ∀ x y, L (x * y) = L x + L y) (x : R) (ys : List R) :
    shan L (ys.map (x * ·)) = x * ys.sum * L x + x * shan L ys := by
  induction ys with
  | nil => simp
  | cons y ys ih => simp only [List.map_cons, shan_cons, List.sum_cons, ih, hL]; ring

/-- The outer-product rule for Shannon's sum. -/
theorem outer_shan (L : R → R) (hL : ∀ x y, L (x * y) = L x + L y) (xs ys : List R) :
    shan L (outer xs ys) = ys.sum * shan L xs + xs.sum * shan L ys := by
  induction xs with
  | nil => simp
  | cons x xs ih => rw [outer_cons, shan_append, shan_map L hL, ih, shan_cons, List.sum_cons]; ring

end Kit

/-! ### Failure chances: the `k = -1` rule -/

section FailureChances

variable [CommRing R]

/-- Union rule for independent events: `P (A or B) = p + q - p * q = wAdd (-1) p q`. This is
`MultAdd.wAdd_neg_one`. -/
theorem union_rule (p q : R) : wAdd (-1) p q = 1 - (1 - p) * (1 - q) := wAdd_neg_one p q

/-- Survival multiplies. -/
theorem survival (p q : R) : 1 - wAdd (-1) p q = (1 - p) * (1 - q) := by
  rw [wAdd_neg_one]; ring

/-- The `k = -1` rule is mult-add read on the negatives. This is
`MultAdd.wAdd_neg_one_eq_neg_multAdd`. -/
theorem union_is_flipped_madd (p q : R) : wAdd (-1) p q = -(-p ✳ -q) :=
  wAdd_neg_one_eq_neg_multAdd p q

/-- Inclusion–exclusion for three independent events falls out of associativity: the signed
corners of the mult-add cube. -/
theorem incl_excl3 (p q r : R) :
    wAdd (-1) (wAdd (-1) p q) r = p + q + r - p * q - p * r - q * r + p * q * r := by
  unfold wAdd; ring

/-- A series system of `n` parts, each failing with chance `p`, fails with chance
`1 - (1 - p) ^ n`. This is `MultAdd.wpow_neg_one`. -/
theorem series_system (p : R) (n : ℕ) : wpow (-1) p n = 1 - (1 - p) ^ n := wpow_neg_one p n

/-- Survival multiplies, so the cumulative hazard `-log (1 - p)` adds. -/
theorem hazards_add {G : Type*} [Add G] (L : R → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (p q : R) : L (1 - wAdd (-1) p q) = L (1 - p) + L (1 - q) := by
  rw [survival, hL]

/-- Adding an independent part never lowers the chance of a failure. -/
theorem union_grows [LinearOrder R] [IsStrictOrderedRing R] (p q : R) (hp : p ≤ 1)
    (hq : 0 ≤ q) : p ≤ wAdd (-1) p q := by
  have h := mul_nonneg hq (sub_nonneg.2 hp)
  unfold wAdd; linarith

end FailureChances

/-! ### Entropy -/

section Entropy

/-- Hartley: typical-set sizes multiply, so bits add. This is `pow_add`. -/
theorem hartley [Semiring R] (m n : ℕ) : (2 : R) ^ (m + n) = 2 ^ m * 2 ^ n := pow_add 2 m n

section CommSemiring

variable [CommSemiring R]

/-- Shannon additivity for independent sources, `H (X, Y) = H X + H Y`, proved for any `L` with
`L (x * y) = L x + L y` and both distributions summing to `1`. -/
theorem shannon_additive (L : R → R) (hL : ∀ x y, L (x * y) = L x + L y) (ps rs : List R)
    (hp : ps.sum = 1) (hr : rs.sum = 1) :
    shan L (outer ps rs) = shan L ps + shan L rs := by
  rw [outer_shan L hL, hp, hr, one_mul, one_mul]

/-- The joint distribution of independent sources is a distribution. -/
theorem joint_sums_to_one (ps rs : List R) (hp : ps.sum = 1) (hr : rs.sum = 1) :
    (outer ps rs).sum = 1 := by
  rw [outer_sum, hp, hr, one_mul]

/-- Power sums multiply. With `xᵢ = pᵢ ^ q` and `yⱼ = rⱼ ^ q` this is
`∑ (pᵢ * rⱼ) ^ q = (∑ pᵢ ^ q) * (∑ rⱼ ^ q)`. -/
theorem power_sums_multiply (ws vs : List R) : (outer ws vs).sum = ws.sum * vs.sum :=
  outer_sum ws vs

/-- Rényi additivity: Rényi entropy is `L (1 + k * S) / k`, so it adds where Tsallis mult-adds.
Stated with the `1 / k` cleared. -/
theorem renyi_additive {G : Type*} [Add G] (L : R → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (k SA SB : R) : L (1 + k * wAdd k SA SB) = L (1 + k * SA) + L (1 + k * SB) :=
  linearize L hL k SA SB

/-- Landmark `q = 1` (`k = 0`): Shannon, plain addition. This is `MultAdd.wAdd_weight_zero`. -/
theorem shannon_weight (a b : R) : wAdd 0 a b = a + b := wAdd_weight_zero a b

end CommSemiring

variable [CommRing R]

/-- Tsallis pseudo-additivity, `S_q (A, B) = S_A + S_B + (1 - q) * S_A * S_B`. With `k = 1 - q`,
Tsallis entropy is `S = (P - 1) / k` where `P` is the power sum. Stated with the `k` cleared: if
`k * S_A = P - 1`, `k * S_B = R - 1` and `k * S_AB = P * R - 1`, then
`k * S_AB = k * wAdd k S_A S_B`. -/
theorem tsallis (k P Q SA SB SAB : R) (hA : k * SA = P - 1) (hB : k * SB = Q - 1)
    (hAB : k * SAB = P * Q - 1) : k * SAB = k * wAdd k SA SB := by
  have h := one_add_mul_wAdd k SA SB
  rw [hA, hB] at h
  linear_combination hAB - h

/-- Landmark `q = 2` (`k = -1`): Gini–Simpson `G = 1 - ∑ p ^ 2` is the chance two draws differ.
Collision sums multiply, so `G` composes by the union rule. -/
theorem gini_simpson (QA QB : R) : 1 - QA * QB = wAdd (-1) (1 - QA) (1 - QB) := by
  unfold wAdd; ring

/-- Landmark `q = 0` (`k = 1`): `S₀ = N - 1`, support sizes multiply, so `N - 1` mult-adds: the
`N × M` grid minus its unit corner. -/
theorem support_count (N M : R) : N * M - 1 = (N - 1) ✳ (M - 1) := by
  unfold multAdd; ring

/-- For `N` equally likely outcomes and `k = 1 - q` a natural number `m`, `1 + k * S = N ^ m`, so
`k * S` is the mult-add power of `N - 1`. -/
theorem uniform_tsallis_is_fused_power (N : R) (m : ℕ) : wpow 1 (N - 1) m = N ^ m - 1 := by
  rw [wpow_one, mpow_eq, add_sub_cancel]

end Entropy

/-! ### Compression -/

section Compression

variable [CommRing R]

/-- Kraft sums of two prefix codes multiply when codewords are concatenated, and the slack
`1 - K` composes like failure chances (`k = -1`). -/
theorem kraft_slack (KA KB : R) : 1 - KA * KB = wAdd (-1) (1 - KA) (1 - KB) :=
  gini_simpson KA KB

/-- The concatenated code still satisfies Kraft. -/
theorem kraft_product [LinearOrder R] [IsStrictOrderedRing R] (KA KB : R) (h0 : 0 ≤ KA)
    (hA : KA ≤ 1) (hB : KB ≤ 1) : KA * KB ≤ 1 := by
  have h := mul_nonneg h0 (sub_nonneg.2 hB)
  linarith

/-- Campbell's exponential cost `M = ∑ pᵢ * D ^ (t * ℓᵢ)` multiplies for concatenated independent
sources, so `M - 1` composes by plain mult-add. -/
theorem campbell (MA MB : R) : MA * MB - 1 = (MA - 1) ✳ (MB - 1) := support_count MA MB

/-- Expected length adds for concatenated codes (the `k = 0` shadow of Campbell). -/
theorem lengths_add (EA EB : R) : wAdd 0 EA EB = EA + EB := wAdd_weight_zero EA EB

end Compression

/-! ### Channels: the alphabet ladder `k = -(1 + c)` -/

section Channels

variable [CommRing R]

/-- Binary symmetric channel cascade: the off-diagonal entry of the product of
`!![1 - p, p; p, 1 - p]` and `!![1 - q, q; q, 1 - q]` is `wAdd (-2) p q`. -/
theorem bsc_cascade (p q : R) : (1 - p) * q + p * (1 - q) = wAdd (-2) p q := by
  unfold wAdd; ring

theorem bsc_cascade_diag (p q : R) : (1 - p) * (1 - q) + p * q = 1 - wAdd (-2) p q := by
  unfold wAdd; ring

/-- Polarization `1 - 2 * p` multiplies along a cascade. -/
theorem bsc_polarization (p q : R) : 1 - 2 * wAdd (-2) p q = (1 - 2 * p) * (1 - 2 * q) := by
  unfold wAdd; ring

/-- `n` identical flips: `1 - 2 * P (odd) = (1 - 2 * p) ^ n`. -/
theorem n_flips (p : R) (n : ℕ) : 1 - 2 * wpow (-2) p n = (1 - 2 * p) ^ n := by
  have h := one_add_mul_wpow (-2) p n
  have e : (1 : R) + -2 * p = 1 - 2 * p := by ring
  rw [e] at h
  rw [← h]; ring

/-- The `q`-ary symmetric channel. An error moves a symbol to one of the other letters at random;
a second error lands back on the right letter with chance `c = 1 / (alphabet - 1)`. The error
composes with weight `k = -(1 + c)`. -/
theorem qary_cascade (e₁ e₂ c : R) :
    1 - ((1 - e₁) * (1 - e₂) + c * e₁ * e₂) = wAdd (-(1 + c)) e₁ e₂ := by
  unfold wAdd; ring

/-- `c = 0`, letters never cancel (huge alphabet, or erasures): the union rule. -/
theorem ladder_never_cancel (e₁ e₂ : R) : wAdd (-(1 + 0)) e₁ e₂ = wAdd (-1) e₁ e₂ := by
  rw [add_zero]

/-- `c = 1`, two letters: a second flip always undoes the first. -/
theorem ladder_binary (e₁ e₂ : R) : wAdd (-(1 + 1)) e₁ e₂ = wAdd (-2) e₁ e₂ := by
  rw [one_add_one_eq_two]

/-- Binary erasure channel: erasures compose by the union rule and capacity `1 - e` multiplies
exactly. -/
theorem bec_capacity (e₁ e₂ : R) : 1 - wAdd (-1) e₁ e₂ = (1 - e₁) * (1 - e₂) := survival e₁ e₂

end Channels

section ChannelsField

variable {K : Type*} [Field K] [CharZero K]

/-- Matsui's piling-up lemma: biases `ε = 1 / 2 - p` combine as `2 * ε₁ * ε₂`. -/
theorem piling_up (p q : K) : 1 / 2 - wAdd (-2) p q = 2 * (1 / 2 - p) * (1 / 2 - q) := by
  unfold wAdd; ring

/-- `c = 1 / 3`, four letters: the depolarizing qubit channel's `-4 / 3`. -/
theorem ladder_four (e₁ e₂ : K) : wAdd (-(1 + 1 / 3)) e₁ e₂ = wAdd (-4 / 3) e₁ e₂ := by
  norm_num

end ChannelsField

/-! ### Data processing: a cascade is never better than its worst link -/

section DataProcessing

/-- BEC: the capacity of the cascade is at most the capacity of each link. This is
`mul_le_of_le_one_right`. -/
theorem bec_dpi [MulOneClass R] [Zero R] [Preorder R] [PosMulMono R] (C₁ C₂ : R) (h₁ : 0 ≤ C₁)
    (h₂ : C₂ ≤ 1) : C₁ * C₂ ≤ C₁ :=
  mul_le_of_le_one_right h₁ h₂

/-- Strong data processing: the contraction coefficient `(1 - 2 * p) ^ 2` multiplies. -/
theorem sdpi_multiplies [CommRing R] (p q : R) :
    (1 - 2 * wAdd (-2) p q) ^ 2 = (1 - 2 * p) ^ 2 * (1 - 2 * q) ^ 2 := by
  rw [bsc_polarization, mul_pow]

/-- Parallel channels: distinguishable messages multiply, capacities add (Hartley). -/
theorem parallel_capacity [Semiring R] (C₁ C₂ : ℕ) : (2 : R) ^ (C₁ + C₂) = 2 ^ C₁ * 2 ^ C₂ :=
  pow_add 2 C₁ C₂

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- BSC (crossovers taken `≤ 1 / 2`): the cascade's crossover is at least each link's. -/
theorem bsc_dpi (p q : K) (hp : p ≤ 1 / 2) (hq : 0 ≤ q) : p ≤ wAdd (-2) p q := by
  have h := mul_nonneg hq (show 0 ≤ 1 - 2 * p by linarith)
  unfold wAdd; linarith

/-- The cascade's crossover is still `≤ 1 / 2`, so it stays on the side where the capacity
`1 - h p` falls as `p` rises. -/
theorem bsc_stays_below_half (p q : K) (hp : p ≤ 1 / 2) (hq : q ≤ 1 / 2) :
    wAdd (-2) p q ≤ 1 / 2 := by
  have h := mul_nonneg (sub_nonneg.2 hp) (sub_nonneg.2 hq)
  have e := piling_up p q
  linarith

end DataProcessing

/-! ### Error-correcting codes -/

section Codes

variable [CommRing R]

/-- Parity of three noisy bits: the odd layers of the mult-add cube (`1, 3, 3, 1`). -/
theorem parity3 (p : R) : wpow (-2) p 3 = 3 * p * (1 - p) ^ 2 + p ^ 3 := by
  simp only [wpow, wAdd]; ring

/-- Every corner of the cube `{0, 1} ^ n`: `(1 + 1) ^ n - 1` nonempty corners. -/
theorem cube_corners (n : ℕ) : wpow 1 (1 : R) n = 2 ^ n - 1 := by
  rw [wpow_one, mpow_eq, one_add_one_eq_two]

/-- Hamming codes are perfect: with `n = m + r` and `n + 1 = 2 ^ r`, the `2 ^ m` balls of radius
`1` (a corner plus its `n` neighbours: layers `0` and `1`) tile the `n`-cube. -/
theorem hamming_perfect (r m : ℕ) (h : m + r + 1 = 2 ^ r) :
    2 ^ m * (1 + (m + r)) = 2 ^ (m + r) := by
  rw [add_comm 1, h, ← pow_add]

end Codes

/-! ### Decoding: belief propagation is mult-add in two coordinates

A bit's belief: chance `p` of a `1`, odds `o = (1 - p) / p`, `ρ = 1 - 2 * p = tanh (LLR / 2)`.
Check node (parity of two independent bits): `ρ` multiplies, `-2 * p` mult-adds. Variable node
(two independent looks at one bit): odds multiply, `o - 1` mult-adds, and in `ρ` that is
Einstein's velocity addition `(ρ₁ + ρ₂) / (1 + ρ₁ * ρ₂)`. -/

section Decoding

variable [CommRing R]

/-- Check node: `-2 * wAdd (-2) p q = (-2 * p) ✳ (-2 * q)`, plain mult-add. This is
`MultAdd.mul_wAdd` at `k = -2`. -/
theorem check_node (p q : R) : -2 * wAdd (-2) p q = (-2 * p) ✳ (-2 * q) := mul_wAdd (-2) p q

/-- Variable node: posterior odds multiply, so `o - 1` mult-adds. -/
theorem variable_node (o₁ o₂ : R) : o₁ * o₂ - 1 = (o₁ - 1) ✳ (o₂ - 1) := support_count o₁ o₂

/-- Odds product in `ρ` coordinates (`ρ = (o - 1) / (o + 1)`), with denominators cleared:
numerator and denominator of the velocity sum. -/
theorem odds_to_velocity (o₁ o₂ : R) :
    (o₁ - 1) * (o₂ + 1) + (o₂ - 1) * (o₁ + 1) = 2 * (o₁ * o₂ - 1) ∧
    (o₁ + 1) * (o₂ + 1) + (o₁ - 1) * (o₂ - 1) = 2 * (o₁ * o₂ + 1) := by
  constructor <;> ring

/-- Velocity addition lifts to mult-add on both corners, so the odds `(1 + v) / (1 - v)`
multiply. -/
theorem velocity_corners (a b : R) :
    (1 + a * b) + (a + b) = (1 + a) * (1 + b) ∧ (1 + a * b) - (a + b) = (1 - a) * (1 - b) := by
  constructor <;> ring

/-- Velocity addition on (numerator, denominator) pairs: `(a / b, c / d) ↦ (a * d + c * b) /
(b * d + a * c)`. -/
def vadd (x y : R × R) : R × R := (x.1 * y.2 + y.1 * x.2, x.2 * y.2 + x.1 * y.1)

/-- Velocity addition is associative. -/
theorem vadd_assoc (x y z : R × R) : vadd (vadd x y) z = vadd x (vadd y z) := by
  simp only [vadd, Prod.mk.injEq]; constructor <;> ring

end Decoding

end InfoFused
