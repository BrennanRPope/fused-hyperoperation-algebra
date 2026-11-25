module

public import FusedAlgebra.Cube
public import FusedAlgebra.Topics.E

/-!
# Statistics as weighted mult-add

Weighted mult-add is `wAdd k a b = a + b + k * a * b`. This file reads elementary probability
through it. It is the mathlib form of `FusedAlgebra.Topics.Statistics`.

## Contents

1. **Indicators.** On idempotents (`IsIdempotentElem`, the `0`/`1` numbers) the corner `a * b` is
   AND, weight `-1` is OR (`or_ind`), weight `-2` is XOR (`xor_ind`, `sq_diff_is_xor`), and the
   complement `x ↦ 1 - x` turns OR into AND (`de_morgan`). Inclusion-exclusion for `n` events is
   the weight `-1` fold `unionL`, the signed corners of the `n`-cube (`inclusion_exclusion`,
   `union_corners`).
2. **Finite probability spaces** (`expect`). For events that need not be independent,
   `P(A ∪ B) = wAdd (-lift) P(A) P(B)` with `lift = P(A ∩ B) / (P(A) * P(B))` (`union_lift`):
   disjoint is weight `0`, independent is weight `-1`. Inclusion-exclusion for `n` events with
   any dependence is `inclusion_exclusion_general`.
3. **The generating-function dial.** `gen d t = E[t ^ X]` multiplies under independent sums
   (`gen_nconv`), so with `t = 1 + k` the `k`-moment `a`, `1 + k * a = gen d (1 + k)`, composes
   by `wAdd k` (`kmoment_law`). A coin's `k`-moment is `p` at every `k` (`bern_kmoment`). `t = 0`
   gives `P(none)` and `t = -1` gives `P(even) - P(odd)` (`heads_none`, `heads_parity`). The
   binomial distribution is the weighted Pascal layers (`bp_layers`).
4. **Moments and cumulants.** Moments of `X + Y` follow the Pascal rows (`moment_conv`, with the
   rows `1 2 1`, `1 3 3 1`, `1 4 6 4 1` spelled out), and the cumulants `k2`, `k3`, `k4` add.
   For `n` copies the shape numbers fall exactly as `1 / n` (`clt_skew`, `clt_kurt`).
5. **Bayes.** Posterior odds are prior odds times the likelihood ratio (`bayes_odds`), so
   `odds - 1` mult-adds; the complement is the mult-add inverse; pooling two opinions is velocity
   addition in `ρ = 2 * p - 1` (`pool_velocity`).
6. **The dependence line.** Union under Ali-Mikhail-Haq dependence `θ` is the lens law at weight
   `-(1 + θ)` and bend `θ` (`amh_union_is_lens`), and a lens law lets certainty absorb exactly
   when `k + q = -1` (`certainty_absorbs_iff`).

All of this is known mathematics (Boolean algebra of idempotents, inclusion-exclusion,
probability generating functions, cumulants, Bayes' rule in odds form, the AMH copula). What the
file records is that each one is a statement about weighted mult-add. Mathlib has the measure
theoretic versions (`MeasureTheory.measure_union_add_inter`, `ProbabilityTheory.IndepFun`,
`PMF.binomial`, `ProbabilityTheory.variance`, …); here a distribution is a finite list of
weighted outcomes, as in the core-Lean original, and no positivity or normalisation of the
weights is assumed unless stated. The original worked over `ℚ`; here parts 1 to 4 hold in any
commutative ring and the parts that divide hold in any field.

Not here (needs limits): the normal density, the limit in the central limit theorem, the Poisson
law.
-/

@[expose] public section

namespace StatFused

open MultAdd Finset

variable {R : Type*} [CommRing R]

/-! ### 1. Indicators: Boolean logic is weighted mult-add on idempotents

An indicator is an idempotent, `x * x = x`: mathlib's `IsIdempotentElem`. In a field these are
exactly `0` and `1` (`IsIdempotentElem.iff_eq_zero_or_one`). AND is `IsIdempotentElem.mul` and
NOT is `IsIdempotentElem.one_sub`. -/

/-- OR of two indicators is an indicator: weight `-1`. -/
theorem or_ind {x y : R} (hx : IsIdempotentElem x) (hy : IsIdempotentElem y) :
    IsIdempotentElem (wAdd (-1) x y) := by
  have e : wAdd (-1) x y = x + y - x * y := by unfold wAdd; ring
  rw [e]; exact hx.add_sub_mul hy

/-- XOR of two indicators is an indicator: weight `-2`. -/
theorem xor_ind {x y : R} (hx : IsIdempotentElem x) (hy : IsIdempotentElem y) :
    IsIdempotentElem (wAdd (-2) x y) := by
  have hx' : x * x = x := hx
  have hy' : y * y = y := hy
  unfold IsIdempotentElem wAdd
  linear_combination hx' + (1 + 4 * x ^ 2 - 4 * x) * hy'

/-- The truth table of OR. -/
theorem or_table : wAdd (-1) (0 : R) 0 = 0 ∧ wAdd (-1) (0 : R) 1 = 1 ∧
    wAdd (-1) (1 : R) 0 = 1 ∧ wAdd (-1) (1 : R) 1 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (unfold wAdd; ring)

/-- The truth table of XOR. -/
theorem xor_table : wAdd (-2) (0 : R) 0 = 0 ∧ wAdd (-2) (0 : R) 1 = 1 ∧
    wAdd (-2) (1 : R) 0 = 1 ∧ wAdd (-2) (1 : R) 1 = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (unfold wAdd; ring)

/-- Weight `0` counts. -/
theorem count_table : wAdd 0 (1 : R) 1 = 2 := by
  rw [wAdd_weight_zero, one_add_one_eq_two]

/-- The squared difference of two indicators is their XOR (weight `-2`). -/
theorem sq_diff_is_xor {x y : R} (hx : IsIdempotentElem x) (hy : IsIdempotentElem y) :
    (x - y) * (x - y) = wAdd (-2) x y := by
  have hx' : x * x = x := hx
  have hy' : y * y = y := hy
  unfold wAdd
  linear_combination hx' + hy'

/-- De Morgan: "not (A or B)" is "not A and not B". The complement `x ↦ 1 - x` turns weight `-1`
into multiplication. -/
theorem de_morgan (x y : R) : 1 - wAdd (-1) x y = (1 - x) * (1 - y) := by
  unfold wAdd; ring

/-- The complement conjugates OR into AND. -/
theorem de_morgan_and (x y : R) : 1 - wAdd (-1) (1 - x) (1 - y) = x * y := by
  unfold wAdd; ring

/-- XOR is self-dual under complement of both inputs. -/
theorem xor_complement (x y : R) : wAdd (-2) (1 - x) (1 - y) = wAdd (-2) x y := by
  unfold wAdd; ring

/-! #### `n` events: the weight `-1` fold and its signed corners -/

/-- Union of a list: `x₁`, `x₂`, … folded with weight `-1`. -/
def unionL (l : List R) : R := l.foldr (wAdd (-1)) 0

@[simp]
theorem unionL_nil : unionL ([] : List R) = 0 := rfl

@[simp]
theorem unionL_cons (a : R) (l : List R) : unionL (a :: l) = wAdd (-1) a (unionL l) := rfl

/-- Nothing happens: `1 - union` is the product of the misses. -/
theorem one_sub_unionL (l : List R) : 1 - unionL l = (l.map (1 - ·)).prod := by
  induction l with
  | nil => simp
  | cons a l ih => rw [unionL_cons, de_morgan, ih, List.map_cons, List.prod_cons]

/-- The weight `-1` fold is mult-add read on the negatives. -/
theorem unionL_eq_neg_multAddList (l : List R) : unionL l = -multAddList (l.map Neg.neg) := by
  induction l with
  | nil => simp
  | cons a l ih =>
    rw [unionL_cons, List.map_cons, multAddList_cons, ih, wAdd_neg_one_eq_neg_multAdd, neg_neg]

/-- Inclusion-exclusion for `n` events: `1 - P(union)` is the sum over every corner `g` of the
`n`-cube of the product of the `-xᵢ` for `i` in `g`. -/
theorem inclusion_exclusion (l : List R) :
    1 - unionL l = ((l.map Neg.neg).sublists'.map List.prod).sum := by
  rw [unionL_eq_neg_multAddList, sub_neg_eq_add, one_add_multAddList_eq_sum_sublists]

/-- Inclusion-exclusion with the signs pulled out: corner `g` carries the sign `(-1) ^ |g|`. -/
theorem union_corners (l : List R) :
    1 - unionL l = (l.sublists'.map fun g ↦ (-1) ^ g.length * g.prod).sum := by
  rw [inclusion_exclusion, List.sublists'_map, List.map_map]
  congr 1
  exact List.map_congr_left fun g _ ↦ List.prod_map_neg g

theorem unionL_pair (p q : R) : unionL [p, q] = p + q - p * q := by
  simp [wAdd]; ring

/-! ### 2. Finite probability spaces: dependence is the weight -/

section Spaces

variable {α : Type*}

/-- A finite probability space is a list of outcomes with weights; `expect ws f` is the
expectation of `f`. -/
def expect (ws : List (α × R)) (f : α → R) : R := (ws.map fun p ↦ p.2 * f p.1).sum

@[simp]
theorem expect_nil (f : α → R) : expect ([] : List (α × R)) f = 0 := rfl

@[simp]
theorem expect_cons (p : α × R) (ws : List (α × R)) (f : α → R) :
    expect (p :: ws) f = p.2 * f p.1 + expect ws f := rfl

theorem expect_add (ws : List (α × R)) (f g : α → R) :
    expect ws (fun a ↦ f a + g a) = expect ws f + expect ws g := by
  induction ws with
  | nil => simp
  | cons p ws ih => simp only [expect_cons, ih]; ring

theorem expect_sub (ws : List (α × R)) (f g : α → R) :
    expect ws (fun a ↦ f a - g a) = expect ws f - expect ws g := by
  induction ws with
  | nil => simp
  | cons p ws ih => simp only [expect_cons, ih]; ring

theorem expect_const_mul (ws : List (α × R)) (c : R) (f : α → R) :
    expect ws (fun a ↦ c * f a) = c * expect ws f := by
  induction ws with
  | nil => simp
  | cons p ws ih => simp only [expect_cons, ih]; ring

theorem expect_const (ws : List (α × R)) (c : R) :
    expect ws (fun _ ↦ c) = c * expect ws (fun _ ↦ 1) := by
  simpa using expect_const_mul ws c (fun _ ↦ 1)

/-- Linearity over a finite list of functions. -/
theorem expect_list_sum {β : Type*} (ws : List (α × R)) (gs : List β) (F : β → α → R) :
    expect ws (fun a ↦ (gs.map fun g ↦ F g a).sum) =
      (gs.map fun g ↦ expect ws (F g)).sum := by
  induction gs with
  | nil => simpa using expect_const ws 0
  | cons g gs ih => simp only [List.map_cons, List.sum_cons]; rw [expect_add, ih]

/-! #### Events that need not be independent -/

/-- Union of two events, dependent or not: `P(A ∪ B) = P(A) + P(B) - P(A ∩ B)`. -/
theorem union_general (ws : List (α × R)) (A B : α → R) :
    expect ws (fun a ↦ wAdd (-1) (A a) (B a)) =
      expect ws A + expect ws B - expect ws (fun a ↦ A a * B a) := by
  rw [← expect_add, ← expect_sub]
  congr 1; funext a; unfold wAdd; ring

/-- Dependence is the weight. If the overlap is `lift` times what independence predicts,
`P(A ∩ B) = lift * P(A) * P(B)`, then `P(A ∪ B) = wAdd (-lift) P(A) P(B)`. -/
theorem union_lift (ws : List (α × R)) (A B : α → R) (lift : R)
    (h : expect ws (fun a ↦ A a * B a) = lift * expect ws A * expect ws B) :
    expect ws (fun a ↦ wAdd (-1) (A a) (B a)) = wAdd (-lift) (expect ws A) (expect ws B) := by
  rw [union_general, h]; unfold wAdd; ring

/-- Disjoint events (lift `0`): weight `0`, plain addition, Kolmogorov's axiom. -/
theorem union_disjoint (ws : List (α × R)) (A B : α → R)
    (h : expect ws (fun a ↦ A a * B a) = 0) :
    expect ws (fun a ↦ wAdd (-1) (A a) (B a)) = wAdd 0 (expect ws A) (expect ws B) := by
  rw [union_general, h, wAdd_weight_zero, sub_zero]

/-- Exclusive or with dependence: weight `-2 * lift`. -/
theorem xor_lift (ws : List (α × R)) (A B : α → R) (lift : R)
    (h : expect ws (fun a ↦ A a * B a) = lift * expect ws A * expect ws B) :
    expect ws (fun a ↦ wAdd (-2) (A a) (B a)) =
      wAdd (-2 * lift) (expect ws A) (expect ws B) := by
  have e : (fun a ↦ wAdd (-2) (A a) (B a)) = fun a ↦ A a + B a - 2 * (A a * B a) := by
    funext a; unfold wAdd; ring
  rw [e, expect_sub, expect_add, expect_const_mul, h]; unfold wAdd; ring

/-- Covariance of two indicators is the corner's excess over the independent corner. -/
theorem cov_is_corner_excess (ws : List (α × R)) (A B : α → R) (lift : R)
    (h : expect ws (fun a ↦ A a * B a) = lift * expect ws A * expect ws B) :
    expect ws (fun a ↦ A a * B a) - expect ws A * expect ws B =
      (lift - 1) * expect ws A * expect ws B := by
  rw [h]; ring

/-! #### Inclusion-exclusion for `n` events, with any dependence -/

/-- The indicator of "every event in `g` happens". -/
def allOf (g : List (α → R)) (a : α) : R := (g.map fun A ↦ A a).prod

/-- Inclusion-exclusion, `n` events, any dependence: `P(no event happens)` is the sum over every
corner `g` of the cube of `(-1) ^ |g| * P(all of g)`. -/
theorem inclusion_exclusion_general (ws : List (α × R)) (As : List (α → R)) :
    expect ws (fun a ↦ 1 - unionL (As.map fun A ↦ A a)) =
      (As.sublists'.map fun g ↦ (-1) ^ g.length * expect ws (allOf g)).sum := by
  have pt : (fun a ↦ 1 - unionL (As.map fun A ↦ A a)) =
      fun a ↦ (As.sublists'.map fun g ↦ (-1) ^ g.length * allOf g a).sum := by
    funext a
    rw [union_corners, List.sublists'_map, List.map_map]
    simp [Function.comp_def, allOf]
  rw [pt, expect_list_sum]
  simp only [expect_const_mul]

end Spaces

/-! ### 3. The generating-function dial `t = 1 + k` -/

/-- A distribution on counts: a list of (value, weight). -/
abbrev NDist (R : Type*) := List (ℕ × R)

/-- The generating function `gen d t = E[t ^ X]`. -/
def gen (d : NDist R) (t : R) : R := (d.map fun p ↦ p.2 * t ^ p.1).sum

@[simp]
theorem gen_nil (t : R) : gen ([] : NDist R) t = 0 := rfl

@[simp]
theorem gen_cons (p : ℕ × R) (d : NDist R) (t : R) :
    gen (p :: d) t = p.2 * t ^ p.1 + gen d t := rfl

/-- Shift every value by `x` and scale every weight by `w`. -/
def nshift (x : ℕ) (w : R) (e : NDist R) : NDist R := e.map fun p ↦ (x + p.1, w * p.2)

/-- The distribution of `X + Y` for independent counts `X`, `Y`. -/
def nconv (d e : NDist R) : NDist R := d.flatMap fun p ↦ nshift p.1 p.2 e

@[simp]
theorem nconv_nil (e : NDist R) : nconv [] e = [] := rfl

@[simp]
theorem nconv_cons (p : ℕ × R) (d e : NDist R) :
    nconv (p :: d) e = nshift p.1 p.2 e ++ nconv d e := rfl

theorem gen_append (d e : NDist R) (t : R) : gen (d ++ e) t = gen d t + gen e t := by
  simp [gen]

theorem gen_nshift (x : ℕ) (w t : R) (e : NDist R) :
    gen (nshift x w e) t = w * t ^ x * gen e t := by
  induction e with
  | nil => simp [nshift]
  | cons p e ih =>
    simp only [nshift, List.map_cons, gen_cons] at ih ⊢
    rw [ih, pow_add]; ring

/-- Generating functions multiply under independent sums: the volume reading
`(1 + a) * (1 + b)`. -/
theorem gen_nconv (d e : NDist R) (t : R) : gen (nconv d e) t = gen d t * gen e t := by
  induction d with
  | nil => simp
  | cons p d ih => rw [nconv_cons, gen_append, gen_nshift, ih, gen_cons]; ring

/-- The dial. If `a` and `b` are the `k`-moments of `X` and `Y`, read through
`1 + k * a = gen d (1 + k)`, then the `k`-moment of `X + Y` is `wAdd k a b`. -/
theorem kmoment_law (k a b : R) (d e : NDist R) (ha : 1 + k * a = gen d (1 + k))
    (hb : 1 + k * b = gen e (1 + k)) : 1 + k * wAdd k a b = gen (nconv d e) (1 + k) := by
  rw [gen_nconv, ← ha, ← hb, one_add_mul_wAdd]

/-- A coin on counts: `1` with weight `p`, `0` with weight `1 - p`. -/
def nbern (p : R) : NDist R := [(0, 1 - p), (1, p)]

/-- A coin's `k`-moment is `p` at every `k`: one number, every weight. -/
theorem bern_kmoment (k p : R) : gen (nbern p) (1 + k) = 1 + k * p := by
  simp [nbern]; ring

/-- Coins with different chances: the number of heads. -/
def heads (ps : List R) : NDist R := ps.foldr (fun p d ↦ nconv (nbern p) d) [(0, 1)]

/-- The weight `k` fold of a list. -/
def foldW (k : R) (ps : List R) : R := ps.foldr (wAdd k) 0

theorem foldW_neg_one (ps : List R) : foldW (-1) ps = unionL ps := rfl

theorem heads_kmoment (k : R) (ps : List R) :
    gen (heads ps) (1 + k) = 1 + k * foldW k ps := by
  induction ps with
  | nil => simp [heads, foldW]
  | cons p ps ih =>
    have h : heads (p :: ps) = nconv (nbern p) (heads ps) := rfl
    have h' : foldW k (p :: ps) = wAdd k p (foldW k ps) := rfl
    rw [h, h', gen_nconv, bern_kmoment, ih, one_add_mul_wAdd]

/-- `t = 0` (weight `-1`): `gen d 0 = P(no heads)`, the union rule's survival. -/
theorem heads_none (ps : List R) : gen (heads ps) 0 = 1 - foldW (-1) ps := by
  have h := heads_kmoment (-1) ps
  rw [add_neg_cancel] at h
  rw [h]; ring

/-- `t = -1` (weight `-2`): `gen d (-1) = P(even) - P(odd)`, piling-up for unequal coins. -/
theorem heads_parity (ps : List R) : gen (heads ps) (-1) = 1 - 2 * foldW (-2) ps := by
  have h := heads_kmoment (-2) ps
  have e : (1 : R) + -2 = -1 := by ring
  rw [e] at h
  rw [h]; ring

/-- The weight on even counts. -/
def pEven (d : NDist R) : R := (d.map fun p ↦ if Even p.1 then p.2 else 0).sum

/-- The weight on odd counts. -/
def pOdd (d : NDist R) : R := (d.map fun p ↦ if Even p.1 then 0 else p.2).sum

theorem gen_neg_one (d : NDist R) : gen d (-1) = pEven d - pOdd d := by
  induction d with
  | nil => simp [pEven, pOdd]
  | cons p d ih =>
    simp only [pEven, pOdd, List.map_cons, List.sum_cons, gen_cons] at ih ⊢
    rw [ih]
    by_cases h : Even p.1
    · simp only [h, h.neg_one_pow, ite_true]; ring
    · simp only [h, (Nat.not_even_iff_odd.1 h).neg_one_pow, ite_false]; ring

/-- The weight on the count `0`. -/
def pZero (d : NDist R) : R := (d.map fun p ↦ if p.1 = 0 then p.2 else 0).sum

theorem gen_zero (d : NDist R) : gen d 0 = pZero d := by
  induction d with
  | nil => simp [pZero]
  | cons p d ih =>
    simp only [pZero, List.map_cons, List.sum_cons, gen_cons] at ih ⊢
    rw [ih, zero_pow_eq, mul_ite, mul_one, mul_zero]

/-- `n` equal coins. -/
def binom (p : R) : ℕ → NDist R
  | 0 => [(0, 1)]
  | n + 1 => nconv (nbern p) (binom p n)

/-- `n` equal coins: the `k`-moment climbs, `gen (binom p n) (1 + k) = 1 + k * wpow k p n`,
which is `(1 + k * p) ^ n`. -/
theorem binom_kmoment (k p : R) (n : ℕ) : gen (binom p n) (1 + k) = 1 + k * wpow k p n := by
  induction n with
  | zero => simp [binom]
  | succ n ih =>
    rw [binom, gen_nconv, bern_kmoment, ih, wpow_succ, one_add_mul_wAdd, mul_comm]

/-- The `k`-moment counts lit corners: `(1 + k) ^ x = ∑ j, x.choose j * k ^ j`, so `gen d (1 + k)`
is the sum over outcomes of weight times `∑ j, x.choose j * k ^ j` (factorial moments). -/
def litCorners (d : NDist R) (k : R) : R :=
  (d.map fun p ↦ p.2 * ∑ j ∈ range (p.1 + 1), (p.1.choose j : R) * k ^ j).sum

theorem gen_eq_litCorners (d : NDist R) (k : R) : gen d (1 + k) = litCorners d k := by
  simp only [gen, litCorners, E.binomial]

/-! #### The binomial distribution is the weighted Pascal layers -/

/-- `bp p q n j` is the chance that `n` coins show `j` heads, built by adding one coin at a
time. -/
def bp (p q : R) : ℕ → ℕ → R
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, 0 => q * bp p q n 0
  | n + 1, j + 1 => q * bp p q n (j + 1) + p * bp p q n j

/-- Binomial is Pascal layers: `P(j heads of n) = n.choose j * p ^ j * q ^ (n - j)`. -/
theorem bp_layers (p q : R) : ∀ n j, bp p q n j = (n.choose j : R) * p ^ j * q ^ (n - j)
  | 0, 0 => by simp [bp]
  | 0, _ + 1 => by simp [bp]
  | n + 1, 0 => by simp [bp, bp_layers p q n 0, pow_succ']
  | n + 1, j + 1 => by
    rw [bp, bp_layers p q n (j + 1), bp_layers p q n j, Nat.choose_succ_succ',
      Nat.add_sub_add_right]
    push_cast
    rcases Nat.lt_or_ge n (j + 1) with h | h
    · rw [Nat.choose_eq_zero_of_lt h]; simp; ring
    · have e : n - j = n - (j + 1) + 1 := by lia
      rw [e]; ring

theorem bp_zero_col (p q : R) (n : ℕ) : bp p q n 0 = q ^ n := by
  simp [bp_layers]

theorem bp_above (p q : R) (n j : ℕ) (h : n < j) : bp p q n j = 0 := by
  simp [bp_layers, Nat.choose_eq_zero_of_lt h]

/-! ### 4. Moments and cumulants: Pascal rows, then the additive coordinates -/

/-- A distribution: a list of (value, weight). Independent sum is convolution. -/
abbrev Dist (R : Type*) := List (R × R)

/-- The `r`-th moment `E[X ^ r]`. `moment 0` is the total weight, `moment 1` the mean. -/
def moment (r : ℕ) (d : Dist R) : R := (d.map fun p ↦ p.2 * p.1 ^ r).sum

@[simp]
theorem moment_nil (r : ℕ) : moment r ([] : Dist R) = 0 := rfl

@[simp]
theorem moment_cons (r : ℕ) (p : R × R) (d : Dist R) :
    moment r (p :: d) = p.2 * p.1 ^ r + moment r d := rfl

/-- Shift every value by `x` and scale every weight by `w`. -/
def shift (x w : R) (e : Dist R) : Dist R := e.map fun p ↦ (x + p.1, w * p.2)

/-- The distribution of `X + Y` for independent `X`, `Y`. -/
def conv (d e : Dist R) : Dist R := d.flatMap fun p ↦ shift p.1 p.2 e

@[simp]
theorem conv_nil (e : Dist R) : conv [] e = [] := rfl

@[simp]
theorem conv_cons (p : R × R) (d e : Dist R) :
    conv (p :: d) e = shift p.1 p.2 e ++ conv d e := rfl

theorem moment_append (r : ℕ) (d e : Dist R) :
    moment r (d ++ e) = moment r d + moment r e := by
  simp [moment]

/-- Moments of a shifted distribution: the binomial theorem, row `r` of Pascal's triangle. -/
theorem moment_shift (r : ℕ) (x w : R) (e : Dist R) :
    moment r (shift x w e) =
      w * ∑ j ∈ range (r + 1), (r.choose j : R) * x ^ (r - j) * moment j e := by
  induction e with
  | nil => simp [shift]
  | cons p e ih =>
    simp only [shift, List.map_cons, moment_cons] at ih ⊢
    rw [ih, add_comm x p.1, add_pow]
    simp only [mul_add, sum_add_distrib, mul_sum]
    congr 1
    exact sum_congr rfl fun j _ ↦ by ring

/-- Moments of a sum of independent variables follow the Pascal rows:
`E[(X + Y) ^ r] = ∑ j, r.choose j * E[X ^ (r - j)] * E[Y ^ j]`. -/
theorem moment_conv (r : ℕ) (d e : Dist R) :
    moment r (conv d e) =
      ∑ j ∈ range (r + 1), (r.choose j : R) * moment (r - j) d * moment j e := by
  induction d with
  | nil => simp
  | cons p d ih =>
    rw [conv_cons, moment_append, moment_shift, ih]
    simp only [moment_cons, mul_add, add_mul, sum_add_distrib, mul_sum]
    congr 1
    exact sum_congr rfl fun j _ ↦ by ring

/-- Total weight multiplies: the grid of outcomes. -/
theorem moment_zero_conv (d e : Dist R) :
    moment 0 (conv d e) = moment 0 d * moment 0 e := by
  simp [moment_conv]

/-- The Pascal row `1 1`. -/
theorem moment_one_conv (d e : Dist R) :
    moment 1 (conv d e) = moment 1 d * moment 0 e + moment 0 d * moment 1 e := by
  simp [moment_conv, sum_range_succ]

/-- The Pascal row `1 2 1`. -/
theorem moment_two_conv (d e : Dist R) :
    moment 2 (conv d e) =
      moment 2 d * moment 0 e + 2 * moment 1 d * moment 1 e + moment 0 d * moment 2 e := by
  simp [moment_conv, sum_range_succ]

/-- The Pascal row `1 3 3 1`. -/
theorem moment_three_conv (d e : Dist R) :
    moment 3 (conv d e) = moment 3 d * moment 0 e + 3 * moment 2 d * moment 1 e
      + 3 * moment 1 d * moment 2 e + moment 0 d * moment 3 e := by
  norm_num [moment_conv, sum_range_succ, Nat.choose]

/-- The Pascal row `1 4 6 4 1`. -/
theorem moment_four_conv (d e : Dist R) :
    moment 4 (conv d e) = moment 4 d * moment 0 e + 4 * moment 3 d * moment 1 e
      + 6 * moment 2 d * moment 2 e + 4 * moment 1 d * moment 3 e + moment 0 d * moment 4 e := by
  norm_num [moment_conv, sum_range_succ, Nat.choose]

/-! Cumulants: the additive coordinates. -/

/-- The second cumulant, the variance. -/
def k2 (d : Dist R) : R := moment 2 d - moment 1 d ^ 2

/-- The third cumulant. -/
def k3 (d : Dist R) : R := moment 3 d - 3 * moment 2 d * moment 1 d + 2 * moment 1 d ^ 3

/-- The fourth cumulant. -/
def k4 (d : Dist R) : R :=
  moment 4 d - 4 * moment 3 d * moment 1 d - 3 * moment 2 d ^ 2
    + 12 * moment 2 d * moment 1 d ^ 2 - 6 * moment 1 d ^ 4

/-- Means add. -/
theorem mean_add (d e : Dist R) (hd : moment 0 d = 1) (he : moment 0 e = 1) :
    moment 1 (conv d e) = moment 1 d + moment 1 e := by
  rw [moment_one_conv, hd, he]; ring

/-- Variances add. -/
theorem var_add (d e : Dist R) (hd : moment 0 d = 1) (he : moment 0 e = 1) :
    k2 (conv d e) = k2 d + k2 e := by
  unfold k2; rw [moment_two_conv, moment_one_conv, hd, he]; ring

/-- Third cumulants add. -/
theorem k3_add (d e : Dist R) (hd : moment 0 d = 1) (he : moment 0 e = 1) :
    k3 (conv d e) = k3 d + k3 e := by
  unfold k3; rw [moment_three_conv, moment_two_conv, moment_one_conv, hd, he]; ring

/-- Fourth cumulants add. -/
theorem k4_add (d e : Dist R) (hd : moment 0 d = 1) (he : moment 0 e = 1) :
    k4 (conv d e) = k4 d + k4 e := by
  unfold k4
  rw [moment_four_conv, moment_three_conv, moment_two_conv, moment_one_conv, hd, he]; ring

/-! #### `n` independent copies and the shrinking shape terms -/

/-- The sum of `n` independent copies of `d`. -/
def convPow (d : Dist R) : ℕ → Dist R
  | 0 => [(0, 1)]
  | n + 1 => conv (convPow d n) d

theorem convPow_mass (d : Dist R) (hd : moment 0 d = 1) (n : ℕ) :
    moment 0 (convPow d n) = 1 := by
  induction n with
  | zero => simp [convPow]
  | succ n ih => rw [convPow, moment_zero_conv, ih, hd, one_mul]

theorem convPow_mean (d : Dist R) (hd : moment 0 d = 1) (n : ℕ) :
    moment 1 (convPow d n) = n * moment 1 d := by
  induction n with
  | zero => simp [convPow]
  | succ n ih => rw [convPow, mean_add _ _ (convPow_mass d hd n) hd, ih]; push_cast; ring

theorem convPow_k2 (d : Dist R) (hd : moment 0 d = 1) (n : ℕ) :
    k2 (convPow d n) = n * k2 d := by
  induction n with
  | zero => simp [convPow, k2]
  | succ n ih => rw [convPow, var_add _ _ (convPow_mass d hd n) hd, ih]; push_cast; ring

theorem convPow_k3 (d : Dist R) (hd : moment 0 d = 1) (n : ℕ) :
    k3 (convPow d n) = n * k3 d := by
  induction n with
  | zero => simp [convPow, k3]
  | succ n ih => rw [convPow, k3_add _ _ (convPow_mass d hd n) hd, ih]; push_cast; ring

theorem convPow_k4 (d : Dist R) (hd : moment 0 d = 1) (n : ℕ) :
    k4 (convPow d n) = n * k4 d := by
  induction n with
  | zero => simp [convPow, k4]
  | succ n ih => rw [convPow, k4_add _ _ (convPow_mass d hd n) hd, ih]; push_cast; ring

/-- Rescaling `X ↦ c * X`. -/
def scale (c : R) (d : Dist R) : Dist R := d.map fun p ↦ (c * p.1, p.2)

/-- The `r`-th moment scales as `c ^ r`. -/
theorem moment_scale (r : ℕ) (c : R) (d : Dist R) :
    moment r (scale c d) = c ^ r * moment r d := by
  induction d with
  | nil => simp [scale]
  | cons p d ih =>
    simp only [scale, List.map_cons, moment_cons] at ih ⊢
    rw [ih]; ring

/-- Cumulant `r` scales as `c ^ r`: the reason the shape terms die after standardising. -/
theorem k2_scale (c : R) (d : Dist R) : k2 (scale c d) = c ^ 2 * k2 d := by
  simp only [k2, moment_scale]; ring

theorem k3_scale (c : R) (d : Dist R) : k3 (scale c d) = c ^ 3 * k3 d := by
  simp only [k3, moment_scale]; ring

theorem k4_scale (c : R) (d : Dist R) : k4 (scale c d) = c ^ 4 * k4 d := by
  simp only [k4, moment_scale]; ring

/-! #### The coin -/

/-- A coin: `1` with weight `p`, `0` with weight `1 - p`. -/
def bern (p : R) : Dist R := [(0, 1 - p), (1, p)]

theorem moment_zero_bern (p : R) : moment 0 (bern p) = 1 := by
  simp [bern]

theorem moment_one_bern (p : R) : moment 1 (bern p) = p := by
  simp [bern]

theorem k2_bern (p : R) : k2 (bern p) = p * (1 - p) := by
  simp [bern, k2]; ring

theorem k3_bern (p : R) : k3 (bern p) = p * (1 - p) * (1 - 2 * p) := by
  simp [bern, k3]; ring

theorem k4_bern (p : R) : k4 (bern p) = p * (1 - p) * (1 - 6 * p + 6 * p * p) := by
  simp [bern, k4]; ring

/-- Binomial mean `n * p` and variance `n * (p * (1 - p))`, from additivity alone. -/
theorem binom_mean_var (p : R) (n : ℕ) :
    moment 1 (convPow (bern p) n) = n * p ∧ k2 (convPow (bern p) n) = n * (p * (1 - p)) := by
  rw [convPow_mean _ (moment_zero_bern p), convPow_k2 _ (moment_zero_bern p), moment_one_bern,
    k2_bern]
  exact ⟨rfl, rfl⟩

/-! ### The parts that divide -/

variable {K : Type*} [Field K]

/-- A coin's variance is half its self-parity: the chance two copies differ, halved. -/
theorem bern_var_parity [NeZero (2 : K)] (p : K) : k2 (bern p) = wAdd (-2) p p / 2 := by
  rw [k2_bern, eq_div_iff two_ne_zero]; unfold wAdd; ring

/-- Skewness squared, a scale-free shape number. -/
def skewSq (d : Dist K) : K := k3 d ^ 2 / k2 d ^ 3

/-- Excess kurtosis, a scale-free shape number. -/
def kurt (d : Dist K) : K := k4 d / k2 d ^ 2

theorem skewSq_scale (c : K) (hc : c ≠ 0) (d : Dist K) : skewSq (scale c d) = skewSq d := by
  rw [skewSq, skewSq, k2_scale, k3_scale]
  by_cases hk : k2 d = 0
  · simp [hk]
  · field_simp

theorem kurt_scale (c : K) (hc : c ≠ 0) (d : Dist K) : kurt (scale c d) = kurt d := by
  rw [kurt, kurt, k2_scale, k4_scale]
  by_cases hk : k2 d = 0
  · simp [hk]
  · field_simp

/-- The exact part of the central limit theorem, for skewness: for a sum of `n` independent
copies, `skewSq` falls as `1 / n`. -/
theorem clt_skew [CharZero K] (d : Dist K) (hd : moment 0 d = 1) (n : ℕ) (hn : n ≠ 0)
    (hk : k2 d ≠ 0) : skewSq (convPow d n) * n = skewSq d := by
  have hn' : (n : K) ≠ 0 := Nat.cast_ne_zero.2 hn
  rw [skewSq, skewSq, convPow_k2 d hd n, convPow_k3 d hd n]
  field_simp

/-- The exact part of the central limit theorem, for kurtosis: for a sum of `n` independent
copies, `kurt` falls as `1 / n`. Only the mean and variance survive, which is the normal law. -/
theorem clt_kurt [CharZero K] (d : Dist K) (hd : moment 0 d = 1) (n : ℕ) (hn : n ≠ 0)
    (hk : k2 d ≠ 0) : kurt (convPow d n) * n = kurt d := by
  have hn' : (n : K) ≠ 0 := Nat.cast_ne_zero.2 hn
  rw [kurt, kurt, convPow_k2 d hd n, convPow_k4 d hd n]
  field_simp

/-- Fair coins: kurtosis `-2 / n`, falling to the normal's `0`. -/
theorem fair_coin_kurt [CharZero K] (n : ℕ) (hn : n ≠ 0) :
    kurt (convPow (bern (1 / 2 : K)) n) * n = -2 := by
  have h2 : k2 (bern (1 / 2 : K)) = 1 / 4 := by rw [k2_bern]; norm_num
  rw [clt_kurt _ (moment_zero_bern _) n hn (by rw [h2]; norm_num), kurt, h2, k4_bern]
  norm_num

/-! ### 5. Bayes: odds multiply, so `odds - 1` mult-adds -/

/-- The odds of a probability: `p / (1 - p)`. -/
def odds (p : K) : K := p / (1 - p)

/-- Bayes' rule in odds form: posterior odds are prior odds times the likelihood ratio. -/
theorem bayes_odds (p lH lN : K) (hp1 : p ≠ 1) (hN : lN ≠ 0)
    (hden : p * lH + (1 - p) * lN ≠ 0) :
    odds (p * lH / (p * lH + (1 - p) * lN)) = odds p * (lH / lN) := by
  have h1 : 1 - p ≠ 0 := sub_ne_zero.2 hp1.symm
  have e : 1 - p * lH / (p * lH + (1 - p) * lN) = (1 - p) * lN / (p * lH + (1 - p) * lN) := by
    field_simp; ring
  rw [odds, odds, e]
  field_simp

/-- Two independent pieces of evidence:
`o * l₁ * l₂ - 1 = (o - 1) ✳ (l₁ - 1) ✳ (l₂ - 1)`. -/
theorem evidence_multAdd (o l1 l2 : R) : o * l1 * l2 - 1 = (o - 1) ✳ (l1 - 1) ✳ (l2 - 1) := by
  unfold multAdd; ring

/-- The same evidence `n` times: `o * l ^ n - 1 = (o - 1) ✳ mpow (l - 1) n`. -/
theorem repeated_evidence (o l : R) (n : ℕ) : o * l ^ n - 1 = (o - 1) ✳ mpow (l - 1) n := by
  rw [mpow_eq]; unfold multAdd; ring

/-- The complement is the way back: `odds (1 - p) - 1` is the mult-add inverse of
`odds p - 1`. -/
theorem complement_is_inverse (p : K) (hp0 : p ≠ 0) (hp1 : p ≠ 1) :
    odds (1 - p) - 1 = multAddInv (odds p - 1) := by
  have h1 : 1 - p ≠ 0 := sub_ne_zero.2 hp1.symm
  unfold odds multAddInv
  have e : 1 + (p / (1 - p) - 1) = p / (1 - p) := by ring
  rw [sub_sub_cancel, e]
  field_simp; ring

theorem complement_cancels (p : K) (hp0 : p ≠ 0) (hp1 : p ≠ 1) :
    (odds p - 1) ✳ (odds (1 - p) - 1) = 0 := by
  have h1 : 1 - p ≠ 0 := sub_ne_zero.2 hp1.symm
  rw [complement_is_inverse p hp0 hp1]
  apply multAdd_multAddInv
  rw [add_sub_cancel, odds]
  exact div_ne_zero hp0 h1

/-- Pooling two independent opinions from a flat prior. -/
def pool (p q : K) : K := p * q / (p * q + (1 - p) * (1 - q))

/-- In the coordinate `ρ = 2 * p - 1` pooling is velocity addition, the lens law with walls at
`±1` (certainty). -/
theorem pool_velocity [NeZero (2 : K)] (p q : K) (h : p * q + (1 - p) * (1 - q) ≠ 0) :
    2 * pool p q - 1 = ((2 * p - 1) + (2 * q - 1)) / (1 + (2 * p - 1) * (2 * q - 1)) := by
  have e : 1 + (2 * p - 1) * (2 * q - 1) = 2 * (p * q + (1 - p) * (1 - q)) := by ring
  have h2 : (2 : K) ≠ 0 := two_ne_zero
  rw [pool, e]
  field_simp; ring

/-- Pooling multiplies odds: the same law in the odds coordinate. -/
theorem pool_odds (p q : K) (hp : p ≠ 1) (hq : q ≠ 1) (h : p * q + (1 - p) * (1 - q) ≠ 0) :
    odds (pool p q) = odds p * odds q := by
  have hp' : 1 - p ≠ 0 := sub_ne_zero.2 hp.symm
  have hq' : 1 - q ≠ 0 := sub_ne_zero.2 hq.symm
  have e : 1 - p * q / (p * q + (1 - p) * (1 - q)) =
      (1 - p) * (1 - q) / (p * q + (1 - p) * (1 - q)) := by
    field_simp; ring
  rw [odds, odds, odds, pool, e]
  field_simp

/-! ### 6. The dependence line: lens laws with a wall at certainty -/

/-- The Ali-Mikhail-Haq copula: `P(both)` for marginals `u`, `v` and dependence `θ`. -/
def amh (θ u v : K) : K := u * v / (1 - θ * (1 - u) * (1 - v))

/-- The lens law at weight `k` and bend `q`. -/
def lens (k q a b : K) : K := (a + b + k * a * b) / (1 - q * a * b)

/-- Union under AMH dependence is a lens law at `k = -(1 + θ)`, `q = θ`: "neither happens"
follows the AMH copula, so "at least one happens" is
`(p + q - (1 + θ) * p * q) / (1 - θ * p * q)`. -/
theorem amh_union_is_lens (θ p q : K) (h : 1 - θ * p * q ≠ 0) :
    1 - amh θ (1 - p) (1 - q) = lens (-(1 + θ)) θ p q := by
  unfold amh lens
  rw [sub_sub_cancel, sub_sub_cancel, eq_div_iff h, sub_mul, div_mul_cancel₀ _ h]
  ring

/-- Its walls: `1 + k * x + q * x ^ 2 = (1 - x) * (1 - θ * x)`. One wall sits at certainty,
`x = 1`; the other at `1 / θ`. -/
theorem amh_walls (θ x : R) : 1 + (-(1 + θ)) * x + θ * (x * x) = (1 - x) * (1 - θ * x) := by
  ring

/-- The lens `x ↦ (1 - x) / (1 - θ * x)` turns the law into multiplication, so it is
associative: two "survivals" multiply in this coordinate. -/
theorem amh_lens_multiplies (θ a b : K) (h : 1 - θ * a * b ≠ 0) :
    (1 - lens (-(1 + θ)) θ a b) * ((1 - θ * a) * (1 - θ * b)) =
      (1 - a) * (1 - b) * (1 - θ * lens (-(1 + θ)) θ a b) := by
  have hL : lens (-(1 + θ)) θ a b * (1 - θ * a * b) = a + b + -(1 + θ) * a * b :=
    div_mul_cancel₀ _ h
  apply mul_right_cancel₀ h
  linear_combination (-((1 - θ * a) * (1 - θ * b) - θ * ((1 - a) * (1 - b)))) * hL

/-- `θ = 0`: independence, the union rule `wAdd (-1) a b`. -/
theorem amh_independent (a b : K) : lens (-(1 + 0)) 0 a b = wAdd (-1) a b := by
  simp [lens, wAdd]

/-- `θ = -1`: walls at `±1`, velocity addition, the same law as pooling opinions. -/
theorem amh_velocity (a b : K) : lens (-(1 + -1)) (-1) a b = (a + b) / (1 + a * b) := by
  simp [lens]

/-- `θ = 1`: a double wall at certainty, the border parabola `k ^ 2 = 4 * q` (the cusp). -/
theorem amh_cusp : (-(1 + (1 : R))) * (-(1 + 1)) = 4 * 1 := by ring

/-- Every `θ` lands on the line `k + q = -1`. -/
theorem amh_line (θ : R) : (-(1 + θ)) + θ = -1 := by ring

/-- Why the line. An "or" law must let certainty swallow everything: `lens k q a 1 = 1`. A lens
law does that exactly when `k + q = -1`, that is, when it has a wall at `1`. -/
theorem certainty_absorbs_iff (k q a : K) (ha : a ≠ 0) (hd : 1 - q * a * 1 ≠ 0) :
    lens k q a 1 = 1 ↔ k + q = -1 := by
  rw [lens, div_eq_one_iff_eq hd]
  constructor
  · intro h
    have h3 : a * (1 + k + q) = 0 := by linear_combination h
    rcases mul_eq_zero.1 h3 with h4 | h4
    · exact absurd h4 ha
    · linear_combination h4
  · intro h
    linear_combination a * h

/-- The AMH copula at `θ = 0` is independence: `P(both) = u * v`. -/
theorem amh_zero (u v : K) : amh 0 u v = u * v := by
  simp [amh]

end StatFused
