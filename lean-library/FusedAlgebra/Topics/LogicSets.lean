/-
  Logic and set theory through mult-add
  =====================================

  Companion page and checks: ../index.html, ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base and FusedAlgebra.Rungs.
  Intended module path: FusedAlgebra/Topics/LogicSets.lean.

  Builds on the statistics page (truth tables, De Morgan, inclusion–exclusion are in
  Topics/Statistics and are not repeated here). What is new here:

    1. Why only OR and XOR: for idempotents a, b (a² = a, b² = b) the dial
       a ⊞ₖ b is idempotent again up to the defect (k + 1)(k + 2)·ab. The two roots are
       k = −1 (OR, the Boolean lattice) and k = −2 (XOR, the Boolean ring). In q-number
       terms q = k + 1 is 0 or −1, and n trues fold to [n]_0 = 1 or [n]_{−1} = n mod 2.
    2. Distributivity and absorption fail by a defect k·abc(1 − a) and ab(1 − a): both
       vanish exactly on idempotents. So the Boolean laws are the dial restricted to 0, 1.
    3. Walls as truth values: OR's wall is 1 (true); XOR's wall is 1/2, the fair coin
       (the one-time pad).
    4. Fuzzy logic: Hamacher's t-conorms are the landscape line p + q = −1, and in the
       odds coordinate x/(1 − x) every one of them is the dial with k = γ. The t-norms are
       their De Morgan mirrors.
    5. Logic of the rungs: the fused exponent alone defines mult-add, shifted add and 0
       (so it carries all of arithmetic), while mult-add alone is multiplication in
       disguise (decidable, Skolem/Mostowski).
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/LogicSets.lean; merged unchanged; no exact twins.
-/
import FusedAlgebra.Base
import FusedAlgebra.Rungs

namespace LogicSets
open CircleOp Lean.Grind

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The dial on idempotents: only k = −1 and k = −2 keep truth values -/

section Ring
variable {R : Type} [CommRing R]

/-- The dial's table on 0 and 1, for every weight k. Only the corner moves: 1 ⊞ₖ 1 = 2 + k. -/
theorem dial_table (k : R) :
    circK k (0:R) 0 = 0 ∧ circK k (0:R) 1 = 1 ∧ circK k (1:R) 0 = 1 ∧ circK k (1:R) 1 = 2 + k := by
  unfold circK; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- **The idempotent defect.** For a² = a and b² = b,
    (a ⊞ₖ b)² − (a ⊞ₖ b) = (k + 1)(k + 2)·ab. -/
theorem idem_defect (k a b : R) (ha : a * a = a) (hb : b * b = b) :
    circK k a b * circK k a b - circK k a b = (k + 1) * (k + 2) * (a * b) := by
  unfold circK
  have e : (a + b + k * a * b) * (a + b + k * a * b)
      = a * a + b * b + 2 * (a * b) + 2 * k * (a * a) * b + 2 * k * a * (b * b)
        + k * k * (a * a) * (b * b) := by grind
  rw [e, ha, hb]; grind

/-- At k = −1 (OR) idempotents stay idempotent. -/
theorem or_keeps_idem (a b : R) (ha : a * a = a) (hb : b * b = b) :
    circK (-1) a b * circK (-1) a b = circK (-1) a b := by
  have h := idem_defect (-1) a b ha hb; grind

/-- At k = −2 (XOR) idempotents stay idempotent. -/
theorem xor_keeps_idem (a b : R) (ha : a * a = a) (hb : b * b = b) :
    circK (-2) a b * circK (-2) a b = circK (-2) a b := by
  have h := idem_defect (-2) a b ha hb; grind

/-- AND (the product, the corner cell) keeps idempotents for free. -/
theorem and_keeps_idem (a b : R) (ha : a * a = a) (hb : b * b = b) :
    (a * b) * (a * b) = a * b := by
  have : (a * b) * (a * b) = (a * a) * (b * b) := by grind
  rw [this, ha, hb]

/-- The distributivity defect: a(b ⊞ₖ c) − (ab ⊞ₖ ac) = k·abc(1 − a). -/
theorem dist_defect (k a b c : R) :
    a * circK k b c - circK k (a * b) (a * c) = k * a * b * c * (1 - a) := by
  unfold circK; grind

/-- On an idempotent, AND distributes over every dial value (so over OR and over XOR). -/
theorem dist_on_idem (k a b c : R) (ha : a * a = a) :
    a * circK k b c = circK k (a * b) (a * c) := by
  have h := dist_defect k a b c
  have : k * a * b * c * (1 - a) = k * b * c * (a - a * a) := by grind
  rw [ha] at this; grind

/-- The absorption defect: a ⊞₋₁ (ab) = a + ab(1 − a). -/
theorem absorb_defect (a b : R) : circK (-1) a (a * b) = a + a * b * (1 - a) := by
  unfold circK; grind

/-- Absorption, both ways, on idempotents: a ∨ (a ∧ b) = a and a ∧ (a ∨ b) = a. -/
theorem absorb_on_idem (a b : R) (ha : a * a = a) :
    circK (-1) a (a * b) = a ∧ a * circK (-1) a b = a := by
  constructor
  · rw [absorb_defect]
    have : a * b * (1 - a) = b * (a - a * a) := by grind
    rw [this, ha]; grind
  · unfold circK
    have : a * (a + b + -1 * a * b) = a * a + b * (a - a * a) := by grind
    rw [this, ha]; grind

/-- OR is idempotent exactly on idempotents: a ⊞₋₁ a − a = a − a². -/
theorem or_self (a : R) : circK (-1) a a - a = a - a * a := by
  unfold circK; grind

/-- XOR of a thing with itself: a ⊞₋₂ a = 2(a − a²), zero on idempotents (self-inverse). -/
theorem xor_self (a : R) : circK (-2) a a = 2 * (a - a * a) := by
  unfold circK; grind

theorem xor_self_idem (a : R) (ha : a * a = a) : circK (-2) a a = 0 := by
  rw [xor_self, ha]; grind

/-- OR is XOR corrected by AND: a ∨ b = (a ⊕ b) ⊕ (a ∧ b) on idempotents. -/
theorem or_from_xor_and (a b : R) (ha : a * a = a) (hb : b * b = b) :
    circK (-1) a b = circK (-2) (circK (-2) a b) (a * b) := by
  unfold circK
  have e : a + b + -2 * a * b + a * b + -2 * (a + b + -2 * a * b) * (a * b)
      = a + b - a * b - 2 * ((a * a) * b) - 2 * (a * (b * b)) + 4 * ((a * a) * (b * b)) := by grind
  rw [e, ha, hb]; grind

/-- XOR from OR and AND: a ⊕ b = (a ∨ b)·(1 − ab) on idempotents. -/
theorem xor_from_or_and (a b : R) (ha : a * a = a) (hb : b * b = b) :
    circK (-2) a b = circK (-1) a b * (1 - a * b) := by
  unfold circK
  have e : (a + b + -1 * a * b) * (1 - a * b)
      = a + b - a * b - (a * a) * b - a * (b * b) + (a * a) * (b * b) := by grind
  rw [e, ha, hb]; grind

/-- OR's wall is "true": 1 absorbs everything under ⊞₋₁. -/
theorem or_wall (b : R) : circK (-1) 1 b = 1 := by unfold circK; grind

/-- The order of the Boolean lattice, read on the dial: on idempotents
    a ⊞₋₁ b = b exactly when a·b = a. -/
theorem order_iff (a b : R) :
    circK (-1) a b - b = a - a * b := by
  unfold circK; grind

end Ring

/-- On the integers: 1 ⊞ₖ 1 is a truth value exactly at k = −1 and k = −2. -/
theorem bits_closed_iff (k : Int) :
    (circK k (1:Int) 1 = 0 ∨ circK k (1:Int) 1 = 1) ↔ (k = -2 ∨ k = -1) := by
  unfold circK; constructor <;> intro h <;> omega

/-- n trues folded with OR give 1: [n]_q at q = 0. -/
theorem fold_or_trues (n : Nat) : wpow (-1 : Int) 1 (n + 1) = 1 := by
  induction n with
  | zero => simp [wpow, circK]
  | succ n ih => simp only [wpow] at *; rw [ih]; unfold circK; omega

/-- n trues folded with XOR give n mod 2: [n]_q at q = −1. -/
theorem fold_xor_trues (n : Nat) : wpow (-2 : Int) 1 n = (n % 2 : Nat) := by
  induction n with
  | zero => simp [wpow]
  | succ n ih =>
    simp only [wpow] at *; rw [ih]; unfold circK
    rcases Nat.mod_two_eq_zero_or_one n with h | h
    · have : (n + 1) % 2 = 1 := by omega
      rw [h, this]; decide
    · have : (n + 1) % 2 = 0 := by omega
      rw [h, this]; decide

/-! ## 2. XOR's wall is the fair coin: the one-time pad -/

/-- XOR's wall sits at 1/2: a fair key bit hides any message bit (Shannon's one-time pad).
    This is the library's `circK_wall` at k = −2, read as a truth value. -/
theorem one_time_pad (m : Rat) : circK (-2) m (1 / 2) = 1 / 2 := by
  have h := circK_wall (-2 : Rat) m (by decide)
  have e : (-1 : Rat) / (-2) = 1 / 2 := by grind
  rw [e] at h; exact h

/-! ## 3. Fuzzy logic: Hamacher's family on the landscape -/

/-- A landscape law F_{p,q}(a, b) = (a + b + p·ab)/(1 − q·ab). -/
def lens (p q a b : Rat) : Rat := (a + b + p * a * b) / (1 - q * a * b)

/-- Hamacher's t-conorm with parameter γ (Hamacher 1978). -/
def hamS (γ a b : Rat) : Rat := (a + b + (γ - 2) * a * b) / (1 + (γ - 1) * a * b)

/-- Hamacher's t-norm with parameter γ. -/
def hamT (γ a b : Rat) : Rat := a * b / (γ + (1 - γ) * (a + b - a * b))

/-- Odds: x/(1 − x). -/
def odds (x : Rat) : Rat := x / (1 - x)

/-- Every Hamacher t-conorm is a landscape law, at p = γ − 2, q = 1 − γ. -/
theorem hamS_is_lens (γ a b : Rat) : hamS γ a b = lens (γ - 2) (1 - γ) a b := by
  unfold hamS lens
  have : 1 + (γ - 1) * a * b = 1 - (1 - γ) * a * b := by grind
  rw [this]

/-- And they all sit on the line p + q = −1. -/
theorem hamacher_line (γ : Rat) : (γ - 2) + (1 - γ) = -1 := by grind

/-- γ = 1 is the probabilistic sum, the dial at k = −1 (no bend). -/
theorem ham_one (a b : Rat) : hamS 1 a b = circK (-1) a b := by
  unfold hamS circK; grind

/-- γ = 2 is the Einstein sum (a + b)/(1 + ab): relativity's point (0, −1). -/
theorem ham_two (a b : Rat) : hamS 2 a b = lens 0 (-1) a b := by
  unfold hamS lens; grind

/-- γ = 0 sits at (−2, 1), on the border parabola p² = 4q (a double wall). -/
theorem ham_zero_cusp : hamS 0 = lens (-2) 1 ∧ (-2 : Rat) * (-2) = 4 * 1 := by
  refine ⟨?_, by grind⟩
  funext a b; unfold hamS lens; grind

/-- The miss of a Hamacher t-conorm does not depend on γ:
    1 − S_γ(a, b) = (1 − a)(1 − b)/(1 + (γ − 1)ab). -/
theorem hamS_miss (γ a b : Rat) (h : 1 + (γ - 1) * a * b ≠ 0) :
    1 - hamS γ a b = (1 - a) * (1 - b) / (1 + (γ - 1) * a * b) := by
  unfold hamS; grind

/-- **Hamacher is the dial on odds.** odds(S_γ(a, b)) = odds a ⊞_γ odds b. -/
theorem hamS_odds (γ a b : Rat) (ha : 1 - a ≠ 0) (hb : 1 - b ≠ 0)
    (h : 1 + (γ - 1) * a * b ≠ 0) :
    odds (hamS γ a b) = circK γ (odds a) (odds b) := by
  unfold odds
  rw [hamS_miss γ a b h]
  have hab : (1 - a) * (1 - b) ≠ 0 := by
    intro e
    have : (1 - a) * (1 - b) / (1 - b) = 0 := by rw [e]; grind
    have : (1 - a) * (1 - b) / (1 - b) = 1 - a := by grind
    grind
  have s1 : hamS γ a b / ((1 - a) * (1 - b) / (1 + (γ - 1) * a * b))
      = (a + b + (γ - 2) * a * b) / ((1 - a) * (1 - b)) := by
    unfold hamS; grind
  have s2 : circK γ (a / (1 - a)) (b / (1 - b))
      = (a + b + (γ - 2) * a * b) / ((1 - a) * (1 - b)) := by
    unfold circK
    have e1 : a / (1 - a) = a * (1 - b) / ((1 - a) * (1 - b)) := by grind
    have e2 : b / (1 - b) = b * (1 - a) / ((1 - a) * (1 - b)) := by grind
    rw [e1, e2]; grind
  rw [s1, s2]

/-- De Morgan: Hamacher's t-norm is the mirror of its t-conorm, T(a, b) = 1 − S(1 − a, 1 − b). -/
theorem hamT_mirror (γ a b : Rat) (h : γ + (1 - γ) * (a + b - a * b) ≠ 0) :
    hamT γ a b = 1 - hamS γ (1 - a) (1 - b) := by
  unfold hamT hamS
  have e : 1 + (γ - 1) * (1 - a) * (1 - b) = γ + (1 - γ) * (a + b - a * b) := by grind
  rw [e]; grind

/-- On truth values every Hamacher t-conorm with γ ≠ 0 is OR. -/
theorem hamS_bits (γ : Rat) (hγ : γ ≠ 0) :
    hamS γ 0 0 = 0 ∧ hamS γ 0 1 = 1 ∧ hamS γ 1 0 = 1 ∧ hamS γ 1 1 = 1 := by
  unfold hamS
  refine ⟨by grind, by grind, by grind, ?_⟩
  have : (1 : Rat) + 1 + (γ - 2) * 1 * 1 = γ := by grind
  have h2 : (1 : Rat) + (γ - 1) * 1 * 1 = γ := by grind
  rw [this, h2]; grind

/-! ## 4. Logic of the rungs: the fused exponent carries arithmetic -/

/-- The fused exponent at base 1 is 2^(n+1) − 1, so it tells exponents apart. -/
theorem fexp_one_inj (m n : Nat) : fexp 1 m = fexp 1 n ↔ m = n := by
  constructor
  · intro h
    have hm := fexp_lift 1 m
    have hn := fexp_lift 1 n
    have : (1 + 1) ^ (m + 1) = (1 + 1) ^ (n + 1) := by omega
    have := (Nat.pow_right_inj (a := 2) (by decide)).mp this
    omega
  · intro h; rw [h]

/-- Mult-add is definable from the fused exponent alone:
    b ⊞ c = d exactly when (a ⊛ b) ⊛ c = a ⊛ d for every a (the tower law). -/
theorem circ_definable (b c d : Nat) :
    b ◯ c = d ↔ ∀ a, fexp (fexp a b) c = fexp a d := by
  constructor
  · intro h a; rw [fexp_tower, h]
  · intro h
    have := h 1
    rw [fexp_tower] at this
    exact (fexp_one_inj _ _).mp this

/-- Shifted add (rung 0) is definable from the fused exponent and mult-add:
    b + c + 1 = d exactly when (a ⊛ b) ⊞ (a ⊛ c) = a ⊛ d for every a. -/
theorem sadd_definable (b c d : Nat) :
    sadd b c = d ↔ ∀ a, (fexp a b) ◯ (fexp a c) = fexp a d := by
  constructor
  · intro h a; rw [← fexp_sadd, h]
  · intro h
    have := h 1
    rw [← fexp_sadd] at this
    exact (fexp_one_inj _ _).mp this

/-- Zero is definable from the fused exponent: the only right identity. -/
theorem zero_definable (y : Nat) : (∀ a, fexp a y = a) ↔ y = 0 := by
  constructor
  · intro h
    have h1 := h 1
    rw [← fexp_zero_right 1] at h1
    exact (fexp_one_inj _ _).mp h1
  · intro h a; rw [h, fexp_zero_right]

/-- With 0 in hand the successor is shifted add with 0, and plain + is shifted add, shifted back. -/
theorem succ_and_add (x y z : Nat) : sadd x 0 = x + 1 ∧ (x + y = z ↔ sadd x y = z + 1) := by
  unfold sadd; constructor <;> omega

/-- Mult-add alone is multiplication in disguise: the lift x ↦ x + 1 carries (Nat, ⊞, 0)
    onto the positive numbers under ×, one-to-one. -/
theorem lift_iso (a b : Nat) : (a ◯ b) + 1 = (a + 1) * (b + 1) ∧ (a + 1 = b + 1 → a = b) :=
  ⟨circ_lift' a b, fun h => by omega⟩

end LogicSets
