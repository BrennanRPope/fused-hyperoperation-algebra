/-
  Lee–Yang for p-balls: the shared core, refactored out of the p = 1, 2, 3 proofs (DRAFT).
  Sits next to MathlibMap/LeeYangBand.lean and reuses `fieldFn`, `Holds` and `EvenWhole`.

  What every proof so far has in common, written once:
    1. a base: Pólya (1918) gives Holds p a for a ≤ 1, at every p;
    2. a climb: if G at b has only real simple zeros and G at b + 1 alternates in sign at them, then G at b + 1 has
       only real zeros (Krein's partial-fraction theorem for Cartwright-class functions, plus G > 0 on the
       imaginary axis);
    3. a break: for p not even and a > p, some zero leaves the real line (Erdélyi endpoint expansion plus the
       Cartwright–Levinson count).
  Bundled as `Toolkit`, these are cited (not proved here). Everything after `Toolkit` is proved in Lean:
  the induction that turns base + climb into the whole band a ≤ p, the "iff" form of the rule, and each
  known roundness as an instance. What is left per roundness is exactly one statement, `AlternatesBelow p`.

  Instances:
    p = 1: needs nothing beyond the toolkit (the band a ≤ 1 is Pólya's base).
    p = 2: `AlternatesBelow 2` is the classical interlacing of Bessel zeros (Watson 1944, §15.22), cited.
    p = 3: `AlternatesBelow 3` is the computer certificate `p3_certificate.py`, cited.
    p = 4 (whole a) used a different tool, Laguerre's operator 1 − D²; it is not an instance of this core.
-/
import MathlibMap.LeeYangBand
-- The full library, so that `deriv`, `⌈·⌉₊` and the tactics used below are in scope (added in the repo).
import Mathlib

namespace MathlibMap.LeeYang

/-- G at b, read on the real line (it is real there). -/
noncomputable def G (p b t : ℝ) : ℝ := (fieldFn p b (t : ℂ)).re

/-- The per-roundness input: the positive zeros of G at b are simple, and G at b + 1 has opposite signs at any two
consecutive ones. -/
def Alternates (p b : ℝ) : Prop :=
  (∀ t : ℝ, 0 < t → G p b t = 0 → deriv (G p b) t ≠ 0) ∧
  ∀ s t : ℝ, 0 < s → s < t → G p b s = 0 → G p b t = 0 → (∀ u : ℝ, s < u → u < t → G p b u ≠ 0) →
    G p (b + 1) s * G p (b + 1) t < 0

/-- Alternation at every step of the climb that Pólya's base doesn't already cover. -/
def AlternatesBelow (p : ℝ) : Prop :=
  ∀ b : ℝ, 0 < b → b ≤ p - 1 → Alternates p b

/-- The cited analytic facts every roundness shares. -/
structure Toolkit : Prop where
  /-- Pólya 1918: a sine transform of a positive nondecreasing weight on (0, 1) has only real zeros. -/
  polya : ∀ p a : ℝ, 1 ≤ p → 0 ≤ a → a ≤ 1 → Holds p a
  /-- Krein: real simple zeros at b plus alternation of b + 1 there give real (interlacing) zeros at b + 1. -/
  climb : ∀ p b : ℝ, 1 ≤ p → 0 < b → Holds p b → Alternates p b → Holds p (b + 1)
  /-- The break past a = p, for p not an even whole number. -/
  breaks : ∀ p a : ℝ, 1 ≤ p → ¬ EvenWhole p → p < a → ¬ Holds p a

/-- The core theorem: base + climb + alternation below p give Lee–Yang on the whole band 0 ≤ a ≤ p. -/
theorem holds_of_alternatesBelow (T : Toolkit) {p : ℝ} (hp : 1 ≤ p) (halt : AlternatesBelow p) :
    ∀ a : ℝ, 0 ≤ a → a ≤ p → Holds p a := by
  have key : ∀ n : ℕ, ∀ a : ℝ, 0 ≤ a → a ≤ p → a ≤ (n : ℝ) + 1 → Holds p a := by
    intro n
    induction n with
    | zero =>
      intro a ha0 _ ha1
      exact T.polya p a hp ha0 (by simpa using ha1)
    | succ n ih =>
      intro a ha0 hap han
      by_cases h1 : a ≤ 1
      · exact T.polya p a hp ha0 h1
      · push_neg at h1
        push_cast at han
        have hb : Holds p (a - 1) := ih (a - 1) (by linarith) (by linarith) (by linarith)
        have hc := T.climb p (a - 1) hp (by linarith) hb (halt (a - 1) (by linarith) (by linarith))
        simpa using hc
  intro a ha0 hap
  exact key ⌈a⌉₊ a ha0 hap (by have := Nat.le_ceil a; linarith)

/-- The rule at one roundness, in "iff" form, for p not an even whole number. -/
theorem rule_of_alternatesBelow (T : Toolkit) {p : ℝ} (hp : 1 ≤ p) (hne : ¬ EvenWhole p)
    (halt : AlternatesBelow p) : ∀ a : ℝ, 0 ≤ a → (Holds p a ↔ a ≤ p) := by
  intro a ha0
  constructor
  · intro h
    by_contra hlt
    exact T.breaks p a hp hne (lt_of_not_le hlt) h
  · intro hap
    exact holds_of_alternatesBelow T hp halt a ha0 hap

/-- p = 1 needs nothing past the toolkit: there is no climb to make. -/
theorem holds_one (T : Toolkit) : ∀ a : ℝ, 0 ≤ a → a ≤ 1 → Holds 1 a :=
  holds_of_alternatesBelow T le_rfl (fun b hb0 hb1 => by exfalso; linarith)

/-- p = 2: the input is the interlacing of Bessel zeros (G at b is a Bessel function J of order b + 1/2,
divided by a power of s). Cited, Watson 1944 §15.22. -/
def BesselInterlacing : Prop := AlternatesBelow 2

theorem holds_two (T : Toolkit) (hB : BesselInterlacing) : ∀ a : ℝ, 0 ≤ a → a ≤ 2 → Holds 2 a :=
  holds_of_alternatesBelow T (by norm_num) hB

/-- p = 3: the input is the computer certificate `p3_certificate.py` (far zone by interval inequalities from the
exact contour formula, near zone by 200 boxes in b). Cited. -/
def Certificate3 : Prop := AlternatesBelow 3

theorem holds_three (T : Toolkit) (hC : Certificate3) : ∀ a : ℝ, 0 ≤ a → a ≤ 3 → Holds 3 a :=
  holds_of_alternatesBelow T (by norm_num) hC

theorem rule_three (T : Toolkit) (hC : Certificate3) (hne : ¬ EvenWhole 3) :
    ∀ a : ℝ, 0 ≤ a → (Holds 3 a ↔ a ≤ 3) :=
  rule_of_alternatesBelow T (by norm_num) hne hC

/-- 3 is not an even whole number. -/
theorem not_evenWhole_three : ¬ EvenWhole 3 := by
  rintro ⟨m, hm⟩
  have h : (3 : ℝ) = 2 * (m : ℝ) := hm
  rcases Nat.even_or_odd m with ⟨k, hk⟩ | ⟨k, hk⟩
  · have : (3 : ℝ) = 4 * (k : ℝ) := by rw [h, hk]; push_cast; ring
    have h' : (3 : ℕ) = 4 * k := by exact_mod_cast this
    omega
  · have : (3 : ℝ) = 4 * (k : ℝ) + 2 := by rw [h, hk]; push_cast; ring
    have h' : (3 : ℕ) = 4 * k + 2 := by exact_mod_cast this
    omega

/-- The goal of the refactor, as one statement: alternation below p at every roundness p ≥ 1 would give the
holding half of the rule for every p at once. -/
theorem holding_half_of_alternation (T : Toolkit) (hall : ∀ p : ℝ, 1 ≤ p → AlternatesBelow p) :
    ∀ p a : ℝ, 1 ≤ p → 0 ≤ a → a ≤ p → Holds p a :=
  fun p a hp ha0 hap => holds_of_alternatesBelow T hp (hall p hp) a ha0 hap

end MathlibMap.LeeYang
