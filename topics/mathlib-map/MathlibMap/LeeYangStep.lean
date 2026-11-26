/-
  STATEMENTS ONLY (DRAFT, from shape-dial-survey section N3, step identity).
  Sits next to LeeYangBand.lean and reuses its definitions; built by CI with the Mathlib map.
  It defines three propositions and proves nothing. No `sorry`.

  StepDown: if Lee–Yang holds at a + 1 then it holds at a (any p > 0, a ≥ 0).
    Proved on paper in shape-dial-survey section N3, step identity:
    (a+1) p G_a = (s d/ds + (a+1) p + 1) G_(a+1); with u = s², the right side is
    2u φ'(u) + κ φ(u) for G_(a+1)(s) = φ(u); for u off the positive axis,
    ψ/φ = κ + 2 Σ u/(u − u_k) has nonzero imaginary part or is > κ, so ψ has only positive zeros.
  TopStripReduces: StepDown plus the strip p − 1 < a ≤ p gives the whole band a ≤ p
    (an induction on ⌈p − a⌉; purely logical once StepDown is in hand).
-/
import MathlibMap.LeeYangBand

namespace MathlibMap.LeeYang

/-- Lee–Yang passes from a + 1 down to a. -/
def StepDown : Prop :=
  ∀ p a : ℝ, 0 < p → 0 ≤ a → Holds p (a + 1) → Holds p a

/-- The top strip under a = p. -/
def TopStrip (p : ℝ) : Prop :=
  ∀ a : ℝ, p - 1 < a → a ≤ p → Holds p a

/-- With StepDown, the top strip is all that is left of the band at roundness p. -/
def TopStripReduces : Prop :=
  StepDown → ∀ p : ℝ, 1 ≤ p → TopStrip p → ∀ a : ℝ, 0 ≤ a → a ≤ p → Holds p a

end MathlibMap.LeeYang
