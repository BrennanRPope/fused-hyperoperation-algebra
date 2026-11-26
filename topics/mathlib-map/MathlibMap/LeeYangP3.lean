/-
  STATEMENTS ONLY (DRAFT, from shape-dial-survey section "N3, p = 3").
  Sits next to LeeYangStep.lean and reuses `Holds` and `TopStrip`; built by CI with the Mathlib map.
  It defines two propositions and proves nothing. No `sorry`.

  RuleAtThree: Lee–Yang holds for the roundness-3 ball whenever 0 ≤ a ≤ 3 (dimension n = 3a + 1 ≤ 10 = 1 + p²).
    Proved in shape-dial-survey section "N3, p = 3" by: Pólya (a ≤ 1); a climb b → b + 1 for every b in [0, 2],
    where G at b + 1 alternates in sign at the zeros of G at b and Krein's partial-fraction theorem rules out
    complex zeros; the alternation is a computer certificate (p3_certificate.py): explicit inequalities for
    s ≥ 21 and error-bounded Taylor models on 200 boxes for s ≤ 25. Two constants there use floats with a
    safety margin, not ball arithmetic.
  Formalising it would need: Krein's theorem, Pólya's theorem, and the certificate's arithmetic in Lean.
-/
import MathlibMap.LeeYangStep

namespace MathlibMap.LeeYang

/-- The full holding side of the rule at roundness 3. -/
def RuleAtThree : Prop :=
  ∀ a : ℝ, 0 ≤ a → a ≤ 3 → Holds 3 a

/-- The p = 3 certificate covers the top strip 2 < a ≤ 3 (and with it, by `TopStripReduces`, the rest). -/
def StripAtThree : Prop := TopStrip 3

end MathlibMap.LeeYang
