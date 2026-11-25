# FusedAlgebra.Topics.Elliptic — review

Theorems (22): GENERAL 16, EQUIV 5 (frey_discriminant, frey_c4, frey_c4_at_A, frey_c4_at_C, gap6_zero_at_rungs), SAME 1 (gap6_whole). Definitions checked: d1–d5, gap6 SAME.

No WEAKER / DIFFERENT / MISSING / VACUOUS items, no low-confidence items.

Curve conventions checked against mathlib's WeierstrassCurve: Δ = -b₂²b₈ - 8b₄³ - 27b₆² + 9b₂b₄b₆ and c₄ = b₂² - 24b₄, with b₂ = a₁² + 4a₂, b₄ = 2a₄ + a₁a₃. For freyCurve = ⟨0, B-A, 0, -AB, 0⟩ this gives Δ = 16(a₂²a₄² - 4a₄³) and c₄ = 16a₂² - 48a₄, exactly the original's hand-written left sides. Same signs and scale (the original's "16 × cubic discriminant" is mathlib's Δ). The port's frey_* statements are stronger readings: they say the expression is mathlib's real Δ / c₄, where the original only stated the polynomial identity.

Route changes: none that matter. Identities are proved by `ring` instead of `grind`; gap6_ne_zero keeps the original's case split (s ≤ 0, 1, 2, ≥ 3) with nlinarith closing each case.

New port theorems with no original: 3 (nodalCurve_Δ, nodalCurve_c₄, freyCurve_equation; gap6_one is half of gap6_zero_at_rungs). Plus 2 new definitions (nodalCurve, freyCurve).
