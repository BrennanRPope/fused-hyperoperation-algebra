import Mathlib.Tactic.Ring

/-!
# The framework's setup for the plane: the two-dial law and its lift

`planeOp k d x y = (x + y + k·xy)/(1 − d·xy)`. With `k = α + β`, `d = αβ`, the map
`x ↦ (1 + αx)/(1 + βx)` turns it into multiplication (stated cross-multiplied).
-/

namespace HCT

variable {R : Type*} [CommRing R]

theorem plane_lift (α β x y : R) :
    let N := x + y + (α + β) * x * y
    let D := 1 - (α * β) * x * y
    (D + α * N) * ((1 + β * x) * (1 + β * y)) = ((1 + α * x) * (1 + α * y)) * (D + β * N) := by
  intro N D; ring

end HCT
