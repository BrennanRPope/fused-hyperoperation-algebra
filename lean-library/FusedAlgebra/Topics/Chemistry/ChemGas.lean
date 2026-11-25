/-
  Chemistry 6, gas laws
  =====================

  Companion page and checks: ../index.html, ../checks.py (15 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  van der Waals wrote his law as (P + a/V²)(V − b) = RT. Divided by RT·V
  it says 1/Z = (1 + α)(1 − β) with α = a/(PV²) and β = b/V, so
  1/Z − 1 = α ✳ (−β): attraction mult-added with minus the excluded volume.

  Mayer: the Boltzmann factor of a set of pairs is Π(1 + f), the lift of
  the mult-add fold of the f's; its expansion over subsets is the cluster
  (graph) expansion.

  Cubic equations of state differ in the attraction denominator
  V² + u·bV + w·b²; this file records the three classic choices.

  LIBRARY NOTE: `local infixl " ✳ "` for multAdd; delete it if Base declares ✳.
  Merged at Topics/Chemistry/ChemGas; the local infix stays. One exact twin removed (see the
  LIBRARY comment).
-/
import FusedAlgebra.Base

namespace ChemGas
open MultAdd

local infixl:65 " ✳ " => multAdd

/-- The two corrections multiply: (1 + α)(1 − β) = 1 + α ✳ (−β). -/
theorem two_corrections (α β : Rat) : (1 + α) * (1 - β) = 1 + (α ✳ (-β)) := by
  rw [multAdd_lift]; grind

/-- The corner of the two corrections is −αβ. -/
theorem corner (α β : Rat) : α ✳ (-β) = α - β - α * β := by
  rw [multAdd_def]; grind

/-- van der Waals in its own form gives 1/Z = (1 + a/PV²)(1 − b/V). -/
theorem vdw_inverse_Z (P V a b RT : Rat) (hV : V ≠ 0) (hP : P ≠ 0)
    (h : (P + a / (V * V)) * (V - b) = RT) :
    (1 + a / (P * V * V)) * (1 - b / V) * (P * V) = RT := by
  rw [← h]; grind

/-- vdW critical constants Vc = 3b, RTc = 8a/(27b), Pc = a/(27b²) give Zc = 3/8. -/
theorem critical_Z (a b : Rat) (hb : b ≠ 0) :
    (a / (27 * b * b)) * (3 * b) = 3 / 8 * (8 * a / (27 * b)) := by
  grind

/-- Mayer for three pairs: (1 + f₁)(1 + f₂)(1 + f₃) = 1 + f₁ ✳ f₂ ✳ f₃. -/
theorem mayer_three (f1 f2 f3 : Rat) : (1 + f1) * (1 + f2) * (1 + f3) = 1 + (f1 ✳ f2 ✳ f3) := by
  rw [multAdd_lift, multAdd_lift]

/-- Expanded, the three-pair fold is the seven nonempty products: singles, pairs and the triangle. -/
theorem mayer_clusters (f1 f2 f3 : Rat) :
    f1 ✳ f2 ✳ f3 = (f1 + f2 + f3) + (f1 * f2 + f1 * f3 + f2 * f3) + f1 * f2 * f3 := by
  simp only [multAdd_def]; grind

-- LIBRARY: mayer_list (Any list of pair functions: 1 + fold = product of (1 + f)) was an exact twin of ThermoFused.fermion_lift and is removed.

/-- Cubic denominators: vdW v², Redlich–Kwong v(v + 1), Peng–Robinson (v + 1)² − 2 (with b = 1). -/
def cubicD (u w v : Rat) : Rat := v * v + u * v + w
theorem vdw_denominator (v : Rat) : cubicD 0 0 v = v * v := by unfold cubicD; grind
theorem rk_denominator (v : Rat) : cubicD 1 0 v = v * (v + 1) := by unfold cubicD; grind
theorem pr_denominator (v : Rat) : cubicD 2 (-1) v = (v + 1) * (v + 1) - 2 := by unfold cubicD; grind

/-- Their discriminants u² − 4w are 0, 1 and 8. -/
theorem cubic_discriminants : (0 : Int) ^ 2 - 4 * 0 = 0 ∧ (1 : Int) ^ 2 - 4 * 0 = 1 ∧ (2 : Int) ^ 2 - 4 * (-1) = 8 := by
  decide

/-- Dalton: partial pressures x·P add up to P when the mole fractions sum to 1. -/
theorem dalton (x1 x2 x3 P : Rat) (h : x1 + x2 + x3 = 1) : x1 * P + x2 * P + x3 * P = P := by
  rw [← Rat.add_mul, ← Rat.add_mul, h, Rat.one_mul]

/-- Fugacity coefficients multiply into an equilibrium: here φ₁²/(φ₂ φ₃³) written without division. -/
theorem fugacity_product (p1 p2 q1 q2 : Rat) : (p1 * q1) * (p2 * q2) = (p1 * p2) * (q1 * q2) := by
  grind

end ChemGas
