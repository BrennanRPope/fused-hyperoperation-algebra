/-
DRAFT (Nov 26, 2026). First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.SpinningRiver.
The spinning river: the cheap algebra behind spinning-river/README.md.

Legs of a cell, in the order (time, radial, azimuth). The map from the flat lattice's legs to the
cell's legs adds v·(time leg) and −v·α·(azimuth leg) to the radial leg: two shears, the dial's
flat middle. Builds on FusedAlgebra.Topics.JetFields (M3, M3.mul).
-/
import FusedAlgebra.Topics.JetFields

namespace CircleOp.SpinningRiver

open CircleOp
open CircleOp.JetFields
open Lean.Grind

section
variable {R : Type} [CommRing R]

/-- The cell map minus the identity: the radial row picks up v (from time) and −v·α (from azimuth). -/
def shear (v α : R) : M3 R := ⟨0, 0, 0, v, 0, -(v * α), 0, 0, 0⟩

/-- The cell map itself, I + shear. -/
def cellMap (v α : R) : M3 R := ⟨1, 0, 0, v, 1, -(v * α), 0, 0, 1⟩

/-- **The cell map is a dual number:** its shear part squares to zero. -/
theorem shear_sq (v α : R) : (shear v α).mul (shear v α) = ⟨0, 0, 0, 0, 0, 0, 0, 0, 0⟩ := by
  simp only [shear, M3.mul, M3.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- **Two cell maps compose by adding their shears** (boosts and turns at the flat middle add):
    (I + S₁)(I + S₂) = I + S₁ + S₂. -/
theorem cellMap_mul (v α w β : R) :
    (cellMap v α).mul (cellMap w β) = ⟨1, 0, 0, v + w, 1, -(v * α) - w * β, 0, 0, 1⟩ := by
  simp only [cellMap, M3.mul, M3.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- The cell map keeps volume: no layer but the 1 survives (P = Q = R = 0 for the shear). -/
theorem shear_layers (v α : R) :
    (shear v α).tr = 0 ∧ (shear v α).e2 = 0 ∧ (shear v α).det = 0 := by
  simp only [shear, M3.tr, M3.e2, M3.det, det9]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- **Kerr's circular-orbit period, cleared of roots.** With r = s², M = m²:
    m·(s⁶ − m²a²) = (m s³ − m² a)(s³ + m a), so T = 2π (s⁶ − m²a²)/(m s³ − m² a) = 2π (s³/m + a). -/
theorem kerr_period (s m a : R) :
    m * (s^6 - m^2 * a^2) = (m * s^3 - m^2 * a) * (s^3 + m * a) := by grind

end

end CircleOp.SpinningRiver
