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

/-! ### Round 7: deriving the turn (checks D1–D6 of spinning-river/checks.py) -/

/-- **Area law (D1a).** A shell turning at rate Ω has turn f = Ω (r²+b²) sin²θ; its clock gap 2π f is
    2 Ω times the area π (r²+b²) sin²θ of the latitude circle (π kept as a symbol). -/
theorem area_law (pi Ω r b s2 : R) :
    2 * pi * (Ω * (r^2 + b^2) * s2) = 2 * Ω * (pi * (r^2 + b^2) * s2) := by grind

/-- **Rim at light speed (D4).** The innermost shell turns at Ω₀ with Ω₀ b² = κ; if its rim speed Ω₀ b
    is 1, then b = κ. -/
theorem rim_at_light (κ b Ω0 : R) (hΩ : Ω0 * b^2 = κ) (hrim : Ω0 * b = 1) : b = κ := by grind

/-- **The t–φ block (D6a).** With V = v², V ρ² = 2 M r, ρ² = r² + a² c², c² + s² = 1:
    g_tφ² − g_tt g_φφ = s² (r² − 2 M r + a²), where g_tt = V − 1, g_tφ = −V a s², g_φφ = (r²+a²) s² + V a² s⁴. -/
theorem tphi_block (V ρ2 M r a c2 s2 : R)
    (hv : V * ρ2 = 2 * M * r) (hρ : ρ2 = r^2 + a^2 * c2) (hc : c2 + s2 = 1) :
    (V * a * s2)^2 - (V - 1) * ((r^2 + a^2) * s2 + V * a^2 * s2^2) = s2 * (r^2 - 2 * M * r + a^2) := by
  grind

/-- **The horizon turns at the shell rate (D6b).** Where r² − 2 M r + a² = 0, the dragging rate
    −g_tφ / g_φφ equals a / (r² + a²), written without division:
    (−g_tφ) (r² + a²) = a g_φφ. -/
theorem horizon_rate (V ρ2 M r a c2 s2 : R)
    (hv : V * ρ2 = 2 * M * r) (hρ : ρ2 = r^2 + a^2 * c2) (hc : c2 + s2 = 1)
    (hor : r^2 - 2 * M * r + a^2 = 0) :
    (V * a * s2) * (r^2 + a^2) = a * ((r^2 + a^2) * s2 + V * a^2 * s2^2) := by
  grind

end

end CircleOp.SpinningRiver
