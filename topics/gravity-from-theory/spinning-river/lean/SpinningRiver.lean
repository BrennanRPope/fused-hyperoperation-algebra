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

/-! ### Round 8: deriving the speed law (checks S1–S3) -/

/-- **The complex distance (S1a).** With an element i, i·i = −1, and a lattice point
    x² + y² = (r² + a²) s², z = r c, c² + s² = 1: the squared distance to the source moved up the
    axis by i·a is (r − i a c)². -/
theorem complex_distance (i r a c s2 : R) (hi : i * i = -1) (hc : c^2 + s2 = 1) :
    (r^2 + a^2) * s2 + (r * c - i * a)^2 = (r - i * a * c)^2 := by grind

/-- **Its real part gives the speed law (S2a).** (r − i a c)(r + i a c) = r² + a² c² = ρ², so
    M/(r − i a c) = M (r + i a c)/ρ², whose real part is M r/ρ². -/
theorem conj_product (i r a c : R) (hi : i * i = -1) :
    (r - i * a * c) * (r + i * a * c) = r^2 + a^2 * c^2 := by grind

/-! ### Round 9: the shifted source's rays (check T1) -/

/-- **The ray has unit length (T1).** For real p, q with k = p + q × p (before dividing by |p|²):
    |p + q × p|² = (p·p)(1 + q·q) − (p·q)². With p·p − q·q = 1 and p·q = 0 (E = p + iq a complex unit
    vector) this is (p·p)², so k/|p|² is a unit vector. -/
theorem ray_norm (p0 p1 p2 q0 q1 q2 : R) :
    (p0 + (q1 * p2 - q2 * p1))^2 + (p1 + (q2 * p0 - q0 * p2))^2 + (p2 + (q0 * p1 - q1 * p0))^2
      = (p0^2 + p1^2 + p2^2) * (1 + q0^2 + q1^2 + q2^2) - (p0 * q0 + p1 * q1 + p2 * q2)^2 := by grind

/-! ### Round 11: two bodies at rest (checks B3, B4) -/

/-- **Two bodies combine by velocity addition (B3).** Cross-multiplied form of
    (1 − u)/(1 + u) = ((1 − u₁)/(1 + u₁))·((1 − u₂)/(1 + u₂)) for u = (u₁ + u₂)/(1 + u₁u₂):
    on the axis each body's clock factor (1 − u)/(1 + u) multiplies, and the u's add like speeds. -/
theorem velocity_addition (u1 u2 : R) :
    ((1 + u1 * u2) - (u1 + u2)) * ((1 + u1) * (1 + u2))
      = ((1 + u1 * u2) + (u1 + u2)) * ((1 - u1) * (1 - u2)) := by grind

/-- **The strut's push (B4).** With e^γ = (ℓ² − (M₁+M₂)²)/(ℓ² − (M₁−M₂)²) on the strut,
    (e^−γ − 1)/4 = M₁M₂/(ℓ² − (M₁+M₂)²): the numerators differ by exactly 4 M₁M₂. -/
theorem strut_push (l M1 M2 : R) :
    (l^2 - (M1 - M2)^2) - (l^2 - (M1 + M2)^2) = 4 * (M1 * M2) := by grind


/-! ### Round 12: two moving bodies (check O1) -/

/-- **Clock rates multiply (O1, rule B).** Each body alone slows clocks by 1 − 2U + 2U² at slow-orbit order.
    Multiplying the two bodies' factors gives the one-body form with U = U₁ + U₂, cross term 4U₁U₂ included,
    plus only terms of higher order (U₁U₂ times another U). -/
theorem clock_multiply (a b : R) :
    (1 - 2 * a + 2 * a^2) * (1 - 2 * b + 2 * b^2)
      = 1 - 2 * (a + b) + 2 * (a + b)^2 - 4 * (a * b) * (a + b) + 4 * (a^2 * b^2) := by grind

/-- **Round 13: the 2PN miss, two ways of writing it** (times 6, to stay in a ring). With total
mass 1, body 1 of mass a and body 2 of mass 1 − a, the three rules' redshift coefficient of x³ is
off by a(1 − a)(23a² − 50a + 36)/3. In the usual variables ν = a(1 − a) and Δ = (1 − a) − a that is
ν(15/2 + 9/2·Δ − 23ν/3). -/
theorem two_pn_miss (a : R) :
    2 * a * (1 - a) * (23 * a^2 - 50 * a + 36)
      = a * (1 - a) * (45 + 27 * ((1 - a) - a) - 46 * (a * (1 - a))) := by grind

/-- **Round 14: adding gravity's own field energy (check O10)**, times 3 to stay in a ring. With body 1 of
mass a and body 2 of mass 1 − a, the closed-form field-energy pieces leave the miss a(1 − a)(5 − 10a + 3a²).
The four static field-on-field numbers α = −1/2, β = 1, γ = 0, δ = −2 change body 1's x³ coefficient by the
Newtonian-order map in the second bracket. What is left is a(1 − a)²(3 − a)/3. -/
theorem fourth_rule_residual (a : R) :
    3 * (a * (1 - a) * (5 - 10 * a + 3 * a^2))
      + (8 * (1 - a) * (a * (1 - a)) + 16 * (1 - a)^2 * (a * (1 - a)) - 24 * a * (1 - a)^3
         + 6 * a^2 * (1 - a) - 12 * a * (1 - a)^2)
      = a * (1 - a)^2 * (3 - a) := by grind

/-- **The leftover is one point-like term (check O10).** Raising β by 1/4 and lowering δ by 1/4 (an extra
X̂ ⊃ m_a m_b²/(4 r₁₂² r_b) at each body) cancels the residual exactly, for every mass ratio:
the map at Δβ = 1/4, Δδ = −1/4, times 3, is minus the residual above. -/
theorem fourth_rule_point_term (a : R) :
    ((1 - a) * (a * (1 - a)) - (1 - a)^2 * (a * (1 - a))
      - 3 * a * (1 - a)^2) = -(a * (1 - a)^2 * (3 - a)) := by grind

/-- **Round 15: the leftover derived (check O11)**, everything times 4 to stay in a ring. A body's potential
is M/r plus, from the delay, −(M/2) s with s = a·n its own acceleration along the direction to the field point;
write u = 1/r, so 2 × V = 2Mu − Ms. The cross term of g₀₀'s −2V² is then −4 times the acceleration piece already
in X̂, Xacc = −(M²/2) s u (4 Xacc = −2M² s u). So g₀₀'s self part is 8 Xacc − 4 Xacc, which is the
"quarter point term" of round 14. -/
theorem own_field_squared (M s u : R) :
    -2 * ((2 * M * u - M * s)^2 - (2 * M * u)^2 - (M * s)^2) = -4 * (-2 * M^2 * s * u) := by grind

end

end CircleOp.SpinningRiver
