/-
  Does the theory predict an unknown force? (DRAFT, Nov 26 2026)
  Uncompiled here: no Lean toolchain in this container.
  Every statement is checked true by check_statements.py.
  Pure algebra over ℝ (ring, field_simp, norm_num).
-/
import Mathlib

namespace FusedAlgebra.Topics.NewForce

/-- The rung tower, every rung coupled equally, sums to coth:
    with q = e^{−r/R} (0 ≤ q < 1), 1 + 2q/(1−q) = (1+q)/(1−q), and (1+q)/(1−q) = coth(r/2R). -/
theorem tower_sum (q : ℝ) (hq : q ≠ 1) :
    1 + 2 * q / (1 - q) = (1 + q) / (1 - q) := by
  have h : (1 - q) ≠ 0 := sub_ne_zero.mpr (Ne.symm hq)
  field_simp
  ring

/-- Far away the tower is one Yukawa term of strength 2: (1+q)/(1−q) − 1 − 2q = 2q²/(1−q). -/
theorem tower_leading_yukawa (q : ℝ) (hq : q ≠ 1) :
    (1 + q) / (1 - q) - 1 - 2 * q = 2 * q ^ 2 / (1 - q) := by
  have h : (1 - q) ≠ 0 := sub_ne_zero.mpr (Ne.symm hq)
  field_simp
  ring

/-- Flux exponent s = 3/2 in D-language is field exponent q = s/(s−1) = 3. -/
theorem dual_exponent_mond : (3 / 2 : ℝ) / (3 / 2 - 1) = 3 := by norm_num

/-- Hölder pairing: 1/s + 1/q = 1 for s = 3/2, q = 3. -/
theorem holder_pair : (1 : ℝ) / (3 / 2) + 1 / 3 = 1 := by norm_num

/-- Asymptotic freedom count for SU(3): 11·3/3 − 2·16/3 > 0 and 11·3/3 − 2·17/3 < 0. -/
theorem su3_flavour_bound :
    (0 : ℝ) < 11 * 3 / 3 - 2 * 16 / 3 ∧ 11 * 3 / 3 - 2 * 17 / 3 < (0 : ℝ) := by
  constructor <;> norm_num

/-- Asymptotic freedom count for SU(2): 11·2/3 − 2·10/3 > 0, and at 11 flavours it is exactly 0. -/
theorem su2_flavour_bound :
    (0 : ℝ) < 11 * 2 / 3 - 2 * 10 / 3 ∧ 11 * 2 / 3 - 2 * 11 / 3 = (0 : ℝ) := by
  constructor <;> norm_num

/-! ### Part G: strength from a closed direction (electromagnetism as gravity's circle) -/

/-- Copies split into 4D pieces: tensor 1/2, vector 2, scalar 3/2 (per E₁E₂).
    The vector is 4 times the tensor and the scalar 3 times. -/
theorem vector_over_tensor : (2 : ℝ) / (1 / 2) = 4 := by norm_num
theorem scalar_over_tensor : (3 / 2 : ℝ) / (1 / 2) = 3 := by norm_num

/-- Like circle charges: tensor + scalar − vector = 0. Opposite: tensor + scalar + vector = 4. -/
theorem like_charges_cancel : (1 / 2 : ℝ) + 3 / 2 - 2 = 0 := by norm_num
theorem opposite_charges_add : (1 / 2 : ℝ) + 3 / 2 + 2 = 4 := by norm_num

/-- Klein's relation in our reading: with G = 1 and circle radius R = 2/√α,
    4/R² = α. Stated with R² = 4/α. -/
theorem klein_strength (α : ℝ) (hα : α ≠ 0) : 4 / (4 / α) = α := by
  field_simp

/-- n circles give a scalar force of strength a = n/(n+2); one circle gives 1/3. -/
theorem one_circle_third : (1 : ℝ) / (1 + 2) = 1 / 3 := by norm_num

/-- A scalar of strength a gives γ = (1 − a)/(1 + a); with a = n/(n+2) that is 1/(n+1). -/
theorem gamma_from_strength (n : ℝ) (hn : 0 ≤ n) :
    (1 - n / (n + 2)) / (1 + n / (n + 2)) = 1 / (n + 1) := by
  have h2 : n + 2 ≠ 0 := by positivity
  have h1 : n + 1 ≠ 0 := by positivity
  have h3 : 2 * n + 2 ≠ 0 := by positivity
  field_simp
  ring

/-- One circle: γ − 1 = −1/2. -/
theorem one_circle_gamma : (1 : ℝ) / (1 + 1) - 1 = -1 / 2 := by norm_num

/-! ### Part H: holding the circle (vacuum energy plus a 5D vacuum energy) -/

/-- With V = a/R + C R*²/R⁶ and a = −6C/R*³, R* is a turning point: V′(R*) = 0. -/
theorem hold_turning_point (C R : ℝ) (hR : R ≠ 0) :
    -(-6 * C / R ^ 3) / R ^ 2 - 6 * C * R ^ 2 / R ^ 7 = 0 := by
  field_simp
  ring

/-- Its curvature there: V″(R*) = 30C/R*⁶ (a minimum when C > 0). -/
theorem hold_curvature (C R : ℝ) (hR : R ≠ 0) :
    2 * (-6 * C / R ^ 3) / R ^ 3 + 42 * C * R ^ 2 / R ^ 8 = 30 * C / R ^ 6 := by
  field_simp
  ring

/-- The vacuum energy left over: V(R*) = −5C/R*⁴. -/
theorem hold_leftover (C R : ℝ) (hR : R ≠ 0) :
    (-6 * C / R ^ 3) / R + C * R ^ 2 / R ^ 6 = -5 * C / R ^ 4 := by
  field_simp
  ring

/-- Radion mass: (2/3)·30 = 20, so m² = 20C/(M̄²R⁴). -/
theorem hold_mass_factor : (2 / 3 : ℝ) * 30 = 20 := by norm_num

/-! ### Part S: strength from shape and material -/

/-- Sphere: ⟨|K|²⟩ = 2R²/3, so the shape factor 4R²/⟨|K|²⟩ is 6 (a circle gives 4). -/
theorem sphere_shape_factor : (4 : ℝ) / (2 / 3) = 6 := by norm_num

/-- Strong force with 5 flavours: b₀ = 11·3/3 − 2·5/3 = 23/3. -/
theorem strong_b0_five : (11 * 3 / 3 : ℝ) - 2 * 5 / 3 = 23 / 3 := by norm_num

/-! ### Part I: the new force up close -/

/-- Each rung between static masses: 5D weight 2/3 over 4D weight 1/2 is 4/3. -/
theorem every_rung_four_thirds : (2 / 3 : ℝ) / (1 / 2) = 4 / 3 := by norm_num

/-- The zero rung's 4/3 is gravity (1) plus the new force (1/3). -/
theorem zero_rung_split : (1 : ℝ) + 1 / 3 = 4 / 3 := by norm_num

/-- Copies sit on the extreme-black-hole line: charge²/(G mass²) = 1 + a² = 1 + 3 = 4. -/
theorem copies_extreme : (1 : ℝ) + 3 = 4 := by norm_num

/-- Near a static mass the circle swells: h₅₅/h₀₀ = (1/3)/(2/3) = 1/2. -/
theorem circle_swell_ratio : (1 / 3 : ℝ) / (2 / 3) = 1 / 2 := by norm_num

/-- Electric strength dips twice as much: −2 · (2/3) = −4/3. -/
theorem alpha_dip : (-2 : ℝ) * (2 / 3) = -4 / 3 := by norm_num

/-- The circle's potential tail e^{−√(2/3) σ} has slope parameter (1/2)(2/3) = 1/3. -/
theorem inflaton_tail_eps : (1 / 2 : ℝ) * (2 / 3) = 1 / 3 := by norm_num

/-! ### Part J: the inverse, and passing for gravity -/

/-- Beyond the size-carrier's range two like copies feel gravity (+1) and the vector (−4): a push of 3. -/
theorem like_copies_far_repel : (1 : ℝ) + 0 - 4 = -3 := by norm_num

/-- Matter (size charge² 1/3) and a copy (size charge² 3): the cross weight is (1/3)·3 = 1. -/
theorem matter_copy_size : (1 / 3 : ℝ) * 3 = 1 := by norm_num

/-- At the cap R₂ = R₁³/(4π²) the first circle equals the 6D fundamental length: 4π²R₁R₂ = R₁⁴. -/
theorem second_circle_cap (R₁ : ℝ) :
    4 * Real.pi ^ 2 * R₁ * (R₁ ^ 3 / (4 * Real.pi ^ 2)) = R₁ ^ 4 := by
  have h : (4 * Real.pi ^ 2 : ℝ) ≠ 0 := by positivity
  field_simp

/-- The cap as a ratio of the two circles: R₂/R₁ ≤ R₁²/(4π²). -/
theorem cap_ratio (R₁ : ℝ) (hR : R₁ ≠ 0) :
    (R₁ ^ 3 / (4 * Real.pi ^ 2)) / R₁ = R₁ ^ 2 / (4 * Real.pi ^ 2) := by
  have h : (4 * Real.pi ^ 2 : ℝ) ≠ 0 := by positivity
  field_simp

/-- A 1/3 force on matter but not light: light bends by (1 + γ)/2 = 3/4 of the expected amount. -/
theorem light_bends_three_quarters :
    (1 / 2 : ℝ) * (1 + (1 - 1 / 3) / (1 + 1 / 3)) = 3 / 4 := by norm_num

/-! ### Part K: the other unknown forces -/

/-- SU(3) with 16 flavours: b₀ = 11 − 32/3 = 1/3, still asymptotically free. -/
theorem bz_b0_sixteen : (11 * 3 / 3 : ℝ) - 2 * 16 / 3 = 1 / 3 := by norm_num

/-- SU(3) with 16 flavours: b₁ = 102 − 38·16/3 = −302/3 < 0, so a fixed point exists. -/
theorem bz_b1_sixteen : (102 : ℝ) - 38 * 16 / 3 = -302 / 3 := by norm_num

/-- The frozen strength: α* = −4π b₀/b₁ = 4π/302 at 16 flavours. -/
theorem bz_alpha_sixteen : -4 * Real.pi * (1 / 3) / (-302 / 3 : ℝ) = 4 * Real.pi / 302 := by
  field_simp

/-- A square torus does not mix the two photons: ε = cos(π/2) = 0. -/
theorem square_torus_no_mixing : Real.cos (Real.pi / 2) = 0 := Real.cos_pi_div_two

/-! ### Parts L and M: chasing the two possible leads -/

/-- Our floor 1/R1 = MP/√(4·137) and the magnetic weak-gravity ceiling g·M̄P = MP/√(2·137) differ by √2: squared, 2. -/
theorem mwgc_ratio_sq : (4 * 137 : ℝ) / (2 * 137) = 2 := by norm_num

/-- On a sphere the rotation arrow's squared length ρ² sin²θ averages 2ρ²/3: ∫ sin³θ dθ = 4/3, times 2π, over 4π. -/
theorem sphere_arrow_average : 2 * Real.pi * (4 / 3) = (2 / 3 : ℝ) * (4 * Real.pi) := by ring

/-- So a sphere's force has α = 4G/(2ρ²/3) = 6G/ρ². -/
theorem sphere_strength_six : (4 : ℝ) / (2 / 3) = 6 := by norm_num

/-- Circle edge: x = 4/R2², R2³ = 4π²R1, R1² = 4/α give x³·(4π²)²·4 = 64α, i.e. x³ = α/π⁴. -/
theorem circle_edge_cubed (a : ℝ) : (4 * Real.pi ^ 2) ^ 2 * 4 * a = 64 * (a * Real.pi ^ 4) := by ring

/-- Sphere edge: x = 6/ρ², ρ³ = 8π²R1 give x³ = 216α/(256π⁴), i.e. x = 6(α/(256π⁴))^(1/3). -/
theorem sphere_edge_cubed (a : ℝ) :
    (8 * Real.pi ^ 2) ^ 2 * 4 * (216 * a) = 216 * (256 * a * Real.pi ^ 4) := by ring

/-! ### Part N: does the hold put a sphere at the fundamental length? -/

/-- The sphere's hold U = −2/ρ⁴ + K/ρ⁹ turns at ρ⁵ = 9K/8: there 8/ρ⁵ = 9K/ρ¹⁰, i.e. 8·(9K/8) = 9K. -/
theorem sphere_hold_rest (K : ℝ) : 8 * (9 * K / 8) = 9 * K := by ring

/-- Its curvature there: ρ⁶U'' = −40 + 90·(8/9) = 40 > 0, a minimum. -/
theorem sphere_hold_minimum : (-40 : ℝ) + 90 * (8 / 9) = 40 := by norm_num

/-- The leftover 4D vacuum energy: ρ⁴U = −2 + 8/9 = −10/9. -/
theorem sphere_hold_leftover : (-2 : ℝ) + 8 / 9 = -10 / 9 := by norm_num

/-- Strength against the vacuum-energy size: α ∝ s^(−10/3) and s ∝ K^(1/5) give α ∝ K^(−2/3). -/
theorem sphere_strength_exponent : (-10 / 3 : ℝ) * (1 / 5) = -2 / 3 := by norm_num

end FusedAlgebra.Topics.NewForce
