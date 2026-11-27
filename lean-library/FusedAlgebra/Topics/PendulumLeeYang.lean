/-
  Fused Hyperoperation Algebra: the pendulum ideas checked against Lee–Yang zeros (DRAFT, Nov 27 2026).
  Core Lean 4 only (no mathlib, no `sorry`). Polynomial identities behind ly_checks.py, over Int.
  k is the dial (the dial gas's edge weight), m the pendulum's modulus, k_O Onsager's modulus.
  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.PendulumLeeYang
  (namespace CircleOp.PendulumLeeYang, the library's convention; no proof changed).
-/
namespace CircleOp.PendulumLeeYang

/-- Lattice gas to Ising, one edge: with σ = 2n − 1, 4·n_i n_j = 1 + σ_i + σ_j + σ_i σ_j. -/
theorem lattice_gas_edge (a b : Int) :
    4 * (a * b) = 1 + (2 * a - 1) + (2 * b - 1) + (2 * a - 1) * (2 * b - 1) := by grind

/-- Onsager's modulus k_O = 4k/(k − 1)² equals 1 exactly when k² − 6k + 1 = 0, i.e. k = 3 ± 2√2. -/
theorem onsager_separatrix (k : Int) : ((k - 1) * (k - 1) = 4 * k) ↔ (k * k - 6 * k + 1 = 0) := by
  constructor <;> intro h <;> grind

/-- The ring's edge as a pendulum amplitude: if the field-angle edge has k·sin²(φ/2) = 1 and
    sin² + cos² = 1, then the pendulum angle θ₀ = π − φ (whose sin(θ₀/2) is cos(φ/2)) has k·sin²(θ₀/2) = k − 1,
    i.e. m = sin²(θ₀/2) = 1 − 1/k. -/
theorem ring_edge_amplitude (k s c : Int) (h1 : s * s + c * c = 1) (h2 : k * (s * s) = 1) :
    k * (c * c) = k - 1 := by grind

/-- The ring's zeros, cleared: 2cos φ + 2 = c_j(1 − 1/k) with c_j = 4w², 2cos φ + 2 = 4v² (v = cos(φ/2))
    gives k·v² = (k − 1)·w², i.e. sin(θ/2) = √m · cos(ψ_j/2). -/
theorem ring_zero_pendulum (k v w : Int) (h : k * (4 * v * v) = 4 * w * w * (k - 1)) :
    k * (v * v) = (k - 1) * (w * w) := by grind

/-- Density: with y = sin(θ/2) = √m · sin τ, the zeros (uniform in τ) have
    dτ/dθ = cos(θ/2)/(2√(m − y²)): the identity m·cos²τ = m − y² behind it. -/
theorem density_identity (m sτ cτ : Int) (h : sτ * sτ + cτ * cτ = 1) :
    m * (cτ * cτ) = m - m * (sτ * sτ) := by grind

end CircleOp.PendulumLeeYang
