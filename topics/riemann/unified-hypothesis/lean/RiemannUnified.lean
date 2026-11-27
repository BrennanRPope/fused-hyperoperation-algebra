import FusedAlgebra.Base

/-!
# Riemann: unified hypothesis and the real-rooted parent hunt (DRAFT)

Small exact facts behind the parent hunt, the two-spin picture and the bath.
-/

namespace RiemannUnified

/-- Rolle bookkeeping: n real zeros of a parent leave n - 1 gaps, each holding a zero of the child. -/
theorem rolle_gaps (n : Nat) (h : 1 ≤ n) : (n - 1) + 1 = n := by omega

/-- Antiderivative parent: if Jensen's formula asks for at least 2N zeros in a disc
but at most N + 1 of them are real, then at least N - 1 are off the line. -/
theorem antiderivative_offline (total real N : Nat) (ht : 2 * N ≤ total) (hr : real ≤ N + 1)
    (hle : real ≤ total) : N - 1 ≤ total - real := by omega

/-- The control sequence 1/(k+1) is strictly log-convex: k(k+2) < (k+1)², so every Turán test fails. -/
theorem turan_control (k : Nat) : k * (k + 2) < (k + 1) * (k + 1) := by grind

/-- Two spins near a zero of G: Z ≈ G² − J·G′². At the zero, Z = −J·G′², negative for J > 0
(two real zeros split off) and positive for J < 0 (none). -/
theorem two_spin_ferro (J g : Int) (hJ : 0 < J) (hg : 0 < g * g) : -J * (g * g) < 0 := by
  have h := Int.mul_pos hJ hg
  rw [Int.neg_mul]
  omega

theorem two_spin_antiferro (J g : Int) (hJ : J < 0) (hg : 0 < g * g) : 0 < -J * (g * g) :=
  Int.mul_pos (by omega) hg

/-- Two-level bath: k′ = −(1 + E) with E = e^{β(Δ−μ−ΩJ)}; half filling (E = 1) is exactly k′ = −2. -/
theorem half_filling (E : Int) : -(1 + E) = -2 ↔ E = 1 := by omega

/-- Three-level symmetric bath with weights (w, 1−2w, w), w = m/N with 0 < 2m < N:
the factor 1 − 2w + 2w·cosh is a Pólya universal factor iff (N − 2m) ≤ 2m, i.e. w ≥ 1/4. -/
theorem three_level_threshold (m N : Int) (h0 : 0 < m) (h1 : 2 * m < N) :
    (N - 2 * m ≤ 2 * m) ↔ N ≤ 4 * m := by omega

end RiemannUnified
