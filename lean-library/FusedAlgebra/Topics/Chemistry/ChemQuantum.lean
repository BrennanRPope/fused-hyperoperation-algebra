/-
  Chemistry 1, quantum chemistry: independent fragments compose by mult-add
  ==========================================================================

  Companion page and checks: ../index.html, ../checks.py (25 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`. Imports the shared base.

  Write a molecule's wavefunction in intermediate normalisation, Ψ = Φ(1 + c),
  where c is everything beyond the reference determinant. Two molecules that
  do not interact have Ψ_AB = Ψ_A Ψ_B, so

      1 + c_AB = (1 + c_A)(1 + c_B),   c_AB = c_A ◯ c_B.

  The corner c_A·c_B is "both molecules excited at once". Truncated CI (CISD)
  keeps only c_A + c_B: it is the k = 0 law, and its size-consistency error is
  exactly the dropped corner. Coupled cluster writes 1 + c = exp T, the log of
  the lift, so T adds and the corner is kept automatically.

  The N-molecule model of ../checks.py: each molecule has a reference and a
  doubly excited state 2Δ above it, coupled by g. CISD's energy E solves
  E(E − 2Δ) = N g² (its 2-level secular equation); the exact energy is N e₁.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged at Topics/Chemistry/ChemQuantum; it now uses Base's ◯.
  Four exact twins of Base removed (see the comment in the Ring section).
-/
import FusedAlgebra.Base

namespace ChemQuantum
open CircleOp

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- Two independent fragments: the excess amplitude mult-adds. -/
theorem fragments_circ (cA cB : R) : (1 + cA) * (1 + cB) = 1 + (cA ◯ cB) := by
  unfold circ; grind

/-- What truncated CI drops when two fragments are joined: exactly the corner. -/
theorem cisd_drops_corner (cA cB : R) : (cA ◯ cB) - (cA + cB) = cA * cB := by
  unfold circ; grind

-- LIBRARY: four exact twins of Base theorems removed; the chemistry reading stays here.
--   cisd_is_weight_zero (k = 0 is CISD)        = CircleOp.circK_weight_zero
--   exact_is_weight_one (k = 1 is exact)       = CircleOp.circK_one
--   n_fragments (1 + c_N = (1 + c)^N)          = CircleOp.mpow_lift
--   fragments_assoc (join N then M = N + M)    = CircleOp.mpow_add

/-- Three fragments: the corner expansion has pairs and one triple, nothing else. -/
theorem three_fragments (a b c : R) :
    (a ◯ b) ◯ c = a + b + c + (a * b + a * c + b * c) + a * b * c := by
  unfold circ; grind

/-- What CISD drops from three fragments: every pair and the triple. -/
theorem cisd_drops_three (a b c : R) :
    ((a ◯ b) ◯ c) - (a + b + c) = a * b + a * c + b * c + a * b * c := by
  unfold circ; grind

end Ring

/-! ## The N-molecule CISD model (over Rat) -/

/-- The CISD secular equation of N molecules: E(E − 2Δ) = N g². -/
def secular (D g E N : Rat) : Prop := E * (E - 2 * D) = N * g * g

/-- If e₁ solves the one-molecule problem, N·e₁ misses the N-molecule CISD
    equation by exactly N(N − 1)e₁², the pairs of molecules both excited. -/
theorem cisd_gap (D g e N : Rat) (h1 : secular D g e 1) :
    (N * e) * (N * e - 2 * D) - N * g * g = N * (N - 1) * (e * e) := by
  unfold secular at h1; grind

/-- So N·e₁ solves N-molecule CISD only when N = 0, N = 1 or e₁ = 0:
    CISD is not size-extensive. -/
theorem cisd_not_extensive (D g e N : Rat) (hN0 : N ≠ 0) (hN1 : N ≠ 1) (he : e ≠ 0)
    (h1 : secular D g e 1) : ¬ secular D g (N * e) N := by
  intro hN'
  unfold secular at hN'
  have hgap := cisd_gap D g e N h1
  have hz : N * (N - 1) * (e * e) = 0 := by grind
  rcases Rat.mul_eq_zero.mp hz with h | h
  · rcases Rat.mul_eq_zero.mp h with h' | h'
    · exact hN0 h'
    · exact hN1 (by grind)
  · rcases Rat.mul_eq_zero.mp h with h' | h' <;> exact he h'

/-! ## Spin coupling counts states by mult-add -/

/-- Multiplicities multiply, so 2s (multiplicity − 1) mult-adds. -/
theorem spin_multiplicity (s1 s2 : Nat) : (s1 + 1) * (s2 + 1) = (s1 ◯ s2) + 1 := by
  unfold circ; grind

/-- The Clebsch–Gordan series for two equal spins j ⊗ j (2j = n): the
    multiplicities 2S+1 for S = 0..n sum to (n+1)². -/
def oddSum : Nat → Nat
  | 0 => 1
  | n + 1 => oddSum n + (2 * (n + 1) + 1)

theorem clebsch_gordan_square (n : Nat) : oddSum n = (n + 1) * (n + 1) := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [oddSum, ih]; grind

end ChemQuantum
