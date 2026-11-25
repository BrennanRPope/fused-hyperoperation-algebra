/-
  Chemistry 3, the periodic law: shells are squares, periods are diagonals
  =======================================================================

  Companion page and checks: ../index.html, ../checks.py (12 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Hydrogen's shell n has subshells l = 0 … n − 1 with 2l + 1 orbitals each.
  They are the L-shaped gnomons of an n × n square, so the shell holds n²
  orbitals: n² = 1 + (n − 1) ✳ (n − 1). The square is SO(4) = SU(2) × SU(2):
  two equal spins j with 2j = n − 1, multiplicities multiplying.

  Madelung fills subshells along the diagonals n + l = s. Diagonal s holds
  2·⌈s/2⌉² electrons: 2, 2, 8, 8, 18, 18, 32, 32, which are Janet's
  left-step periods; the usual table merges the first two, giving
  2, 8, 8, 18, 18, 32, 32 and noble gases at 2, 10, 18, 36, 54, 86, 118.

  The dial: the q-integer [n]_q is n copies of 1 under weighted mult-add
  with k = q − 1, so a q-deformed shell would hold [n]_q² orbitals.

  LIBRARY NOTE: `local infixl " ✳ "` for multAdd; delete it if Base declares ✳.
-/
import FusedAlgebra.Base

namespace ChemPeriodic
open MultAdd

local infixl:65 " ✳ " => multAdd

/-- The shell square: (m + 1)² = 1 + m ✳ m, with m = n − 1. -/
theorem shell_square (m : Nat) : (m + 1) * (m + 1) = 1 + (m ✳ m) := by
  unfold multAdd; grind

/-- Orbitals in subshells l = 0 … m. -/
def shellOrbitals : Nat → Nat
  | 0 => 1
  | m + 1 => shellOrbitals m + (2 * (m + 1) + 1)

/-- The gnomons fill the square: shell n = m + 1 holds n² orbitals. -/
theorem gnomons_fill (m : Nat) : shellOrbitals m = (m + 1) * (m + 1) := by
  induction m with
  | zero => rfl
  | succ m ih => simp only [shellOrbitals, ih]; grind

/-- Each gnomon is corner + two edges of the next square: 2l + 1 = 1 + l + l. -/
theorem gnomon (l : Nat) : (l + 1) * (l + 1) - l * l = 1 + l + l := by
  have : (l + 1) * (l + 1) = l * l + (1 + l + l) := by grind
  omega

/-- With spin, shell n holds 2n² electrons. -/
theorem shell_capacity (m : Nat) : 2 * shellOrbitals m = 2 * ((m + 1) * (m + 1)) := by
  rw [gnomons_fill]

/-- SO(4): two independent spins with multiplicity n each give n · n states,
    and (n − 1) mult-adds with itself. -/
theorem two_spins (a b : Nat) : (a + 1) * (b + 1) = 1 + (a ✳ b) := by
  unfold multAdd; grind

/-! ## Madelung diagonals and the periods -/

/-- Electrons on the Madelung diagonal n + l = s: subshells l = 0 … ⌊(s−1)/2⌋,
    each holding 2(2l + 1). -/
def diagCap (s : Nat) : Nat := 2 * shellOrbitals ((s - 1) / 2)

theorem diagonals : (List.range 8).map (fun i => diagCap (i + 1)) = [2, 2, 8, 8, 18, 18, 32, 32] := by
  decide

/-- Diagonal s holds 2 · ⌈s/2⌉² electrons, in closed form. -/
theorem diagCap_closed (s : Nat) : diagCap s = 2 * (((s - 1) / 2 + 1) * ((s - 1) / 2 + 1)) := by
  unfold diagCap; rw [gnomons_fill]

/-- The usual table merges the first two diagonals into period 1. -/
def periodLen (p : Nat) : Nat := 2 * ((p / 2 + 1) * (p / 2 + 1))

theorem periods : (List.range 7).map (fun i => periodLen (i + 1)) = [2, 8, 8, 18, 18, 32, 32] := by
  decide

/-- Noble gases: running totals of the period lengths. -/
def nobleZ : Nat → Nat
  | 0 => 0
  | p + 1 => nobleZ p + periodLen (p + 1)

theorem noble_gases : (List.range 7).map (fun i => nobleZ (i + 1)) = [2, 10, 18, 36, 54, 86, 118] := by
  decide

/-- Period p (≥ 2) and the diagonal pair it covers agree: periods 2k and 2k+1
    have the length of diagonals 2k+1 and 2k+2. -/
theorem period_is_diagonal (k : Nat) : periodLen (2 * k + 1) = diagCap (2 * k + 2) := by
  unfold periodLen diagCap; rw [gnomons_fill]
  have h1 : (2 * k + 1) / 2 = k := by omega
  have h2 : (2 * k + 2 - 1) / 2 = k := by omega
  rw [h1, h2]

/-! ## The dial: q-integers -/

section Dial
variable {R : Type} [Lean.Grind.CommRing R]

/-- n copies of 1 under ⊞ₖ lift to (1 + k)ⁿ: with q = k + 1, that is qⁿ,
    so the weighted power is the q-integer (qⁿ − 1)/(q − 1). -/
theorem q_integer (k : R) (n : Nat) : 1 + k * wpow k 1 n = (1 + k) ^ n := by
  have := wpow_lift k 1 n
  rw [show k * (1 : R) = k by grind] at this
  exact this

end Dial

/-- At k = 0 (q = 1) the weighted power of 1 is n itself: the ordinary shell. -/
theorem q_integer_classical (n : Nat) : wpow (0 : Rat) 1 n = (n : Rat) := by
  induction n with
  | zero => simp only [wpow]; grind
  | succ n ih => simp only [wpow, wAdd, ih]; grind

end ChemPeriodic
