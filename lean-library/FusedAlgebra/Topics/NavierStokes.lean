/-
  Navier–Stokes, the exact statements behind the page (Burgers' corner,
  the helical shell model, the triad-sign cube).

  From ns-rungs/lean/FusedRungs.lean; the rung ladder itself moved to
  FusedAlgebra.Rungs. Core Lean only, no `sorry`.
-/
import FusedAlgebra.Rungs

namespace FusedRungs

open MultAdd

/-! ### Burgers: why the top corner of the Cole–Hopf cube must be empty
  A term e^(a x + b t) solves the heat equation φ_t = ν φ_xx iff b = ν a².
  The corner term of e1·e2 has slope a1+a2 and rate b1+b2; its heat defect is -2 ν a1 a2. -/
theorem burgers_corner_defect (ν a1 a2 : Int) :
    (ν * a1 ^ 2 + ν * a2 ^ 2) - ν * (a1 + a2) ^ 2 = -2 * ν * a1 * a2 := by grind
theorem burgers_corner_ok_iff (ν a1 a2 : Int) (hν : ν ≠ 0) :
    (ν * a1 ^ 2 + ν * a2 ^ 2) - ν * (a1 + a2) ^ 2 = 0 ↔ a1 = 0 ∨ a2 = 0 := by
  rw [burgers_corner_defect]
  constructor
  · intro h
    have h2 : ν * (a1 * a2) = 0 := by grind
    rcases Int.mul_eq_zero.mp h2 with h3 | h3
    · exact absurd h3 hν
    · exact Int.mul_eq_zero.mp h3
  · rintro (h | h) <;> subst h <;> simp

/-! ### Helical shell model: every corner of the triad-sign cube conserves energy and helicity
  c0 = (s1 l - s2 l²) G,  c1 = (s2 l² - s0) G,  c2 = c0 + c1,  G = s0 + s1 l + s2 l².
  Energy change ∝ c0 + c1 - c2; helicity change ∝ s0 c0 + s1 l c1 - s2 l² c2.  -/
def c0 (s0 s1 s2 l : Int) : Int := (s1 * l - s2 * l ^ 2) * (s0 + s1 * l + s2 * l ^ 2)
def c1 (s0 s1 s2 l : Int) : Int := (s2 * l ^ 2 - s0) * (s0 + s1 * l + s2 * l ^ 2)
def c2 (s0 s1 s2 l : Int) : Int := c0 s0 s1 s2 l + c1 s0 s1 s2 l
theorem shell_energy (s0 s1 s2 l : Int) : c0 s0 s1 s2 l + c1 s0 s1 s2 l - c2 s0 s1 s2 l = 0 := by
  unfold c2; omega
/-- Helicity is conserved for every sign corner and every shell ratio l, as a polynomial identity
    (it does not even need s² = 1). -/
theorem shell_helicity (s0 s1 s2 l : Int) :
    s0 * c0 s0 s1 s2 l + s1 * l * c1 s0 s1 s2 l - s2 * l ^ 2 * c2 s0 s1 s2 l = 0 := by
  unfold c2 c0 c1; grind

/-! ### The triad-sign cube has layers 1, 3, 3, 1 (by number of minus signs) -/
def signs : List (Bool × Bool × Bool) :=
  [false, true].flatMap fun a => [false, true].flatMap fun b => [false, true].map fun c => (a, b, c)
def minuses (t : Bool × Bool × Bool) : Nat :=
  (if t.1 then 1 else 0) + (if t.2.1 then 1 else 0) + (if t.2.2 then 1 else 0)
theorem sign_cube_layers :
    (List.range 4).map (fun j => (signs.filter (fun t => minuses t == j)).length) = [1, 3, 3, 1] := by
  decide
/-- The homochiral corners (+++) and (---) are exactly layers 0 and 3. -/
theorem homochiral_layers (t : Bool × Bool × Bool) (ht : t.1 = t.2.1 ∧ t.2.1 = t.2.2) :
    minuses t = 0 ∨ minuses t = 3 := by
  obtain ⟨a, b, c⟩ := t; simp at ht; obtain ⟨rfl, rfl⟩ := ht; cases a <;> simp [minuses]


end FusedRungs
