module

public import FusedAlgebra.Rungs

/-!
# Navier–Stokes: Burgers' corner, the helical shell model, the triad-sign cube

The exact statements behind the Navier–Stokes page.

## Main results

* `burgers_corner_defect`, `burgers_corner_ok_iff`: a term `exp (a * x + b * t)` solves the heat
  equation `φ_t = ν * φ_xx` iff `b = ν * a ^ 2`. The corner term of a product of two such terms
  has slope `a₁ + a₂` and rate `b₁ + b₂`; its heat defect is `-2 * ν * a₁ * a₂`, which vanishes
  only when one slope is zero. This is why the top corner of the Cole–Hopf cube must be empty.
* `shell_energy`, `shell_helicity`: in the helical shell model every corner of the triad-sign
  cube conserves energy and helicity, as polynomial identities.
* `sign_cube_layers`, `homochiral_layers`: the triad-sign cube has layers `1, 3, 3, 1` by number
  of minus signs, and the homochiral corners are layers `0` and `3`.

The identities were stated over `ℤ` in the core-Lean original; here they hold in any commutative
ring (a domain where `2 ≠ 0` for `burgers_corner_ok_iff`). The rung ladder itself is in
`FusedAlgebra.Rungs`.
-/

@[expose] public section

namespace FusedRungs

open Finset

variable {R : Type*} [CommRing R]

/-! ### Burgers: why the top corner of the Cole–Hopf cube must be empty -/

/-- The heat defect of the corner term is `-2 * ν * a₁ * a₂`. -/
theorem burgers_corner_defect (ν a₁ a₂ : R) :
    (ν * a₁ ^ 2 + ν * a₂ ^ 2) - ν * (a₁ + a₂) ^ 2 = -2 * ν * a₁ * a₂ := by
  ring

/-- The corner term solves the heat equation only when one of the two slopes is zero. -/
theorem burgers_corner_ok_iff [IsDomain R] [NeZero (2 : R)] (ν a₁ a₂ : R) (hν : ν ≠ 0) :
    (ν * a₁ ^ 2 + ν * a₂ ^ 2) - ν * (a₁ + a₂) ^ 2 = 0 ↔ a₁ = 0 ∨ a₂ = 0 := by
  simp [burgers_corner_defect, hν, (two_ne_zero : (2 : R) ≠ 0)]

/-! ### Helical shell model: every corner of the triad-sign cube conserves energy and helicity

With `G = s₀ + s₁ * l + s₂ * l ^ 2`: `c0 = (s₁ * l - s₂ * l ^ 2) * G`,
`c1 = (s₂ * l ^ 2 - s₀) * G`, `c2 = c0 + c1`. The energy change is proportional to
`c0 + c1 - c2` and the helicity change to `s₀ * c0 + s₁ * l * c1 - s₂ * l ^ 2 * c2`. -/

/-- First coupling coefficient of the helical shell model. -/
def c0 (s₀ s₁ s₂ l : R) : R := (s₁ * l - s₂ * l ^ 2) * (s₀ + s₁ * l + s₂ * l ^ 2)

/-- Second coupling coefficient of the helical shell model. -/
def c1 (s₀ s₁ s₂ l : R) : R := (s₂ * l ^ 2 - s₀) * (s₀ + s₁ * l + s₂ * l ^ 2)

/-- Third coupling coefficient of the helical shell model, fixed by energy conservation. -/
def c2 (s₀ s₁ s₂ l : R) : R := c0 s₀ s₁ s₂ l + c1 s₀ s₁ s₂ l

/-- Energy is conserved. -/
theorem shell_energy (s₀ s₁ s₂ l : R) : c0 s₀ s₁ s₂ l + c1 s₀ s₁ s₂ l - c2 s₀ s₁ s₂ l = 0 :=
  sub_self _

/-- Helicity is conserved for every sign corner and every shell ratio `l`, as a polynomial
identity (it does not even need `s ^ 2 = 1`). -/
theorem shell_helicity (s₀ s₁ s₂ l : R) :
    s₀ * c0 s₀ s₁ s₂ l + s₁ * l * c1 s₀ s₁ s₂ l - s₂ * l ^ 2 * c2 s₀ s₁ s₂ l = 0 := by
  unfold c2 c0 c1; ring

/-! ### The triad-sign cube has layers `1, 3, 3, 1` (by number of minus signs) -/

/-- The number of minus signs of a sign triple (`true` is minus). -/
def minuses (t : Bool × Bool × Bool) : ℕ :=
  (if t.1 then 1 else 0) + (if t.2.1 then 1 else 0) + (if t.2.2 then 1 else 0)

/-- The triad-sign cube has layers `1, 3, 3, 1`. -/
theorem sign_cube_layers :
    (List.range 4).map (fun j ↦ #{t : Bool × Bool × Bool | minuses t = j}) = [1, 3, 3, 1] := by
  decide

/-- Layer `j` of the triad-sign cube has `Nat.choose 3 j` corners. -/
theorem card_minuses_eq (j : ℕ) : #{t : Bool × Bool × Bool | minuses t = j} = Nat.choose 3 j := by
  obtain hj | hj := lt_or_ge j 4
  · interval_cases j <;> decide
  · rw [Nat.choose_eq_zero_of_lt (by lia), card_eq_zero, filter_eq_empty_iff]
    rintro ⟨a, b, c⟩ -
    cases a <;> cases b <;> cases c <;> simp [minuses] <;> lia

/-- The homochiral corners `(+++)` and `(---)` are exactly layers `0` and `3`. -/
theorem homochiral_layers (t : Bool × Bool × Bool) (ht : t.1 = t.2.1 ∧ t.2.1 = t.2.2) :
    minuses t = 0 ∨ minuses t = 3 := by
  obtain ⟨a, b, c⟩ := t
  obtain ⟨rfl, rfl⟩ : a = b ∧ b = c := ht
  cases a <;> simp [minuses]

end FusedRungs
