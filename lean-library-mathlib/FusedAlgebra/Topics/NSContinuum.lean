module

public import FusedAlgebra.Cube

/-!
# The dimension continuum and the 3D Navier–Stokes gap

The exact, algebraic statements behind the "dimension continuum" picture.

1. The dyadic shell is the mult-add cube: doubling a box in `n` dimensions adds
   `2 ^ (n * j) * mpow 1 n` new cells (`annulus_cells`), and those cells fall into the cube's
   layers `1, …, n` (`annulus_graded`).
2. The fan: every scaling threshold is a ray from the pivot `(n, α) = (-2, 0)`; lying on the same
   ray is reflexive, symmetric and transitive (`onRay_refl`, `onRay_symm`, `onRay_trans`), so
   sliding the dimension and sliding the dissipation are one dial. 3D Navier–Stokes sits on the
   ray of 2D flow with dissipation power `4 / 5`.
3. The K41 edge of the dyadic model: the constant-flux margin `6 * α - (n + 2)` vanishes at
   `n = 4` when `α = 1`.
4. The dyadic energy identity (`energy_telescopes`): the nonlinearity only passes energy upward.
5. Planar flow has an empty stretching corner (`planar_stretching_zero`).

This is the mathlib form of `FusedAlgebra.Topics.NSContinuum`. The corner count
`mpow 1 n + 1 = 2 ^ n` is `MultAdd.Cube.mpow_one_add_one`, the hand-made sums are `Finset.sum`s,
and vortex stretching is `Matrix.mulVec`. Nothing here is about actual Navier–Stokes solutions.
-/

@[expose] public section

namespace MultAdd.NSContinuum

open MultAdd MultAdd.Cube Finset Matrix

/-! ### The dyadic shell is the mult-add cube -/

/-- Doubling the side of an `n`-dimensional box of side `2 ^ j` multiplies its cells by `2 ^ n`. -/
theorem box_doubles (n j : ℕ) : 2 ^ (n * (j + 1)) = 2 ^ (n * j) * 2 ^ n := by
  rw [mul_add_one, pow_add]

/-- The dyadic annulus is a mult-add power: the new cells added when a box of side `2 ^ j`
doubles number `2 ^ (n * j) * mpow 1 n`, the old box times the `n`-fold mult-add of `1`. -/
theorem annulus_cells (n j : ℕ) : 2 ^ (n * (j + 1)) - 2 ^ (n * j) = 2 ^ (n * j) * mpow 1 n := by
  rw [box_doubles, ← mpow_one_add_one n, mul_add_one, Nat.add_sub_cancel]

/-- The annulus is graded by the cube's layers: `2 ^ n - 1 = ∑ k = 1..n, n.choose k`. The new
cells of a doubled box split by how many of the `n` coordinates landed in the upper half. -/
theorem annulus_graded (n : ℕ) : (2 : ℤ) ^ n - 1 = ∑ k ∈ range n, (n.choose (k + 1) : ℤ) := by
  have h := mpow_eq_sum_choose (1 : ℤ) n
  rw [mpow_eq, one_add_one_eq_two] at h
  simpa using h

/-! ### The fan of scaling rays

Dissipation is measured in twentieths, `A = 20 * α`, so 3D Navier–Stokes is `(n, A) = (3, 20)`,
Lions' threshold in 3D is `(3, 25)` and 2D Navier–Stokes is `(2, 20)`. Two points lie on the same
scaling ray when `A * (n' + 2) = A' * (n + 2)`: the ratio `2 * α / (n + 2)` is the same. -/

section Fan

variable {R : Type*} [CommRing R]

/-- Two points `(n, A)` and `(n', A')` lie on the same scaling ray through the pivot
`(n, A) = (-2, 0)`. -/
def onRay (n A n' A' : R) : Prop := A * (n' + 2) = A' * (n + 2)

instance (n A n' A' : ℤ) : Decidable (onRay n A n' A') :=
  inferInstanceAs (Decidable (_ = _))

theorem onRay_refl (n A : R) : onRay n A n A := rfl

theorem onRay_symm {n A n' A' : R} (h : onRay n A n' A') : onRay n' A' n A :=
  Eq.symm h

/-- Lying on one ray is transitive (away from the pivot dimension `n = -2`), so the rays split
the plane: sliding the dimension and sliding the dissipation are one dial. -/
theorem onRay_trans [IsDomain R] {n A n' A' n'' A'' : R} (hne : n' + 2 ≠ 0)
    (h1 : onRay n A n' A') (h2 : onRay n' A' n'' A'') : onRay n A n'' A'' := by
  unfold onRay at *
  apply mul_right_cancel₀ hne
  linear_combination (n'' + 2) * h1 + (n + 2) * h2

end Fan

/-- Lions' 3D threshold `α = 5 / 4` is on the same ray as 2D Navier–Stokes (`α = 1`): both are
critical. -/
theorem lions_on_2d_ray : onRay (3 : ℤ) 25 2 20 := by decide

/-- 3D Navier–Stokes (`α = 1`) is on the same ray as 2D flow with dissipation power
`α = 4 / 5`. -/
theorem ns3d_on_ray_of_2d_four_fifths : onRay (3 : ℤ) 20 2 16 := by decide

/-- 3D Navier–Stokes is off the critical ray: the scaling gap is a quarter power. -/
theorem ns3d_supercritical : ¬ onRay (3 : ℤ) 20 3 25 := by decide

/-- The three known thresholds are rays: critical `α = (n + 2) / 4`, the dyadic regularity line
`α = (n + 2) / 5` (Barbato–Morandin–Romito) and the dyadic blow-up line `α = (n + 2) / 6`
(Cheskidov). In twentieths: `5 * (n + 2)`, `4 * (n + 2)`, `(10 / 3) * (n + 2)`; all three vanish
at the pivot `n = -2`. -/
theorem fan_pivot :
    5 * ((-2 : ℤ) + 2) = 0 ∧ 4 * ((-2 : ℤ) + 2) = 0 ∧ 10 * ((-2 : ℤ) + 2) = 0 := by
  decide

/-- 3D Navier–Stokes lies exactly on the dyadic regularity line `A = 4 * (n + 2)`. -/
theorem ns3d_on_bmr_line : (20 : ℤ) = 4 * (3 + 2) := by decide

/-! ### The K41 edge of the dyadic model

Constant energy flux `λ ^ j * a j ^ 3 = const` gives `a j ∝ λ ^ (-j / 3)`. Dissipation beats
transfer at small scales iff `6 * α > n + 2`. -/

/-- The constant-flux margin `6 * α - (n + 2)`, at `α = 1`. -/
def k41Margin (n : ℤ) : ℤ := 6 - (n + 2)

theorem k41_edge_at_four : k41Margin 4 = 0 := by decide

theorem k41_3d_positive : 0 < k41Margin 3 := by decide

theorem k41_5d_negative : k41Margin 5 < 0 := by decide

/-! ### The dyadic energy identity -/

section Energy

variable {R : Type*} [CommRing R]

/-- Energy flux from shell `j` to shell `j + 1` in the dyadic model
`a j' = c (j - 1) * a (j - 1) ^ 2 - c j * a j * a (j + 1)` (times `a j`, for the energy). -/
def flux (c a : ℕ → R) (j : ℕ) : R := c j * a j * a j * a (j + 1)

/-- Shell `j`'s energy gain from the nonlinearity. -/
def gain (c a : ℕ → R) : ℕ → R
  | 0 => -flux c a 0
  | j + 1 => flux c a j - flux c a (j + 1)

/-- The gain is `a j` times the model's nonlinear term. -/
theorem gain_is_model (c a : ℕ → R) (j : ℕ) :
    gain c a (j + 1) = a (j + 1) * (c j * a j ^ 2 - c (j + 1) * a (j + 1) * a (j + 2)) := by
  simp only [gain, flux]; ring

/-- Energy only moves: summed over shells `0, …, N`, the nonlinear gains telescope to minus the
flux out of the top shell. -/
theorem energy_telescopes (c a : ℕ → R) (N : ℕ) :
    ∑ j ∈ range (N + 1), gain c a j = -flux c a N := by
  induction N with
  | zero => simp [gain]
  | succ N ih => rw [sum_range_succ, ih, gain]; ring

end Energy

/-! ### The empty stretching corner of planar flow

Vortex stretching is `(ω · ∇) u`, with `i`-th component `∑ j, ω j * ∂ j (u i)`: the velocity
gradient matrix `G i j = ∂ j (u i)` applied to the vorticity, `G *ᵥ ω`. For planar flow
`u = (u₁ (x, y), u₂ (x, y), 0)` has `ω = (0, 0, w)` and nothing depends on `z`, so the only term
left, `w * ∂ z (u i)`, is zero: the corner that feeds 3D growth is empty. -/

/-- Planar flow has no vortex stretching: if the velocity gradient `G` has zero `z`-column and
the vorticity points along `z`, then `G *ᵥ ω = 0`. -/
theorem planar_stretching_zero {R : Type*} [CommRing R] (G : Matrix (Fin 3) (Fin 3) R) (w : R)
    (hz : ∀ i, G i 2 = 0) : G *ᵥ Pi.single 2 w = 0 := by
  ext i
  simp [hz]

end MultAdd.NSContinuum
-- touched
