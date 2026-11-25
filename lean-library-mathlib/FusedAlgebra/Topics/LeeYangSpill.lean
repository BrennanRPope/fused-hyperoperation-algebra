module

public import FusedAlgebra.Basic

/-!
# Fill, linger, spill at the landscape's border

Repeat one landscape law with a fixed step `x`:
`a ↦ (a + x + k * a * x) / (1 - q * a * x)`. The step is a Möbius map with trace `2 + k * x` and
determinant `1 + k * x + q * x ^ 2`. Its discriminant is `x ^ 2 * (k ^ 2 - 4 * q)`
(`step_disc`), so for every step size it is

* hyperbolic (two walls, the orbit settles) when `4 * q < k ^ 2`,
* parabolic (one double wall) on the border `4 * q = k ^ 2`,
* a turn (the orbit lingers, spills through infinity, comes round) when `4 * q > k ^ 2`.

How far `a` moves in one step is `x * (1 + k * a + q * a ^ 2) / (1 - q * a * x)` (`step_gap`), and
`4 * q * (1 + k * a + q * a ^ 2) = (2 * q * a + k) ^ 2 + (4 * q - k ^ 2)` (`wall_square`): on the
turning side the slowest step is set by the distance to the border (`slowest_step`,
`turning_never_stalls`).

All of this is the standard classification of Möbius maps by trace squared over determinant and
completing the square. `-(k ^ 2 - 4 * q)`, up to the factor `x ^ 2`, is mathlib's `discrim q k 1`
of the rate `1 + k * a + q * a ^ 2`. The core-Lean original proved the two inequalities over `ℤ`;
here they hold in any linearly ordered commutative ring.
-/

@[expose] public section

namespace LeeYangSpill

open Matrix

section Ring

variable {R : Type*} [CommRing R]

/-- The step's discriminant: `trace ^ 2 - 4 * det = x ^ 2 * (k ^ 2 - 4 * q)`, the same sign for
every step size. -/
theorem step_disc (k q x : R) :
    (2 + k * x) ^ 2 - 4 * (1 + k * x + q * x ^ 2) = x ^ 2 * (k ^ 2 - 4 * q) := by
  ring

/-- One step moves `a` by `x * (1 + k * a + q * a ^ 2)`, over the denominator `1 - q * a * x`. -/
theorem step_gap (k q x a : R) :
    (a + x + k * a * x) - a * (1 - q * a * x) = x * (1 + k * a + q * a ^ 2) := by
  ring

/-- Completing the square: the step's size is a square plus the distance to the border. -/
theorem wall_square (k q a : R) :
    4 * q * (1 + k * a + q * a ^ 2) = (2 * q * a + k) ^ 2 + (4 * q - k ^ 2) := by
  ring

/-- Two walls: when `k ^ 2 - 4 * q = m ^ 2`, the step size factors through the two walls, where
it is `0`. -/
theorem walls_factor (k q a m : R) (h : k ^ 2 - 4 * q = m ^ 2) :
    4 * q * (1 + k * a + q * a ^ 2) = (2 * q * a + k - m) * (2 * q * a + k + m) := by
  linear_combination -h

/-- The border (`k = 2 * u`, `q = u ^ 2`): `n` steps of `x` from `0` land at
`n * x / (1 - (n - 1) * u * x)`. Stated for a general `t` in place of `n`, cross-multiplied: the
step takes the `t`-th point to the `(t + 1)`-th. The denominator is linear in `t`: it passes
infinity once, then fills toward the double wall `-1 / u` like `1 / t`, forever. -/
theorem border_orbit (u x t : R) :
    ((t * x) + x * (1 - (t - 1) * u * x) + (2 * u) * (t * x) * x) * (1 - t * u * x)
      = ((t + 1) * x) * ((1 - (t - 1) * u * x) - u ^ 2 * (t * x) * x) := by
  ring

end Ring

section Ordered

variable {S : Type*} [CommRing S] [LinearOrder S] [IsStrictOrderedRing S]

/-- On the turning side the slowest step is bounded below by the distance to the border:
`4 * q * (1 + k * a + q * a ^ 2) ≥ 4 * q - k ^ 2` for every `a`. -/
theorem slowest_step (k q a : S) : 4 * q - k ^ 2 ≤ 4 * q * (1 + k * a + q * a ^ 2) := by
  rw [wall_square]
  exact le_add_of_nonneg_left (sq_nonneg _)

/-- On the turning side the step never stalls: if `q > 0` and `4 * q > k ^ 2`, the step size is
positive. -/
theorem turning_never_stalls (k q a : S) (hq : 0 < q) (hd : k ^ 2 < 4 * q) :
    0 < 1 + k * a + q * a ^ 2 := by
  have h : 0 < 4 * q * (1 + k * a + q * a ^ 2) :=
    (sub_pos.2 hd).trans_le (slowest_step k q a)
  exact (pos_iff_pos_of_mul_pos h).1 (by positivity)

end Ordered

/-! ### Tangent addition with step `1`

Tangent addition (`k = 0`, `q = 1`) with step `1` is `a ↦ (a + 1) / (1 - a)`, the Möbius map of
the matrix `!![1, 1; -1, 1]`. Its fourth power is the scalar `-4`: every orbit spills through
infinity and returns to its start in exactly four steps (`tan (θ + π / 4)`). -/

/-- The tangent step squared is `!![0, 2; -2, 0]`. -/
theorem tan_square : (!![1, 1; -1, 1] : Matrix (Fin 2) (Fin 2) ℤ) ^ 2 = !![0, 2; -2, 0] := by
  simp [pow_two]

/-- The tangent step to the fourth is the scalar `-4`. -/
theorem tan_fourth :
    (!![1, 1; -1, 1] : Matrix (Fin 2) (Fin 2) ℤ) ^ 4
      = (-4 : ℤ) • (1 : Matrix (Fin 2) (Fin 2) ℤ) := by
  rw [show 4 = 2 * 2 from rfl, pow_mul, tan_square]
  simp [pow_two, one_fin_two]

end LeeYangSpill
