/-
  Fused Hyperoperation Algebra: developmental biology (growth curves, allometry).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Development.lean

  1. The Richards family of growth curves: size/asymptote = 1 + mpow(x, p) with the excess x
     decaying exponentially. Logistic growth is p = −1, the mult-add inverse (saturation);
     von Bertalanffy is p = 3.
  2. Gompertz in discrete time iterates one mult-add power; powers compose by multiplying
     exponents, so n steps are a single power.
  3. Allometry y ∝ xᵃ: the growth of y is the mult-add power of the growth of x.
  LIBRARY NOTE: first compiled in the shared library. 3 exact twins removed (star_def,
  mpow_mpow, logistic); each is marked where it stood. pow_mul' is a twin of
  Dynamics.pow_mul_R but stays, because later proofs here use it. No proof step or
  statement changed.
-/
import FusedAlgebra.Base

namespace Development
open MultAdd Lean.Grind

set_option linter.unusedSectionVars false

/-- Brennan's notation for mult-add in this file: a ✳ b = a + b + ab. -/
local infixl:65 " ✳ " => multAdd

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

-- LIBRARY: pow_mul' is an exact twin of Dynamics.pow_mul_R, kept because later proofs here use it and Dynamics is not imported.
theorem pow_mul' (a : R) (m n : Nat) : a ^ (m * n) = (a ^ m) ^ n := by
  induction n with
  | zero => rw [Nat.mul_zero, Semiring.pow_zero, Semiring.pow_zero]
  | succ n ih => rw [Nat.mul_succ, pow_add', ih, Semiring.pow_succ]

-- LIBRARY: mpow_mpow (**Powers compose.** mpow(mpow(x, a), b) = mpow(x, a·b)) was an exact twin of Dynamics.mpow_mpow and is removed.

/-- Iterating one mult-add power n times. -/
def iter (l : Nat) : Nat → R → R
  | 0, x => x
  | n + 1, x => mpow (iter l n x) l

/-- **Gompertz steps.** n steps of x ↦ mpow(x, λ) give the single power λⁿ. -/
theorem gompertz_steps (l n : Nat) (x : R) : 1 + iter l n x = (1 + x) ^ (l ^ n) := by
  induction n with
  | zero => simp only [iter]; rw [Nat.pow_zero, Semiring.pow_one]
  | succ n ih =>
    simp only [iter]
    rw [mpow_lift, ih, ← pow_mul', Nat.pow_succ]

/-- Von Bertalanffy growth (power 3): size/asymptote = (1 + x)³ = 1 + mpow(x, 3). -/
theorem bertalanffy (x : R) : (1 + x) ^ 3 = 1 + mpow x 3 := (mpow_lift x 3).symm

/-- Allometry: if x grows by the factor 1 + g, then y ∝ xᵃ grows by 1 + mpow(g, a). -/
theorem allometry (b g : R) (a : Nat) (x : R) :
    b * ((1 + g) * x) ^ a = (1 + mpow g a) * (b * x ^ a) := by
  rw [mpow_lift, mul_pow']; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

-- LIBRARY: logistic (**Logistic growth is saturation by the inverse.** With excess x = Q·e^(−rt), size/asymptote = 1/(1 + x) = 1 + inv(x), the Michaelis–Menten form again) was an exact twin of Metabolism.saturation and is removed.

/-- The logistic shortfall below the asymptote is −inv(x) = x/(1 + x). -/
theorem logistic_shortfall (x : K) (_h : 1 + x ≠ 0) : 1 - 1 / (1 + x) = -multAddInv x := by
  unfold multAddInv; grind

end Field

end Development
