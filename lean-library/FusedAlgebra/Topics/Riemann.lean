/-
  Fused Hyperoperation Algebra: the Riemann hypothesis, first look (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib), no `sorry`. First compiled in the shared library on Nov 26,
  as FusedAlgebra.Topics.Riemann (old names wAdd and foldMA updated to circK and foldCirc).

  What is here (all finite, algebraic pieces of the write-up; nothing here touches the
  zeros of zeta, which need analysis that core Lean does not have):

  1. The prime cube fold.  Folding ⊞ₖ over the weights x_p = p^(-s) lifts to the Euler
     product: 1 + k·(x₁ ⊞ₖ … ⊞ₖ xₙ) = ∏ (1 + k·xᵢ).  k = 1 is ζ(s)/ζ(2s) on the squarefree
     corners, k = −1 is 1/ζ(s) (Möbius signs).
  2. Boson cube times fermion cube.  ∏(1 + x)·∏(1 − x) = ∏(1 − x²): the k = 1 and k = −1
     cubes multiply to the k = −1 cube at doubled s (finite form of ζ/ζ(2s) · 1/ζ = 1/ζ(2s)).
  3. The local factor as a dial point.  1 + k·u + d·u² : its roots come in reciprocal pairs
     u ↦ 1/(d·u) (the functional-equation symmetry), and on the elliptic side k² < 4d it has
     no integer root at all (finite shadow of "both roots on the circle |u| = d^(-1/2)").
  4. The step identity's kernel side and the heat step: the two algebraic facts behind
     sections 4 and 5 of the README.
-/
import FusedAlgebra.Base

namespace CircleOp.Riemann
open CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The prime cube fold -/

section Semiring
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- Fold ⊞ₖ over a list of corner weights. -/
def cubeFold (k : R) : List R → R
  | [] => 0
  | a :: l => circK k a (cubeFold k l)

/-- The fold lifts to the Euler product: 1 + k·fold = ∏ (1 + k·x). -/
theorem cubeFold_lift (k : R) (l : List R) :
    1 + k * cubeFold k l = prodL (l.map (fun x => 1 + k * x)) := by
  induction l with
  | nil => simp only [cubeFold, List.map, prodL]; grind
  | cons a l ih => simp only [cubeFold, List.map, prodL]; rw [← ih, circK_lift]

/-- At k = 1 the fold is the plain circle operation. -/
theorem cubeFold_one (l : List R) : cubeFold 1 l = foldCirc l := by
  induction l with
  | nil => rfl
  | cons a l ih => simp only [cubeFold, foldCirc]; rw [ih, circK_one]

/-- Products of lists multiply pointwise. -/
theorem prodL_map_mul {α : Type} (f g : α → R) (l : List α) :
    prodL (l.map f) * prodL (l.map g) = prodL (l.map (fun x => f x * g x)) := by
  induction l with
  | nil => simp only [List.map, prodL]; grind
  | cons a l ih => simp only [List.map, prodL]; rw [← ih]; grind

end Semiring

/-! ## 2. Boson cube times fermion cube -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- ∏(1 + x) · ∏(1 − x) = ∏(1 − x²), written with the two cube folds:
    (1 + fold₁ l)·(1 − fold₋₁ l) = 1 − fold₋₁ (l²). -/
theorem boson_times_fermion (l : List R) :
    (1 + 1 * cubeFold 1 l) * (1 + (-1) * cubeFold (-1) l)
      = 1 + (-1) * cubeFold (-1) (l.map (fun x => x * x)) := by
  rw [cubeFold_lift, cubeFold_lift, cubeFold_lift, prodL_map_mul, List.map_map]
  congr 1
  congr 1
  funext x
  simp only [Function.comp]
  grind

/-- Kernel side of the Lee–Yang step identity, with y = xᵖ and c = (a+1)p:
    c·y·(1 − y)^a + c·(1 − y)^(a+1) = c·(1 − y)^a.
    (−x·d/dx of (1 − xᵖ)^(a+1) is the first term; this is the whole identity.) -/
theorem step_identity_kernel (c y : R) (a : Nat) :
    c * y * (1 - y) ^ a + c * (1 - y) ^ (a + 1) = c * (1 - y) ^ a := by
  rw [Lean.Grind.Semiring.pow_succ]; grind

/-- One heat step: multiplying a kernel by 1 + b²u² acts on its transform as
    1 − b²D² = (1 − bD)(1 + bD), a product of two Laguerre steps. Algebraic core. -/
theorem heat_step_factor (b X : R) : (1 - b * X) * (1 + b * X) = 1 - b * b * (X * X) := by
  grind

end Ring

/-! ## 3. The local factor as a dial point -/

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- Roots of 1 + k·u + d·u² pair up as u and 1/(d·u): the reflection that fixes the
    circle |u| = d^(-1/2) swaps them. -/
theorem localFactor_reciprocal (k d u : K) (hd : d ≠ 0)
    (h : 1 + k * u + d * u * u = 0) : 1 + k * (1 / (d * u)) + d * (1 / (d * u)) * (1 / (d * u)) = 0 := by
  have hu : u ≠ 0 := by
    intro hu0; rw [hu0] at h; grind
  have hdu : d * u ≠ 0 := by
    intro h0; grind
  have key : (d * u * u) * (1 + k * (1 / (d * u)) + d * (1 / (d * u)) * (1 / (d * u)))
      = 1 + k * u + d * u * u := by
    grind
  rw [h] at key
  have hduu : d * u * u ≠ 0 := by
    intro h0; grind
  grind

end Field

/-- Squares of integers are nonnegative (local helper). -/
theorem int_sq_nonneg (a : Int) : 0 ≤ a * a := by
  rcases Int.le_total 0 a with h | h
  · exact Int.mul_nonneg h h
  · have h' : 0 ≤ -a := Int.neg_nonneg_of_nonpos h
    have := Int.mul_nonneg h' h'
    rw [Int.neg_mul_neg] at this
    exact this

/-- Elliptic side, no integer root: if 0 < d and k² < 4d then 1 + k·u + d·u² ≠ 0
    for every integer u.  (4d·(1 + ku + du²) = (2du + k)² + (4d − k²) > 0.) -/
theorem elliptic_no_int_root (k d : Int) (_hd : 0 < d) (hk : k * k < 4 * d) (u : Int) :
    1 + k * u + d * u * u ≠ 0 := by
  intro h
  have key : 4 * d * (1 + k * u + d * u * u)
      = (2 * d * u + k) * (2 * d * u + k) + (4 * d - k * k) := by
    grind
  rw [h, Int.mul_zero] at key
  have hs := int_sq_nonneg (2 * d * u + k)
  omega

/-- The genus-1 curve y² = x³ + x + 1 at p = 5: a = p + 1 − #E = −3 (#E = 9),
    dial point (k, d) = (3, 5), Δ = 9 − 20 < 0: elliptic.  (Worked example, by decide.) -/
theorem genus1_p5_elliptic : (3 : Int) * 3 < 4 * 5 := by decide

end CircleOp.Riemann
