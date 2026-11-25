/-
  BooleanLattice: the circle operation lives on the corners of a cube (Nov 25, 2026)
  =============================================================================
  A chain of n circle operations is a sum over the 2ⁿ − 1 non-empty corners of the
  n-cube (the Boolean lattice Bₙ). The dial k weights a corner of rank m by k^(m−1);
  the plane's two dials (k, d) weight it by the Lucas numbers U_m(k, d). At k = −1
  the chain is inclusion–exclusion. In a poset's incidence algebra, the way back of
  the strict zeta function is the Möbius function (Rota, Hall).

  Companion: math-theory/boolean-lattice/README.md, checks.py (92 exact checks).
  Core Lean 4.37, imports the library's Base for `circK`.
-/
import FusedAlgebra.Base

namespace BooleanLattice
open CircleOp

section Comm
variable {R : Type} [Lean.Grind.CommRing R]

/-- Four steps: every non-empty corner of the 4-cube (the tesseract) once,
    a rank-m corner weighted k^(m−1). -/
theorem chain4_corners (k a b c d : R) :
    circK k (circK k (circK k a b) c) d =
      (a + b + c + d) + k * (a*b + a*c + a*d + b*c + b*d + c*d)
        + k * k * (a*b*c + a*b*d + a*c*d + b*c*d) + k * k * k * (a*b*c*d) := by
  unfold circK; grind

/-- k = −1 is inclusion–exclusion: the chance that at least one of three
    independent events happens. -/
theorem incl_excl3 (p q r : R) :
    circK (-1) (circK (-1) p q) r = 1 - (1 - p) * (1 - q) * (1 - r) := by
  unfold circK; grind

/-- The plane's law (x + y + k·xy)/(1 − d·xy) with k = α + β and d = αβ turns into
    multiplication under x ↦ (1 + αx)/(1 + βx). Stated cross-multiplied, so it
    holds in any commutative ring. -/
theorem plane_lift (α β x y : R) :
    let N := x + y + (α + β) * x * y
    let D := 1 - (α * β) * x * y
    (D + α * N) * ((1 + β * x) * (1 + β * y)) = ((1 + α * x) * (1 + α * y)) * (D + β * N) := by
  intro N D; grind

/-- Three steps on the plane, numerator: corners weighted by the Lucas numbers
    U₁ = 1, U₂ = k, U₃ = k² − d. (The two-step law is N/D; one more step with z
    has numerator N + D·z + k·N·z.) -/
theorem plane_chain3_num (k d x y z : R) :
    let N := x + y + k * x * y
    let D := 1 - d * x * y
    N + D * z + k * N * z = (x + y + z) + k * (x*y + x*z + y*z) + (k * k - d) * (x*y*z) := by
  intro N D; grind

/-- Three steps on the plane, denominator: 1 − d·(U₁·e₂ + U₂·e₃). -/
theorem plane_chain3_den (k d x y z : R) :
    let N := x + y + k * x * y
    let D := 1 - d * x * y
    D - d * N * z = 1 - d * ((x*y + x*z + y*z) + k * (x*y*z)) := by
  intro N D; grind

/-- The tangent of a sum of three angles, (k, d) = (0, 1): odd corners with
    alternating signs over even corners. -/
theorem tangent3 (x y z : R) :
    let N := x + y
    let D := 1 - x * y
    (N + D * z) = (x + y + z) - x*y*z ∧ (D - N * z) = 1 - (x*y + x*z + y*z) := by
  intro N D; constructor <;> grind

/-- Relativistic velocity addition, (k, d) = (0, −1): odd corners over even corners. -/
theorem velocity3 (x y z : R) :
    let N := x + y
    let D := 1 + x * y
    (N + D * z) = (x + y + z) + x*y*z ∧ (D + N * z) = 1 + (x*y + x*z + y*z) := by
  intro N D; constructor <;> grind

end Comm

section NonComm
variable {R : Type} [Lean.Grind.Ring R]

/-- The circle operation with dial k in a ring where order matters (k written on the left). -/
def ncirc (k a b : R) : R := a + b + k * a * b

/-- Three steps when order matters: every non-empty corner of the cube once, its
    factors in chain order. -/
theorem nchain3_corners (k a b c : R) :
    ncirc k (ncirc k a b) c = (a + b + c) + k * (a * b + a * c + b * c) + k * k * (a * b * c) := by
  unfold ncirc; grind

/-- Rota and Hall in one line: if no chain is longer than 2 steps (η³ = 0, as for the
    strict zeta function of a poset of rank 2), the way back of η under the circle
    operation is −η + η², the alternating count of chains; so 1 + (−η + η²) is the
    Möbius function. (For dial k, replace η by k·η.) -/
theorem way_back_of_strict_zeta (η : R) (h : η * η * η = 0) :
    ncirc 1 η (-η + η * η) = 0 := by
  unfold ncirc
  -- grind is given the identity with η generalised, since the hypothesis h distracts it
  have e : ∀ x : R, x + (-x + x * x) + 1 * x * (-x + x * x) = x * x * x := by intro x; grind
  rw [e, h]

end NonComm
end BooleanLattice
