/-
DRAFT (Nov 26, 2026). First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.JetRiver.
The river as a jet field: the cheap algebra behind jet-river/README.md.

The river cell's velocity gradient, in principal axes, is diag(a, b, b) with a = v' (radial)
and b = v/r (sideways). Its second layer Q = e₂ is what the field equations of the river
reduce to (checked symbolically in checks.py, J2). Here: the layer identities that need no
calculus. Builds on FusedAlgebra.Topics.JetFields (M3, M3.e2, diag3).
-/
import FusedAlgebra.Topics.JetFields

namespace CircleOp.JetRiver

open CircleOp
open CircleOp.JetFields
open Lean.Grind

section
variable {R : Type} [CommRing R]

/-- The river cell's second layer: Q = 2ab + b². -/
theorem river_e2 (a b : R) : (diag3 a b b).e2 = 2 * a * b + b * b := by
  simp only [diag3, M3.e2]; grind

/-- Inflow or outflow gives the same Q (the layer is quadratic). -/
theorem river_e2_neg (a b : R) : (diag3 (-a) (-b) (-b)).e2 = (diag3 a b b).e2 := by
  simp only [diag3, M3.e2]; grind

/-- **Q is the derivative of r·v², over r².** With b·r = v (sideways rate = v/r) and a = v':
    r²·Q = v² + 2 r v v' = (r v²)'. So Q = 0 is r v² = constant. -/
theorem river_e2_gauss (a b r v : R) (hb : b * r = v) :
    r * r * (diag3 a b b).e2 = v * v + 2 * r * v * a := by
  subst hb; simp only [diag3, M3.e2]; grind

/-- **Empty space:** if the radial rate is minus half the sideways rate (v ∝ r^(-1/2),
    that is v² = 2M/r), the cell's second layer vanishes. -/
theorem vacuum_e2 (a b : R) (h : 2 * a = -b) : (diag3 a b b).e2 = 0 := by
  simp only [diag3, M3.e2]
  have h2 : a * b * 2 = -(b * b) := by
    have : (2 * a) * b = (-b) * b := by rw [h]
    grind
  grind

/-- **A uniform ball:** every rate equal (v = H r), Q = 3H². With Q = 8πGρ this is
    Friedmann's equation H² = 8πGρ/3. -/
theorem friedmann_e2 (H : R) : (diag3 H H H).e2 = 3 * (H * H) := by
  simp only [diag3, M3.e2]; grind

/-- X + t·Y, entrywise. -/
def M3.addSmul (X Y : M3 R) (t : R) : M3 R :=
  ⟨X.m11 + t*Y.m11, X.m12 + t*Y.m12, X.m13 + t*Y.m13, X.m21 + t*Y.m21, X.m22 + t*Y.m22,
   X.m23 + t*Y.m23, X.m31 + t*Y.m31, X.m32 + t*Y.m32, X.m33 + t*Y.m33⟩

/-- **The gradient of the Q layer.** Q(X + tY) = Q(X) + t·(tr X · tr Y − tr(XY)) + t²·Q(Y).
    The linear term is ⟨(tr X)I − Xᵀ, Y⟩: the vector law of the river (momentum constraint)
    is the divergence of (tr S)I − S. -/
theorem e2_addSmul (X Y : M3 R) (t : R) :
    (M3.addSmul X Y t).e2 = X.e2 + t * (X.tr * Y.tr - (X.mul Y).tr) + t^2 * Y.e2 := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  obtain ⟨a', b', c', d', e', f', g', k', i'⟩ := Y
  simp only [M3.addSmul, M3.e2, M3.tr, M3.mul]; grind

end

end CircleOp.JetRiver
