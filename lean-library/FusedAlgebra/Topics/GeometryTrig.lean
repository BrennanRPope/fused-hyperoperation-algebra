/-
  Fused Hyperoperation Algebra: geometry and trigonometry out of mult-add.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the library and NumberPlane.

  Geometry:  the corner square ab of (1 + a)(1 + b) is Euclid's gnomon (Elements II.4);
             mult-add of a with itself is the L-shaped strip that grows a square,
             mult-add three times is the cube's layers 3, 3, 1; Galileo's falling-body
             distances 1, 3, 5, 7 are exactly those gnomons.
  Pythagoras: 1 + x·i is the right triangle with legs 1 and x. Mult-adding two of them adds
             their angles; the slopes combine by tangent addition (the landscape point (0, 1))
             and the squared hypotenuses multiply. Squaring one gives Euclid's triples.
  Trigonometry: the sine and cosine addition laws are the multiplication at (0, 1);
             cubing 1 + x·i sends the cube's layers 1, 3, 3, 1 alternately to the two axes
             (triple-angle formula); the law of cosines is the norm on the line q = 1, and
             real triangles exist exactly in the turning country of that line.

  LIBRARY NOTE: copied from simple-laws/lean/GeometryTrig.lean (Nov 25) as Topics/GeometryTrig.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.NumberPlane

namespace SimpleLaws
open CircleOp

/-! ## Geometry -/

/-- the rectangle: (1 + a)(1 + b) = the unit square + two strips + the corner square -/
theorem rectangle_corner (a b : Int) : (1 + a) * (1 + b) = 1 + a + b + a * b := by grind

/-- Euclid II.4, the gnomon: mult-add of a with itself is the L-shaped strip 2a + a² -/
theorem gnomon (a : Int) : circ a a = 2 * a + a * a := by unfold circ; grind

/-- Euclid II.4 in full: (a + b)² = a² + 2ab + b² -/
theorem euclid_II_4 (a b : Int) : (a + b) * (a + b) = a * a + 2 * a * b + b * b := by grind

/-- the cube: mult-add three times gives the layers 3, 3, 1 (with the unit cube, 1, 3, 3, 1) -/
theorem cube_layers (a : Int) : circ (circ a a) a = 3 * a + 3 * (a * a) + a * a * a := by
  unfold circ; grind

/-- the box: three different sides give the corners a+b+c, ab+bc+ca, abc -/
theorem box_corners (a b c : Int) :
    circ (circ a b) c = (a + b + c) + (a * b + b * c + c * a) + a * b * c := by
  unfold circ; grind

/-- Galileo's square-cube law (1638): grow every length by a; area grows by A = a ⊞ a and
    volume by A ⊞ a, so weight per unit of strength grows by the factor 1 + a -/
theorem square_cube (a : Int) : 1 + circ (circ a a) a = (1 + circ a a) * (1 + a) :=
  circ_lift _ _

/-- Galileo's odd-number rule: the distances fallen in successive equal times -/
def fall : Nat → Nat
  | 0 => 0
  | n + 1 => 2 * n + 1

/-- total distance after n equal times -/
def fallen : Nat → Nat
  | 0 => 0
  | n + 1 => fallen n + fall (n + 1)

/-- each step's distance is the gnomon that grows an n-square into an (n+1)-square -/
theorem fall_is_gnomon (n : Nat) : (n + 1) * (n + 1) = n * n + fall (n + 1) := by
  simp only [fall]; grind

/-- so the distance fallen is the square of the time: d ∝ t² -/
theorem fallen_square (n : Nat) : fallen n = n * n := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [fallen, fall]; rw [ih]; grind

/-- the gnomon step is mult-add's lift on Nat: (n ⊞ n) + 1 = (n + 1)² -/
theorem gnomon_lift (n : Nat) : circ n n + 1 = (n + 1) * (n + 1) := circ_lift' n n

/-- Pythagoras by the "bride's chair": the square on a + b holds the tilted square on c and
    four right triangles of area ab/2. Given that area count, c² = a² + b². -/
theorem pythagoras_from_areas (a b c : Int) (h : (a + b) * (a + b) = c * c + 2 * a * b) :
    c * c = a * a + b * b := by grind

/-- the right triangle 1 + x·i: mult-adding two of them multiplies squared hypotenuses,
    (1 + x²)(1 + y²) = (1 − xy)² + (x + y)²; its new slope is (x + y)/(1 − xy) -/
theorem triangles_compose (x y : Int) :
    (1 + x * x) * (1 + y * y) = (1 - x * y) * (1 - x * y) + (x + y) * (x + y) := by grind

theorem triangles_compose_norm (x y : Int) :
    norm 0 1 (mre 1 1 x 1 y) (mu 0 1 x 1 y) = norm 0 1 1 x * norm 0 1 1 y := norm_mul 0 1 1 x 1 y

/-- Euclid's triples: square the triangle q + p·i (mult-add power 2 in the complex lift) -/
theorem euclid_triple_parts (p q : Int) : mre 1 q p q p = q * q - p * p ∧ mu 0 q p q p = 2 * p * q := by
  unfold mre mu; constructor <;> grind

theorem euclid_triple (p q : Int) :
    (q * q - p * p) * (q * q - p * p) + (2 * p * q) * (2 * p * q) = (q * q + p * p) * (q * q + p * p) := by
  grind

example : 3 * 3 + 4 * 4 = 5 * 5 ∧ 5 * 5 + 12 * 12 = 13 * 13 := by decide

/-! ## Trigonometry -/

/-- angle addition: (cos α + i sin α)(cos β + i sin β) at the complex point (0, 1) -/
theorem angle_addition (c1 s1 c2 s2 : Int) :
    mre 1 c1 s1 c2 s2 = c1 * c2 - s1 * s2 ∧ mu 0 c1 s1 c2 s2 = c1 * s2 + s1 * c2 := by
  unfold mre mu; constructor <;> grind

/-- cos² + sin² = 1 survives angle addition: the unit circle is closed under the product -/
theorem unit_circle_closed (c1 s1 c2 s2 : Int)
    (h1 : c1 * c1 + s1 * s1 = 1) (h2 : c2 * c2 + s2 * s2 = 1) :
    (c1 * c2 - s1 * s2) * (c1 * c2 - s1 * s2) + (c1 * s2 + s1 * c2) * (c1 * s2 + s1 * c2) = 1 := by
  have h := norm_mul 0 1 c1 s1 c2 s2
  unfold norm mre mu at h
  grind

/-- tangent addition is the landscape law at (0, 1): numerator x + y, denominator 1 − xy -/
theorem tangent_addition (x y : Int) :
    LandscapeMap.num 0 x y = x + y ∧ LandscapeMap.den 1 x y = 1 - x * y := by
  unfold LandscapeMap.num LandscapeMap.den circK; constructor <;> grind

/-- double angle: (1 + x i)² = (1 − x²) + 2x·i, so tan 2θ = 2t/(1 − t²) -/
theorem double_angle (x : Int) : mre 1 1 x 1 x = 1 - x * x ∧ mu 0 1 x 1 x = 2 * x := by
  unfold mre mu; constructor <;> grind

/-- triple angle: (1 + x i)³ = (1 − 3x²) + (3x − x³)·i: the cube's layers 1, 3, 3, 1
    land alternately on the real and the i axis, with the signs of i², i³ -/
theorem triple_angle (x : Int) :
    mre 1 (mre 1 1 x 1 x) (mu 0 1 x 1 x) 1 x = 1 - 3 * (x * x) ∧
    mu 0 (mre 1 1 x 1 x) (mu 0 1 x 1 x) 1 x = 3 * x - x * x * x := by
  unfold mre mu; constructor <;> grind

/-- the law of cosines is the norm on the line q = 1, with k = −2 cos C:
    c² = a² + k·ab + b². At k = 0 it is Pythagoras. -/
theorem law_of_cosines_norm (k a b : Int) : norm k 1 a b = a * a + k * a * b + b * b := by
  unfold norm; grind

theorem law_of_cosines_right_angle (a b : Int) : norm 0 1 a b = a * a + b * b := norm_complex a b

/-- on the border (k = ±2, angle 180° or 0°) the triangle is flat: c = a + b or c = a − b -/
theorem flat_triangle_open (a b : Int) : norm 2 1 a b = (a + b) * (a + b) := by unfold norm; grind
theorem flat_triangle_shut (a b : Int) : norm (-2) 1 a b = (a - b) * (a - b) := by unfold norm; grind

/-- inside the turning country (k² ≤ 4) every side-squared is a true length, never negative -/
theorem triangle_lengths_real (k a b : Int) (h : k * k ≤ 4) : 0 ≤ 4 * norm k 1 a b :=
  norm_nonneg_turning k 1 a b (by omega)

/-- past the border (k² > 4, "cos C" beyond ±1) a side-squared can go negative:
    the triangles of relativity, where cosine becomes cosh -/
example : norm 3 1 1 (-1) = -1 := by decide

end SimpleLaws
