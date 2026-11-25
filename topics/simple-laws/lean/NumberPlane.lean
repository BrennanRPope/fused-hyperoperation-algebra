/-
  Fused Hyperoperation Algebra: every point of the landscape is mult-add in a number system.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library.

  Pick a new unit u with u² = k·u − q (so u + ū = k and u·ū = q). Numbers are a + b·u.
  Mult-add the two "tilted ones" 1 + x·u and 1 + y·u:
      (1 + x u)(1 + y u) = (1 − q x y) + (x + y + k x y)·u
  Read the result as a slope (u-part over real part) and you get exactly the landscape law
  at (k, q):  (x + y + k x y) / (1 − q x y).
  The norm N(a + b u) = a² + k a b + q b² multiplies, so each point carries its own geometry:
    turning country (k² < 4q)  : N is a sum of squares in disguise  -> Euclid (circles, rotation)
    border          (k² = 4q)  : N is one square                    -> Galileo (no length one way)
    two walls       (k² > 4q)  : N splits into two factors          -> Einstein (light cone = walls)
  Classical: Clifford (1873), Study, Yaglom "A Simple Non-Euclidean Geometry" (1979).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Landscape

namespace SimpleLaws
open MultAdd

/-! ### 1. The number system at (k, q): a + b·u with u² = k·u − q -/

/-- real part of (a + b u)(c + d u) -/
def mre (q a b c d : Int) : Int := a * c - q * b * d
/-- u-part of (a + b u)(c + d u) -/
def mu (k a b c d : Int) : Int := a * d + b * c + k * b * d

/-- the multiplication really is "u² = k u − q": expand (a + b u)(c + d u) and replace u² -/
theorem mul_rule (k q a b c d u : Int) (hu : u * u = k * u - q) :
    (a + b * u) * (c + d * u) = mre q a b c d + mu k a b c d * u := by
  unfold mre mu
  have : (a + b * u) * (c + d * u) = a * c + (a * d + b * c) * u + b * d * (u * u) := by grind
  rw [this, hu]; grind

theorem mul_comm_re (q a b c d : Int) : mre q a b c d = mre q c d a b := by unfold mre; grind
theorem mul_comm_u (k a b c d : Int) : mu k a b c d = mu k c d a b := by unfold mu; grind

/-- associativity, real part and u-part -/
theorem mul_assoc_re (k q a b c d e f : Int) :
    mre q (mre q a b c d) (mu k a b c d) e f = mre q a b (mre q c d e f) (mu k c d e f) := by
  unfold mre mu; grind
theorem mul_assoc_u (k q a b c d e f : Int) :
    mu k (mre q a b c d) (mu k a b c d) e f = mu k a b (mre q c d e f) (mu k c d e f) := by
  unfold mre mu; grind

/-! ### 2. The main fact: mult-add of tilted ones is the landscape law -/

/-- (1 + x u)(1 + y u): its real part is the landscape denominator -/
theorem tilted_re (q x y : Int) : mre q 1 x 1 y = LandscapeMap.den q x y := by
  unfold mre LandscapeMap.den; grind

/-- ... and its u-part is the landscape numerator, the weighted mult-add of x and y -/
theorem tilted_u (k x y : Int) : mu k 1 x 1 y = LandscapeMap.num k x y := by
  unfold mu LandscapeMap.num wAdd; grind

/-- at mult-add's own point (k, q) = (1, 0) the unit is idempotent (u² = u) and the
    u-part is plain mult-add, with real part 1: nothing to divide by -/
theorem multadd_point_u (x y : Int) : mu 1 1 x 1 y = multAdd x y := by
  unfold mu multAdd; grind
theorem multadd_point_re (x y : Int) : mre 0 1 x 1 y = 1 := by unfold mre; grind

/-- the complex corner (0, 1), u = i: slopes add like tangents, (x + y)/(1 − x y) -/
theorem complex_point (x y : Int) : mu 0 1 x 1 y = x + y ∧ mre 1 1 x 1 y = 1 - x * y := by
  unfold mu mre; constructor <;> grind

/-- the split corner (0, −1), u² = 1: velocity addition, (x + y)/(1 + x y) -/
theorem split_point (x y : Int) : mu 0 1 x 1 y = x + y ∧ mre (-1) 1 x 1 y = 1 + x * y := by
  unfold mu mre; constructor <;> grind

/-- the dual corner (0, 0), u² = 0: Galileo's velocity addition, plain x + y over 1 -/
theorem dual_point (x y : Int) : mu 0 1 x 1 y = x + y ∧ mre 0 1 x 1 y = 1 := by
  unfold mu mre; constructor <;> grind

/-! ### 3. The norm: each point's own Pythagoras -/

def norm (k q a b : Int) : Int := a * a + k * a * b + q * b * b

/-- the norm multiplies (the generalized Brahmagupta–Fibonacci identity) -/
theorem norm_mul (k q a b c d : Int) :
    norm k q (mre q a b c d) (mu k a b c d) = norm k q a b * norm k q c d := by
  unfold norm mre mu; grind

/-- the norm is the number times its conjugate: with ū = k − u, (a + b u)(a + b ū) = N -/
theorem norm_conj (k q a b u : Int) (hu : u * u = k * u - q) :
    (a + b * u) * (a + b * (k - u)) = norm k q a b := by
  unfold norm
  have : (a + b * u) * (a + b * (k - u)) = a * a + k * a * b + b * b * (k * u - u * u) := by grind
  rw [this, hu]; grind

/-- 4·N splits as a square plus (4q − k²)·b²: the sign of 4q − k² picks the country -/
theorem four_norm (k q a b : Int) :
    4 * norm k q a b = (2 * a + k * b) * (2 * a + k * b) + (4 * q - k * k) * (b * b) := by
  unfold norm; grind

/-- a square is never negative (core Lean has no lemma by this name, so we prove it) -/
theorem sq_nonneg (x : Int) : 0 ≤ x * x := by
  rcases Int.le_total 0 x with h | h
  · exact Int.mul_nonneg h h
  · exact Int.mul_nonneg_of_nonpos_of_nonpos h h

/-- turning country: when 4q ≥ k² the norm is never negative (a true length) -/
theorem norm_nonneg_turning (k q a b : Int) (h : k * k ≤ 4 * q) : 0 ≤ 4 * norm k q a b := by
  rw [four_norm]
  have h1 : 0 ≤ (2 * a + k * b) * (2 * a + k * b) := sq_nonneg _
  have h2 : 0 ≤ 4 * q - k * k := by omega
  have h3 : 0 ≤ (4 * q - k * k) * (b * b) := Int.mul_nonneg h2 (sq_nonneg b)
  omega

/-- border: when k² = 4q the norm is one square, (2a + k b)²/4 (Galileo: lengths only one way) -/
theorem norm_border (k q a b : Int) (h : k * k = 4 * q) :
    4 * norm k q a b = (2 * a + k * b) * (2 * a + k * b) := by
  rw [four_norm, ← h]; grind

/-- two walls: write k = −(r + s), q = r s; the norm splits into two light-cone factors -/
theorem norm_walls (r s a b : Int) : norm (-(r + s)) (r * s) a b = (a - r * b) * (a - s * b) := by
  unfold norm; grind

/-- Pythagoras at the complex corner: a² + b² -/
theorem norm_complex (a b : Int) : norm 0 1 a b = a * a + b * b := by unfold norm; grind
/-- the spacetime interval at the split corner: a² − b² -/
theorem norm_split (a b : Int) : norm 0 (-1) a b = a * a - b * b := by unfold norm; grind
/-- Galileo at the dual corner: only a² survives -/
theorem norm_dual (a b : Int) : norm 0 0 a b = a * a := by unfold norm; grind

/-! ### 4. Mult-add's own geometry: its unit circle is Boyle's hyperbola -/

/-- at (1, 0) the norm is a·(a + b): in the coordinates X = a, Y = a + b it is X·Y -/
theorem norm_multadd (a b : Int) : norm 1 0 a b = a * (a + b) := by unfold norm; grind

/-- the tilted one 1 + x u at mult-add's point has norm 1·(1 + x): its lift -/
theorem norm_multadd_tilted (x : Int) : norm 1 0 1 x = 1 + x := by unfold norm; grind

/-- the product of two tilted ones has norm (1 + x)(1 + y) = 1 + (x ⊞ y) -/
theorem norm_multadd_lift (x y : Int) :
    norm 1 0 (mre 0 1 x 1 y) (mu 1 1 x 1 y) = 1 + multAdd x y := by
  unfold norm mre mu multAdd; grind

/-- the unit circle N = 1 at mult-add's point is the curve x ⊞ y = 0, i.e. (1 + x)(1 + y) = 1:
    a hyperbola whose asymptotes are the walls x = −1 and y = −1 -/
theorem unit_circle_is_hyperbola (x y : Int) : multAdd x y = 0 ↔ (1 + x) * (1 + y) = 1 := by
  rw [← multAdd_lift]; omega

end SimpleLaws
