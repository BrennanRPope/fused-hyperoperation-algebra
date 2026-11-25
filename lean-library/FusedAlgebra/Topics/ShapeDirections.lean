/-
  The shape's directions: which dials of Fused Hyperoperation Algebra are truly separate
  ======================================================================================
  Companion page and checks: ../index.html, ../checks.py (29 checks, sympy + mpmath).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports only FusedAlgebra.Base.
  Intended path: FusedAlgebra/Topics/ShapeDirections.lean.
  LIBRARY NOTE: first compiled in the shared library, unchanged. No exact twins. Int
  specializations kept beside their general forms: star_def (CircleOp.circ_def),
  weight_is_units (LeeYangCluster.circK_units, ring-general).
  Not compiled in the session that wrote it (no Lean there); written in the library's
  `unfold …; grind` style for its first compile. Every identity was checked with sympy.

  Three levels of "the same law":
    units  (a ↦ λa: change the size of the measuring stick),
    rulers (any change of coordinate near 0: a lens),
    global (the whole group: how many periods its log has).

  What the file proves:
    A. The weight is units: λ·(a ◯[λk] b) = (λa) ◯[k] (λb). The q of the q-integers is k + 1,
       the same dial. On a graph (three sites in a row) the weight can no longer be scaled away.
    B. The landscape bend: a projective ruler x ↦ x/(1+cx) slides (p, q) along the parabola of
       constant discriminant D = p² − 4q; units multiply D by λ². So after units and rulers only
       the sign of D is left (two real walls: a line; complex walls: a circle).
    C. Rung height: the ladder E^s ∘ E = E ∘ E^s forces the jets of the fractional rung
       (c₃ = s(3s−1)/12, c₄ = s(2s−1)(3s−1)/48). From them the rung-s law has corner p = s,
       bend q = s(s−1)/4, so D = s (rung height = the bend's discriminant through third order),
       and a fourth-order coefficient −s(s−1)/12 that no point of the landscape plane has
       (zero only at rungs 0 and 1). In jet space the rung adds a direction; the weight does not.
    D. The nome: gluing a ↦ a ◯ (q − 1) commutes with the law for every q and the law descends
       to every quotient, so nothing near 0 can see q. Torsion through the 3-division polynomial:
       the cusp (addition) has no 3-step loop, the node (mult-add) has two plus the start.
    E. Several waves: a pair weight survives every separate rescaling of the waves. Riemann's
       period matrices (g(g+1)/2 numbers) outnumber curves (3g − 3) exactly from four waves on.
    F. Whole numbers: the p-fold mult-add power of a + b + k·ab is ≡ 0 mod p when p ∣ k
       (additive there) and ≡ x^p otherwise (multiplicative), shown for p = 2, 3, 5.
-/
import FusedAlgebra.Base

namespace ShapeDirections
open CircleOp

open scoped CircleOp

theorem star_def (a b : Int) : a ◯ b = a + b + a * b := rfl

/-! ## A. The weight is the size of the measuring stick -/

/-- Zooming the inputs by λ and the weight by 1/λ gives the same law: the weight is units. -/
theorem weight_is_units (k l a b : Int) : l * circK (l * k) a b = circK k (l * a) (l * b) := by
  unfold circK; grind

/- (Every nonzero weight is mult-add after a zoom: that is the base's `circK_scale`.) -/

/-- The q of the q-integers (and of the Potts colours) is the weight plus one. -/
theorem qcount_is_weight (q a b : Int) : circK (q - 1) a b = a + b + q * a * b - a * b := by
  unfold circK; grind

/-- Three sites in a row, each site with fugacity z and each occupied bond with weight k:
    Z = 1 + 3z + (2k+1)z² + k²z³. -/
def path3 (k z : Int) : Int := 1 + 3 * z + (2 * k + 1) * z * z + k * k * z * z * z

/-- The list of all eight occupation patterns adds up to the closed form. -/
theorem path3_sum (k z : Int) :
    path3 k z = 1 + z + z + z + (z * z * k) + (z * z * k) + (z * z) + (z * z * z * k * k) := by
  unfold path3; grind

/-- On the graph the weight survives: if a rescaling of the sites (z ↦ λz) and of the whole
    sum (by m) turns the coefficients for weight k into those for weight k', then k = k'. -/
theorem graph_weight_survives (k k' l m : Rat)
    (h0 : m * 1 = 1) (h1 : m * (3 * l) = 3) (h2 : m * ((2 * k + 1) * (l * l)) = 2 * k' + 1) :
    k = k' := by
  have hm : m = 1 := by grind
  subst hm
  have hl : l = 1 := by grind
  subst hl
  grind

/-! ## B. The landscape bend: rulers slide along parabolas, units zoom them -/

/-- Cross-multiplied numerator of F_{p,q}(N a, N b) with N(x) = x/(1+cx): it is the numerator
    of F_{p+2c, ·}(a, b). -/
theorem ruler_num (p c a b : Int) :
    a * (1 + c * b) + b * (1 + c * a) + p * a * b = circK (p + 2 * c) a b := by
  unfold circK; grind

/-- Cross-multiplied denominator: (1+ca)(1+cb) − q·ab is the denominator of N applied to
    F_{p+2c, q+pc+c²}(a, b). Together with ruler_num: N⁻¹(F_{p,q}(N a, N b)) = F_{p+2c, q+pc+c²}(a,b). -/
theorem ruler_den (p q c a b : Int) :
    (1 + c * a) * (1 + c * b) - q * a * b
      = (1 - (q + p * c + c * c) * a * b) + c * circK (p + 2 * c) a b := by
  unfold circK; grind

/-- The ruler keeps the discriminant: (p+2c)² − 4(q+pc+c²) = p² − 4q. The ruler's orbits are
    the parabolas q = (p² − D)/4, one for each D. -/
theorem ruler_keeps_disc (p q c : Int) :
    (p + 2 * c) * (p + 2 * c) - 4 * (q + p * c + c * c) = p * p - 4 * q := by grind

/-- Units zoom the discriminant: (λp)² − 4(λ²q) = λ²(p² − 4q). Only its sign is free of units. -/
theorem zoom_disc (p q l : Int) : (l * p) * (l * p) - 4 * (l * l * q) = l * l * (p * p - 4 * q) := by
  grind

/-! ## C. Rung height: the ladder forces the jets, the jets give D = s and a fourth-order gap

  A jet x + c₂x² + c₃x³ + c₄x⁴ + c₅x⁵ is stored as its four coefficients. -/

structure Jet where
  c2 : Rat
  c3 : Rat
  c4 : Rat
  c5 : Rat

/-- Composition f ∘ g of two jets, cut after x⁵ (orders 4 and 5 shown). -/
def comp4 (f g : Jet) : Rat := f.c2 * (g.c2 * g.c2 + 2 * g.c3) + 3 * f.c3 * g.c2 + f.c4 + g.c4
def comp5 (f g : Jet) : Rat :=
  f.c2 * (2 * g.c2 * g.c3 + 2 * g.c4) + f.c3 * (3 * g.c2 * g.c2 + 3 * g.c3) + 4 * f.c4 * g.c2 + f.c5 + g.c5

/-- E(x) = eˣ − 1 to fifth order. -/
def E : Jet := ⟨1/2, 1/6, 1/24, 1/120⟩

/-- The ladder at fourth order forces c₃ of any rung with c₂ = s/2. -/
theorem ladder_forces_c3 (s c3 c4 c5 : Rat)
    (h4 : comp4 E ⟨s/2, c3, c4, c5⟩ = comp4 ⟨s/2, c3, c4, c5⟩ E) :
    c3 = s * (3 * s - 1) / 12 := by
  unfold comp4 E at h4; simp only at h4; grind

/-- The ladder at fifth order then forces c₄. -/
theorem ladder_forces_c4 (s c3 c4 c5 : Rat)
    (h4 : comp4 E ⟨s/2, c3, c4, c5⟩ = comp4 ⟨s/2, c3, c4, c5⟩ E)
    (h5 : comp5 E ⟨s/2, c3, c4, c5⟩ = comp5 ⟨s/2, c3, c4, c5⟩ E) :
    c4 = s * (2 * s - 1) * (3 * s - 1) / 48 := by
  have h3 := ladder_forces_c3 s c3 c4 c5 h4
  unfold comp5 E at h5; simp only at h5
  rw [h3] at h5
  grind

/-- ...and the forced jets do satisfy the ladder (so the rung exists to this order). -/
theorem ladder_jets_commute (s c5 : Rat) :
    comp4 E ⟨s/2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩
      = comp4 ⟨s/2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩ E ∧
    comp5 E ⟨s/2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩
      = comp5 ⟨s/2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩ E := by
  unfold comp4 comp5 E; simp only; constructor <;> grind

/-- The rung-s law is plain addition in the ruler ℓ = E^{−s}, whose jets are the forced ones at −s:
    ℓ₂ = −s/2, ℓ₃ = s(3s+1)/12, ℓ₄ = −s(2s+1)(3s+1)/48. Its rate R = 1/ℓ′ has coefficients
    r₁ = s, r₂ = s(s−1)/4, r₃ = −s(s−1)/12: these satisfy ℓ′·R = 1 through a³. -/
theorem rung_rate_inverts (s : Rat) :
    s + 2 * (-s / 2) = 0 ∧
    s * (s - 1) / 4 + 2 * (-s / 2) * s + 3 * (s * (3 * s + 1) / 12) = 0 ∧
    -(s * (s - 1)) / 12 + 2 * (-s / 2) * (s * (s - 1) / 4) + 3 * (s * (3 * s + 1) / 12) * s
      + 4 * (-(s * (2 * s + 1) * (3 * s + 1)) / 48) = 0 := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- The jets of E^{−s} are the forced jets with s replaced by −s. -/
theorem lens_jets_at_minus (s : Rat) :
    (-s) * (3 * (-s) - 1) / 12 = s * (3 * s + 1) / 12 ∧
    (-s) * (2 * (-s) - 1) * (3 * (-s) - 1) / 48 = -(s * (2 * s + 1) * (3 * s + 1)) / 48 := by
  constructor <;> grind

/-- Through third order the rung's discriminant is its height: r₁² − 4r₂ = s.
    So rung s looks, to third order, like the law with one period 2πi/√s
    (rung 1: 2πi, mult-add; rung 0: no period; rung −1: real period 2π, a circle). -/
theorem rung_disc_is_height (s : Rat) : s * s - 4 * (s * (s - 1) / 4) = s := by grind

/-- Every point of the landscape plane has rate 1 + p a + q a² exactly, so r₃ = 0 there;
    the rung's r₃ = −s(s−1)/12 differs from 0 at s = 2, s = −1 and s = 1/2. -/
theorem rung_leaves_plane :
    -((2 : Rat) * (2 - 1)) / 12 ≠ 0 ∧ -((-1 : Rat) * (-1 - 1)) / 12 ≠ 0 ∧
    -((1/2 : Rat) * (1/2 - 1)) / 12 ≠ 0 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide +kernel

/-- The rung crosses the plane at rungs 0 and 1 (r₃ = 0 there), crossing, not touching:
    dr₃/ds = −(2s−1)/12 is +1/12 at s = 0 and −1/12 at s = 1. -/
theorem rung_crossings :
    -((0 : Rat) * (0 - 1)) / 12 = 0 ∧ -((1 : Rat) * (1 - 1)) / 12 = 0 ∧
    -(2 * (0 : Rat) - 1) / 12 = 1/12 ∧ -(2 * (1 : Rat) - 1) / 12 = -1/12 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-- Jet-space rank at mult-add, coordinates (r₁, r₂, r₃), scaled by 12:
    zoom (12, 0, 0), weight (12, 0, 0), bend (0, 12, 0), rung (12, 3, −1).
    The weight is the zoom (no new direction); zoom, bend, rung have determinant −144 ≠ 0. -/
def det3 (a1 a2 a3 b1 b2 b3 c1 c2 c3 : Int) : Int :=
  a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)

theorem jet_rank :
    det3 12 0 0 0 12 0 12 3 (-1) = -144 ∧ det3 12 0 0 12 0 0 0 12 0 = 0 := by decide

/-! ## D. The nome: a gluing the law never sees -/

/-- Gluing a ↦ a ◯ (q − 1) commutes with the law, for every q. -/
theorem glue_commutes (q a b : Int) : (a ◯ (q - 1)) ◯ b = (a ◯ b) ◯ (q - 1) := by
  unfold circ; grind

/-- Two gluings are one gluing by the product: the gluings form the group q^ℤ (in mult-add form). -/
theorem glue_twice (q q' a : Int) : (a ◯ (q - 1)) ◯ (q' - 1) = a ◯ (q * q' - 1) := by
  unfold circ; grind

/-- The law descends to the quotient: glued inputs give a glued output, whatever the gluings.
    So the torus ℂ*/q^ℤ carries mult-add for every nome, and the nome is invisible near 0. -/
theorem glue_descends (x y a b : Int) : (a ◯ (x - 1)) ◯ (b ◯ (y - 1)) = (a ◯ b) ◯ (x * y - 1) := by
  unfold circ; grind

/-- The 3-division polynomial of y² + a₁xy = x³ + a₂x² (a₃ = a₄ = a₆ = 0):
    ψ₃ = 3x⁴ + b₂x³ with b₂ = a₁² + 4a₂ (b₄ = b₆ = b₈ = 0). -/
def psi3 (a1 a2 x : Rat) : Rat := 3 * x ^ 4 + (a1 * a1 + 4 * a2) * x ^ 3

/-- Cusp y² = x³ (addition): ψ₃ = 3x⁴, whose only root is the cusp itself: no point other than
    the start comes back in three steps. -/
theorem cusp_psi3 (x : Rat) : psi3 0 0 x = 3 * x ^ 4 := by unfold psi3; grind

/-- Node y² + xy = x³ (mult-add): ψ₃ = x³(3x + 1). Besides the node (x = 0) the root x = −1/3. -/
theorem node_psi3 (x : Rat) : psi3 1 0 x = x ^ 3 * (3 * x + 1) := by unfold psi3; grind

/-- At x = −1/3 the node's equation y² + xy − x³ = 0 reads y² − y/3 + 1/27 = 0, whose
    discriminant 1/9 − 4/27 = −1/27 is not 0: two points. With the start, three: mult-add's
    three cube roots of 1. (The Tate torus has nine: see checks.py, 5a.) -/
theorem node_three_points : (-1/3 : Rat) * (-1/3) - 4 * (-(-1/3 : Rat) ^ 3) = -1/27 ∧ (-1/27 : Rat) ≠ 0 := by
  constructor <;> decide +kernel

/-! ## E. Several waves -/

/-- Two waves with pair weight A: τ = 1 + w₁ + w₂ + A·w₁w₂. -/
def tau2 (A w1 w2 : Int) : Int := 1 + w1 + w2 + A * w1 * w2

/-- No separate rescaling of the waves changes the pair weight: if τ_A(λ₁w₁, λ₂w₂) = τ_{A'}(w₁, w₂)
    for all w₁, w₂, then λ₁ = λ₂ = 1 and A = A'. (A single law's weight is units; a pair weight is not.) -/
theorem pair_weight_survives (A A' l1 l2 : Int)
    (h : ∀ w1 w2 : Int, tau2 A (l1 * w1) (l2 * w2) = tau2 A' w1 w2) : A = A' := by
  have h10 := h 1 0
  have h01 := h 0 1
  have h11 := h 1 1
  unfold tau2 at h10 h01 h11
  have e1 : l1 = 1 := by grind
  have e2 : l2 = 1 := by grind
  subst e1; subst e2
  grind

/-- Period matrices for g waves have g(g+1)/2 entries; curves of genus g have 3g − 3 moduli.
    The first exceeds the second exactly from g = 4 on (checked for 2 ≤ g < 40): Schottky's problem. -/
theorem periods_outnumber_curves :
    ∀ g, g < 40 → 2 ≤ g → (6 * g - 6 < g * (g + 1) ↔ 4 ≤ g) := by decide

/-! ## F. Whole numbers: the weight keeps its primes -/

theorem pseries_two (k x : Int) : wpow k x 2 = 2 * x + k * x * x := by
  rw [show wpow k x 2 = circK k (circK k 0 x) x from rfl]; unfold circK; grind

theorem pseries_three (k x : Int) : wpow k x 3 = 3 * x + 3 * k * x * x + k * k * x * x * x := by
  rw [show wpow k x 3 = circK k (circK k (circK k 0 x) x) x from rfl]; unfold circK; grind

theorem pseries_five (k x : Int) :
    wpow k x 5 = 5 * x + 10 * k * x * x + 10 * k * k * x * x * x + 5 * k * k * k * x * x * x * x
      + k * k * k * k * x * x * x * x * x := by
  rw [show wpow k x 5 = circK k (circK k (circK k (circK k (circK k 0 x) x) x) x) x from rfl]
  unfold circK; grind

/-- p ∣ k: the p-fold power is p times a whole number (the law is additive at p). -/
theorem additive_at_two (m x : Int) : wpow (2 * m) x 2 = 2 * (x + m * x * x) := by
  rw [pseries_two]; grind
theorem additive_at_three (m x : Int) :
    wpow (3 * m) x 3 = 3 * (x + 3 * m * x * x + 3 * m * m * x * x * x) := by
  rw [pseries_three]; grind
theorem additive_at_five (m x : Int) :
    wpow (5 * m) x 5 = 5 * (x + 10 * m * x * x + 50 * m * m * x * x * x
      + 125 * m * m * m * x * x * x * x + 125 * m * m * m * m * x * x * x * x * x) := by
  rw [pseries_five]; grind

/-- p ∤ k: the p-fold power is k^{p−1}x^p plus p times a whole number, and k^{p−1} is 1 more than
    a multiple of p (shown for p = 3; the law is multiplicative at p, height 1). -/
theorem mult_at_three (k x : Int) : wpow k x 3 - k * k * x * x * x = 3 * (x + k * x * x) := by
  rw [pseries_three]; grind
theorem fermat_three_a (m : Int) : (3 * m + 1) * (3 * m + 1) - 1 = 3 * (3 * m * m + 2 * m) := by grind
theorem fermat_three_b (m : Int) : (3 * m + 2) * (3 * m + 2) - 1 = 3 * (3 * m * m + 4 * m + 1) := by grind
theorem mult_at_two (k x : Int) : wpow k x 2 - k * x * x = 2 * x := by
  rw [pseries_two]; grind

end ShapeDirections
