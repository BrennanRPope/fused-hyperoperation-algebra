/-
  Quantum action  (DRAFT, Nov 26 2026)  -- first compiled in the shared library on Nov 26;
  `mul_one` and `one_mul` needed a closing `grind` after `simp`.
  ======================================================================================

  Companion: ../README.md and ../checks.py. Written for core Lean 4.37.0-rc1 (no mathlib), no `sorry`,
  in the style of Topics/Precession.lean (rationals, `ext`, `simp`, `grind`). Every identity below was
  also checked symbolically with sympy (checks.py section B1 and the script lean_identities_check.py).

    1. Phase-space moves (the Heisenberg group over the rationals). Composing two moves adds them and
       adds half their cross product to the phase. Swapping the order changes the phase by exactly
       the parallelogram area q₁p₂ − p₁q₂: the swap rule in phase space.
    2. A closed triangle of moves returns to the start with phase = its area, and a closed
       quadrilateral with phase = its shoelace area: action is enclosed area.
    3. One dial for Fermat, Newton and Einstein: the boost (t, x) ↦ (t + d·u·x, x − u·t) scales the
       form t² + d·x² by (1 + d·u²), and two boosts compose by the plane law w = (u + v)/(1 − d·u·v).
    4. The k part of the full law (x + y + kxy)/(1 − dxy) is only a change of ruler:
       x ↦ x/(1 + (k/2)x) turns it into the k = 0 law with d' = d − k²/4 (cleared of denominators).
    5. Regge in 2D: with equilateral triangles the deficit at a vertex is (6 − degree)·π/3, the sum of
       (6 − degree) is 6·(V − E + F) whenever 2E = 3F, and an edge flip leaves it unchanged.
-/
import FusedAlgebra.Base

namespace QuantumAction

/-! ## 1. Moves in phase space -/

/-- A move in phase space: shift position by `q`, momentum by `p`, and carry a phase `s`
    (the phase is measured in units of ħ, so `s` is an action). -/
structure Move where
  q : Rat
  p : Rat
  s : Rat

@[ext] theorem Move.ext' {a b : Move} (h0 : a.q = b.q) (h1 : a.p = b.p) (h2 : a.s = b.s) : a = b := by
  cases a; cases b; simp_all

/-- Composing moves: positions and momenta add, the phase picks up half the cross product. -/
def Move.mul (a b : Move) : Move :=
  ⟨a.q + b.q, a.p + b.p, a.s + b.s + (a.q * b.p - a.p * b.q) / 2⟩

def Move.one : Move := ⟨0, 0, 0⟩
def Move.inv (a : Move) : Move := ⟨-a.q, -a.p, -a.s⟩
/-- A pure move with no phase. -/
def mv (q p : Rat) : Move := ⟨q, p, 0⟩
/-- The parallelogram area spanned by two moves. -/
def cross (a b : Move) : Rat := a.q * b.p - a.p * b.q

theorem mul_assoc (a b c : Move) : Move.mul (Move.mul a b) c = Move.mul a (Move.mul b c) := by
  ext <;> simp [Move.mul] <;> grind

theorem mul_one (a : Move) : Move.mul a Move.one = a := by
  ext <;> simp [Move.mul, Move.one] <;> grind

theorem one_mul (a : Move) : Move.mul Move.one a = a := by
  ext <;> simp [Move.mul, Move.one] <;> grind

theorem mul_inv (a : Move) : Move.mul a (Move.inv a) = Move.one := by
  ext <;> simp [Move.mul, Move.inv, Move.one] <;> grind

/-- Swap rule in phase space: the two orders land in the same place... -/
theorem swap_same_place (a b : Move) :
    (Move.mul a b).q = (Move.mul b a).q ∧ (Move.mul a b).p = (Move.mul b a).p := by
  constructor <;> simp [Move.mul] <;> grind

/-- ...and differ in phase by exactly the area they span. -/
theorem swap_costs_area (a b : Move) : (Move.mul a b).s - (Move.mul b a).s = cross a b := by
  simp [Move.mul, cross]; grind

/-- One line, no swap cost: parallel moves commute. -/
theorem parallel_no_cost (q p r : Rat) : cross (mv q p) (mv (r * q) (r * p)) = 0 := by
  simp [cross, mv]; grind

/-! ## 2. Closed loops: phase = enclosed area -/

/-- Go along `a`, then `b`, then straight back: a closed triangle. -/
theorem triangle_loop (q₁ p₁ q₂ p₂ : Rat) :
    Move.mul (Move.mul (mv q₁ p₁) (mv q₂ p₂)) (mv (-(q₁ + q₂)) (-(p₁ + p₂)))
      = ⟨0, 0, (q₁ * p₂ - p₁ * q₂) / 2⟩ := by
  ext <;> simp [Move.mul, mv] <;> grind

/-- Shoelace area of a quadrilateral with corners 0, A, B, C. -/
def quadArea (a b c d e f : Rat) : Rat := ((a * d - b * c) + (c * f - d * e)) / 2

/-- Going round the quadrilateral 0 → (a,b) → (c,d) → (e,f) → 0 returns with phase = its area. -/
theorem quad_loop (a b c d e f : Rat) :
    Move.mul (Move.mul (Move.mul (mv a b) (mv (c - a) (d - b))) (mv (e - c) (f - d))) (mv (-e) (-f))
      = ⟨0, 0, quadArea a b c d e f⟩ := by
  ext <;> simp [Move.mul, mv, quadArea] <;> grind

/-- Reversing a loop reverses its phase. -/
theorem triangle_loop_reversed (q₁ p₁ q₂ p₂ : Rat) :
    Move.mul (Move.mul (mv q₂ p₂) (mv q₁ p₁)) (mv (-(q₁ + q₂)) (-(p₁ + p₂)))
      = ⟨0, 0, -((q₁ * p₂ - p₁ * q₂) / 2)⟩ := by
  ext <;> simp [Move.mul, mv] <;> grind

/-! ## 3. One dial: Fermat (d > 0), Newton (d = 0), Einstein (d < 0) -/

/-- The unnormalised boost with velocity `u` in the plane with dial `d`, as a 2×2 matrix
    acting on (t, x): t' = t + d·u·x, x' = −u·t + x. -/
structure M2 where
  a : Rat
  b : Rat
  c : Rat
  e : Rat

@[ext] theorem M2.ext' {m n : M2} (h0 : m.a = n.a) (h1 : m.b = n.b) (h2 : m.c = n.c) (h3 : m.e = n.e) :
    m = n := by
  cases m; cases n; simp_all

def M2.mul (m n : M2) : M2 := ⟨m.a * n.a + m.b * n.c, m.a * n.b + m.b * n.e, m.c * n.a + m.e * n.c, m.c * n.b + m.e * n.e⟩
def M2.smul (r : Rat) (m : M2) : M2 := ⟨r * m.a, r * m.b, r * m.c, r * m.e⟩
def boost (d u : Rat) : M2 := ⟨1, d * u, -u, 1⟩

/-- The boost keeps the plane's length form t² + d·x², up to the factor 1 + d·u². -/
theorem boost_keeps_form (d u t x : Rat) :
    (t + d * u * x) * (t + d * u * x) + d * ((x - u * t) * (x - u * t))
      = (1 + d * u * u) * (t * t + d * (x * x)) := by
  grind

/-- Two boosts make one boost, with velocities composed by the plane law w(1 − d·u·v) = u + v. -/
theorem boost_compose (d u v w : Rat) (hw : w * (1 - d * u * v) = u + v) :
    M2.mul (boost d u) (boost d v) = M2.smul (1 - d * u * v) (boost d w) := by
  ext <;> simp [M2.mul, M2.smul, boost] <;> grind

/-- Galilean end of the dial: at d = 0 velocities simply add. -/
theorem galilean_compose (u v : Rat) :
    M2.mul (boost 0 u) (boost 0 v) = boost 0 (u + v) := by
  ext <;> simp [M2.mul, boost] <;> grind

/-- The length element √(1 + d·v²) expands as 1 + d·v²/2 + O(d²): the Newtonian term is ½v² for
    every d. Stated as the exact identity (1 + d v²/2)² = 1 + d v² + (d v²/2)². -/
theorem newton_from_dial (d v : Rat) :
    (1 + d * v * v / 2) * (1 + d * v * v / 2) = 1 + d * v * v + (d * v * v / 2) * (d * v * v / 2) := by
  grind

/-! ## 4. The k part of the full law is a change of ruler (d' = d − k²/4) -/

/-- Cleared form of: f(x ∘ y) = f(x) ⊕_{d'} f(y), with x ∘ y = (x + y + kxy)/(1 − dxy),
    f(x) = x/(1 + a·x), a = k/2, d' = d − k²/4, and ⊕_{d'} the k = 0 plane law.
    Left: numerator of f(N/D) is N, denominator D + aN.  Right: numerator x(1+ay) + y(1+ax),
    denominator (1+ax)(1+ay) − d'·x·y.  Cross-multiplied they agree. -/
theorem k_is_a_slide (k d x y : Rat) :
    (x + y + k * x * y) * ((1 + k / 2 * x) * (1 + k / 2 * y) - (d - k * k / 4) * x * y)
      = (x * (1 + k / 2 * y) + y * (1 + k / 2 * x)) * ((1 - d * x * y) + k / 2 * (x + y + k * x * y)) := by
  grind

/-! ## 5. Regge in two dimensions: the action does not move -/

/-- With equilateral triangles every vertex's deficit is (6 − degree)·π/3. Summed over the vertices,
    Σ(6 − deg) = 6V − 2E (each edge has two ends), and when every face is a triangle (2E = 3F)
    this is 6·(V − E + F): six times the Euler number, whatever the triangulation. -/
theorem regge2d_total (V E F : Int) (htri : 2 * E = 3 * F) : 6 * V - 2 * E = 6 * (V - E + F) := by
  grind

/-- An edge flip lowers the degrees of the edge's two ends by one and raises the two opposite
    corners by one, so the total deficit (6 − d₁) + (6 − d₂) + (6 − d₃) + (6 − d₄) is unchanged. -/
theorem regge2d_flip (d₁ d₂ d₃ d₄ : Int) :
    (6 - (d₁ - 1)) + (6 - (d₂ - 1)) + (6 - (d₃ + 1)) + (6 - (d₄ + 1))
      = (6 - d₁) + (6 - d₂) + (6 - d₃) + (6 - d₄) := by
  grind

/-- On a sphere (V − E + F = 2) the deficits add to 12 units of π/3, i.e. 4π (checks E1a, E1c). -/
theorem regge2d_sphere (V E F : Int) (htri : 2 * E = 3 * F) (hχ : V - E + F = 2) :
    6 * V - 2 * E = 12 := by
  grind

end QuantumAction
