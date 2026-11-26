/-
  Cut-independence from the swap rule's q weights  (DRAFT, Nov 26 2026). First compiled in the
  shared library on Nov 26, as FusedAlgebra.Topics.CutIndependence. Every statement was also checked in Python
  (../lean_identities_check.py), and the numbers behind them in ../checks.py.

  Companion: ../README.md. Written for core Lean 4 (no mathlib), no `sorry`, in the style of the
  library's other topic files (`grind` for ring identities, `decide` for small integer facts).
  The q-binomial theorem and the U_q(sl2) commutator are already in the library
  (FusedAlgebra.Topics.QAcross.QCounting: `q_binomial_theorem`, `uq_sl2_commutator`); they are not
  repeated here.

    1. The swap rule y·x = q·x·y makes monomials multiply as (xᵃyᵇ)(xᶜyᵈ) = q^(b·c) x^(a+c) y^(b+d).
       The q-exponents agree for either bracketing of three monomials, so the plane is associative.
    2. A symmetry acting by E y = x, F x = y, K x = p·x, extended by the twisted product rule,
       leaves (p − q)·x² behind on y·x − q·x·y (and (p − q)·y² for F). It respects the swap rule
       exactly when p = q: the swap rule fixes the symmetry's parameter.
    3. The dial's own sequence u(n+2) = k·u(n+1) − u(n) (d = 1) obeys Cassini's identity
       u(n+1)² − u(n)·u(n+2) = 1 for every n: the step matrix [[k, −1], [1, 0]] has determinant 1.
    4. The same recursion with k = [2] = q + 1/q generates the quantum dimensions [n] directly.
       At the two roots of unity where [2] is an integer: [2] = 1 (q = e^(iπ/3)) and
       [2] = 0 (q = e^(iπ/2)), the quantum dimensions [n] first vanish at n = 3 and n = 2, the vertex
       weight w² = Σ [2j+1]² is 2 and 1, and the step matrix closes: M³ = −1 and M² = −1.
    5. d − k²/4 changes sign at k² = 4d: 4 − k² = (2 − k)(2 + k) with d = 1.
-/
import FusedAlgebra.Base

namespace CutIndependence

/-! ## 1. The swap rule's plane is associative -/

/-- The power of q picked up when (xᵃyᵇ) is multiplied by (xᶜyᵈ): each of the b copies of y
    passes each of the c copies of x once. -/
def swapExp (b c : Int) : Int := b * c

/-- ((xᵃyᵇ)(xᶜyᵈ))(xᵉyᶠ) and (xᵃyᵇ)((xᶜyᵈ)(xᵉyᶠ)) pick up the same power of q. -/
theorem swap_assoc (b c d e : Int) :
    swapExp b c + swapExp (b + d) e = swapExp d e + swapExp b (c + e) := by
  unfold swapExp; grind

/-! ## 2. The symmetry must carry the swap rule's q -/

/-- Coefficient of x² in E(y·x) − q·E(x·y), with E(y·x) = E(y)K(x) + y E(x) = p·x² and
    E(x·y) = E(x)K(y) + x E(y) = x². -/
def leftoverE (p q : Rat) : Rat := p * 1 - q * 1

/-- Coefficient of y² in F(y·x) − q·F(x·y), with F(y·x) = F(y)x + K⁻¹(y)F(x) = p·y² and
    F(x·y) = F(x)y + K⁻¹(x)F(y) = y². -/
def leftoverF (p q : Rat) : Rat := p * 1 - q * 1

theorem leftoverE_zero_iff (p q : Rat) : leftoverE p q = 0 ↔ p = q := by
  unfold leftoverE; constructor <;> intro h <;> grind

theorem leftoverF_zero_iff (p q : Rat) : leftoverF p q = 0 ↔ p = q := by
  unfold leftoverF; constructor <;> intro h <;> grind

/-! ## 3. The dial's sequence and Cassini's identity -/

/-- The dial's sequence with d = 1: u 0 = 0, u 1 = 1, u (n+2) = k·u (n+1) − u n. -/
def u (k : Int) : Nat → Int
  | 0 => 0
  | 1 => 1
  | n + 2 => k * u k (n + 1) - u k n

theorem u_step (k : Int) (n : Nat) : u k (n + 2) = k * u k (n + 1) - u k n := rfl

theorem u_two (k : Int) : u k 2 = k := by simp [u]

/-- Cassini: u(n+1)² − u(n)·u(n+2) = 1 = dⁿ, the determinant of the step matrix to the n. -/
theorem cassini (k : Int) (n : Nat) : u k (n + 1) * u k (n + 1) - u k n * u k (n + 2) = 1 := by
  induction n with
  | zero => simp [u]
  | succ n ih =>
    have h1 := u_step k n
    have h2 := u_step k (n + 1)
    grind

/-! ## 4. Roots of unity with an integer [2]  (here u is run with k = [2] = q + 1/q) -/

/-- q = e^(iπ/3): [2] = q + 1/q = 1, the dimensions run 0, 1, 1, 0 and first vanish at n = 3. -/
theorem dims_r3 : u 1 1 = 1 ∧ u 1 2 = 1 ∧ u 1 3 = 0 := by decide

/-- q = e^(iπ/2): [2] = 0, the dimensions run 0, 1, 0 and first vanish at n = 2. -/
theorem dims_r2 : u 0 1 = 1 ∧ u 0 2 = 0 := by decide

/-- The vertex weight w² = Σ_{2j+1 < r} [2j+1]²: 2 at r = 3 and 1 at r = 2, both equal to
    r / (2 sin²(π/r)). -/
theorem w2_r3 : u 1 1 * u 1 1 + u 1 2 * u 1 2 = 2 := by decide
theorem w2_r2 : u 0 1 * u 0 1 = 1 := by decide

/-- 2 × 2 integer matrices as quadruples (a, b, c, d) = [[a, b], [c, d]]. -/
def mm (m n : Int × Int × Int × Int) : Int × Int × Int × Int :=
  let (a, b, c, d) := m; let (e, f, g, h) := n
  (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)

/-- The step matrix [[k, −1], [1, 0]]. -/
def step (k : Int) : Int × Int × Int × Int := (k, -1, 1, 0)

/-- It closes up at a root of unity: M³ = −1 for [2] = 1 and M² = −1 for [2] = 0. -/
theorem step_closes_r3 : mm (step 1) (mm (step 1) (step 1)) = (-1, 0, 0, -1) := by decide
theorem step_closes_r2 : mm (step 0) (step 0) = (-1, 0, 0, -1) := by decide

/-- At q = 1 ([2] = 2) the sequence is u n = n: it never vanishes, so nothing truncates. -/
theorem u_at_two (n : Nat) : u 2 n = n := by
  induction n using Nat.strongRecOn with
  | ind n ih =>
    match n with
    | 0 => rfl
    | 1 => rfl
    | m + 2 =>
      rw [u_step, ih (m + 1) (by omega), ih m (by omega)]
      grind

/-! ## 5. Where d − k²/4 changes sign -/

theorem disc_factor (k : Rat) : 4 - k * k = (2 - k) * (2 + k) := by grind

end CutIndependence
