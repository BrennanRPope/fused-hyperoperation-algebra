/-
  Fused Hyperoperation Algebra: elliptic counting.
  Core Lean 4 only (no mathlib), no `sorry`. NOT YET COMPILED where it was written (no Lean in
  that workspace). Every statement was also checked numerically or exactly in
  ../elliptic_counting.py (14 of 14 checks pass).

  To compile: copy to lean-library/FusedAlgebra/Topics/EllipticCounting.lean, add it after
  FusedAlgebra.Topics.Patterns in build.sh and FusedAlgebra.lean, run ./build.sh.

  The story in one line: a landscape law counts with a Lucas sequence U(n); an elliptic law
  counts with an elliptic divisibility sequence W(n) (Ward 1948). The Lucas sequence IS an
  elliptic divisibility sequence, up to the factor Q^(n−1) (removed by a gauge): the landscape
  is the pinched edge of the elliptic counting.

  Known mathematics throughout: Ward (1948), Somos (1989), Tate normal form, Mazur (1977).
  LIBRARY NOTE: first compiled in the shared library, after Patterns. Changes: lucas_shift_one's
  step case states the two recurrence steps as hypotheses for grind; lucas_eds_unit uses
  MultAdd.one_pow'; disc gets a `;` after `let b8 := …`, since the next line's leading minus was
  being read as part of b8; order_five_point's last conjunct read w3 * (w3 * w3 * w3), which is
  false (it is −2); it now reads w3 * w3 * w3, the W3³ its docstring names. No exact twins.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Patterns

namespace EllipticCounting
open Patterns

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ### 1. Lucas sequences are elliptic divisibility sequences (up to Q^(n−1)) -/

/-- one-step shift: U(a+1)U(a+e) − U(a)U(a+e+1) = Q^a · U(e) -/
theorem lucas_shift_one (P Q : R) (e a : Nat) :
    lucas P Q (a + 1) * lucas P Q (a + e) - lucas P Q a * lucas P Q (a + e + 1)
      = Q ^ a * lucas P Q e := by
  induction a with
  | zero =>
    rw [Nat.zero_add, Nat.zero_add, lucas_one, lucas_zero, Lean.Grind.Semiring.pow_zero]; grind
  | succ a ih =>
    have h1 := lucas_succ_succ P Q a
    have h2 := lucas_succ_succ P Q (a + e)
    rw [show a + 1 + e = a + e + 1 by omega, show a + e + 1 + 1 = a + e + 2 by omega,
      show a + 1 + 1 = a + 2 by omega, Lean.Grind.Semiring.pow_succ]
    grind

/-- the general exchange identity: U(a+c)U(a+e) − U(a)U(a+e+c) = Q^a · U(c) · U(e) -/
theorem lucas_exchange_pair (P Q : R) (a e c : Nat) :
    (lucas P Q (a + c) * lucas P Q (a + e) - lucas P Q a * lucas P Q (a + e + c)
        = Q ^ a * lucas P Q c * lucas P Q e) ∧
    (lucas P Q (a + (c + 1)) * lucas P Q (a + e) - lucas P Q a * lucas P Q (a + e + (c + 1))
        = Q ^ a * lucas P Q (c + 1) * lucas P Q e) := by
  induction c with
  | zero =>
    constructor
    · rw [Nat.add_zero, Nat.add_zero, lucas_zero]; grind
    · rw [Nat.zero_add, lucas_one, lucas_shift_one]; grind
  | succ c ih =>
    obtain ⟨h0, h1⟩ := ih
    refine ⟨h1, ?_⟩
    rw [show a + (c + 1 + 1) = (a + c) + 2 by omega, lucas_succ_succ P Q (a + c),
      show a + e + (c + 1 + 1) = (a + e + c) + 2 by omega, lucas_succ_succ P Q (a + e + c),
      show c + 1 + 1 = c + 2 by omega, lucas_succ_succ P Q c,
      show a + c + 1 = a + (c + 1) by omega, show a + e + c + 1 = a + e + (c + 1) by omega]
    grind

theorem lucas_exchange (P Q : R) (a e c : Nat) :
    lucas P Q (a + c) * lucas P Q (a + e) - lucas P Q a * lucas P Q (a + e + c)
      = Q ^ a * lucas P Q c * lucas P Q e :=
  (lucas_exchange_pair P Q a e c).1

/-- Catalan's identity in the form used below: U(n+d)² − U(d)·U(2n+d) = Q^d · U(n)² -/
theorem lucas_catalan (P Q : R) (n d : Nat) :
    lucas P Q (d + n) * lucas P Q (d + n) - lucas P Q d * lucas P Q (d + n + n)
      = Q ^ d * lucas P Q n * lucas P Q n :=
  lucas_exchange P Q d n n

/-- **Lucas counting is elliptic counting, weighted by Q^(n−1).**
    Writing n = f + 1 and m = f + e + 2 (so m > n ≥ 1):
    Q^(n−1)·U(m+n)·U(m−n) = U(m+1)U(m−1)U(n)² − U(n+1)U(n−1)U(m)².
    This is the elliptic divisibility recurrence (with U(1) = 1), with the extra factor Q^(n−1). -/
theorem lucas_eds (P Q : R) (f e : Nat) :
    Q ^ f * (lucas P Q (2 * f + e + 3) * lucas P Q (e + 1))
      = lucas P Q (f + e + 3) * lucas P Q (f + e + 1) * (lucas P Q (f + 1) * lucas P Q (f + 1))
        - lucas P Q (f + 2) * lucas P Q f * (lucas P Q (f + e + 2) * lucas P Q (f + e + 2)) := by
  have hM := cassini P Q (f + e + 1)
  have hN := cassini P Q f
  have hC := lucas_catalan P Q (f + 1) (e + 1)
  rw [show f + e + 1 + 1 = f + e + 2 by omega, show f + e + 1 + 2 = f + e + 3 by omega] at hM
  rw [show f + 1 + 1 = f + 2 by omega] at hN
  rw [show e + 1 + (f + 1) = f + e + 2 by omega,
    show f + e + 2 + (f + 1) = 2 * f + e + 3 by omega] at hC
  have hp : Q ^ (f + e + 1) = Q ^ f * Q ^ (e + 1) := by
    rw [show f + e + 1 = f + (e + 1) by omega, MultAdd.pow_add']
  grind

/-- with Q = 1 the Lucas sequence satisfies the elliptic divisibility recurrence exactly -/
theorem lucas_eds_unit (P : R) (f e : Nat) :
    lucas P 1 (2 * f + e + 3) * lucas P 1 (e + 1)
      = lucas P 1 (f + e + 3) * lucas P 1 (f + e + 1) * (lucas P 1 (f + 1) * lucas P 1 (f + 1))
        - lucas P 1 (f + 2) * lucas P 1 f * (lucas P 1 (f + e + 2) * lucas P 1 (f + e + 2)) := by
  have h := lucas_eds P 1 f e
  rw [MultAdd.one_pow' f] at h
  grind

/-- the Somos-4 edge (n = 2): Q·U(m+4)U(m) = P²·U(m+3)U(m+1) − (P² − Q)·U(m+2)² -/
theorem lucas_somos4 (P Q : R) (m : Nat) :
    Q * (lucas P Q (m + 4) * lucas P Q m)
      = P * P * (lucas P Q (m + 3) * lucas P Q (m + 1))
        - (P * P - Q) * (lucas P Q (m + 2) * lucas P Q (m + 2)) := by
  cases m with
  | zero =>
    simp only [lucas, Nat.zero_add]; grind
  | succ e =>
    have h := lucas_eds P Q 1 e
    have h2 : lucas P Q 2 = P := by simp only [lucas]; grind
    have h3 : lucas P Q 3 = P * P - Q := by simp only [lucas]; grind
    rw [show 2 * 1 + e + 3 = e + 1 + 4 by omega, show 1 + e + 3 = e + 1 + 3 by omega,
      show 1 + e + 1 = e + 1 + 1 by omega, show 1 + e + 2 = e + 1 + 2 by omega,
      show (1 : Nat) + 1 = 2 by rfl, show (1 : Nat) + 2 = 3 by rfl, h2, h3, lucas_one,
      Lean.Grind.Semiring.pow_succ, Lean.Grind.Semiring.pow_zero] at h
    grind

/-! ### 2. The family with a point of order 5, and its golden pinches -/

/-- the standard b-invariants and discriminant of y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6 -/
def disc (a1 a2 a3 a4 a6 : R) : R :=
  let b2 := a1 * a1 + 4 * a2
  let b4 := 2 * a4 + a1 * a3
  let b6 := a3 * a3 + 4 * a6
  let b8 := a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
  -(b2 * b2 * b8) - 8 * (b4 * b4 * b4) - 27 * (b6 * b6) + 9 * b2 * b4 * b6

/-- Tate normal form with a point of order 5: y² + (1−b)xy − by = x³ − bx².
    Its discriminant is b⁵(b² − 11b − 1), so it pinches only at b = 0 and b² = 11b + 1. -/
theorem tate5_disc (b : R) :
    disc (1 - b) (-b) (-b) 0 0 = b * b * b * b * b * (b * b - 11 * b - 1) := by
  simp only [disc]; grind

/-- the pinch points b² − 11b − 1 = 0 are φ⁵ = 5φ + 3 and −φ⁻⁵ = 8 − 5φ, for any φ with φ² = φ + 1 -/
theorem tate5_golden_pinch (φ : R) (h : φ * φ = φ + 1) :
    (5 * φ + 3) * (5 * φ + 3) - 11 * (5 * φ + 3) - 1 = 0 ∧
    (8 - 5 * φ) * (8 - 5 * φ) - 11 * (8 - 5 * φ) - 1 = 0 := by
  constructor <;> grind

/-- and 5φ + 3 really is φ⁵ -/
theorem golden_fifth_power (φ : R) (h : φ * φ = φ + 1) : φ ^ 5 = 5 * φ + 3 := by
  simp only [Lean.Grind.Semiring.pow_succ, Lean.Grind.Semiring.pow_zero]; grind

end Ring

/-! ### 3. A concrete elliptic counting sequence: y² + y = x³ − x at P = (0, 0) -/

/-- division-polynomial seeds at a point with x = 0: W2 = 2y + a1 x + a3, W3 = b8, W4 = W2·(b4 b8 − b6²) -/
def seedW2 (a3 y : Int) : Int := 2 * y + a3
def seedW3 (a1 a2 a3 a4 a6 : Int) : Int :=
  a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
def seedW4 (a1 a2 a3 a4 a6 y : Int) : Int :=
  let b4 := 2 * a4 + a1 * a3
  let b6 := a3 * a3 + 4 * a6
  seedW2 a3 y * (b4 * seedW3 a1 a2 a3 a4 a6 - b6 * b6)

/-- the counting sequence of y² + y = x³ − x, W(0) .. W(25) -/
def w37List : List Int :=
  [0, 1, 1, -1, 1, 2, -1, -3, -5, 7, -4, -23, 29, 59, 129, -314, -65, 1529, -3689, -8209,
   -16264, 83313, 113689, -620297, 2382785, 7869898]

def w37 (n : Nat) : Int := w37List.getD n 0

/-- the seeds are the curve's division polynomials at (0, 0): a = (0, 0, 1, −1, 0) -/
theorem w37_seeds :
    w37 1 = 1 ∧ w37 2 = seedW2 1 0 ∧ w37 3 = seedW3 0 0 1 (-1) 0 ∧ w37 4 = seedW4 0 0 1 (-1) 0 0 := by
  decide

/-- (0, 0) is on y² + y = x³ − x -/
theorem w37_point : (0 : Int) * 0 + 0 = 0 * 0 * 0 - 0 := by decide

/-- Somos-4: W(n+4)W(n) = W(n+3)W(n+1) + W(n+2)² -/
theorem w37_somos4 : ∀ n, n < 22 → w37 (n + 4) * w37 n = w37 (n + 3) * w37 (n + 1) + w37 (n + 2) * w37 (n + 2) := by
  decide

/-- the elliptic divisibility recurrence, for 1 ≤ n < m, m + n ≤ 24 -/
theorem w37_eds : ∀ m, m < 13 → ∀ n, n < m → 1 ≤ n →
    w37 (m + n) * w37 (m - n) = w37 (m + 1) * w37 (m - 1) * (w37 n * w37 n)
      - w37 (n + 1) * w37 (n - 1) * (w37 m * w37 m) := by
  decide

/-- strong divisibility, like Fibonacci: gcd(W m, W n) = |W gcd(m, n)| -/
theorem w37_strong_div : ∀ m, m < 21 → ∀ n, n < 21 →
    Int.gcd (w37 m) (w37 n) = (w37 (Nat.gcd m n)).natAbs := by
  decide

/-- its odd terms are the classic Somos-4 numbers 1, 1, 2, 3, 7, 23, 59, 314, 1529, 8209 (up to sign) -/
theorem w37_odd_somos :
    [1, 3, 5, 7, 9, 11, 13, 15, 17, 19].map (fun n => (w37 n).natAbs)
      = [1, 1, 2, 3, 7, 23, 59, 314, 1529, 8209] := by
  decide

/-! ### 4. A rational point of order 5, which no line law can have -/

/-- y² + y = x³ − x², a = (0, −1, 1, 0, 0), at (0, 0): W2 = 1, W3 = −1, W4 = −1, W5 = W4·W2³ − W3³ = 0 -/
theorem order_five_point :
    let w2 := seedW2 1 0
    let w3 := seedW3 0 (-1) 1 0 0
    let w4 := seedW4 0 (-1) 1 0 0 0
    w2 = 1 ∧ w3 = -1 ∧ w4 = -1 ∧ w4 * (w2 * w2 * w2) - w3 * w3 * w3 = 0 := by
  decide

/-- its multiples, all on the curve: P = (0,0), 2P = (1,−1), 3P = (1,0), 4P = (0,−1),
    and 3P = −2P, 4P = −P (negation is y ↦ −y − 1), so 5P = O -/
theorem order_five_points :
    (∀ p ∈ [((0 : Int), (0 : Int)), (1, -1), (1, 0), (0, -1)],
        p.2 * p.2 + p.2 = p.1 * p.1 * p.1 - p.1 * p.1) ∧
    (1 : Int) = 1 ∧ (0 : Int) = -(-1) - 1 ∧ (-1 : Int) = -0 - 1 := by
  decide

/-- why: the tangent at P is y = 0, cutting x³ − x² = x²(x − 1), so −2P = (1, 0);
    the chord through P and 2P is y = −x, cutting x(x − 1)², so P + 2P = −(1, −1) = (1, 0) -/
theorem order_five_lines (x : Int) :
    (0 * 0 + 0) - (x * x * x - x * x) = -(x * x * (x - 1)) ∧
    ((-x) * (-x) + (-x)) - (x * x * x - x * x) = -(x * (x - 1) * (x - 1)) := by
  constructor <;> grind

end EllipticCounting
