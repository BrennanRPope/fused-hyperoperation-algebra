/-
  Fused Hyperoperation Algebra: e, i and π from one dial (topics/why-i/e-i-pi).
  Core Lean 4 only (no mathlib), no `sorry`, no imports.
  Library path: FusedAlgebra/Topics/EIPi.lean (compiled, 17 theorems, standard axioms only).

  Only exact integer identities live here. The limits (e as the chain of n tiny steps,
  π as the half-turn) are numerical and are checked in ../checks.py (A1–A11).
  1. the lift: 1 + k·(a ∘ₖ b) = (1 + k·a)(1 + k·b), for two and three steps;
  2. a dial is a direction: k·(a ∘ₖ b) = (k·a) ∘₁ (k·b);
  3. a discrete Euler on the Gaussian integers: q = −1 + i lifts to i, q ∘ q = −2 lifts to −1,
     four copies give 0, which lifts to 1 (check A10);
  4. local tomography as mult-add of parameter counts: (K_A − 1) ∘₁ (K_B − 1) = K_A·K_B − 1,
     true for complex counts (n²), false for real (n(n+1)/2) and quaternionic (n(2n−1))
     counts (checks B1–B4).
-/
namespace EIPi

/-! ## 1. The mult-add dial and its lift -/

/-- Mult-add at dial k: a ∘ₖ b = a + b + k·a·b. -/
def madd (k a b : Int) : Int := a + b + k * a * b

/-- The lift 1 + k·x turns the dial into ordinary multiplication. -/
theorem lift_two (k a b : Int) : 1 + k * madd k a b = (1 + k * a) * (1 + k * b) := by
  unfold madd; grind

theorem lift_three (k a b c : Int) :
    1 + k * madd k (madd k a b) c = (1 + k * a) * (1 + k * b) * (1 + k * c) := by
  unfold madd; grind

/-- Chaining is associative, so "n steps" has one meaning. -/
theorem madd_assoc (k a b c : Int) : madd k (madd k a b) c = madd k a (madd k b c) := by
  unfold madd; grind

theorem madd_zero (k a : Int) : madd k a 0 = a := by
  unfold madd; grind

/-! ## 2. Choosing a dial is choosing a direction for the inputs -/

theorem dial_is_direction (k a b : Int) : k * madd k a b = madd 1 (k * a) (k * b) := by
  unfold madd; grind

/-! ## 3. Discrete Euler on the Gaussian integers -/

/-- Gaussian integer re + im·i. -/
structure G where
  re : Int
  im : Int
deriving DecidableEq, Repr

def G.add (x y : G) : G := ⟨x.re + y.re, x.im + y.im⟩
def G.mul (x y : G) : G := ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re⟩
def G.one : G := ⟨1, 0⟩
def G.i : G := ⟨0, 1⟩

/-- Mult-add at dial 1 on Gaussian integers. -/
def G.madd (x y : G) : G := G.add (G.add x y) (G.mul x y)

/-- The lift at dial 1. -/
def G.lift (x : G) : G := G.add G.one x

def q : G := ⟨-1, 1⟩

theorem q_lifts_to_i : G.lift q = G.i := by decide

theorem i_squared : G.mul G.i G.i = ⟨-1, 0⟩ := by decide

/-- Two quarter turns: q ∘ q = −2, the mirror point −2/k at k = 1 ... -/
theorem two_quarter_turns : G.madd q q = ⟨-2, 0⟩ := by decide

/-- ... whose lift is −1: the discrete e^{iπ} = −1. -/
theorem half_turn_lifts_to_minus_one : G.lift (G.madd q q) = ⟨-1, 0⟩ := by decide

/-- Four quarter turns come back to the identity 0 (lift 1). -/
theorem four_quarter_turns : G.madd (G.madd q q) (G.madd q q) = ⟨0, 0⟩ := by decide

theorem full_turn_lifts_to_one : G.lift (G.madd (G.madd q q) (G.madd q q)) = G.one := by decide

/-- The lift turns Gaussian mult-add into Gaussian multiplication. -/
theorem G_lift_hom (x y : G) : G.lift (G.madd x y) = G.mul (G.lift x) (G.lift y) := by
  cases x; cases y; simp only [G.lift, G.madd, G.add, G.mul, G.one, G.mk.injEq]; grind

/-! ## 4. Local tomography is mult-add of parameter counts -/

/-- Hardy's K_AB = K_A·K_B, read as mult-add of the counts minus one. -/
theorem counts_mult_add (KA KB : Int) : KA * KB - 1 = madd 1 (KA - 1) (KB - 1) := by
  unfold madd; grind

/-- Complex counts K(n) = n² always satisfy it: (a²)(b²) − 1 = (a² − 1) ∘₁ (b² − 1). -/
theorem complex_counts (a b : Int) :
    a * a * (b * b) - 1 = (a * a - 1) + (b * b - 1) + (a * a - 1) * (b * b - 1) := by
  grind

def Kc (n : Nat) : Nat := n * n
def Kr (n : Nat) : Nat := n * (n + 1) / 2
def Kq (n : Nat) : Nat := n * (2 * n - 1)

/-- Two qubits: 16 = 4·4, and 15 = 3 + 3 + 9 (the Bloch / correlation split). -/
theorem two_qubits : Kc 4 = Kc 2 * Kc 2 ∧ Kc 4 - 1 = 3 + 3 + 3 * 3 := by decide

/-- Two real bits: 10 parameters, but local data gives only 3·3 = 9 (one hidden: σy⊗σy). -/
theorem two_rebits : Kr 4 = 10 ∧ Kr 2 * Kr 2 = 9 := by decide

/-- Two quaternionic bits: 28 parameters, but local products give 6·6 = 36 cells. -/
theorem two_quaterbits : Kq 4 = 28 ∧ Kq 2 * Kq 2 = 36 := by decide

end EIPi
