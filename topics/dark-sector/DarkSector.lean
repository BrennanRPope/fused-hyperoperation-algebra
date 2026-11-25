/-
  Fused Hyperoperation Algebra and the dark sector: the exact statements behind the page.
  Core Lean 4 only (no mathlib), no `sorry`.

  circ a b = a + b + a*b = (1+a)(1+b) - 1        (rung 2, mult-add)
  fexp a b = (1+a)^(1+b) - 1                      (rung 3, fused exponent)
-/
namespace DarkSector

def circ (a b : Nat) : Nat := a + b + a * b
def fexp (a b : Nat) : Nat := (a + 1) ^ (b + 1) - 1
/-- n-fold mult-add power of z:  z ⊕ z ⊕ … ⊕ z  (n copies), with 0 copies = 0. -/
def mpow (z : Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => circ (mpow z n) z

theorem circ_lift (a b : Nat) : circ a b + 1 = (a + 1) * (b + 1) := by unfold circ; grind

/-! ## 1. Redshifts compose by mult-add (exact, no model needed)
Light stretched by 1+z₁ and then by 1+z₂ is stretched by (1+z₁)(1+z₂).  So the total redshift is z₁ ⊕ z₂.
Geometrically: the (1+z₁)×(1+z₂) rectangle = unit square + two strips z₁, z₂ + corner z₁z₂. -/
theorem redshift_compose (z₁ z₂ : Nat) : (1 + z₁) * (1 + z₂) = 1 + circ z₁ z₂ := by unfold circ; grind
/-- n equal stretches: (1+z)^n = 1 + (z ⊕ z ⊕ … ⊕ z). -/
theorem redshift_power (z n : Nat) : (1 + z) ^ n = 1 + mpow z n := by
  induction n with
  | zero => simp [mpow]
  | succ n ih => rw [Nat.pow_succ, ih, mpow]; unfold circ; grind
/-- The same with integers, so blueshifts (negative z) are covered too. -/
theorem redshift_compose_int (z₁ z₂ : Int) : (1 + z₁) * (1 + z₂) = 1 + (z₁ + z₂ + z₁ * z₂) := by grind

/-! ## 2. Matter dilutes as the mult-add cube of redshift
(1+z)³ = 1 + 3z + 3z² + z³ : unit cube, three slabs, three bars, one small cube (Pascal layers 1,3,3,1). -/
theorem cube_layers (z : Nat) : (1 + z) ^ 3 = 1 + 3 * z + 3 * z ^ 2 + z ^ 3 := by grind
theorem cube_is_mpow3 (z : Nat) : mpow z 3 = 3 * z + 3 * z ^ 2 + z ^ 3 := by
  have h := redshift_power z 3; rw [cube_layers] at h; omega
/-- Flat ΛCDM, exactly:  Ωm(1+z)³ + (1-Ωm) = 1 + Ωm·(z ⊕ z ⊕ z).
The cosmological constant is the unit cube; matter is the weighted shell of the (1+z) cube around it. -/
theorem lcdm_is_weighted_cube (Om z : Rat) :
    Om * (1 + z) ^ 3 + (1 - Om) = 1 + Om * (3 * z + 3 * z ^ 2 + z ^ 3) := by grind

/-! ## 3. Constant expansion rate = growth fractions that compose by mult-add (rung 2)
If the fractional growth A(t) over a time t satisfies A(t+s) = A(t) ⊕ A(s), the scale factor 1+A grows
exponentially: 1 + A(t) = (1 + A(1))^t.  That is de Sitter expansion, the cosmological constant. -/
theorem circ_growth_is_exponential (A : Nat → Nat) (hA : ∀ t s, A (t + s) = circ (A t) (A s)) :
    ∀ t, A t + 1 = (A 1 + 1) ^ t := by
  have h0 : A 0 = 0 := by
    have := hA 0 0; simp [circ] at this; grind
  intro t
  induction t with
  | zero => simp [h0]
  | succ n ih => rw [hA n 1, circ_lift, ih, Nat.pow_succ]

/-! ## 4. No-go one rung up
Expansion cannot compose at the fused-exponent rung.  If growth over t+s were A(t) ⊛ A(s), in either order,
and nothing has grown at t = 0, then nothing ever grows.  So the ladder tops out at rung 2 for expansion:
the cosmological constant is the highest rung that can carry it. -/
theorem fexp_growth_trivial (A : Nat → Nat) (h0 : A 0 = 0) (hA : ∀ t s, A (t + s) = fexp (A t) (A s)) :
    ∀ t, A t = 0 := by
  intro t
  have := hA 0 t; rw [Nat.zero_add, h0] at this; unfold fexp at this; simpa using this
theorem fexp_growth_trivial' (A : Nat → Nat) (h0 : A 0 = 0) (hA : ∀ t s, A (t + s) = fexp (A s) (A t)) :
    ∀ t, A t = 0 := by
  intro t
  have := hA t 0; rw [Nat.add_zero, h0] at this; unfold fexp at this; simpa using this
/-- The same in the unshifted form: a(t+s) = a(t)^a(s) with a(0) = 1 forces a ≡ 1. -/
theorem pow_growth_trivial (a : Nat → Nat) (h0 : a 0 = 1) (ha : ∀ t s, a (t + s) = a t ^ a s) :
    ∀ t, a t = 1 := by
  intro t; have := ha 0 t; rw [Nat.zero_add, h0, Nat.one_pow] at this; exact this

/-! ## 5. Why mult-add of gravity points the wrong way
Combining non-negative pulls by mult-add never gives less than their sum, and the extra is exactly the
product: a corner rectangle a·b.  For many pieces the extra only grows. -/
theorem circ_excess (a b : Nat) : circ a b = (a + b) + a * b := rfl
theorem circ_ge_add (a b : Nat) : a + b ≤ circ a b := by unfold circ; omega
def mfold : List Nat → Nat
  | [] => 0
  | x :: xs => circ x (mfold xs)
theorem mfold_ge_sum : ∀ l : List Nat, l.sum ≤ mfold l
  | [] => by simp [mfold]
  | x :: xs => by
      have ih := mfold_ge_sum xs
      simp only [List.sum_cons, mfold]; unfold circ; omega
/-- The corner is at most a quarter of the square on the sum:  4ab ≤ (a+b)².
In units of the scale s, the extra pull relative to the plain sum is at most (a+b)/4 of it, so it fades as
the fields weaken.  Dark matter shows up in exactly the opposite place: where fields are weakest. -/
theorem corner_le_quarter (a b : Int) : 4 * (a * b) ≤ (a + b) ^ 2 := by
  have h : 0 ≤ (a - b) * (a - b) := by
    rcases Int.le_total 0 (a - b) with h | h
    · exact Int.mul_nonneg h h
    · exact Int.mul_nonneg_of_nonpos_of_nonpos h h
  have e : (a + b) ^ 2 = (a - b) * (a - b) + 4 * (a * b) := by grind
  omega

/-! ## 6. MOND's simple interpolation is the mult-add inverse (up to sign)
The mult-add inverse of x is y = -x/(1+x):  x ⊕ y = 0.  MOND's "simple" μ(x) = x/(1+x) is exactly -y.
Stated without division: if x ⊕ y = 0 then (1+x)(1+y) = 1 and (1+x)·(-y) = x. -/
theorem minv_lift (x y : Int) (h : x + y + x * y = 0) : (1 + x) * (1 + y) = 1 := by grind
theorem mu_is_neg_minv (x y : Int) (h : x + y + x * y = 0) : (1 + x) * (-y) = x := by grind
/-- Simple-MOND law  g_N = g·μ(g/a₀)  written without division is  g·g = g_N·g + g_N·a₀ :
the square on the observed pull equals two rectangles, g_N by g and g_N by a₀. -/
theorem simple_mond_square (g gN a0 : Int) (h : gN * (g + a0) = g * g) : g * g = gN * g + gN * a0 := by grind
/-- Deep-MOND limit: the observed pull is the side of the square whose area is the rectangle g_N × a₀. -/
theorem deep_mond_square (g gN a0 : Nat) (h : g * g = gN * a0) : g ^ 2 = gN * a0 := by rw [Nat.pow_two]; exact h

/-! ## Small kernel checks -/
example : (1 + 2) * (1 + 3) = 1 + circ 2 3 := by decide
example : mpow 1 3 = 7 := by decide
example : mpow 2 3 = 26 := by decide

end DarkSector
