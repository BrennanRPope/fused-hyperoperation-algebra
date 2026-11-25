import FusedAlgebra.Base

/-! # Closing rate of the fused natural constants (core Lean 4, no Mathlib)

What is here: the structural ladder lemma behind the closed form on [1,2],
a growth template for the escape count, and the exact integer facts behind
the numbers ln 2 < 0.7, ln 4 > 1.38, 2(ln 2)^2 < 1.  The real-analysis steps
(fixed points of s, the escape-time estimate, the two sandwich bounds) are
argued on the page and checked numerically in Python, not proved here. -/

namespace FusedRate

/-- Plain iteration, defined here so nothing depends on library lemmas. -/
def iter {α : Type} (s : α → α) : Nat → α → α
  | 0, x => x
  | n + 1, x => iter s n (s x)

/-- Ladder lemma: if each rung is the one below it run after one step of `s`,
    then rung j is the bottom rung run after j steps of `s`.
    On [1,2] the fused-linear rungs satisfy this with f x = b^x and s x = b^(x-1),
    which gives H_k(b,x) = b^(s^(k-3)(x)). -/
theorem ladder {α : Type} (f s : α → α) (G : Nat → α → α)
    (h0 : ∀ x, G 0 x = f x) (hs : ∀ j x, G (j + 1) x = G j (s x)) :
    ∀ j x, G j x = f (iter s j x)
  | 0, x => by rw [h0]; rfl
  | j + 1, x => by rw [hs, ladder f s G h0 hs j (s x)]; rfl

/-- Growth template: if each step multiplies the gap by at least p/q,
    then after j steps it is at least (p/q)^j times the start. -/
theorem grow (d : Nat → Nat) (p q : Nat) (h : ∀ j, p * d j ≤ q * d (j + 1)) :
    ∀ j, p ^ j * d 0 ≤ q ^ j * d j
  | 0 => by simp
  | j + 1 => by
      have ih := grow d p q h j
      have hj := h j
      calc p ^ (j + 1) * d 0 = p * (p ^ j * d 0) := by rw [Nat.pow_succ]; rw [Nat.mul_comm (p ^ j) p, Nat.mul_assoc]
        _ ≤ p * (q ^ j * d j) := Nat.mul_le_mul_left p ih
        _ = q ^ j * (p * d j) := by rw [← Nat.mul_assoc, Nat.mul_comm p (q ^ j), Nat.mul_assoc]
        _ ≤ q ^ j * (q * d (j + 1)) := Nat.mul_le_mul_left _ hj
        _ = q ^ (j + 1) * d (j + 1) := by rw [Nat.pow_succ, Nat.mul_assoc]

/-! ### Integer facts behind the numbers in the rate

Each line is an exact integer inequality checked by the kernel. What they mean
uses two textbook facts not proved here: (1 + 1/n)^n < e for every n, and
e < 2.7182819.

* `two_lt` : 2·1000^700 < 1001^700, so 2 < (1.001)^700 < e^0.7, so ln 2 < 0.7.
* `e_pow_lt` : 27182819^69 < 2^100·10^483, so e^69 < 2^100, so ln 2 > 0.69.
Together: 1.38 < ln 4 < 1.4, so 0.714 < 1/ln 4 < 0.7247, and
2(ln 2)^2 < 0.98 < 1 (the fact that makes the upper bound work). -/

set_option exponentiation.threshold 1000 in
theorem two_lt : 2 * 1000 ^ 700 < 1001 ^ 700 := by decide
set_option exponentiation.threshold 1000 in
theorem e_pow_lt : 27182819 ^ 69 < 2 ^ 100 * 10 ^ 483 := by decide
theorem sq_bound : 2 * 70 * 70 < 100 * 100 := by decide        -- 2·(0.7)^2 < 1
-- 1/1.4 > 0.714 and 1/1.38 < 0.7247
theorem inv_lo : 714 * 14 < 10000 := by decide
theorem inv_hi : 1000000 < 7247 * 138 := by decide
end FusedRate

