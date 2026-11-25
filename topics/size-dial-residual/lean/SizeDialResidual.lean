/-
  Size weighting as a dial, and what the circle operation (called mult-add before 1.0) leaves unexplained
  ==============================================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  math-theory/size-dial-residual/ (index.html, checks.py).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base only.

  What is proved (the elementary identities behind the page):

  1. A size weight is a per-element twist. Weighing an element by w and
     combining by the circle operation equals the twisted (weighted) circle operation at k = w, scaled by w
     (`weight_is_twist`, Base's circK_scale). Union of two weighted dials differs
     from weighing the union dial by exactly t·s·w·(w − 1) (`weighted_merge_gap`),
     so the two commute when w·w = w, the sharp count (`sharp_weight_commutes`).
     Under a twist k the half dial h, with (1 + kh)² = 1 + k, is a half:
     k·(h ◯ₖ h) = k (`twisted_half`).
  2. The residual is repetition at doubled size. One element used any number of
     times (1/(1 − x)) is the merge of "at most once" factors at x, x², x⁴, …;
     the base fact is (1 + x)(1 − x) = 1 − x², and over a list of elements the
     subsets times "everything at doubled size" give the multisets
     (`distinct_times_doubled`). The finite η–ζ law: the alternating sum is the
     plain sum minus twice the even terms (`alt_eq_sum_sub_evens`), and when
     doubling an element scales its weight by c (= 2^(−s)) the even terms are c
     times the whole (`alt_doubling`), which is the 1 − 2^(−β) share of the page.
     At most m copies: (1 − x)(1 + x + … + x^m) = 1 − x^(m+1) (`at_most_m`).
  3. The class-group residual is a sign (character) dial. Over a list of
     elements with dials x, the products ∏(1 + x) and ∏(1 − x) are the even part
     ± the odd part (`prod_plus_even_odd`, `prod_minus_even_odd`), so the even
     part (products of an even number of non-principal primes: the principal
     ones when the class group has order 2) is half their sum (`even_half`).
     Two elements: (1 + x)(1 + y) + (1 − x)(1 − y) = 2(1 + xy) (`two_sign_dials`).
  The polylogarithm limits (Li_{β+1}(−t)/Li_{β+1}(−1)) are computed in checks.py,
  not proved here.
-/
import FusedAlgebra.Base

namespace SizeDialResidual

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. Size weight = per-element twist -/

section Weight
variable {R : Type} [Lean.Grind.CommRing R]

/-- **A size weight is a twist.** Weighing two dials by w and combining by the circle operation is the
    twisted circle operation at k = w, scaled by w. -/
theorem weight_is_twist (w t s : R) : circ (w * t) (w * s) = w * circK w t s :=
  (circK_scale w t s).symm

/-- **Weighing does not commute with merging.** Union of the weighted dials vs the
    weighted union dial: the gap is t·s·w·(w − 1). -/
theorem weighted_merge_gap (t s w : R) :
    circ (t * w) (s * w) = circ t s * w + t * s * w * (w - 1) := by
  simp only [circ]; grind

/-- **Only a sharp count commutes.** If w·w = w (w is 0 or 1 in a domain), weighing
    and merging commute. -/
theorem sharp_weight_commutes (t s w : R) (h : w * w = w) :
    circ (t * w) (s * w) = circ t s * w := by
  rw [weighted_merge_gap]
  have : t * s * w * (w - 1) = t * s * (w * w - w) := by grind
  rw [this, h]; grind

/-- Weight 1 is the plain count. -/
theorem weight_one (t s : R) : circ (t * 1) (s * 1) = circ t s * 1 :=
  sharp_weight_commutes t s 1 (by grind)

/-- **The half under a twist.** If (1 + k·h)² = 1 + k then k·(h ◯ₖ h) = k. -/
theorem twisted_half (k h : R) (hh : (1 + k * h) * (1 + k * h) = 1 + k) :
    k * circK k h h = k := by
  have l := circK_lift k h h
  grind

/-- Dial 1/2 under weight w is the dial w/2: the uniform dial and the weight multiply. -/
theorem stacked_dials (t w : R) : 1 + t * w = 1 + w * t := by grind

end Weight

/-! ## 2. The residual: repetition at doubled size -/

section Residual
variable {R : Type} [Lean.Grind.CommRing R]

/-- Subsets times multisets-at-doubled-size: element by element,
    ∏(1 + x) · ∏(1 − x²)⁻¹… in the finite form ∏(1 + x) · ∏(1 − x) = ∏(1 − x²). -/
theorem distinct_times_doubled (l : List R) :
    prodL (l.map (fun x => 1 + x)) * prodL (l.map (fun x => 1 - x))
      = prodL (l.map (fun x => 1 - x * x)) := by
  induction l with
  | nil => simp only [List.map, prodL]; grind
  | cons a l ih => simp only [List.map, prodL]; rw [← ih]; grind

/-- Sum f(1) + … + f(n). -/
def S (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => S f n + f (n + 1)

/-- Alternating sum over 1..2N: f(1) − f(2) + f(3) − … − f(2N). -/
def A (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => A f n + f (2 * n + 1) - f (2 * n + 2)

/-- Plain sum over 1..2N. -/
def S2 (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => S2 f n + f (2 * n + 1) + f (2 * n + 2)

/-- Even terms f(2) + f(4) + … + f(2N). -/
def E (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => E f n + f (2 * n + 2)

theorem S2_eq_S (f : Nat → R) (N : Nat) : S2 f N = S f (2 * N) := by
  induction N with
  | zero => rfl
  | succ n ih =>
    have h : 2 * (n + 1) = (2 * n + 1) + 1 := by omega
    rw [h]; simp only [S2, S]; rw [ih]

/-- **The finite η–ζ law.** Alternating sum = plain sum − 2 × even terms. -/
theorem alt_eq_sum_sub_evens (f : Nat → R) (N : Nat) :
    A f N = S2 f N - 2 * E f N := by
  induction N with
  | zero => simp only [A, S2, E]; grind
  | succ n ih => simp only [A, S2, E]; rw [ih]; grind

/-- If doubling scales a weight by c, the even terms are c times the first half. -/
theorem evens_scale (f : Nat → R) (c : R) (hc : ∀ n, f (2 * n) = c * f n) (N : Nat) :
    E f N = c * S f N := by
  induction N with
  | zero => simp only [E, S]; grind
  | succ n ih =>
    simp only [E, S]; rw [ih]
    have : f (2 * n + 2) = c * f (n + 1) := by
      have := hc (n + 1); rwa [show 2 * (n + 1) = 2 * n + 2 by omega] at this
    rw [this]; grind

/-- **Doubling law (η = (1 − 2c)·ζ in finite form).** With c = 2^(−s): the
    "at most once" share is the whole minus twice c times the first half. -/
theorem alt_doubling (f : Nat → R) (c : R) (hc : ∀ n, f (2 * n) = c * f n) (N : Nat) :
    A f N = S f (2 * N) - 2 * c * S f N := by
  rw [alt_eq_sum_sub_evens, S2_eq_S, evens_scale f c hc]; grind

/-- Geometric sum 1 + x + … + x^m. -/
def geo (x : R) : Nat → R
  | 0 => 1
  | m + 1 => geo x m + x ^ (m + 1)

/-- **At most m copies.** (1 − x)(1 + x + … + x^m) = 1 − x^(m+1); m = 1 is the circle operation's
    "at most once", m → ∞ is full repetition. -/
theorem at_most_m (x : R) (m : Nat) : (1 - x) * geo x m = 1 - x ^ (m + 1) := by
  induction m with
  | zero => simp only [geo]; rw [Semiring.pow_succ, Semiring.pow_zero]; grind
  | succ m ih =>
    simp only [geo]
    have e : x ^ (m + 1 + 1) = x ^ (m + 1) * x := Semiring.pow_succ x (m + 1)
    rw [e]; grind

/-- m = 1 is the circle-operation factor. -/
theorem at_most_one (x : R) : geo x 1 = 1 + x := by
  simp only [geo]; rw [Semiring.pow_succ, Semiring.pow_zero]; grind

end Residual

/-! ## 3. The class-group residual: sign dials -/

section Signs
variable {R : Type} [Lean.Grind.CommRing R]

mutual
/-- Products of an even number of the dials (sum over even-size subsets). -/
def ev : List R → R
  | [] => 1
  | x :: l => ev l + x * od l
/-- Products of an odd number of the dials. -/
def od : List R → R
  | [] => 0
  | x :: l => od l + x * ev l
end

theorem prod_plus_even_odd (l : List R) :
    prodL (l.map (fun x => 1 + x)) = ev l + od l := by
  induction l with
  | nil => simp only [List.map, prodL, ev, od]; grind
  | cons a l ih => simp only [List.map, prodL, ev, od]; rw [ih]; grind

theorem prod_minus_even_odd (l : List R) :
    prodL (l.map (fun x => 1 - x)) = ev l - od l := by
  induction l with
  | nil => simp only [List.map, prodL, ev, od]; grind
  | cons a l ih => simp only [List.map, prodL, ev, od]; rw [ih]; grind

/-- **Class group of order 2.** The principal part (even number of non-principal
    primes) is half of ∏(1 + x) + ∏(1 − x): the sign dial is the missing piece. -/
theorem even_half (l : List R) :
    prodL (l.map (fun x => 1 + x)) + prodL (l.map (fun x => 1 - x)) = 2 * ev l := by
  rw [prod_plus_even_odd, prod_minus_even_odd]; grind

theorem two_sign_dials (x y : R) : (1 + x) * (1 + y) + (1 - x) * (1 - y) = 2 * (1 + x * y) := by
  grind

/-- Two non-principal primes make a principal product (the xy term), one alone does not. -/
theorem two_sign_dials_example (x y : R) : ev [x, y] = 1 + x * y := by
  simp only [ev, od]; grind

end Signs

end SizeDialResidual
