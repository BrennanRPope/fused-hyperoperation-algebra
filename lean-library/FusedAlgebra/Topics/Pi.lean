/-
  Fused Hyperoperation Algebra: pi from the circle operation (called mult-add before 1.0).
  ===============================================

  Core Lean 4 only (no mathlib), no `sorry`, only Lean's standard axioms.
  COMPILED with Lean 4.37.0-rc1 against the shared library's FusedAlgebra.Base
  (Nov 25, 2026). Companion page and checks: ../pi_checks.py, ../constants_hunt.py.
  To fold into the library: copy to lean-library/FusedAlgebra/Topics/Pi.lean,
  add `FusedAlgebra.Topics.Pi` to build.sh and FusedAlgebra.lean, run ./build.sh.

  The picture. The circle operation a ◯ b = a + b + ab = (1+a)(1+b) − 1, over the Gaussian
  numbers. The turning law (landscape point p = 0, q = 1) is slope addition
  x ⊕ y = (x + y)/(1 − xy). The lens 1 + a = (1 + ix)/(1 − ix) turns the turning
  law into the circle operation, and it puts every real slope on the circle |1 + a| = 1:
  the circle of radius 1 around the wall −1. That circle meets the real line at
  exactly two points: 0 (the identity) and −2 (the parity point, slope = ∞).
  pi is how far the turning flow runs from 0 to −2, measured in the law's own
  units (half of the circle operation's period 2πi).

  Proved here (all exact, no analysis):
  * The gaussian circle operation: lift, commutative, associative, identity 0.
  * The lens: slope addition on pairs (num, den) IS multiplication of den + i·num;
    the conjugate lift multiplies too, so 1 + a(x ⊕ y) = (1 + a(x))(1 + a(y)).
  * The circle |1 + a| = 1 is closed under the circle operation (norms multiply).
  * On the real line, the circle is exactly {0, −2}; and the slope ∞ lands on −2.
  * The quarter step a = −1 + i: a ◯ a = −2 and four copies make 0
    (slope 1 ⊕ 1 = 2/0 = ∞, four slope-1 steps come back to 0).
  * Machin (1706) and the other three two-term identities Størmer (1899) proved are
    the only ones, as exact Gaussian circle-operation identities: four steps of i/5 and one
    step back of i/239 land exactly on the 45° ray.
  * Step counting: 32 circle-operation steps of i/10 pass the far side (−2 side), 31 do not;
    315 steps of i/100. So the counts 3.2, 3.15 head for pi from the circle operation alone.
  * Gregory sums nest: for 0 ≤ x ≤ 1, the even partial sums climb, the odd ones
    fall, and each odd one sits above each even one (nested intervals).
  * Kernel check: the Machin combination of Gregory sums traps a number between
    3.14159265358979323846264338327950288419716939937510 and …511.

  NOT proved here (core Lean has no reals, limits or trigonometry):
  * that these nested intervals close on pi = the half-period of the turning law,
    i.e. that Im log(1 + ix) = arctan x and arctan 1 = π/4. That is classical
    (Gregory 1671, Leibniz 1674, Machin 1706) and checked numerically to 1000
    digits in ../pi_checks.py.

  LIBRARY NOTE: copied from pi/lean/PiFromMultAdd.lean (Nov 25) as Topics/Pi.
-/
import FusedAlgebra.Base

namespace PiCircle
open CircleOp

/-! ## 1. Gaussian integers and the circle operation on them -/

/-- a Gaussian integer re + i·im -/
structure G where
  re : Int
  im : Int
deriving DecidableEq, Repr

@[ext] theorem G.ext' {a b : G} (h1 : a.re = b.re) (h2 : a.im = b.im) : a = b := by
  cases a; cases b; simp_all

def G.add (a b : G) : G := ⟨a.re + b.re, a.im + b.im⟩
def G.mul (a b : G) : G := ⟨a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re⟩
def G.one : G := ⟨1, 0⟩
def G.zero : G := ⟨0, 0⟩
def G.conj (a : G) : G := ⟨a.re, -a.im⟩
/-- squared length |a|² -/
def G.norm (a : G) : Int := a.re * a.re + a.im * a.im

/-- The gaussian circle operation a ◯ b = a + b + ab -/
def gma (a b : G) : G := G.add (G.add a b) (G.mul a b)

/-- circle power: n copies of a combined together by the circle operation -/
def gpow (a : G) : Nat → G
  | 0 => G.zero
  | n + 1 => gma (gpow a n) a

/-- ordinary power of a Gaussian integer -/
def G.pow (a : G) : Nat → G
  | 0 => G.one
  | n + 1 => G.mul (G.pow a n) a

theorem gma_lift (a b : G) : G.add G.one (gma a b) = G.mul (G.add G.one a) (G.add G.one b) := by
  ext <;> simp [gma, G.add, G.mul, G.one] <;> grind

theorem gma_comm (a b : G) : gma a b = gma b a := by
  ext <;> simp [gma, G.add, G.mul] <;> grind

theorem gma_assoc (a b c : G) : gma (gma a b) c = gma a (gma b c) := by
  ext <;> simp [gma, G.add, G.mul] <;> grind

theorem gma_zero (a : G) : gma a G.zero = a := by
  ext <;> simp [gma, G.add, G.mul, G.zero]

theorem mul_assoc' (a b c : G) : G.mul (G.mul a b) c = G.mul a (G.mul b c) := by
  ext <;> simp [G.mul] <;> grind

/-- the lift of a circle power is an ordinary power: 1 + a^(◯n) = (1 + a)ⁿ -/
theorem gpow_lift (a : G) (n : Nat) : G.add G.one (gpow a n) = G.pow (G.add G.one a) n := by
  induction n with
  | zero => ext <;> simp [gpow, G.pow, G.add, G.one, G.zero]
  | succ n ih => rw [gpow, G.pow, gma_lift, ih]

/-! ## 2. The circle around the wall -/

/-- lengths multiply (Brahmagupta–Fibonacci) -/
theorem norm_mul (a b : G) : G.norm (G.mul a b) = G.norm a * G.norm b := by
  simp [G.norm, G.mul]; grind

theorem norm_conj (a : G) : G.norm (G.conj a) = G.norm a := by
  simp [G.norm, G.conj]; grind

/-- the circle |1 + a| = 1 (centre: the wall −1, radius 1) is closed under the circle operation -/
theorem circle_closed (a b : G) (ha : G.norm (G.add G.one a) = 1) (hb : G.norm (G.add G.one b) = 1) :
    G.norm (G.add G.one (gma a b)) = 1 := by
  rw [gma_lift, norm_mul, ha, hb]; rfl

/-- On the real line (Rat), the circle (1 + x)² = 1 is exactly the two points 0 and −2. -/
theorem circle_real_points (x : Rat) : (1 + x) * (1 + x) = 1 ↔ x = 0 ∨ x = -2 := by
  constructor
  · intro h
    have h2 : x * (x + 2) = 0 := by grind
    rcases Rat.mul_eq_zero.mp h2 with h3 | h3
    · exact Or.inl h3
    · exact Or.inr (by grind)
  · rintro (h | h) <;> subst h <;> decide +kernel

/-- in circle-operation units: x ◯ x = 0 (a half turn twice is a full turn) exactly at 0 and −2 -/
theorem circ_self_zero (x : Rat) : circ x x = 0 ↔ x = 0 ∨ x = -2 := by
  rw [← circle_real_points]; unfold circ; constructor <;> intro h <;> grind

/-! ## 3. The lens: the turning law is the circle operation

  A slope x = num/den is kept as a pair, so ∞ = (n, 0) is an ordinary pair.
  Slope addition (x + y)/(1 − xy) on pairs: (p, q) ⊕ (r, s) = (ps + rq, qs − pr). -/

structure Slope where
  num : Int
  den : Int
deriving DecidableEq, Repr

def tl (x y : Slope) : Slope := ⟨x.num * y.den + y.num * x.den, x.den * y.den - x.num * y.num⟩

/-- the Gaussian number den + i·num that a slope stands for (its lift is this over its conjugate) -/
def toG (x : Slope) : G := ⟨x.den, x.num⟩

/-- THE LENS. Slope addition is exactly Gaussian multiplication. -/
theorem lens (x y : Slope) : toG (tl x y) = G.mul (toG x) (toG y) := by
  ext <;> simp [toG, tl, G.mul] <;> grind

/-- and the conjugates multiply the same way, so the lift W/W̄ of a sum is the product of
    the lifts: 1 + a(x ⊕ y) = (1 + a(x))·(1 + a(y)), i.e. a(x ⊕ y) = a(x) ◯ a(y). -/
theorem lens_conj (x y : Slope) : G.conj (toG (tl x y)) = G.mul (G.conj (toG x)) (G.conj (toG y)) := by
  ext <;> simp [toG, tl, G.mul, G.conj] <;> grind

/-- every real slope lands on the circle: W and W̄ have the same length, so |W/W̄| = 1 -/
theorem lens_on_circle (x : Slope) : G.norm (toG x) = G.norm (G.conj (toG x)) := (norm_conj _).symm

/-- slope 0 lands on the identity: W = W̄, lift 1, a = 0 -/
theorem lens_zero (d : Int) : toG ⟨0, d⟩ = G.conj (toG ⟨0, d⟩) := by
  ext <;> simp [toG, G.conj]

/-- slope ∞ lands on the parity point: W = −W̄, lift −1, a = −2 -/
theorem lens_infinity (n : Int) : toG ⟨n, 0⟩ = G.mul ⟨-1, 0⟩ (G.conj (toG ⟨n, 0⟩)) := by
  ext <;> simp [toG, G.conj, G.mul]

/-! ## 4. The quarter step -/

/-- a = −1 + i, the lens image of slope 1 (a quarter of the way round) -/
def q4 : G := ⟨-1, 1⟩

theorem quarter_twice : gma q4 q4 = ⟨-2, 0⟩ := by decide
theorem quarter_four : gpow q4 4 = G.zero := by decide
theorem quarter_on_circle : G.norm (G.add G.one q4) = 1 := by decide
/-- the same in slopes: 1 ⊕ 1 = 2/0 = ∞, and four steps of slope 1 come back to slope 0 -/
theorem slope_one_twice : tl ⟨1, 1⟩ ⟨1, 1⟩ = ⟨2, 0⟩ := by decide
theorem slope_one_four : tl (tl (tl ⟨1, 1⟩ ⟨1, 1⟩) ⟨1, 1⟩) ⟨1, 1⟩ = ⟨0, -4⟩ := by decide
/-- the parity point is its own way back: (−2) ◯ (−2) = 0 -/
theorem parity_twice : circ (-2 : Rat) (-2) = 0 := by decide +kernel

/-! ## 5. Steps that land exactly on the eighth turn (the 45° ray, slope 1)

  Circle-operation steps of i/m are, after scaling by m, the Gaussian numbers m + i.
  A mix of steps lands on the 45° ray when the product is k·(1 + i), k > 0. -/

/-- Machin (1706): four steps of i/5, one step back of i/239 -/
theorem machin : G.mul (G.pow ⟨5, 1⟩ 4) ⟨239, -1⟩ = ⟨114244, 114244⟩ := by decide
/-- Euler: one step of i/2, one of i/3 -/
theorem euler : G.mul ⟨2, 1⟩ ⟨3, 1⟩ = ⟨5, 5⟩ := by decide
/-- Hermann: two steps of i/2, one back of i/7 -/
theorem hermann : G.mul (G.pow ⟨2, 1⟩ 2) ⟨7, -1⟩ = ⟨25, 25⟩ := by decide
/-- Hutton: two steps of i/3, one of i/7 -/
theorem hutton : G.mul (G.pow ⟨3, 1⟩ 2) ⟨7, 1⟩ = ⟨50, 50⟩ := by decide

/-- Machin again, in the turning law itself: 1/5 ⊕ 1/5 ⊕ 1/5 ⊕ 1/5 ⊕ (−1/239) = 1 -/
theorem machin_slopes :
    tl (tl (tl (tl ⟨1, 5⟩ ⟨1, 5⟩) ⟨1, 5⟩) ⟨1, 5⟩) ⟨-1, 239⟩ = ⟨114244, 114244⟩ := by decide

/-! ## 6. Counting tiny turning steps until we pass the far side

  (1 + i/n)^m = (n + i)^m / n^m, so the sign of the imaginary part of (n + i)^m tells
  which side of the real line m circle-operation steps of i/n have reached. -/

/-- the first m ≤ fuel at which (n + i)^m has negative imaginary part (0 if none) -/
def firstCross (n : Int) : Nat → Nat → G → Nat
  | 0, _, _ => 0
  | fuel + 1, m, z =>
      let z' := G.mul z ⟨n, 1⟩
      if z'.im < 0 then m + 1 else firstCross n fuel (m + 1) z'

/-- 32 steps of i/10 pass −2; 31 do not -/
theorem cross_10 : firstCross 10 100 0 G.one = 32 := by decide +kernel
/-- 315 steps of i/100 -/
theorem cross_100 : firstCross 100 1000 0 G.one = 315 := by decide +kernel
/-- 3142 steps of i/1000 -/
theorem cross_1000 : firstCross 1000 4000 0 G.one = 3142 := by decide +kernel

/-! ## 7. The circle operation's logarithm, Gregory's sums, and nested intervals

  The way back from the circle operation to addition is L(a) = a − a²/2 + a³/3 − …  At a = i·x its
  imaginary part is x − x³/3 + x⁵/5 − … (Gregory's series). -/

/-- Gregory's k-th term x^(2k+1)/(2k+1) -/
def gt (x : Rat) (k : Nat) : Rat := x ^ (2 * k + 1) / ((2 * k + 1 : Nat) : Rat)

/-- Gregory's partial sums: Σ_{k<n} (−1)^k x^(2k+1)/(2k+1) -/
def greg (x : Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => greg x n + (-1) ^ n * gt x n

/-- a'/d ≤ a/b when 0 < b ≤ d and 0 ≤ a' ≤ a -/
theorem div_le_div_of {a a' b d : Rat} (hb : 0 < b) (hbd : b ≤ d) (ha' : 0 ≤ a') (haa : a' ≤ a) :
    a' / d ≤ a / b := by
  have hd : 0 < d := by grind
  have e1 : a' / d * d = a' := Rat.div_mul_cancel (by grind)
  have e2 : a / b * b = a := Rat.div_mul_cancel (by grind)
  apply Classical.byContradiction
  intro h
  have h' : a / b < a' / d := by grind
  have hq : 0 ≤ a' / d := by
    apply Classical.byContradiction; intro hn
    have : a' / d < 0 := by grind
    have := Rat.mul_lt_mul_of_pos_left this hd
    grind
  have s1 := Rat.mul_lt_mul_of_pos_left h' hb
  have s2 := Rat.mul_le_mul_of_nonneg_left hbd hq
  grind

theorem neg_one_pow_even (m : Nat) : ((-1 : Rat)) ^ (2 * m) = 1 := by
  induction m with
  | zero => decide +kernel
  | succ m ih =>
    rw [show 2 * (m + 1) = 2 * m + 1 + 1 by omega, Rat.pow_succ, Rat.pow_succ, ih]; decide +kernel

theorem neg_one_pow_odd (m : Nat) : ((-1 : Rat)) ^ (2 * m + 1) = -1 := by
  rw [Rat.pow_succ, neg_one_pow_even]; decide +kernel

theorem gt_nonneg (x : Rat) (h0 : 0 ≤ x) (k : Nat) : 0 ≤ gt x k := by
  unfold gt
  have hp : 0 ≤ x ^ (2 * k + 1) := Rat.pow_nonneg h0
  have hd : (0 : Rat) < ((2 * k + 1 : Nat) : Rat) := by exact_mod_cast (by omega : 0 < 2 * k + 1)
  have := div_le_div_of hd (Rat.le_refl) (Rat.le_refl) hp
  have e : (0 : Rat) / ((2 * k + 1 : Nat) : Rat) = 0 := by grind
  grind

/-- the terms shrink for 0 ≤ x ≤ 1 -/
theorem term_shrink (x : Rat) (h0 : 0 ≤ x) (h1 : x ≤ 1) (k : Nat) : gt x (k + 1) ≤ gt x k := by
  unfold gt
  have hp : 0 ≤ x ^ (2 * k + 1) := Rat.pow_nonneg h0
  have hx2 : x * x ≤ 1 := by
    have := Rat.mul_le_mul_of_nonneg_left h1 h0
    grind
  have hpow : x ^ (2 * (k + 1) + 1) = x ^ (2 * k + 1) * (x * x) := by
    rw [show 2 * (k + 1) + 1 = 2 * k + 1 + 1 + 1 by omega, Rat.pow_succ, Rat.pow_succ]; grind
  have hle : x ^ (2 * (k + 1) + 1) ≤ x ^ (2 * k + 1) := by
    rw [hpow]
    have := Rat.mul_le_mul_of_nonneg_left hx2 hp
    grind
  have hb : (0 : Rat) < ((2 * k + 1 : Nat) : Rat) := by exact_mod_cast (by omega : 0 < 2 * k + 1)
  have hbd : ((2 * k + 1 : Nat) : Rat) ≤ ((2 * (k + 1) + 1 : Nat) : Rat) := by
    exact_mod_cast (by omega : 2 * k + 1 ≤ 2 * (k + 1) + 1)
  exact div_le_div_of hb hbd (Rat.pow_nonneg h0) hle

theorem greg_two (x : Rat) (m : Nat) :
    greg x (2 * m + 2) = greg x (2 * m) + gt x (2 * m) - gt x (2 * m + 1) := by
  rw [show 2 * m + 2 = (2 * m + 1) + 1 by omega, greg, greg, neg_one_pow_even, neg_one_pow_odd]
  grind

/-- even partial sums climb -/
theorem greg_even_up (x : Rat) (h0 : 0 ≤ x) (h1 : x ≤ 1) (m : Nat) :
    greg x (2 * m) ≤ greg x (2 * m + 2) := by
  rw [greg_two]
  have := term_shrink x h0 h1 (2 * m)
  grind

/-- each odd partial sum sits above the even one before it -/
theorem greg_odd_above (x : Rat) (h0 : 0 ≤ x) (m : Nat) :
    greg x (2 * m) ≤ greg x (2 * m + 1) := by
  rw [greg, neg_one_pow_even]
  have := gt_nonneg x h0 (2 * m)
  grind

/-- odd partial sums fall -/
theorem greg_odd_down (x : Rat) (h0 : 0 ≤ x) (h1 : x ≤ 1) (m : Nat) :
    greg x (2 * m + 3) ≤ greg x (2 * m + 1) := by
  rw [show 2 * m + 3 = (2 * m + 2) + 1 by omega, greg, greg_two, greg,
    show 2 * m + 2 = 2 * (m + 1) by omega, neg_one_pow_even, neg_one_pow_even]
  have := term_shrink x h0 h1 (2 * m + 1)
  rw [show 2 * m + 1 + 1 = 2 * (m + 1) by omega] at this
  grind

/-- so [even sum, odd sum] is a nested interval: every even sum sits below every odd sum -/
theorem greg_even_le_odd (x : Rat) (h0 : 0 ≤ x) (h1 : x ≤ 1) (m n : Nat) :
    greg x (2 * m) ≤ greg x (2 * n + 1) := by
  have up : ∀ a b, a ≤ b → greg x (2 * a) ≤ greg x (2 * b) := by
    intro a b hab
    induction b with
    | zero => have : a = 0 := by omega
              subst this; exact Rat.le_refl
    | succ b ih =>
      rcases Nat.lt_or_ge a (b + 1) with h | h
      · have := ih (by omega)
        have := greg_even_up x h0 h1 b
        rw [show 2 * (b + 1) = 2 * b + 2 by omega]; exact Rat.le_trans ‹_› this
      · have : a = b + 1 := by omega
        subst this; exact Rat.le_refl
  have down : ∀ a b, a ≤ b → greg x (2 * b + 1) ≤ greg x (2 * a + 1) := by
    intro a b hab
    induction b with
    | zero => have : a = 0 := by omega
              subst this; exact Rat.le_refl
    | succ b ih =>
      rcases Nat.lt_or_ge a (b + 1) with h | h
      · have := ih (by omega)
        have := greg_odd_down x h0 h1 b
        rw [show 2 * (b + 1) + 1 = 2 * b + 3 by omega]; exact Rat.le_trans this ‹_›
      · have : a = b + 1 := by omega
        subst this; exact Rat.le_refl
  rcases Nat.lt_or_ge m n with h | h
  · exact Rat.le_trans (up m n (by omega)) (greg_odd_above x h0 n)
  · exact Rat.le_trans (greg_odd_above x h0 m) (down n m h)

/-! ## 8. The Machin interval, computed by the kernel

  Lower: 16·(even sum at 1/5) − 4·(odd sum at 1/239); upper: the other way round.
  With 44 and 45 terms the two ends agree to 61 digits. -/

def machinLo (N : Nat) : Rat := 16 * greg (1 / 5) (2 * N) - 4 * greg (1 / 239) (2 * N + 1)
def machinHi (N : Nat) : Rat := 16 * greg (1 / 5) (2 * N + 1) - 4 * greg (1 / 239) (2 * N)

theorem machin_lo_le_hi (N : Nat) : machinLo N ≤ machinHi N := by
  unfold machinLo machinHi
  have a := greg_odd_above (1 / 5) (by decide +kernel) N
  have b := greg_odd_above (1 / 239) (by decide +kernel) N
  grind

/-- the Machin intervals nest: the lower end climbs and the upper end falls with N -/
theorem machin_nested (N : Nat) : machinLo N ≤ machinLo (N + 1) ∧ machinHi (N + 1) ≤ machinHi N := by
  unfold machinLo machinHi
  have h5 : (0 : Rat) ≤ 1 / 5 := by decide +kernel
  have h5' : (1 / 5 : Rat) ≤ 1 := by decide +kernel
  have h2 : (0 : Rat) ≤ 1 / 239 := by decide +kernel
  have h2' : (1 / 239 : Rat) ≤ 1 := by decide +kernel
  have a := greg_even_up (1 / 5) h5 h5' N
  have b := greg_odd_down (1 / 239) h2 h2' N
  have c := greg_odd_down (1 / 5) h5 h5' N
  have d := greg_even_up (1 / 239) h2 h2' N
  rw [show 2 * (N + 1) = 2 * N + 2 by omega, show 2 * N + 2 + 1 = 2 * N + 3 by omega]
  constructor <;> grind

/-- 50 decimals: 3.14159265358979323846264338327950288419716939937510 < lo ≤ hi < …511 -/
theorem machin_50_digits :
    (314159265358979323846264338327950288419716939937510 / 10 ^ 50 : Rat) < machinLo 22 ∧
    machinHi 22 < (314159265358979323846264338327950288419716939937511 / 10 ^ 50 : Rat) := by
  decide +kernel

end PiCircle
