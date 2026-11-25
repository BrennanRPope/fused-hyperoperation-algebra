/-
  Information theory through mult-add: the full pass
  ==================================================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`.  Intended module path:
  FusedAlgebra/Topics/InformationFull.lean (imports Base and Topics.InformationTheory).

  Builds on InfoFused (union rule, Kraft, Campbell, BSC cascades, Tsallis via
  power sums) and adds:
  1. The q-logarithm and q-exponential as the dial's own log and exp
     (stated for whole-number weights, where core Lean has the powers).
  2. Tsallis entropy as the expectation of a q-surprise s with
     s(xy) = s(x) ✳ₖ s(y): the weighted outer-log lemma, which contains
     Shannon (k = 0) and Tsallis (any k) at once.
  3. The duality q ↔ 2 − q is the dial's sign flip; Sharma–Mittal entropies.
  4. Binary trees: grafting composes internal nodes by mult-add; leaves are the lift.
  5. Any binary channel: the crossover sum a + b composes by the union rule.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/InformationFull.lean; two exact twins of InfoFused removed (see the LIBRARY comments).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.InformationTheory

namespace InformationFull

open MultAdd
open InfoFused (lsum outer shan lsum_append lsum_map shan_append)

set_option linter.unusedSectionVars false

local infixl:65 " ✳ " => multAdd

/-! # 1. The q-logarithm and q-exponential -/

/-- ln_q at whole-number weight k = m: (x^m − 1)/m. -/
def qlogN (m : Nat) (x : Rat) : Rat := (x ^ m - 1) / m

/-- The q-logarithm turns products into the dial with k = m. -/
theorem qlogN_mul (m : Nat) (hm : m ≠ 0) (x y : Rat) :
    qlogN m (x * y) = wAdd (m : Rat) (qlogN m x) (qlogN m y) := by
  unfold qlogN wAdd
  have h : (m : Rat) ≠ 0 := by
    intro e; exact hm (by exact_mod_cast e)
  rw [mul_pow']; grind

/-- Its lift is the plain power: 1 + m·ln(x) = x^m. -/
theorem qlogN_lift (m : Nat) (hm : m ≠ 0) (x : Rat) : 1 + (m : Rat) * qlogN m x = x ^ m := by
  unfold qlogN
  have h : (m : Rat) ≠ 0 := by
    intro e; exact hm (by exact_mod_cast e)
  grind

theorem qlogN_one (m : Nat) : qlogN m 1 = 0 := by
  unfold qlogN; rw [one_pow']; grind

/-- k = 1 (Tsallis q = 0): ln x = x − 1, the "count minus one". -/
theorem qlogN_one_weight (x : Rat) : qlogN 1 x = x - 1 := by
  unfold qlogN; grind

/-- The q-exponential inverts it: if y^m = 1 + m·x then ln y = x. -/
theorem qlog_of_root (m : Nat) (hm : m ≠ 0) (x y : Rat) (h : y ^ m = 1 + (m : Rat) * x) :
    qlogN m y = x := by
  unfold qlogN
  have h' : (m : Rat) ≠ 0 := by
    intro e; exact hm (by exact_mod_cast e)
  rw [h]; grind

/-- k = −1 (Tsallis q = 2): ln x = 1 − 1/x, and it obeys the union rule. -/
theorem qlog_neg_one (x y : Rat) (hx : x ≠ 0) (hy : y ≠ 0) :
    1 - 1 / (x * y) = wAdd (-1) (1 - 1 / x) (1 - 1 / y) := by
  unfold wAdd
  have : x * y ≠ 0 := by grind
  grind

/-- The dual index 2 − q: weight −m, ln(x) = ((1/x)^m − 1)/(−m). -/
def qlogNeg (m : Nat) (x : Rat) : Rat := (x⁻¹ ^ m - 1) / (-(m : Rat))

/-- Duality: ln_{2−q}(x) = −ln_q(1/x). -/
theorem qlog_dual (m : Nat) (x : Rat) : qlogNeg m x = -(qlogN m x⁻¹) := by
  unfold qlogNeg qlogN; grind

/-- The dual logarithm turns products into the dial with the opposite weight. -/
theorem qlogNeg_mul (m : Nat) (hm : m ≠ 0) (x y : Rat) (hx : x ≠ 0) (hy : y ≠ 0) :
    qlogNeg m (x * y) = wAdd (-(m : Rat)) (qlogNeg m x) (qlogNeg m y) := by
  rw [qlog_dual, qlog_dual, qlog_dual]
  have e : (x * y)⁻¹ = x⁻¹ * y⁻¹ := by grind
  rw [e, qlogN_mul m hm]
  unfold wAdd; grind

/-- The q-exponential at weight 1/n: (1 + x/n)^n. -/
def qexpN (n : Nat) (x : Rat) : Rat := (1 + x / n) ^ n

/-- It turns the dial into multiplication: e(x)·e(y) = e(x ✳_{1/n} y). -/
theorem qexpN_mul (n : Nat) (hn : n ≠ 0) (x y : Rat) :
    qexpN n x * qexpN n y = qexpN n (wAdd (1 / (n : Rat)) x y) := by
  unfold qexpN
  rw [← mul_pow']
  congr 1
  have h : (n : Rat) ≠ 0 := by
    intro e; exact hn (by exact_mod_cast e)
  unfold wAdd; grind

/-- Compound interest: e_q(x) − 1 is n mult-add steps of x/n (rung 2). -/
theorem qexp_climb (n : Nat) (x : Rat) : qexpN n x - 1 = mpow (x / n) n := by
  unfold qexpN; rw [mpow_eq]

theorem qexpN_zero (n : Nat) : qexpN n 0 = 1 := by
  unfold qexpN
  have : (1 : Rat) + 0 / n = 1 := by grind
  rw [this, one_pow']

/-! # 2. Tsallis entropy as an expectation: the weighted outer-log lemma -/

theorem shan_map_w (s : Rat → Rat) (k : Rat) (hs : ∀ x y, s (x * y) = wAdd k (s x) (s y))
    (x : Rat) (ys : List Rat) :
    shan s (ys.map (fun y => x * y))
      = x * lsum ys * s x + x * shan s ys + k * (x * s x) * shan s ys := by
  induction ys with
  | nil => simp only [List.map_nil, shan, lsum]; grind
  | cons y ys ih => simp only [List.map_cons, shan, lsum]; rw [ih, hs]; unfold wAdd; grind

/-- For any surprise s that turns products into ⊞ₖ, the joint "entropy" Σ p·s(p)
    of independent sources is
      (Σr)·S_A + (Σp)·S_B + k·S_A·S_B.
    k = 0 is the information page's outer-log lemma (Shannon); any k is Tsallis. -/
theorem outer_shan_w (s : Rat → Rat) (k : Rat) (hs : ∀ x y, s (x * y) = wAdd k (s x) (s y))
    (xs ys : List Rat) :
    shan s (outer xs ys) = lsum ys * shan s xs + lsum xs * shan s ys + k * shan s xs * shan s ys := by
  induction xs with
  | nil => simp only [outer, shan, lsum]; grind
  | cons x xs ih =>
    simp only [outer, shan, lsum]; rw [shan_append, shan_map_w s k hs, ih]; grind

/-- Tsallis' pseudo-additivity, proved as "expectation passes through the dial":
    S(A × B) = S(A) ✳ₖ S(B) for normalized independent sources. -/
theorem tsallis_engine (s : Rat → Rat) (k : Rat) (hs : ∀ x y, s (x * y) = wAdd k (s x) (s y))
    (xs ys : List Rat) (hx : lsum xs = 1) (hy : lsum ys = 1) :
    shan s (outer xs ys) = wAdd k (shan s xs) (shan s ys) := by
  rw [outer_shan_w s k hs, hx, hy]; unfold wAdd; grind

-- LIBRARY: shannon_engine (Shannon entropy is additive, the k = 0 corner) was an exact twin of InfoFused.shannon_additive and is removed.

-- LIBRARY: joint_mass (the joint law of independent sources sums to 1) was an exact twin of InfoFused.joint_sums_to_one and is removed.

/-! # 3. Sharma–Mittal: an additive entropy pushed through the dial's exponential -/

/-- If R adds over independent sources and E turns + into ⊞ₖ, then E∘R composes
    by ⊞ₖ.  With R = Rényi of order α and E(x) = (e^{kx} − 1)/k this is
    Sharma–Mittal H_{α,β}, k = 1 − β. -/
theorem sm_compose (E : Rat → Rat) (k : Rat) (hE : ∀ x y, E (x + y) = wAdd k (E x) (E y))
    (RA RB RAB : Rat) (hR : RAB = RA + RB) : E RAB = wAdd k (E RA) (E RB) := by
  rw [hR, hE]

/-- Flipping the weight: the dial at −k is the dial at k read on negatives. -/
theorem flip_weight (k a b : Rat) : wAdd (-k) a b = -(wAdd k (-a) (-b)) := by
  unfold wAdd; grind

/-! # 4. Binary trees: grafting is mult-add -/

inductive BT where
  | leaf : BT
  | node : BT → BT → BT

def leaves : BT → Nat
  | .leaf => 1
  | .node l r => leaves l + leaves r

def inner : BT → Nat
  | .leaf => 0
  | .node l r => inner l + inner r + 1

/-- Graft a copy of b onto every leaf of a. -/
def graft : BT → BT → BT
  | .leaf, b => b
  | .node l r, b => .node (graft l b) (graft r b)

/-- The lift: leaves = 1 + inner nodes. -/
theorem leaves_lift (t : BT) : leaves t = 1 + inner t := by
  induction t with
  | leaf => rfl
  | node l r ihl ihr => simp only [leaves, inner]; omega

theorem graft_leaves (a b : BT) : leaves (graft a b) = leaves a * leaves b := by
  induction a with
  | leaf => simp only [graft, leaves]; omega
  | node l r ihl ihr => simp only [graft, leaves]; rw [ihl, ihr, Nat.add_mul]

theorem graft_inner (a b : BT) : inner (graft a b) = inner a + leaves a * inner b := by
  induction a with
  | leaf => simp only [graft, inner, leaves]; omega
  | node l r ihl ihr =>
    simp only [graft, inner, leaves]; rw [ihl, ihr, Nat.add_mul]; omega

/-- Grafting composes internal nodes (decisions) by mult-add. -/
theorem graft_multAdd (a b : BT) : inner (graft a b) = inner a ✳ inner b := by
  rw [graft_inner, leaves_lift]; unfold multAdd
  rw [Nat.add_mul, Nat.one_mul]; omega

/-- The complete tree of depth n. -/
def full : Nat → BT
  | 0 => .leaf
  | n + 1 => .node (full n) (full n)

/-- It has 2^n − 1 = the mult-add power of 1 internal nodes: the q-number [n]₂. -/
theorem full_inner (n : Nat) : inner (full n) = mpow 1 n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [full, inner, mpow]; rw [ih]; unfold multAdd; omega

/-- Kraft: weight 1/2 per level; every full binary tree has Kraft sum exactly 1. -/
def kraft : BT → Rat
  | .leaf => 1
  | .node l r => (kraft l + kraft r) / 2

theorem kraft_one (t : BT) : kraft t = 1 := by
  induction t with
  | leaf => rfl
  | node l r ihl ihr => simp only [kraft]; rw [ihl, ihr]; grind

/-- D-ary trees by counts: graft (leaves 1 + c·I) — internal nodes compose by the dial c = D − 1. -/
theorem dary_graft (c IA IB : Nat) : IA + (1 + c * IA) * IB = wAdd c IA IB := by
  unfold wAdd; rw [Nat.add_mul, Nat.one_mul, Nat.mul_assoc]; omega

/-- ... and leaves multiply: the lift of the dial. -/
theorem dary_leaves (c IA IB : Nat) : (1 + c * IA) * (1 + c * IB) = 1 + c * wAdd c IA IB :=
  (wAdd_lift c IA IB).symm

/-- The complete D-ary tree of depth n: c·[n]_D + 1 = D^n, with D = 1 + c. -/
theorem dary_full (c : Nat) (n : Nat) : 1 + c * wpow c 1 n = (1 + c * 1) ^ n := wpow_lift c 1 n

/-! # 5. Any binary channel: the crossover sum composes by the union rule -/

/-- Cascade of [[1−a, a], [b, 1−b]] and [[1−a₂, a₂], [b₂, 1−b₂]]:
    the new crossover sum is (a + b) ✳₋₁ (a₂ + b₂). -/
theorem channel_cascade (a b a2 b2 : Rat) :
    ((1 - a) * a2 + a * (1 - b2)) + (b * (1 - a2) + (1 - b) * b2)
      = wAdd (-1) (a + b) (a2 + b2) := by
  unfold wAdd; grind

/-- The determinant 1 − s multiplies. -/
theorem channel_det (a b a2 b2 : Rat) :
    ((1 - a) * (1 - a2) + a * b2) * (b * a2 + (1 - b) * (1 - b2))
      - ((1 - a) * a2 + a * (1 - b2)) * (b * (1 - a2) + (1 - b) * b2)
      = (1 - a - b) * (1 - a2 - b2) := by grind

/-- The binary symmetric channel (s = 2p): the union rule on 2p is the parity rule on p. -/
theorem bsc_from_union (p q : Rat) : wAdd (-1) (2 * p) (2 * q) = 2 * wAdd (-2) p q := by
  unfold wAdd; grind

/-- The Z-channel (b = 0): its crossover composes by the union rule. -/
theorem z_cascade (a a2 : Rat) : (1 - a) * a2 + a * (1 - 0) = wAdd (-1) a a2 := by
  unfold wAdd; grind

end InformationFull
