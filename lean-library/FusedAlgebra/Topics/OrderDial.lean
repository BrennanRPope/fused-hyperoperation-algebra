/-
  Fused Hyperoperation Algebra: order, composition and the collection laws.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base.
  Intended library path: FusedAlgebra/Topics/OrderDial.lean (imports only FusedAlgebra.Base).

  1. Where the circle operation (called mult-add before 1.0) sits on the cube of laws: its only idempotents are 0 and the wall.
  2. Composition first: one compounding period x ↦ (1 + r)x + d. In circle-operation coordinates the
     stretches compose by ◯ and the shifts pick up order; the order gap is the wedge
     r₁d₂ − r₂d₁; maps commute exactly when the wedge is 0; the k-family is one such slice.
  3. Order inside the formula: a + b + αab + βba. Its commutator is (α − β)[a, b]; on
     2×2 matrices its grouping error is αβ times a fixed nonzero matrix, so it is
     associative there only when αβ = 0.
  4. Anticommuting forces "at most once"; Pauli X and Z.
  5. The discrete Weyl relation: shift and multiply-by-λⁿ commute up to the factor λ^a.
  6. Kleisli composition: the k-sum is a semiring sum only at k = 0 or on 0/1 values,
     where k = −1 is OR (relations).
  7. Dictionary merge (last write wins) is a right regular band, not commutative.
  8. Arrangements of x x y y by out-of-order pairs: 1, 1, 2, 1, 1.
  9. No functor half power set: the size skeleton of H ∘ H ≅ P fails at sizes 0, 1, 2.
-/
import FusedAlgebra.Base

namespace OrderDial
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. The circle operation on the cube of laws -/

/-- The grouping-free, order-free, but repetition-counting corner: a ◯ₖ a = a only at 0 and the wall. -/
theorem circK_idem_iff (k a : Int) : circK k a a = a ↔ a = 0 ∨ 1 + k * a = 0 := by
  have e : circK k a a - a = a * (1 + k * a) := by simp only [circK]; grind
  constructor
  · intro h
    have : a * (1 + k * a) = 0 := by rw [← e, h]; grind
    exact Int.mul_eq_zero.mp this
  · intro h
    have : a * (1 + k * a) = 0 := by
      rcases h with h | h <;> simp [h]
    grind

/-- So the circle operation counts repeats: 1 ◯ 1 = 3. -/
theorem one_circ_one : (1 : Int) ◯ 1 = 3 := by decide

/-- The wall −1 absorbs everything (the circle operation's zero, not its unit). -/
theorem wall_absorbs (b : Int) : (-1 : Int) ◯ b = -1 := by simp only [circ]; grind

/-! ## 2. Composition first: stretch-and-shift maps -/

/-- A stretch-and-shift in circle-operation coordinates: (a, c) is x ↦ (1 + a)x + c. -/
def app (p : Int × Int) (x : Int) : Int := (1 + p.1) * x + p.2

/-- Composition of two stretch-and-shifts: stretches combine by the circle operation, shifts pick up order. -/
def comp (p q : Int × Int) : Int × Int := (p.1 ◯ q.1, p.2 + (1 + p.1) * q.2)

/-- `comp` really is composition (functoriality of the one-object category). -/
theorem app_comp (p q : Int × Int) (x : Int) : app (comp p q) x = app p (app q x) := by
  simp only [app, comp, circ]; grind

theorem comp_id (p : Int × Int) : comp p (0, 0) = p ∧ comp (0, 0) p = p := by
  simp only [comp, circ]; constructor <;> (apply Prod.ext <;> simp <;> grind)

theorem comp_assoc (p q r : Int × Int) : comp (comp p q) r = comp p (comp q r) := by
  simp only [comp, circ]; apply Prod.ext <;> simp <;> grind

/-- The stretch parts never care about order. -/
theorem stretch_order_free (p q : Int × Int) : (comp p q).1 = (comp q p).1 := circ_comm _ _

/-- The order gap is a pure shift, equal to the wedge a·d − b·c. -/
theorem order_gap (p q : Int × Int) : (comp p q).2 - (comp q p).2 = p.1 * q.2 - q.1 * p.2 := by
  simp only [comp]; grind

/-- Two stretch-and-shifts commute exactly when their wedge vanishes. -/
theorem commute_iff (p q : Int × Int) : comp p q = comp q p ↔ p.1 * q.2 = q.1 * p.2 := by
  constructor
  · intro h; have := order_gap p q; rw [h] at this; grind
  · intro h; apply Prod.ext
    · exact stretch_order_free p q
    · have := order_gap p q; grind

/-- The k-family (stretch k·a, shift a; all about the wall −1/k) is a commuting slice, and it is ◯ₖ. -/
theorem k_slice (k a b : Int) : comp (k * a, a) (k * b, b) = (k * circK k a b, circK k a b) := by
  simp only [comp, circ, circK]; apply Prod.ext <;> simp <;> grind

theorem k_slice_commutes (k a b : Int) : comp (k * a, a) (k * b, b) = comp (k * b, b) (k * a, a) := by
  rw [commute_iff]; simp; grind

/-- Two stretches by 2 about different centres (0 and 1) do not commute. -/
theorem different_centres : comp (1, 0) (1, -1) ≠ comp (1, -1) (1, 0) := by decide

/-- Compounding: deposit-then-grow minus grow-then-deposit is the corner r·d. -/
theorem compounding_corner (B r d : Int) : (B + d) * (1 + r) - (B * (1 + r) + d) = r * d := by grind

/-- A period's centre: if r·x* = −d then x* is fixed (the interest-only balance). -/
theorem centre_fixed (r d x : Int) (h : r * x = -d) : app (r, d) x = x := by
  simp only [app]; grind

/-! ## 3. Order inside the formula costs grouping -/

/-- 2×2 integer matrices (local copy, rows (a b) (c d)). -/
structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq, Repr

namespace M2
def add (x y : M2) : M2 := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def sub (x y : M2) : M2 := ⟨x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d⟩
def mul (x y : M2) : M2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def smul (s : Int) (x : M2) : M2 := ⟨s * x.a, s * x.b, s * x.c, s * x.d⟩
def neg (x : M2) : M2 := ⟨-x.a, -x.b, -x.c, -x.d⟩
theorem ext' {x y : M2} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c) (h4 : x.d = y.d) : x = y := by
  cases x; cases y; simp_all
end M2

open M2 in
/-- The order-weighted circle operation: a + b + α·ab + β·ba. -/
def mutAdd (al be : Int) (x y : M2) : M2 := add (add x y) (add (smul al (mul x y)) (smul be (mul y x)))

open M2 in
/-- Its order error is (α − β) times the commutator: the order part is a Lie bracket. -/
theorem mut_comm (al be : Int) (x y : M2) :
    sub (mutAdd al be x y) (mutAdd al be y x) = smul (al - be) (sub (mul x y) (mul y x)) := by
  apply ext' <;> simp only [mutAdd, add, sub, smul, mul] <;> grind

open M2 in
/-- Every weighting stays flexible: (a ◯ b) ◯ a = a ◯ (b ◯ a). -/
theorem mut_flexible (al be : Int) (x y : M2) :
    mutAdd al be (mutAdd al be x y) x = mutAdd al be x (mutAdd al be y x) := by
  apply ext' <;> simp only [mutAdd, add, smul, mul] <;> grind

def E11 : M2 := ⟨1, 0, 0, 0⟩
def E12 : M2 := ⟨0, 1, 0, 0⟩
def E21 : M2 := ⟨0, 0, 1, 0⟩

open M2 in
/-- On E11, E12, E21 the grouping error is αβ times the fixed matrix [[1,0],[0,−1]]. -/
theorem mut_assoc_error (al be : Int) :
    sub (mutAdd al be (mutAdd al be E12 E21) E11) (mutAdd al be E12 (mutAdd al be E21 E11)) = smul (al * be) ⟨1, 0, 0, -1⟩ := by
  apply ext' <;> simp only [mutAdd, add, sub, smul, mul, E11, E12, E21] <;> grind

open M2 in
/-- So order-in-the-formula is associative on these three only at the ends αβ = 0. -/
theorem mut_assoc_only_at_ends (al be : Int) :
    mutAdd al be (mutAdd al be E12 E21) E11 = mutAdd al be E12 (mutAdd al be E21 E11) ↔ al * be = 0 := by
  constructor
  · intro h
    have e := mut_assoc_error al be
    rw [h] at e
    have := congrArg M2.a e
    simp only [sub, smul] at this; grind
  · intro h
    have e := mut_assoc_error al be
    rw [h] at e
    apply ext'
    · have := congrArg M2.a e; simp only [sub, smul] at this; grind
    · have := congrArg M2.b e; simp only [sub, smul] at this; grind
    · have := congrArg M2.c e; simp only [sub, smul] at this; grind
    · have := congrArg M2.d e; simp only [sub, smul] at this; grind

/-! ## 4. Anticommuting forces "at most once" -/

open M2 in
/-- If swapping two copies of x costs a sign (x·x = −x·x), then x·x = 0. -/
theorem anticommute_self_zero (x : M2) (h : mul x x = neg (mul x x)) : mul x x = ⟨0, 0, 0, 0⟩ := by
  have h1 := congrArg M2.a h; have h2 := congrArg M2.b h
  have h3 := congrArg M2.c h; have h4 := congrArg M2.d h
  simp only [neg] at h1 h2 h3 h4
  apply ext' <;> simp only <;> omega

def PX : M2 := ⟨0, 1, 1, 0⟩
def PZ : M2 := ⟨1, 0, 0, -1⟩

/-- Pauli X and Z: each squares to 1 and they anticommute (order dial q = −1). -/
theorem pauli : M2.mul PX PX = ⟨1, 0, 0, 1⟩ ∧ M2.mul PZ PZ = ⟨1, 0, 0, 1⟩ ∧
    M2.mul PZ PX = M2.neg (M2.mul PX PZ) := by decide

/-! ## 5. The discrete Weyl relation -/

/-- Shift a sequence: (shift a f) n = f (n + a). -/
def shift (a : Nat) (f : Nat → Int) : Nat → Int := fun n => f (n + a)
/-- Multiply a sequence by λⁿ. -/
def mulPow (l : Int) (f : Nat → Int) : Nat → Int := fun n => l ^ n * f n

theorem pow_add' (l : Int) (n a : Nat) : l ^ (n + a) = l ^ n * l ^ a := by
  induction a with
  | zero => simp
  | succ a ih => rw [← Nat.add_assoc, Int.pow_succ, ih, Int.pow_succ]; grind

/-- Shifting after multiplying = λ^a times multiplying after shifting: the q-commutation yx = q·xy with q = λ^a. -/
theorem weyl (a : Nat) (l : Int) (f : Nat → Int) (n : Nat) :
    shift a (mulPow l f) n = l ^ a * mulPow l (shift a f) n := by
  simp only [shift, mulPow, pow_add']; grind

/-! ## 6. Kleisli composition and the k-sum -/

/-- The k-sum distributes only where k·c·(c − 1)·a·b = 0. -/
theorem ksum_distrib_gap (k a b c : Int) : circK k (c * a) (c * b) - c * circK k a b = k * a * b * c * (c - 1) := by
  simp only [circK]; grind

/-- At k = 0 the k-sum is plain addition (bags: count paths). -/
theorem ksum_zero (a b : Int) : circK 0 a b = a + b := by simp only [circK]; grind

/-- At k = −1 on 0/1 values the k-sum is OR (relations: is there a path). -/
theorem ksum_or (a b : Int) (ha : a = 0 ∨ a = 1) (hb : b = 0 ∨ b = 1) : circK (-1) a b = max a b := by
  rcases ha with rfl | rfl <;> rcases hb with rfl | rfl <;> decide

/-- Off 0/1 the k = −1 sum fails to distribute: c = 2, a = b = 1. -/
theorem ksum_not_distrib : circK (-1 : Int) (2 * 1) (2 * 1) ≠ 2 * circK (-1) 1 1 := by decide

/-! ## 7. Dictionary merge -/

/-- Merge two dictionaries; the right one wins on shared keys. -/
def over (m₁ m₂ : Nat → Option Nat) : Nat → Option Nat := fun key =>
  match m₂ key with
  | some v => some v
  | none => m₁ key

theorem over_assoc (x y z : Nat → Option Nat) : over (over x y) z = over x (over y z) := by
  funext key; simp only [over]; cases z key <;> cases y key <;> rfl

theorem over_idem (x : Nat → Option Nat) : over x x = x := by
  funext key; simp only [over]; cases x key <;> simp

/-- Right regular band: x y x = y x. -/
theorem over_right_regular (x y : Nat → Option Nat) : over (over x y) x = over y x := by
  funext key; simp only [over]; cases x key <;> cases y key <;> rfl

/-- Merge order matters on a shared key. -/
theorem over_not_comm :
    over (fun _ => some 1) (fun _ => some 2) 0 ≠ over (fun _ => some 2) (fun _ => some 1) 0 := by decide

/-! ## 8. The order dial on x x y y -/

/-- Out-of-order pairs (a larger letter before a smaller one). -/
def inv : List Nat → Nat
  | [] => 0
  | x :: xs => (xs.filter (fun y => decide (y < x))).length + inv xs

/-- The six arrangements of x x y y (x = 0, y = 1). -/
def xxyy : List (List Nat) := [[0,0,1,1],[0,1,0,1],[0,1,1,0],[1,0,0,1],[1,0,1,0],[1,1,0,0]]

/-- Their out-of-order counts give 1 + q + 2q² + q³ + q⁴ (the Gaussian binomial [4,2]_q). -/
theorem xxyy_inversions : xxyy.map inv = [0, 1, 2, 2, 3, 4] := by decide

/-! ## 9. No functor half power set (size skeleton) -/

/-- Sizes of a would-be functor H with H(H X) ≅ P(X): h n = size of H applied to an n-element set.
    A functor keeps split injections (a nonempty set's injections have retractions), so h is
    monotone from 1 up; and H(H n) has 2ⁿ elements. Already at n ≤ 2 this is impossible.
    (The case H(1) infinite is ruled out before this skeleton: H(1) injects into H(H 1), which has 2 elements.) -/
theorem no_functor_half (h : Nat → Nat) (mono : ∀ n m, 1 ≤ n → n ≤ m → h n ≤ h m)
    (hh : ∀ n, n ≤ 2 → h (h n) = 2 ^ n) : False := by
  have h0 := hh 0 (by omega); have h1 := hh 1 (by omega); have h2 := hh 2 (by omega)
  simp at h0 h1 h2
  by_cases a0 : h 1 = 0
  · rw [a0] at h1
    have : h 2 = 1 := by rw [h1] at h0; exact h0
    rw [this, a0] at h2; omega
  · have m := mono 1 (h 1) (Nat.le_refl 1) (by omega)
    rw [h1] at m
    have : h 1 = 1 ∨ h 1 = 2 := by omega
    rcases this with e | e
    · rw [e, e] at h1; omega
    · rw [e] at h1; rw [h1, h1] at h2; omega

end OrderDial
