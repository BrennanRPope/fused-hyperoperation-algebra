/-
  Fused Hyperoperation Algebra: chaos as repeated self-composition.
  Core Lean 4.37 only (no mathlib), no `sorry`.
  Intended library path: FusedAlgebra/Topics/ChaosComposition.lean (imports Dynamics and OrderDial).

  1. Iterating one map is a functor from (ℕ, +) to (maps, ∘); a map h with h ∘ f = g ∘ h
     (a morphism of dynamical systems) carries orbits to orbits.
  2. Every dial square x ◯ₖ x has the dynamics of x ◯ x, orbit for orbit, through x ↦ kx.
  3. The logistic map at r = 4 is the fold of circle-operation squaring, at every step of the orbit
     (Ulam and von Neumann 1947; the fold through the circle inverse is the repo's reading).
  4. Order is invisible to iteration: g ∘ (f ∘ g)ⁿ = (g ∘ f)ⁿ ∘ g, so f ∘ g and g ∘ f have the same
     orbits and periodic points (shift equivalence, Williams 1973). For affine maps the second
     map carries one order's centre onto the other's.
  5. Order is what makes chaos: two branches that share a centre (they commute) give composites
     that depend only on the bag of letters; two branches with different centres, x ↦ 2x and
     x ↦ 2x + 1, give a different composite for every word (the list corner).
  6. Doubling is dropping the leading binary digit: the shift on digit strings (the list `tail`).
-/
import FusedAlgebra.Topics.Dynamics
import FusedAlgebra.Topics.OrderDial

namespace ChaosComposition
open CircleOp
open Dynamics (iter msq fold)
open Lean.Grind

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. Iteration is a functor; morphisms carry orbits to orbits -/

/-- fᵐ⁺ⁿ = fᵐ ∘ fⁿ: iterating sends + in ℕ to ∘ on maps. -/
theorem iter_add {α : Type} (f : α → α) (m n : Nat) (x : α) :
    iter f (m + n) x = iter f m (iter f n x) := by
  induction m with
  | zero => rw [Nat.zero_add]; rfl
  | succ m ih =>
    rw [Nat.succ_add]
    show f (iter f (m + n) x) = f (iter f m (iter f n x))
    rw [ih]

/-- fⁿ⁺¹ = fⁿ ∘ f as well as f ∘ fⁿ. -/
theorem iter_succ' {α : Type} (f : α → α) (n : Nat) (x : α) :
    iter f (n + 1) x = iter f n (f x) := by
  induction n with
  | zero => rfl
  | succ n ih => show f (iter f (n + 1) x) = f (iter f n (f x)); rw [ih]

/-- A morphism of dynamical systems (h ∘ f = g ∘ h, a natural transformation between the two
    functors out of ℕ) carries every orbit of f onto an orbit of g. -/
theorem iter_semiconj {α β : Type} (f : α → α) (g : β → β) (h : α → β)
    (hc : ∀ x, h (f x) = g (h x)) (n : Nat) (x : α) : h (iter f n x) = iter g n (h x) := by
  induction n with
  | zero => rfl
  | succ n ih => show h (f (iter f n x)) = g (iter g n (h x)); rw [hc, ih]

/-- So it carries periodic points to periodic points of the same period. -/
theorem periodic_semiconj {α β : Type} (f : α → α) (g : β → β) (h : α → β)
    (hc : ∀ x, h (f x) = g (h x)) (p : Nat) (x : α) (hx : iter f p x = x) :
    iter g p (h x) = h x := by
  rw [← iter_semiconj f g h hc, hx]

/-! ## 2. The dial does not change the dynamics -/

section Ring
variable {R : Type} [CommRing R]

/-- Orbit for orbit, x ◯ₖ x is x ◯ x in the units x ↦ kx. -/
theorem dial_iter_units (k x : R) (n : Nat) :
    k * iter (fun y => circK k y y) n x = iter msq n (k * x) :=
  iter_semiconj (fun y => circK k y y) msq (fun y => k * y)
    (fun y => Dynamics.dial_square_units k y) n x

/-- The quadratic family z² + c in lift coordinates: square by the circle operation, then shift by c. -/
theorem quadratic_lift (y c : R) : 1 + (msq y + c) = (1 + y) * (1 + y) + c := by
  unfold msq circ; grind

/-- Lyapunov through the lift, n steps at once: the slope of the n-th iterate of x ◯ x is
    2ⁿ (1 + x)^(2ⁿ − 1), written without subtraction as slope·(1 + x) = 2ⁿ (1 + x)^(2ⁿ). -/
theorem iter_slope (p : R × R) (n : Nat) :
    (iter Dynamics.msqJet n p).2 * (1 + p.1) = 2 ^ n * p.2 * (1 + p.1) ^ (2 ^ n) := by
  rw [Dynamics.iter_jet_ruler, Dynamics.iter_jet_fst, Dynamics.iter_msq_lift]

end Ring

/-! ## 3. The logistic map at r = 4 is folded squaring, along the whole orbit -/

section Field
variable {K : Type} [Field K]

/-- The logistic map at r = 4. -/
def logistic4 (x : K) : K := 4 * x * (1 - x)

/-- Squaring never reaches the wall from off it. -/
theorem msq_ne_wall (w : K) (h : 1 + w ≠ 0) : 1 + msq w ≠ 0 := by
  rw [Dynamics.msq_lift]; intro e
  have : (1 + w) * (1 + w) * ((1 + w)⁻¹ * (1 + w)⁻¹) = 0 := by rw [e]; grind
  grind

theorem iter_msq_ne_wall (w : K) (h : 1 + w ≠ 0) (n : Nat) : 1 + iter msq n w ≠ 0 := by
  induction n with
  | zero => exact h
  | succ n ih => exact msq_ne_wall _ ih

/-- n steps of the logistic map at r = 4 are n steps of x ◯ x, folded. -/
theorem logistic_iter_fold (w : K) (h : 1 + w ≠ 0) (n : Nat) :
    fold (iter msq n w) = iter logistic4 n (fold w) := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show fold (msq (iter msq n w)) = logistic4 (iter logistic4 n (fold w))
    rw [← ih, ← Dynamics.logistic_fold _ (iter_msq_ne_wall w h n)]
    rfl

end Field

/-! ## 4. Order is invisible to iteration -/

/-- Shift equivalence: g ∘ (f ∘ g)ⁿ = (g ∘ f)ⁿ ∘ g. The two orders are the same system, seen
    through g. -/
theorem order_invisible {α β : Type} (f : β → α) (g : α → β) (n : Nat) (x : α) :
    g (iter (fun y => f (g y)) n x) = iter (fun z => g (f z)) n (g x) :=
  iter_semiconj _ _ g (fun _ => rfl) n x

/-- A cycle of f ∘ g is carried by g to a cycle of g ∘ f with the same period. -/
theorem order_periodic {α β : Type} (f : β → α) (g : α → β) (p : Nat) (x : α)
    (hx : iter (fun y => f (g y)) p x = x) : iter (fun z => g (f z)) p (g x) = g x := by
  rw [← order_invisible, hx]

/-- …and back again by f, so periodic points of the two orders match up. -/
theorem order_periodic_back {α β : Type} (f : β → α) (g : α → β) (p : Nat) (z : β)
    (hz : iter (fun z => g (f z)) p z = z) : iter (fun y => f (g y)) p (f z) = f z :=
  order_periodic g f p z hz

/-- Affine maps: doing q carries the centre of p ∘ q onto the centre of q ∘ p. -/
theorem affine_centre_transfer (p q : Int × Int) (x : Int)
    (hx : OrderDial.app (OrderDial.comp p q) x = x) :
    OrderDial.app (OrderDial.comp q p) (OrderDial.app q x) = OrderDial.app q x := by
  rw [OrderDial.app_comp]; rw [OrderDial.app_comp] at hx; rw [hx]

/-- The two orders stretch by the same factor (the multiplier, hence the Lyapunov exponent, is
    order-free); only the shift, and so the centre, feels the order, by the wedge. -/
theorem affine_orders (p q : Int × Int) :
    (OrderDial.comp p q).1 = (OrderDial.comp q p).1 ∧
    (OrderDial.comp p q).2 - (OrderDial.comp q p).2 = p.1 * q.2 - q.1 * p.2 :=
  ⟨OrderDial.stretch_order_free p q, OrderDial.order_gap p q⟩

/-- The cost of order is one shift: p ∘ q = (shift by the wedge) ∘ (q ∘ p). In the affine group
    every commutator is a pure translation by the wedge, which slopes (and so Lyapunov exponents)
    cannot see: Lyapunov reads the abelian shadow, the wedge is the rest. -/
theorem order_is_a_shift (p q : Int × Int) :
    OrderDial.comp p q = OrderDial.comp (0, p.1 * q.2 - q.1 * p.2) (OrderDial.comp q p) := by
  simp only [OrderDial.comp, circ]
  apply Prod.ext <;> simp only <;> grind

/-! ## 5. Order is what makes chaos -/

/-- A 0/1 digit. -/
def bit (b : Bool) : Nat := if b then 1 else 0

/-- Run a word of the two branches x ↦ 2x + b (b = 0, 1), the first letter applied last.
    Their centres are 0 and −1, so they do not commute. -/
def runB : List Bool → Nat → Nat
  | [], x => x
  | b :: w, x => 2 * runB w x + bit b

/-- Every composite is x ↦ 2ⁿx + (its value at 0). -/
theorem runB_affine : ∀ (w : List Bool) (x : Nat), runB w x = 2 ^ w.length * x + runB w 0
  | [], x => by simp [runB]
  | b :: w, x => by
    simp only [runB, List.length_cons]
    rw [runB_affine w x, Nat.pow_succ]
    grind

/-- The list corner: different words of the same length give different composites (the branch
    monoid is free), so there are 2ⁿ composites of length n: growth rate log 2. -/
theorem runB_free : ∀ (w v : List Bool), w.length = v.length → runB w 0 = runB v 0 → w = v
  | [], [], _, _ => rfl
  | [], _ :: _, h, _ => by simp at h
  | _ :: _, [], h, _ => by simp at h
  | b :: w, c :: v, h, e => by
    simp only [runB] at e
    have hl : w.length = v.length := by simpa using h
    have hb : b = c := by cases b <;> cases c <;> simp [bit] at e ⊢ <;> omega
    subst hb
    have : runB w 0 = runB v 0 := by cases b <;> simp [bit] at e <;> omega
    rw [runB_free w v hl this]

/-- Two branches that share a centre: shifts x ↦ x + b (the k = 0 slice, centre at infinity). -/
def runT : List Bool → Nat → Nat
  | [], x => x
  | b :: w, x => runT w x + bit b

/-- The bag corner: a composite only counts its 1s, so length n gives just n + 1 composites:
    growth rate 0. -/
theorem runT_count : ∀ (w : List Bool) (x : Nat), runT w x = x + w.count true
  | [], x => by simp [runT]
  | b :: w, x => by
    cases b <;> simp [runT, bit, runT_count w x] <;> omega

theorem runT_perm {w v : List Bool} (h : List.Perm w v) (x : Nat) : runT w x = runT v x := by
  rw [runT_count, runT_count, h.count_eq]

section Stretch
variable {R : Type} [CommRing R]

/-- A word of stretches x ↦ m + (1 + a)(x − m), all about one centre m. -/
def stretchRun (m : R) : List R → R → R
  | [], x => x
  | a :: l, x => m + (1 + a) * (stretchRun m l x - m)

/-- The composite stretches about m by 1 + (the circle-operation fold of the letters). -/
theorem stretchRun_eq (m : R) : ∀ (l : List R) (x : R),
    stretchRun m l x = m + (1 + foldCirc l) * (x - m)
  | [], x => by simp only [stretchRun, foldCirc]; grind
  | a :: l, x => by
    simp only [stretchRun, foldCirc]; rw [stretchRun_eq m l x]; unfold circ; grind

/-- Shared centre ⇒ bag corner: the composite depends only on the bag of letters. -/
theorem stretchRun_perm (m : R) {l l' : List R} (h : List.Perm l l') (x : R) :
    stretchRun m l x = stretchRun m l' x := by
  rw [stretchRun_eq, stretchRun_eq, foldCirc_perm h]

end Stretch

/-! ## 6. Doubling is the shift on digit strings -/

/-- On N-digit binary fractions x/2ᴺ, doubling mod 1 moves every digit up one place: digit j + 1
    of D(x) is digit j of x (digits counted from the least significant). Read from the most
    significant end, D drops the first digit: the list `tail`. -/
theorem doubling_shift (N j x : Nat) (hj : j + 1 < N) :
    Nat.testBit ((x <<< 1) % 2 ^ N) (j + 1) = Nat.testBit x j := by
  rw [Nat.testBit_mod_two_pow, Nat.testBit_shiftLeft]
  simp [hj]

end ChaosComposition
