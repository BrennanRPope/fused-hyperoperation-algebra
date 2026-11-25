/-
  Measure theory and probability through the circle operation (called mult-add before 1.0)
  ===============================================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`.  Intended module path:
  FusedAlgebra/Topics/MeasureProbability.lean (imports Base and Topics.Statistics).

  1. Sugeno λ-measures (Sugeno 1974).  A finite set of atoms with densities gᵢ;
     μ(A) is the dial fold of A's densities with weight k = λ.  Disjoint union is
     ◯ₖ, the lift 1 + λμ multiplies, and the Möbius masses are the corners of the
     circle-operation cube, λ^(|S|−1)·∏g_S.
  2. Inclusion–exclusion becomes multiplication; duality μ*(A) = 1 − μ(Aᶜ) is a
     λ*-measure with λ* = the circle inverse of λ; m equal atoms give a q-number.
  3. The wall: values below the wall −1/k stay below it, so for k ≤ −1 no set
     of such atoms can have measure 1 (this is why Sugeno needs λ > −1).
  4. Odds of an independent "or" combine by the circle operation; odds-against of an "and" combine by the circle operation.
  5. Expectation passes through the dial for independent variables; for
     dependent ones the gap is k·Cov.  Multiplicative martingales.
  6. Yor's formula, discrete: stochastic exponentials multiply step by step
     through Δx ◯ Δy.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/MeasureProbability.lean; merged unchanged; no exact twins.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Statistics

namespace MeasureProbability

open CircleOp
open StatFused (Ex Ex_add Ex_smul Ex_const)

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! # 1. The λ-measure of a list of atoms -/

/-- μ of a set of atoms with densities `l`, under weight k: the dial fold. -/
def lam (k : Rat) : List Rat → Rat
  | [] => 0
  | g :: l => circK k g (lam k l)

/-- Disjoint union: μ(A ∪ B) = μ(A) ◯ₖ μ(B). -/
theorem lam_append (k : Rat) (A B : List Rat) : lam k (A ++ B) = circK k (lam k A) (lam k B) := by
  induction A with
  | nil => simp only [List.nil_append, lam]; unfold circK; grind
  | cons g A ih => simp only [List.cons_append, lam]; rw [ih, circK_assoc]

/-- The lift: 1 + λ·μ(A) = ∏ (1 + λ·gᵢ). -/
theorem lam_lift (k : Rat) (l : List Rat) : 1 + k * lam k l = prodL (l.map (fun g => 1 + k * g)) := by
  induction l with
  | nil => simp only [lam, List.map, prodL]; grind
  | cons g l ih => simp only [lam, List.map, prodL]; rw [circK_lift, ih]

/-- λ·μ is the plain circle operation of the scaled densities λ·gᵢ. -/
theorem lam_scale (k : Rat) (l : List Rat) : k * lam k l = foldCirc (l.map (fun g => k * g)) := by
  induction l with
  | nil => simp only [lam, List.map, foldCirc]; grind
  | cons g l ih => simp only [lam, List.map, foldCirc]; rw [circK_scale, ih]

/-- Möbius masses are cube corners: 1 + λμ(A) = Σ over sub-lists S of ∏_{i∈S} λgᵢ. -/
theorem lam_corners (k : Rat) (l : List Rat) :
    1 + k * lam k l = sumL ((groups (l.map (fun g => k * g))).map prodL) := by
  rw [lam_scale, one_add_foldCirc_eq_corners]

/-- The three-atom expansion: layer 1, layer 2 weighted λ, layer 3 weighted λ². -/
theorem lam_three (k a b c : Rat) :
    lam k [a, b, c] = (a + b + c) + k * (a * b + a * c + b * c) + k * k * (a * b * c) := by
  simp only [lam]; unfold circK; grind

/-- λ = 0 is an ordinary additive measure. -/
theorem lam_additive (l : List Rat) : lam 0 l = sumL l := by
  induction l with
  | nil => rfl
  | cons g l ih => simp only [lam, sumL]; rw [circK_weight_zero, ih]

/-- λ = −1: μ(A) = 1 − ∏(1 − gᵢ), the chance at least one atom fires. -/
theorem lam_neg_one (l : List Rat) : lam (-1) l = 1 - prodL (l.map (fun g => 1 - g)) := by
  induction l with
  | nil => simp only [lam, List.map, prodL]; grind
  | cons g l ih => simp only [lam, List.map, prodL]; rw [ih]; unfold circK; grind

/-- Order of the atoms never matters. -/
theorem lam_comm_two (k : Rat) (A B : List Rat) : lam k (A ++ B) = lam k (B ++ A) := by
  rw [lam_append, lam_append, circK_comm]

/-! # 2. Normalization, inclusion–exclusion, duality, q-numbers -/

/-- Two atoms: λ = (1 − g₁ − g₂)/(g₁g₂) solves (1 + λg₁)(1 + λg₂) = 1 + λ. -/
theorem norm_two (g1 g2 : Rat) (h1 : g1 ≠ 0) (h2 : g2 ≠ 0) :
    (1 + (1 - g1 - g2) / (g1 * g2) * g1) * (1 + (1 - g1 - g2) / (g1 * g2) * g2)
      = 1 + (1 - g1 - g2) / (g1 * g2) := by
  have h : g1 * g2 ≠ 0 := by grind
  grind

/-- Any nonzero root of the two-atom normalization is that one. -/
theorem norm_two_root (k g1 g2 : Rat) (hk : k ≠ 0)
    (h : (1 + k * g1) * (1 + k * g2) = 1 + k) : k * (g1 * g2) = 1 - g1 - g2 := by
  have e : k * (k * (g1 * g2) - (1 - g1 - g2)) = 0 := by grind
  rcases Rat.mul_eq_zero.mp e with h0 | h0
  · exact absurd h0 hk
  · grind

/-- Two atoms with positive densities: λ > 0 exactly when g₁ + g₂ < 1. -/
theorem norm_two_sign (k g1 g2 : Rat) (hk : k ≠ 0) (h1 : 0 < g1) (h2 : 0 < g2)
    (h : (1 + k * g1) * (1 + k * g2) = 1 + k) : 0 < k ↔ g1 + g2 < 1 := by
  have r := norm_two_root k g1 g2 hk h
  have p : 0 < g1 * g2 := Rat.mul_pos h1 h2
  constructor
  · intro hk0; have := Rat.mul_pos hk0 p; grind
  · intro hs
    by_cases hlt : 0 < k
    · exact hlt
    · have : k * (g1 * g2) ≤ 0 := by
        have := Rat.mul_nonneg (show 0 ≤ -k by grind) (Rat.le_of_lt p); grind
      grind

/-- Inclusion–exclusion is multiplication: with atoms x (only A), y (only B),
    z (both), (1 + λμ(A∪B))(1 + λμ(A∩B)) = (1 + λμ(A))(1 + λμ(B)). -/
theorem lam_incl_excl (k : Rat) (X Y Z : List Rat) :
    (1 + k * lam k (X ++ Y ++ Z)) * (1 + k * lam k Z)
      = (1 + k * lam k (X ++ Z)) * (1 + k * lam k (Y ++ Z)) := by
  simp only [lam_append, circK_lift]; grind

/-- Duality: if the whole space a, b, c is normalized, the dual measure
    μ*(S) = 1 − μ(complement) obeys the dial with weight circInv λ. -/
theorem dual_weight (k a b c : Rat) (hk : 1 + k ≠ 0) (hk0 : k ≠ 0)
    (hn : (1 + k * a) * (1 + k * b) * (1 + k * c) = 1 + k) :
    1 - c = circK (circInv k) (1 - circK k b c) (1 - circK k a c) := by
  -- the normalization says N = 0 for the polynomial N below (k ≠ 0)
  have hN : k * (a*b*c*k*k + a*b*k + a*c*k + a + b*c*k + b + c - 1) = 0 := by grind
  have N0 : a*b*c*k*k + a*b*k + a*c*k + a + b*c*k + b + c - 1 = 0 := by
    rcases Rat.mul_eq_zero.mp hN with h | h
    · exact absurd h hk0
    · exact h
  have hkk : (1 + k) * circInv k = -k := by unfold circInv; grind
  have e : (1 + k) * circK (circInv k) (1 - circK k b c) (1 - circK k a c)
      = (1 + k) * ((1 - circK k b c) + (1 - circK k a c))
        + ((1 + k) * circInv k) * ((1 - circK k b c) * (1 - circK k a c)) := by
    unfold circK; grind
  have E0 : (1 + k) * (1 - c) - (1 + k) * ((1 - circK k b c) + (1 - circK k a c))
      + k * ((1 - circK k b c) * (1 - circK k a c))
      = (c * k + 1) * (a*b*c*k*k + a*b*k + a*c*k + a + b*c*k + b + c - 1) := by
    unfold circK; grind
  rw [N0] at E0
  have := Rat.mul_eq_zero.mp (show (1 + k) * ((1 - c) - circK (circInv k) (1 - circK k b c) (1 - circK k a c)) = 0 by
    have d : (1 + k) * ((1 - c) - circK (circInv k) (1 - circK k b c) (1 - circK k a c))
        = (1 + k) * (1 - c) - (1 + k) * circK (circInv k) (1 - circK k b c) (1 - circK k a c) := by grind
    rw [d, e, hkk]; generalize circK k b c = u at E0 ⊢; generalize circK k a c = v at E0 ⊢; grind)
  rcases this with h0 | h0
  · exact absurd h0 hk
  · grind

/-- The dual of the dual weight is the weight again. -/
theorem dual_involution (k : Rat) (hk : 1 + k ≠ 0) : circInv (circInv k) = k := by
  unfold circInv
  have e : 1 + -k / (1 + k) = 1 / (1 + k) := by grind
  rw [e]; grind

/-- In q = 1 + λ terms the dual is the reciprocal: q* = 1/q. -/
theorem dual_q (k : Rat) (hk : 1 + k ≠ 0) : 1 + circInv k = 1 / (1 + k) := by
  unfold circInv; grind

/-- q-number [m]_q = 1 + q + … + q^(m−1). -/
def qnum (q : Rat) : Nat → Rat
  | 0 => 0
  | m + 1 => qnum q m + pw q m

/-- m equal atoms of size g have measure g·[m]_q with q = 1 + λg. -/
theorem lam_equal_atoms (k g : Rat) (m : Nat) :
    lam k (List.replicate m g) = g * qnum (1 + k * g) m := by
  have key : ∀ m, (1 + k * g) * pw (1 + k * g) m = pw (1 + k * g) m * (1 + k * g) := by
    intro m; grind
  induction m with
  | zero => simp only [List.replicate, lam, qnum]; grind
  | succ m ih =>
    simp only [List.replicate, lam]; rw [ih]
    unfold circK
    have : qnum (1 + k * g) (m + 1) = qnum (1 + k * g) m + pw (1 + k * g) m := rfl
    rw [this]
    have hq : ∀ n, g * k * qnum (1 + k * g) n + 1 = pw (1 + k * g) n := by
      intro n; induction n with
      | zero => simp only [qnum, pw]; grind
      | succ n ihn => simp only [qnum, pw]; grind
    have := hq m; grind

/-! # 3. The wall must stay outside [0, 1] -/

/-- Below the wall stays below the wall: positive lifts multiply to a positive lift. -/
theorem wall_blocks (k a b : Rat) (ha : 0 < 1 + k * a) (hb : 0 < 1 + k * b) :
    0 < 1 + k * circK k a b := by
  rw [circK_lift]; exact Rat.mul_pos ha hb

/-- For any list of atoms below the wall, the measure stays below the wall. -/
theorem wall_list (k : Rat) (l : List Rat) (h : ∀ g ∈ l, 0 < 1 + k * g) :
    0 < 1 + k * lam k l := by
  induction l with
  | nil => simp only [lam]; grind
  | cons g l ih =>
    simp only [lam]
    exact wall_blocks k g (lam k l) (h g (by simp)) (ih (fun x hx => h x (by simp [hx])))

/-- So when the wall sits in (0, 1] (k ≤ −1), no such set has measure 1:
    a normalized λ-measure needs λ > −1. -/
theorem no_norm_past_wall (k : Rat) (hk : k ≤ -1) (l : List Rat)
    (h : ∀ g ∈ l, 0 < 1 + k * g) : lam k l ≠ 1 := by
  intro e
  have := wall_list k l h
  rw [e] at this; grind

/-- At k = −1 the only way to reach 1 is a sure atom: the product of misses vanishes. -/
theorem sure_atom (l : List Rat) (h : lam (-1) l = 1) : prodL (l.map (fun g => 1 - g)) = 0 := by
  rw [lam_neg_one] at h; grind

theorem prodL_eq_zero (l : List Rat) (h : prodL l = 0) : (0 : Rat) ∈ l := by
  induction l with
  | nil => simp only [prodL] at h; grind
  | cons x l ih =>
    simp only [prodL] at h
    rcases Rat.mul_eq_zero.mp h with h0 | h0
    · simp [h0]
    · simp [ih h0]

/-- k = −1, measure 1 ⇒ some atom has density exactly 1. -/
theorem sure_atom_exists (l : List Rat) (h : lam (-1) l = 1) : ∃ g ∈ l, g = 1 := by
  have h0 := prodL_eq_zero _ (sure_atom l h)
  rcases List.mem_map.mp h0 with ⟨g, hg, e⟩
  exact ⟨g, hg, by grind⟩

/-! # 4. Odds -/

def odds (p : Rat) : Rat := p / (1 - p)
def against (p : Rat) : Rat := (1 - p) / p

/-- The lift of the odds is 1/(1 − p). -/
theorem odds_lift (p : Rat) (hp : p ≠ 1) : 1 + odds p = 1 / (1 - p) := by
  unfold odds; have : (1 : Rat) - p ≠ 0 := by grind
  grind

/-- Odds of an independent "or" combine by the circle operation: odds(p ◯₋₁ q) = odds p ◯ odds q. -/
theorem odds_or (p q : Rat) (hp : p ≠ 1) (hq : q ≠ 1) :
    odds (circK (-1) p q) = odds p ◯ odds q := by
  unfold odds circ circK
  have h1 : (1 : Rat) - p ≠ 0 := by grind
  have h2 : (1 : Rat) - q ≠ 0 := by grind
  have e : 1 - (p + q + -1 * p * q) = (1 - p) * (1 - q) := by grind
  rw [e]
  have h3 : (1 - p) * (1 - q) ≠ 0 := by grind
  grind

/-- Odds against an independent "and" combine by the circle operation. -/
theorem against_and (p q : Rat) (hp : p ≠ 0) (hq : q ≠ 0) :
    against (p * q) = against p ◯ against q := by
  unfold against circ
  have h3 : p * q ≠ 0 := by grind
  grind

/-- Three independent events: the odds of "any" is the three-way circle operation. -/
theorem odds_or_three (p q r : Rat) (hp : p ≠ 1) (hq : q ≠ 1) (hr : r ≠ 1) :
    odds (circK (-1) (circK (-1) p q) r) = odds p ◯ odds q ◯ odds r := by
  have hpq : circK (-1) p q ≠ 1 := by
    unfold circK; intro e
    have : (1 - p) * (1 - q) = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h | h <;> grind
  rw [odds_or _ _ hpq hr, odds_or _ _ hp hq]

/-! # 5. Expectation and independence -/

/-- The product of two finite distributions (independence). -/
def prodD {α β : Type} : List (α × Rat) → List (β × Rat) → List ((α × β) × Rat)
  | [], _ => []
  | (a, w) :: ws, vs => vs.map (fun bv => ((a, bv.1), w * bv.2)) ++ prodD ws vs

def mass {α : Type} : List (α × Rat) → Rat
  | [] => 0
  | (_, w) :: ws => w + mass ws

theorem Ex_append {α : Type} (l m : List (α × Rat)) (f : α → Rat) :
    Ex (l ++ m) f = Ex l f + Ex m f := by
  induction l with
  | nil => simp only [List.nil_append, Ex]; grind
  | cons p l ih => obtain ⟨a, w⟩ := p; simp only [List.cons_append, Ex]; rw [ih]; grind

theorem Ex_row {α β : Type} (a : α) (w : Rat) (vs : List (β × Rat)) (F : α × β → Rat) :
    Ex (vs.map (fun bv => ((a, bv.1), w * bv.2))) F = w * Ex vs (fun b => F (a, b)) := by
  induction vs with
  | nil => simp only [List.map, Ex]; grind
  | cons p vs ih => obtain ⟨b, v⟩ := p; simp only [List.map, Ex] at ih ⊢; rw [ih]; grind

/-- Independent product: E[f(X)·g(Y)] = E f(X) · E g(Y). -/
theorem Ex_prod {α β : Type} (ws : List (α × Rat)) (vs : List (β × Rat)) (f : α → Rat) (g : β → Rat) :
    Ex (prodD ws vs) (fun ab => f ab.1 * g ab.2) = Ex ws f * Ex vs g := by
  induction ws with
  | nil => simp only [prodD, Ex]; grind
  | cons p ws ih =>
    obtain ⟨a, w⟩ := p
    simp only [prodD, Ex]; rw [Ex_append, Ex_row, ih]
    have : ∀ (vs : List (β × Rat)), Ex vs (fun b => f a * g b) = f a * Ex vs g := by
      intro vs; induction vs with
      | nil => simp only [Ex]; grind
      | cons p vs ih2 => obtain ⟨b, v⟩ := p; simp only [Ex]; rw [ih2]; grind
    rw [this]; grind

theorem Ex_one {α : Type} (ws : List (α × Rat)) : Ex ws (fun _ => 1) = mass ws := by
  induction ws with
  | nil => rfl
  | cons p ws ih => obtain ⟨a, w⟩ := p; simp only [Ex, mass]; rw [ih]; grind

/-- The engine: for independent X and Y (total mass 1 each),
    E[X ◯ₖ Y] = E X ◯ₖ E Y. -/
theorem indep_dial {α β : Type} (k : Rat) (ws : List (α × Rat)) (vs : List (β × Rat))
    (hw : mass ws = 1) (hv : mass vs = 1) (X : α → Rat) (Y : β → Rat) :
    Ex (prodD ws vs) (fun ab => circK k (X ab.1) (Y ab.2)) = circK k (Ex ws X) (Ex vs Y) := by
  have e1 := Ex_prod ws vs X (fun _ => 1)
  have e2 := Ex_prod ws vs (fun _ => 1) Y
  have e3 := Ex_prod ws vs X Y
  rw [Ex_one, hv] at e1; rw [Ex_one, hw] at e2
  have split : Ex (prodD ws vs) (fun ab => circK k (X ab.1) (Y ab.2))
      = Ex (prodD ws vs) (fun ab => X ab.1 * 1) + Ex (prodD ws vs) (fun ab => 1 * Y ab.2)
        + k * Ex (prodD ws vs) (fun ab => X ab.1 * Y ab.2) := by
    unfold circK
    rw [← Ex_smul, ← Ex_add, ← Ex_add]
    apply StatFused.Ex_congr; intro ab; grind
  rw [split, e1, e2, e3]; unfold circK; grind

/-- For any joint law (dependent or not), the gap is k times the covariance:
    the corner's excess. -/
theorem dial_gap_cov {α : Type} (k : Rat) (ws : List (α × Rat)) (X Y : α → Rat) :
    Ex ws (fun a => circK k (X a) (Y a)) - circK k (Ex ws X) (Ex ws Y)
      = k * (Ex ws (fun a => X a * Y a) - Ex ws X * Ex ws Y) := by
  have : Ex ws (fun a => circK k (X a) (Y a))
      = Ex ws X + Ex ws Y + k * Ex ws (fun a => X a * Y a) := by
    unfold circK
    rw [← Ex_smul, ← Ex_add, ← Ex_add]
    apply StatFused.Ex_congr; intro a; grind
  rw [this]; unfold circK; grind

/-- Multiplicative martingale step: independent mean-zero returns keep mean zero. -/
theorem martingale_step {α β : Type} (ws : List (α × Rat)) (vs : List (β × Rat))
    (hw : mass ws = 1) (hv : mass vs = 1) (X : α → Rat) (Y : β → Rat)
    (h0 : Ex ws X = 0) (h1 : Ex vs Y = 0) :
    Ex (prodD ws vs) (fun ab => X ab.1 ◯ Y ab.2) = 0 := by
  have := indep_dial 1 ws vs hw hv X Y
  simp only [circK_one] at this
  rw [this, h0, h1]; unfold circ; grind

/-! # 6. Yor's formula, discrete -/

/-- The discrete stochastic exponential of increments: ℰ = ∏(1 + Δx) − 1. -/
def sexp (l : List Rat) : Rat := foldCirc l

/-- Yor's formula step by step: ℰ(x) ◯ ℰ(y) = ℰ(x ◯ y) increment-wise,
    where Δx ◯ Δy = Δx + Δy + ΔxΔy carries the bracket [X, Y]. -/
theorem yor (xs ys : List Rat) (h : xs.length = ys.length) :
    sexp xs ◯ sexp ys = sexp (List.zipWith circ xs ys) := by
  induction xs generalizing ys with
  | nil =>
    cases ys with
    | nil => simp only [sexp, foldCirc, List.zipWith]; rw [zero_circ]
    | cons _ _ => simp at h
  | cons x xs ih =>
    cases ys with
    | nil => simp at h
    | cons y ys =>
      simp only [List.length_cons] at h
      simp only [sexp, foldCirc, List.zipWith] at ih ⊢
      rw [← ih ys (by omega)]
      unfold circ; grind

/-- One step: the increment of X + Y + [X, Y] is exactly Δx ◯ Δy. -/
theorem yor_step (dx dy : Rat) : (1 + dx) * (1 + dy) = 1 + (dx + dy + dx * dy) := by grind

/-- The self-bracket: one step of ℰ(x) against its own circle inverse is flat. -/
theorem sexp_inverse (dx : Rat) (h : 1 + dx ≠ 0) : dx ◯ circInv dx = 0 :=
  circ_circInv dx h

end MeasureProbability
