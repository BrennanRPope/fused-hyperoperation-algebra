/-
  Category theory through mult-add
  ================================

  Companion page and checks: ../index.html, ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports only FusedAlgebra.Base.
  Intended module path: FusedAlgebra/Topics/CategoryTheory.lean.

  A ✳ B := A + B + A×B is the type `These A B` (Haskell's `These`, Rust's `EitherOrBoth`).

    1. Option turns ✳ into ×: Option (These A B) ≃ Option A × Option B, naturally.
       Unit (These Empty B ≃ B), symmetry and associativity: (FinSet, ✳, ∅) is symmetric monoidal.
       The dial with k tags: 1 + K×(A ✳ₖ B) ≃ (1 + K×A) × (1 + K×B).
    2. ✳ is the categorical product in Par, the category of sets and partial functions
       (functions X → Option Y): projections, pairing, and uniqueness of the pairing.
       Empty is a zero object there.
    3. Counting: the list of all values of These has length a ✳ b; pointed product
       = wedge (a + b) + smash (ab).
    4. align / zipLongest: pairing two lists into These loses nothing; both inputs come back,
       the length is the max, and the zip is the `both` part.
    5. Decategorification: cell counts of CW complexes multiply under products (as polynomials
       in t); the reduced count P − 1 mult-adds, at every t. At t = −1 this is the reduced Euler
       characteristic; the circle (χ̃ = −1) is the wall, spheres sit at ±1, suspension multiplies by t.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/CategoryTheory.lean; merged unchanged; no exact twins.
-/
import FusedAlgebra.Base

namespace CategoryTheory
open MultAdd

local infixl:65 " ✳ " => multAdd

set_option linter.unusedSectionVars false

/-! ## 1. The type A ✳ B and Option -/

/-- A ✳ B = A + B + A×B: left only, right only, or both. -/
inductive These (α β : Type) where
  | this : α → These α β
  | that : β → These α β
  | both : α → β → These α β

open These

/-- A type isomorphism: maps both ways that undo each other. -/
structure Iso (α β : Type) where
  to : α → β
  inv : β → α
  left : ∀ x, inv (to x) = x
  right : ∀ y, to (inv y) = y

variable {α β γ δ : Type}

/-- Option of These, split into a pair of Options. -/
def toPair : Option (These α β) → Option α × Option β
  | none => (none, none)
  | some (this a) => (some a, none)
  | some (that b) => (none, some b)
  | some (both a b) => (some a, some b)

/-- A pair of Options, fused back into an Option of These. -/
def ofPair : Option α × Option β → Option (These α β)
  | (none, none) => none
  | (some a, none) => some (this a)
  | (none, some b) => some (that b)
  | (some a, some b) => some (both a b)

theorem ofPair_toPair (x : Option (These α β)) : ofPair (toPair x) = x := by
  rcases x with _ | (a | b | ⟨a, b⟩) <;> rfl

theorem toPair_ofPair (p : Option α × Option β) : toPair (ofPair p) = p := by
  rcases p with ⟨_ | a, _ | b⟩ <;> rfl

/-- **Option turns ✳ into ×**: Option (A ✳ B) ≅ Option A × Option B. -/
def optionThese : Iso (Option (These α β)) (Option α × Option β) :=
  ⟨toPair, ofPair, ofPair_toPair, toPair_ofPair⟩

/-- Option of the empty type is the one-point type: the unit goes to the unit. -/
theorem option_empty (x : Option Empty) : x = none := by
  rcases x with _ | e
  · rfl
  · exact nomatch e

/-- Functoriality: These is a functor in both slots. -/
def These.map (f : α → γ) (g : β → δ) : These α β → These γ δ
  | this a => this (f a)
  | that b => that (g b)
  | both a b => both (f a) (g b)

theorem map_id (t : These α β) : t.map id id = t := by
  rcases t with a | b | ⟨a, b⟩ <;> rfl

theorem map_comp {ε ζ : Type} (f : α → γ) (g : β → δ) (f' : γ → ε) (g' : δ → ζ) (t : These α β) :
    (t.map f g).map f' g' = t.map (f' ∘ f) (g' ∘ g) := by
  rcases t with a | b | ⟨a, b⟩ <;> rfl

/-- Naturality: the split commutes with mapping. -/
theorem toPair_natural (f : α → γ) (g : β → δ) (x : Option (These α β)) :
    toPair (x.map (These.map f g)) = ((toPair x).1.map f, (toPair x).2.map g) := by
  rcases x with _ | (a | b | ⟨a, b⟩) <;> rfl

/-- Unit: ∅ ✳ B ≅ B. -/
def unitL : Iso (These Empty β) β where
  to := fun t => match t with
    | this e => nomatch e
    | that b => b
    | both e _ => nomatch e
  inv := that
  left := by intro t; rcases t with e | b | ⟨e, _⟩ <;> first | exact nomatch e | rfl
  right := by intro b; rfl

/-- Symmetry: A ✳ B ≅ B ✳ A. -/
def swap : These α β → These β α
  | this a => that a
  | that b => this b
  | both a b => both b a

theorem swap_swap (t : These α β) : swap (swap t) = t := by
  rcases t with a | b | ⟨a, b⟩ <;> rfl

/-- Associativity, left to right: (A ✳ B) ✳ C → A ✳ (B ✳ C). Seven corners each side. -/
def assocTo : These (These α β) γ → These α (These β γ)
  | this (this a) => this a
  | this (that b) => that (this b)
  | this (both a b) => both a (this b)
  | that c => that (that c)
  | both (this a) c => both a (that c)
  | both (that b) c => that (both b c)
  | both (both a b) c => both a (both b c)

def assocInv : These α (These β γ) → These (These α β) γ
  | this a => this (this a)
  | that (this b) => this (that b)
  | that (that c) => that c
  | that (both b c) => both (that b) c
  | both a (this b) => this (both a b)
  | both a (that c) => both (this a) c
  | both a (both b c) => both (both a b) c

/-- **Associativity**: (A ✳ B) ✳ C ≅ A ✳ (B ✳ C), the mult-add cube's seven corners. -/
def assoc : Iso (These (These α β) γ) (These α (These β γ)) where
  to := assocTo
  inv := assocInv
  left := by intro t; rcases t with (a | b | ⟨a, b⟩) | c | ⟨(a | b | ⟨a, b⟩), c⟩ <;> rfl
  right := by intro t; rcases t with a | (b | c | ⟨b, c⟩) | ⟨a, (b | c | ⟨b, c⟩)⟩ <;> rfl


/-! ## 1b. The dial with k tags: A ✳ₖ B = A + B + K×A×B -/

/-- The dial for a set K of k tags: the `both` corner carries a tag. -/
inductive TheseK (κ α β : Type) where
  | this : α → TheseK κ α β
  | that : β → TheseK κ α β
  | both : κ → α → β → TheseK κ α β

/-- The tagged lift A ↦ 1 + K×A (an Option with a tag) splits the tagged dial into a pair. -/
def toPairK {κ : Type} : Option (κ × TheseK κ α β) → Option (κ × α) × Option (κ × β)
  | none => (none, none)
  | some (i, .this a) => (some (i, a), none)
  | some (j, .that b) => (none, some (j, b))
  | some (i, .both j a b) => (some (i, a), some (j, b))

def ofPairK {κ : Type} : Option (κ × α) × Option (κ × β) → Option (κ × TheseK κ α β)
  | (none, none) => none
  | (some (i, a), none) => some (i, .this a)
  | (none, some (j, b)) => some (j, .that b)
  | (some (i, a), some (j, b)) => some (i, .both j a b)

theorem ofPairK_toPairK {κ : Type} (x : Option (κ × TheseK κ α β)) : ofPairK (toPairK x) = x := by
  rcases x with _ | ⟨i, a | b | ⟨j, a, b⟩⟩ <;> rfl

theorem toPairK_ofPairK {κ : Type} (p : Option (κ × α) × Option (κ × β)) : toPairK (ofPairK p) = p := by
  rcases p with ⟨_ | ⟨i, a⟩, _ | ⟨j, b⟩⟩ <;> rfl

/-- **The dial categorified**: 1 + K(A ✳ₖ B) ≅ (1 + KA) × (1 + KB), the lift 1 + k·x on sets. -/
def optionTheseK {κ : Type} : Iso (Option (κ × TheseK κ α β)) (Option (κ × α) × Option (κ × β)) :=
  ⟨toPairK, ofPairK, ofPairK_toPairK, toPairK_ofPairK⟩

/- The count 1 + k(a + b + k·ab) = (1 + ka)(1 + kb) is the library's `wAdd_lift`, on Nat. -/

/-! ## 2. ✳ is the product in Par (sets and partial functions) -/

/-- Composition of partial functions (Kleisli composition for Option). -/
def kc {X Y Z : Type} (h : X → Option Y) (k : Y → Option Z) : X → Option Z := fun x => (h x).bind k

/-- The two projections out of A ✳ B: partial, defined where that side is present. -/
def fst : These α β → Option α
  | this a => some a
  | that _ => none
  | both a _ => some a

def snd : These α β → Option β
  | this _ => none
  | that b => some b
  | both _ b => some b

/-- Pairing two partial functions. -/
def pair {X : Type} (f : X → Option α) (g : X → Option β) : X → Option (These α β) :=
  fun x => ofPair (f x, g x)

theorem fst_pair {X : Type} (f : X → Option α) (g : X → Option β) :
    kc (pair f g) fst = f := by
  funext x; unfold kc pair
  cases f x <;> cases g x <;> rfl

theorem snd_pair {X : Type} (f : X → Option α) (g : X → Option β) :
    kc (pair f g) snd = g := by
  funext x; unfold kc pair
  cases f x <;> cases g x <;> rfl

/-- **Uniqueness**: any partial map whose projections are f and g is the pairing.
    So A ✳ B is the categorical product of A and B in Par. -/
theorem pair_unique {X : Type} (f : X → Option α) (g : X → Option β) (h : X → Option (These α β))
    (h1 : kc h fst = f) (h2 : kc h snd = g) : h = pair f g := by
  funext x
  have e1 := congrFun h1 x
  have e2 := congrFun h2 x
  unfold kc at e1 e2; unfold pair
  rw [← e1, ← e2]
  rcases h x with _ | (a | b | ⟨a, b⟩) <;> rfl

/-- The empty type is terminal in Par: the only partial map into it is "undefined". -/
theorem empty_terminal {X : Type} (h : X → Option Empty) : h = fun _ => none := by
  funext x; exact option_empty (h x)

/-- And initial: there is exactly one map out of it. -/
theorem empty_initial {X : Type} (h k : Empty → Option X) : h = k := by
  funext e; exact nomatch e

/-- The terminal object is the unit of ✳ (checked above: `unitL`), as it must be for a product. -/
theorem empty_pair {X : Type} (g : X → Option β) :
    kc (pair (fun _ => (none : Option Empty)) g) snd = g := snd_pair _ g

/-! ## 3. Counting -/

def bothL : List α → List β → List (These α β)
  | [], _ => []
  | a :: as, bs => bs.map (both a) ++ bothL as bs

/-- Every value of A ✳ B, listed from lists of A and of B. -/
def allThese (la : List α) (lb : List β) : List (These α β) :=
  la.map this ++ lb.map that ++ bothL la lb

theorem bothL_length (la : List α) (lb : List β) : (bothL la lb).length = la.length * lb.length := by
  induction la with
  | nil => simp [bothL]
  | cons a as ih => simp [bothL, ih, Nat.succ_mul]; omega

/-- **The count**: |A ✳ B| = |A| ✳ |B|. -/
theorem allThese_length (la : List α) (lb : List β) :
    (allThese la lb).length = la.length ✳ lb.length := by
  simp [allThese, bothL_length, multAdd]; omega

theorem bothL_complete (la : List α) (lb : List β) (a : α) (b : β) (ha : a ∈ la) (hb : b ∈ lb) :
    both a b ∈ bothL la lb := by
  induction la with
  | nil => simp at ha
  | cons x xs ih =>
    simp only [bothL, List.mem_append, List.mem_map]
    rcases List.mem_cons.mp ha with rfl | h
    · exact Or.inl ⟨b, hb, rfl⟩
    · exact Or.inr (ih h)

/-- Nothing is missed: every value whose parts are listed is in the list. -/
theorem allThese_complete (la : List α) (lb : List β) (t : These α β)
    (h1 : ∀ a, a ∈ la) (h2 : ∀ b, b ∈ lb) : t ∈ allThese la lb := by
  unfold allThese
  rcases t with a | b | ⟨a, b⟩
  · simp [h1]
  · simp [h2]
  · simp only [List.mem_append]; exact Or.inr (bothL_complete la lb a b (h1 a) (h2 b))

/-- Pointed sets: the product of A₊ and B₊ has (a+1)(b+1) points; take away the base point and
    what is left is the wedge (a + b) plus the smash (ab). -/
theorem wedge_smash (a b : Nat) : (a + 1) * (b + 1) = 1 + (a + b) + a * b ∧ (a + 1) * (b + 1) = (a ✳ b) + 1 := by
  refine ⟨by simp [Nat.add_mul, Nat.mul_add]; omega, (multAdd_lift' a b).symm⟩

/-! ## 4. align and zipLongest: the software pattern -/

/-- Pair two lists position by position, keeping the overhang (Haskell's `align`). -/
def align : List α → List β → List (These α β)
  | [], ys => ys.map that
  | x :: xs, [] => (x :: xs).map this
  | x :: xs, y :: ys => both x y :: align xs ys

theorem align_length (xs : List α) (ys : List β) : (align xs ys).length = max xs.length ys.length := by
  induction xs generalizing ys with
  | nil => simp [align]
  | cons x xs ih =>
    cases ys with
    | nil => simp [align]
    | cons y ys => simp [align, ih]

theorem fm_fst_this (xs : List α) : xs.filterMap (fst ∘ (this : α → These α β)) = xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [List.filterMap_cons, Function.comp, fst]; exact congrArg (x :: ·) ih

theorem fm_snd_that (ys : List β) : ys.filterMap (snd ∘ (that : β → These α β)) = ys := by
  induction ys with
  | nil => rfl
  | cons y ys ih => simp only [List.filterMap_cons, Function.comp, snd]; exact congrArg (y :: ·) ih

/-- **Nothing lost**: the left list comes back out of the alignment. -/
theorem align_fst (xs : List α) (ys : List β) : (align xs ys).filterMap fst = xs := by
  induction xs generalizing ys with
  | nil => simp [align, List.filterMap_map, fst]
  | cons x xs ih =>
    cases ys with
    | nil => simp only [align, List.map, List.filterMap_cons, fst, List.filterMap_map, fm_fst_this]
    | cons y ys => simp [align, fst, ih]

/-- And so does the right list. -/
theorem align_snd (xs : List α) (ys : List β) : (align xs ys).filterMap snd = ys := by
  induction xs generalizing ys with
  | nil => simp only [align, List.filterMap_map, fm_snd_that]
  | cons x xs ih =>
    cases ys with
    | nil => simp [align, snd]
    | cons y ys => simp [align, snd, ih]

/-- The `both` corner. -/
def bothPart : These α β → Option (α × β)
  | both a b => some (a, b)
  | _ => none

/-- The ordinary zip is the corner cell of the alignment (the ab term). -/
theorem align_zip (xs : List α) (ys : List β) : (align xs ys).filterMap bothPart = xs.zip ys := by
  induction xs generalizing ys with
  | nil => simp [align, List.filterMap_map, bothPart]
  | cons x xs ih =>
    cases ys with
    | nil => simp [align, bothPart]
    | cons y ys => simp [align, bothPart, ih]

/-! ## 5. Decategorification: cell counts and the reduced Euler characteristic -/

/-- A cell count c₀ + c₁t + c₂t² + … evaluated at t. -/
def ev (t : Int) : List Int → Int
  | [] => 0
  | c :: cs => c + t * ev t cs

def padd : List Int → List Int → List Int
  | [], ys => ys
  | x :: xs, [] => x :: xs
  | x :: xs, y :: ys => (x + y) :: padd xs ys

/-- Cells of a product: an i-cell times a j-cell is an (i+j)-cell (convolution). -/
def pmul : List Int → List Int → List Int
  | [], _ => []
  | x :: xs, ys => padd (ys.map (x * ·)) (0 :: pmul xs ys)

theorem ev_padd (t : Int) (xs ys : List Int) : ev t (padd xs ys) = ev t xs + ev t ys := by
  induction xs generalizing ys with
  | nil => simp [padd, ev]
  | cons x xs ih =>
    cases ys with
    | nil => simp [padd, ev]
    | cons y ys => simp only [padd, ev, ih]; grind

theorem ev_smul (t c : Int) (ys : List Int) : ev t (ys.map (c * ·)) = c * ev t ys := by
  induction ys with
  | nil => simp [ev]
  | cons y ys ih => simp only [List.map, ev, ih]; grind

/-- Counting cells is multiplicative: P(X × Y) = P(X)·P(Y) at every t. -/
theorem ev_pmul (t : Int) (xs ys : List Int) : ev t (pmul xs ys) = ev t xs * ev t ys := by
  induction xs with
  | nil => simp [pmul, ev]
  | cons x xs ih => simp only [pmul, ev_padd, ev_smul, ev, ih]; grind

/-- **The reduced count mult-adds**: P(X×Y) − 1 = (P(X) − 1) ✳ (P(Y) − 1), at every t.
    At t = 1: non-base cells of a pointed product = wedge + smash. At t = −1: χ̃(X×Y) = χ̃X ✳ χ̃Y. -/
theorem reduced_product (t : Int) (xs ys : List Int) :
    ev t (pmul xs ys) - 1 = (ev t xs - 1) ✳ (ev t ys - 1) := by
  rw [ev_pmul]; unfold multAdd; grind

/-- The sphere Sⁿ: a base point and one n-cell (S⁰ is two points). -/
def sphere : Nat → List Int
  | 0 => [2]
  | n + 1 => 1 :: (List.replicate n 0 ++ [1])

theorem ev_top (t : Int) (n : Nat) : ev t (List.replicate n 0 ++ [1]) = t ^ n := by
  induction n with
  | zero => simp [ev]
  | succ n ih => simp only [List.replicate_succ, List.cons_append, ev, ih]; grind

/-- Spheres sit at ±1: χ̃(Sⁿ) = (−1)ⁿ. -/
theorem sphere_rchi (n : Nat) : ev (-1) (sphere n) - 1 = (-1) ^ n := by
  cases n with
  | zero => simp [sphere, ev]
  | succ n => simp only [sphere, ev, ev_top]; grind

/-- **The circle is the wall**: χ̃(S¹ × Y) = −1 for every Y. -/
theorem circle_wall (ys : List Int) : ev (-1) (pmul (sphere 1) ys) - 1 = -1 := by
  rw [reduced_product]
  have : ev (-1) (sphere 1) - 1 = -1 := by simp [sphere, ev]
  rw [this]; exact neg_one_absorbs _

/-- The torus S¹ × S¹: one 0-cell, two 1-cells, one 2-cell. -/
theorem torus_cells : pmul (sphere 1) (sphere 1) = [1, 2, 1] := by decide

/-- Even spheres sit at mult-add's 1: χ̃(S² × S²) = 1 ✳ 1 = 3. -/
theorem s2_s2 : ev (-1) (pmul (sphere 2) (sphere 2)) - 1 = 3 := by decide

/-- Suspension multiplies the reduced count by t, so the splitting
    Σ(X×Y) ≃ ΣX ∨ ΣY ∨ Σ(X∧Y) has the count t·(a ✳ b) = ta + tb + t·ab. -/
theorem suspension_split (t : Int) (xs ys : List Int) :
    t * (ev t (pmul xs ys) - 1) = t * (ev t xs - 1) + t * (ev t ys - 1) + t * ((ev t xs - 1) * (ev t ys - 1)) := by
  rw [reduced_product]; unfold multAdd; grind

end CategoryTheory
