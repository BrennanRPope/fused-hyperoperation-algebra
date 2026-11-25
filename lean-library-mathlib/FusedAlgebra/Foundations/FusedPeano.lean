module

public import FusedAlgebra.Rungs

/-!
# Fused Peano arithmetic: the natural numbers with mult-add as a primitive

Peano's axioms with the fused ladder in place of `+` and `*`. A model (`MAFound.FusedPeano`) is a
type `N` with a starting point `zero`, a next-step map `succ`, rung 1 (`sadd`, which will be
`a + b + 1`) and rung 2 (`madd`, which will be `a + b + a * b`), such that

* P1: `succ n ≠ zero`; P2: `succ` is injective;
* R1: `sadd a zero = succ a` and `sadd a (succ b) = succ (sadd a b)`;
* R2: `madd a zero = a` and `madd a (succ b) = sadd (madd a b) a`;
* I: induction.

Rung 1 is needed in the axioms for the same reason ordinary Peano arithmetic needs `+` before it
can define `*`: the recursion for rung 2 is written with rung 1.

## Main results

* From the axioms, by induction alone: `sadd` and `madd` are commutative and associative, `zero`
  is the identity of `madd`, and `madd` distributes over `sadd` exactly (`madd_sadd`).
* Categoricity: `equivNat : N ≃ ℕ`, the inverse of counting `toN : ℕ → N`, carries `sadd` to
  `MultAdd.sadd` and `madd` to `✳` (`toN_sadd`, `toN_madd`).
* `instCommSemiring`, `ringEquivNat : N ≃+* ℕ`: every model is a commutative semiring isomorphic
  to mathlib's `ℕ`; in it `succ x = x + 1`, `sadd x y = x + y + 1` and `madd x y = x ✳ y`
  (`succ_eq_add_one`, `sadd_eq`, `multAdd_eq`).
* `recN`: recursion on any model, and `ladN`, the whole fused ladder on any model from one
  schema; counting carries `MultAdd.rung` to it (`toN_rung`).
* `mul_definable`, `add_definable`: on `ℕ`, multiplication is definable from successor and `✳`,
  and addition from successor and multiplication (Julia Robinson's identity, 1949).

This is the mathlib form of `FusedAlgebra.Foundations.FusedPeano`. The mathematics is known:
categoricity of second-order Peano arithmetic (Dedekind), here for a different choice of
primitives, and Robinson's definability of `+` from successor and `*`. The fused ladder on `ℕ`
(`lad` in the original, with its lift `lad_lift` to the ordinary hyperoperations) is
`MultAdd.rung` and `MultAdd.rung_add_one` in `FusedAlgebra.Rungs`, built on mathlib's
`hyperoperation`.
-/

@[expose] public section

namespace MAFound

open scoped MultAdd

/-- A model of fused Peano arithmetic: zero, successor, rung 1 (`sadd`) and rung 2 (`madd`), with
the two recursions and induction. -/
class FusedPeano (N : Type*) where
  /-- The starting point. -/
  zero : N
  /-- The next-step map. -/
  succ : N → N
  /-- Rung 1, the fused sum: `a + b + 1`. -/
  sadd : N → N → N
  /-- Rung 2, mult-add: `a + b + a * b`. -/
  madd : N → N → N
  /-- P1: `zero` is where counting starts. -/
  succ_ne_zero : ∀ n, succ n ≠ zero
  /-- P2: different numbers have different next steps. -/
  succ_inj : ∀ m n, succ m = succ n → m = n
  /-- R1: joining a rod to the unit rod. -/
  sadd_zero : ∀ a, sadd a zero = succ a
  /-- R1: the recursion for rung 1. -/
  sadd_succ : ∀ a b, sadd a (succ b) = succ (sadd a b)
  /-- R2: a box one cell wide is just the row. -/
  madd_zero : ∀ a, madd a zero = a
  /-- R2: a box grows one row at a time. -/
  madd_succ : ∀ a b, madd a (succ b) = sadd (madd a b) a
  /-- Induction: what holds at `zero` and passes from `n` to `succ n` holds everywhere. -/
  ind : ∀ (P : N → Prop), P zero → (∀ n, P n → P (succ n)) → ∀ n, P n

namespace FusedPeano

variable {N : Type*} [FusedPeano N]

/-! ### Rung 1 from the axioms -/

theorem succ_sadd (a b : N) : sadd (succ a) b = succ (sadd a b) := by
  refine ind (fun b ↦ sadd (succ a) b = succ (sadd a b)) ?_ ?_ b
  · rw [sadd_zero, sadd_zero]
  · intro b ih; rw [sadd_succ, sadd_succ, ih]

theorem zero_sadd (b : N) : sadd zero b = succ b := by
  refine ind (fun b ↦ sadd zero b = succ b) ?_ ?_ b
  · rw [sadd_zero]
  · intro b ih; rw [sadd_succ, ih]

/-- Rung 1 is commutative: two rods joined either way round. -/
theorem sadd_comm (a b : N) : sadd a b = sadd b a := by
  refine ind (fun b ↦ sadd a b = sadd b a) ?_ ?_ b
  · rw [sadd_zero, zero_sadd]
  · intro b ih; rw [sadd_succ, succ_sadd, ih]

/-- Rung 1 is associative. -/
theorem sadd_assoc (a b c : N) : sadd (sadd a b) c = sadd a (sadd b c) := by
  refine ind (fun c ↦ sadd (sadd a b) c = sadd a (sadd b c)) ?_ ?_ c
  · rw [sadd_zero, sadd_zero, sadd_succ]
  · intro c ih; rw [sadd_succ, sadd_succ, sadd_succ, ih]

theorem sadd_right_cancel {a b c : N} (h : sadd a c = sadd b c) : a = b := by
  revert h
  refine ind (fun c ↦ sadd a c = sadd b c → a = b) ?_ ?_ c
  · intro h; rw [sadd_zero, sadd_zero] at h; exact succ_inj _ _ h
  · intro c ih h; rw [sadd_succ, sadd_succ] at h; exact ih (succ_inj _ _ h)

theorem sadd_ne_zero (a b : N) : sadd a b ≠ zero := by
  refine ind (fun b ↦ sadd a b ≠ zero) ?_ ?_ b
  · rw [sadd_zero]; exact succ_ne_zero a
  · intro b _; rw [sadd_succ]; exact succ_ne_zero _

/-! ### Rung 2 (mult-add) from the axioms -/

/-- `zero` is a left identity too. -/
theorem zero_madd (b : N) : madd zero b = b := by
  refine ind (fun b ↦ madd zero b = b) ?_ ?_ b
  · rw [madd_zero]
  · intro b ih; rw [madd_succ, ih, sadd_zero]

/-- Growing the box along its other side. -/
theorem succ_madd (a b : N) : madd (succ a) b = sadd (madd a b) b := by
  refine ind (fun b ↦ madd (succ a) b = sadd (madd a b) b) ?_ ?_ b
  · rw [madd_zero, madd_zero, sadd_zero]
  · intro b ih
    rw [madd_succ, ih, madd_succ, sadd_succ, sadd_succ, sadd_assoc, sadd_assoc, sadd_comm b a]

/-- Mult-add is commutative: the box turned on its side. -/
theorem madd_comm (a b : N) : madd a b = madd b a := by
  refine ind (fun b ↦ madd a b = madd b a) ?_ ?_ b
  · rw [madd_zero, zero_madd]
  · intro b ih; rw [madd_succ, succ_madd, ih]

/-- Rung 2 distributes over rung 1 exactly, with no correction term:
`madd a (sadd b c) = sadd (madd a b) (madd a c)`. Two boxes of the same width, stacked. -/
theorem madd_sadd (a b c : N) : madd a (sadd b c) = sadd (madd a b) (madd a c) := by
  refine ind (fun c ↦ madd a (sadd b c) = sadd (madd a b) (madd a c)) ?_ ?_ c
  · rw [sadd_zero, madd_succ, madd_zero]
  · intro c ih; rw [sadd_succ, madd_succ, ih, madd_succ, sadd_assoc]

/-- Mult-add is associative: a 3-D box built two ways. -/
theorem madd_assoc (a b c : N) : madd (madd a b) c = madd a (madd b c) := by
  refine ind (fun c ↦ madd (madd a b) c = madd a (madd b c)) ?_ ?_ c
  · rw [madd_zero, madd_zero]
  · intro c ih; rw [madd_succ, ih, madd_succ, madd_sadd]

/-! ### Categoricity: every model is a copy of `ℕ` -/

/-- Counting up in the model: `k ↦ succ (succ (… zero))`. -/
def toN : ℕ → N
  | 0 => zero
  | k + 1 => succ (toN k)

@[simp] theorem toN_zero : toN (N := N) 0 = zero := rfl

theorem toN_succ (k : ℕ) : toN (N := N) (k + 1) = succ (toN k) := rfl

theorem toN_injective : Function.Injective (toN (N := N))
  | 0, 0, _ => rfl
  | 0, _ + 1, h => absurd h.symm (succ_ne_zero _)
  | _ + 1, 0, h => absurd h (succ_ne_zero _)
  | _ + 1, _ + 1, h => congrArg Nat.succ (toN_injective (succ_inj _ _ h))

theorem toN_surjective : Function.Surjective (toN (N := N)) := fun n ↦ by
  refine ind (fun n ↦ ∃ k, toN k = n) ⟨0, rfl⟩ ?_ n
  rintro n ⟨k, rfl⟩
  exact ⟨k + 1, rfl⟩

theorem toN_bijective : Function.Bijective (toN (N := N)) := ⟨toN_injective, toN_surjective⟩

/-- Counting carries rung 1 to rung 1: `a + b + 1 ↦ sadd a b`. -/
theorem toN_sadd (a b : ℕ) : toN (N := N) (MultAdd.sadd a b) = sadd (toN a) (toN b) := by
  induction b with
  | zero => exact (sadd_zero _).symm
  | succ b ih => rw [toN_succ, sadd_succ, ← ih]; rfl

/-- Counting carries mult-add to mult-add: `a ✳ b ↦ madd a b`. -/
theorem toN_madd (a b : ℕ) : toN (N := N) (a ✳ b) = madd (toN a) (toN b) := by
  induction b with
  | zero => rw [MultAdd.multAdd_zero]; exact (madd_zero _).symm
  | succ b ih => rw [MultAdd.multAdd_succ, toN_sadd, ih, toN_succ, madd_succ, sadd_comm]

/-- Every model of the fused Peano axioms is in bijection with `ℕ`; the inverse map is counting,
`toN`. This is the way back from the model to `ℕ` (`fromN` in the original). -/
noncomputable def equivNat : N ≃ ℕ := (Equiv.ofBijective toN toN_bijective).symm

@[simp] theorem equivNat_symm_apply (k : ℕ) : (equivNat (N := N)).symm k = toN k := rfl

@[simp] theorem toN_equivNat (n : N) : toN (equivNat n) = n := equivNat.symm_apply_apply n

@[simp] theorem equivNat_toN (k : ℕ) : equivNat (toN (N := N) k) = k := equivNat.apply_symm_apply k

/-! ### The commutative semiring on any model -/

/-- Any model of the fused Peano axioms is a commutative semiring: the structure of `ℕ` carried
along `equivNat`. Its `+` and `*` are pinned down by the primitives: `succ (x + y) = sadd x y`
(`succ_add_eq_sadd`) and `madd x y = x + y + x * y` (`madd_eq`). -/
noncomputable instance instCommSemiring : CommSemiring N := equivNat.commSemiring

/-- Every model of the fused Peano axioms is isomorphic to `ℕ` as a semiring. -/
noncomputable def ringEquivNat : N ≃+* ℕ := equivNat.ringEquiv

@[simp] theorem ringEquivNat_apply (n : N) : ringEquivNat n = equivNat n := rfl

/-- Addition on any model cancels on the right (`add_right_cancel'` in the original). It is
proved here for the model itself, by carrying the equation to `ℕ` along `ringEquivNat`. -/
protected theorem add_right_cancel {a b c : N} (h : a + c = b + c) : a = b :=
  ringEquivNat.injective <| Nat.add_right_cancel (m := ringEquivNat c) <| by
    rw [← map_add, ← map_add, h]

/-- Addition on any model cancels on the left. -/
protected theorem add_left_cancel {a b c : N} (h : a + b = a + c) : b = c :=
  FusedPeano.add_right_cancel (c := a) <| by rw [add_comm b a, add_comm c a, h]

/-- Every model is a cancellative additive monoid, so mathlib's `add_right_cancel`,
`add_left_cancel` and their `iff` forms apply to it directly. -/
instance instIsCancelAdd : IsCancelAdd N where
  add_left_cancel _ _ _ := FusedPeano.add_left_cancel
  add_right_cancel _ _ _ := FusedPeano.add_right_cancel

@[simp] theorem ringEquivNat_symm_apply (k : ℕ) : (ringEquivNat (N := N)).symm k = toN k := rfl

theorem zero_eq : (0 : N) = zero := (map_zero (ringEquivNat (N := N)).symm).symm

theorem toN_one : toN (N := N) 1 = 1 := map_one (ringEquivNat (N := N)).symm

theorem toN_add (a b : ℕ) : toN (N := N) (a + b) = toN a + toN b :=
  map_add (ringEquivNat (N := N)).symm a b

theorem toN_mul (a b : ℕ) : toN (N := N) (a * b) = toN a * toN b :=
  map_mul (ringEquivNat (N := N)).symm a b

/-- Counting is the canonical map `ℕ → N` of the semiring. -/
theorem natCast_eq_toN (k : ℕ) : (k : N) = toN k :=
  (map_natCast (ringEquivNat (N := N)).symm k).symm

theorem succ_eq_add_one (x : N) : succ x = x + 1 := by
  obtain ⟨a, rfl⟩ := toN_surjective x
  rw [← toN_succ, toN_add, toN_one]

/-- Rung 1 is `x + y + 1`. -/
theorem sadd_eq (x y : N) : sadd x y = x + y + 1 := by
  obtain ⟨a, rfl⟩ := toN_surjective x
  obtain ⟨b, rfl⟩ := toN_surjective y
  rw [← toN_sadd, MultAdd.sadd, toN_add, toN_add, toN_one]

/-- Addition from the primitives: `x + y` is the number one step below `sadd x y`. -/
theorem succ_add_eq_sadd (x y : N) : succ (x + y) = sadd x y := by
  rw [succ_eq_add_one, sadd_eq]

/-- In the derived semiring, the library's mult-add `a ✳ b = a + b + a * b` is the primitive
`madd`. -/
theorem multAdd_eq (x y : N) : x ✳ y = madd x y := by
  obtain ⟨a, rfl⟩ := toN_surjective x
  obtain ⟨b, rfl⟩ := toN_surjective y
  rw [← toN_madd, MultAdd.multAdd_def, MultAdd.multAdd_def, toN_add, toN_add, toN_mul]

/-- Multiplication from the primitives: `x * y` is what is left of `madd x y` after taking off
`x` and `y`. -/
theorem madd_eq (x y : N) : madd x y = x + y + x * y := (multAdd_eq x y).symm

/-! ### Recursion on any model, from categoricity -/

/-- Define a function on the model by "value at `zero`" and "how to take a step". -/
noncomputable def recN {α : Type*} (z : α) (s : N → α → α) (n : N) : α :=
  Nat.rec (motive := fun _ ↦ α) z (fun k acc ↦ s (toN k) acc) (equivNat n)

theorem recN_zero {α : Type*} (z : α) (s : N → α → α) : recN z s (zero : N) = z := by
  rw [recN, ← toN_zero, equivNat_toN]; rfl

theorem recN_succ {α : Type*} (z : α) (s : N → α → α) (n : N) :
    recN z s (succ n) = s n (recN z s n) := by
  obtain ⟨k, rfl⟩ := toN_surjective n
  rw [recN, recN, ← toN_succ, equivNat_toN, equivNat_toN]

/-! ### The fused ladder, from one schema

Rung 0 is the step, and every rung above it is the rung below, iterated:
`F 0 a b = succ b`, `F 1 a zero = succ a`, `F (n + 2) a zero = a` (from rung 2 up, `zero` is the
right identity), and `F (n + 1) a (succ b) = F n a (F (n + 1) a b)`.

On `ℕ` this schema is `MultAdd.rung` (`rung_zero`, `rung_one_zero_right`,
`rung_add_two_zero_right` below and `MultAdd.rung_succ_succ`). -/

/-- Rung 0 of the fused ladder is the step. -/
theorem rung_zero (a b : ℕ) : MultAdd.rung 0 a b = b + 1 := by
  simp [MultAdd.rung, hyperoperation_zero]

/-- Rung 1 starts at `a + 1`. -/
theorem rung_one_zero_right (a : ℕ) : MultAdd.rung 1 a 0 = a + 1 := by
  simp [MultAdd.rung_one, MultAdd.sadd]

/-- From rung 2 up, `0` is the right identity. -/
theorem rung_add_two_zero_right (n a : ℕ) : MultAdd.rung (n + 2) a 0 = a := by
  simp [MultAdd.rung, hyperoperation_ge_two_eq_self]

/-- The fused ladder on any model, defined by the schema with `recN`. -/
noncomputable def ladN : ℕ → N → N → N
  | 0, _, b => succ b
  | 1, a, b => recN (succ a) (fun _ acc ↦ succ acc) b
  | n + 2, a, b => recN a (fun _ acc ↦ ladN (n + 1) a acc) b

theorem ladN_zero (a b : N) : ladN 0 a b = succ b := rfl

theorem ladN_zero_right (n : ℕ) (a : N) :
    ladN (n + 1) a zero = if n = 0 then succ a else a := by
  cases n <;> simp [ladN, recN_zero]

theorem ladN_succ (n : ℕ) (a b : N) : ladN (n + 1) a (succ b) = ladN n a (ladN (n + 1) a b) := by
  cases n <;> simp only [ladN, recN_succ]

/-- The schema's rung 1 is the primitive `sadd`. -/
theorem ladN_one (a b : N) : ladN 1 a b = sadd a b := by
  refine ind (fun b ↦ ladN 1 a b = sadd a b) ?_ ?_ b
  · rw [ladN_zero_right, sadd_zero]; rfl
  · intro b ih; rw [ladN_succ, ih, sadd_succ]; rfl

/-- The schema's rung 2 is the primitive `madd`. -/
theorem ladN_two (a b : N) : ladN 2 a b = madd a b := by
  refine ind (fun b ↦ ladN 2 a b = madd a b) ?_ ?_ b
  · rw [ladN_zero_right, madd_zero]; rfl
  · intro b ih; rw [ladN_succ, ih, ladN_one, madd_succ, sadd_comm]

/-- Counting carries the whole ladder over: every model has every rung, and rung `n` of the
model is `MultAdd.rung n`, the hyperoperation `n` seen through `x ↦ x + 1`. -/
theorem toN_rung (n a b : ℕ) :
    toN (N := N) (MultAdd.rung n a b) = ladN n (toN a) (toN b) := by
  induction n generalizing b with
  | zero => rw [rung_zero, ladN_zero, toN_succ]
  | succ n ihn =>
    induction b with
    | zero =>
      rw [toN_zero, ladN_zero_right]
      cases n with
      | zero => rw [rung_one_zero_right, toN_succ]; rfl
      | succ m => rw [rung_add_two_zero_right]; rfl
    | succ b ihb => rw [MultAdd.rung_succ_succ, ihn, ihb, toN_succ, ladN_succ]

/-! ### Strength: `+` and `*` can be defined from successor and mult-add

Mult-add alone on `ℕ` is multiplication in disguise (`a ✳ b + 1 = (a + 1) * (b + 1)`), and the
first-order theory of multiplication alone is decidable (Skolem 1930, Mostowski 1952): too weak
to carry arithmetic. With the successor added, `+` and `*` become definable, so the structure
`(ℕ, 0, succ, ✳)` carries all of first-order arithmetic. -/

/-- `*` from successor and `✳`: `m * n = k` exactly when one side is `0` and `k = 0`, or
`m = m' + 1`, `n = n' + 1`, `k = k' + 1` and `m' ✳ n' = k'`. -/
theorem mul_definable (m n k : ℕ) :
    m * n = k ↔ (m = 0 ∧ k = 0) ∨ (n = 0 ∧ k = 0) ∨
      ∃ m' n' k', m = m' + 1 ∧ n = n' + 1 ∧ k = k' + 1 ∧ m' ✳ n' = k' := by
  constructor
  · rintro rfl
    rcases m with _ | m'
    · simp
    rcases n with _ | n'
    · simp
    · exact .inr <| .inr ⟨m', n', m' ✳ n', rfl, rfl, (MultAdd.multAdd_add_one m' n').symm, rfl⟩
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨m', n', k', rfl, rfl, rfl, rfl⟩)
    · simp
    · simp
    · exact (MultAdd.multAdd_add_one m' n').symm

/-- `+` from successor and `*`: Julia Robinson's identity (1949). For `z ≠ 0`, `x + y = z` exactly
when `(x * z + 1) * (y * z + 1) = z * z * (x * y + 1) + 1`. -/
theorem add_definable (x y z : ℕ) :
    x + y = z ↔ (z = 0 ∧ x = 0 ∧ y = 0) ∨
      (z ≠ 0 ∧ (x * z + 1) * (y * z + 1) = z * z * (x * y + 1) + 1) := by
  constructor
  · rintro rfl
    rcases Nat.eq_zero_or_pos (x + y) with h | h
    · exact .inl ⟨h, by lia, by lia⟩
    · exact .inr ⟨h.ne', by ring⟩
  · rintro (⟨rfl, rfl, rfl⟩ | ⟨hz, h⟩)
    · rfl
    · have e : z * (x + y) = z * z := by linarith
      exact Nat.eq_of_mul_eq_mul_left (Nat.pos_of_ne_zero hz) e

end FusedPeano

end MAFound
