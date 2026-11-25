module

public import FusedAlgebra.Basic

/-!
# Infinities between countable and the continuum, through mult-add

## Main results

* `sumOnes_eq`: `n` ones added make `n` (the countable side).
* `maddOnes_succ`: `n` ones mult-added make `2 ^ n - 1` (the corner side).
* `layers_eval`: the mult-add cube's layers `n.choose k`, weighted by `t ^ k`, fill
  `(1 + t) ^ n`.
* `wOnes_spec`: the zoomed mult-add `a + b + t * a * b`, folded `n` times on `1`, gives
  `t * W + 1 = (1 + t) ^ n`; `t = 0` is plain addition (`wOnes_zero`).
* `cantor`: no `α`-indexed list of 0/1 sequences contains them all: the infinite cube has more
  corners than axes. `single_inj`: but every axis gives its own corner.
* `half_iter_upper`, `half_iter_lower`, `half_not_poly`: a "half rung" `f` with `f ∘ f = g` is
  squeezed: if `f` were polynomially bounded, `g` would be.

All of this is known: `layers_eval` is the binomial theorem (`add_pow`), `cantor` is Cantor's
theorem (`Function.cantor_surjective`), the rest are short inductions. The core-Lean original
defined its own folds `sumOnes`, `maddOnes`, `wOnes` and Pascal rows `layers`; here they are
`Finset.sum`, `MultAdd.mpow`, `MultAdd.wpow` and `Nat.choose`.
-/

@[expose] public section

namespace Infinities

open MultAdd Finset

/-! ### Sum of ones versus mult-add of ones -/

/-- `n` ones added make `n`. -/
theorem sumOnes_eq (n : ℕ) : ∑ _i ∈ range n, 1 = n := by simp

/-- `n` ones mult-added make `2 ^ n - 1`. -/
theorem maddOnes_succ (n : ℕ) : mpow 1 n + 1 = 2 ^ n := by
  rw [add_comm, one_add_mpow]

/-! ### The mult-add cube's layers (Pascal), weighted by `t` -/

/-- Layer `k` of the `n`-cube is the `n.choose k` corners with `k` switches on. Weighted by
`t ^ k`, the layers fill `(1 + t) ^ n`. This is the binomial theorem, `add_pow`. -/
theorem layers_eval {R : Type*} [CommSemiring R] (t : R) (n : ℕ) :
    ∑ k ∈ range (n + 1), (n.choose k : R) * t ^ k = (1 + t) ^ n := by
  rw [add_comm 1 t, add_pow]
  simp [mul_comm]

example : (List.range 4).map (Nat.choose 3) = [1, 3, 3, 1] := by decide

/-! ### Zoomed mult-add: `a + b + t * a * b` (`t = 0` addition, `t = 1` mult-add) -/

/-- The zoomed mult-add folded `n` times on `1` lifts to `(1 + t) ^ n`. -/
theorem wOnes_spec {R : Type*} [CommSemiring R] (t : R) (n : ℕ) :
    t * wpow t 1 n + 1 = (1 + t) ^ n := by
  simpa [add_comm] using one_add_mul_wpow t 1 n

/-- At `t = 0` the fold is plain counting. -/
theorem wOnes_zero {R : Type*} [CommSemiring R] (n : ℕ) : wpow (0 : R) 1 n = n := by
  induction n with
  | zero => simp
  | succ n ih => simp [wpow_succ, ih]

/-! ### The infinite cube: more corners than axes (Cantor) -/

/-- Cantor's diagonal argument: no `α`-indexed family of 0/1 sequences contains them all. In
mathlib this is `Function.cantor_surjective` (stated with `Set α`). -/
theorem cantor {α : Type*} (f : α → (α → Bool)) : ∃ g : α → Bool, ∀ a, f a ≠ g := by
  refine ⟨fun a ↦ !(f a a), fun a h ↦ ?_⟩
  have h2 := congrFun h a
  simp at h2

/-- The corner of the infinite cube with only switch `n` on. -/
def single (n : ℕ) : ℕ → Bool := fun k ↦ k == n

/-- Every axis gives its own corner, so there are at least as many corners as axes. -/
theorem single_inj : Function.Injective single := by
  intro m n h
  have h2 := congrFun h m
  simpa [single] using h2.symm

/-! ### Half rungs are squeezed between their square and the whole rung -/

section Preorder

variable {α : Type*} [Preorder α]

theorem half_iter_upper (f g h : α → α) (N : α)
    (hfg : ∀ n, f (f n) = g n)
    (hgrow : ∀ n, n ≤ f n)
    (hmono : Monotone h)
    (hbound : ∀ n, N ≤ n → f n ≤ h n) :
    ∀ n, N ≤ n → g n ≤ h (h n) := by
  intro n hn
  rw [← hfg n]
  exact (hbound (f n) (hn.trans (hgrow n))).trans (hmono (hbound n hn))

theorem half_iter_lower (f g h : α → α) (N : α)
    (hfg : ∀ n, f (f n) = g n)
    (hfmono : Monotone f)
    (hh : ∀ n, N ≤ n → N ≤ h n)
    (hbound : ∀ n, N ≤ n → h n ≤ f n) :
    ∀ n, N ≤ n → h (h n) ≤ g n := by
  intro n hn
  rw [← hfg n]
  exact (hbound (h n) (hh n hn)).trans (hfmono (hbound n hn))

end Preorder

/-- If a half rung stayed under `n ^ k`, the whole rung would stay under `n ^ (k * k)`. Since the
whole rung (exponential) outgrows every polynomial, no half rung is polynomial. -/
theorem half_not_poly (f g : ℕ → ℕ) (N k : ℕ)
    (hfg : ∀ n, f (f n) = g n)
    (hgrow : ∀ n, n ≤ f n)
    (hbound : ∀ n, N ≤ n → f n ≤ n ^ k) :
    ∀ n, N ≤ n → g n ≤ n ^ (k * k) := by
  intro n hn
  rw [pow_mul]
  exact half_iter_upper f g (fun m ↦ m ^ k) N hfg hgrow
    (fun _ _ hab ↦ Nat.pow_le_pow_left hab k) hbound n hn

end Infinities
