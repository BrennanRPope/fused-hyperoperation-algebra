/-
  Fractional power sets: can there be a dial between X and P(X)?
  ===============================================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  ../index.html, ../checks.py (18 of 18 pass).

  STATUS: core Lean 4 only (no mathlib, no imports), no `sorry`.
  UNCOMPILED here: this workspace has no Lean toolchain. Written for the
  shared library (core Lean 4.37.0-rc1); please compile and merge as
  Topics/FractionalPowerset.

  What is proved:
  * no_monotone_half_pow2   the finite wall: no non-decreasing f : ℕ → ℕ has
                            f (f n) = 2ⁿ. Only n = 0 and n = 1 are used.
  * strictMono_add          a strictly increasing f gains at least k over k steps.
  * half_shift_even         if f is strictly increasing and f (f n) = n + d,
                            then d = 2 · f 0. So "jump d successors" has a
                            monotone square root only when d is even.
  * no_half_shift_odd       the odd case, as a contradiction.
  * no_descending           no strictly decreasing sequence of naturals.
  * no_descending_wf        the same for any well-founded relation.
  * dial_inflationary       in a dyadic dial (each level the square root of the
                            one above, all monotone, the top strictly growing),
                            every level is strictly growing.
  * no_dyadic_dial          so on a well-ordered scale of sizes the dial is
                            impossible: the values P^(1/2^k)(0) would fall forever.
                            This is the ZFC obstruction with ℕ standing in for
                            the well-ordered alephs.
  * flat_mul_flat_add, flat_add_flat_succ, flat_circ_iff_flat_mul
                            the rung chain for "which rung is flat on infinite
                            sizes", from rules that hold for cardinals in ZF
                            (Cantor–Bernstein, monotone + and ·, m·2 = m + m).
-/

namespace FractionalPowerset

/-! ## 1. The finite wall -/

theorem no_monotone_half_pow2 (f : Nat → Nat)
    (mono : ∀ m n, m ≤ n → f m ≤ f n) (h : ∀ n, f (f n) = 2 ^ n) : False := by
  have h0 : f (f 0) = 1 := (h 0).trans (by decide)
  have h1 : f (f 1) = 2 := (h 1).trans (by decide)
  have hle : f 0 ≤ f (f 0) := mono 0 (f 0) (Nat.zero_le _)
  generalize ha : f 0 = a at h0 hle
  have ha1 : a = 1 := by
    cases a with
    | zero => omega
    | succ k => omega
  subst ha1
  rw [h0] at h1
  omega

/-! ## 2. Square roots of a shift: only even jumps have one -/

theorem strictMono_add (f : Nat → Nat) (hs : ∀ n, f n < f (n + 1)) :
    ∀ m k, f m + k ≤ f (m + k) := by
  intro m k
  induction k with
  | zero => simp
  | succ k ih =>
    have := hs (m + k)
    show f m + (k + 1) ≤ f (m + k + 1)
    omega

theorem half_shift_even (f : Nat → Nat) (d : Nat)
    (hs : ∀ n, f n < f (n + 1)) (h : ∀ n, f (f n) = n + d) : d = 2 * f 0 := by
  have ha : f (f 0) = d := by have := h 0; omega
  have hd : f d = f 0 + d := by
    have := h (f 0)
    rw [ha] at this
    exact this
  have l1 := strictMono_add f hs 0 (f 0)
  simp only [Nat.zero_add] at l1
  have hle : f 0 ≤ d := by omega
  have l2 := strictMono_add f hs (f 0) (d - f 0)
  rw [Nat.add_sub_of_le hle] at l2
  omega

theorem no_half_shift_odd (f : Nat → Nat) (d : Nat) (hodd : d % 2 = 1)
    (hs : ∀ n, f n < f (n + 1)) (h : ∀ n, f (f n) = n + d) : False := by
  have := half_shift_even f d hs h
  omega

/-! ## 3. Well-ordered sizes cannot carry a dyadic dial -/

theorem no_descending (s : Nat → Nat) (h : ∀ k, s (k + 1) < s k) : False := by
  have key : ∀ k, s k + k ≤ s 0 := by
    intro k
    induction k with
    | zero => simp
    | succ k ih => have := h k; omega
  have := key (s 0 + 1)
  omega

theorem no_descending_wf {α : Type} (lt : α → α → Prop) (wf : WellFounded lt)
    (s : Nat → α) (h : ∀ k, lt (s (k + 1)) (s k)) : False := by
  have key : ∀ a, ∀ k, s k ≠ a := fun a =>
    wf.induction (C := fun x => ∀ k, s k ≠ x) a
      (fun x ih k hk => ih (s (k + 1)) (hk ▸ h k) (k + 1) rfl)
  exact key (s 0) 0 rfl

/-- `F k` plays the role of `P^(1/2^k)`: each level composed with itself
    gives the level above, all levels are monotone, and the top level
    strictly grows (Cantor: `X < P(X)`). -/
theorem dial_inflationary (F : Nat → Nat → Nat)
    (mono : ∀ k x y, x ≤ y → F k x ≤ F k y)
    (half : ∀ k x, F (k + 1) (F (k + 1) x) = F k x)
    (grow : ∀ x, x < F 0 x) : ∀ k x, x < F k x := by
  intro k
  induction k with
  | zero => exact grow
  | succ k ih =>
    intro x
    rcases Nat.lt_or_ge x (F (k + 1) x) with hlt | hge
    · exact hlt
    · exfalso
      have h1 : F (k + 1) (F (k + 1) x) ≤ F (k + 1) x := mono (k + 1) _ _ hge
      have h2 := ih x
      rw [← half k x] at h2
      omega

theorem no_dyadic_dial (F : Nat → Nat → Nat)
    (mono : ∀ k x y, x ≤ y → F k x ≤ F k y)
    (half : ∀ k x, F (k + 1) (F (k + 1) x) = F k x)
    (grow : ∀ x, x < F 0 x) : False := by
  have infl := dial_inflationary F mono half grow
  apply no_descending (fun k => F k 0)
  intro k
  have := infl (k + 1) (F (k + 1) 0)
  rw [half k 0] at this
  exact this

/-! ## 4. Which rung is flat on infinite sizes

  Abstract cardinal rules that hold in ZF: `le` is antisymmetric
  (Cantor–Bernstein), `+` and `·` are monotone, `m·2 = m+m`,
  `a ≤ a + b`. Flat multiplication forces flat addition, which forces a flat
  successor. (Tarski: flat multiplication on every infinite size is the axiom
  of choice; flat successor says every infinite set is Dedekind-infinite.) -/

theorem flat_mul_flat_add {C : Type} (add mul : C → C → C) (le : C → C → Prop)
    (antisymm : ∀ a b, le a b → le b a → a = b)
    (two : C) (two_def : ∀ m, mul m two = add m m)
    (mul_mono : ∀ m a b, le a b → le (mul m a) (mul m b))
    (self_le_add : ∀ a b, le a (add a b))
    (m : C) (h2 : le two m) (hm : mul m m = m) : add m m = m := by
  apply antisymm
  · have := mul_mono m two m h2
    rw [two_def, hm] at this
    exact this
  · exact self_le_add m m

theorem flat_add_flat_succ {C : Type} (add : C → C → C) (le : C → C → Prop)
    (antisymm : ∀ a b, le a b → le b a → a = b)
    (one : C) (add_mono : ∀ m a b, le a b → le (add m a) (add m b))
    (self_le_add : ∀ a b, le a (add a b))
    (m : C) (h1 : le one m) (hm : add m m = m) : add m one = m := by
  apply antisymm
  · have := add_mono m one m h1
    rw [hm] at this
    exact this
  · exact self_le_add m one

/-- Mult-add `m ◯ m = (m + m) + m·m` is flat exactly when multiplication is. -/
theorem flat_circ_iff_flat_mul {C : Type} (add mul : C → C → C) (le : C → C → Prop)
    (antisymm : ∀ a b, le a b → le b a → a = b)
    (two : C) (two_def : ∀ m, mul m two = add m m)
    (mul_mono : ∀ m a b, le a b → le (mul m a) (mul m b))
    (self_le_add : ∀ a b, le a (add a b))
    (add_comm : ∀ a b, add a b = add b a)
    (m : C) (h2 : le two m) (self_le_sq : le m (mul m m)) :
    add (add m m) (mul m m) = m ↔ mul m m = m := by
  constructor
  · intro hm
    apply antisymm
    · have := self_le_add (mul m m) (add m m)
      rw [add_comm (mul m m) (add m m), hm] at this
      exact this
    · exact self_le_sq
  · intro hm
    have hadd := flat_mul_flat_add add mul le antisymm two two_def mul_mono self_le_add m h2 hm
    rw [hm, hadd, hadd]

end FractionalPowerset
