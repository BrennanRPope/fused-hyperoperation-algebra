/-
  FusedAlgebra.Topics.HyperLadder: the fused hyperoperation ladder, made complete  (DRAFT, Nov 26 2026)
  ======================================================================================================

  Companion: topics/hyperoperation-ladder/README.md and checks.py. Core Lean 4.37.0-rc1 only
  (no mathlib), no `sorry`, standard axioms only.

  The repo uses the word "rung" for two different ladders. This file proves what each one is.

  Part A. The whole-rung ladder (`lad`, from Foundations.FusedPeano). Rung n is Goodstein's
  hyperoperation H n read through x ↦ 1 + x (`lad_lift`, already proved there for every n).
  Proved here, for every n at once:
    * one definition behind the named rungs: lad 1 = sadd, 2 = circ, 3 = fexp, 4 = ftet, 5 = fpent;
    * uniqueness: the fused recursion and its three starting values pin the ladder down;
    * Knuth's arrows: lad (m+3) a b + 1 = (a+1) ↑^(m+1) (b+1);
    * 0 is a right identity from rung 2 up and a left absorber from rung 3 up;
    * the fixed corner lad n 1 1 = 3 for every n ≥ 1 (Goodstein's 2 # 2 = 4, fused);
    * rung n+1 at 1 is rung n on the diagonal: lad (n+2) a 1 = lad (n+1) a a;
    * exact characterizations: commutative ⟺ n ∈ {1, 2}; associative ⟺ n ∈ {1, 2};
      a two-sided identity ⟺ n = 2; right identity 0 ⟺ n ≥ 2; left absorber 0 ⟺ n ≥ 3;
    * distributivity of rung n + 1 over rung n, exactly: on the left ⟺ n ≤ 1, on the right
      ⟺ n ≤ 2. The giant rungs (n ≥ 4) are settled by remainders mod 7, since their values
      cannot be computed (`lad_right_dist_iff`, `lad_left_dist_iff`);
    * no rung other than 1 and 2 is the transport of any commutative operation, by any maps.

  Part B. The height ladder (B_s, the operation a +⁽ˢ⁾ b of fractional-rungs/), stated for any
  flow E_s with E_s ∘ E_t = E_(s+t) and E_s 0 = 0. Its existence for E = eˣ − 1 on [0, ∞) is
  Szekeres's regular iteration (1958), cited, not proved here. Given a flow:
    * every height is commutative, associative, with identity 0, and B_0 = +;
    * heights climb by transport: B_(s+t) = (E_s)₊(B_t);
    * if E_1 carries + to the circle operation (eᵘ⁺ᵛ − 1 = (eᵘ − 1) ◯ (eᵛ − 1)), then B_1 = ◯,
      the hypernatural log turns ◯ into the height-(1 − s) operation, heights compose, and
      B_(s+1) distributes over its partner (E_s)₊(rung 1);
    * B_1 does not distribute over B_0 (◯ over +);
    * so the only rung the height ladder shares with the whole-rung ladder is rung 2 (◯).

  Part C. The dial k ≠ 0: rung 1 at dial k is a + b + 1/k, rung 2 is ◯ₖ, rung 3 is
  ((1 + ka)^m − 1)/k at 1 + kb = m. Each is the k = 1 rung rescaled, and rung 3 is rung 2
  iterated at every k ≠ 0. At k = 0 the lift 1 + kx is not one-to-one, so there is no ladder.
-/
import FusedAlgebra.Rungs
import FusedAlgebra.Foundations.FusedPeano

namespace HyperLadder

open CircleOp
open CircFound.FusedPeano (lad hyp lad_lift hyp_one hyp_add hyp_mul hyp_pow lad_one lad_two)

set_option linter.unusedSectionVars false

/-! ## Part A. The whole-rung ladder -/

/-! ### A1. The recursion, as equations -/

theorem lad_zero (a b : Nat) : lad 0 a b = b + 1 := rfl
theorem lad_one_zero (a : Nat) : lad 1 a 0 = a + 1 := rfl
/-- From rung 2 up, 0 is a right identity. -/
theorem lad_zero_right (n a : Nat) : lad (n + 2) a 0 = a := rfl
/-- Each rung is the rung below, iterated. -/
theorem lad_succ (n a b : Nat) : lad (n + 1) a (b + 1) = lad n a (lad (n + 1) a b) := rfl

theorem hyp_succ (n a b : Nat) : hyp (n + 1) a (b + 1) = hyp n a (hyp (n + 1) a b) := rfl

/-- Uniqueness: any ladder with the fused recursion and the same three starting values
    (step, a + 1 at rung 1, a from rung 2 up) is `lad`. -/
theorem lad_unique (G : Nat → Nat → Nat → Nat)
    (h0 : ∀ a b, G 0 a b = b + 1) (h1 : ∀ a, G 1 a 0 = a + 1)
    (h2 : ∀ n a, G (n + 2) a 0 = a)
    (hs : ∀ n a b, G (n + 1) a (b + 1) = G n a (G (n + 1) a b)) :
    ∀ n a b, G n a b = lad n a b := by
  intro n
  induction n with
  | zero => intro a b; rw [h0]; rfl
  | succ n ih =>
    intro a b
    induction b with
    | zero =>
      cases n with
      | zero => rw [h1]; rfl
      | succ m => rw [h2]; rfl
    | succ b ihb => rw [hs, ihb, ih]; rfl

/-! ### A2. One definition behind the named rungs -/

theorem hyp_four (x n : Nat) : hyp 4 x n = tet x n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [hyp_succ, ih, hyp_pow]; rfl

theorem hyp_five (x n : Nat) : hyp 5 x n = pent x n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [hyp_succ, ih, hyp_four]; rfl

theorem lad_eq_sadd (a b : Nat) : lad 1 a b = sadd a b := by rw [lad_one]; rfl
theorem lad_eq_circ (a b : Nat) : lad 2 a b = circ a b := by rw [lad_two]; rfl
theorem lad_eq_fexp (a b : Nat) : lad 3 a b = fexp a b := by
  apply lift_inj; rw [lad_lift, hyp_pow, fexp_lift]
theorem lad_eq_ftet (a b : Nat) : lad 4 a b = ftet a b := by
  apply lift_inj; rw [lad_lift, hyp_four, ftet_lift]
theorem lad_eq_fpent (a b : Nat) : lad 5 a b = fpent a b := by
  apply lift_inj; rw [lad_lift, hyp_five, fpent_lift]

/-! ### A3. Knuth's up-arrows -/

/-- `arrow m a b` is Knuth's a ↑^(m+1) b: one arrow is a power, and each further arrow
    iterates the one before, starting from 1. -/
def arrow : Nat → Nat → Nat → Nat
  | 0 => fun a b => a ^ b
  | m + 1 => fun a b => Nat.rec (motive := fun _ => Nat) 1 (fun _ acc => arrow m a acc) b

theorem hyp_arrow (m a b : Nat) : hyp (m + 3) a b = arrow m a b := by
  induction m generalizing b with
  | zero => rw [hyp_pow]; rfl
  | succ m ih =>
    induction b with
    | zero => rfl
    | succ b ihb =>
      show hyp (m + 3) a (hyp (m + 4) a b) = arrow m a (arrow (m + 1) a b)
      rw [ih, ihb]

/-- From rung 3 up the fused ladder is Knuth's arrows seen through x ↦ 1 + x. -/
theorem lad_arrow (m a b : Nat) : lad (m + 3) a b + 1 = arrow m (a + 1) (b + 1) := by
  rw [lad_lift, hyp_arrow]

theorem arrow_values : arrow 0 2 3 = 8 ∧ arrow 1 2 3 = 16 ∧ arrow 1 3 2 = 27 ∧ arrow 2 2 3 = 65536 := by
  decide

/-! ### A4. Laws at every rung -/

/-- From rung 3 up, 1 is a left absorber of the ordinary ladder: 1 ↑^m b = 1. -/
theorem hyp_one_left : ∀ n b, hyp (n + 3) 1 b = 1
  | 0, b => by rw [hyp_pow, Nat.one_pow]
  | n + 1, b => by
    induction b with
    | zero => rfl
    | succ b ihb => rw [hyp_succ, ihb, hyp_one_left n 1]

/-- From rung 3 up, 0 is a left absorber of the fused ladder. -/
theorem lad_zero_left (n b : Nat) : lad (n + 3) 0 b = 0 := by
  have h := lad_lift (n + 3) 0 b; rw [hyp_one_left] at h; omega

/-- Goodstein's 2 # 2 = 4 on every rung from addition up. -/
theorem hyp_two_two : ∀ n, hyp (n + 1) 2 2 = 4
  | 0 => by rw [hyp_add]
  | n + 1 => by
    show hyp (n + 1) 2 (hyp (n + 2) 2 1) = 4
    rw [hyp_one n 2]; exact hyp_two_two n

/-- The fixed corner: 1 ⋆ 1 = 3 on every fused rung from 1 up. -/
theorem lad_one_one (n : Nat) : lad (n + 1) 1 1 = 3 := by
  have h := lad_lift (n + 1) 1 1; rw [hyp_two_two] at h; omega

/-- One step up at 1 is the rung below on the diagonal: a ⋆ₙ₊₁ 1 = a ⋆ₙ a. -/
theorem lad_diag (n a : Nat) : lad (n + 2) a 1 = lad (n + 1) a a := rfl

theorem lad_zero_one_one : lad 0 1 1 = 2 := rfl

/-! ### A5. Exact characterizations of the laws -/

/-- Commutative. -/
def IsComm (op : Nat → Nat → Nat) : Prop := ∀ a b, op a b = op b a
/-- Associative. -/
def IsAssoc (op : Nat → Nat → Nat) : Prop := ∀ a b c, op (op a b) c = op a (op b c)
/-- Has a two-sided identity. -/
def HasIdentity (op : Nat → Nat → Nat) : Prop := ∃ e, ∀ a, op e a = a ∧ op a e = a

theorem sadd_comm' (a b : Nat) : sadd a b = sadd b a := by unfold sadd; omega
theorem sadd_assoc' (a b c : Nat) : sadd (sadd a b) c = sadd a (sadd b c) := by unfold sadd; omega

/-- Only rungs 1 and 2 are commutative. -/
theorem lad_comm_iff (n : Nat) : IsComm (lad n) ↔ n = 1 ∨ n = 2 := by
  constructor
  · intro h
    match n with
    | 0 => have := h 0 1; simp [lad_zero] at this
    | 1 => exact Or.inl rfl
    | 2 => exact Or.inr rfl
    | n + 3 =>
      have := h 0 1
      rw [lad_zero_left, lad_zero_right (n + 1) 1] at this
      exact absurd this (by decide)
  · rintro (rfl | rfl)
    · intro a b; rw [lad_eq_sadd, lad_eq_sadd]; exact sadd_comm' a b
    · intro a b; rw [lad_eq_circ, lad_eq_circ]; exact circ_comm a b

/-- Only rungs 1 and 2 are associative. -/
theorem lad_assoc_iff (n : Nat) : IsAssoc (lad n) ↔ n = 1 ∨ n = 2 := by
  constructor
  · intro h
    match n with
    | 0 => have := h 0 0 0; simp [lad_zero] at this
    | 1 => exact Or.inl rfl
    | 2 => exact Or.inr rfl
    | n + 3 =>
      have := h 1 0 1
      rw [lad_zero_right (n + 1) 1, lad_zero_left, lad_zero_right (n + 1) 1,
        lad_one_one (n + 2)] at this
      exact absurd this (by decide)
  · rintro (rfl | rfl)
    · intro a b c; simp only [lad_eq_sadd]; exact sadd_assoc' a b c
    · intro a b c; simp only [lad_eq_circ]; exact circ_assoc a b c

/-- Only rung 2 has a two-sided identity in the whole numbers (it is 0). Rung 1's identity
    is −1, outside Nat. -/
theorem lad_identity_iff (n : Nat) : HasIdentity (lad n) ↔ n = 2 := by
  constructor
  · rintro ⟨e, he⟩
    match n with
    | 0 => have := (he 0).2; simp [lad_zero] at this
    | 1 => have := (he 0).1; rw [lad_one] at this; omega
    | 2 => rfl
    | n + 3 =>
      have h0 : e = 0 := (lad_zero_right (n + 1) e).symm.trans (he 0).1
      have h1 := (he 1).1
      rw [h0, lad_zero_left] at h1
      exact absurd h1 (by decide)
  · rintro rfl
    exact ⟨0, fun a => ⟨by rw [lad_eq_circ]; exact zero_circ a, by rw [lad_eq_circ]; exact circ_zero a⟩⟩

/-- 0 is a right identity exactly from rung 2 up. -/
theorem lad_right_id_iff (n : Nat) : (∀ a, lad n a 0 = a) ↔ 2 ≤ n := by
  constructor
  · intro h
    match n with
    | 0 => have := h 0; simp [lad_zero] at this
    | 1 => have := h 0; simp [lad_one_zero] at this
    | n + 2 => omega
  · intro hn a
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
    rfl

/-- 0 is a left absorber exactly from rung 3 up. -/
theorem lad_left_absorb_iff (n : Nat) : (∀ b, lad n 0 b = 0) ↔ 3 ≤ n := by
  constructor
  · intro h
    match n with
    | 0 => have := h 0; simp [lad_zero] at this
    | 1 => have := h 0; simp [lad_one_zero] at this
    | 2 => have := h 1; rw [lad_eq_circ] at this; unfold circ at this; omega
    | n + 3 => omega
  · intro hn b
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 3 := ⟨n - 3, by omega⟩
    exact lad_zero_left m b

/-! ### A6. Distributivity between neighbours -/

/-- Rung 2 distributes over rung 1 on the left … -/
theorem circ_sadd_left (a b c : Nat) : circ a (sadd b c) = sadd (circ a b) (circ a c) := by
  unfold circ sadd; grind
/-- … and on the right. -/
theorem circ_sadd_right (a b c : Nat) : circ (sadd a b) c = sadd (circ a c) (circ b c) := by
  unfold circ sadd; grind

/-- Rung 3 distributes over rung 2 on the right ((xy)^z = x^z y^z, fused). -/
theorem fexp_circ_right' (a b c : Nat) : fexp (circ a b) c = circ (fexp a c) (fexp b c) :=
  fexp_circ_right a b c

/-- … but not on the left (x^(yz) ≠ x^y · x^z). -/
theorem fexp_circ_left_fails : fexp 1 (circ 1 2) ≠ circ (fexp 1 1) (fexp 1 2) := by decide

/-- Rung 4 distributes over rung 3 on neither side. -/
theorem ftet_fexp_right_fails : ftet (fexp 1 2) 1 ≠ fexp (ftet 1 1) (ftet 2 1) := by decide
theorem ftet_fexp_left_fails : ftet 1 (fexp 1 1) ≠ fexp (ftet 1 1) (ftet 1 1) := by decide

/-- The obvious witnesses for rung n + 1 over rung n tell nothing: at 1, 1, 1 both sides agree on
    every rung (each side is the corner 3 pushed up one rung). -/
theorem lad_dist_right_at_ones (n : Nat) :
    lad (n + 3) (lad (n + 2) 1 1) 1 = lad (n + 2) (lad (n + 3) 1 1) (lad (n + 3) 1 1) := by
  rw [lad_one_one, lad_one_one, lad_diag]

/-- Left distributivity of rung n + 1 over rung n holds exactly at the bottom two steps
    (rung 1 over successor, ◯ over rung 1). From rung 3 over ◯ upward it fails, and one
    witness serves every rung: 1 ⋆ (0 ⋆' 0) = 1 but (1 ⋆ 0) ⋆' (1 ⋆ 0) = 3. -/
theorem lad_left_dist_iff (n : Nat) :
    (∀ a b c, lad (n + 1) a (lad n b c) = lad n (lad (n + 1) a b) (lad (n + 1) a c)) ↔ n ≤ 1 := by
  constructor
  · intro h
    match n with
    | 0 => omega
    | 1 => omega
    | m + 2 =>
      have := h 1 0 0
      rw [lad_zero_right m 0, lad_zero_right (m + 1) 1, lad_one_one (m + 1)] at this
      exact absurd this (by decide)
  · intro hn
    match n, hn with
    | 0, _ => intro a b c; rfl
    | 1, _ =>
      intro a b c
      show lad 2 a (lad 1 b c) = lad 1 (lad 2 a b) (lad 2 a c)
      simp only [lad_eq_circ, lad_eq_sadd]; exact circ_sadd_left a b c

/-- Right distributivity of rung n + 1 over rung n holds at the bottom three steps. -/
theorem lad_right_dist_le_two (n : Nat) (hn : n ≤ 2) :
    ∀ a b c, lad (n + 1) (lad n a b) c = lad n (lad (n + 1) a c) (lad (n + 1) b c) := by
  match n, hn with
  | 0, _ =>
    intro a b c
    show lad 1 (lad 0 a b) c = lad 0 (lad 1 a c) (lad 1 b c)
    simp only [lad_eq_sadd, lad_zero]; unfold sadd; omega
  | 1, _ =>
    intro a b c
    show lad 2 (lad 1 a b) c = lad 1 (lad 2 a c) (lad 2 b c)
    simp only [lad_eq_circ, lad_eq_sadd]; exact circ_sadd_right a b c
  | 2, _ =>
    intro a b c
    show lad 3 (lad 2 a b) c = lad 2 (lad 3 a c) (lad 3 b c)
    simp only [lad_eq_fexp, lad_eq_circ]; exact fexp_circ_right' a b c

/-- …and fails at the fourth (rung 4 over rung 3). Whether it fails at every higher rung is
    open in Lean; the checks show it failing at rung 5 over rung 4 by tower height. -/
theorem lad_right_dist_three_fails :
    ¬ ∀ a b c, lad 4 (lad 3 a b) c = lad 3 (lad 4 a c) (lad 4 b c) := by
  intro h; have := h 1 2 1
  rw [lad_eq_ftet, lad_eq_fexp, lad_eq_ftet, lad_eq_ftet, lad_eq_fexp] at this
  exact ftet_fexp_right_fails this

/-! ### A6b. Right distributivity fails on every rung from 4 up: a mod-7 argument

Rung 4 over rung 3 fails by direct computation. For the giant rungs the numbers cannot be
computed, but their remainders can. With the witness a, b, c = 1, 2, 1 (Goodstein's 2, 3, 2) the
two sides become H_n(P, P) with P = H_n(2, 3), and H_n(4, H_n(3, 3)). Every H_n(2, y) with y ≥ 3
leaves 16 on division by 42, so P leaves 2 on division by 7. A base that leaves 1 on division by
3 and whose cube leaves 1 on division by 7 passes its remainder mod 7 up every rung from 4. So
the left side leaves 2 and the right side leaves 4: they differ. -/

theorem two_pow_mod42 : ∀ k, 2 ^ (4 + 6 * k) % 42 = 16
  | 0 => by decide
  | k + 1 => by
    rw [show 4 + 6 * (k + 1) = (4 + 6 * k) + 6 by omega, Nat.pow_add, Nat.mul_mod,
      two_pow_mod42 k]

theorem two_pow_mod42' (z : Nat) (h : z % 6 = 4) : 2 ^ z % 42 = 16 := by
  obtain ⟨k, rfl⟩ : ∃ k, z = 4 + 6 * k := ⟨z / 6, by omega⟩
  exact two_pow_mod42 k

/-- 2 ↑↑ (z + 3) leaves 16 on division by 42. -/
theorem tet_two_mod42 : ∀ z, hyp 4 2 (z + 3) % 42 = 16
  | 0 => by decide
  | z + 1 => by
    have ih := tet_two_mod42 z
    show hyp 3 2 (hyp 4 2 (z + 3)) % 42 = 16
    rw [hyp_pow]; exact two_pow_mod42' _ (by omega)

/-- On every rung from 4 up, H_n(2, y) leaves 16 on division by 42 once y ≥ 3. -/
theorem hyp_two_mod42 : ∀ m z, hyp (m + 4) 2 (z + 3) % 42 = 16
  | 0, z => tet_two_mod42 z
  | m + 1, 0 => by
    show hyp (m + 4) 2 (hyp (m + 5) 2 2) % 42 = 16
    rw [hyp_two_two (m + 4)]; exact hyp_two_mod42 m 1
  | m + 1, z + 1 => by
    have ihz : hyp (m + 5) 2 (z + 3) % 42 = 16 := hyp_two_mod42 (m + 1) z
    show hyp (m + 4) 2 (hyp (m + 5) 2 (z + 3)) % 42 = 16
    obtain ⟨w, hw⟩ : ∃ w, hyp (m + 5) 2 (z + 3) = w + 3 := ⟨hyp (m + 5) 2 (z + 3) - 3, by omega⟩
    rw [hw]; exact hyp_two_mod42 m w

/-- A base that leaves 1 mod 3 keeps leaving 1 mod 3 on every rung from 3 up. -/
theorem hyp_mod3 : ∀ m x y, x % 3 = 1 → hyp (m + 3) x y % 3 = 1
  | 0, x, y, hx => by rw [hyp_pow, Nat.pow_mod, hx, Nat.one_pow]
  | m + 1, x, 0, _ => rfl
  | m + 1, x, y + 1, hx => by
    show hyp (m + 3) x (hyp (m + 4) x y) % 3 = 1
    exact hyp_mod3 m x _ hx

theorem pow_mod7_cycle (x : Nat) (h3 : x ^ 3 % 7 = 1) : ∀ q, x ^ (3 * q + 1) % 7 = x % 7
  | 0 => by rw [Nat.mul_zero, Nat.zero_add, Nat.pow_one]
  | q + 1 => by
    rw [show 3 * (q + 1) + 1 = (3 * q + 1) + 3 by omega, Nat.pow_add, Nat.mul_mod,
      pow_mod7_cycle x h3 q, h3, Nat.mul_one, Nat.mod_mod]

/-- A base x with x ≡ 1 (mod 3) and x³ ≡ 1 (mod 7) passes its remainder mod 7 up every rung from 4. -/
theorem hyp_mod7 : ∀ m x y, x % 3 = 1 → x ^ 3 % 7 = 1 → hyp (m + 4) x (y + 1) % 7 = x % 7
  | 0, x, y, hx, h3 => by
    show hyp 3 x (hyp 4 x y) % 7 = x % 7
    rw [hyp_pow]
    have ht : hyp 4 x y % 3 = 1 := hyp_mod3 1 x y hx
    obtain ⟨q, hq⟩ : ∃ q, hyp 4 x y = 3 * q + 1 := ⟨hyp 4 x y / 3, by omega⟩
    rw [hq]; exact pow_mod7_cycle x h3 q
  | m + 1, x, 0, _, _ => by
    show hyp (m + 4) x (hyp (m + 5) x 0) % 7 = x % 7
    rw [show hyp (m + 5) x 0 = 1 from rfl, hyp_one (m + 2) x]
  | m + 1, x, y + 1, hx, h3 => by
    have ihy := hyp_mod7 (m + 1) x y hx h3
    have hx0 : x % 7 ≠ 0 := by
      intro h0; rw [Nat.pow_mod, h0] at h3; exact absurd h3 (by decide)
    show hyp (m + 4) x (hyp (m + 5) x (y + 1)) % 7 = x % 7
    obtain ⟨w, hw⟩ : ∃ w, hyp (m + 5) x (y + 1) = w + 1 :=
      ⟨hyp (m + 5) x (y + 1) - 1, by
        have : hyp (m + 1 + 4) x (y + 1) % 7 = x % 7 := ihy
        have : hyp (m + 5) x (y + 1) ≠ 0 := by intro hz; simp_all
        omega⟩
    rw [hw]; exact hyp_mod7 m x w hx h3

theorem hyp_pos : ∀ m x y, 1 ≤ x → 1 ≤ hyp (m + 3) x y
  | 0, x, y, hx => by rw [hyp_pow]; exact Nat.pow_pos hx
  | m + 1, x, 0, _ => Nat.le_refl 1
  | m + 1, x, y + 1, hx => by
    show 1 ≤ hyp (m + 3) x (hyp (m + 4) x y)
    exact hyp_pos m x _ hx

/-- Goodstein's ladder: from rung 5 over rung 4 upward, right distributivity fails at 2, 3, 2. -/
theorem hyp_right_dist_fails (m : Nat) :
    hyp (m + 5) (hyp (m + 4) 2 3) 2 ≠ hyp (m + 4) (hyp (m + 5) 2 2) (hyp (m + 5) 3 2) := by
  intro h
  have hd : ∀ x, hyp (m + 5) x 2 = hyp (m + 4) x x := fun x => by
    show hyp (m + 4) x (hyp (m + 5) x 1) = _
    rw [hyp_one (m + 3) x]
  rw [hd, hd, hd, hyp_two_two (m + 3)] at h
  have hP : hyp (m + 4) 2 3 % 42 = 16 := hyp_two_mod42 m 0
  obtain ⟨p, hp⟩ : ∃ p, hyp (m + 4) 2 3 = p + 1 := ⟨hyp (m + 4) 2 3 - 1, by omega⟩
  have hQ : 1 ≤ hyp (m + 4) 3 3 := hyp_pos (m + 1) 3 3 (by decide)
  obtain ⟨q, hq⟩ : ∃ q, hyp (m + 4) 3 3 = q + 1 := ⟨hyp (m + 4) 3 3 - 1, by omega⟩
  have hP3 : hyp (m + 4) 2 3 % 3 = 1 := by omega
  have hP7 : hyp (m + 4) 2 3 ^ 3 % 7 = 1 := by
    rw [Nat.pow_mod, show hyp (m + 4) 2 3 % 7 = 2 by omega]
  have L := hyp_mod7 m _ p hP3 hP7
  rw [← hp] at L
  have R := hyp_mod7 m 4 q (by decide) (by decide)
  rw [← hq] at R
  rw [h] at L
  omega

/-- The same, read on the fused ladder: rung n + 1 over rung n, right side, fails at 1, 2, 1 for n ≥ 4. -/
theorem lad_right_dist_fails_high (m : Nat) :
    lad (m + 5) (lad (m + 4) 1 2) 1 ≠ lad (m + 4) (lad (m + 5) 1 1) (lad (m + 5) 2 1) := by
  intro h
  apply hyp_right_dist_fails m
  have e1 := lad_lift (m + 5) (lad (m + 4) 1 2) 1
  rw [lad_lift (m + 4) 1 2] at e1
  have e2 := lad_lift (m + 4) (lad (m + 5) 1 1) (lad (m + 5) 2 1)
  rw [lad_lift (m + 5) 1 1, lad_lift (m + 5) 2 1] at e2
  rw [← e1, ← e2, h]

/-- Right distributivity of rung n + 1 over rung n holds exactly for n = 0, 1, 2. -/
theorem lad_right_dist_iff (n : Nat) :
    (∀ a b c, lad (n + 1) (lad n a b) c = lad n (lad (n + 1) a c) (lad (n + 1) b c)) ↔ n ≤ 2 := by
  constructor
  · intro h
    match n with
    | 0 => omega
    | 1 => omega
    | 2 => omega
    | 3 => exact absurd h lad_right_dist_three_fails
    | m + 4 => exact absurd (h 1 2 1) (lad_right_dist_fails_high m)
  · exact lad_right_dist_le_two n

/-! ### A7. Transport, and what it rules out -/

/-- Transport of an operation along a pair of maps: a ⋆ b = g(f a ∘ f b). -/
def transport {β : Type} (f : Nat → β) (g : β → Nat) (op : β → β → β) (a b : Nat) : Nat :=
  g (op (f a) (f b))

/-- Every rung is the ordinary hyperoperation transported by x ↦ x + 1 and back. -/
theorem lad_eq_transport (n : Nat) : lad n = transport (· + 1) (· - 1) (hyp n) := by
  funext a b; show lad n a b = hyp n (a + 1) (b + 1) - 1; have := lad_lift n a b; omega

/-- Transport keeps commutativity, whatever the maps. -/
theorem transport_comm {β : Type} (f : Nat → β) (g : β → Nat) (op : β → β → β)
    (h : ∀ x y, op x y = op y x) : IsComm (transport f g op) := by
  intro a b; unfold transport; rw [h]

/-- Transport keeps associativity when the maps undo each other on the values that occur. -/
theorem transport_assoc {β : Type} (f : Nat → β) (g : β → Nat) (op : β → β → β)
    (hfg : ∀ a b, f (g (op (f a) (f b))) = op (f a) (f b))
    (h : ∀ x y z, op (op x y) z = op x (op y z)) : IsAssoc (transport f g op) := by
  intro a b c; unfold transport; rw [hfg, hfg, h]

/-- No rung other than 1 and 2 is a transport of a commutative operation, by any maps at all.
    In particular no height of a commutative ladder (Bennett's, or B_s) is fused exponentiation
    or anything above it. -/
theorem not_transport_of_comm (n : Nat) (hn : n = 0 ∨ 3 ≤ n) {β : Type} (f : Nat → β)
    (g : β → Nat) (op : β → β → β) (h : ∀ x y, op x y = op y x) :
    lad n ≠ transport f g op := by
  intro heq
  have hc : IsComm (lad n) := by rw [heq]; exact transport_comm f g op h
  rcases (lad_comm_iff n).1 hc with h1 | h2 <;> omega

/-! ## Part B. The height ladder, for any flow -/

section Heights

open Lean.Grind

variable {S X : Type} [CommRing S] [Add X] [Zero X]

/-- A height flow: E_0 = id, E_s ∘ E_t = E_(s+t), and 0 is fixed. For E = eˣ − 1 on [0, ∞)
    such a flow exists by Szekeres (1958); this file only uses the three laws. -/
structure Flow (S X : Type) [CommRing S] [Add X] [Zero X] where
  E : S → X → X
  E_zero : ∀ x, E 0 x = x
  E_add : ∀ s t x, E s (E t x) = E (s + t) x
  E_fix : ∀ s, E s 0 = 0

namespace Flow

variable (F : Flow S X)

/-- The height-s operation a +⁽ˢ⁾ b = E_s(E_(−s) a + E_(−s) b) (B_s in fractional-rungs/). -/
def B (s : S) (a b : X) : X := F.E s (F.E (-s) a + F.E (-s) b)

theorem E_cancel (s : S) (x : X) : F.E s (F.E (-s) x) = x := by
  rw [F.E_add, show s + -s = (0 : S) by grind, F.E_zero]

theorem E_cancel' (s : S) (x : X) : F.E (-s) (F.E s x) = x := by
  rw [F.E_add, show -s + s = (0 : S) by grind, F.E_zero]

/-- Height 0 is addition. -/
theorem B_zero_height (a b : X) : F.B 0 a b = a + b := by
  unfold B; rw [show (-0 : S) = 0 by grind, F.E_zero, F.E_zero, F.E_zero]

/-- Every height is commutative. -/
theorem B_comm (hc : ∀ x y : X, x + y = y + x) (s : S) (a b : X) : F.B s a b = F.B s b a := by
  unfold B; rw [hc]

/-- Every height is associative. -/
theorem B_assoc (ha : ∀ x y z : X, x + y + z = x + (y + z)) (s : S) (a b c : X) :
    F.B s (F.B s a b) c = F.B s a (F.B s b c) := by
  unfold B; rw [F.E_cancel', F.E_cancel', ha]

/-- 0 is the identity at every height. -/
theorem B_id (hz : ∀ x : X, x + 0 = x) (s : S) (a : X) : F.B s a 0 = a := by
  unfold B; rw [F.E_fix, hz, F.E_cancel]

/-- Heights climb by transport: B_(s+t) = (E_s)₊(B_t). -/
theorem B_climb (s t : S) (a b : X) : F.B (s + t) a b = F.E s (F.B t (F.E (-s) a) (F.E (-s) b)) := by
  unfold B
  rw [F.E_add, F.E_add, F.E_add, show -t + -s = -(s + t) by grind]

/-! ### The hypernatural pair, in shifted coordinates

  ln_s(1 + x) = E_(−s)(x) and exp_s(t) = 1 + E_s(t). Here `lnh s x` is ln_s(1 + x) and
  `exph s t` is exp_s(t) − 1, so the statements need no subtraction on X. -/

/-- The hypernatural logarithm, ln_s(1 + x). -/
def lnh (s : S) (x : X) : X := F.E (-s) x
/-- The hypernatural exponential, exp_s(t) − 1. -/
def exph (s : S) (t : X) : X := F.E s t

theorem lnh_exph (s : S) (t : X) : F.lnh s (F.exph s t) = t := F.E_cancel' s t
theorem exph_lnh (s : S) (x : X) : F.exph s (F.lnh s x) = x := F.E_cancel s x
/-- ln_s(1) = 0. -/
theorem lnh_zero (s : S) : F.lnh s 0 = 0 := F.E_fix (-s)
/-- Heights compose: ln_(s+t)(X) = ln_s(1 + ln_t(X)). -/
theorem lnh_compose (s t : S) (x : X) : F.lnh (s + t) x = F.lnh s (F.lnh t x) := by
  unfold lnh; rw [F.E_add, show -s + -t = -(s + t) by grind]

/-- If E_1 carries + to a circle operation ◯ (as eᵘ⁺ᵛ − 1 = (eᵘ − 1) ◯ (eᵛ − 1)),
    then height 1 is ◯. -/
theorem B_one (circ' : X → X → X) (hexp : ∀ p q, circ' (F.E 1 p) (F.E 1 q) = F.E 1 (p + q))
    (a b : X) : F.B 1 a b = circ' a b := by
  conv => rhs; rw [← F.E_cancel 1 a, ← F.E_cancel 1 b]
  rw [hexp]; rfl

/-- The hypernatural log turns ◯ (ordinary × in lifted coordinates) into the height-(1 − s)
    operation: ln_s(XY) = ln_s X +⁽¹⁻ˢ⁾ ln_s Y. -/
theorem lnh_circ (circ' : X → X → X) (hexp : ∀ p q, circ' (F.E 1 p) (F.E 1 q) = F.E 1 (p + q))
    (s : S) (x y : X) : F.lnh s (circ' x y) = F.B (1 - s) (F.lnh s x) (F.lnh s y) := by
  rw [← F.B_one circ' hexp x y]
  unfold lnh B
  rw [F.E_add, F.E_add, F.E_add, show -(1 - s) + -s = (-1 : S) by grind,
    show -s + 1 = (1 - s : S) by grind]

/-- Height s + 1 distributes over its partner (E_s)₊(r), where r is the operation ◯ distributes
    over (rung 1, a + b + 1, for the circle operation). Neighbouring heights do not distribute
    (`circ_not_dist_add`), so this partner is the right second operation at every height. -/
theorem B_partner_dist (circ' r : X → X → X)
    (hexp : ∀ p q, circ' (F.E 1 p) (F.E 1 q) = F.E 1 (p + q))
    (hd : ∀ a b c, circ' a (r b c) = r (circ' a b) (circ' a c)) (s : S) (a b c : X) :
    F.B (s + 1) a (F.E s (r (F.E (-s) b) (F.E (-s) c))) =
      F.E s (r (F.E (-s) (F.B (s + 1) a b)) (F.E (-s) (F.B (s + 1) a c))) := by
  simp only [F.B_climb s 1, F.B_one circ' hexp, F.E_cancel', hd]

end Flow

/-- Neighbouring heights do not distribute: ◯ (height 1) over + (height 0) fails. -/
theorem circ_not_dist_add : ¬ ∀ a b c : Nat, circ a (b + c) = circ a b + circ a c := by
  intro h; have := h 1 1 1; unfold circ at this; omega

/-- Rung 1 (a + b + 1) is no height of any flow: every height has identity 0, and rung 1 does not. -/
theorem sadd_not_height {S X : Type} [CommRing S] [CommRing X] (F : Flow S X)
    (h10 : (1 : X) ≠ 0) (s : S) : ¬ ∀ a b : X, F.B s a b = a + b + 1 := by
  intro h
  have h0 := h 0 0
  rw [F.B_id (Lean.Grind.Semiring.add_zero) s 0] at h0
  rw [Lean.Grind.Semiring.add_zero, Lean.Grind.Semiring.add_comm,
    Lean.Grind.Semiring.add_zero] at h0
  exact h10 h0.symm

end Heights

/-! ## Part C. The dial k -/

section Dial

open Lean.Grind

variable {K : Type} [Field K]

/-- Rung 1 at dial k: transport of + by the lift 1 + kx gives a + b + 1/k (Pap's ⊕ for g = 1 + kx). -/
def saddK (k a b : K) : K := a + b + 1 / k

/-- Its lift: 1 + k(a ⊕ₖ b) = (1 + ka) + (1 + kb). -/
theorem saddK_lift (k a b : K) (hk : k ≠ 0) : 1 + k * saddK k a b = (1 + k * a) + (1 + k * b) := by
  unfold saddK; grind

/-- Rung 1 at dial k is rung 1 rescaled: k(a ⊕ₖ b) = (ka) + (kb) + 1. -/
theorem saddK_scale (k a b : K) (hk : k ≠ 0) : k * saddK k a b = k * a + k * b + 1 := by
  unfold saddK; grind

/-- Rung 2 at dial k distributes over rung 1 at dial k, for every k ≠ 0 (Pap's generated semiring). -/
theorem circK_saddK (k a b c : K) (hk : k ≠ 0) :
    circK k a (saddK k b c) = saddK k (circK k a b) (circK k a c) := by
  unfold circK saddK; grind

/-- Rung 3 at dial k, at the whole-number lift 1 + kb = m: ((1 + ka)^m − 1)/k. -/
def fexpK (k a : K) (m : Nat) : K := ((1 + k * a) ^ m - 1) / k

theorem fexpK_lift (k a : K) (m : Nat) (hk : k ≠ 0) : 1 + k * fexpK k a m = (1 + k * a) ^ m := by
  unfold fexpK; grind

/-- Rung 3 is rung 2 iterated, at every dial k ≠ 0. -/
theorem fexpK_succ (k a : K) (m : Nat) (hk : k ≠ 0) :
    fexpK k a (m + 1) = circK k a (fexpK k a m) := by
  unfold fexpK circK; rw [Lean.Grind.Semiring.pow_succ]; grind

/-- Rung 3 at dial k is the k = 1 rung rescaled. -/
theorem fexpK_scale (k a : K) (m : Nat) (hk : k ≠ 0) :
    k * fexpK k a m = (1 + k * a) ^ m - 1 := by
  unfold fexpK; grind

/-- At k = 0 the lift 1 + kx sends everything to 1, so nothing can be transported along it. -/
theorem lift_zero_not_injective : ¬ ∀ x y : K, (1 + 0 * x = 1 + 0 * y → x = y) := by
  intro h
  have := h 0 1 (by grind)
  exact (show (0 : K) ≠ 1 by grind) this

end Dial

end HyperLadder
