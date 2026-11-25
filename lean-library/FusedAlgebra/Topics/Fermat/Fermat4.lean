/-
  Fermat's Last Theorem for exponent 4
  ====================================

  What is proved (core Lean 4, no mathlib, no `sorry`, standard axioms only):

  * `no_x4_y4_z2` : if x and y are nonzero naturals, x⁴ + y⁴ is never a perfect square.
  * `fermat4`     : so x⁴ + y⁴ = z⁴ has no solution in nonzero naturals (z⁴ is the square of z²).
  * `fourth_power_ne_mpow_four` : the same fact said with mult-add. Since 1 + s^(⊞4) = (1 + s)⁴,
    "q⁴ is the mult-add 4th power of s" means 1 + q⁴ = (1 + s)⁴. No nonzero rational q and
    rational s do this, i.e. no nonzero rational fourth power is a mult-add fourth power.

  The route is Fermat's own infinite descent, by strong induction on z:

  1. If x, y share a factor d ≥ 2, divide it out: d⁴ | z² gives d² | z, and we get a smaller z.
  2. Otherwise one of x, y is odd (they are not both even); call it x.
     (x², y², z) is then a primitive Pythagorean triple, so (by `pyth_param`)
     x² = m² − n², y² = 2mn, z = m² + n² with m, n coprime.
  3. (x, n, m) is again a primitive Pythagorean triple: x = r² − s², n = 2rs, m = r² + s².
  4. Now y² = 4·m·r·s with m, r, s pairwise coprime, so each is a square
     (`sq_of_coprime_mul`, proved without prime factorization by dividing out gcd(a, c)):
     m = c², r = a², s = b².
  5. Then c² = m = r² + s² = a⁴ + b⁴ is a new solution with 0 < c ≤ m < m² + n² = z.
     A smallest solution cannot exist, so there is none.

  Parity facts come from "a square is 0 or 1 mod 4" (`sq_mod_four`); the ring arithmetic is
  done by `grind` and the linear bookkeeping by `omega`.

  LIBRARY NOTE: copied from fermat/lean/Fermat4.lean (Nov 25) as Topics/Fermat/Fermat4.
-/
import FusedAlgebra.Base

open MultAdd

namespace Fermat4

/-- Squaring is one-to-one on the natural numbers: a² = b² gives a = b. -/
theorem sq_inj {a b : Nat} (h : a^2 = b^2) : a = b :=
  Nat.le_antisymm ((Nat.pow_le_pow_iff_left (by decide)).1 (Nat.le_of_eq h))
    ((Nat.pow_le_pow_iff_left (by decide)).1 (Nat.le_of_eq h.symm))

/-- Every square is 0 or 1 mod 4: 0 for an even number, 1 for an odd one. -/
theorem sq_mod_four (n : Nat) : (n % 2 = 0 ∧ n^2 % 4 = 0) ∨ (n % 2 = 1 ∧ n^2 % 4 = 1) := by
  obtain ⟨k, hk | hk⟩ : ∃ k, n = 2*k ∨ n = 2*k+1 := ⟨n/2, by omega⟩
  · left; subst hk
    have : (2*k)^2 = 4*(k*k) := by grind
    rw [this]; omega
  · right; subst hk
    have : (2*k+1)^2 = 4*(k*k+k)+1 := by grind
    rw [this]; omega

/-- A number divides its own square. -/
theorem dvd_sq_self (n : Nat) : n ∣ n^2 := ⟨n, Nat.pow_two n⟩

/-- If a·b is a square and a, b share no factor, then a is itself a square. -/
theorem sq_of_coprime_mul {a b c : Nat} (hab : Nat.Coprime a b) (h : a * b = c^2) :
    ∃ k, a = k^2 := by
  rcases Nat.eq_zero_or_pos a with ha | ha
  · exact ⟨0, by simp [ha]⟩
  obtain ⟨a', c', hco, ha', hc'⟩ := Nat.exists_coprime a c
  have hg : 0 < Nat.gcd a c := Nat.gcd_pos_of_pos_left c ha
  generalize Nat.gcd a c = g at ha' hc' hg
  subst ha' hc'
  have h1 : a' * b = c'^2 * g := by
    apply Nat.eq_of_mul_eq_mul_left hg
    grind
  have h2 : a' ∣ g := Nat.Coprime.dvd_of_dvd_mul_left (hco.pow_right 2) ⟨b, h1.symm⟩
  obtain ⟨t, rfl⟩ := h2
  have ha'pos : 0 < a' := Nat.pos_of_ne_zero (by intro h0; subst h0; simp at ha)
  have hb : b = c'^2 * t := by
    apply Nat.eq_of_mul_eq_mul_left ha'pos
    grind
  have ht : t = 1 := by
    have d1 : t ∣ a' * (a' * t) := ⟨a' * a', by grind⟩
    have d2 : t ∣ b := ⟨c'^2, by rw [hb]; grind⟩
    have := Nat.dvd_gcd d1 d2
    rw [hab.gcd_eq_one] at this
    exact Nat.eq_one_of_dvd_one this
  subst ht
  exact ⟨a', by grind⟩

/-- Primitive Pythagorean triples: a² + b² = c² with a odd and gcd(a,b) = 1 gives coprime m, n with
    a = m² − n² (written a + n² = m²), b = 2mn, c = m² + n². -/
theorem pyth_param {a b c : Nat} (h : a^2 + b^2 = c^2) (ha : a % 2 = 1) (hab : Nat.Coprime a b) :
    ∃ m n, a + n^2 = m^2 ∧ c = m^2 + n^2 ∧ b = 2*m*n ∧ Nat.Coprime m n := by
  have hb : b % 2 = 0 ∧ c % 2 = 1 := by
    rcases sq_mod_four a with ⟨_, h1⟩ | ⟨_, h1⟩ <;> rcases sq_mod_four b with ⟨_, h2⟩ | ⟨_, h2⟩ <;>
      rcases sq_mod_four c with ⟨_, h3⟩ | ⟨_, h3⟩ <;> omega
  have hac : a ≤ c := (Nat.pow_le_pow_iff_left (n := 2) (by decide)).1 (by omega)
  obtain ⟨u, hu⟩ : ∃ u, c = a + 2*u := ⟨(c-a)/2, by omega⟩
  obtain ⟨w, hw⟩ : ∃ w, b = 2*w := ⟨b/2, by omega⟩
  subst hu hw
  have huv : u * (a + u) = w^2 := by
    have : 4 * (u * (a+u)) = 4 * w^2 := by grind
    omega
  have hco : Nat.Coprime u (a+u) := by
    show Nat.gcd u (a+u) = 1
    have d1 : Nat.gcd u (a+u) ∣ a := by
      have := Nat.dvd_sub (Nat.gcd_dvd_right u (a+u)) (Nat.gcd_dvd_left u (a+u))
      rwa [Nat.add_sub_cancel] at this
    have d2 : Nat.gcd u (a+u) ∣ (2*w)^2 := by
      have e : (2*w)^2 = (u * (a+u)) * 4 := by grind
      rw [e]
      exact Nat.dvd_mul_right_of_dvd (Nat.dvd_mul_right_of_dvd (Nat.gcd_dvd_left _ _) _) _
    have := Nat.dvd_gcd d1 d2
    rw [(hab.pow_right 2).gcd_eq_one] at this
    exact Nat.eq_one_of_dvd_one this
  obtain ⟨n, hn⟩ := sq_of_coprime_mul hco huv
  obtain ⟨m, hm⟩ := sq_of_coprime_mul hco.symm (by rw [Nat.mul_comm]; exact huv)
  refine ⟨m, n, by omega, by omega, ?_, ?_⟩
  · have : w = m * n := sq_inj (by rw [← huv, hm, hn]; grind)
    rw [this]; grind
  · rw [hm, hn] at hco
    exact (hco.symm.coprime_dvd_left (dvd_sq_self m)).coprime_dvd_right (dvd_sq_self n)

/-- One step of Fermat's descent: a coprime solution of x⁴ + y⁴ = z² (x odd, y ≠ 0) yields a
    nonzero solution a⁴ + b⁴ = c² with c < z. -/
theorem descent_step {x y z : Nat} (hx : x % 2 = 1) (hy : y ≠ 0) (hxy : Nat.Coprime x y)
    (h : x^4 + y^4 = z^2) : ∃ a b c, a ≠ 0 ∧ b ≠ 0 ∧ c < z ∧ a^4 + b^4 = c^2 := by
  have hx2 : (x^2) % 2 = 1 := by rcases sq_mod_four x with ⟨_, _⟩ | ⟨_, _⟩ <;> omega
  obtain ⟨m, n, h1, h2, h3, hmn⟩ :=
    pyth_param (a := x^2) (b := y^2) (c := z) (by grind) hx2 (Nat.Coprime.pow 2 2 hxy)
  have hxn : Nat.Coprime x n := by
    show Nat.gcd x n = 1
    have d1 : Nat.gcd x n ∣ m^2 := by
      rw [← h1]
      exact Nat.dvd_add (Nat.dvd_trans (Nat.gcd_dvd_left _ _) (dvd_sq_self x))
        (Nat.dvd_trans (Nat.gcd_dvd_right _ _) (dvd_sq_self n))
    have := Nat.dvd_gcd (Nat.gcd_dvd_right x n) d1
    rw [(hmn.symm.pow_right 2).gcd_eq_one] at this
    exact Nat.eq_one_of_dvd_one this
  obtain ⟨r, s, h4, h5, h6, hrs⟩ := pyth_param (a := x) (b := n) (c := m) h1 hx hxn
  subst h6
  have hy2 : y^2 = 4 * (m * (r * s)) := by rw [h3]; grind
  have hyev : y % 2 = 0 := by rcases sq_mod_four y with ⟨_, _⟩ | ⟨_, _⟩ <;> omega
  obtain ⟨v, rfl⟩ : ∃ v, y = 2 * v := ⟨y/2, by omega⟩
  have hv : m * (r * s) = v^2 := by
    have : 4 * (m * (r * s)) = 4 * v^2 := by grind
    omega
  have hmrs : Nat.Coprime m (r * s) := hmn.coprime_dvd_right (show r * s ∣ 2 * r * s from ⟨2, by grind⟩)
  obtain ⟨c, hc⟩ := sq_of_coprime_mul hmrs hv
  obtain ⟨e, he⟩ := sq_of_coprime_mul hmrs.symm (by rw [Nat.mul_comm]; exact hv)
  obtain ⟨a, ha⟩ := sq_of_coprime_mul hrs he
  obtain ⟨b, hb⟩ := sq_of_coprime_mul hrs.symm (by rw [Nat.mul_comm]; exact he)
  have hr0 : r ≠ 0 := by
    intro h0; subst h0; apply hy
    have : v = 0 := sq_inj (by rw [← hv]; simp)
    omega
  have hs0 : s ≠ 0 := by
    intro h0; subst h0; apply hy
    have : v = 0 := sq_inj (by rw [← hv]; simp)
    omega
  refine ⟨a, b, c, ?_, ?_, ?_, ?_⟩
  · intro h0; subst h0; simp at ha; exact hr0 ha
  · intro h0; subst h0; simp at hb; exact hs0 hb
  · have hn : 0 < (2 * r * s)^2 := Nat.pow_pos (by
      have := Nat.pos_of_ne_zero hr0; have := Nat.pos_of_ne_zero hs0
      exact Nat.mul_pos (Nat.mul_pos (by decide) ‹0 < r›) ‹0 < s›)
    have c1 : c ≤ m := by rw [hc, Nat.pow_two]; exact Nat.le_mul_self c
    have c2 : m ≤ m^2 := by rw [Nat.pow_two]; exact Nat.le_mul_self m
    omega
  · rw [← hc, h5, ha, hb]; grind

end Fermat4

open Fermat4

/-- No sum of two nonzero fourth powers is a square. -/
theorem no_x4_y4_z2 (x y z : Nat) (hx : x ≠ 0) (hy : y ≠ 0) : x^4 + y^4 ≠ z^2 := by
  induction z using Nat.strongRecOn generalizing x y with
  | _ z ih =>
  intro h
  by_cases hco : Nat.Coprime x y
  · have hodd : x % 2 = 1 ∨ y % 2 = 1 := by
      refine Classical.byContradiction fun hc => ?_
      have h2 : 2 ∣ Nat.gcd x y :=
        Nat.dvd_gcd (Nat.dvd_of_mod_eq_zero (by omega)) (Nat.dvd_of_mod_eq_zero (by omega))
      rw [hco.gcd_eq_one] at h2
      exact absurd (Nat.le_of_dvd (by decide) h2) (by decide)
    rcases hodd with hx1 | hy1
    · obtain ⟨a, b, c, ha, hb, hc, habc⟩ := descent_step hx1 hy hco h
      exact ih c hc a b ha hb habc
    · obtain ⟨a, b, c, ha, hb, hc, habc⟩ := descent_step hy1 hx hco.symm (by rw [← h]; grind)
      exact ih c hc a b ha hb habc
  · obtain ⟨x', y', hco', hx', hy'⟩ := Nat.exists_coprime x y
    have hd1 : Nat.gcd x y ≠ 1 := hco
    have hd0 : Nat.gcd x y ≠ 0 := by
      intro h0; rw [h0] at hx'; exact hx (by simpa using hx')
    generalize Nat.gcd x y = d at hx' hy' hd1 hd0
    subst hx' hy'
    have hdz : d^2 ∣ z :=
      (Nat.pow_dvd_pow_iff (n := 2) (by decide)).1 ⟨x'^4 + y'^4, by rw [← h]; grind⟩
    obtain ⟨z', rfl⟩ := hdz
    have hd4 : 0 < d^4 := Nat.pow_pos (Nat.pos_of_ne_zero hd0)
    have hx'0 : x' ≠ 0 := by intro h0; subst h0; simp at hx
    have hy'0 : y' ≠ 0 := by intro h0; subst h0; simp at hy
    have heq : x'^4 + y'^4 = z'^2 := by
      apply Nat.eq_of_mul_eq_mul_left hd4
      grind
    have hz'0 : 0 < z' := by
      rcases Nat.eq_zero_or_pos z' with h0 | h0
      · subst h0
        have : 0 < x'^4 := Nat.pow_pos (Nat.pos_of_ne_zero hx'0)
        simp at heq; omega
      · exact h0
    have hlt : z' < d^2 * z' := by
      have : 4 ≤ d^2 := by
        rw [Nat.pow_two]; exact Nat.mul_le_mul (n₁ := 2) (m₁ := 2) (by omega) (by omega)
      have := Nat.mul_le_mul_right z' this
      omega
    exact ih z' hlt x' y' hx'0 hy'0 heq

/-- Fermat's Last Theorem for exponent 4. -/
theorem fermat4 (x y z : Nat) (hx : x ≠ 0) (hy : y ≠ 0) : x^4 + y^4 ≠ z^4 := by
  intro h
  exact no_x4_y4_z2 x y (z^2) hx hy (by rw [h, ← Nat.pow_mul])



/-- Mult-add form of FLT(4): for rationals q ≠ 0 and s, q⁴ is never the mult-add 4th power
    s^(⊞4) = (1 + s)⁴ − 1. -/
theorem fourth_power_ne_mpow_four (q s : Rat) (hq : q ≠ 0) : q^4 ≠ mpow s 4 := by
  intro h
  have hlift : 1 + q^4 = (1 + s)^4 := by rw [h]; exact mpow_lift s 4
  generalize 1 + s = r at hlift
  have den_ne (t : Rat) : (t.den : Rat) ≠ 0 := by
    intro h0
    have h1 : ((t.den : Int) : Rat) = ((0 : Int) : Rat) := by
      rw [Rat.intCast_natCast]; exact h0
    have := Rat.intCast_inj.1 h1
    have := t.den_nz
    omega
  have eq_q : q = (q.num : Rat) / (q.den : Rat) := by rw [← Rat.mkRat_eq_div, Rat.mkRat_self]
  have eq_r : r = (r.num : Rat) / (r.den : Rat) := by rw [← Rat.mkRat_eq_div, Rat.mkRat_self]
  have hB := den_ne q
  have hD := den_ne r
  -- clear denominators in the rationals
  have hR : ((q.den : Rat) * r.den)^4 + ((q.num : Rat) * r.den)^4 = ((r.num : Rat) * q.den)^4 := by
    rw [eq_q, eq_r] at hlift
    grind
  have hA0 : q.num ≠ 0 := fun h0 => hq (Rat.num_eq_zero.1 h0)
  generalize q.num = A at hR hA0
  generalize q.den = B at hR hB
  generalize r.num = C at hR
  generalize r.den = D at hR hD
  -- move to the integers
  have hI : ((B : Int) * D)^4 + (A * D)^4 = (C * B)^4 := by
    apply Rat.intCast_inj.1
    simp only [Rat.intCast_add, Rat.intCast_mul, Rat.intCast_pow, Rat.intCast_natCast]
    exact hR
  -- and to the naturals, through absolute values
  have h4 (t : Int) : (t.natAbs : Int)^4 = t^4 := by
    have := Int.natAbs_pow_two t
    calc (t.natAbs : Int)^4 = ((t.natAbs : Int)^2)^2 := by grind
      _ = (t^2)^2 := by rw [this]
      _ = t^4 := by grind
  have hN : (B * D)^4 + (A.natAbs * D)^4 = (C.natAbs * B)^4 := by
    apply Int.ofNat_inj.1
    simp only [Int.natCast_add, Int.natCast_mul, Int.natCast_pow]
    have hA := h4 A
    have hC := h4 C
    grind
  have hB0 : B ≠ 0 := by intro h0; subst h0; exact hB (by simp)
  have hD0 : D ≠ 0 := by intro h0; subst h0; exact hD (by simp)
  have hA' : A.natAbs ≠ 0 := fun h0 => hA0 (Int.natAbs_eq_zero.1 h0)
  exact fermat4 (B * D) (A.natAbs * D) (C.natAbs * B)
    (Nat.mul_ne_zero hB0 hD0) (Nat.mul_ne_zero hA' hD0) hN


