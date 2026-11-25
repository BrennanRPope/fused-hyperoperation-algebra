module

public import FusedAlgebra.Topics.Elliptic

/-!
# The elementary pieces of the Wiles route to Fermat's Last Theorem

Fermat's Last Theorem for the exponent `n` is mathlib's `FermatLastTheoremFor n`: no nonzero
`a b c : ℕ` with `a ^ n + b ^ n = c ^ n`. This file collects the elementary steps around the Frey
curve. **The deep steps of the route are not proved here**: modularity of semistable elliptic
curves (Wiles, Taylor–Wiles), irreducibility of the mod `p` representation (Mazur) and level
lowering (Ribet) appear nowhere in this file, and neither does the genus formula for `X₀(N)`.

This is the mathlib form of `FusedAlgebra.Topics.Fermat.FreyPieces`. The original built its own
`FLT`, primality, Weierstrass invariants and congruences; here they are mathlib's
`FermatLastTheoremFor`, `Nat.Prime`, `WeierstrassCurve` (with its `c₄` and `Δ`) and `Int.ModEq`.

## Main results

1. Exponent reduction (known; in mathlib as `FermatLastTheoremFor.mono` and
   `FermatLastTheorem.of_odd_primes`): `flt_reduction`, `flt_reduction_frey`.
2. Normalisation: `coprime_reduction`, `flt_of_coprime`, `coprime_parity`, and
   `frey_normalize`: a coprime integer solution for odd `p ≥ 5` can be re-signed and reordered so
   that `a ≡ 3 mod 4`, `b` is even and `32 ∣ b ^ p`.
3. The Frey curve `y² = x(x - A)(x + B)` with `A + B = C` (`EllipticLead.freyCurve`):
   its 2-torsion (`frey_two_torsion`), the model `Y² + XY = X³ + a₂X² + a₄X` (`freyModel`,
   `frey_minimal_model`, `variableChange_freyCurve`), the invariants `c₄ = A² + AB + B²` and
   `256Δ = (ABC)²` (`frey_c4`, `frey_disc`), semistability (`frey_semistable`, `frey_at_two`),
   reduction to the node whose chord law is mult-add (`freyModel_map_zmod`,
   `node_chord_multAdd`), and the fact that every odd prime's exact exponent in `Δ` is a multiple
   of `p` (`frey_disc_valuation`).
4. Mult-add restatement (`flt_multAdd`): `x ^ n + y ^ n = z ^ n ↔ (y / x) ^ n = mpow (z / x - 1) n`.
5. Level 2: the counts `μ = 3`, `ν₂ = 1`, `ν₃ = 0`, `ν∞ = 2` for `Γ₀(2)` (`level2_counts`),
   which put the genus of `X₀(2)` at `0` *through the genus formula, which is cited as a
   hypothesis and not proved* (`genus_X0_2`).
-/

@[expose] public section

namespace FreyPieces

open MultAdd EllipticLead

/-! ### 1. Exponent reduction -/

/-- Mathlib's `FermatLastTheoremFor n` asks for `c ≠ 0` as well; the original library's `FLT n`
did not. The two say the same thing. -/
theorem fermatLastTheoremFor_iff {n : ℕ} :
    FermatLastTheoremFor n ↔ ∀ x y z : ℕ, x ≠ 0 → y ≠ 0 → x ^ n + y ^ n ≠ z ^ n := by
  refine ⟨fun h x y z hx hy hxyz ↦ h x y z hx hy ?_ hxyz, fun h x y z hx hy _ ↦ h x y z hx hy⟩
  rintro rfl
  rcases eq_or_ne n 0 with rfl | hn
  · simp at hxyz
  · rw [zero_pow hn] at hxyz
    exact pow_ne_zero n hx (Nat.eq_zero_of_add_eq_zero_right hxyz)

/-- Fermat's Last Theorem for `4` and for every odd prime gives it for every `n ≥ 3`. This is
known; mathlib's `FermatLastTheorem.of_odd_primes` is the same reduction with exponent `4`
already discharged. -/
theorem flt_reduction (h4 : FermatLastTheoremFor 4)
    (hp : ∀ p : ℕ, p.Prime → p ≠ 2 → FermatLastTheoremFor p) : FermatLastTheorem := by
  intro n hn
  obtain h | ⟨p, hpp, hpn, hodd⟩ := Nat.four_dvd_or_exists_odd_prime_and_dvd_of_two_lt hn
  · exact h4.mono h
  · exact (hp p hpp fun h ↦ by subst h; exact absurd hodd (by decide)).mono hpn

/-- Fermat's Last Theorem for `4`, for `3` and for every prime `p ≥ 5` (the Frey range) gives it
for every `n ≥ 3`. -/
theorem flt_reduction_frey (h4 : FermatLastTheoremFor 4) (h3 : FermatLastTheoremFor 3)
    (hp : ∀ p : ℕ, p.Prime → 5 ≤ p → FermatLastTheoremFor p) : FermatLastTheorem := by
  refine flt_reduction h4 fun p hpp hp2 ↦ ?_
  obtain rfl | h : p = 3 ∨ 5 ≤ p := by
    have := hpp.two_le
    have : p ≠ 4 := by rintro rfl; norm_num at hpp
    lia
  · exact h3
  · exact hp p hpp h

/-! ### 2. Normalisation -/

/-- A solution with `x, y ≠ 0` shrinks (divide by `gcd x y`) to a coprime one. -/
theorem coprime_reduction {n x y z : ℕ} (hn : n ≠ 0) (hx : x ≠ 0) (hy : y ≠ 0)
    (h : x ^ n + y ^ n = z ^ n) :
    ∃ x' y' z', x' ≠ 0 ∧ y' ≠ 0 ∧ Nat.gcd x' y' = 1 ∧ x' ^ n + y' ^ n = z' ^ n ∧
      x' ∣ x ∧ y' ∣ y := by
  obtain ⟨x', y', hcop, hx', hy'⟩ := Nat.exists_coprime x y
  have hg : x.gcd y ≠ 0 := Nat.gcd_ne_zero_left hx
  generalize x.gcd y = g at hx' hy' hg
  have hsum : z ^ n = (x' ^ n + y' ^ n) * g ^ n := by rw [← h, hx', hy']; ring
  obtain ⟨z', rfl⟩ : g ∣ z := (Nat.pow_dvd_pow_iff hn).1 (hsum ▸ Dvd.intro_left _ rfl)
  refine ⟨x', y', z', ?_, ?_, hcop, ?_, Dvd.intro _ hx'.symm, Dvd.intro _ hy'.symm⟩
  · rintro rfl; exact hx (by simpa using hx')
  · rintro rfl; exact hy (by simpa using hy')
  · rw [mul_pow, mul_comm] at hsum
    exact (Nat.eq_of_mul_eq_mul_right (pow_pos (Nat.pos_of_ne_zero hg) n) hsum).symm

/-- So Fermat's Last Theorem for `n` holds once coprime solutions are ruled out. -/
theorem flt_of_coprime {n : ℕ} (hn : n ≠ 0)
    (h : ∀ x y z : ℕ, x ≠ 0 → y ≠ 0 → Nat.gcd x y = 1 → x ^ n + y ^ n ≠ z ^ n) :
    FermatLastTheoremFor n := by
  refine fermatLastTheoremFor_iff.2 fun x y z hx hy hxyz ↦ ?_
  obtain ⟨x', y', z', hx', hy', hc, he, -, -⟩ := coprime_reduction hn hx hy hxyz
  exact h x' y' z' hx' hy' hc he

/-- In a coprime solution exactly one of `x`, `y`, `z` is even. -/
theorem coprime_parity {n x y z : ℕ} (hn : n ≠ 0) (hc : Nat.gcd x y = 1)
    (h : x ^ n + y ^ n = z ^ n) :
    (x % 2 = 0 ∧ y % 2 = 1 ∧ z % 2 = 1) ∨ (x % 2 = 1 ∧ y % 2 = 0 ∧ z % 2 = 1) ∨
      (x % 2 = 1 ∧ y % 2 = 1 ∧ z % 2 = 0) := by
  have hpow : ∀ w : ℕ, w ^ n % 2 = w % 2 := fun w ↦ by
    rw [Nat.pow_mod]
    rcases Nat.mod_two_eq_zero_or_one w with hw | hw <;> simp [hw, hn]
  have hm : (x ^ n + y ^ n) % 2 = z ^ n % 2 := by rw [h]
  rw [Nat.add_mod, hpow, hpow, hpow] at hm
  have hnot : ¬(x % 2 = 0 ∧ y % 2 = 0) := fun ⟨h1, h2⟩ ↦ by
    have := Nat.dvd_gcd (Nat.dvd_of_mod_eq_zero h1) (Nat.dvd_of_mod_eq_zero h2)
    rw [hc] at this
    lia
  lia

/-- Over `ℤ`, `x ^ n` has the parity of `x` for `n ≠ 0`. -/
theorem int_pow_emod_two (x : ℤ) {n : ℕ} (hn : n ≠ 0) : x ^ n % 2 = x % 2 := by
  have h : x ^ n ≡ (x % 2) ^ n [ZMOD 2] := (Int.mod_modEq x 2).symm.pow n
  rw [Int.ModEq] at h
  rcases Int.emod_two_eq_zero_or_one x with hx | hx <;> simpa [hx, hn] using h

/-- From `o₁ ^ p + e ^ p + o₂ ^ p = 0` with `o₁` odd and `e` even: re-sign to `a' ≡ 3 mod 4` and
`b'` even, keeping the absolute values. -/
theorem normalize_core {p : ℕ} (hp : Odd p) {o₁ e o₂ : ℤ} (h1 : o₁ % 2 = 1) (he : e % 2 = 0)
    (h : o₁ ^ p + e ^ p + o₂ ^ p = 0) :
    ∃ a' b' c' : ℤ, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0 ∧
      a'.natAbs = o₁.natAbs ∧ b'.natAbs = e.natAbs ∧ c'.natAbs = o₂.natAbs := by
  obtain h4 | h4 : o₁ % 4 = 3 ∨ o₁ % 4 = 1 := by lia
  · refine ⟨o₁, e, -o₂, ?_, h4, he, rfl, rfl, Int.natAbs_neg o₂⟩
    rw [hp.neg_pow]; linear_combination h
  · refine ⟨-o₁, -e, o₂, ?_, by lia, by lia, Int.natAbs_neg o₁, Int.natAbs_neg e, rfl⟩
    rw [hp.neg_pow, hp.neg_pow]; linear_combination -h

/-- `b` even and `p ≥ 5` give `32 ∣ b ^ p`. -/
theorem two_pow_dvd {p : ℕ} (h5 : 5 ≤ p) {b : ℤ} (hb : b % 2 = 0) : (32 : ℤ) ∣ b ^ p :=
  calc (32 : ℤ) = 2 ^ 5 := by norm_num
    _ ∣ b ^ 5 := pow_dvd_pow_of_dvd (Int.dvd_of_emod_eq_zero hb) 5
    _ ∣ b ^ p := pow_dvd_pow b h5

/-- In a solution `a ^ p + b ^ p = c ^ p` with `gcd a b = 1`, also `gcd a c = 1`. -/
theorem gcd_eq_one_of_add_pow_eq {p : ℕ} {a b c : ℤ} (hp : p ≠ 0) (hcop : Int.gcd a b = 1)
    (h : a ^ p + b ^ p = c ^ p) : Int.gcd a c = 1 := by
  rw [← Int.isCoprime_iff_gcd_eq_one] at hcop ⊢
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hp
  have h' := ((IsCoprime.pow_right_iff m.succ_pos).2 hcop).add_mul_left_right (a ^ m)
  rwa [← pow_succ', add_comm, h, IsCoprime.pow_right_iff m.succ_pos] at h'

/-- Frey normalisation: a coprime integer solution for odd `p ≥ 5` becomes one with
`a' ≡ 3 mod 4`, `b'` even and `32 ∣ b' ^ p`, still coprime and with the same `|abc|`. -/
theorem frey_normalize {p : ℕ} (hp : Odd p) (h5 : 5 ≤ p) {a b c : ℤ} (hcop : Int.gcd a b = 1)
    (habc : a * b * c ≠ 0) (h : a ^ p + b ^ p = c ^ p) :
    ∃ a' b' c' : ℤ, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0 ∧
      (a' * b' * c').natAbs = (a * b * c).natAbs ∧ Int.gcd a' b' = 1 ∧ a' * b' * c' ≠ 0 ∧
      (32 : ℤ) ∣ b' ^ p := by
  have hp0 : p ≠ 0 := by lia
  have hac : Int.gcd a c = 1 := gcd_eq_one_of_add_pow_eq hp0 hcop h
  have hbc : Int.gcd b c = 1 :=
    gcd_eq_one_of_add_pow_eq hp0 (Int.gcd_comm a b ▸ hcop) (add_comm (a ^ p) _ ▸ h)
  have ha := int_pow_emod_two a hp0
  have hb := int_pow_emod_two b hp0
  have hc := int_pow_emod_two c hp0
  have hnot : ¬(a % 2 = 0 ∧ b % 2 = 0) := fun ⟨h1, h2⟩ ↦ by
    have := Int.dvd_coe_gcd (Int.dvd_of_emod_eq_zero h1) (Int.dvd_of_emod_eq_zero h2)
    rw [hcop] at this
    norm_num at this
  have hm : (a ^ p + b ^ p) % 2 = c ^ p % 2 := by rw [h]
  rw [Int.add_emod, ha, hb, hc] at hm
  have hnc := hp.neg_pow c
  suffices H : ∃ a' b' c' : ℤ, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0 ∧
      (a' * b' * c').natAbs = (a * b * c).natAbs ∧ Int.gcd a' b' = 1 by
    obtain ⟨a', b', c', h1, h2, h3, h4, hg⟩ := H
    refine ⟨a', b', c', h1, h2, h3, h4, hg, fun h0 ↦ habc ?_, two_pow_dvd h5 h3⟩
    rw [h0] at h4
    exact Int.natAbs_eq_zero.1 h4.symm
  obtain ⟨h1, h2, -⟩ | ⟨h1, h2⟩ | ⟨h1, h2, h3⟩ :
      (a % 2 = 0 ∧ b % 2 = 1 ∧ c % 2 = 1) ∨ (a % 2 = 1 ∧ b % 2 = 0) ∨
        (a % 2 = 1 ∧ b % 2 = 1 ∧ c % 2 = 0) := by lia
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core hp h2 h1 (o₂ := -c) (by rw [hnc]; linear_combination h)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · simp only [Int.natAbs_mul, e4, e5, e6, Int.natAbs_neg]; ring
    · rw [Int.gcd_def, e4, e5, ← Int.gcd_def, Int.gcd_comm]; exact hcop
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core hp h1 h2 (o₂ := -c) (by rw [hnc]; linear_combination h)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · simp only [Int.natAbs_mul, e4, e5, e6, Int.natAbs_neg]
    · rw [Int.gcd_def, e4, e5, ← Int.gcd_def]; exact hcop
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core hp h1 (show -c % 2 = 0 by lia) (o₂ := b) (by rw [hnc]; linear_combination h)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · simp only [Int.natAbs_mul, e4, e5, e6, Int.natAbs_neg]; ring
    · rw [Int.gcd_def, e4, e5, Int.natAbs_neg, ← Int.gcd_def]; exact hac

/-! ### 3. The Frey curve `y² = x(x - A)(x + B)` -/

section CommRing

variable {R : Type*} [CommRing R]

/-- The 2-torsion points `(0, 0)`, `(A, 0)`, `(-B, 0)` lie on the Frey curve. -/
theorem frey_two_torsion (A B : R) :
    (freyCurve A B).toAffine.Equation 0 0 ∧ (freyCurve A B).toAffine.Equation A 0 ∧
      (freyCurve A B).toAffine.Equation (-B) 0 := by
  simp only [freyCurve_equation]
  refine ⟨?_, ?_, ?_⟩ <;> ring

/-- `16` times the squared product of the root differences of `x(x - A)(x + B)` is `16(ABC)²`
when `A + B = C`. -/
theorem frey_disc_roots (A B C : R) (hC : A + B = C) :
    16 * ((0 - A) * (A - -B) * (0 - -B)) ^ 2 = 16 * (A * B * C) ^ 2 := by
  subst hC; ring

/-- The identity behind the change of variables `x = 4X`, `y = 8Y + 4X`. -/
theorem frey_change_of_vars (A B X Y : R) :
    (8 * Y + 4 * X) ^ 2 - 4 * X * (4 * X - A) * (4 * X + B) =
      64 * (Y ^ 2 + X * Y - X ^ 3) - 16 * (B - A - 1) * X ^ 2 + 4 * (A * B) * X := by
  ring

/-- The model `Y² + XY = X³ + a₂X² + a₄X` (`a₁ = 1`, `a₃ = a₆ = 0`). -/
def freyModel (a₂ a₄ : R) : WeierstrassCurve R := ⟨1, a₂, 0, a₄, 0⟩

/-- `c₄` of the Frey model is `A² + AB + B²`. -/
theorem frey_c4 {A B a₂ a₄ : R} (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) :
    (freyModel a₂ a₄).c₄ = A ^ 2 + A * B + B ^ 2 := by
  simp only [freyModel, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  linear_combination (4 * a₂ + B - A + 1) * h2 - 3 * h4

/-- `256Δ` of the Frey model is `(ABC)²`. -/
theorem frey_disc {A B C a₂ a₄ : R} (hC : A + B = C) (h2 : 4 * a₂ = B - A - 1)
    (h4 : 16 * a₄ = -(A * B)) : 256 * (freyModel a₂ a₄).Δ = (A * B * C) ^ 2 := by
  subst hC
  have hb : 1 + 4 * a₂ = B - A := by linear_combination h2
  calc 256 * (freyModel a₂ a₄).Δ = (16 * a₄) ^ 2 * ((1 + 4 * a₂) ^ 2 - 4 * (16 * a₄)) := by
        simp only [freyModel, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
          WeierstrassCurve.b₆, WeierstrassCurve.b₈]
        ring
    _ = (A * B * (A + B)) ^ 2 := by rw [hb, h4]; ring

/-- The point `(4X, 8Y + 4X)` is on the Frey curve exactly when `(X, Y)` is on the model, as
soon as `64` can be cancelled. -/
theorem frey_model_equation [IsDomain R] (h64 : (64 : R) ≠ 0) {A B a₂ a₄ : R}
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) (X Y : R) :
    (freyCurve A B).toAffine.Equation (4 * X) (8 * Y + 4 * X) ↔
      (freyModel a₂ a₄).toAffine.Equation X Y := by
  rw [freyCurve_equation, WeierstrassCurve.Affine.equation_iff]
  simp only [freyModel]
  have key : (8 * Y + 4 * X) ^ 2 - 4 * X * (4 * X - A) * (4 * X + B) =
      64 * (Y ^ 2 + 1 * X * Y + 0 * Y - (X ^ 3 + a₂ * X ^ 2 + a₄ * X + 0)) := by
    linear_combination (16 * X ^ 2) * h2 + 4 * X * h4
  rw [← sub_eq_zero, key, mul_eq_zero, or_iff_right h64, sub_eq_zero]

end CommRing

/-- Over a field in which `2 ≠ 0`, the Frey model is the Frey curve after mathlib's change of
variables `x = 4X`, `y = 8Y + 4X` (that is `u = 2`, `r = 0`, `s = 1`, `t = 0`). -/
theorem variableChange_freyCurve {K : Type*} [Field K] (h2 : (2 : K) ≠ 0) (A B : K) :
    (⟨Units.mk0 2 h2, 0, 1, 0⟩ : WeierstrassCurve.VariableChange K) • freyCurve A B =
      freyModel ((B - A - 1) / 2 ^ 2) (-(A * B) / 2 ^ 4) := by
  ext <;> simp [freyCurve, freyModel, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₃,
    WeierstrassCurve.variableChange_a₄, WeierstrassCurve.variableChange_a₆] <;> field_simp

/-- `A ≡ 3 mod 4` and `16 ∣ B` give `4 ∣ B - A - 1` and `16 ∣ A * B`. -/
theorem frey_minimal_divisibility {A B : ℤ} (hA : A % 4 = 3) (hB : 16 ∣ B) :
    4 ∣ B - A - 1 ∧ 16 ∣ A * B :=
  ⟨by lia, hB.mul_left A⟩

/-- The coefficient `a₂ = (B - A - 1) / 4` of the minimal model. -/
def a2of (A B : ℤ) : ℤ := (B - A - 1) / 4

/-- The coefficient `a₄ = -AB / 16` of the minimal model. -/
def a4of (A B : ℤ) : ℤ := -(A * B / 16)

/-- With these integers, `(4X, 8Y + 4X)` is on the Frey curve exactly when
`Y² + XY = X³ + a₂X² + a₄X`. -/
theorem frey_minimal_model {A B : ℤ} (hA : A % 4 = 3) (hB : 16 ∣ B) :
    4 * a2of A B = B - A - 1 ∧ 16 * a4of A B = -(A * B) ∧
      ∀ X Y : ℤ, (freyCurve A B).toAffine.Equation (4 * X) (8 * Y + 4 * X) ↔
        (freyModel (a2of A B) (a4of A B)).toAffine.Equation X Y := by
  obtain ⟨h4, h16⟩ := frey_minimal_divisibility hA hB
  have h1 : 4 * a2of A B = B - A - 1 := Int.mul_ediv_cancel' h4
  have h2 : 16 * a4of A B = -(A * B) := by rw [a4of, mul_neg, Int.mul_ediv_cancel' h16]
  exact ⟨h1, h2, frey_model_equation (by norm_num) h1 h2⟩

/-! ### Semistability: only nodes -/

/-- A prime cannot divide two coprime integers. -/
theorem prime_dvd_gcd_false {q : ℕ} (hq : q.Prime) {A B : ℤ} (hcop : Int.gcd A B = 1)
    (hA : (q : ℤ) ∣ A) (hB : (q : ℤ) ∣ B) : False := by
  have h := Int.dvd_coe_gcd hA hB
  rw [hcop] at h
  exact hq.one_lt.ne' (Nat.dvd_one.1 (Int.natCast_dvd_natCast.1 h))

/-- If `gcd A B = 1`, a prime dividing `ABC` does not divide `A² + AB + B²`. -/
theorem frey_c4_not_dvd {A B C : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1) {q : ℕ}
    (hq : q.Prime) (hd : (q : ℤ) ∣ A * B * C) : ¬(q : ℤ) ∣ A ^ 2 + A * B + B ^ 2 := by
  have hq' : Prime (q : ℤ) := Nat.prime_iff_prime_int.1 hq
  subst hC
  intro hc
  obtain (hA | hB) | hAB := (hq'.dvd_or_dvd hd).imp_left hq'.dvd_or_dvd
  · refine prime_dvd_gcd_false hq hcop hA (hq'.dvd_of_dvd_pow (n := 2) ?_)
    have e : B ^ 2 = A ^ 2 + A * B + B ^ 2 - A * (A + B) := by ring
    rw [e]; exact dvd_sub hc (hA.mul_right _)
  · refine prime_dvd_gcd_false hq hcop (hq'.dvd_of_dvd_pow (n := 2) ?_) hB
    have e : A ^ 2 = A ^ 2 + A * B + B ^ 2 - B * (A + B) := by ring
    rw [e]; exact dvd_sub hc (hB.mul_right _)
  · have hA : (q : ℤ) ∣ A := by
      refine hq'.dvd_of_dvd_pow (n := 2) ?_
      have e : A ^ 2 = A ^ 2 + A * B + B ^ 2 - B * (A + B) := by ring
      rw [e]; exact dvd_sub hc (hAB.mul_left _)
    exact prime_dvd_gcd_false hq hcop hA ((dvd_add_right hA).1 hAB)

/-- Semistability: a prime dividing `Δ` never divides `c₄` (a node, not a cusp). -/
theorem frey_semistable {A B C a₂ a₄ : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1)
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) {q : ℕ} (hq : q.Prime)
    (hd : (q : ℤ) ∣ (freyModel a₂ a₄).Δ) : ¬(q : ℤ) ∣ (freyModel a₂ a₄).c₄ := by
  rw [frey_c4 h2 h4]
  refine frey_c4_not_dvd hC hcop hq ((Nat.prime_iff_prime_int.1 hq).dvd_of_dvd_pow (n := 2) ?_)
  rw [← frey_disc hC h2 h4]
  exact hd.mul_left 256

/-- At `2`: `A` odd and `32 ∣ B` make `c₄` odd and `4 ∣ Δ`, so `2` is a bad prime with a node. -/
theorem frey_at_two {A B C a₂ a₄ : ℤ} (hC : A + B = C) (hA : A % 2 = 1) (hB : 32 ∣ B)
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) :
    (freyModel a₂ a₄).c₄ % 2 = 1 ∧ 4 ∣ (freyModel a₂ a₄).Δ := by
  obtain ⟨m, rfl⟩ := hB
  constructor
  · rw [frey_c4 h2 h4]
    obtain ⟨k, rfl⟩ : ∃ k, A = 2 * k + 1 := ⟨A / 2, by lia⟩
    have e : (2 * k + 1) ^ 2 + (2 * k + 1) * (32 * m) + (32 * m) ^ 2 =
        2 * (2 * k ^ 2 + 2 * k + 32 * k * m + 16 * m + 512 * m ^ 2) + 1 := by ring
    rw [e]; lia
  · have hD := frey_disc hC h2 h4
    have e : (A * (32 * m) * C) ^ 2 = 1024 * (A * m * C) ^ 2 := by ring
    rw [e] at hD
    exact ⟨(A * m * C) ^ 2, by lia⟩

/-! ### Node shape and the mult-add chord law -/

/-- For an odd prime `q` dividing `A`, `q` divides `a₄`. -/
theorem frey_a4_mod_q {A B a₄ : ℤ} (h4 : 16 * a₄ = -(A * B)) {q : ℕ} (hq : q.Prime)
    (h2 : q ≠ 2) (hA : (q : ℤ) ∣ A) : (q : ℤ) ∣ a₄ := by
  have hq' : Prime (q : ℤ) := Nat.prime_iff_prime_int.1 hq
  have h16 : (q : ℤ) ∣ 2 ^ 4 * a₄ := by
    rw [show (2 : ℤ) ^ 4 = 16 by norm_num, h4]; exact (hA.mul_right B).neg_right
  refine (hq'.dvd_or_dvd h16).resolve_left fun h ↦ h2 ?_
  have h' : q ∣ 2 := by exact_mod_cast hq'.dvd_of_dvd_pow h
  exact (Nat.prime_dvd_prime_iff_eq hq Nat.prime_two).1 h'

/-- Then modulo `q` the model is the node family `Y² + XY = X³ + a₂X²`, equation by equation. -/
theorem frey_reduces_to_node {a₂ a₄ : ℤ} {q : ℕ} (h : (q : ℤ) ∣ a₄) (X Y : ℤ) :
    Y ^ 2 + X * Y - (X ^ 3 + a₂ * X ^ 2 + a₄ * X) ≡
      Y ^ 2 + 1 * X * Y - (X ^ 3 + a₂ * X ^ 2) [ZMOD q] := by
  rw [Int.modEq_iff_dvd]
  have e : Y ^ 2 + 1 * X * Y - (X ^ 3 + a₂ * X ^ 2) -
      (Y ^ 2 + X * Y - (X ^ 3 + a₂ * X ^ 2 + a₄ * X)) = a₄ * X := by ring
  rw [e]; exact h.mul_right X

/-- The same as an equality of curves: the reduction of the Frey model modulo `q` is the nodal
cubic with `a₁ = 1`. -/
theorem freyModel_map_zmod {a₂ a₄ : ℤ} {q : ℕ} (h : (q : ℤ) ∣ a₄) :
    (freyModel a₂ a₄).map (Int.castRingHom (ZMod q)) = nodalCurve 1 (a₂ : ZMod q) := by
  ext <;> simp [freyModel, nodalCurve, (ZMod.intCast_zmod_eq_zero_iff_dvd a₄ q).2 h]

section Chord

variable {R : Type*} [CommRing R]

/-- The `x`-coordinate of the second meeting point of the line `y = tx` with the node of slopes
`m₁`, `m₂`. -/
def nodeX (m₁ m₂ t : R) : R := (t - m₁) * (t - m₂)

/-- The collinearity determinant of three points. -/
def chordDet (x₁ y₁ x₂ y₂ x₃ y₃ : R) : R :=
  x₁ * (y₂ - y₃) - y₁ * (x₂ - x₃) + (x₂ * y₃ - x₃ * y₂)

/-- `(m₂ - m₁)` times the collinearity determinant of three points of the node is the
Vandermonde of the slopes times `∏ (tᵢ - m₁) - ∏ (tᵢ - m₂)`. -/
theorem node_chord_identity (m₁ m₂ t₁ t₂ t₃ : R) :
    (m₂ - m₁) * chordDet (nodeX m₁ m₂ t₁) (t₁ * nodeX m₁ m₂ t₁) (nodeX m₁ m₂ t₂)
        (t₂ * nodeX m₁ m₂ t₂) (nodeX m₁ m₂ t₃) (t₃ * nodeX m₁ m₂ t₃) =
      (t₁ - t₂) * (t₂ - t₃) * (t₃ - t₁) *
        ((t₁ - m₁) * (t₂ - m₁) * (t₃ - m₁) - (t₁ - m₂) * (t₂ - m₂) * (t₃ - m₂)) := by
  unfold chordDet nodeX; ring

/-- The point `(x, tx)` with `x = nodeX m₁ m₂ t` lies on the nodal cubic with `a₁ = -(m₁ + m₂)`
and `a₂ = -(m₁m₂)`. -/
theorem node_point (m₁ m₂ t : R) :
    (nodalCurve (-(m₁ + m₂)) (-(m₁ * m₂))).toAffine.Equation (nodeX m₁ m₂ t)
      (t * nodeX m₁ m₂ t) := by
  have h := node_param (-(m₁ + m₂)) (-(m₁ * m₂)) t
  rwa [slopes_factor] at h

end Chord

/-- Three points of the node are collinear exactly when `s₁ ✳ s₂ ✳ s₃ = 0`, where
`sᵢ = (m₂ - m₁) / (tᵢ - m₂)`: the chord law of the node is mult-add. -/
theorem node_chord_multAdd {K : Type*} [Field K] (m₁ m₂ t₁ t₂ t₃ : K) (hm : m₁ ≠ m₂)
    (h12 : t₁ ≠ t₂) (h23 : t₂ ≠ t₃) (h31 : t₃ ≠ t₁) (h1 : t₁ ≠ m₂) (h2 : t₂ ≠ m₂)
    (h3 : t₃ ≠ m₂) :
    chordDet (nodeX m₁ m₂ t₁) (t₁ * nodeX m₁ m₂ t₁) (nodeX m₁ m₂ t₂)
        (t₂ * nodeX m₁ m₂ t₂) (nodeX m₁ m₂ t₃) (t₃ * nodeX m₁ m₂ t₃) = 0 ↔
      (m₂ - m₁) / (t₁ - m₂) ✳ (m₂ - m₁) / (t₂ - m₂) ✳ (m₂ - m₁) / (t₃ - m₂) = 0 := by
  have e1 := sub_ne_zero.2 h1
  have e2 := sub_ne_zero.2 h2
  have e3 := sub_ne_zero.2 h3
  have hV : (t₁ - t₂) * (t₂ - t₃) * (t₃ - t₁) ≠ 0 :=
    mul_ne_zero (mul_ne_zero (sub_ne_zero.2 h12) (sub_ne_zero.2 h23)) (sub_ne_zero.2 h31)
  have hB : (t₁ - m₂) * (t₂ - m₂) * (t₃ - m₂) ≠ 0 := mul_ne_zero (mul_ne_zero e1 e2) e3
  have hd : m₂ - m₁ ≠ 0 := sub_ne_zero.2 hm.symm
  have hlift : (1 + ((m₂ - m₁) / (t₁ - m₂) ✳ (m₂ - m₁) / (t₂ - m₂) ✳ (m₂ - m₁) / (t₃ - m₂))) *
      ((t₁ - m₂) * (t₂ - m₂) * (t₃ - m₂)) = (t₁ - m₁) * (t₂ - m₁) * (t₃ - m₁) := by
    rw [one_add_multAdd, one_add_multAdd]; field_simp; ring
  have key := node_chord_identity m₁ m₂ t₁ t₂ t₃
  rw [← hlift, add_mul, one_mul, add_sub_cancel_left, ← mul_assoc] at key
  constructor <;> intro h
  · rw [h, mul_zero] at key
    exact (mul_eq_zero.1 ((mul_eq_zero.1 key.symm).resolve_right hB)).resolve_left hV
  · rw [h, mul_zero, zero_mul] at key
    exact (mul_eq_zero.1 key).resolve_left hd

/-! ### The Tate parameter is a `p`-th power -/

/-- With `A = a ^ p`, `B = b ^ p`, `C = c ^ p`: `256Δ = ((abc)²) ^ p`. -/
theorem frey_disc_pth_power {R : Type*} [CommRing R] {a b c Δ : R} {p : ℕ}
    (hΔ : 256 * Δ = (a ^ p * b ^ p * c ^ p) ^ 2) : 256 * Δ = ((a * b * c) ^ 2) ^ p := by
  rw [hΔ]; ring

/-- The exact exponent of a prime in a `p`-th power is a multiple of `p`. -/
theorem exact_exponent_pth_power {q M e : ℕ} (p : ℕ) (hq : q.Prime) (hM : M ≠ 0)
    (he : q ^ e ∣ M ^ p) (he1 : ¬q ^ (e + 1) ∣ M ^ p) : p ∣ e := by
  rw [hq.pow_dvd_iff_le_factorization (pow_ne_zero p hM), Nat.factorization_pow] at he he1
  exact ⟨M.factorization q, by simp only [Finsupp.smul_apply, smul_eq_mul] at he he1; lia⟩

/-- For every odd prime `q`, the exact exponent of `q` in the discriminant `Δ` of the Frey model
of a solution `a ^ p + b ^ p = c ^ p` is a multiple of `p`. -/
theorem frey_disc_valuation {a b c Δ : ℤ} {p : ℕ} (hΔ : 256 * Δ = (a ^ p * b ^ p * c ^ p) ^ 2)
    (habc : a * b * c ≠ 0) {q : ℕ} (hq : q.Prime) (h2 : q ≠ 2) {e : ℕ} (he : (q : ℤ) ^ e ∣ Δ)
    (he1 : ¬(q : ℤ) ^ (e + 1) ∣ Δ) : p ∣ e := by
  have hN : 2 ^ 8 * Δ.natAbs = ((a * b * c).natAbs ^ 2) ^ p := by
    simpa [Int.natAbs_mul, Int.natAbs_pow] using congrArg Int.natAbs (frey_disc_pth_power hΔ)
  have hcop : ∀ f, Nat.Coprime (q ^ f) (2 ^ 8) := fun f ↦
    Nat.Coprime.pow f 8 ((Nat.coprime_primes hq Nat.prime_two).2 h2)
  have tr : ∀ f, (q : ℤ) ^ f ∣ Δ ↔ q ^ f ∣ ((a * b * c).natAbs ^ 2) ^ p := fun f ↦ by
    rw [← Int.natAbs_dvd_natAbs, Int.natAbs_pow, Int.natAbs_natCast, ← hN]
    exact ⟨fun h ↦ h.mul_left _, (hcop f).dvd_of_dvd_mul_left⟩
  exact exact_exponent_pth_power p hq (pow_ne_zero 2 (Int.natAbs_ne_zero.2 habc)) ((tr e).1 he)
    fun h ↦ he1 ((tr (e + 1)).2 h)

/-! ### 4. Mult-add restatement -/

/-- `x ^ n + y ^ n = z ^ n` exactly when `(y / x) ^ n = mpow (z / x - 1) n`: an `n`-th power equals
a mult-add `n`-th power. -/
theorem flt_multAdd {K : Type*} [Field K] (x y z : K) (hx : x ≠ 0) (n : ℕ) :
    x ^ n + y ^ n = z ^ n ↔ (y / x) ^ n = mpow (z / x - 1) n := by
  have hxn : x ^ n ≠ 0 := pow_ne_zero n hx
  rw [mpow_eq, add_sub_cancel, div_pow, div_pow, eq_sub_iff_add_eq', ← div_self hxn, ← add_div,
    div_left_inj₀ hxn]

/-! ### 5. Level 2 target: `X₀(2)` has genus `0` -/

/-- `μ = [SL₂(ℤ) : Γ₀(2)] = #ℙ¹(ℤ/2) = 3`, counted as the nonzero vectors of `(ℤ/2)²` (the only
unit of `ℤ/2` is `1`). -/
def mu2 : ℕ := Finset.card {v : ZMod 2 × ZMod 2 | v ≠ 0}

/-- `ν₂`, the number of `x` modulo `2` with `x² + 1 = 0`. -/
def nu2 : ℕ := Finset.card {x : ZMod 2 | x ^ 2 + 1 = 0}

/-- `ν₃`, the number of `x` modulo `2` with `x² + x + 1 = 0`. -/
def nu3 : ℕ := Finset.card {x : ZMod 2 | x ^ 2 + x + 1 = 0}

/-- `ν∞ = ∑ d ∣ 2, φ(gcd(d, 2 / d))`, the number of cusps of `Γ₀(2)`. -/
def nuInf : ℕ := ∑ d ∈ Nat.divisors 2, Nat.totient (Nat.gcd d (2 / d))

/-- `μ = 3`, `ν₂ = 1`, `ν₃ = 0`, `ν∞ = 2` for `Γ₀(2)`, by computation. -/
theorem level2_counts : mu2 = 3 ∧ nu2 = 1 ∧ nu3 = 0 ∧ nuInf = 2 := by
  decide

/-- The genus formula `12g = 12 + μ - 3ν₂ - 4ν₃ - 6ν∞` gives `g(X₀(2)) = 0`. The genus formula
is a hypothesis here: it is cited, not proved. -/
theorem genus_X0_2 (g : ℤ) (hformula : 12 * g = 12 + (mu2 : ℤ) - 3 * nu2 - 4 * nu3 - 6 * nuInf) :
    g = 0 := by
  obtain ⟨h1, h2, h3, h4⟩ := level2_counts
  rw [h1, h2, h3, h4] at hformula
  lia

end FreyPieces
