module

public import FusedAlgebra.Basic

/-!
# Fermat's Last Theorem for exponent 3, and its mult-add form

The core-Lean library proved by hand, with Euler's infinite descent in the Eisenstein integers,
that `x ^ 3 + y ^ 3 = z ^ 3` has no solution in nonzero integers. That is known mathematics and
is already in mathlib: `fermatLastTheoremThree : FermatLastTheoremFor 3`
(`Mathlib.NumberTheory.FLT.Three`). So
this file does not redo the descent. The headline theorems are stated in the library's form as
corollaries of the mathlib theorem, and the parts of the original that are statements in their
own right are kept, on mathlib's structures.

## Main results

Corollaries of mathlib's `fermatLastTheoremThree`:

* `Fermat3.flt3`, `Fermat3.flt3_nat`: no solutions in nonzero integers, or naturals.
* `Fermat3.flt3_multAdd`: for nonzero integers, `(y / x) ^ 3` is never the mult-add cube of
  `z / x - 1`.
* `Fermat3.no_sol`: `x ^ 3 + y ^ 3 + z ^ 3 = 0` has no primitive nontrivial solution.

The mult-add angle:

* `Fermat3.mpow_three`: `mpow s 3 = s ^ 3 + 3 * (s + s ^ 2)`, the middle layers of the cube carry
  the factor `3`.
* `Fermat3.flt3_multAdd_iff`: in a field, `x ^ 3 + y ^ 3 = z ^ 3` iff
  `(y / x) ^ 3 = mpow (z / x - 1) 3`.
* `Fermat3.mpow_three_mod3`, `Fermat3.sum_mod3`: mult-add cubes and plain cubes agree mod `3`,
  so every solution has `x + y ≡ z` mod `3`.

Elementary facts about cubes, proved here: `case1` (the first case, by cubes mod `9`; mathlib
has it too, but private to its file), `cube_mod9`, `cube_mod9_all`, `cube_ne_two`, `cube_mod2`,
`cube_mod3`, `cube_of_mul3`, `sum_cubes`, `second_factor`, `gcd_factors`, `cop_of_cubes`,
`exactly_one_three`.

The Eisenstein integers and Euler's lemma (this part is not in mathlib in this form):

* `Fermat3.Eis`: the Eisenstein integers `ℤ[ω]`, `ω ^ 2 = -1 - ω`, as mathlib's
  `QuadraticAlgebra ℤ (-1) (-1)`. The ring structure, norm and conjugation (`star`) come from
  mathlib.
* `Fermat3.Eis.instEuclideanDomain`: `ℤ[ω]` is a Euclidean domain for the norm
  (`Fermat3.Eis.norm_mod_lt` is division with remainder). Bezout and unique factorisation then
  come from mathlib.
* `Fermat3.Eis.isUnit_iff_norm_eq_one`, `Fermat3.Eis.units_list`: the six units.
* `Fermat3.key_cube`: Euler's lemma. If `p ^ 2 + 3 * q ^ 2` is a cube, `p` and `q` are coprime
  and `p + q` is odd, then `p = a ^ 3 - 9 * a * b ^ 2` and `q = 3 * a ^ 2 * b - 3 * b ^ 3`, that
  is, `p + q √-3 = (a + b √-3) ^ 3`.

Mathlib's proof of the exponent `3` case works in `𝓞 (CyclotomicField 3 ℚ)` with `λ = ζ - 1` and
does not go through Euler's lemma, so `key_cube` and the Euclidean structure are kept here with
proofs; the descent step built on them in the original is not ported.

Coprimality of integers is stated with `IsCoprime` where the original used `Int.gcd _ _ = 1`
(`Int.isCoprime_iff_gcd_eq_one`), and parity with `Even`/`Odd` where the original used `% 2`.
-/

@[expose] public section

namespace Fermat3

open MultAdd

/-! ### The first case: cubes mod 9 -/

/-- A cube of a number prime to `3` is `1` or `8` mod `9`. -/
theorem cube_mod9 (x : ℤ) (h : ¬ (3 : ℤ) ∣ x) : x ^ 3 % 9 = 1 ∨ x ^ 3 % 9 = 8 := by
  have h9 : x ^ 3 % 9 = (x % 9) ^ 3 % 9 := (Int.ModEq.pow 3 (Int.mod_modEq x 9)).symm
  have : x % 9 = 1 ∨ x % 9 = 2 ∨ x % 9 = 4 ∨ x % 9 = 5 ∨ x % 9 = 7 ∨ x % 9 = 8 := by lia
  rcases this with e | e | e | e | e | e <;> rw [h9, e] <;> decide

/-- The first case: if `3` does not divide `x * y * z`, then `x ^ 3 + y ^ 3 ≠ z ^ 3`, because
cubes of numbers prime to `3` are `1` or `8` mod `9`. Mathlib has this as
`fermatLastTheoremThree_case_1`, which is private to `Mathlib.NumberTheory.FLT.Three`. -/
theorem case1 (x y z : ℤ) (h : ¬ (3 : ℤ) ∣ x * y * z) : x ^ 3 + y ^ 3 ≠ z ^ 3 := by
  have hx := cube_mod9 x fun d ↦ h ((d.mul_right y).mul_right z)
  have hy := cube_mod9 y fun d ↦ h ((d.mul_left x).mul_right z)
  have hz := cube_mod9 z fun d ↦ h (d.mul_left (x * y))
  lia

/-- The mult-add cube: `mpow s 3 = s ^ 3 + 3 * (s + s ^ 2)`; the middle layers carry the `3`. -/
theorem mpow_three {R : Type*} [CommRing R] (s : R) : mpow s 3 = s ^ 3 + 3 * (s + s ^ 2) := by
  rw [mpow_eq]; ring

/-- In a field, with `x ≠ 0`: `x ^ 3 + y ^ 3 = z ^ 3` iff `(y / x) ^ 3` is the mult-add cube of
`z / x - 1`. -/
theorem flt3_multAdd_iff {K : Type*} [Field K] (x y z : K) (hx : x ≠ 0) :
    x ^ 3 + y ^ 3 = z ^ 3 ↔ (y / x) ^ 3 = mpow (z / x - 1) 3 := by
  rw [mpow_eq, add_sub_cancel, div_pow, div_pow, eq_sub_iff_add_eq, div_add_one (pow_ne_zero 3 hx),
    div_left_inj₀ (pow_ne_zero 3 hx), add_comm]

/-- Fermat mod `3`: `x ^ 3 ≡ x`. -/
theorem cube_mod3 (x : ℤ) : (3 : ℤ) ∣ x ^ 3 - x := by
  have h : ∀ a : ZMod 3, a ^ 3 - a = 0 := by decide
  have h3 : ((x ^ 3 - x : ℤ) : ZMod 3) = 0 := by push_cast; exact h x
  exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd _ 3).1 h3

/-- Every solution of `x ^ 3 + y ^ 3 = z ^ 3` has `x + y ≡ z` mod `3`. -/
theorem sum_mod3 (x y z : ℤ) (h : x ^ 3 + y ^ 3 = z ^ 3) : (3 : ℤ) ∣ x + y - z := by
  have e : x + y - z = (z ^ 3 - z) - (x ^ 3 - x) - (y ^ 3 - y) := by linear_combination h
  rw [e]
  exact dvd_sub (dvd_sub (cube_mod3 z) (cube_mod3 x)) (cube_mod3 y)

/-- Mult-add cubes and plain cubes agree mod `3`: `(1 + s) ^ 3 - 1 ≡ s ^ 3 ≡ s`. -/
theorem mpow_three_mod3 (s : ℤ) : (3 : ℤ) ∣ mpow s 3 - s := by
  have e : mpow s 3 - s = (s ^ 3 - s) + 3 * (s + s ^ 2) := by rw [mpow_three]; ring
  rw [e]
  exact dvd_add (cube_mod3 s) (dvd_mul_right 3 _)

/-! ### Cubes and coprime factors -/

/-- Three pairwise coprime factors of a cube are cubes. (For two factors this is mathlib's
`Int.eq_pow_of_mul_eq_pow_odd`.) -/
theorem cube_of_mul3 {a b c w : ℤ} (h : a * b * c = w ^ 3) (hab : IsCoprime a b)
    (hac : IsCoprime a c) (hbc : IsCoprime b c) :
    (∃ d, a = d ^ 3) ∧ (∃ d, b = d ^ 3) ∧ (∃ d, c = d ^ 3) := by
  have h' : a * (b * c) = w ^ 3 := by rw [← h, mul_assoc]
  have hodd : Odd 3 := by decide
  obtain ⟨e, he⟩ := Int.eq_pow_of_mul_eq_pow_odd_right (hab.mul_right hac) hodd h'
  exact ⟨Int.eq_pow_of_mul_eq_pow_odd_left (hab.mul_right hac) hodd h',
    Int.eq_pow_of_mul_eq_pow_odd_left hbc hodd he, Int.eq_pow_of_mul_eq_pow_odd_right hbc hodd he⟩

/-- Cubes mod `9` are `0`, `1` or `8`. -/
theorem cube_mod9_all (x : ℤ) : x ^ 3 % 9 = 0 ∨ x ^ 3 % 9 = 1 ∨ x ^ 3 % 9 = 8 := by
  have h9 : x ^ 3 % 9 = (x % 9) ^ 3 % 9 := (Int.ModEq.pow 3 (Int.mod_modEq x 9)).symm
  have : x % 9 = 0 ∨ x % 9 = 1 ∨ x % 9 = 2 ∨ x % 9 = 3 ∨ x % 9 = 4 ∨ x % 9 = 5 ∨
    x % 9 = 6 ∨ x % 9 = 7 ∨ x % 9 = 8 := by lia
  rcases this with e | e | e | e | e | e | e | e | e <;> rw [h9, e] <;> decide

/-- `2` and `-2` are not cubes (cubes mod `9` are `0`, `1`, `8`). -/
theorem cube_ne_two (x : ℤ) : x ^ 3 ≠ 2 ∧ x ^ 3 ≠ -2 := by
  have := cube_mod9_all x
  constructor <;> intro e <;> rw [e] at this <;> revert this <;> decide

/-- A cube has the parity of its root. -/
theorem cube_mod2 (x : ℤ) : x ^ 3 % 2 = x % 2 := by
  have h2 : x ^ 3 % 2 = (x % 2) ^ 3 % 2 := (Int.ModEq.pow 3 (Int.mod_modEq x 2)).symm
  rcases Int.emod_two_eq_zero_or_one x with e | e <;> rw [h2, e] <;> decide

/-- In a solution of `x ^ 3 + y ^ 3 + z ^ 3 = 0`, if `x` and `y` are coprime then so are `y` and
`z`. -/
theorem cop_of_cubes {x y z : ℤ} (heq : x ^ 3 + y ^ 3 + z ^ 3 = 0) (hxy : IsCoprime x y) :
    IsCoprime y z := by
  have h1 : IsCoprime y (-x ^ 3 + y * (-y ^ 2)) :=
    (hxy.symm.pow_right (n := 3)).neg_right.add_mul_left_right _
  have e : -x ^ 3 + y * (-y ^ 2) = z ^ 3 := by linear_combination (-1 : ℤ) * heq
  rw [e] at h1
  exact (IsCoprime.pow_right_iff (by decide)).1 h1

/-- The sum of two cubes factors: `x ^ 3 + y ^ 3 = (x + y) * (x ^ 2 - x * y + y ^ 2)`. -/
theorem sum_cubes {R : Type*} [CommRing R] (x y : R) :
    x ^ 3 + y ^ 3 = (x + y) * (x ^ 2 - x * y + y ^ 2) := by ring

/-- The second factor is `(x + y) ^ 2 - 3 * (x * y)`. -/
theorem second_factor {R : Type*} [CommRing R] (x y : R) :
    x ^ 2 - x * y + y ^ 2 = (x + y) ^ 2 - 3 * (x * y) := by ring

/-- For coprime `x`, `y` the two factors of `x ^ 3 + y ^ 3` share only `1` or `3`. -/
theorem gcd_factors {x y : ℤ} (h : IsCoprime x y) :
    Int.gcd (x + y) (x ^ 2 - x * y + y ^ 2) = 1 ∨ Int.gcd (x + y) (x ^ 2 - x * y + y ^ 2) = 3 := by
  obtain ⟨u, v, huv⟩ := h
  have c1 : IsCoprime (x + y) x := ⟨v, u - v, by linear_combination huv⟩
  have c2 : IsCoprime (x + y) y := ⟨u, v - u, by linear_combination huv⟩
  have hd1 := Int.gcd_dvd_left (x + y) (x ^ 2 - x * y + y ^ 2)
  have hd2 := Int.gcd_dvd_right (x + y) (x ^ 2 - x * y + y ^ 2)
  generalize Int.gcd (x + y) (x ^ 2 - x * y + y ^ 2) = d at hd1 hd2
  have hcop : IsCoprime (d : ℤ) (x * y) := (c1.mul_right c2).of_isCoprime_of_dvd_left hd1
  have h3 : (d : ℤ) ∣ 3 * (x * y) := by
    have e : 3 * (x * y) = (x + y) * (x + y) - (x ^ 2 - x * y + y ^ 2) := by ring
    rw [e]
    exact dvd_sub (hd1.mul_right _) hd2
  have hd3 : d ∣ 3 := Int.natCast_dvd_natCast.1 (hcop.dvd_of_dvd_mul_right h3)
  exact (Nat.dvd_prime Nat.prime_three).1 hd3

/-- The second case: in a coprime solution with `3 ∣ x * y * z`, exactly one of `x`, `y`, `z` is
divisible by `3`. -/
theorem exactly_one_three {x y z : ℤ} (h : x ^ 3 + y ^ 3 = z ^ 3) (hxy : IsCoprime x y)
    (h3 : (3 : ℤ) ∣ x * y * z) :
    ((3 : ℤ) ∣ x ∧ ¬ (3 : ℤ) ∣ y ∧ ¬ (3 : ℤ) ∣ z) ∨
    (¬ (3 : ℤ) ∣ x ∧ (3 : ℤ) ∣ y ∧ ¬ (3 : ℤ) ∣ z) ∨
    (¬ (3 : ℤ) ∣ x ∧ ¬ (3 : ℤ) ∣ y ∧ (3 : ℤ) ∣ z) := by
  have hyz : IsCoprime y z :=
    (IsCoprime.neg_right_iff y z).1 (cop_of_cubes (by linear_combination h) hxy)
  have hxz : IsCoprime x z :=
    (IsCoprime.neg_right_iff x z).1 (cop_of_cubes (by linear_combination h) hxy.symm)
  have no2 : ∀ a b : ℤ, IsCoprime a b → (3 : ℤ) ∣ a → (3 : ℤ) ∣ b → False := by
    intro a b hab ha hb
    have := Int.isUnit_iff.1 (hab.isUnit_of_dvd' ha hb)
    lia
  rcases Int.prime_three.dvd_mul.1 h3 with h1 | h1
  · rcases Int.prime_three.dvd_mul.1 h1 with h2 | h2
    · exact Or.inl ⟨h2, no2 x y hxy h2, no2 x z hxz h2⟩
    · exact Or.inr (Or.inl ⟨fun hx ↦ no2 x y hxy hx h2, h2, no2 y z hyz h2⟩)
  · exact Or.inr (Or.inr ⟨fun hx ↦ no2 x z hxz hx h1, fun hy ↦ no2 y z hyz hy h1, h1⟩)

/-! ### The Eisenstein integers `ℤ[ω]`, `ω ^ 2 = -1 - ω` -/

/-- The Eisenstein integers `a + b * ω` with `ω ^ 2 = -1 - ω`: mathlib's quadratic algebra
`QuadraticAlgebra ℤ (-1) (-1)`. The norm is `QuadraticAlgebra.norm` and the conjugate is `star`. -/
abbrev Eis : Type := QuadraticAlgebra ℤ (-1) (-1)

namespace Eis

open QuadraticAlgebra

/-- The norm of `a + b * ω` is `a ^ 2 - a * b + b ^ 2`. -/
theorem norm_eq (x : Eis) : norm x = x.re * x.re - x.re * x.im + x.im * x.im := by
  rw [norm_def]; ring

/-- Four times the norm is a square plus three times a square. -/
theorem four_mul_norm (x : Eis) : 4 * norm x = (2 * x.re - x.im) ^ 2 + 3 * x.im ^ 2 := by
  rw [norm_eq]; ring

theorem norm_nonneg (x : Eis) : 0 ≤ norm x := by
  nlinarith [four_mul_norm x, sq_nonneg (2 * x.re - x.im), sq_nonneg x.im]

theorem norm_eq_zero_iff {x : Eis} : norm x = 0 ↔ x = 0 := by
  refine ⟨fun h ↦ ?_, fun h ↦ by rw [h, QuadraticAlgebra.norm_zero]⟩
  have h4 := four_mul_norm x
  have him : x.im = 0 := by
    have : x.im ^ 2 = 0 := by nlinarith [sq_nonneg (2 * x.re - x.im), sq_nonneg x.im]
    exact pow_eq_zero_iff two_ne_zero |>.1 this
  have hre : x.re = 0 := by
    have : x.re ^ 2 = 0 := by nlinarith
    exact pow_eq_zero_iff two_ne_zero |>.1 this
  ext <;> simp [him, hre]

theorem norm_pos {x : Eis} (h : x ≠ 0) : 0 < norm x :=
  lt_of_le_of_ne (norm_nonneg x) (Ne.symm (mt norm_eq_zero_iff.1 h))

/-- The Eisenstein integers are an integral domain. -/
instance : IsDomain Eis :=
  QuadraticAlgebra.Int.isDomain_iff.2 fun ⟨r, hr⟩ ↦ by
    rw [discr_def] at hr
    nlinarith [mul_self_nonneg r]

/-- The units of `ℤ[ω]` are the elements of norm `1`. -/
theorem isUnit_iff_norm_eq_one {u : Eis} : IsUnit u ↔ norm u = 1 := by
  rw [isUnit_iff_norm_isUnit, Int.isUnit_iff]
  exact ⟨fun h ↦ h.resolve_right fun h' ↦ by linarith [norm_nonneg u], Or.inl⟩

/-- The six units `±1`, `±ω`, `±(1 + ω)` (and `1 + ω = -ω ^ 2`). -/
theorem units_list {u : Eis} (h : norm u = 1) :
    u = 1 ∨ u = -1 ∨ u = ω ∨ u = -ω ∨ u = 1 + ω ∨ u = -(1 + ω) := by
  have h4 := four_mul_norm u
  obtain ⟨a, b⟩ := u
  rw [h] at h4
  dsimp only at h4
  have hb : -1 ≤ b ∧ b ≤ 1 := by constructor <;> nlinarith [sq_nonneg (2 * a - b)]
  have ha : -1 ≤ a ∧ a ≤ 1 := by constructor <;> nlinarith [sq_nonneg (2 * a - b), sq_nonneg b]
  obtain ⟨hb1, hb2⟩ := hb
  obtain ⟨ha1, ha2⟩ := ha
  interval_cases a <;> interval_cases b <;> simp at h4 <;>
    simp [QuadraticAlgebra.ext_iff]

/-- Division in `ℤ[ω]`, rounding `a * star b / norm b` to the nearest Eisenstein integer. -/
instance : Div Eis :=
  ⟨fun a b ↦ ⟨(2 * (a * star b).re + norm b) / (2 * norm b),
    (2 * (a * star b).im + norm b) / (2 * norm b)⟩⟩

/-- The remainder of the rounding division in `ℤ[ω]`. -/
instance : Mod Eis := ⟨fun a b ↦ a - b * (a / b)⟩

theorem div_def (a b : Eis) :
    a / b = ⟨(2 * (a * star b).re + norm b) / (2 * norm b),
      (2 * (a * star b).im + norm b) / (2 * norm b)⟩ := rfl

theorem mod_def (a b : Eis) : a % b = a - b * (a / b) := rfl

/-- The algebra behind division with remainder: the norm of `a - b * q`, times `4 * norm b`, in
terms of the distance from `a * star b` to `norm b * q`. -/
theorem four_mul_norm_sub_mul (a b : Eis) (k₁ k₂ : ℤ) :
    4 * (norm (a - b * ⟨k₁, k₂⟩) * norm b) =
      (2 * (a * star b).re - 2 * norm b * k₁) ^ 2 -
        (2 * (a * star b).re - 2 * norm b * k₁) * (2 * (a * star b).im - 2 * norm b * k₂) +
          (2 * (a * star b).im - 2 * norm b * k₂) ^ 2 := by
  simp only [norm_def, re_sub, im_sub, re_mul, im_mul, re_star, im_star]
  ring

/-- Division with remainder: the Eisenstein integers are Euclidean for the norm. -/
theorem norm_mod_lt (a : Eis) {b : Eis} (hb : b ≠ 0) : norm (a % b) < norm b := by
  have hn := norm_pos hb
  have key := four_mul_norm_sub_mul a b ((2 * (a * star b).re + norm b) / (2 * norm b))
    ((2 * (a * star b).im + norm b) / (2 * norm b))
  rw [← div_def, ← mod_def] at key
  have hn2 : 0 < 2 * norm b := by linarith
  have e1 := Int.emod_add_mul_ediv (2 * (a * star b).re + norm b) (2 * norm b)
  have e2 := Int.emod_add_mul_ediv (2 * (a * star b).im + norm b) (2 * norm b)
  have m1 := Int.emod_nonneg (2 * (a * star b).re + norm b) hn2.ne'
  have m2 := Int.emod_nonneg (2 * (a * star b).im + norm b) hn2.ne'
  have l1 := Int.emod_lt_of_pos (2 * (a * star b).re + norm b) hn2
  have l2 := Int.emod_lt_of_pos (2 * (a * star b).im + norm b) hn2
  generalize (2 * (a * star b).re + norm b) / (2 * norm b) = k₁ at *
  generalize (2 * (a * star b).im + norm b) / (2 * norm b) = k₂ at *
  generalize (2 * (a * star b).re + norm b) % (2 * norm b) = M₁ at *
  generalize (2 * (a * star b).im + norm b) % (2 * norm b) = M₂ at *
  generalize hS : 2 * (a * star b).re - 2 * norm b * k₁ = S at key
  generalize hT : 2 * (a * star b).im - 2 * norm b * k₂ = T at key
  generalize norm b = n at *
  generalize norm (a % b) = r at *
  have hS1 : -n ≤ S := by linarith
  have hS2 : S < n := by linarith
  have hT1 : -n ≤ T := by linarith
  have hT2 : T < n := by linarith
  have p1 : 0 ≤ (n - S) * (n + S) := mul_nonneg (by linarith) (by linarith)
  have p2 : 0 ≤ (n - T) * (n + T) := mul_nonneg (by linarith) (by linarith)
  have p3 : 0 ≤ (n - S) * (n - T) := mul_nonneg (by linarith) (by linarith)
  have p4 : 0 ≤ (n + S) * (n + T) := mul_nonneg (by linarith) (by linarith)
  have hle : 4 * (r * n) ≤ 3 * (n * n) := by nlinarith
  by_contra hcon
  have : n * n ≤ r * n := mul_le_mul_of_nonneg_right (not_lt.1 hcon) hn.le
  nlinarith

/-- The Eisenstein integers are a Euclidean domain for the norm. Mathlib then supplies Bezout's
identity, greatest common divisors and unique factorisation. -/
instance instEuclideanDomain : EuclideanDomain Eis :=
  { (inferInstance : CommRing Eis), (inferInstance : Nontrivial Eis) with
    quotient := (· / ·)
    remainder := (· % ·)
    quotient_zero := fun a ↦ by ext <;> simp [div_def]
    quotient_mul_add_remainder_eq := fun a b ↦ by rw [mod_def]; ring
    r := fun a b ↦ (norm a).natAbs < (norm b).natAbs
    r_wellFounded := (measure fun a : Eis ↦ (norm a).natAbs).wf
    remainder_lt := fun a b hb ↦
      Int.natAbs_lt_natAbs_of_nonneg_of_lt (norm_nonneg _) (norm_mod_lt a hb)
    mul_left_not_lt := fun a b hb ↦ not_lt_of_ge <| by
      rw [map_mul, Int.natAbs_mul]
      exact Nat.le_mul_of_pos_right _ (Int.natAbs_pos.2 (norm_pos hb).ne') }

/-- `ω` is a cube root of unity. -/
theorem omega_pow_three : (ω : Eis) ^ 3 = 1 := by
  ext <;> simp [pow_succ]

/-- Every cube is the cube of something with even `ω`-part (something in `ℤ[√-3]`). -/
theorem even_cube_root (δ : Eis) : ∃ δ' : Eis, δ' ^ 3 = δ ^ 3 ∧ Even δ'.im := by
  have hw : ∀ z : Eis, (ω * z) ^ 3 = z ^ 3 := fun z ↦ by
    rw [mul_pow, omega_pow_three, one_mul]
  by_cases h1 : Even δ.im
  · exact ⟨δ, rfl, h1⟩
  by_cases h2 : Even δ.re
  · refine ⟨ω * (ω * δ), by rw [hw, hw], ?_⟩
    simpa using h2
  · refine ⟨ω * δ, hw δ, ?_⟩
    rw [Int.not_even_iff_odd] at h1 h2
    simpa [sub_eq_add_neg] using h2.sub_odd h1

/-- The cube of `x + 2 * k * ω = (x - k) + k * √-3`, written out. -/
theorem cube_sqrt3 (x k : ℤ) :
    (⟨x, 2 * k⟩ : Eis) ^ 3 =
      ⟨((x - k) ^ 3 - 9 * (x - k) * k ^ 2) + (3 * (x - k) ^ 2 * k - 3 * k ^ 3),
        2 * (3 * (x - k) ^ 2 * k - 3 * k ^ 3)⟩ := by
  ext <;> simp [pow_succ] <;> ring

end Eis

/-! ### Euler's lemma: when `p ^ 2 + 3 * q ^ 2` is a cube, `p + q √-3` is a cube -/

theorem odd_norm {p q : ℤ} (hpar : Odd (p + q)) : Odd (p ^ 2 + 3 * q ^ 2) := by
  have e : p ^ 2 + 3 * q ^ 2 = (p + q) ^ 2 + 2 * (q ^ 2 - p * q) := by ring
  rw [e]
  exact hpar.pow.add_even (even_two_mul _)

/-- The integer heart of the coprimality of `p + q √-3` and `p - q √-3`: a common divisor of
`p ^ 2 + 3 * q ^ 2`, `4 * p ^ 2` and `12 * q ^ 2` is `1`. -/
theorem eq_one_of_dvd_sq_add_three_mul_sq {n p q c : ℤ} (hn : 0 ≤ n)
    (h : p ^ 2 + 3 * q ^ 2 = c ^ 3) (hpq : IsCoprime p q) (hpar : Odd (p + q))
    (h1 : n ∣ p ^ 2 + 3 * q ^ 2) (h2 : n ∣ 4 * p ^ 2) (h3 : n ∣ 12 * q ^ 2) : n = 1 := by
  have hodd : Odd n := by
    rw [← Int.not_even_iff_odd]
    exact fun he ↦ Int.not_even_iff_odd.2 (odd_norm hpar) (he.trans_dvd h1)
  obtain ⟨m, hm⟩ := hodd
  have hn4 : IsCoprime n 4 := ⟨1 - 2 * m, m * m, by rw [hm]; ring⟩
  have hp2 : n ∣ p ^ 2 := hn4.dvd_of_dvd_mul_left h2
  have hq2 : IsCoprime n (q ^ 2) := (hpq.pow (m := 2) (n := 2)).of_isCoprime_of_dvd_left hp2
  have h12 : n ∣ 12 := hq2.dvd_of_dvd_mul_right h3
  have h3' : n ∣ 3 := hn4.dvd_of_dvd_mul_left (by simpa using h12)
  have hle : n ≤ 3 := Int.le_of_dvd (by decide) h3'
  interval_cases n
  · simp at h3'
  · rfl
  · lia
  · exfalso
    obtain ⟨p', rfl⟩ := Int.prime_three.dvd_of_dvd_pow hp2
    obtain ⟨c', rfl⟩ := Int.prime_three.dvd_of_dvd_pow (n := 3) (a := c)
      ⟨p' ^ 2 * 3 + q ^ 2, by rw [← h]; ring⟩
    have hq : (3 : ℤ) ∣ q := Int.prime_three.dvd_of_dvd_pow (n := 2)
      ⟨3 * c' ^ 3 - p' ^ 2, by linarith⟩
    have := Int.isUnit_iff.1 (hpq.isUnit_of_dvd' (dvd_mul_right 3 p') hq)
    lia

open QuadraticAlgebra in
/-- Euler's lemma: if `p ^ 2 + 3 * q ^ 2` is a cube, `p` and `q` are coprime and `p + q` is odd,
then `p + q √-3` is a cube: `p = a ^ 3 - 9 * a * b ^ 2` and `q = 3 * a ^ 2 * b - 3 * b ^ 3`. The
proof works in the Eisenstein integers, where `p + q √-3 = (p + q) + 2 * q * ω` and its conjugate
are coprime with product a cube, so each is a unit times a cube. -/
theorem key_cube {p q c : ℤ} (h : p ^ 2 + 3 * q ^ 2 = c ^ 3) (hpq : IsCoprime p q)
    (hpar : Odd (p + q)) :
    ∃ a b : ℤ, p = a ^ 3 - 9 * a * b ^ 2 ∧ q = 3 * a ^ 2 * b - 3 * b ^ 3 := by
  obtain ⟨α, hα⟩ : ∃ α : Eis, α = ⟨p + q, 2 * q⟩ := ⟨_, rfl⟩
  have hN : norm α = p ^ 2 + 3 * q ^ 2 := by rw [hα, Eis.norm_eq]; ring
  have hprod : α * star α = (c : Eis) ^ 3 := by
    rw [← algebraMap_norm_eq_mul_star, hN, h, map_pow, eq_intCast]
  have hcop : IsCoprime α (star α) := by
    refine isCoprime_of_dvd α (star α) ?_ ?_
    · rintro ⟨h0, -⟩
      have h1 := congrArg QuadraticAlgebra.re h0
      rw [hα] at h1
      obtain ⟨m, hm⟩ := hpar
      simp at h1
      lia
    · intro d hd _ hd1 hd2
      apply hd
      refine Eis.isUnit_iff_norm_eq_one.2 (eq_one_of_dvd_sq_add_three_mul_sq (Eis.norm_nonneg d)
        h hpq hpar ?_ ?_ ?_)
      · rw [← hN]; exact map_dvd QuadraticAlgebra.norm hd1
      · have e : norm (α + star α) = 4 * p ^ 2 := by
          rw [hα, Eis.norm_eq]; simp; ring
        rw [← e]; exact map_dvd QuadraticAlgebra.norm (dvd_add hd1 hd2)
      · have e : norm (α - star α) = 12 * q ^ 2 := by
          rw [hα, Eis.norm_eq]; simp; ring
        rw [← e]; exact map_dvd QuadraticAlgebra.norm (dvd_sub hd1 hd2)
  obtain ⟨δ, u, hu⟩ := exists_associated_pow_of_mul_eq_pow' hcop hprod
  obtain ⟨δ', hδ', k, hk⟩ := Eis.even_cube_root δ
  have hδ'' : δ' = ⟨δ'.re, 2 * k⟩ := by ext <;> simp [hk, two_mul]
  generalize δ'.re = x at hδ''
  rw [← hδ', hδ'', Eis.cube_sqrt3, hα] at hu
  generalize hP : (x - k) ^ 3 - 9 * (x - k) * k ^ 2 = P at hu
  generalize hQ : 3 * (x - k) ^ 2 * k - 3 * k ^ 3 = Q at hu
  have hre := congrArg QuadraticAlgebra.re hu
  have him := congrArg QuadraticAlgebra.im hu
  have eP : (-(x - k)) ^ 3 - 9 * (-(x - k)) * (-k) ^ 2 = -P := by rw [← hP]; ring
  have eQ : 3 * (-(x - k)) ^ 2 * (-k) - 3 * (-k) ^ 3 = -Q := by rw [← hQ]; ring
  obtain ⟨m, hm⟩ := hpar
  rcases Eis.units_list (Eis.isUnit_iff_norm_eq_one.1 u.isUnit) with e | e | e | e | e | e <;>
    rw [e] at hre him <;> simp at hre him
  · exact ⟨x - k, k, by lia, by lia⟩
  · exact ⟨-(x - k), -k, by lia, by lia⟩
  all_goals lia

/-- If `p + q` is odd then `a + b` is odd in Euler's parametrisation. -/
theorem par_ab {a b : ℤ} (h : Odd (a ^ 3 - 9 * a * b ^ 2 + (3 * a ^ 2 * b - 3 * b ^ 3))) :
    Odd (a + b) := by
  by_contra hc
  rw [Int.not_odd_iff_even] at hc
  have e : a ^ 3 - 9 * a * b ^ 2 + (3 * a ^ 2 * b - 3 * b ^ 3) =
      (a + b) ^ 3 + 2 * (-6 * a * b ^ 2 - 2 * b ^ 3) := by ring
  rw [e] at h
  exact Int.not_even_iff_odd.2 h ((hc.pow_of_ne_zero three_ne_zero).add (even_two_mul _))

/-- In Euler's parametrisation, coprime `p`, `q` come from coprime `a`, `b`. -/
theorem cop_of_param {a b p q : ℤ} (hpq : IsCoprime p q) (hp : p = a ^ 3 - 9 * a * b ^ 2)
    (hq : q = 3 * a ^ 2 * b - 3 * b ^ 3) : IsCoprime a b := by
  have h1 : p = a * (a ^ 2 - 9 * b ^ 2) := by rw [hp]; ring
  have h2 : q = b * (3 * a ^ 2 - 3 * b ^ 2) := by rw [hq]; ring
  rw [h1, h2] at hpq
  exact hpq.of_mul_left_left.of_mul_right_left

/-! ### Fermat's Last Theorem for exponent 3, from mathlib -/

/-- Fermat's Last Theorem for exponent 3, over the integers: no nonzero `x`, `y`, `z` with
`x ^ 3 + y ^ 3 = z ^ 3`. This is mathlib's `fermatLastTheoremThree`. -/
theorem flt3 (x y z : ℤ) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0) : x ^ 3 + y ^ 3 ≠ z ^ 3 :=
  fermatLastTheoremFor_iff_int.1 fermatLastTheoremThree x y z hx hy hz

/-- Fermat's Last Theorem for exponent 3, over the naturals: `x ^ 3 + y ^ 3 ≠ z ^ 3` when
`x, y ≠ 0`. This is mathlib's `fermatLastTheoremThree`. -/
theorem flt3_nat (x y z : ℕ) (hx : x ≠ 0) (hy : y ≠ 0) : x ^ 3 + y ^ 3 ≠ z ^ 3 := by
  intro h
  refine fermatLastTheoremThree x y z hx hy ?_ h
  rintro rfl
  have : 0 < x ^ 3 := by positivity
  lia

/-- In the mult-add picture: for nonzero integers `x`, `y`, `z`, the cube `(y / x) ^ 3` is never
the mult-add cube of `z / x - 1`. A corollary of mathlib's `fermatLastTheoremThree`. -/
theorem flt3_multAdd (x y z : ℤ) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0) :
    ((y : ℚ) / x) ^ 3 ≠ mpow ((z : ℚ) / x - 1) 3 := by
  intro h
  have := (flt3_multAdd_iff (x : ℚ) y z (mod_cast hx)).2 h
  exact flt3 x y z hx hy hz (mod_cast this)

/-- A primitive nontrivial solution of `x ^ 3 + y ^ 3 + z ^ 3 = 0`. -/
def Sol (x y z : ℤ) : Prop :=
  x ≠ 0 ∧ y ≠ 0 ∧ z ≠ 0 ∧ x ^ 3 + y ^ 3 + z ^ 3 = 0 ∧ IsCoprime x y

/-- No primitive nontrivial solution of `x ^ 3 + y ^ 3 + z ^ 3 = 0` exists. A corollary of
mathlib's `fermatLastTheoremThree`. -/
theorem no_sol (x y z : ℤ) : ¬ Sol x y z := fun ⟨hx, hy, hz, h, _⟩ ↦
  flt3 x y (-z) hx hy (neg_ne_zero.2 hz) (by linear_combination h)

end Fermat3
