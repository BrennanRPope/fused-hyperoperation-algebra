/-
  Fermat's Last Theorem for exponent 3, in core Lean 4 (no mathlib, no `sorry`)
  ===========================================================================

  WHAT IS PROVED (everything below; nothing is left as a hypothesis):

  * `flt3`      : for integers x, y, z, all nonzero, x³ + y³ ≠ z³.
  * `flt3_nat`  : for natural numbers with x, y ≠ 0, x³ + y³ ≠ z³.
  * `case1`     : the first case on its own: if 3 ∤ xyz then x³ + y³ ≠ z³,
                  because cubes of numbers prime to 3 are 1 or 8 mod 9.

  Supporting lemmas, each fully proved:
  * `primitive_of_sol`  : any solution can be divided down to one with gcd(x, y) = 1.
  * `exactly_one_three` : in a coprime solution with 3 ∣ xyz, exactly one of x, y, z
                          is divisible by 3.
  * `sum_cubes`, `gcd_factors` : x³ + y³ = (x + y)(x² − xy + y²), and for coprime
                          x, y the two factors have gcd 1 or 3.
  * `key_cube` (Euler's lemma): if p² + 3q² is a cube, gcd(p, q) = 1 and p + q is odd,
                          then p = a³ − 9ab², q = 3a²b − 3b³, i.e. p + q√−3 = (a + b√−3)³.
                          Proved in the Eisenstein integers ℤ[ω], which are built here
                          from scratch: ring, norm, Euclidean division, Bezout, units,
                          and "a coprime product that is a cube is a unit times a cube".
  * `step`, `no_sol`    : Euler's infinite descent on |xyz| for x³ + y³ + z³ = 0.

  The mult-add angle (a ⊞ b = a + b + ab, mpow s n = (1 + s)ⁿ − 1, from FusedAlgebra.Base):
  * `flt3_circ_iff` : over ℚ, with x ≠ 0, x³ + y³ = z³ iff (y/x)³ = mpow (z/x − 1) 3.
  * `flt3_circ`     : so for nonzero integers, (y/x)³ is never the mult-add cube of z/x − 1.
  * `mpow_three`       : mpow s 3 = s³ + 3(s + s²): the middle layers of the cube
                         (Pascal 1, 3, 3, 1) carry the factor 3,
  * `mpow_three_mod3`  : so mult-add cubes and plain cubes agree mod 3 (both ≡ s),
  * `sum_mod3`         : which is why every solution has x + y ≡ z (mod 3).

  Route of the proof of `flt3`: reduce to gcd(x, y) = 1; write the equation as
  x³ + y³ + z³ = 0 and rotate so z is even; set x = p + q, y = p − q, so that
  2p(p² + 3q²) = (−z)³. If 3 ∤ p the two factors are coprime cubes, Euler's lemma gives
  p = a(a − 3b)(a + 3b), and 2a, a − 3b, a + 3b are coprime cubes with
  (a − 3b) + (a + 3b) − 2a = 0: a smaller solution. If p = 3r the same happens with
  18r(q² + 3r²), giving 2b, a − b, a + b. The descent measure |xyz| drops each time.

  Only Lean's standard axioms are used (propext, Classical.choice, Quot.sound); the
  `#print axioms` lines at the end check this.

  LIBRARY NOTE: copied from fermat/lean/Fermat3.lean (Nov 25) as Topics/Fermat/Fermat3.
-/
import FusedAlgebra.Base
open CircleOp

namespace Fermat3

/-! ## 1. The first case: cubes mod 9 -/

theorem pow_emod_of (x m : Int) (n : Nat) : x^n % m = ((x % m)^n) % m := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Int.pow_succ, Int.pow_succ, Int.mul_emod, ih, Int.mul_emod ((x % m)^n), Int.emod_emod]

/-- A cube of a number prime to 3 is 1 or 8 mod 9. -/
theorem cube_mod9 (x : Int) (h : ¬ (3 : Int) ∣ x) : x^3 % 9 = 1 ∨ x^3 % 9 = 8 := by
  have h3 : x % 3 ≠ 0 := fun e => h (Int.dvd_of_emod_eq_zero e)
  have : x % 9 = 1 ∨ x % 9 = 2 ∨ x % 9 = 4 ∨ x % 9 = 5 ∨ x % 9 = 7 ∨ x % 9 = 8 := by omega
  rw [pow_emod_of]
  rcases this with e|e|e|e|e|e <;> rw [e] <;> decide

/-- **First case**: if 3 does not divide xyz, then x³ + y³ ≠ z³. -/
theorem case1 (x y z : Int) (h : ¬ (3 : Int) ∣ x * y * z) : x^3 + y^3 ≠ z^3 := by
  have hx : ¬ (3 : Int) ∣ x := fun d => h (Int.dvd_mul_of_dvd_left (Int.dvd_mul_of_dvd_left d))
  have hy : ¬ (3 : Int) ∣ y := fun d => h (Int.dvd_mul_of_dvd_left (Int.dvd_mul_of_dvd_right d))
  have hz : ¬ (3 : Int) ∣ z := fun d => h (Int.dvd_mul_of_dvd_right d)
  have a := cube_mod9 x hx
  have b := cube_mod9 y hy
  have c := cube_mod9 z hz
  intro e
  generalize x^3 = X at *
  generalize y^3 = Y at *
  generalize z^3 = Z at *
  omega

/-- The mult-add cube: mpow s 3 = s³ + 3(s + s²); the middle layers carry the 3. -/
theorem mpow_three (s : Int) : mpow s 3 = s^3 + 3*(s + s^2) := by
  rw [mpow_eq]; grind

/-- Over ℚ with x ≠ 0: x³ + y³ = z³ iff (y/x)³ is the mult-add cube of z/x − 1. -/
theorem flt3_circ_iff (x y z : Rat) (hx : x ≠ 0) :
    x^3 + y^3 = z^3 ↔ (y/x)^3 = mpow (z/x - 1) 3 := by
  rw [mpow_eq]
  constructor
  · intro h; grind
  · intro h; grind

/-- Fermat mod 3: x³ ≡ x. -/
theorem cube_mod3 (x : Int) : (3 : Int) ∣ x^3 - x := by
  have e : x^3 - x = (x - 1) * x * (x + 1) := by grind
  rw [e]
  have : x % 3 = 0 ∨ x % 3 = 1 ∨ x % 3 = 2 := by omega
  rcases this with h|h|h
  · exact Int.dvd_mul_of_dvd_left (Int.dvd_mul_of_dvd_right (Int.dvd_of_emod_eq_zero h))
  · exact Int.dvd_mul_of_dvd_left (Int.dvd_mul_of_dvd_left (Int.dvd_of_emod_eq_zero (by omega)))
  · exact Int.dvd_mul_of_dvd_right (Int.dvd_of_emod_eq_zero (by omega))

/-- Every solution of x³ + y³ = z³ has x + y ≡ z (mod 3). -/
theorem sum_mod3 (x y z : Int) (h : x^3 + y^3 = z^3) : (3 : Int) ∣ x + y - z := by
  have a := cube_mod3 x
  have b := cube_mod3 y
  have c := cube_mod3 z
  have : x + y - z = (z^3 - z) - (x^3 - x) - (y^3 - y) := by grind
  rw [this]
  exact Int.dvd_sub (Int.dvd_sub c a) b

/-! ## 2. Integer gcd tools and "a coprime product of cubes is made of cubes" -/

theorem dvd_of_dvd_mul_cop {k m n : Int} (hc : Int.gcd k m = 1) (h : k ∣ m * n) : k ∣ n := by
  have := (Int.dvd_mul_gcd_iff_dvd_mul (k := k) (n := n) (m := m)).2 (by rwa [Int.mul_comm])
  rwa [hc, Int.natCast_one, Int.mul_one] at this

theorem cop_symm {a b : Int} (h : Int.gcd a b = 1) : Int.gcd b a = 1 := by
  rwa [Int.gcd_comm]

theorem cop_of_dvd_left {a b d : Int} (h : Int.gcd a b = 1) (hd : d ∣ a) : Int.gcd d b = 1 := by
  rw [Int.gcd_eq_one_iff] at h ⊢
  exact fun c hc hb => h c (Int.dvd_trans hc hd) hb

theorem cop_of_dvd_right {a b d : Int} (h : Int.gcd a b = 1) (hd : d ∣ b) : Int.gcd a d = 1 :=
  cop_symm (cop_of_dvd_left (cop_symm h) hd)

theorem cop_mul_right {a b c : Int} (h1 : Int.gcd a b = 1) (h2 : Int.gcd a c = 1) :
    Int.gcd a (b * c) = 1 := by
  rw [Int.gcd_mul_right_right_of_gcd_eq_one h1, h2]

theorem cop_mul_left {a b c : Int} (h1 : Int.gcd a c = 1) (h2 : Int.gcd b c = 1) :
    Int.gcd (a * b) c = 1 :=
  cop_symm (cop_mul_right (cop_symm h1) (cop_symm h2))

theorem cop_of_pow {a b : Int} {n : Nat} (hn : n ≠ 0) (h : Int.gcd (a^n) (b^n) = 1) :
    Int.gcd a b = 1 := by
  rw [Int.pow_gcd_pow] at h
  rcases Nat.pow_eq_one.1 h with h | h
  · exact h
  · exact absurd h hn

theorem cop_two_of_odd {a : Int} (h : a % 2 = 1) : Int.gcd 2 a = 1 := by
  rw [Int.gcd_eq_one_iff]
  intro c h2 ha
  have e : (1 : Int) = a - 2 * (a / 2) := by omega
  rw [e]; exact Int.dvd_sub ha (Int.dvd_mul_of_dvd_left h2)

theorem cop_three {a : Int} (h : a % 3 ≠ 0) : Int.gcd 3 a = 1 := by
  rw [Int.gcd_eq_one_iff]
  intro c h3 ha
  have : a % 3 = 1 ∨ a % 3 = 2 := by omega
  rcases this with e | e
  · have e' : (1 : Int) = a - 3 * (a / 3) := by omega
    rw [e']; exact Int.dvd_sub ha (Int.dvd_mul_of_dvd_left h3)
  · have e' : (1 : Int) = 3 - (a - 3 * (a / 3)) := by omega
    rw [e']; exact Int.dvd_sub h3 (Int.dvd_sub ha (Int.dvd_mul_of_dvd_left h3))

theorem three_dvd_of_sq {x : Int} (h : (3 : Int) ∣ x^2) : (3 : Int) ∣ x := by
  have h' := Int.emod_eq_zero_of_dvd h
  rw [pow_emod_of] at h'
  apply Int.dvd_of_emod_eq_zero
  have : x % 3 = 0 ∨ x % 3 = 1 ∨ x % 3 = 2 := by omega
  rcases this with e|e|e <;> rw [e] at h' ⊢ <;> revert h' <;> decide

theorem three_dvd_of_cube {x : Int} (h : (3 : Int) ∣ x^3) : (3 : Int) ∣ x := by
  have h' := Int.emod_eq_zero_of_dvd h
  rw [pow_emod_of] at h'
  apply Int.dvd_of_emod_eq_zero
  have : x % 3 = 0 ∨ x % 3 = 1 ∨ x % 3 = 2 := by omega
  rcases this with e|e|e <;> rw [e] at h' ⊢ <;> revert h' <;> decide

/-- If a coprime product a·b is a cube, then a is a cube. -/
theorem cube_of_mul_cop {a b c : Int} (h : a * b = c^3) (hc : Int.gcd a b = 1) :
    ∃ d : Int, a = d^3 := by
  have h1 : Int.gcd (a^3) (c^3) = (Int.gcd a c)^3 := Int.pow_gcd_pow
  have h2 : Int.gcd (a^3) (c^3) = a.natAbs := by
    have e1 : a^3 = a * a^2 := by grind
    rw [← h, e1, Int.gcd_mul_left, Int.gcd_pow_left_of_gcd_eq_one hc, Nat.mul_one]
  rw [h2] at h1
  have h3 : (a.natAbs : Int) = (Int.gcd a c : Int)^3 := by rw [h1, Int.natCast_pow]
  rcases Int.natAbs_eq a with e | e
  · rw [h3] at e; exact ⟨_, e⟩
  · rw [h3] at e; exact ⟨-(Int.gcd a c : Int), by rw [e]; grind⟩

theorem cube_of_mul_cop' {a b c : Int} (h : a * b = c^3) (hc : Int.gcd a b = 1) :
    ∃ d : Int, b = d^3 :=
  cube_of_mul_cop (by rw [Int.mul_comm]; exact h) (cop_symm hc)

/-- Three pairwise coprime factors of a cube are cubes. -/
theorem cube_of_mul3 {a b c w : Int} (h : a * b * c = w^3) (hab : Int.gcd a b = 1)
    (hac : Int.gcd a c = 1) (hbc : Int.gcd b c = 1) :
    (∃ d, a = d^3) ∧ (∃ d, b = d^3) ∧ (∃ d, c = d^3) := by
  have h' : a * (b * c) = w^3 := by rw [← h, Int.mul_assoc]
  have hcop : Int.gcd a (b * c) = 1 := cop_mul_right hab hac
  obtain ⟨e, he⟩ := cube_of_mul_cop' h' hcop
  exact ⟨cube_of_mul_cop h' hcop, cube_of_mul_cop he.symm.symm hbc |>.elim (fun d hd => ⟨d, hd⟩),
    cube_of_mul_cop' he hbc⟩

/-! ## 3. The Eisenstein integers ℤ[ω], ω² = −1 − ω -/

/-- Eisenstein integers a + b·ω with ω² = −1 − ω. -/
structure Eis where
  re : Int
  im : Int

namespace Eis

@[ext] theorem ext {x y : Eis} (h1 : x.re = y.re) (h2 : x.im = y.im) : x = y := by
  cases x; cases y; simp only at h1 h2; subst h1; subst h2; rfl

instance : Add Eis := ⟨fun x y => ⟨x.re + y.re, x.im + y.im⟩⟩
instance : Neg Eis := ⟨fun x => ⟨-x.re, -x.im⟩⟩
instance : Sub Eis := ⟨fun x y => ⟨x.re - y.re, x.im - y.im⟩⟩
instance : Mul Eis := ⟨fun x y => ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re - x.im * y.im⟩⟩
instance : NatCast Eis := ⟨fun n => ⟨n, 0⟩⟩
instance : IntCast Eis := ⟨fun n => ⟨n, 0⟩⟩
instance (n : Nat) : OfNat Eis n := ⟨⟨n, 0⟩⟩
instance : SMul Nat Eis := ⟨fun n x => ⟨n * x.re, n * x.im⟩⟩
instance : SMul Int Eis := ⟨fun n x => ⟨n * x.re, n * x.im⟩⟩

def npow (x : Eis) : Nat → Eis
  | 0 => ⟨1, 0⟩
  | n + 1 => npow x n * x

instance : HPow Eis Nat Eis := ⟨npow⟩

@[simp] theorem add_re (x y : Eis) : (x + y).re = x.re + y.re := rfl
@[simp] theorem add_im (x y : Eis) : (x + y).im = x.im + y.im := rfl
@[simp] theorem neg_re (x : Eis) : (-x).re = -x.re := rfl
@[simp] theorem neg_im (x : Eis) : (-x).im = -x.im := rfl
@[simp] theorem sub_re (x y : Eis) : (x - y).re = x.re - y.re := rfl
@[simp] theorem sub_im (x y : Eis) : (x - y).im = x.im - y.im := rfl
@[simp] theorem mul_re (x y : Eis) : (x * y).re = x.re * y.re - x.im * y.im := rfl
@[simp] theorem mul_im (x y : Eis) : (x * y).im = x.re * y.im + x.im * y.re - x.im * y.im := rfl
@[simp] theorem ofNat_re (n : Nat) : (OfNat.ofNat n : Eis).re = n := rfl
@[simp] theorem ofNat_im (n : Nat) : (OfNat.ofNat n : Eis).im = 0 := rfl
@[simp] theorem natCast_re (n : Nat) : (n : Eis).re = n := rfl
@[simp] theorem natCast_im (n : Nat) : (n : Eis).im = 0 := rfl
@[simp] theorem intCast_re (n : Int) : (n : Eis).re = n := rfl
@[simp] theorem intCast_im (n : Int) : (n : Eis).im = 0 := rfl
@[simp] theorem nsmul_re (n : Nat) (x : Eis) : (n • x).re = n * x.re := rfl
@[simp] theorem nsmul_im (n : Nat) (x : Eis) : (n • x).im = n * x.im := rfl
@[simp] theorem zsmul_re (n : Int) (x : Eis) : (n • x).re = n * x.re := rfl
@[simp] theorem zsmul_im (n : Int) (x : Eis) : (n • x).im = n * x.im := rfl
@[simp] theorem zero_re : (0 : Eis).re = 0 := rfl
@[simp] theorem zero_im : (0 : Eis).im = 0 := rfl
@[simp] theorem one_re : (1 : Eis).re = 1 := rfl
@[simp] theorem one_im : (1 : Eis).im = 0 := rfl
theorem pow_zero' (x : Eis) : x ^ 0 = ⟨1, 0⟩ := rfl
theorem pow_succ' (x : Eis) (n : Nat) : x ^ (n + 1) = x ^ n * x := rfl

instance : Lean.Grind.CommRing Eis where
  add_zero := by intro a; ext <;> simp
  add_comm := by intro a b; ext <;> simp <;> omega
  add_assoc := by intro a b c; ext <;> simp <;> omega
  mul_assoc := by intro a b c; ext <;> simp <;> grind
  mul_one := by intro a; ext <;> simp
  one_mul := by intro a; ext <;> simp
  left_distrib := by intro a b c; ext <;> simp <;> grind
  right_distrib := by intro a b c; ext <;> simp <;> grind
  zero_mul := by intro a; ext <;> simp
  mul_zero := by intro a; ext <;> simp
  pow_zero := by intro a; rfl
  pow_succ := by intro a n; rfl
  ofNat_succ := by intro n; ext <;> simp
  ofNat_eq_natCast := by intro n; rfl
  nsmul_eq_natCast_mul := by intro n a; ext <;> simp
  neg_add_cancel := by intro a; ext <;> simp <;> omega
  sub_eq_add_neg := by intro a b; ext <;> simp <;> omega
  neg_zsmul := by intro i a; ext <;> simp <;> grind
  zsmul_natCast_eq_nsmul := by intro n a; ext <;> simp
  intCast_ofNat := by intro n; rfl
  intCast_neg := by intro i; ext <;> simp
  mul_comm := by intro a b; ext <;> simp <;> grind

instance : Dvd Eis := ⟨fun a b => ∃ c, b = a * c⟩

theorem dvd_def {a b : Eis} : a ∣ b ↔ ∃ c, b = a * c := Iff.rfl

theorem dvd_trans {a b c : Eis} (h1 : a ∣ b) (h2 : b ∣ c) : a ∣ c := by
  obtain ⟨x, hx⟩ := h1; obtain ⟨y, hy⟩ := h2
  exact ⟨x * y, by rw [hy, hx]; grind⟩

theorem dvd_add {a b c : Eis} (h1 : a ∣ b) (h2 : a ∣ c) : a ∣ b + c := by
  obtain ⟨x, hx⟩ := h1; obtain ⟨y, hy⟩ := h2
  exact ⟨x + y, by rw [hy, hx]; grind⟩

theorem dvd_sub {a b c : Eis} (h1 : a ∣ b) (h2 : a ∣ c) : a ∣ b - c := by
  obtain ⟨x, hx⟩ := h1; obtain ⟨y, hy⟩ := h2
  exact ⟨x - y, by rw [hy, hx]; grind⟩

theorem dvd_mul_left_of {a b : Eis} (c : Eis) (h : a ∣ b) : a ∣ c * b := by
  obtain ⟨x, hx⟩ := h
  exact ⟨c * x, by rw [hx]; grind⟩

theorem dvd_refl (a : Eis) : a ∣ a := ⟨1, by grind⟩

/-- The norm N(a + bω) = a² − ab + b². -/
def N (x : Eis) : Int := x.re * x.re - x.re * x.im + x.im * x.im

/-- The conjugate of a + bω is (a − b) − bω. -/
def conj (x : Eis) : Eis := ⟨x.re - x.im, -x.im⟩

theorem N_mul (x y : Eis) : N (x * y) = N x * N y := by
  simp only [N, mul_re, mul_im]; grind

theorem N_conj (x : Eis) : N (conj x) = N x := by
  simp only [N, conj]; grind

theorem mul_conj (x : Eis) : x * conj x = ⟨N x, 0⟩ := by
  ext <;> simp [conj, N] <;> grind

theorem sq_nonneg' (x : Int) : 0 ≤ x * x := by
  rcases Int.le_total 0 x with h | h
  · exact Int.mul_nonneg h h
  · exact Int.mul_nonneg_of_nonpos_of_nonpos h h

theorem sq_pos' {x : Int} (h : x ≠ 0) : 0 < x * x := by
  rcases Int.lt_or_gt_of_ne h with h | h
  · exact Int.mul_pos_of_neg_of_neg h h
  · exact Int.mul_pos h h

theorem four_N (x : Eis) : 4 * N x = (2 * x.re - x.im) * (2 * x.re - x.im) + 3 * (x.im * x.im) := by
  simp only [N]; grind

theorem N_nonneg (x : Eis) : 0 ≤ N x := by
  have h1 := sq_nonneg' (2 * x.re - x.im)
  have h2 := sq_nonneg' x.im
  have := four_N x
  omega

theorem eq_zero_of_N {x : Eis} (h : N x = 0) : x = 0 := by
  have h4 := four_N x
  have h1 := sq_nonneg' (2 * x.re - x.im)
  have him : x.im = 0 := by
    apply Classical.byContradiction; intro hne
    have := sq_pos' hne
    omega
  have hre : x.re * x.re = 0 := by simp only [N, him] at h; grind
  rcases Int.mul_eq_zero.1 hre with e | e <;> ext <;> simp_all

theorem N_zero : N 0 = 0 := rfl

theorem mul_eq_zero {x y : Eis} (h : x * y = 0) : x = 0 ∨ y = 0 := by
  have := N_mul x y
  rw [h, N_zero] at this
  rcases Int.mul_eq_zero.1 this.symm with e | e
  · exact Or.inl (eq_zero_of_N e)
  · exact Or.inr (eq_zero_of_N e)

theorem N_pos {x : Eis} (h : x ≠ 0) : 0 < N x := by
  have := N_nonneg x
  have : N x ≠ 0 := fun e => h (eq_zero_of_N e)
  omega

theorem N_dvd {a b : Eis} (h : a ∣ b) : N a ∣ N b := by
  obtain ⟨c, hc⟩ := h
  exact ⟨N c, by rw [hc, N_mul]⟩

/-- A unit: something with an inverse. -/
def IsUnit (u : Eis) : Prop := ∃ v, u * v = 1

theorem isUnit_of_N {u : Eis} (h : N u = 1) : IsUnit u :=
  ⟨conj u, by rw [mul_conj, h]; rfl⟩

theorem N_of_isUnit {u : Eis} (h : IsUnit u) : N u = 1 := by
  obtain ⟨v, hv⟩ := h
  have := N_mul u v
  rw [hv] at this
  exact Int.eq_one_of_mul_eq_one_right (N_nonneg u) this.symm

theorem small_of_sq_le_one {b : Int} (h : b * b ≤ 1) : b = -1 ∨ b = 0 ∨ b = 1 := by
  have h1 := sq_nonneg' (b - 2)
  have h2 := sq_nonneg' (b + 2)
  have : -1 ≤ b ∧ b ≤ 1 := by grind
  omega

/-- The six units ±1, ±ω, ±ω². -/
theorem units_list {u : Eis} (h : N u = 1) :
    u = ⟨1, 0⟩ ∨ u = ⟨-1, 0⟩ ∨ u = ⟨0, 1⟩ ∨ u = ⟨0, -1⟩ ∨ u = ⟨1, 1⟩ ∨ u = ⟨-1, -1⟩ := by
  have h4 := four_N u
  have h1 := sq_nonneg' (2 * u.re - u.im)
  have hb : u.im * u.im ≤ 1 := by omega
  obtain ⟨a, b⟩ := u
  simp only [N] at h h4 hb
  rcases small_of_sq_le_one hb with e | e | e <;> subst e
  · have : (a + 1) * a = 0 := by grind
    rcases Int.mul_eq_zero.1 this with e | e
    · right; right; right; right; right; congr; omega
    · right; right; right; left; congr
  · have : (a + 1) * (a - 1) = 0 := by grind
    rcases Int.mul_eq_zero.1 this with e | e
    · right; left; congr; omega
    · left; congr; omega
  · have : a * (a - 1) = 0 := by grind
    rcases Int.mul_eq_zero.1 this with e | e
    · right; right; left; congr
    · right; right; right; right; left; congr; omega

/-- Division with remainder: the Eisenstein integers are Euclidean for the norm. -/
theorem div_rem (a b : Eis) (hb : b ≠ 0) : ∃ q r, a = q * b + r ∧ N r < N b := by
  have hn := N_pos hb
  generalize hnd : N b = n at hn
  generalize hw : a * conj b = w
  have hbb : b * conj b = ⟨n, 0⟩ := by rw [mul_conj, hnd]
  have hn2 : 0 < 2 * n := by omega
  have e1 := Int.emod_add_mul_ediv (2 * w.re + n) (2 * n)
  have e2 := Int.emod_add_mul_ediv (2 * w.im + n) (2 * n)
  have m1 := Int.emod_nonneg (2 * w.re + n) (Int.ne_of_gt hn2)
  have m2 := Int.emod_nonneg (2 * w.im + n) (Int.ne_of_gt hn2)
  have l1 := Int.emod_lt_of_pos (2 * w.re + n) hn2
  have l2 := Int.emod_lt_of_pos (2 * w.im + n) hn2
  generalize (2 * w.re + n) / (2 * n) = k1 at *
  generalize (2 * w.im + n) / (2 * n) = k2 at *
  generalize (2 * w.re + n) % (2 * n) = M1 at *
  generalize (2 * w.im + n) % (2 * n) = M2 at *
  refine ⟨⟨k1, k2⟩, a - ⟨k1, k2⟩ * b, by grind, ?_⟩
  generalize hr : a - ⟨k1, k2⟩ * b = r
  have hrc : r * conj b = ⟨w.re - n * k1, w.im - n * k2⟩ := by
    have : r * conj b = a * conj b - ⟨k1, k2⟩ * (b * conj b) := by rw [← hr]; grind
    rw [this, hw, hbb]; ext <;> simp <;> grind
  have hN : N r * n = N ⟨w.re - n * k1, w.im - n * k2⟩ := by
    rw [← hrc, N_mul, N_conj, hnd]
  generalize hS : 2 * (w.re - n * k1) = S
  generalize hT : 2 * (w.im - n * k2) = T
  have hS1 : -n ≤ S := by grind
  have hS2 : S < n := by grind
  have hT1 : -n ≤ T := by grind
  have hT2 : T < n := by grind
  have p1 : 0 ≤ (n - S) * (n + S) := Int.mul_nonneg (by omega) (by omega)
  have p2 : 0 ≤ (n - T) * (n + T) := Int.mul_nonneg (by omega) (by omega)
  have p3 : 0 ≤ (n - S) * (n - T) := Int.mul_nonneg (by omega) (by omega)
  have p4 : 0 ≤ (n + S) * (n + T) := Int.mul_nonneg (by omega) (by omega)
  have p5 : 0 < n * n := Int.mul_pos hn hn
  have h4 : 4 * (N r * n) = S * S - S * T + T * T := by
    rw [hN, ← hS, ← hT]; simp only [N]; grind
  have hle : 4 * (N r * n) ≤ 3 * (n * n) := by grind
  apply Classical.byContradiction; intro hcon
  have : n * n ≤ N r * n := Int.mul_le_mul_of_nonneg_right (by omega) (by omega)
  omega

/-- Bezout: every pair has a common divisor that is a combination of the two. -/
theorem bezout_aux (k : Nat) : ∀ a b : Eis, (N b).toNat < k →
    ∃ g s t : Eis, g ∣ a ∧ g ∣ b ∧ g = s * a + t * b := by
  induction k with
  | zero => intro a b h; omega
  | succ k ih =>
    intro a b hk
    by_cases hb : b = 0
    · subst hb
      exact ⟨a, 1, 0, dvd_refl a, ⟨0, by grind⟩, by grind⟩
    · obtain ⟨q, r, hqr, hr⟩ := div_rem a b hb
      have h0 := N_nonneg r
      obtain ⟨g, s, t, hgb, hgr, hg⟩ := ih b r (by omega)
      refine ⟨g, t, s - t * q, ?_, hgb, by rw [hg, hqr]; grind⟩
      rw [hqr]; exact dvd_add (dvd_mul_left_of q hgb) hgr

/-- Bezout in ℤ[ω]. -/
theorem bezout (a b : Eis) : ∃ g s t : Eis, g ∣ a ∧ g ∣ b ∧ g = s * a + t * b :=
  bezout_aux _ a b (Nat.lt_succ_self _)

/-- Coprime: every common divisor is a unit. -/
def Cop (a b : Eis) : Prop := ∀ d, d ∣ a → d ∣ b → IsUnit d

theorem cop_bezout {a b : Eis} (h : Cop a b) : ∃ s t : Eis, s * a + t * b = 1 := by
  obtain ⟨g, s, t, hga, hgb, hg⟩ := bezout a b
  obtain ⟨v, hv⟩ := h g hga hgb
  exact ⟨v * s, v * t, by rw [← hv, hg]; grind⟩

/-- If a coprime product α·β is a cube, then α is a unit times a cube. -/
theorem cube_of_cop {α β γ : Eis} (hcop : Cop α β) (h : α * β = γ^3) :
    ∃ u δ : Eis, IsUnit u ∧ α = u * δ^3 := by
  by_cases hα : α = 0
  · exact ⟨1, 0, ⟨1, by grind⟩, by rw [hα]; grind⟩
  obtain ⟨g, s, t, hga, hgc, hg⟩ := bezout α γ
  have hg3 : g^3 = α * (s^3*α^2 + 3*s^2*α*t*γ + 3*s*t^2*γ^2 + t^3*β) := by
    rw [hg]; grind
  have hgb : Cop g β := fun d hd hb => hcop d (dvd_trans hd hga) hb
  obtain ⟨A, B, hAB⟩ := cop_bezout hgb
  obtain ⟨γ', hγ'⟩ := hgc
  have hα' : α = α * (A * g + B * β)^3 := by rw [hAB]; grind
  have he : α = g^3 * (α * A^3 + γ'^3 * (3*A^2*g^2*B + 3*A*g*B^2*β + B^3*β^2)) := by
    grind
  generalize (α * A^3 + γ'^3 * (3*A^2*g^2*B + 3*A*g*B^2*β + B^3*β^2)) = e at he
  generalize (s^3*α^2 + 3*s^2*α*t*γ + 3*s*t^2*γ^2 + t^3*β) = m at hg3
  have h0 : α * (1 - m * e) = 0 := by grind
  rcases mul_eq_zero h0 with h1 | h1
  · exact absurd h1 hα
  · exact ⟨e, g, ⟨m, by grind⟩, by grind⟩

end Eis

/-! ## 4. Euler's lemma: when p² + 3q² is a cube, p + q√−3 is a cube -/

theorem odd_norm {p q : Int} (hpar : (p + q) % 2 = 1) : (p^2 + 3*q^2) % 2 = 1 := by
  have hk : p + q = 2 * ((p + q) / 2) + 1 := by omega
  generalize (p + q) / 2 = k at hk
  have e : p^2 + 3*q^2 = 2 * (2*k^2 + 2*k + q^2 - p*q) + 1 := by
    have : p = 2*k + 1 - q := by omega
    subst this; grind
  rw [e]; omega

/-- The integer heart of the coprimality of p + q√−3 and p − q√−3. -/
theorem norm_one {n p q c : Int} (hn : 0 ≤ n) (h : p^2 + 3*q^2 = c^3)
    (hpq : Int.gcd p q = 1) (hpar : (p + q) % 2 = 1)
    (h1 : n ∣ p^2 + 3*q^2) (h2 : n ∣ 4 * p^2) (h3 : n ∣ 12 * q^2) : n = 1 := by
  have hodd := odd_norm hpar
  have hn2 : Int.gcd n 2 = 1 := cop_symm (cop_of_dvd_right (cop_two_of_odd hodd) h1)
  have hp2 : n ∣ p^2 := by
    have : 4 * p^2 = 2 * (2 * p^2) := by grind
    rw [this] at h2
    exact dvd_of_dvd_mul_cop hn2 (dvd_of_dvd_mul_cop hn2 h2)
  have hq2 : Int.gcd n (q^2) = 1 :=
    cop_of_dvd_left (Int.pow_gcd_pow_of_gcd_eq_one (k := 2) (l := 2) hpq) hp2
  have h12 : n ∣ 12 := by
    rw [Int.mul_comm] at h3; exact dvd_of_dvd_mul_cop hq2 h3
  have h3' : n ∣ 3 := by
    have : (12 : Int) = 2 * (2 * 3) := by decide
    rw [this] at h12
    exact dvd_of_dvd_mul_cop hn2 (dvd_of_dvd_mul_cop hn2 h12)
  have hle : n ≤ 3 := Int.le_of_dvd (by decide) h3'
  have : n = 0 ∨ n = 1 ∨ n = 2 ∨ n = 3 := by omega
  rcases this with e | e | e | e
  · subst e; have := Int.zero_dvd.1 h3'; omega
  · exact e
  · subst e; exact absurd h3' (by decide)
  · subst e
    exfalso
    have hp : (3 : Int) ∣ p := three_dvd_of_sq hp2
    obtain ⟨p', hp'⟩ := hp
    have hc : (3 : Int) ∣ c := three_dvd_of_cube ⟨p'^2 * 3 + q^2, by rw [← h, hp']; grind⟩
    obtain ⟨c', hc'⟩ := hc
    have hq : (3 : Int) ∣ q := by
      apply three_dvd_of_sq
      exact ⟨3 * c'^3 - p'^2, by
        have : (3*p')^2 + 3*q^2 = (3*c')^3 := by rw [← hp', ← hc']; exact h
        grind⟩
    have : (3 : Nat) ∣ Int.gcd p q := Int.dvd_gcd ⟨p', hp'⟩ hq
    rw [hpq] at this
    exact absurd this (by decide)

namespace Eis

theorem cube_eq (z : Eis) : z^3 = z * z * z := by grind

theorem cube_int (c : Int) : (⟨c, 0⟩ : Eis)^3 = ⟨c^3, 0⟩ := by
  rw [cube_eq]; ext <;> simp only [mul_re, mul_im] <;> grind

theorem omega_cube : (⟨0, 1⟩ : Eis)^3 = 1 := by
  rw [cube_eq]; ext <;> simp

/-- Every cube is the cube of something with even ω-part (something in ℤ[√−3]). -/
theorem even_cube_root (δ : Eis) : ∃ δ' : Eis, δ'^3 = δ^3 ∧ δ'.im % 2 = 0 := by
  obtain ⟨x, y⟩ := δ
  have hw : ∀ z : Eis, (⟨0, 1⟩ * z)^3 = z^3 := by
    intro z
    have : (⟨0, 1⟩ * z)^3 = (⟨0, 1⟩ : Eis)^3 * z^3 := by grind
    rw [this, omega_cube]; grind
  by_cases h1 : y % 2 = 0
  · exact ⟨⟨x, y⟩, rfl, h1⟩
  by_cases h2 : (x - y) % 2 = 0
  · refine ⟨⟨-y, x - y⟩, ?_, h2⟩
    have e : (⟨-y, x - y⟩ : Eis) = ⟨0, 1⟩ * ⟨x, y⟩ := by ext <;> simp
    rw [e, hw]
  · refine ⟨⟨y - x, -x⟩, ?_, by simp only; omega⟩
    have e : (⟨y - x, -x⟩ : Eis) = ⟨0, 1⟩ * (⟨0, 1⟩ * ⟨x, y⟩) := by ext <;> simp <;> omega
    rw [e, hw, hw]

theorem cube_sqrt3 (x k : Int) :
    (⟨x, 2 * k⟩ : Eis)^3 =
      ⟨((x-k)^3 - 9*(x-k)*k^2) + (3*(x-k)^2*k - 3*k^3), 2 * (3*(x-k)^2*k - 3*k^3)⟩ := by
  rw [cube_eq]; ext <;> simp only [mul_re, mul_im] <;> grind

end Eis

/-- Euler's lemma: if p² + 3q² is a cube, gcd(p, q) = 1 and p + q is odd, then
    p + q√−3 is a cube: p = a³ − 9ab², q = 3a²b − 3b³. -/
theorem key_cube {p q c : Int} (h : p^2 + 3*q^2 = c^3) (hpq : Int.gcd p q = 1)
    (hpar : (p + q) % 2 = 1) :
    ∃ a b : Int, p = a^3 - 9*a*b^2 ∧ q = 3*a^2*b - 3*b^3 := by
  let α : Eis := ⟨p + q, 2 * q⟩
  have hN : Eis.N α = p^2 + 3*q^2 := by simp only [Eis.N, α]; grind
  have hprod : α * Eis.conj α = (⟨c, 0⟩ : Eis)^3 := by
    rw [Eis.mul_conj, Eis.cube_int, hN, h]
  have hcop : Eis.Cop α (Eis.conj α) := by
    intro d hd hd'
    apply Eis.isUnit_of_N
    apply norm_one (Eis.N_nonneg d) h hpq hpar
    · rw [← hN]; exact Eis.N_dvd hd
    · have := Eis.N_dvd (Eis.dvd_add hd hd')
      have e : Eis.N (α + Eis.conj α) = 4 * p^2 := by simp only [Eis.N, Eis.conj, α, Eis.add_re, Eis.add_im]; grind
      rwa [e] at this
    · have := Eis.N_dvd (Eis.dvd_sub hd hd')
      have e : Eis.N (α - Eis.conj α) = 12 * q^2 := by simp only [Eis.N, Eis.conj, α, Eis.sub_re, Eis.sub_im]; grind
      rwa [e] at this
  obtain ⟨u, δ, hu, hα⟩ := Eis.cube_of_cop hcop hprod
  obtain ⟨δ', hδ', hev⟩ := Eis.even_cube_root δ
  obtain ⟨x, y⟩ := δ'
  simp only at hev
  have hy : y = 2 * (y / 2) := by omega
  generalize y / 2 = k at hy
  subst hy
  rw [← hδ', Eis.cube_sqrt3] at hα
  have hpre : p + q = (α).re := rfl
  have hpim : 2 * q = (α).im := rfl
  rw [hα] at hpre hpim
  generalize hP : (x-k)^3 - 9*(x-k)*k^2 = P at hpre hpim
  generalize hQ : 3*(x-k)^2*k - 3*k^3 = Q at hpre hpim
  rcases Eis.units_list (Eis.N_of_isUnit hu) with e|e|e|e|e|e <;> subst e <;>
    simp at hpre hpim
  · exact ⟨x - k, k, by omega, by omega⟩
  · exact ⟨-(x - k), -k, by grind, by grind⟩
  all_goals omega

/-! ## 5. Euler's descent for the second case -/

theorem cube_mod9_all (x : Int) : x^3 % 9 = 0 ∨ x^3 % 9 = 1 ∨ x^3 % 9 = 8 := by
  rw [pow_emod_of]
  have : x % 9 = 0 ∨ x % 9 = 1 ∨ x % 9 = 2 ∨ x % 9 = 3 ∨ x % 9 = 4 ∨ x % 9 = 5 ∨
    x % 9 = 6 ∨ x % 9 = 7 ∨ x % 9 = 8 := by omega
  rcases this with e|e|e|e|e|e|e|e|e <;> rw [e] <;> decide

/-- 2 and −2 are not cubes (cubes mod 9 are 0, 1, 8). -/
theorem cube_ne_two (x : Int) : x^3 ≠ 2 ∧ x^3 ≠ -2 := by
  have := cube_mod9_all x
  constructor <;> intro e <;> rw [e] at this <;> revert this <;> decide

theorem cube_mod2 (x : Int) : x^3 % 2 = x % 2 := by
  rw [pow_emod_of]
  have : x % 2 = 0 ∨ x % 2 = 1 := by omega
  rcases this with e|e <;> rw [e] <;> decide

theorem eq_zero_of_cube {w : Int} (h : w^3 = 0) : w = 0 := by
  have : w * (w * w) = 0 := by rw [← h]; grind
  rcases Int.mul_eq_zero.1 this with e | e
  · exact e
  · rcases Int.mul_eq_zero.1 e with e | e <;> exact e

theorem cube_ne_zero {w : Int} (h : w ≠ 0) : w^3 ≠ 0 := fun e => h (eq_zero_of_cube e)

/-- In a solution of x³ + y³ + z³ = 0, gcd(x, y) = 1 forces gcd(y, z) = 1. -/
theorem cop_of_cubes {x y z : Int} (heq : x^3 + y^3 + z^3 = 0) (hxy : Int.gcd x y = 1) :
    Int.gcd y z = 1 := by
  have h1 : ((Int.gcd y z : Int))^3 ∣ x^3 := by
    have a := Int.pow_dvd_pow_iff (n := 3) (by decide) |>.2 (Int.gcd_dvd_left y z)
    have b := Int.pow_dvd_pow_iff (n := 3) (by decide) |>.2 (Int.gcd_dvd_right y z)
    have : x^3 = -(y^3 + z^3) := by grind
    rw [this]; exact Int.dvd_neg.2 (Int.dvd_add a b)
  have h2 := (Int.pow_dvd_pow_iff (by decide)).1 h1
  have : Int.gcd y z ∣ Int.gcd x y := Int.dvd_gcd h2 (Int.gcd_dvd_left y z)
  rw [hxy] at this
  exact Nat.dvd_one.1 this

/-- A primitive nontrivial solution of x³ + y³ + z³ = 0. -/
def Sol (x y z : Int) : Prop := x ≠ 0 ∧ y ≠ 0 ∧ z ≠ 0 ∧ x^3 + y^3 + z^3 = 0 ∧ Int.gcd x y = 1

theorem sol_cop {x y z : Int} (hs : Sol x y z) : Int.gcd y z = 1 :=
  cop_of_cubes hs.2.2.2.1 hs.2.2.2.2

theorem sol_swap {x y z : Int} (hs : Sol x y z) : Sol y x z := by
  obtain ⟨a, b, c, heq, hxy⟩ := hs
  exact ⟨b, a, c, by grind, cop_symm hxy⟩

theorem sol_rot {x y z : Int} (hs : Sol x y z) : Sol y z x := by
  have h := sol_cop hs
  obtain ⟨a, b, c, heq, hxy⟩ := hs
  exact ⟨b, c, a, by grind, h⟩

theorem odd_of_cop_even {x z : Int} (h : Int.gcd x z = 1) (hz : z % 2 = 0) : x % 2 = 1 := by
  apply Classical.byContradiction; intro hx
  have : (2 : Nat) ∣ Int.gcd x z :=
    Int.dvd_gcd (Int.dvd_of_emod_eq_zero (by omega)) (Int.dvd_of_emod_eq_zero hz)
  rw [h] at this; exact absurd this (by decide)

/-- Some entry of a primitive solution is even. -/
theorem sol_even {x y z : Int} (hs : Sol x y z) : x % 2 = 0 ∨ y % 2 = 0 ∨ z % 2 = 0 := by
  obtain ⟨-, -, -, heq, -⟩ := hs
  have a := cube_mod2 x
  have b := cube_mod2 y
  have c := cube_mod2 z
  generalize x^3 = X at *
  generalize y^3 = Y at *
  generalize z^3 = Z at *
  omega

theorem measure_lt {x y z : Int} {N1 : Nat} (hx : x ≠ 0) (hy : y ≠ 0)
    (h : N1^3 < (z^3).natAbs) : N1 < (x * y * z).natAbs := by
  apply Classical.byContradiction; intro hc
  have h1 : (x * y * z).natAbs ^ 3 ≤ N1^3 := Nat.pow_le_pow_left (by omega) 3
  have h2 : (x * y * z).natAbs ^ 3 = (x^3).natAbs * (y^3).natAbs * (z^3).natAbs := by
    rw [← Int.natAbs_pow, ← Int.natAbs_mul, ← Int.natAbs_mul]; congr 1; grind
  have hx3 : 0 < (x^3).natAbs := Int.natAbs_pos.2 (cube_ne_zero hx)
  have hy3 : 0 < (y^3).natAbs := Int.natAbs_pos.2 (cube_ne_zero hy)
  have h3 : (z^3).natAbs ≤ (x^3).natAbs * (y^3).natAbs * (z^3).natAbs :=
    Nat.le_mul_of_pos_left _ (Nat.mul_pos hx3 hy3)
  omega

/-- (p + q) odd forces (a + b) odd in Euler's parametrization. -/
theorem par_ab {a b : Int} (h : (a^3 - 9*a*b^2 + (3*a^2*b - 3*b^3)) % 2 = 1) : (a + b) % 2 = 1 := by
  apply Classical.byContradiction; intro hc
  have hk : a = b + 2 * ((a - b) / 2) := by omega
  generalize (a - b) / 2 = k at hk
  have e : a^3 - 9*a*b^2 + (3*a^2*b - 3*b^3) = 2 * (-4*b^3 + 12*b*k^2 + 4*k^3) := by
    subst hk; grind
  rw [e] at h; omega

theorem cop_of_param {a b p q : Int} (hpq : Int.gcd p q = 1) (hp : p = a^3 - 9*a*b^2)
    (hq : q = 3*a^2*b - 3*b^3) : Int.gcd a b = 1 := by
  rw [Int.gcd_eq_one_iff] at hpq ⊢
  intro c ha hb
  apply hpq c
  · have : p = a * (a^2 - 9*b^2) := by rw [hp]; grind
    rw [this]; exact Int.dvd_mul_of_dvd_left ha
  · have : q = b * (3*a^2 - 3*b^2) := by rw [hq]; grind
    rw [this]; exact Int.dvd_mul_of_dvd_left hb

/-- The descent step (Euler): from a primitive solution with z even, build a smaller one. -/
theorem step {x y z : Int} (hs : Sol x y z) (hz : z % 2 = 0) :
    ∃ x' y' z', Sol x' y' z' ∧ (x' * y' * z').natAbs < (x * y * z).natAbs := by
  have hyz := sol_cop hs
  have hxz := sol_cop (sol_rot hs)
  obtain ⟨hx0, hy0, hz0, heq, hxy⟩ := hs
  have hxz' : Int.gcd x z = 1 := cop_symm hxz
  have hxo := odd_of_cop_even hxz' hz
  have hyo := odd_of_cop_even hyz hz
  have hpx : x + y = 2 * ((x + y) / 2) := by omega
  have hqx : x - y = 2 * ((x - y) / 2) := by omega
  generalize (x + y) / 2 = p at hpx
  generalize (x - y) / 2 = q at hqx
  have hx : x = p + q := by omega
  have hy : y = p - q := by omega
  have hsum : 2 * p * (p^2 + 3*q^2) = (-z)^3 := by
    have : x^3 + y^3 = 2 * p * (p^2 + 3*q^2) := by rw [hx, hy]; grind
    grind
  have hz3 : z^3 = -(2 * p * (p^2 + 3*q^2)) := by grind
  have hpq : Int.gcd p q = 1 := by
    rw [Int.gcd_eq_one_iff] at hxy ⊢
    intro c hc1 hc2
    exact hxy c (by rw [hx]; exact Int.dvd_add hc1 hc2) (by rw [hy]; exact Int.dvd_sub hc1 hc2)
  have hpar : (p + q) % 2 = 1 := by omega
  have hp0 : p ≠ 0 := by
    intro e; subst e
    have : (-z)^3 = 0 := by rw [← hsum]; grind
    exact hz0 (by have := eq_zero_of_cube this; omega)
  have hodd := odd_norm hpar
  have hN1 : 1 ≤ p^2 + 3*q^2 := by
    have := Eis.sq_pos' hp0
    have := Eis.sq_nonneg' q
    grind
  by_cases h3 : (3 : Int) ∣ p
  · -- Case B: 3 ∣ p
    obtain ⟨r, hr⟩ := h3
    have hq3 : ¬ (3 : Int) ∣ q := by
      intro hq
      have : (3 : Nat) ∣ Int.gcd p q := Int.dvd_gcd ⟨r, hr⟩ hq
      rw [hpq] at this; exact absurd this (by decide)
    have hsum' : 18 * r * (q^2 + 3*r^2) = (-z)^3 := by rw [← hsum, hr]; grind
    have hqr : Int.gcd q r = 1 := by
      rw [Int.gcd_eq_one_iff] at hpq ⊢
      intro c hc1 hc2
      exact hpq c (by rw [hr]; exact Int.dvd_mul_of_dvd_right hc2) hc1
    have hpar' : (q + r) % 2 = 1 := by omega
    have hr0 : r ≠ 0 := by intro e; apply hp0; rw [hr, e]; rfl
    have hodd' := odd_norm hpar'
    have hcop : Int.gcd (18 * r) (q^2 + 3*r^2) = 1 := by
      have c2 := cop_two_of_odd hodd'
      have c3 : Int.gcd 3 (q^2 + 3*r^2) = 1 := by
        apply cop_three
        have : ¬ (3 : Int) ∣ q^2 := fun h => hq3 (three_dvd_of_sq h)
        have : q^2 % 3 ≠ 0 := fun e => this (Int.dvd_of_emod_eq_zero e)
        generalize q^2 = Q at *
        generalize r^2 = R at *
        omega
      have cr : Int.gcd r (q^2 + 3*r^2) = 1 := by
        have : q^2 + 3*r^2 = q^2 + r * (3 * r) := by grind
        rw [this, Int.gcd_add_mul_left_right]
        exact Int.gcd_pow_right_of_gcd_eq_one (cop_symm hqr)
      have : (18 : Int) * r = 2 * 3 * 3 * r := by grind
      rw [this]
      exact cop_mul_left (cop_mul_left (cop_mul_left c2 c3) c3) cr
    obtain ⟨c, hc⟩ := cube_of_mul_cop' hsum' hcop
    obtain ⟨w, hw⟩ := cube_of_mul_cop hsum' hcop
    obtain ⟨a, b, hqa, hrb⟩ := key_cube hc hqr hpar'
    have hab := cop_of_param hqr hqa hrb
    have habpar : (a + b) % 2 = 1 := par_ab (by rw [← hqa, ← hrb]; exact hpar')
    have hb0 : b ≠ 0 := by intro e; apply hr0; rw [hrb, e]; grind
    have h27 : w^3 = 27 * (2 * b * (a - b) * (a + b)) := by rw [← hw, hrb]; grind
    obtain ⟨w', hw'⟩ := three_dvd_of_cube ⟨9 * (2 * b * (a - b) * (a + b)), by rw [h27]; grind⟩
    have hm : 2 * b * (a - b) * (a + b) = w'^3 := by
      rw [hw'] at h27; grind
    have o1 : (a - b) % 2 = 1 := by omega
    have o2 : (a + b) % 2 = 1 := by omega
    have g1 : Int.gcd b (a - b) = 1 := by rw [Int.gcd_sub_self_right]; exact cop_symm hab
    have g2 : Int.gcd b (a + b) = 1 := by rw [Int.gcd_add_self_right]; exact cop_symm hab
    have c1 : Int.gcd (2 * b) (a - b) = 1 := cop_mul_left (cop_two_of_odd o1) g1
    have c2 : Int.gcd (2 * b) (a + b) = 1 := cop_mul_left (cop_two_of_odd o2) g2
    have c3 : Int.gcd (a - b) (a + b) = 1 := by
      have : a + b = (a - b) + 2 * b := by omega
      rw [this, Int.gcd_self_add_right]; exact cop_symm c1
    obtain ⟨⟨G, hG⟩, ⟨A, hA⟩, ⟨B, hB⟩⟩ := cube_of_mul3 hm c1 c2 c3
    refine ⟨B, -A, -G, ⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
    · intro e; rw [e] at hB; omega
    · intro e; have : A = 0 := by omega
      rw [this] at hA; omega
    · intro e; have : G = 0 := by omega
      rw [this] at hG; exact hb0 (by omega)
    · grind
    · have hBA : Int.gcd B A = 1 := cop_of_pow (n := 3) (by decide) (by rw [← hA, ← hB]; exact cop_symm c3)
      rw [Int.gcd_eq_one_iff] at hBA ⊢
      intro d h1 h2; exact hBA d h1 (Int.dvd_neg.1 h2)
    · apply measure_lt hx0 hy0
      have e1 : (B * -A * -G).natAbs ^ 3 = (2 * b * (a - b) * (a + b)).natAbs := by
        rw [← Int.natAbs_pow]; congr 1; rw [hG, hA, hB]; grind
      have e2 : (z^3).natAbs = (18 * r).natAbs * (q^2 + 3*r^2).natAbs := by
        rw [← Int.natAbs_mul]
        have : z^3 = -(18 * r * (q^2 + 3*r^2)) := by grind
        rw [this, Int.natAbs_neg]
      have e3 : (18 * r).natAbs = 27 * (2 * b * (a - b) * (a + b)).natAbs := by
        have : 18 * r = 27 * (2 * b * (a - b) * (a + b)) := by rw [hrb]; grind
        rw [this, Int.natAbs_mul]; rfl
      have hm0 : 0 < (2 * b * (a - b) * (a + b)).natAbs := by
        apply Int.natAbs_pos.2
        intro e
        rcases Int.mul_eq_zero.1 e with e | e
        · rcases Int.mul_eq_zero.1 e with e | e
          · rcases Int.mul_eq_zero.1 e with e | e
            · omega
            · exact hb0 e
          · omega
        · omega
      have hQ : 1 ≤ (q^2 + 3*r^2).natAbs := by
        have := Eis.sq_nonneg' q
        have := Eis.sq_pos' hr0
        have : 1 ≤ q^2 + 3*r^2 := by grind
        omega
      rw [e1, e2, e3]
      have := Nat.le_mul_of_pos_right (27 * (2 * b * (a - b) * (a + b)).natAbs) hQ
      omega
  · -- Case A: 3 ∤ p
    have hcop : Int.gcd (2 * p) (p^2 + 3*q^2) = 1 := by
      have c2 := cop_two_of_odd hodd
      have cp : Int.gcd p (p^2 + 3*q^2) = 1 := by
        have : p^2 + 3*q^2 = 3*q^2 + p * p := by grind
        rw [this, Int.gcd_add_mul_left_right]
        apply cop_mul_right
        · exact cop_symm (cop_three (fun e => h3 (Int.dvd_of_emod_eq_zero e)))
        · exact Int.gcd_pow_right_of_gcd_eq_one hpq
      exact cop_mul_left c2 cp
    obtain ⟨c, hc⟩ := cube_of_mul_cop' hsum hcop
    obtain ⟨w, hw⟩ := cube_of_mul_cop hsum hcop
    obtain ⟨a, b, hpa, hqb⟩ := key_cube hc hpq hpar
    have hab := cop_of_param hpq hpa hqb
    have habpar : (a + b) % 2 = 1 := par_ab (by rw [← hpa, ← hqb]; exact hpar)
    have ha3 : ¬ (3 : Int) ∣ a := by
      intro h; apply h3
      have : p = a * (a^2 - 9*b^2) := by rw [hpa]; grind
      rw [this]; exact Int.dvd_mul_of_dvd_left h
    have ha0 : a ≠ 0 := by intro e; apply hp0; rw [hpa, e]; grind
    have hm : 2 * a * (a - 3 * b) * (a + 3 * b) = w^3 := by rw [← hw, hpa]; grind
    have o1 : (a - 3 * b) % 2 = 1 := by omega
    have o2 : (a + 3 * b) % 2 = 1 := by omega
    have ga : Int.gcd a (3 * b) = 1 :=
      cop_mul_right (cop_symm (cop_three (fun e => ha3 (Int.dvd_of_emod_eq_zero e)))) hab
    have g1 : Int.gcd a (a - 3 * b) = 1 := by rw [Int.gcd_self_sub_right]; exact ga
    have g2 : Int.gcd a (a + 3 * b) = 1 := by rw [Int.gcd_self_add_right]; exact ga
    have c1 : Int.gcd (2 * a) (a - 3 * b) = 1 := cop_mul_left (cop_two_of_odd o1) g1
    have c2 : Int.gcd (2 * a) (a + 3 * b) = 1 := cop_mul_left (cop_two_of_odd o2) g2
    have c3 : Int.gcd (a - 3 * b) (a + 3 * b) = 1 := by
      have : a + 3 * b = (a - 3 * b) + 2 * (3 * b) := by omega
      rw [this, Int.gcd_self_add_right]
      apply cop_mul_right (cop_symm (cop_two_of_odd o1))
      rw [Int.gcd_sub_self_left]; exact ga
    obtain ⟨⟨G, hG⟩, ⟨A, hA⟩, ⟨B, hB⟩⟩ := cube_of_mul3 hm c1 c2 c3
    refine ⟨A, B, -G, ⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
    · intro e; rw [e] at hA; omega
    · intro e; rw [e] at hB; omega
    · intro e; have : G = 0 := by omega
      rw [this] at hG; exact ha0 (by omega)
    · grind
    · exact cop_of_pow (n := 3) (by decide) (by rw [← hA, ← hB]; exact c3)
    · apply measure_lt hx0 hy0
      have e1 : (A * B * -G).natAbs ^ 3 = (2 * p).natAbs := by
        rw [← Int.natAbs_pow]
        have : (A * B * -G)^3 = -(2 * p) := by
          have : (A * B * -G)^3 = -(A^3 * B^3 * G^3) := by grind
          rw [this, ← hA, ← hB, ← hG, hpa]; grind
        rw [this, Int.natAbs_neg]
      have e2 : (z^3).natAbs = (2 * p).natAbs * (p^2 + 3*q^2).natAbs := by
        rw [← Int.natAbs_mul, hz3, Int.natAbs_neg]
      have hp2 : 0 < (2 * p).natAbs := Int.natAbs_pos.2 (by omega)
      rw [e1, e2]
      by_cases hbig : 2 ≤ p^2 + 3*q^2
      · have : 2 ≤ (p^2 + 3*q^2).natAbs := by omega
        have := Nat.mul_le_mul_left (2 * p).natAbs this
        omega
      · exfalso
        have hq0 : q = 0 := by
          apply Classical.byContradiction; intro hq
          have := Eis.sq_pos' hq
          have := Eis.sq_pos' hp0
          grind
        subst hq0
        have hp1 : p * p ≤ 1 := by grind
        rcases Eis.small_of_sq_le_one hp1 with e | e | e
        · subst e; exact (cube_ne_two z).1 (by grind)
        · exact hp0 e
        · subst e; exact (cube_ne_two z).2 (by grind)

/-- No primitive solution exists (infinite descent on |xyz|). -/
theorem no_sol : ∀ n : Nat, ∀ x y z : Int, (x * y * z).natAbs = n → ¬ Sol x y z := by
  intro n
  induction n using Nat.strongRecOn with
  | _ n ih =>
    intro x y z hn hs
    have key : ∀ x y z : Int, Sol x y z → z % 2 = 0 → (x * y * z).natAbs = n → False := by
      intro x y z hs hz hn
      obtain ⟨x', y', z', hs', hlt⟩ := step hs hz
      exact ih _ (hn ▸ hlt) x' y' z' rfl hs'
    rcases sol_even hs with h | h | h
    · exact key y z x (sol_rot hs) h (by rw [← hn]; congr 1; grind)
    · exact key z x y (sol_rot (sol_rot hs)) h (by rw [← hn]; congr 1; grind)
    · exact key x y z hs h hn

/-! ## 6. The theorem, and the supporting lemmas -/

/-- Reduction to coprime solutions: dividing by gcd(x, y) gives a solution with gcd 1. -/
theorem primitive_of_sol {x y z : Int} (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (h : x^3 + y^3 = z^3) :
    ∃ x' y' z' : Int, x' ≠ 0 ∧ y' ≠ 0 ∧ z' ≠ 0 ∧ x'^3 + y'^3 = z'^3 ∧ Int.gcd x' y' = 1 := by
  have hg : 0 < Int.gcd x y := Int.gcd_pos_of_ne_zero_left y hx
  generalize hgd : Int.gcd x y = g at hg
  have hdx : (g : Int) ∣ x := hgd ▸ Int.gcd_dvd_left x y
  have hdy : (g : Int) ∣ y := hgd ▸ Int.gcd_dvd_right x y
  obtain ⟨x', hx'⟩ := hdx
  obtain ⟨y', hy'⟩ := hdy
  have hg0 : (g : Int) ≠ 0 := by omega
  have hcop : Int.gcd x' y' = 1 := by
    have := Int.gcd_ediv_gcd_ediv_gcd (i := x) (j := y) (by omega)
    rw [hgd] at this
    rwa [hx', hy', Int.mul_ediv_cancel_left _ hg0, Int.mul_ediv_cancel_left _ hg0] at this
  have hz3 : z^3 = (g : Int)^3 * (x'^3 + y'^3) := by rw [← h, hx', hy']; grind
  have hdz : (g : Int) ∣ z :=
    (Int.pow_dvd_pow_iff (n := 3) (by decide)).1 ⟨x'^3 + y'^3, hz3⟩
  obtain ⟨z', hz'⟩ := hdz
  refine ⟨x', y', z', ?_, ?_, ?_, ?_, hcop⟩
  · intro e; apply hx; rw [hx', e]; grind
  · intro e; apply hy; rw [hy', e]; grind
  · intro e; apply hz; rw [hz', e]; grind
  · have hg3 : (g : Int)^3 ≠ 0 := cube_ne_zero hg0
    have : (g : Int)^3 * (x'^3 + y'^3) = (g : Int)^3 * z'^3 := by rw [← hz3, hz']; grind
    exact Int.eq_of_mul_eq_mul_left hg3 this

/-- **Fermat's Last Theorem for exponent 3** (integers): no nonzero x, y, z with x³ + y³ = z³. -/
theorem flt3 (x y z : Int) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0) : x^3 + y^3 ≠ z^3 := by
  intro h
  obtain ⟨a, b, c, ha, hb, hc, h', hab⟩ := primitive_of_sol hx hy hz h
  have hs : Sol a b (-c) := ⟨ha, hb, by omega, by grind, hab⟩
  exact no_sol _ a b (-c) rfl hs

/-- **Fermat's Last Theorem for exponent 3** (natural numbers): x³ + y³ ≠ z³ when x, y ≠ 0. -/
theorem flt3_nat (x y z : Nat) (hx : x ≠ 0) (hy : y ≠ 0) : x^3 + y^3 ≠ z^3 := by
  intro h
  have h' : (x : Int)^3 + (y : Int)^3 = (z : Int)^3 := by exact_mod_cast h
  have hz : (z : Int) ≠ 0 := by
    intro e
    rw [e] at h'
    have h1 : (0 : Int) < (x : Int)^3 := Int.pow_pos (by omega)
    have h2 : (0 : Int) < (y : Int)^3 := Int.pow_pos (by omega)
    grind
  exact flt3 x y z (by omega) (by omega) hz h'

/-- In the mult-add picture: for integers x, y, z ≠ 0, (y/x)³ is never the mult-add cube of z/x − 1. -/
theorem flt3_circ (x y z : Int) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0) :
    ((y : Rat) / x)^3 ≠ mpow ((z : Rat) / x - 1) 3 := by
  intro h
  have hx' : (x : Rat) ≠ 0 := by exact_mod_cast hx
  have := (flt3_circ_iff (x : Rat) y z hx').2 h
  exact flt3 x y z hx hy hz (by exact_mod_cast this)

/-! ## Supporting lemmas -/

/-- The sum of two cubes factors: x³ + y³ = (x + y)(x² − xy + y²). -/
theorem sum_cubes (x y : Int) : x^3 + y^3 = (x + y) * (x^2 - x*y + y^2) := by grind

/-- The second factor is (x + y)² − 3xy. -/
theorem second_factor (x y : Int) : x^2 - x*y + y^2 = (x + y)^2 - 3*(x*y) := by grind

/-- For coprime x, y the two factors of x³ + y³ share only 1 or 3. -/
theorem gcd_factors {x y : Int} (h : Int.gcd x y = 1) :
    Int.gcd (x + y) (x^2 - x*y + y^2) = 1 ∨ Int.gcd (x + y) (x^2 - x*y + y^2) = 3 := by
  have e : x^2 - x*y + y^2 = (-3) * (x*y) + (x + y) * (x + y) := by grind
  have c1 : Int.gcd (x + y) x = 1 := by
    have := Int.gcd_self_add_right x y
    exact cop_symm (this.trans h)
  have c2 : Int.gcd (x + y) y = 1 := by
    have := Int.gcd_add_self_right y x
    exact cop_symm (this.trans (cop_symm h))
  have c3 : Int.gcd (x + y) (x * y) = 1 := cop_mul_right c1 c2
  rw [e, Int.gcd_add_mul_left_right, Int.gcd_mul_left_right_of_gcd_eq_one c3]
  have hd : Int.gcd (x + y) (-3) ∣ 3 := Int.gcd_dvd_natAbs_right _ _
  have hle := Nat.le_of_dvd (by decide) hd
  have : Int.gcd (x + y) (-3) ≠ 0 := Int.gcd_ne_zero_right (by decide)
  have : Int.gcd (x + y) (-3) ≠ 2 := fun e => by rw [e] at hd; exact absurd hd (by decide)
  omega


theorem three_dvd_mul {a b : Int} (h : (3 : Int) ∣ a * b) : (3 : Int) ∣ a ∨ (3 : Int) ∣ b := by
  have h' := Int.emod_eq_zero_of_dvd h
  rw [Int.mul_emod] at h'
  have ha : a % 3 = 0 ∨ a % 3 = 1 ∨ a % 3 = 2 := by omega
  have hb : b % 3 = 0 ∨ b % 3 = 1 ∨ b % 3 = 2 := by omega
  rcases ha with ea|ea|ea
  · exact Or.inl (Int.dvd_of_emod_eq_zero ea)
  all_goals rcases hb with eb|eb|eb
  all_goals first
    | exact Or.inr (Int.dvd_of_emod_eq_zero eb)
    | (rw [ea, eb] at h'; exact absurd h' (by decide))

/-- Second case: in a coprime solution with 3 ∣ xyz, exactly one of x, y, z is divisible by 3. -/
theorem exactly_one_three {x y z : Int} (h : x^3 + y^3 = z^3) (hxy : Int.gcd x y = 1)
    (h3 : (3 : Int) ∣ x * y * z) :
    ((3 : Int) ∣ x ∧ ¬ (3 : Int) ∣ y ∧ ¬ (3 : Int) ∣ z) ∨
    (¬ (3 : Int) ∣ x ∧ (3 : Int) ∣ y ∧ ¬ (3 : Int) ∣ z) ∨
    (¬ (3 : Int) ∣ x ∧ ¬ (3 : Int) ∣ y ∧ (3 : Int) ∣ z) := by
  have heq : x^3 + y^3 + (-z)^3 = 0 := by grind
  have hyz : Int.gcd y (-z) = 1 := cop_of_cubes heq hxy
  have hxz : Int.gcd x (-z) = 1 := by
    have heq' : y^3 + x^3 + (-z)^3 = 0 := by grind
    exact cop_of_cubes heq' (cop_symm hxy)
  have no2 : ∀ a b : Int, Int.gcd a b = 1 → (3 : Int) ∣ a → (3 : Int) ∣ b → False := by
    intro a b hab ha hb
    have : (3 : Nat) ∣ Int.gcd a b := Int.dvd_gcd ha hb
    rw [hab] at this; exact absurd this (by decide)
  have zn : (3 : Int) ∣ z → (3 : Int) ∣ -z := fun h => Int.dvd_neg.2 h
  rcases three_dvd_mul h3 with h1 | h1
  · rcases three_dvd_mul h1 with h2 | h2
    · exact Or.inl ⟨h2, fun hy => no2 x y hxy h2 hy, fun hz => no2 x (-z) hxz h2 (zn hz)⟩
    · exact Or.inr (Or.inl ⟨fun hx => no2 x y hxy hx h2, h2, fun hz => no2 y (-z) hyz h2 (zn hz)⟩)
  · exact Or.inr (Or.inr ⟨fun hx => no2 x (-z) hxz hx (zn h1),
      fun hy => no2 y (-z) hyz hy (zn h1), h1⟩)

/-- Mult-add cubes and plain cubes agree mod 3: (1+s)³ − 1 ≡ s³ ≡ s. -/
theorem mpow_three_mod3 (s : Int) : (3 : Int) ∣ mpow s 3 - s := by
  rw [mpow_three]
  have := cube_mod3 s
  have e : s^3 + 3*(s + s^2) - s = (s^3 - s) + 3*(s + s^2) := by grind
  rw [e]; exact Int.dvd_add this (Int.dvd_mul_right 3 _)

end Fermat3

/-! ## Axiom check: only propext, Classical.choice, Quot.sound may appear. -/

