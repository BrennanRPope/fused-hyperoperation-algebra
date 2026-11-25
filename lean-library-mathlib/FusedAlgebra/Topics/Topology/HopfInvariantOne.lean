module

public import FusedAlgebra.Basic

/-!
# Only 1, 2 and 4: the number step of the Hopf invariant one theorem

The theorem (Adams 1960; the short proof is Adams and Atiyah 1966): a map `S^(4n-1) → S^(2n)` of
odd Hopf invariant exists only for `n = 1, 2, 4`. This file contains the arithmetic of the
Adams–Atiyah proof, not the topology.

## Part 1: numbers

* `MultAdd.HopfOne.level_three_pow`: the exact power of `2` in `3 ^ n - 1`. It is `2` when `n` is
  odd and `2 ^ (k + 2)` when `n = 2 ^ k * odd` with `k ≥ 1`.
* `MultAdd.HopfOne.only_one_two_four`: for `n ≥ 1`, `2 ^ n ∣ 3 ^ n - 1` forces `n = 1, 2, 4`.

## Part 2: commuting Adams operations

* `MultAdd.HopfOne.Elt`, `MultAdd.HopfOne.psi`: the two-generator K-theory of the glued space
  and the Adams operations `ψ^k (a) = k ^ n * a + μ * b`, `ψ^k (b) = k ^ (2 * n) * b` on it.
* `MultAdd.HopfOne.hopf_invariant_one`: if `ψ²` and `ψ³` commute, `ψ² (a) ≡ a²` mod `2` and the
  Hopf invariant `h` is odd, then `n = 1, 2, 4`.

## Part 3: the mult-add cube

* `MultAdd.HopfOne.Trunc`: the polynomials `α * x + β * x ^ 2` in `ℤ[x] / (x ^ 3)`, a
  commutative ring without unit. Mult-add on it is the library's `✳` and the `k`-fold mult-add
  power of `x` is the library's `MultAdd.mpow Trunc.x k`.
* `MultAdd.HopfOne.mpow_formula`: `mpow x k = k * x + k.choose 2 * x ^ 2`.
* `MultAdd.HopfOne.psi_compose`: on `ℂP²` the Adams operations are mult-add powers, and they
  compose, `ψ^j ψ^k = ψ^(j * k)`.
* `MultAdd.HopfOne.cp2_case`: `ℂP²` (the case `n = 1`, `h = 1`) meets every hypothesis of Part 2.

## What is taken as given

The topology is not formalised: that the glued space's K-theory has this two-generator shape,
that the Adams operations exist on it, are additive and compose, that `ψ^k` acts on a `2m`-cell
by `k ^ m`, and that `ψ^p (x) ≡ x ^ p` mod `p`. These standard K-theory facts appear here only as
the hypotheses of `hopf_invariant_one`.

This is the mathlib form of `FusedAlgebra.Topics.Topology.HopfInvariantOne`. Part 1 is the
lifting-the-exponent lemma at the prime `2`, which is in mathlib (`Nat.two_pow_sub_pow`); the
hand-built predicate `OneMod` of the core-Lean file is replaced by `multiplicity 2`.
-/

@[expose] public section

namespace MultAdd.HopfOne

open MultAdd

/-! ### Part 1: the exact power of 2 in `3 ^ n - 1` -/

/-- For odd `n`, `3 ^ n - 1` is `2` times an odd number. -/
theorem multiplicity_two_three_pow_sub_one_of_odd {n : ℕ} (hn : Odd n) :
    multiplicity 2 (3 ^ n - 1) = 1 := by
  obtain ⟨t, rfl⟩ := hn
  have h9 : 9 ^ t % 4 = 1 := by rw [Nat.pow_mod]; simp
  have h3 : 3 ^ (2 * t + 1) % 4 = 3 := by
    rw [pow_succ, pow_mul, Nat.mul_mod]
    norm_num [h9]
  have h1 : 1 ≤ 3 ^ (2 * t + 1) := Nat.one_le_pow _ _ (by decide)
  generalize 3 ^ (2 * t + 1) = N at h3 h1
  refine multiplicity_eq_of_dvd_of_not_dvd ?_ ?_
  · rw [pow_one]; omega
  · norm_num; omega

/-- For even `n ≠ 0`, the power of `2` in `3 ^ n - 1` is two more than the power of `2` in `n`.
This is the lifting-the-exponent lemma `Nat.two_pow_sub_pow`. -/
theorem multiplicity_two_three_pow_sub_one_of_even {n : ℕ} (hn : Even n) (h0 : n ≠ 0) :
    multiplicity 2 (3 ^ n - 1) = multiplicity 2 n + 2 := by
  have h := Nat.two_pow_sub_pow (x := 3) (y := 1) (by decide) (by decide) hn
  have h1 : 1 < 3 ^ n := Nat.one_lt_pow h0 (by decide)
  have hf1 : FiniteMultiplicity 2 (3 ^ n - 1) := Nat.finiteMultiplicity_iff.2 ⟨by decide, by omega⟩
  have hf2 : FiniteMultiplicity 2 n := Nat.finiteMultiplicity_iff.2 ⟨by decide, by omega⟩
  have h4 : emultiplicity 2 (3 + 1) = ((2 : ℕ) : ℕ∞) :=
    emultiplicity_eq_of_dvd_of_not_dvd (by decide) (by decide)
  have h2 : emultiplicity 2 (3 - 1) = ((1 : ℕ) : ℕ∞) :=
    emultiplicity_eq_of_dvd_of_not_dvd (by decide) (by decide)
  rw [one_pow, h4, h2, hf1.emultiplicity_eq_multiplicity, hf2.emultiplicity_eq_multiplicity] at h
  have h' : multiplicity 2 (3 ^ n - 1) + 1 = 2 + 1 + multiplicity 2 n := by exact_mod_cast h
  omega

/-- The power of `2` in `3 ^ n - 1`: it is `1` if `n` is odd, and `k + 2` if `n = 2 ^ k * odd`
with `k ≥ 1`. -/
theorem level_three_pow (k t : ℕ) :
    multiplicity 2 (3 ^ (2 ^ k * (2 * t + 1)) - 1) = if k = 0 then 1 else k + 2 := by
  split_ifs with hk
  · subst hk
    exact multiplicity_two_three_pow_sub_one_of_odd ⟨t, by ring⟩
  · have hmul : multiplicity 2 (2 ^ k * (2 * t + 1)) = k := by
      refine multiplicity_eq_of_dvd_of_not_dvd (Dvd.intro _ rfl) fun h ↦ ?_
      rw [pow_succ, Nat.mul_dvd_mul_iff_left (by positivity)] at h
      omega
    rw [multiplicity_two_three_pow_sub_one_of_even
      ((Nat.even_pow.2 ⟨even_two, hk⟩).mul_right _) (by positivity), hmul]

/-- `k + 3 ≤ 2 ^ k` once `3 ≤ k`. -/
theorem add_three_le_two_pow (k : ℕ) (hk : 3 ≤ k) : k + 3 ≤ 2 ^ k := by
  induction k, hk using Nat.le_induction with
  | base => decide
  | succ k _ ih => rw [pow_succ]; omega

/-- **Part 1.** `2 ^ n ∣ 3 ^ n - 1` forces `n = 1`, `2` or `4` (for `n ≥ 1`). -/
theorem only_one_two_four (n : ℕ) (hn : 1 ≤ n) (h : 2 ^ n ∣ 3 ^ n - 1) :
    n = 1 ∨ n = 2 ∨ n = 4 := by
  have h1 : 1 < 3 ^ n := Nat.one_lt_pow (by omega) (by decide)
  have hf : FiniteMultiplicity 2 (3 ^ n - 1) := Nat.finiteMultiplicity_iff.2 ⟨by decide, by omega⟩
  have hle := hf.pow_dvd_iff_le_multiplicity.1 h
  rcases Nat.even_or_odd n with he | ho
  · rw [multiplicity_two_three_pow_sub_one_of_even he (by omega)] at hle
    have hpow : 2 ^ multiplicity 2 n ≤ n := Nat.le_of_dvd hn (pow_multiplicity_dvd 2 n)
    generalize multiplicity 2 n = m at hle hpow
    have hm : m < 3 := by
      by_contra hm
      have := add_three_le_two_pow m (by omega)
      omega
    obtain ⟨r, rfl⟩ := he
    omega
  · rw [multiplicity_two_three_pow_sub_one_of_odd ho] at hle
    omega

/-- And the three survivors really work. -/
theorem one_two_four_work : 2 ^ 1 ∣ 3 ^ 1 - 1 ∧ 2 ^ 2 ∣ 3 ^ 2 - 1 ∧ 2 ^ 4 ∣ 3 ^ 4 - 1 := by decide

/-! ### Part 2: commuting Adams operations force `2 ^ n ∣ 3 ^ n - 1` -/

/-- An element `α * a + β * b` of the glued space's K-theory (`a` the bottom class, `b` the top
class). -/
@[ext]
structure Elt where
  /-- The coefficient of the bottom class `a`. -/
  α : ℤ
  /-- The coefficient of the top class `b`. -/
  β : ℤ

/-- The Adams operation `ψ^k`: `a ↦ k ^ n * a + μ * b`, `b ↦ k ^ (2 * n) * b`. -/
def psi (n : ℕ) (k μ : ℤ) (e : Elt) : Elt :=
  ⟨e.α * k ^ n, e.α * μ + e.β * (k ^ n * k ^ n)⟩

/-- `ψ² ψ³ = ψ³ ψ²` on `a` gives Adams and Atiyah's equation. -/
theorem commute_equation (n : ℕ) (μ₂ μ₃ : ℤ)
    (hc : psi n 2 μ₂ (psi n 3 μ₃ ⟨1, 0⟩) = psi n 3 μ₃ (psi n 2 μ₂ ⟨1, 0⟩)) :
    2 ^ n * (2 ^ n - 1) * μ₃ = 3 ^ n * (3 ^ n - 1) * μ₂ := by
  have hβ := congrArg Elt.β hc
  simp only [psi] at hβ
  linear_combination hβ

/-- **Part 2.** Commuting Adams operations, the mod-2 rule and an odd Hopf invariant `h` force
`n = 1`, `2` or `4`. -/
theorem hopf_invariant_one (n : ℕ) (hn : 1 ≤ n) (h μ₂ μ₃ : ℤ) (hodd : Odd h)
    (hmod : (2 : ℤ) ∣ μ₂ - h)
    (hc : psi n 2 μ₂ (psi n 3 μ₃ ⟨1, 0⟩) = psi n 3 μ₃ (psi n 2 μ₂ ⟨1, 0⟩)) :
    n = 1 ∨ n = 2 ∨ n = 4 := by
  have heq := commute_equation n μ₂ μ₃ hc
  have hμ : Odd μ₂ := by
    obtain ⟨m, rfl⟩ := hodd
    obtain ⟨c, hc⟩ := hmod
    exact ⟨c + m, by omega⟩
  have hdvd : (2 : ℤ) ^ n ∣ ((3 : ℤ) ^ n - 1) * (3 ^ n * μ₂) :=
    ⟨(2 ^ n - 1) * μ₃, by linear_combination (-1 : ℤ) * heq⟩
  have hodd' : Odd ((3 : ℤ) ^ n * μ₂) := (Odd.pow (by decide)).mul hμ
  have hI : (2 : ℤ) ^ n ∣ (3 : ℤ) ^ n - 1 :=
    (Int.isCoprime_two_left.2 hodd').pow_left.dvd_of_dvd_mul_right hdvd
  have h1 : 1 ≤ 3 ^ n := Nat.one_le_pow _ _ (by decide)
  refine only_one_two_four n hn (Int.natCast_dvd_natCast.1 ?_)
  rw [Nat.cast_sub h1]
  exact_mod_cast hI

/-! ### Part 3: on `ℂP²` the Adams operations are mult-add powers -/

/-- `α * x + β * x ^ 2` in `ℤ[x] / (x ^ 3)`: the polynomials without constant term, truncated
above degree `2`. -/
@[ext]
structure Trunc where
  /-- The coefficient of `x`. -/
  α : ℤ
  /-- The coefficient of `x ^ 2`. -/
  β : ℤ
deriving DecidableEq

namespace Trunc

instance : Zero Trunc := ⟨⟨0, 0⟩⟩

instance : Add Trunc := ⟨fun p q ↦ ⟨p.α + q.α, p.β + q.β⟩⟩

instance : Neg Trunc := ⟨fun p ↦ ⟨-p.α, -p.β⟩⟩

/-- Multiplication with `x ^ 3 = 0`: only the product of the two `x`-coefficients survives. -/
instance : Mul Trunc := ⟨fun p q ↦ ⟨0, p.α * q.α⟩⟩

@[simp] theorem zero_α : (0 : Trunc).α = 0 := rfl

@[simp] theorem zero_β : (0 : Trunc).β = 0 := rfl

@[simp] theorem add_α (p q : Trunc) : (p + q).α = p.α + q.α := rfl

@[simp] theorem add_β (p q : Trunc) : (p + q).β = p.β + q.β := rfl

@[simp] theorem neg_α (p : Trunc) : (-p).α = -p.α := rfl

@[simp] theorem neg_β (p : Trunc) : (-p).β = -p.β := rfl

@[simp] theorem mul_α (p q : Trunc) : (p * q).α = 0 := rfl

@[simp] theorem mul_β (p q : Trunc) : (p * q).β = p.α * q.α := rfl

/-- `Trunc` is a commutative ring without unit: the ideal `(x)` of `ℤ[x] / (x ^ 3)`. -/
instance : NonUnitalCommRing Trunc where
  add_assoc p q r := by ext <;> simp [add_assoc]
  zero_add p := by ext <;> simp
  add_zero p := by ext <;> simp
  add_comm p q := by ext <;> simp [add_comm]
  neg_add_cancel p := by ext <;> simp
  left_distrib p q r := by ext <;> simp [mul_add]
  right_distrib p q r := by ext <;> simp [add_mul]
  zero_mul p := by ext <;> simp
  mul_zero p := by ext <;> simp
  mul_assoc p q r := by ext <;> simp
  mul_comm p q := by ext <;> simp [mul_comm]
  nsmul := nsmulRec
  zsmul := zsmulRec

/-- The generator `x`. -/
def x : Trunc := ⟨1, 0⟩

@[simp] theorem x_α : x.α = 1 := rfl

@[simp] theorem x_β : x.β = 0 := rfl

/-- Mult-add on `Trunc` is `(1 + p) * (1 + q) - 1` with `x ^ 3 = 0`. -/
theorem multAdd_eq (p q : Trunc) : p ✳ q = ⟨p.α + q.α, p.β + q.β + p.α * q.α⟩ := by
  ext <;> simp [multAdd_def]

end Trunc


/-- `(1 + x) ^ k - 1 = k * x + k.choose 2 * x ^ 2`, written without division. -/
theorem mpow_formula (k : ℕ) :
    (mpow Trunc.x k).α = k ∧ 2 * (mpow Trunc.x k).β = k * (k - 1 : ℤ) := by
  induction k with
  | zero => simp
  | succ k ih =>
    obtain ⟨h1, h2⟩ := ih
    rw [mpow_succ, Trunc.multAdd_eq]
    refine ⟨by simp [h1], ?_⟩
    simp only [Trunc.x_α, Trunc.x_β, add_zero, mul_one, h1]
    push_cast
    linear_combination h2

/-- `ψ^j` applied to `α * x + β * x ^ 2`, using `ψ^j (x) = mpow Trunc.x j` and
`ψ^j (x ^ 2) = (mpow Trunc.x j) ^ 2`. -/
def psiT (j : ℕ) (p : Trunc) : Trunc :=
  ⟨p.α * (mpow Trunc.x j).α,
    p.α * (mpow Trunc.x j).β + p.β * ((mpow Trunc.x j).α * (mpow Trunc.x j).α)⟩

/-- Mult-add powers compose: `ψ^j ψ^k = ψ^(j * k)`, for every `j`, `k`. -/
theorem psi_compose (j k : ℕ) : psiT j (mpow Trunc.x k) = mpow Trunc.x (j * k) := by
  obtain ⟨a1, b1⟩ := mpow_formula j
  obtain ⟨a2, b2⟩ := mpow_formula k
  obtain ⟨a3, b3⟩ := mpow_formula (j * k)
  ext
  · simp only [psiT, a1, a2, a3]
    push_cast
    ring
  · refine mul_left_cancel₀ (two_ne_zero (α := ℤ)) ?_
    simp only [psiT, a1, a2, b3]
    push_cast
    linear_combination (k : ℤ) * b1 + (j : ℤ) * j * b2

/-- `ℂP²` meets every hypothesis of Part 2 with `n = 1`, `h = 1`: `μ₂ = 2.choose 2 = 1`,
`μ₃ = 3.choose 2 = 3`, and `ψ² ψ³ = ψ³ ψ²`. -/
theorem cp2_case :
    mpow Trunc.x 2 = ⟨2, 1⟩ ∧ mpow Trunc.x 3 = ⟨3, 3⟩ ∧
    psi 1 2 1 (psi 1 3 3 ⟨1, 0⟩) = psi 1 3 3 (psi 1 2 1 ⟨1, 0⟩) ∧
    (2 : ℤ) ∣ 1 - 1 := by
  refine ⟨?_, ?_, by simp [psi], by simp⟩ <;> simp [mpow_succ, Trunc.multAdd_eq]

end MultAdd.HopfOne
