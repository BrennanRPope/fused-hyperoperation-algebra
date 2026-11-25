/-
  Only 1, 2 and 4: the number step of the Hopf invariant one theorem
  ================================================================

  Part of Fused Mult-Add Algebra. Companion to the topology page and to
  ../../tesseract/ToricCharts.lean. Python checks: ../checks.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`. Run it with
  `lean HopfInvariantOne.lean`.

  The theorem (Adams 1960; this short proof Adams and Atiyah 1966):
  a map S^(4n-1) → S^(2n) of odd Hopf invariant exists only for
  n = 1, 2, 4. So only the spheres S^1, S^3, S^7 can carry a continuous
  multiplication with a unit, and the only real division algebras have
  dimension 1, 2, 4, 8.

  What is PROVED here, from scratch:

  Part 1 (numbers). 2^n divides 3^n - 1 exactly when n = 1, 2 or 4
    (for n ≥ 1). Proof: track the exact power of 2 in 3^n - 1. It is 2
    when n is odd, and 2^(k+2) when n = 2^k · odd with k ≥ 1. That is
    far too small once n ≥ 5.

  Part 2 (mult-add powers). On the K-theory of the glued space, with a
    bottom class a and a top class b = a² / h, the Adams operations act
    by ψ^k(a) = k^n a + μ_k b and ψ^k(b) = k^(2n) b. If ψ² and ψ³
    commute, and ψ²(a) ≡ a² mod 2, and h is odd, then 2^n | 3^n - 1,
    so n is 1, 2 or 4.

  Part 3 (the mult-add cube). On CP², which is the n = 1 case built on
    the Hopf map, the Adams operations ARE truncated mult-add powers:
    the k-fold mult-add sum x ◯ x ◯ … ◯ x in Z[x]/(x³) is
    k x + C(k,2) x². They compose (ψ^j ψ^k = ψ^(jk)), satisfy the mod 2
    rule with h = 1, and so meet every hypothesis of Part 2.

  What is TAKEN AS GIVEN (topology, not formalised here): that the
  glued space's K-theory has this two-generator shape, that the Adams
  operations exist on it, are additive and compose, that ψ^k acts on a
  2m-cell by k^m, and that ψ^p(x) ≡ x^p mod p. These are the standard
  K-theory facts (Atiyah, Hirzebruch, Adams); for line bundles they are
  exactly the statements that ψ^k is the k-th mult-add power and that
  the middle of the Pascal row vanishes mod p.
-/

namespace CircleOp.HopfOne

/-! ## Part 1: the exact power of 2 in 3^n − 1 -/

/-- `OneMod j x`: x is one more than 2^j times an odd number. -/
def OneMod (j x : Nat) : Prop := ∃ u, x = 2 ^ j * (2 * u + 1) + 1

/-- Raising to an odd power keeps the level (for j ≥ 1). -/
theorem oneMod_oddPow {j A : Nat} (hj : 1 ≤ j) (hA : OneMod j A) :
    ∀ t, OneMod j (A ^ (2 * t + 1)) := by
  obtain ⟨u, rfl⟩ := hA
  obtain ⟨i, rfl⟩ : ∃ i, j = i + 1 := ⟨j - 1, by omega⟩
  intro t
  induction t with
  | zero => exact ⟨u, by simp⟩
  | succ t ih =>
    obtain ⟨v, hv⟩ := ih
    refine ⟨(2 * (2 ^ i * (2 * u + 1)) * (2 ^ i * (2 * u + 1)) + 2 * (2 ^ i * (2 * u + 1))) * (2 * v + 1)
        + v + 2 ^ i * (2 * u + 1) * (2 * u + 1) + (2 * u + 1), ?_⟩
    have e : 2 * (t + 1) + 1 = (2 * t + 1) + 2 := by omega
    rw [e, Nat.pow_add, hv, Nat.pow_succ 2 i]
    grind

/-- Squaring raises the level by one (for j ≥ 2). -/
theorem oneMod_sq {j A : Nat} (hj : 2 ≤ j) (hA : OneMod j A) : OneMod (j + 1) (A ^ 2) := by
  obtain ⟨u, rfl⟩ := hA
  obtain ⟨i, rfl⟩ : ∃ i, j = i + 2 := ⟨j - 2, by omega⟩
  refine ⟨u + 2 ^ i * (2 * u + 1) * (2 * u + 1), ?_⟩
  rw [show i + 2 + 1 = i + 3 by omega, Nat.pow_add 2 i 3, Nat.pow_add 2 i 2]
  grind

/-- 3^(2^k) is one more than 2^(k+2) times an odd number, for k ≥ 1. -/
theorem three_pow_two_pow : ∀ k, 1 ≤ k → OneMod (k + 2) (3 ^ (2 ^ k)) := by
  intro k hk
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.eq_zero_or_pos k with rfl | hpos
    · exact ⟨0, by decide⟩
    · have h := oneMod_sq (by omega) (ih hpos)
      rw [show 2 ^ (k + 1) = 2 ^ k * 2 from Nat.pow_succ 2 k, Nat.pow_mul]
      exact h

/-- Every positive n is 2^k times an odd number. -/
theorem two_adic_split : ∀ n, 1 ≤ n → ∃ k t, n = 2 ^ k * (2 * t + 1) := by
  intro n
  induction n using Nat.strongRecOn with
  | _ n ih =>
    intro hn
    rcases Nat.mod_two_eq_zero_or_one n with h0 | h1
    · obtain ⟨k, t, ht⟩ := ih (n / 2) (by omega) (by omega)
      have hm : n = 2 * (n / 2) := by omega
      exact ⟨k + 1, t, by rw [hm, ht, Nat.pow_succ]; grind⟩
    · exact ⟨0, n / 2, by simp; omega⟩

/-- The level of 3^n: 1 if n is odd, k + 2 if n = 2^k · odd with k ≥ 1. -/
theorem level_three_pow (k t : Nat) :
    OneMod (if k = 0 then 1 else k + 2) (3 ^ (2 ^ k * (2 * t + 1))) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simpa using oneMod_oddPow (j := 1) (A := 3) (Nat.le_refl 1) ⟨0, rfl⟩ t
  · have hk0 : k ≠ 0 := by omega
    simp only [hk0, ↓reduceIte]
    rw [Nat.pow_mul]
    exact oneMod_oddPow (by omega) (three_pow_two_pow k hk) t

/-- If 2^n divides 2^e times an odd number, then n ≤ e. -/
theorem le_of_dvd_two_pow_odd {n e u : Nat} (h : 2 ^ n ∣ 2 ^ e * (2 * u + 1)) : n ≤ e := by
  apply Nat.le_of_not_lt
  intro hlt
  have h1 : 2 ^ (e + 1) ∣ 2 ^ e * (2 * u + 1) :=
    Nat.dvd_trans (Nat.pow_dvd_pow 2 hlt) h
  obtain ⟨c, hc⟩ := h1
  rw [Nat.pow_succ, Nat.mul_assoc] at hc
  have := Nat.eq_of_mul_eq_mul_left (Nat.pow_pos (by decide : 0 < 2)) hc
  omega

/-- 2^k ≥ k + 3 once k ≥ 3. -/
theorem pow_ge (k : Nat) (hk : 3 ≤ k) : k + 3 ≤ 2 ^ k := by
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.lt_or_ge k 3 with h | h
    · have : k = 2 := by omega
      subst this; decide
    · have := ih h; rw [Nat.pow_succ]; omega

/-- **Part 1.** 2^n | 3^n − 1 forces n = 1, 2 or 4 (for n ≥ 1). -/
theorem only_one_two_four (n : Nat) (hn : 1 ≤ n) (h : 2 ^ n ∣ 3 ^ n - 1) :
    n = 1 ∨ n = 2 ∨ n = 4 := by
  obtain ⟨k, t, hkt⟩ := two_adic_split n hn
  have hlev := level_three_pow k t
  rw [← hkt] at hlev
  obtain ⟨u, hu⟩ := hlev
  have hle : n ≤ (if k = 0 then 1 else k + 2) := by
    apply le_of_dvd_two_pow_odd (u := u)
    rw [hu, Nat.add_sub_cancel] at h
    exact h
  have hbig : 2 ^ k ≤ n := by
    rw [hkt]; exact Nat.le_mul_of_pos_right _ (by omega)
  have hn4 : n ≤ 4 := by
    by_cases hk : k = 0
    · simp only [hk, ↓reduceIte] at hle; omega
    · simp only [hk, ↓reduceIte] at hle
      rcases Nat.lt_or_ge k 3 with h3 | h3
      · omega
      · have := pow_ge k h3; omega
  have : n = 1 ∨ n = 2 ∨ n = 3 ∨ n = 4 := by omega
  rcases this with rfl | rfl | rfl | rfl
  · exact Or.inl rfl
  · exact Or.inr (Or.inl rfl)
  · exact absurd h (by decide)
  · exact Or.inr (Or.inr rfl)

/-- And the three survivors really work. -/
theorem one_two_four_work : 2 ^ 1 ∣ 3 ^ 1 - 1 ∧ 2 ^ 2 ∣ 3 ^ 2 - 1 ∧ 2 ^ 4 ∣ 3 ^ 4 - 1 := by decide

/-! ## Part 2: commuting mult-add powers force 2^n | 3^n − 1 -/

/-- An element α·a + β·b of the glued space's K-theory (a bottom, b top). -/
structure Elt where
  α : Int
  β : Int

/-- The Adams operation ψ^k: a ↦ k^n a + μ b, b ↦ k^(2n) b. -/
def psi (n : Nat) (k μ : Int) (e : Elt) : Elt :=
  ⟨e.α * k ^ n, e.α * μ + e.β * (k ^ n * k ^ n)⟩

/-- ψ² ψ³ = ψ³ ψ² on a gives Adams and Atiyah's equation. -/
theorem commute_equation (n : Nat) (μ₂ μ₃ : Int)
    (hc : psi n 2 μ₂ (psi n 3 μ₃ ⟨1, 0⟩) = psi n 3 μ₃ (psi n 2 μ₂ ⟨1, 0⟩)) :
    2 ^ n * (2 ^ n - 1) * μ₃ = 3 ^ n * (3 ^ n - 1) * μ₂ := by
  have hβ := congrArg Elt.β hc
  simp only [psi] at hβ
  grind

/-- An odd integer: 2m + 1. -/
def IsOdd (y : Int) : Prop := ∃ m, y = 2 * m + 1

theorem odd_mul {x y : Int} (hx : IsOdd x) (hy : IsOdd y) : IsOdd (x * y) := by
  obtain ⟨a, rfl⟩ := hx; obtain ⟨b, rfl⟩ := hy
  exact ⟨2 * a * b + a + b, by grind⟩

theorem odd_three_pow : ∀ n : Nat, IsOdd ((3 : Int) ^ n) := by
  intro n
  induction n with
  | zero => exact ⟨0, by decide⟩
  | succ n ih => rw [Int.pow_succ]; exact odd_mul ih ⟨1, by decide⟩

/-- Euclid's step for powers of 2: an odd factor can be cancelled. -/
theorem two_pow_dvd_of_mul_odd : ∀ (n : Nat) (x y : Int), IsOdd y → (2 : Int) ^ n ∣ x * y → (2 : Int) ^ n ∣ x := by
  intro n
  induction n with
  | zero => intro x y _ _; exact ⟨x, by simp⟩
  | succ n ih =>
    intro x y hy h
    obtain ⟨m, rfl⟩ := hy
    obtain ⟨c, hc⟩ := h
    have hx2 : x % 2 = 0 := by
      rcases Int.emod_two_eq x with hx | hx
      · exact hx
      · exfalso
        have hx' : x = 2 * (x / 2) + 1 := by omega
        rw [Int.pow_succ] at hc
        have : 2 * (2 * (x / 2) * m + x / 2 + m) + 1 = 2 * (2 ^ n * c) := by
          rw [hx'] at hc; grind
        omega
    have hx' : x = 2 * (x / 2) := by omega
    have : (2 : Int) ^ n ∣ (x / 2) * (2 * m + 1) := by
      refine ⟨c, ?_⟩
      rw [Int.pow_succ] at hc
      have h2 : 2 * ((x / 2) * (2 * m + 1)) = 2 * (2 ^ n * c) := by
        rw [hx'] at hc; grind
      omega
    obtain ⟨d, hd⟩ := ih (x / 2) (2 * m + 1) ⟨m, rfl⟩ this
    exact ⟨d, by rw [hx', hd, Int.pow_succ]; grind⟩

/-- **Part 2.** Commuting Adams operations, the mod-2 rule and an odd Hopf
    invariant h force n = 1, 2 or 4. -/
theorem hopf_invariant_one (n : Nat) (hn : 1 ≤ n) (h μ₂ μ₃ : Int)
    (hodd : IsOdd h)
    (hmod : (2 : Int) ∣ μ₂ - h)
    (hc : psi n 2 μ₂ (psi n 3 μ₃ ⟨1, 0⟩) = psi n 3 μ₃ (psi n 2 μ₂ ⟨1, 0⟩)) :
    n = 1 ∨ n = 2 ∨ n = 4 := by
  have heq := commute_equation n μ₂ μ₃ hc
  have hμ : IsOdd μ₂ := by
    obtain ⟨m, rfl⟩ := hodd; obtain ⟨c, hc⟩ := hmod
    exact ⟨c + m, by omega⟩
  have hdvd : (2 : Int) ^ n ∣ ((3 : Int) ^ n - 1) * (3 ^ n * μ₂) :=
    ⟨(2 ^ n - 1) * μ₃, by rw [← Int.mul_assoc, Int.mul_comm (3 ^ n - 1) (3 ^ n), ← heq]; grind⟩
  have hI : (2 : Int) ^ n ∣ (3 : Int) ^ n - 1 :=
    two_pow_dvd_of_mul_odd n _ _ (odd_mul (odd_three_pow n) hμ) hdvd
  have h1 : 1 ≤ 3 ^ n := Nat.pow_pos (by decide)
  have hN : 2 ^ n ∣ 3 ^ n - 1 := by
    have hcast : ((3 ^ n - 1 : Nat) : Int) = (3 : Int) ^ n - 1 := by
      rw [Int.ofNat_sub h1]; simp
    rw [← hcast, show (2 : Int) ^ n = ((2 ^ n : Nat) : Int) by simp] at hI
    exact Int.ofNat_dvd.mp hI
  exact only_one_two_four n hn hN

/-! ## Part 3: on CP² the Adams operations are mult-add powers -/

/-- αx + βx² in Z[x]/(x³). -/
structure Trunc where
  α : Int
  β : Int
deriving DecidableEq

/-- Mult-add: (1+p)(1+q) − 1 with x³ = 0. -/
def circ (p q : Trunc) : Trunc := ⟨p.α + q.α, p.β + q.β + p.α * q.α⟩

/-- The k-fold mult-add sum x ◯ x ◯ … ◯ x, i.e. (1 + x)^k − 1. -/
def mpow : Nat → Trunc
  | 0 => ⟨0, 0⟩
  | k + 1 => circ (mpow k) ⟨1, 0⟩

/-- (1 + x)^k − 1 = k x + C(k,2) x², written without division. -/
theorem mpow_formula (k : Nat) : (mpow k).α = k ∧ 2 * (mpow k).β = k * (k - 1 : Int) := by
  induction k with
  | zero => simp [mpow]
  | succ k ih =>
    obtain ⟨h1, h2⟩ := ih
    simp only [mpow, circ]
    constructor
    · rw [h1]; simp
    · simp only [Int.add_zero, Int.mul_one]
      rw [h1, Int.mul_add, h2]; simp only [Int.natCast_add, Int.natCast_one]; grind

/-- ψ^j applied to αx + βx², using ψ^j(x) = mpow j and ψ^j(x²) = (mpow j)². -/
def psiT (j : Nat) (p : Trunc) : Trunc :=
  ⟨p.α * (mpow j).α, p.α * (mpow j).β + p.β * ((mpow j).α * (mpow j).α)⟩

/-- Mult-add powers compose: ψ^j ψ^k = ψ^(jk), for every j, k. -/
theorem psi_compose (j k : Nat) : psiT j (mpow k) = mpow (j * k) := by
  obtain ⟨a1, b1⟩ := mpow_formula j
  obtain ⟨a2, b2⟩ := mpow_formula k
  obtain ⟨a3, b3⟩ := mpow_formula (j * k)
  have : 2 * (psiT j (mpow k)).β = 2 * (mpow (j * k)).β := by
    simp only [psiT]
    rw [Int.mul_add, a1, a2, b3]
    have e1 : 2 * ((k : Int) * (mpow j).β) = k * (2 * (mpow j).β) := by grind
    have e2 : 2 * ((mpow k).β * ((j : Int) * j)) = (2 * (mpow k).β) * (j * j) := by grind
    rw [e1, e2, b1, b2]; push_cast; grind
  cases hp : psiT j (mpow k) with
  | mk p q =>
    cases hm : mpow (j * k) with
    | mk r s =>
      rw [hp, hm] at this
      simp only [psiT] at hp
      have hα : p = r := by
        have := congrArg Trunc.α hp; have := congrArg Trunc.α hm
        simp_all; push_cast; grind
      simp at this; subst hα; simp; omega

/-- CP² meets every hypothesis of Part 2 with n = 1, h = 1:
    μ₂ = C(2,2) = 1, μ₃ = C(3,2) = 3, and ψ² ψ³ = ψ³ ψ². -/
theorem cp2_case :
    mpow 2 = ⟨2, 1⟩ ∧ mpow 3 = ⟨3, 3⟩ ∧
    psi 1 2 1 (psi 1 3 3 ⟨1, 0⟩) = psi 1 3 3 (psi 1 2 1 ⟨1, 0⟩) ∧
    (2 : Int) ∣ 1 - 1 := by
  refine ⟨by decide, by decide, by simp [psi], ⟨0, by decide⟩⟩

end CircleOp.HopfOne
