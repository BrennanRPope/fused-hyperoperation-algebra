/-
  New ties, round 2: the union fractal. Core Lean 4.37 only (no mathlib), no `sorry`.
  Self-contained (no imports), so it can sit in the shared library as FusedAlgebra/Topics/NewTies.lean.

  1. A union in a box of size M is a stretch about the point M (the box full).
  2. Stretches about one centre commute; about two centres they differ by (c₁ − c₂)(1 − s₁)(1 − s₂).
  3. So two stretches commute exactly when they share a centre (slopes ≠ 1).
  4. The key identity of the lower bound: (1 − P)(1 − Q) = PQ ↔ P + Q = 1.
  5. With one map per box the dimension equation is Moran's.
  6. A small exact count: three integer stretches, two sharing a centre, give 3, 8, 21 distinct composites
     (the shared-centre count), not 3, 9, 27.
  7. Plane laws: a law's maps fix w exactly when d·w² + k·w + 1 = 0, and that line is tangent to the border.
-/
namespace NewTies

/-- A stretch about the centre c by the factor s: x ↦ c + s(x − c). -/
def hom (c s x : Int) : Int := c + s * (x - c)

/-- Every stretch fixes its centre. -/
theorem hom_centre (c s : Int) : hom c s c = c := by simp only [hom]; grind

/-- 1. A union in a box of size M, scaled by M: M·(a + x − ax/M) = M·M + (M − a)(x − M),
    i.e. the stretch about M (the full box) by the factor 1 − a/M. -/
theorem union_is_stretch (M a x : Int) : M * a + M * x - a * x = M * M + (M - a) * (x - M) := by grind

/-- 2a. Stretches about one centre commute. -/
theorem same_centre_commute (c s t x : Int) : hom c s (hom c t x) = hom c t (hom c s x) := by
  simp only [hom]; grind

/-- 2b. The order gap between stretches about two centres. -/
theorem order_gap (c₁ s₁ c₂ s₂ x : Int) :
    hom c₁ s₁ (hom c₂ s₂ x) - hom c₂ s₂ (hom c₁ s₁ x) = (c₁ - c₂) * ((1 - s₁) * (1 - s₂)) := by
  simp only [hom]; grind

/-- 3. Two stretches with slopes ≠ 1 commute exactly when they share a centre. -/
theorem commute_iff (c₁ s₁ c₂ s₂ : Int) (h₁ : s₁ ≠ 1) (h₂ : s₂ ≠ 1) :
    (∀ x, hom c₁ s₁ (hom c₂ s₂ x) = hom c₂ s₂ (hom c₁ s₁ x)) ↔ c₁ = c₂ := by
  constructor
  · intro h
    have g := order_gap c₁ s₁ c₂ s₂ 0
    rw [h 0, Int.sub_self] at g
    have z : (c₁ - c₂) * ((1 - s₁) * (1 - s₂)) = 0 := g.symm
    rcases Int.mul_eq_zero.mp z with e | e
    · grind
    · rcases Int.mul_eq_zero.mp e with e' | e'
      · exact absurd (by grind) h₁
      · exact absurd (by grind) h₂
  · intro e x; subst e; exact same_centre_commute c₁ s₁ s₂ x

/-- 4. The lower bound's key step: the first-return generators have total weight 1 exactly when
    the shared-centre dimension equation holds, (1/P − 1)(1/Q − 1) = 1 ↔ P + Q = 1, cleared of fractions. -/
theorem lower_bound_key (P Q : Int) : (1 - P) * (1 - Q) = P * Q ↔ P + Q = 1 := by
  constructor <;> intro h <;> grind

/-- 5. One map per box: 1 + Σ[(1 − wᵢ) − 1] = 1 − Σ wᵢ, Moran's equation. -/
theorem moran_special (ws : List Int) :
    1 + (ws.map (fun w => (1 - w) - 1)).sum = 1 - ws.sum := by
  induction ws with
  | nil => simp
  | cons w ws ih => simp only [List.map_cons, List.sum_cons]; grind

/-! ## 6. An exact count -/

/-- A map x ↦ s·x + b as a pair (s, b); composing p after q. -/
def after (p q : Int × Int) : Int × Int := (p.1 * q.1, p.1 * q.2 + p.2)

/-- Three integer stretches: about 0 by 2, about 0 by 3, about 1 by 5 (x ↦ 5x − 4). -/
def gens : List (Int × Int) := [(2, 0), (3, 0), (5, -4)]

/-- All composites of length n (with repeats). -/
def words : Nat → List (Int × Int)
  | 0 => [(1, 0)]
  | n + 1 => (words n).flatMap (fun w => gens.map (fun g => after g w))

/-- Distinct composites of length 1, 2, 3: 3, 8, 21, the shared-centre count (not 3, 9, 27). -/
theorem count_one : (words 1).eraseDups.length = 3 := by decide
theorem count_two : (words 2).eraseDups.length = 8 := by decide
theorem count_three : (words 3).eraseDups.length = 21 := by decide

/-! ## 7. Plane laws: sharing a wall is a tangent line -/

/-- The plane-law map x ↦ (x + a + kax)/(1 − dax) fixes w exactly when d·w² + k·w + 1 = 0 (cleared of the fraction). -/
theorem fixes_wall (k d a w : Int) :
    w * (1 - d * a * w) - (w + a + k * a * w) = -(a * (d * w * w + k * w + 1)) := by grind

/-- On the border 4d = k², the wall line d·w² + k·w + 1 = 0 becomes (k·w + 2)² = 0 (times 4): a double root,
    so the laws sharing the wall w form the tangent line to the border at k = −2/w. -/
theorem wall_line_tangent (k w : Int) : k * k * w * w + 4 * (k * w) + 4 = (k * w + 2) * (k * w + 2) := by grind

/-! ## 8. Off the tangent line: the two-law threshold -/

/-- Two laws with block weights x₁, x₂: the switching matrix [[0, x₁], [x₂, 0]] has spectral radius 1 exactly when
    x₁·x₂ = 1, which is the same as x₁/(1+x₁) + x₂/(1+x₂) = 1 (cleared of fractions). -/
theorem two_law_threshold (x₁ x₂ : Int) :
    (1 + x₁) * (1 + x₂) = x₁ * (1 + x₂) + x₂ * (1 + x₁) ↔ x₁ * x₂ = 1 := by
  constructor <;> intro h <;> grind

/-! ## 9. No accidental collisions: the two-prime ping-pong example -/

/-- Law 1 of the ping-pong example has multipliers 5/18 and 5/24. A ratio of two of its blocks is
    (5/18)^m (5/24)^n with whole numbers m, n; its 2-adic and 3-adic valuations are −m − 3n and −2m − n.
    Both vanish only for m = n = 0, so a law-1 block that is not the identity moves 2-adically or 3-adically. -/
theorem valuations_separate (m n : Int) (h2 : -m - 3 * n = 0) (h3 : -2 * m - n = 0) : m = 0 ∧ n = 0 := by
  omega

/-- XOR with an independent set of density a, in a box of size 1, is union in a box of size 1/2:
    a + x − 2ax (cleared of the fraction 1/2). -/
theorem xor_is_half_box_union (a x : Int) : 2 * (a + x - 2 * a * x) = 2 * a + 2 * x - 4 * a * x := by grind

/-- De Morgan for densities: intersecting with density a is complement ∘ union with density 1 − a ∘ complement. -/
theorem and_from_or (a x : Int) : 1 - ((1 - a) + (1 - x) - (1 - a) * (1 - x)) = a * x := by grind

/-! ## 10. No accidental collisions for the original examples (round 6) -/

/-- 2-adic step for even walls: every power of 5 is 1 more than a multiple of 4, so a run of the letter a alone
    (multiplier 5^(-α)) moves z₁ by a factor that is 1 mod 4 and cannot leave the 2-adic ball the proof uses. -/
theorem five_pow_mod_four (n : Nat) : 5 ^ n % 4 = 1 := by
  induction n with
  | zero => decide
  | succ n ih => rw [Nat.pow_succ, Nat.mul_mod, ih]

/-- For even Q the number (Q + 1)² is odd, so κ = 4Q/(Q + 1)² has 2-adic valuation 2 + v₂(Q) ≥ 3. -/
theorem even_Q_square_odd (Q : Nat) (h : Q % 2 = 0) : (Q + 1) * (Q + 1) % 2 = 1 := by
  rw [Nat.mul_mod, Nat.add_mod, h]

/-- The last step of the 2-adic proof: every remainder lands at one 2-adic level s, and the ratio map shifts the
    level by 2β with β ≥ 1, so the two sides of a supposed collision can't both be at level s. -/
theorem levels_differ (s β : Int) (hβ : 1 ≤ β) : s + 2 * β ≠ s := by omega

/-- 5-adic step at walls ±3: an odd run of the letter c multiplies z₂ by 4^(−k) ≡ −1 (mod 5) ... -/
theorem four_pow_odd_mod_five (j : Nat) : 4 ^ (2 * j + 1) % 5 = 4 := by
  induction j with
  | zero => decide
  | succ j ih =>
    have : 2 * (j + 1) + 1 = (2 * j + 1) + 2 := by omega
    rw [this, Nat.pow_add, Nat.mul_mod, ih]

/-- ... and an even run multiplies it by 4^(−k) ≡ 1 (mod 5). -/
theorem four_pow_even_mod_five (j : Nat) : 4 ^ (2 * j) % 5 = 1 := by
  induction j with
  | zero => decide
  | succ j ih =>
    have : 2 * (j + 1) = 2 * j + 2 := by omega
    rw [this, Nat.pow_add, Nat.mul_mod, ih]

/-- The real-number step at walls ±3: the ratios (3/4)^n, 5·(3/4)^n and (3/4)^n/5 all miss the window
    [12/13, 13/12]. Since each family decreases in n, it is enough that 5·(3/4)^5 > 13/12, 5·(3/4)^6 < 12/13 and
    3/4 < 12/13 (cleared of fractions). -/
theorem window_miss : 12 * 5 * 3 ^ 5 > 13 * 4 ^ 5 ∧ 13 * 5 * 3 ^ 6 < 12 * 4 ^ 6 ∧ 13 * 3 < 12 * 4 := by decide

end NewTies
