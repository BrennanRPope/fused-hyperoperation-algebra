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

/-! ## Round 7: the plane union fractal -/

/-- Gaussian integers mod 3 have no zero divisors (3 stays prime in ℤ[i], so ℤ[i]/3 is the field with 9 elements).
    This is what makes the 3-adic valuation of a Gaussian rational well defined and additive. -/
theorem gauss_mod3_no_zero_divisors :
    ∀ a b c d : Fin 3, (a ≠ 0 ∨ b ≠ 0) → (c ≠ 0 ∨ d ≠ 0) → (a * c - b * d ≠ 0 ∨ a * d + b * c ≠ 0) := by decide

/-- The numbers in the multipliers other than 3 are 3-adic units: 10, 4, and 1 + i (its norm is 2). The centres
    0, 1, i are distinct mod 3 (1 − i has norm 2). -/
theorem units_mod3 : 10 % 3 ≠ 0 ∧ 4 % 3 ≠ 0 ∧ 2 % 3 ≠ 0 := by decide

/-- The modified example's class-2 and class-3 maps, 1 + (3i/10)(z − 1) and i + (3/10)(z − i), send every 3-adic
    integer z to 1 (resp. i) plus a multiple of 3, so to a 3-adic unit; class-1 maps multiply by (1+i)/4 or 3/10,
    keeping 3-adic integers. Then g′(1) has valuation 0, while λ₁^α λ₂^(−β) has valuation −β: -/
theorem three_adic_levels_differ (β : Int) (hβ : 1 ≤ β) : (0 : Int) - β ≠ 0 := by omega

/-- Upper bound for √2 used in every disc radius: (99/70)² > 2. -/
theorem sqrt_two_upper : 99 * 99 > 2 * 70 * 70 := by decide

/-- The size bound that ends the original example's last step for long a-runs: for d ≥ 8 and e ≤ 6,
    |λ₁|^d (10/3)^e · 1.15 ≤ (99/280)^8 (10/3)^6 · 1.15 < 0.59, and 0.59 is below every |w| − r on the class-2/3
    covers, 1.15 above every |w| + r (plane.py, C7). Cleared of fractions: -/
theorem tail_bound : 99 ^ 8 * 10 ^ 6 * 115 < 59 * 280 ^ 8 * 3 ^ 6 := by decide

/-- Linear irreducibility for Hochman's theorem: a 45° turn, matrix [[1, −1], [1, 1]] up to scale, has
    characteristic polynomial t² − 2t + 2, whose discriminant is negative, so it fixes no line. -/
theorem rot45_no_invariant_line : (-2 : Int) * (-2) - 4 * 1 * 2 < 0 := by decide

/-- First-return pieces: with p = 1/P₁ and q = 1/Q (the reciprocal growth series of class 1 and of the rest),
    (1/p − 1)(1/q − 1) = 1 exactly when p + q = 1, i.e. at the root of the trace-monoid formula. Cleared of the
    denominators pq: -/
theorem first_return_threshold (p q : Int) : (1 - p) * (1 - q) = p * q ↔ p + q = 1 := by
  constructor <;> intro h <;> grind


/-! ### Round 8: a general no-collision test (README section 12) -/

/-- Lemma A′, invariance: with ε = a/b, the disc about c₁ of radius R is sent into itself by every map exactly at
    R(1 − ε) = D(1 + ε). Cleared of fractions: -/
theorem lemmaA_invariant (a b D R : Int) : D * (b + a) + a * R = b * R ↔ R * (b - a) = D * (b + a) := by
  constructor <;> intro h <;> grind

/-- Lemma A′, separation: the class pieces sit in discs of radius ε(R + D) = 2εD/(1 − ε) about their centres, so they
    are disjoint when 4εD/(1 − ε) < δ, which is ε(4D + δ) < δ. Cleared of fractions (ε = a/b): -/
theorem lemmaA_key (a b D δ : Int) : 4 * a * D < δ * (b - a) ↔ a * (4 * D + δ) < b * δ := by
  constructor <;> intro h <;> grind

/-- Test (B) for a class of two maps and one good prime, valuations (0, v) with v ≥ 1: every integer relation
    w₁·0 + w₂·v = 0 has w₂ = 0, so the relations are the multiples of (1, 0), one-signed. -/
theorem two_maps_one_prime (w₁ w₂ v : Int) (hv : 1 ≤ v) (h : w₁ * 0 + w₂ * v = 0) : w₂ = 0 := by
  have : w₂ * v = 0 := by grind
  rcases Int.mul_eq_zero.mp this with h | h <;> omega

/-- Corollary P9's recipe: λ_j carries its own prime p_j (valuation ≥ 1) and the others' valuations there are 0, so a
    relation must vanish on the first k − 1 coordinates and is a multiple of the last unit vector, one-signed.
    The case k = 3: -/
theorem diagonal_kernel (w₁ w₂ a b : Int) (ha : 1 ≤ a) (hb : 1 ≤ b) (h₁ : w₁ * a = 0) (h₂ : w₂ * b = 0) :
    w₁ = 0 ∧ w₂ = 0 := by
  rcases Int.mul_eq_zero.mp h₁ with h | h <;> rcases Int.mul_eq_zero.mp h₂ with k | k <;> omega

/-- The good-prime step: a word whose outer map has centre c_j and a multiplier divisible by p sends a p-integral
    point to c_j + p·t, which differs from c_i by c_j − c_i mod p, a unit. -/
theorem good_prime_unit (cj ci p t : Int) : (cj - ci + p * t) % p = (cj - ci) % p := by
  simp [Int.add_mul_emod_self_left]

/-- Centres 0, 1, …, N − 1 are pairwise distinct mod every prime p ≥ N. -/
theorem centres_distinct_mod (a b p : Nat) (hab : a < b) (hbp : b < p) : (b - a) % p ≠ 0 := by
  rw [Nat.mod_eq_of_lt (by omega)]; omega

/-- Density of the family: every residue class mod P has a member in every window of length P, so the fractions
    m/2ⁿ with m in a prescribed class come within P/2ⁿ of any target. -/
theorem exists_near (T r P : Int) (hP : 0 < P) : ∃ m, m % P = r % P ∧ T ≤ m ∧ m < T + P := by
  refine ⟨T + (r - T) % P, ?_, ?_, ?_⟩
  · rw [Int.add_emod, Int.emod_emod_of_dvd _ (Int.dvd_refl P), ← Int.add_emod]; congr 1; omega
  · have := Int.emod_nonneg (r - T) (by omega : P ≠ 0); omega
  · have := Int.emod_lt_of_pos (r - T) hP; omega

/-- The example with two two-map classes (G4): every multiplier is below the Lemma A′ bound 1/9 for centres 0, 1, 2. -/
theorem g4_small : 5 * 9 < 56 ∧ 15 * 9 < 196 ∧ 3 * 9 < 56 := by decide


/-! ### Round 9: which parameters break the formula (README section 13) -/

/-- The concrete break: centres 0, 1, 2 with λ = 19/200, μ = 1/20 about 0, ν = 1/19 about 1, κ = 1/10 about 2.
    The pieces λ∘ν and μ∘κ have the same slope, (19/200)(1/19) = 1/200 = (1/20)(1/10): -/
theorem example_slopes : 19 * 1 * (20 * 10) = 1 * 1 * (200 * 19) := by decide

/-- … and the same intercept: both send z to (18 + z)/200. Scaled by 38000:
    λ(1 + ν(z − 1)) ↦ 190(19 + (z − 1)) and μ(2 + κ(z − 2)) ↦ 190(20 + (z − 2)). -/
theorem example_collision (z : Int) : 190 * (19 + (z - 1)) = 190 * (20 + (z - 2)) := by omega

/-- The break sits inside the separated region: every multiplier is below Lemma A′'s bound 1/9. -/
theorem example_in_box : 19 * 9 < 200 ∧ 1 * 9 < 20 ∧ 1 * 9 < 19 ∧ 1 * 9 < 10 := by decide

/-- Theorem P10's reduction terminates: replacing a factor of a word by a shorter one shortens the word. -/
theorem replace_shortens (a b x y : Nat) (h : x < y) : a + x + b < a + y + b := by omega


/-! ### Round 10: the exact size of a break (README section 14) -/

/-- Every point of the collision surface collides. With κ = k/d, λ = μ(2d − k)/d and ν = k/(2d − k):
    the slopes agree, λν = μκ (cleared of fractions) … -/
theorem surface_slopes (m d k : Int) : (m * (2 * d - k)) * k * d = m * k * (d * (2 * d - k)) := by grind

/-- … and the intercepts agree, λ(1 − ν) = 2μ(1 − κ), i.e. (2d − k)(2d − 2k)·d = 2(d − k)·d·(2d − k). -/
theorem surface_intercept (d k : Int) : (2 * d - k) * (2 * d - 2 * k) * d = 2 * (d - k) * d * (2 * d - k) := by grind

/-- The 3-adic break: μ = 3/80, κ = 3/40, λ = 231/3200, ν = 3/77. Slopes: λν = μκ. -/
theorem threeadic_slopes : 231 * 3 * (80 * 40) = 3 * 3 * (3200 * 77) := by decide

/-- Intercepts: λ(1 − ν) = (231/3200)(74/77) equals 2μ(1 − κ) = 2(3/80)(37/40). -/
theorem threeadic_intercept : 231 * 74 * (80 * 40) = 2 * 3 * 37 * (3200 * 77) := by decide

/-- Each multiplier has 3-adic valuation exactly 1: the numerator has one factor 3, the denominator none. -/
theorem threeadic_valuations :
    231 = 3 * 77 ∧ 77 % 3 ≠ 0 ∧ 3200 % 3 ≠ 0 ∧ 80 % 3 ≠ 0 ∧ 40 % 3 ≠ 0 := by decide

/-- Hypothesis (I) for the 3-adic break: the landing zones [0, 2λ], [1 − ν, 1 + ν], [2 − 2κ, 2] are disjoint:
    2λ < 1 − ν is 462·77 < 74·3200, and 1 + ν < 2 − 2κ is 80·40 < 74·77. -/
theorem threeadic_zones : 462 * 77 < 74 * 3200 ∧ 80 * 40 < 74 * 77 := by decide

/-- Hypothesis (II) for the 3-adic break: (2 − κ)(1 − κ) > 1 is 77·37 > 40·40. -/
theorem threeadic_push : 40 * 40 < 77 * 37 := by decide

/-- Hypothesis (II) holds whenever 0 < κ ≤ 1/3 (κ = k/d): (2d − k)(d − k) > d·d. -/
theorem push_general (d k : Int) (hk : 0 < k) (hd : 3 * k ≤ d) : d * d < (2 * d - k) * (d - k) := by
  have h1 : 0 ≤ d * (d - 3 * k) := Int.mul_nonneg (by omega) (by omega)
  have h2 : 0 < k * k := Int.mul_pos hk hk
  have e : (2 * d - k) * (d - k) = d * d + d * (d - 3 * k) + k * k := by grind
  omega

/-- The p-adic step: if a·v = b·v with v ≥ 1 (the valuations of λ^a and μ^b at the shared prime are equal),
    then a = b. -/
theorem equal_powers (a b v : Nat) (hv : 1 ≤ v) (h : a * v = b * v) : a = b :=
  Nat.eq_of_mul_eq_mul_right (by omega) h

/-- The push only grows with the power: x ≤ x^a for x ≥ 1 and a ≥ 1 (so ρ = (2 − κ)^a ≥ 2 − κ). -/
theorem push_grows (x a : Nat) (ha : 1 ≤ a) : x ≤ x ^ a := Nat.le_self_pow (by omega) x

/-- Pushing a zone clear: if ρ₀·lo > hi (the stretched zone starts beyond the target zone), the same holds for
    every larger stretch ρ ≥ ρ₀ (lo ≥ 0). -/
theorem push_clear (r r0 lo hi : Int) (h0 : r0 ≤ r) (hlo : 0 ≤ lo) (h : hi < r0 * lo) : hi < r * lo := by
  have : r0 * lo ≤ r * lo := Int.mul_le_mul_of_nonneg_right h0 hlo
  omega

/-- The original break has no shared odd prime (both numerators are 1), and the argument first fails at
    a = 13, b = 10: ρ = λ¹³/μ¹⁰ = 19¹³·20¹⁰/200¹³ lies in the window [9/19, 100/171] where ρ·(K-zone) can meet
    the N-zone. -/
theorem original_first_bad :
    9 * 200 ^ 13 ≤ 19 * (19 ^ 13 * 20 ^ 10) ∧ 171 * (19 ^ 13 * 20 ^ 10) ≤ 100 * 200 ^ 13 := by decide

end NewTies
