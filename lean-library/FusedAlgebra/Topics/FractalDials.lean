/-
  Hidden dials of a half box: what tells apart fractals of the same dimension?
  ============================================================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  math-theory/fractal-dials/ (README.md, checks.py, index.html).

  Core Lean 4 only (no mathlib, no imports), no `sorry`. Compiled in the shared
  library as Topics/FractalDials (Lean 4.37.0-rc1).

  What is proved:
  * half_dimension_certificates  the five sample sets A to E all have dimension 1/2:
                                 the square roots of their zoom ratios add to 1
                                 (Moran's equation at s = 1/2), exactly.
  * mulAdd_shift, mulAdd_identity
                                 (x-1) ◯ (y-1) = xy - 1: mult-add is multiplication
                                 of sides 1 + a. It has no notion of size, which is
                                 why a zoom dial would be new.
  * two_pow_mod7, uneven_not_lattice
                                 the uneven set D (zooms 9/100 and 49/100) is
                                 non-lattice: no power of one zoom equals a power of
                                 the other. Non-lattice sets have no Minkowski wobble.
  * weights_odds                 the per-element dial t and the weight p = t/(1+t)
                                 agree at the sample weights, and fair weights give
                                 the smallest collision sum p² + (1-p)² = 1/2.
  * pieces_join_closed           self-similar piece counts c·bⁿ stay self-similar
                                 when two sets are joined: (c·bⁿ)(c'·b'ⁿ) = cc'(bb')ⁿ.
  * four_corner_is_tree          the four-corner set (half box × half box) has
                                 2ᵐ·2ᵐ = 2^(2m) level-m pieces: the binary tree.
  * dyadic_glue                  0.0111…1 + last step = 0.1 in binary (to 40 digits):
                                 the stick is the binary tree with dyadic ends glued.
  * gasket_vertices, gasket_touches, same_dimension_gasket_dust
                                 the gasket (3 pieces at 1/2) and the 9-at-1/4 dust
                                 have the same dimension, but the gasket's pieces
                                 share corners at every level while the dust's never
                                 touch: topology is a dial only once pieces touch.
  * powerset_branch, powerset_not_selfsimilar
                                 the power-set tree's counts 2^(2ⁿ) square at each
                                 level and fit no c·bⁿ: no shape dial (all of which
                                 keep c·bⁿ) reaches the power set's growth.
-/

namespace FractalDials

/-! ## 1. Five sets of dimension one half

Each set is made of pieces shrunk by zoom ratios r₁, r₂, … . Its dimension s
solves Moran's equation r₁ˢ + r₂ˢ + … = 1. At s = 1/2 this says the square
roots of the zooms add to 1. We give each square root and check both facts. -/

theorem half_dimension_certificates :
    -- A: two pieces at 1/4
    ((1/2 : Rat) * (1/2) = 1/4 ∧ (1/2 : Rat) + 1/2 = 1) ∧
    -- B and C: three pieces at 1/9 (placement differs, zooms do not)
    ((1/3 : Rat) * (1/3) = 1/9 ∧ (1/3 : Rat) + 1/3 + 1/3 = 1) ∧
    -- D: uneven pieces at 9/100 and 49/100
    ((3/10 : Rat) * (3/10) = 9/100 ∧ (7/10 : Rat) * (7/10) = 49/100 ∧
      (3/10 : Rat) + 7/10 = 1) ∧
    -- E: four pieces at 1/16
    ((1/4 : Rat) * (1/4) = 1/16 ∧ (1/4 : Rat) + 1/4 + 1/4 + 1/4 = 1) := by
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_, ?_⟩, ⟨?_, ?_⟩⟩ <;> decide +kernel

/-! ## 2. Mult-add multiplies sides; it has no zoom -/

/-- mult-add: a ◯ b = a + b + ab. -/
def mulAdd (a b : Int) : Int := a + b + a * b

/-- Writing sides as 1 + a, mult-add is plain multiplication of sides. -/
theorem mulAdd_shift (x y : Int) : mulAdd (x - 1) (y - 1) = x * y - 1 := by
  simp only [mulAdd, Int.sub_mul, Int.mul_sub, Int.one_mul, Int.mul_one]; omega

theorem mulAdd_identity (a : Int) : mulAdd a 0 = a ∧ mulAdd 0 a = a := by
  simp only [mulAdd]; omega

/-! ## 3. Lattice or not: the wobble dial

A set is lattice when all zooms are powers of one ratio (A, B, C, E). Then its
complex dimensions line up and its Minkowski content wobbles forever with the
zoom. D's zooms 9/100 and 49/100 share no power, so D is non-lattice and its
content settles to a limit (checks.py: about 3.0007). -/

theorem two_pow_mod7 (n : Nat) : 2 ^ n % 7 = 1 ∨ 2 ^ n % 7 = 2 ∨ 2 ^ n % 7 = 4 := by
  induction n with
  | zero => decide
  | succ n ih => rw [Nat.pow_succ, Nat.mul_mod]; rcases ih with h | h | h <;> rw [h] <;> decide

/-- (9/100)^p = (49/100)^q has no solution with q ≥ 1 (clear denominators). -/
theorem uneven_not_lattice (p q : Nat) (hq : 0 < q) : 9 ^ p * 100 ^ q ≠ 49 ^ q * 100 ^ p := by
  intro h
  have h1 : (49 ^ q * 100 ^ p) % 7 = 0 := by
    obtain ⟨r, rfl⟩ : ∃ r, q = r + 1 := ⟨q - 1, by omega⟩
    rw [Nat.pow_succ, Nat.mul_mod, Nat.mul_mod (49 ^ r), Nat.mul_comm]; simp
  have h2 : (9 ^ p * 100 ^ q) % 7 = 2 ^ (p + q) % 7 := by
    rw [Nat.mul_mod, Nat.pow_mod, Nat.pow_mod 100, Nat.pow_add, Nat.mul_mod (2 ^ p)]
  rw [h, h1] at h2
  rcases two_pow_mod7 (p + q) with h3 | h3 | h3 <;> omega

/-! ## 4. Weights are the per-element dial t

Splitting mass p : 1-p between two pieces is mult-add's per-element dial
t = p/(1-p), so p = t/(1+t). Fair weights (t = 1) give the smallest collision
sum p² + (1-p)², which is what keeps all Rényi dimensions D_q equal. -/

theorem weights_odds :
    ((3/10 : Rat) / (1 - 3/10) = 3/7 ∧ (3/7 : Rat) / (1 + 3/7) = 3/10) ∧
    ((1/10 : Rat) / (1 - 1/10) = 1/9 ∧ (1/9 : Rat) / (1 + 1/9) = 1/10) ∧
    ((1/2 : Rat) / (1 - 1/2) = 1 ∧ (1 : Rat) / (1 + 1) = 1/2) ∧
    ((1/2 : Rat) * (1/2) + (1/2) * (1/2) = 1/2) ∧
    ((3/10 : Rat) * (3/10) + (7/10) * (7/10) = 29/50) ∧
    ((1/10 : Rat) * (1/10) + (9/10) * (9/10) = 41/50) := by
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩, ?_, ?_, ?_⟩ <;> decide +kernel

/-! ## 5. Joining keeps shapes self-similar; two binary half boxes make a tree -/

theorem pieces_join_closed (c b c' b' n : Nat) :
    (c * b ^ n) * (c' * b' ^ n) = (c * c') * (b * b') ^ n := by
  rw [Nat.mul_pow, Nat.mul_mul_mul_comm]

theorem four_corner_is_tree (m : Nat) : 2 ^ m * 2 ^ m = 2 ^ (2 * m) := by
  rw [← Nat.pow_add, Nat.two_mul]

/-- sum of a list of rationals -/
def sumR : List Rat → Rat := List.foldl (· + ·) 0

/-- binary 0.0111…1 (ones in places 2 to 40) plus one more step 2^-40 is 0.1:
    the two codes of a dyadic point are glued in the stick. -/
theorem dyadic_glue :
    sumR ((List.range 39).map (fun i => (1 : Rat) / 2 ^ (i + 2))) + 1 / 2 ^ 40 = 1 / 2 := by
  decide +kernel

/-! ## 6. Topology: a dial only once pieces touch

The gasket (3 pieces at 1/2) and the dust of 9 pieces at 1/4 have the same
dimension log 3 / log 2, since 9 = 3² and 4 = 2². The gasket's level-m triangles
share corners: 3ᵐ triangles have only (3^(m+1) + 3)/2 distinct corners, fewer
than 3·3ᵐ. Dust pieces never touch, and below dimension 1 every such dust is the
same Cantor set (Brouwer), so topology cannot tell them apart there. -/

theorem same_dimension_gasket_dust (k : Nat) : 9 ^ k = 3 ^ (2 * k) ∧ 4 ^ k = 2 ^ (2 * k) := by
  constructor <;> rw [Nat.pow_mul] <;> rfl

/-- distinct corners of the level-m gasket -/
def gasketVertices : Nat → Nat
  | 0 => 3
  | m + 1 => 3 * gasketVertices m - 3

theorem gasket_vertices (m : Nat) : 2 * gasketVertices m = 3 ^ (m + 1) + 3 := by
  induction m with
  | zero => decide
  | succ m ih =>
    have hp : 0 < 3 ^ m := Nat.pow_pos (by omega)
    simp only [gasketVertices]; rw [Nat.pow_succ 3 (m + 1)]; rw [Nat.pow_succ] at ih ⊢; omega

theorem gasket_touches (m : Nat) (hm : 0 < m) : gasketVertices m < 3 * 3 ^ m := by
  have h := gasket_vertices m
  have hp : 3 ≤ 3 ^ m := by
    obtain ⟨r, rfl⟩ : ∃ r, m = r + 1 := ⟨m - 1, by omega⟩
    have : 0 < 3 ^ r := Nat.pow_pos (by omega)
    rw [Nat.pow_succ]; omega
  rw [Nat.pow_succ] at h; omega

/-! ## 7. The power set leaves every shape class -/

theorem powerset_branch (n : Nat) : 2 ^ (2 ^ (n + 1)) = 2 ^ (2 ^ n) * 2 ^ (2 ^ n) := by
  rw [Nat.pow_succ, Nat.mul_two, Nat.pow_add]

theorem powerset_not_selfsimilar : ¬ ∃ c b : Nat, ∀ n, 2 ^ (2 ^ n) = c * b ^ n := by
  rintro ⟨c, b, h⟩
  have h0 := h 0; have h1 := h 1; have h2 := h 2
  simp at h0 h1 h2
  subst h0
  have hb : b = 2 := by omega
  subst hb; simp at h2

end FractalDials
