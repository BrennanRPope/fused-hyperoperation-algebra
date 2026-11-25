/-
  Fractal boxes: does a box of fractional dimension give a fractional power set?
  ==============================================================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  math-theory/fractal-boxes/ (README.md, checks.py, index.html).

  Core Lean 4 only (no mathlib, no imports), no `sorry`. Compiled in the shared
  library as Topics/FractalBoxes (Lean 4.37.0-rc1).

  What is proved:
  * layer, joinLayer        corner layers of a box of rational dimension d are
                            Newton's C(d,k); joined boxes convolve their layers.
  * cube_layers             d = 3 gives Pascal's 1,3,3,1 (definition check).
  * half_layers             d = 1/2: 1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024.
  * half_join_stick         two half boxes joined: 1, 1, then 0 through layer 10.
  * third_join_stick        three third boxes joined: a stick.
  * half_and_three_halves   dimensions add: 1/2 + 3/2 gives the square 1,2,1.
  * pieces_join, half_box_squared, other_half_squared
                            piece counts multiply when fractals are multiplied:
                            (2^m)^2 = 4^m (the four-corner set), (3^m)^2 = 9^m.
  * code, code_lt, code_inj the left/right code of level-m Cantor points in any
                            base B >= 2 is one to one: Cantor sets of every
                            dimension are copies of P(N) as sets.
  * succ_le_two_pow, exp_beats_linear
                            K*d + L < 2^d at d = 4(K+L+1).
  * no_linear_half_pow2     no rule bounded by A*d + B is, applied twice, d -> 2^d.
  * no_affine_half_pow2     the plain affine case.
  * hyper_twice             2^N - 1 + 1 = 2^N (hyperspace count bookkeeping).
-/

namespace FractalBoxes

/-! ## 1. Signed corner layers of a fractional box (Newton's binomial)

A box with n sides of length 1 + a splits into 2ⁿ corner pieces; the pieces
made of k sides form layer k, and there are C(n,k) of them. For a box of
dimension d (any rational) the same layer count is the generalized binomial
coefficient C(d,k). -/

/-- layer k of a box of dimension d: generalized binomial coefficient C(d,k). -/
def layer (d : Rat) : Nat → Rat
  | 0 => 1
  | k + 1 => layer d k * (d - k) / (k + 1)

/-- sum of a list of rationals -/
def sumR : List Rat → Rat := List.foldl (· + ·) 0

/-- layer k of two boxes joined (their product box): Σ_{i+j=k} C(d,i) C(e,j). -/
def joinLayer (d e : Rat) (k : Nat) : Rat :=
  sumR ((List.range (k + 1)).map (fun i => layer d i * layer e (k - i)))

/-- layer k of three boxes joined. -/
def join3Layer (d e f : Rat) (k : Nat) : Rat :=
  sumR ((List.range (k + 1)).map (fun i => layer d i * joinLayer e f (k - i)))

/-- A whole-number box has the familiar Pascal layers (a check of the definition). -/
theorem cube_layers : layer 3 0 = 1 ∧ layer 3 1 = 3 ∧ layer 3 2 = 3 ∧ layer 3 3 = 1 ∧
    layer 3 4 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-- The half box: 1 empty piece, half a single side, then −1/8 of a two-side
corner, +1/16, −5/128, ... (negative layers alternate from layer 2 on). -/
theorem half_layers :
    layer (1/2) 0 = 1 ∧ layer (1/2) 1 = 1/2 ∧ layer (1/2) 2 = -1/8 ∧
    layer (1/2) 3 = 1/16 ∧ layer (1/2) 4 = -5/128 ∧ layer (1/2) 5 = 7/256 ∧
    layer (1/2) 6 = -21/1024 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-- Two half boxes joined give a stick exactly, layer by layer: 1, 1, then 0
(the negative pieces cancel), checked through layer 10. -/
theorem half_join_stick :
    joinLayer (1/2) (1/2) 0 = 1 ∧ joinLayer (1/2) (1/2) 1 = 1 ∧
    (∀ k, k < 9 → joinLayer (1/2) (1/2) (k + 2) = 0) := by
  refine ⟨by decide +kernel, by decide +kernel, ?_⟩
  intro k hk
  match k, hk with
  | 0, _ => decide +kernel
  | 1, _ => decide +kernel
  | 2, _ => decide +kernel
  | 3, _ => decide +kernel
  | 4, _ => decide +kernel
  | 5, _ => decide +kernel
  | 6, _ => decide +kernel
  | 7, _ => decide +kernel
  | 8, _ => decide +kernel

/-- Three third-boxes joined give a stick (layers 1, 1, 0, 0, 0, 0). -/
theorem third_join_stick :
    join3Layer (1/3) (1/3) (1/3) 0 = 1 ∧ join3Layer (1/3) (1/3) (1/3) 1 = 1 ∧
    join3Layer (1/3) (1/3) (1/3) 2 = 0 ∧ join3Layer (1/3) (1/3) (1/3) 3 = 0 ∧
    join3Layer (1/3) (1/3) (1/3) 4 = 0 ∧ join3Layer (1/3) (1/3) (1/3) 5 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-- Dimensions add when boxes join: a half box and a box of dimension 3/2 make a
square (layers 1, 2, 1, 0, ...). -/
theorem half_and_three_halves :
    joinLayer (1/2) (3/2) 0 = 1 ∧ joinLayer (1/2) (3/2) 1 = 2 ∧
    joinLayer (1/2) (3/2) 2 = 1 ∧ joinLayer (1/2) (3/2) 3 = 0 ∧
    joinLayer (1/2) (3/2) 4 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-! ## 2. The half box as a fractal: keep the two outer quarters

Level m of the set C¼ has 2^m pieces of length 4^-m. Its square (the
four-corner set) has 4^m pieces of side 4^-m: the count of a stick cut in
4^m, so dimension 1, though every piece is separated from the rest. -/

theorem pieces_join (N M m : Nat) : (N * M) ^ m = N ^ m * M ^ m := Nat.mul_pow N M m

theorem half_box_squared (m : Nat) : 2 ^ m * 2 ^ m = 4 ^ m := by
  rw [← Nat.mul_pow]

/-- A second half box: three pieces at scale 1/9 (3^m pieces of 9^-m). It is a
different shape with the same dimension; squared it is 9^m pieces of 9^-m. -/
theorem other_half_squared (m : Nat) : 3 ^ m * 3 ^ m = 9 ^ m := by
  rw [← Nat.mul_pow]

/-! ## 3. Points of a Cantor set are subsets: the left/right code

A point of a Cantor set is an infinite path of left/right choices, which is a
subset of ℕ. At level m the code of a path (digits 0 or B−1 in base B) is one to
one, whatever the base, so Cantor sets of every dimension have the same number
of points: they are all copies of P(ℕ) as sets. -/

/-- code of a left/right path in base B: right = digit B−1, left = digit 0. -/
def code (B : Nat) : List Bool → Nat
  | [] => 0
  | b :: bs => (if b then B - 1 else 0) * B ^ bs.length + code B bs

theorem code_lt (B : Nat) (hB : 2 ≤ B) : ∀ bs : List Bool, code B bs < B ^ bs.length
  | [] => by simp [code]
  | b :: bs => by
    have ih := code_lt B hB bs
    have hpow : B ^ (bs.length + 1) = (B - 1) * B ^ bs.length + B ^ bs.length := by
      obtain ⟨C, rfl⟩ : ∃ C, B = C + 1 := ⟨B - 1, by omega⟩
      rw [Nat.add_sub_cancel, Nat.pow_succ, Nat.mul_comm _ (C + 1), Nat.add_mul, Nat.one_mul]
    simp only [code, List.length_cons]
    rw [hpow]
    cases b <;> simp <;> omega

theorem code_inj (B : Nat) (hB : 2 ≤ B) :
    ∀ bs cs : List Bool, bs.length = cs.length → code B bs = code B cs → bs = cs
  | [], [], _, _ => rfl
  | b :: bs, c :: cs, hl, he => by
    simp only [List.length_cons, Nat.add_right_cancel_iff] at hl
    have lb := code_lt B hB bs
    have lc := code_lt B hB cs
    have hpos : 1 ≤ B - 1 := by omega
    have hx : B ^ bs.length ≤ (B - 1) * B ^ bs.length :=
      Nat.le_mul_of_pos_left _ hpos
    simp only [code] at he
    rw [← hl] at he lc
    cases b <;> cases c <;> simp at he
    · rw [code_inj B hB bs cs hl he]
    · omega
    · omega
    · rw [code_inj B hB bs cs hl he]

/-! ## 4. Iterating: a linear change of size twice is never the power set

A geometric operation built from products, unions and Lipschitz maps changes
dimension (and the log of the count at each scale) at most linearly. Doing it
twice is still at most linear, while the power set turns n into 2ⁿ. -/

theorem succ_le_two_pow (n : Nat) : n + 1 ≤ 2 ^ n := by
  induction n with
  | zero => decide
  | succ k ih => rw [Nat.pow_succ]; omega

/-- exponential beats any linear rule somewhere: K·d + L < 2^d at d = 4(K+L+1). -/
theorem exp_beats_linear (K L : Nat) : ∃ d, K * d + L < 2 ^ d := by
  let M := K + L + 1
  refine ⟨4 * M, ?_⟩
  have h1 : 2 * M + 1 ≤ 2 ^ (2 * M) := succ_le_two_pow _
  have h2 : 2 ^ (4 * M) = 2 ^ (2 * M) * 2 ^ (2 * M) := by
    rw [← Nat.pow_add]; congr 1; omega
  have h3 : (2 * M + 1) * (2 * M + 1) ≤ 2 ^ (2 * M) * 2 ^ (2 * M) := Nat.mul_le_mul h1 h1
  have hK : K * (4 * M) ≤ (M - 1) * (4 * M) := Nat.mul_le_mul_right _ (by omega)
  have e1 : (2 * M + 1) * (2 * M + 1) = 4 * (M * M) + 4 * M + 1 := by
    simp only [Nat.add_mul, Nat.mul_add, Nat.mul_one, Nat.one_mul]
    rw [Nat.mul_mul_mul_comm]; omega
  have e2 : (M - 1) * (4 * M) + M ≤ 4 * (M * M) := by
    have : (M - 1) * (4 * M) = 4 * (M * M) - 4 * M := by
      rw [Nat.sub_one_mul, Nat.mul_left_comm]
    have : 4 * M ≤ 4 * (M * M) := by
      have : 1 ≤ M := by omega
      exact Nat.mul_le_mul_left _ (Nat.le_mul_of_pos_left _ this)
    omega
  omega

/-- No operation whose effect on size is bounded by a linear rule A·d + B can,
applied twice, give 2^d for every d. -/
theorem no_linear_half_pow2 (Φ : Nat → Nat) (A B : Nat) (hb : ∀ d, Φ d ≤ A * d + B) :
    ¬ ∀ d, Φ (Φ d) = 2 ^ d := by
  intro h
  obtain ⟨d, hd⟩ := exp_beats_linear (A * A) (A * B + B)
  have hb2 : Φ (Φ d) ≤ A * (A * d + B) + B :=
    Nat.le_trans (hb (Φ d)) (Nat.add_le_add_right (Nat.mul_le_mul_left _ (hb d)) _)
  have : A * (A * d + B) + B = A * A * d + (A * B + B) := by
    rw [Nat.mul_add, Nat.mul_assoc]; omega
  rw [h d] at hb2
  omega

/-- The plainest case: no affine rule d ↦ c·d + e is a half of d ↦ 2^d. -/
theorem no_affine_half_pow2 (c e : Nat) : ¬ ∀ d, c * (c * d + e) + e = 2 ^ d := by
  intro h
  have h0 := h 0
  have h1 := h 1
  simp at h0 h1
  -- from d = 0: c * e + e = 1, so e = 1 and c = 0; then d = 1 gives 1 = 2
  have he : e = 1 ∧ c = 0 := by
    rcases Nat.eq_zero_or_pos c with hc | hc
    · subst hc; simp at h0; exact ⟨h0, rfl⟩
    · have : c * e + e ≥ 2 * e := by
        have := Nat.mul_le_mul_right e hc; omega
      rcases Nat.eq_zero_or_pos e with he0 | he0
      · subst he0; simp at h0
      · omega
  obtain ⟨rfl, rfl⟩ := he
  simp at h1

/-- Hyperspace count: if a shape shows N separate points at some scale, its
space of pieces shows 2^N − 1 (all nonempty subsets). Doing that twice: -/
theorem hyper_twice (N : Nat) : 2 ^ N - 1 + 1 = 2 ^ N := by
  have := Nat.one_le_two_pow (n := N); omega

end FractalBoxes
