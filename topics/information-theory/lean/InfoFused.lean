/-
  Information theory in Fused Hyperoperation Algebra
  ==================================================

  Companion page: ../index.html.  Numbers: ../checks.py.  Pictures: ../figures.py.

  STATUS: written for core Lean 4 only (no mathlib), no `sorry`, in the same
  style as the project's compiled files (core `Rat`, `grind` for ring
  arithmetic, `Rat.mul_nonneg` / `Rat.mul_pos` for signs, `omega` for whole
  numbers).  NOT YET COMPILED: Lean is not installed in the workspace where
  this was written.  Every ring identity in this file was checked
  symbolically with sympy and every list lemma numerically (see ../checks.py
  and ../checks_output.txt).  Run `lean InfoFused.lean` to compile it.

  The idea: a small kit of fused lemmas (Part 0) proves a whole family of
  known information-theory laws (Parts 1 to 7).  The kit is

    K1  lift       1 + k·(a ◯[k] b) = (1 + k·a)(1 + k·b)
    K2  climb      1 + k·(a ◯[k] a ◯[k] … ◯[k] a) = (1 + k·a)ⁿ
    K3  flip       (−a) ◯[k] (−b) = −(a ◯[−k] b)
    K4  linearize  any L with L(xy) = L(x) + L(y) turns ◯[k] into +
    K5  outer sum  Σᵢⱼ xᵢ·yⱼ = (Σ xᵢ)(Σ yⱼ)
    K6  outer log  Σᵢⱼ xᵢyⱼ·L(xᵢyⱼ) = (Σy)·Σ x·L(x) + (Σx)·Σ y·L(y)

  where a ◯[k] b = a + b + k·a·b.  The weight k says how the corner cell a·b of
  the mult-add square is counted.

  What Lean does NOT check here: anything that needs real numbers (log, real
  powers p^q, the binary entropy function h(p) and its monotonicity), limits,
  and the probability models themselves (that a channel is a transition matrix
  and a cascade is a matrix product).  Where a law needs log, the proof is
  stated for any L that turns products into sums, which log does.
-/
namespace InfoFused

/-! ## 0. The common kit -/

/-- Weighted mult-add. -/
def wmul (k a b : Rat) : Rat := a + b + k * a * b
/-- Plain mult-add, rung 1. -/
def circ (a b : Rat) : Rat := a + b + a * b

/-- Our own power, so every lemma used is proved here. -/
def pw (x : Rat) : Nat → Rat
  | 0 => 1
  | n + 1 => pw x n * x

/-- The n-fold weighted mult-add of a with itself (the mult-add power). -/
def mpow (k a : Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => wmul k (mpow k a n) a

theorem wmul_comm (k a b : Rat) : wmul k a b = wmul k b a := by unfold wmul; grind
theorem wmul_assoc (k a b c : Rat) : wmul k (wmul k a b) c = wmul k a (wmul k b c) := by
  unfold wmul; grind
theorem wmul_zero (k a : Rat) : wmul k a 0 = a := by unfold wmul; grind
theorem circ_eq (a b : Rat) : circ a b = wmul 1 a b := by unfold circ wmul; grind

/-- K1, lift: the coordinate 1 + k·a multiplies. -/
theorem lift (k a b : Rat) : 1 + k * wmul k a b = (1 + k * a) * (1 + k * b) := by
  unfold wmul; grind

/-- Every nonzero weight is plain mult-add after scaling by k. -/
theorem scale (k a b : Rat) : k * wmul k a b = circ (k * a) (k * b) := by
  unfold wmul circ; grind

/-- K3, flip: reading a and b as their negatives swaps the sign of the weight. -/
theorem flip (k a b : Rat) : wmul k (-a) (-b) = -(wmul (-k) a b) := by unfold wmul; grind

theorem pw_add (x : Rat) (m n : Nat) : pw x (m + n) = pw x m * pw x n := by
  induction n with
  | zero => simp only [Nat.add_zero, pw]; grind
  | succ n ih => rw [Nat.add_succ, pw, pw, ih]; grind

theorem pw_mul_base (x y : Rat) (n : Nat) : pw (x * y) n = pw x n * pw y n := by
  induction n with
  | zero => simp only [pw]; grind
  | succ n ih => rw [pw, pw, pw, ih]; grind

/-- K2, climb: the mult-add power climbs the cube. -/
theorem climb (k a : Rat) (n : Nat) : 1 + k * mpow k a n = pw (1 + k * a) n := by
  induction n with
  | zero => simp only [mpow, pw]; grind
  | succ n ih => rw [mpow, lift, ih, pw]

/-- K4, linearize: any map that turns products into sums turns ◯[k] into +
    in the coordinate 1 + k·a.  Taking L = log gives log(1 + k·a). -/
theorem linearize {G : Type} [Add G] (L : Rat → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (k a b : Rat) : L (1 + k * wmul k a b) = L (1 + k * a) + L (1 + k * b) := by
  rw [lift, hL]

/-! ### Sums over a product of two lists (the joint distribution) -/

def lsum : List Rat → Rat
  | [] => 0
  | x :: xs => x + lsum xs

/-- All products xᵢ·yⱼ, row by row: the cells of the outer-product grid. -/
def outer : List Rat → List Rat → List Rat
  | [], _ => []
  | x :: xs, ys => ys.map (fun y => x * y) ++ outer xs ys

theorem lsum_append (xs ys : List Rat) : lsum (xs ++ ys) = lsum xs + lsum ys := by
  induction xs with
  | nil => simp only [List.nil_append, lsum]; grind
  | cons x xs ih => simp only [List.cons_append, lsum]; rw [ih]; grind

theorem lsum_map (x : Rat) (ys : List Rat) : lsum (ys.map (fun y => x * y)) = x * lsum ys := by
  induction ys with
  | nil => simp only [List.map_nil, lsum]; grind
  | cons y ys ih => simp only [List.map_cons, lsum]; rw [ih]; grind

/-- K5, outer sum: the grid's total is the product of the two edges. -/
theorem outer_sum (xs ys : List Rat) : lsum (outer xs ys) = lsum xs * lsum ys := by
  induction xs with
  | nil => simp only [outer, lsum]; grind
  | cons x xs ih => simp only [outer, lsum]; rw [lsum_append, lsum_map, ih]; grind

/-- Σ x·L(x): Shannon's sum with L in place of log (it is −H when L = log). -/
def shan (L : Rat → Rat) : List Rat → Rat
  | [] => 0
  | x :: xs => x * L x + shan L xs

theorem shan_append (L : Rat → Rat) (xs ys : List Rat) :
    shan L (xs ++ ys) = shan L xs + shan L ys := by
  induction xs with
  | nil => simp only [List.nil_append, shan]; grind
  | cons x xs ih => simp only [List.cons_append, shan]; rw [ih]; grind

theorem shan_map (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y) (x : Rat) (ys : List Rat) :
    shan L (ys.map (fun y => x * y)) = x * lsum ys * L x + x * shan L ys := by
  induction ys with
  | nil => simp only [List.map_nil, shan, lsum]; grind
  | cons y ys ih => simp only [List.map_cons, shan, lsum]; rw [ih, hL]; grind

/-- K6, outer log. -/
theorem outer_shan (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y) (xs ys : List Rat) :
    shan L (outer xs ys) = lsum ys * shan L xs + lsum xs * shan L ys := by
  induction xs with
  | nil => simp only [outer, shan, lsum]; grind
  | cons x xs ih =>
    simp only [outer, shan, lsum]; rw [shan_append, shan_map L hL, ih]; grind

/-! ## 1. Failure chances: the k = −1 rule -/

/-- Union rule for independent events: P(A or B) = p + q − pq = p ◯[−1] q. -/
theorem union_rule (p q : Rat) : wmul (-1) p q = 1 - (1 - p) * (1 - q) := by unfold wmul; grind

/-- Survival multiplies (K1 at k = −1). -/
theorem survival (p q : Rat) : 1 - wmul (-1) p q = (1 - p) * (1 - q) := by
  have h := lift (-1) p q; grind

/-- The k = −1 rule is mult-add read on the negatives (K3). -/
theorem union_is_flipped_circ (p q : Rat) : wmul (-1) p q = -(circ (-p) (-q)) := by
  unfold wmul circ; grind

/-- Inclusion–exclusion for three independent events falls out of associativity:
    the signed corners of the mult-add cube. -/
theorem incl_excl3 (p q r : Rat) :
    wmul (-1) (wmul (-1) p q) r = p + q + r - p*q - p*r - q*r + p*q*r := by
  unfold wmul; grind

/-- A series system of n parts, each failing with chance p, fails with chance 1 − (1 − p)ⁿ (K2). -/
theorem series_system (p : Rat) (n : Nat) : mpow (-1) p n = 1 - pw (1 - p) n := by
  have h := climb (-1) p n
  have e : (1 : Rat) + -1 * p = 1 - p := by grind
  rw [e] at h; grind

/-- Survival multiplies, so cumulative hazard −log(1 − p) adds (K4). -/
theorem hazards_add (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y) (p q : Rat) :
    L (1 - wmul (-1) p q) = L (1 - p) + L (1 - q) := by
  rw [survival, hL]

/-- Adding an independent part never lowers the chance of a failure. -/
theorem union_grows (p q : Rat) (hp : p ≤ 1) (hq : 0 ≤ q) : p ≤ wmul (-1) p q := by
  have h := Rat.mul_nonneg hq (show (0 : Rat) ≤ 1 - p by grind)
  unfold wmul; grind

/-! ## 2. Entropy -/

/-- Hartley: typical-set sizes multiply, so bits add (log₂ of 2ᵐ·2ⁿ is m + n). -/
theorem hartley (m n : Nat) : pw 2 (m + n) = pw 2 m * pw 2 n := pw_add 2 m n

/-- Shannon additivity for independent sources, H(X,Y) = H(X) + H(Y),
    proved for any L with L(xy) = L(x) + L(y) and both distributions summing to 1. -/
theorem shannon_additive (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y)
    (ps rs : List Rat) (hp : lsum ps = 1) (hr : lsum rs = 1) :
    shan L (outer ps rs) = shan L ps + shan L rs := by
  rw [outer_shan L hL, hp, hr]; grind

/-- The joint distribution of independent sources is a distribution (K5). -/
theorem joint_sums_to_one (ps rs : List Rat) (hp : lsum ps = 1) (hr : lsum rs = 1) :
    lsum (outer ps rs) = 1 := by
  rw [outer_sum, hp, hr]; grind

/-- Power sums multiply.  With xᵢ = pᵢ^q and yⱼ = rⱼ^q (and (pr)^q = p^q·r^q for
    real powers), this is Σ (pᵢrⱼ)^q = (Σ pᵢ^q)(Σ rⱼ^q) (K5). -/
theorem power_sums_multiply (ws vs : List Rat) : lsum (outer ws vs) = lsum ws * lsum vs :=
  outer_sum ws vs

/-- Tsallis pseudo-additivity, S_q(A,B) = S_A + S_B + (1 − q)·S_A·S_B.
    With k = 1 − q, Tsallis entropy is S = (P − 1)/k where P is the power sum.
    Stated with the k cleared: if k·S_A = P − 1, k·S_B = R − 1 and
    k·S_AB = P·R − 1, then k·S_AB = k·(S_A ◯[k] S_B). -/
theorem tsallis (k P R SA SB SAB : Rat) (hA : k * SA = P - 1) (hB : k * SB = R - 1)
    (hAB : k * SAB = P * R - 1) : k * SAB = k * wmul k SA SB := by
  have h := lift k SA SB
  rw [hA, hB] at h
  grind

/-- Rényi additivity: Rényi entropy is L(1 + k·S)/k, so it adds where Tsallis
    mult-adds (K1 + K4).  Stated with the 1/k cleared. -/
theorem renyi_additive (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y)
    (k SA SB : Rat) : L (1 + k * wmul k SA SB) = L (1 + k * SA) + L (1 + k * SB) :=
  linearize L hL k SA SB

/-- Landmark q = 2 (k = −1): Gini–Simpson G = 1 − Σp² is the chance two draws differ.
    Collision sums multiply, so G composes by the union rule:
    "the pair differs" = "A differs or B differs". -/
theorem gini_simpson (QA QB : Rat) : 1 - QA * QB = wmul (-1) (1 - QA) (1 - QB) := by
  unfold wmul; grind

/-- Landmark q = 0 (k = 1): S₀ = N − 1, support sizes multiply, so N − 1 mult-adds:
    the N×M grid minus its unit corner. -/
theorem support_count (N M : Rat) : N * M - 1 = circ (N - 1) (M - 1) := by
  unfold circ; grind

/-- Landmark q = 1 (k = 0): Shannon, plain addition. -/
theorem shannon_weight (a b : Rat) : wmul 0 a b = a + b := by unfold wmul; grind

/-- Rung 2: for N equally likely outcomes and k = 1 − q a whole number m,
    1 + k·S = Nᵐ, so k·S is the fused power (N − 1)^◯m (K2 with k = 1). -/
theorem uniform_tsallis_is_fused_power (N : Rat) (m : Nat) : mpow 1 (N - 1) m = pw N m - 1 := by
  have h := climb 1 (N - 1) m
  have e : (1 : Rat) + 1 * (N - 1) = N := by grind
  have e2 : (1 : Rat) * mpow 1 (N - 1) m = mpow 1 (N - 1) m := by grind
  rw [e, e2] at h; grind

/-! ## 3. Compression -/

/-- Kraft sums of two prefix codes multiply when you concatenate codewords (K5),
    and the slack 1 − K composes like failure chances (k = −1). -/
theorem kraft_slack (KA KB : Rat) : 1 - KA * KB = wmul (-1) (1 - KA) (1 - KB) := by
  unfold wmul; grind

/-- The concatenated code still satisfies Kraft. -/
theorem kraft_product (KA KB : Rat) (h0 : 0 ≤ KA) (hA : KA ≤ 1) (hB : KB ≤ 1) : KA * KB ≤ 1 := by
  have h := Rat.mul_nonneg h0 (show (0 : Rat) ≤ 1 - KB by grind)
  grind

/-- Campbell's exponential cost: M = Σ pᵢ·D^(t·ℓᵢ) multiplies for concatenated
    independent sources (K5 on the weights), so M − 1 composes by plain mult-add. -/
theorem campbell (MA MB : Rat) : MA * MB - 1 = circ (MA - 1) (MB - 1) := by
  unfold circ; grind

/-- Expected length adds for concatenated codes (the k = 0 shadow of Campbell):
    Σᵢⱼ pᵢrⱼ(ℓᵢ + mⱼ) = Σ pℓ + Σ rm when both sum to 1.  Stated on the grid sums. -/
theorem lengths_add (EA EB : Rat) : wmul 0 EA EB = EA + EB := shannon_weight EA EB

/-! ## 4. Channels: the alphabet ladder k = −(1 + c) -/

/-- Binary symmetric channel cascade: the off-diagonal entry of the product of
    [[1−p, p],[p, 1−p]] and [[1−q, q],[q, 1−q]] is p ◯[−2] q. -/
theorem bsc_cascade (p q : Rat) : (1 - p) * q + p * (1 - q) = wmul (-2) p q := by
  unfold wmul; grind

theorem bsc_cascade_diag (p q : Rat) : (1 - p) * (1 - q) + p * q = 1 - wmul (-2) p q := by
  unfold wmul; grind

/-- Polarization 1 − 2p multiplies along a cascade (K1 at k = −2). -/
theorem bsc_polarization (p q : Rat) : 1 - 2 * wmul (-2) p q = (1 - 2 * p) * (1 - 2 * q) := by
  have h := lift (-2) p q; grind

/-- Matsui's piling-up lemma: biases ε = 1/2 − p combine as 2·ε₁·ε₂. -/
theorem piling_up (p q : Rat) : 1/2 - wmul (-2) p q = 2 * (1/2 - p) * (1/2 - q) := by
  unfold wmul; grind

/-- n identical flips: 1 − 2·P(odd) = (1 − 2p)ⁿ (K2). -/
theorem n_flips (p : Rat) (n : Nat) : 1 - 2 * mpow (-2) p n = pw (1 - 2 * p) n := by
  have h := climb (-2) p n
  have e : (1 : Rat) + -2 * p = 1 - 2 * p := by grind
  rw [e] at h; grind

/-- The q-ary symmetric channel.  An error moves a symbol to one of the other
    letters at random; a second error lands back on the right letter with chance
    c = 1/(alphabet − 1).  Right after two hops: both right, or wrong then sent
    back.  The error composes with weight k = −(1 + c). -/
theorem qary_cascade (e1 e2 c : Rat) :
    1 - ((1 - e1) * (1 - e2) + c * e1 * e2) = wmul (-(1 + c)) e1 e2 := by
  unfold wmul; grind

/-- c = 0, letters never cancel (huge alphabet, or erasures): the union rule. -/
theorem ladder_never_cancel (e1 e2 : Rat) : wmul (-(1 + 0)) e1 e2 = wmul (-1) e1 e2 := by
  unfold wmul; grind
/-- c = 1, two letters: a second flip always undoes the first. -/
theorem ladder_binary (e1 e2 : Rat) : wmul (-(1 + 1)) e1 e2 = wmul (-2) e1 e2 := by
  unfold wmul; grind
/-- c = 1/3, four letters: the depolarizing qubit channel's −4/3. -/
theorem ladder_four (e1 e2 : Rat) : wmul (-(1 + 1/3)) e1 e2 = wmul (-4/3) e1 e2 := by
  unfold wmul; grind

/-- Binary erasure channel: erasures compose by the union rule and capacity
    1 − e multiplies exactly. -/
theorem bec_capacity (e1 e2 : Rat) : 1 - wmul (-1) e1 e2 = (1 - e1) * (1 - e2) := survival e1 e2

/-! ## 5. Data processing: a cascade is never better than its worst link -/

/-- BEC: capacity of the cascade ≤ capacity of each link. -/
theorem bec_dpi (C1 C2 : Rat) (h1 : 0 ≤ C1) (h2 : C2 ≤ 1) : C1 * C2 ≤ C1 := by
  have h := Rat.mul_nonneg h1 (show (0 : Rat) ≤ 1 - C2 by grind)
  grind

/-- BSC (crossovers taken ≤ 1/2): the cascade's crossover is at least each link's. -/
theorem bsc_dpi (p q : Rat) (hp : p ≤ 1/2) (hq : 0 ≤ q) : p ≤ wmul (-2) p q := by
  have h := Rat.mul_nonneg hq (show (0 : Rat) ≤ 1 - 2 * p by grind)
  unfold wmul; grind

/-- …and still ≤ 1/2, so it stays on the side where capacity 1 − h(p) falls as p
    rises.  (That monotonicity is real analysis and is not checked here.) -/
theorem bsc_stays_below_half (p q : Rat) (hp : p ≤ 1/2) (hq : q ≤ 1/2) :
    wmul (-2) p q ≤ 1/2 := by
  have h := Rat.mul_nonneg (show (0 : Rat) ≤ 1/2 - p by grind) (show (0 : Rat) ≤ 1/2 - q by grind)
  have e := piling_up p q
  grind

/-- Strong data processing: the contraction coefficient (1 − 2p)² multiplies. -/
theorem sdpi_multiplies (p q : Rat) :
    (1 - 2 * wmul (-2) p q) * (1 - 2 * wmul (-2) p q)
      = ((1 - 2 * p) * (1 - 2 * p)) * ((1 - 2 * q) * (1 - 2 * q)) := by
  rw [bsc_polarization]; grind

/-- Parallel channels: distinguishable messages multiply, capacities add (Hartley). -/
theorem parallel_capacity (C1 C2 : Nat) : pw 2 (C1 + C2) = pw 2 C1 * pw 2 C2 := pw_add 2 C1 C2

/-! ## 6. Error-correcting codes -/

/-- Parity of three noisy bits: the odd layers of the mult-add cube (1·3·3·1). -/
theorem parity3 (p : Rat) : mpow (-2) p 3 = 3 * p * (1 - p) * (1 - p) + p * p * p := by
  rw [show mpow (-2) p 3 = wmul (-2) (wmul (-2) (wmul (-2) 0 p) p) p from rfl]
  unfold wmul; grind

/-- Every corner of the cube {0,1}ⁿ: (1 + 1)ⁿ − 1 nonempty corners (K2 at k = 1, a = 1). -/
theorem cube_corners (n : Nat) : mpow 1 1 n = pw 2 n - 1 := by
  have h := uniform_tsallis_is_fused_power 2 n
  have e : (2 : Rat) - 1 = 1 := by grind
  rw [e] at h; exact h

/-- Hamming codes are perfect: with n = m + r and n + 1 = 2ʳ, the 2ᵐ balls of
    radius 1 (a corner plus its n neighbours: layers 0 and 1) tile the n-cube. -/
theorem hamming_perfect (r m : Nat) (h : m + r + 1 = 2 ^ r) :
    2 ^ m * (1 + (m + r)) = 2 ^ (m + r) := by
  rw [show 1 + (m + r) = 2 ^ r by omega, ← Nat.pow_add]

/-! ## 7. Decoding: belief propagation is mult-add in two coordinates

  A bit's belief: chance p of a 1, odds o = (1−p)/p, ρ = 1 − 2p = tanh(LLR/2).
  Check node (parity of two independent bits): ρ multiplies, −2p mult-adds.
  Variable node (two independent looks at one bit): odds multiply, o − 1 mult-adds,
  and in ρ that is Einstein's velocity addition (ρ₁ + ρ₂)/(1 + ρ₁ρ₂).
-/

/-- Check node: −2·(p ◯[−2] q) = (−2p) ◯ (−2q), plain mult-add. -/
theorem check_node (p q : Rat) : -2 * wmul (-2) p q = circ (-2 * p) (-2 * q) := by
  unfold wmul circ; grind

/-- Variable node: posterior odds multiply, so o − 1 mult-adds. -/
theorem variable_node (o1 o2 : Rat) : o1 * o2 - 1 = circ (o1 - 1) (o2 - 1) := by
  unfold circ; grind

/-- Odds product in ρ coordinates (ρ = (o−1)/(o+1)), with denominators cleared:
    numerator and denominator of the velocity sum. -/
theorem odds_to_velocity (o1 o2 : Rat) :
    (o1 - 1) * (o2 + 1) + (o2 - 1) * (o1 + 1) = 2 * (o1 * o2 - 1) ∧
    (o1 + 1) * (o2 + 1) + (o1 - 1) * (o2 - 1) = 2 * (o1 * o2 + 1) := by
  constructor <;> grind

/-- Velocity addition lifts to mult-add on both corners:
    (1 + ab) + (a + b) = (1 + a)(1 + b) and (1 + ab) − (a + b) = (1 − a)(1 − b),
    so (1 + v)/(1 − v), the odds, multiplies. -/
theorem velocity_corners (a b : Rat) :
    (1 + a * b) + (a + b) = (1 + a) * (1 + b) ∧ (1 + a * b) - (a + b) = (1 - a) * (1 - b) := by
  constructor <;> grind

/-- Velocity addition is associative, written on (numerator, denominator) pairs. -/
def vadd (x y : Rat × Rat) : Rat × Rat := (x.1 * y.2 + y.1 * x.2, x.2 * y.2 + x.1 * y.1)

theorem vadd_assoc (x y z : Rat × Rat) : vadd (vadd x y) z = vadd x (vadd y z) := by
  simp only [vadd, Prod.mk.injEq]; constructor <;> grind

end InfoFused
