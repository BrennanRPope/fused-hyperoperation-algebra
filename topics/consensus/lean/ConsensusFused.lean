/-
  Consensus (alternatives to proof of work) in Fused Hyperoperation Algebra
  =========================================================================

  Companion page: ../index.html.  Numbers: ../checks.py.  Pictures: ../figures.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`. Checked with Lean 4.24.0:
      lean ConsensusFused.lean
  Uses core `Rat`, `grind` for ring arithmetic, `omega` for whole numbers and
  `decide +kernel` for exact rational checks.

  Weighted mult-add:  a ⊕ₖ b = a + b + k·a·b.
    k =  1   plain mult-add, 1 + (a ⊕ b) = (1+a)(1+b)
    k = −1   "either one happens":  a ⊕₋₁ b = 1 − (1−a)(1−b)

  What Lean does NOT check: anything with real numbers, logarithms or limits;
  the probability models themselves (that a race is a gambler's ruin, that a
  committee is a binomial draw); and the simulations. Those are in checks.py.
-/
namespace ConsensusFused

def wmul (k a b : Rat) : Rat := a + b + k * a * b
def madd (a b : Rat) : Rat := a + b + a * b

/-- Our own power, so every lemma we use is proved here. -/
def pw (x : Rat) : Nat → Rat
  | 0 => 1
  | n + 1 => pw x n * x

/-- n-fold weighted mult-add of a with itself (the mult-add power). -/
def mpow (k a : Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => wmul k (mpow k a n) a

/-! ## 0. The algebra -/

theorem madd_lift (a b : Rat) : 1 + madd a b = (1 + a) * (1 + b) := by unfold madd; grind
theorem madd_comm (a b : Rat) : madd a b = madd b a := by unfold madd; grind
theorem madd_assoc (a b c : Rat) : madd (madd a b) c = madd a (madd b c) := by unfold madd; grind
theorem madd_zero (a : Rat) : madd a 0 = a := by unfold madd; grind
theorem wmul_lift (k a b : Rat) : 1 + k * wmul k a b = (1 + k * a) * (1 + k * b) := by
  unfold wmul; grind
/-- k = −1 is "at least one of two independent events". -/
theorem either (a b : Rat) : wmul (-1) a b = 1 - (1 - a) * (1 - b) := by unfold wmul; grind

theorem pw_add (x : Rat) (m n : Nat) : pw x (m + n) = pw x m * pw x n := by
  induction n with
  | zero => simp [pw]
  | succ n ih => rw [Nat.add_succ, pw, pw, ih]; grind

theorem pw_mul (x : Rat) (m n : Nat) : pw x (m * n) = pw (pw x m) n := by
  induction n with
  | zero => simp [pw]
  | succ n ih => rw [Nat.mul_succ, pw_add, ih, pw]

theorem pw_one (n : Nat) : pw 1 n = 1 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [pw, ih]

/-- The mult-add power climbs the cube: 1 + k·(a ⊕ₖ … ⊕ₖ a) = (1 + k·a)^n. -/
theorem mpow_lift (k a : Rat) (n : Nat) : 1 + k * mpow k a n = pw (1 + k * a) n := by
  induction n with
  | zero => simp only [mpow, pw]; grind
  | succ n ih => rw [mpow, wmul_lift, ih, pw]

theorem mpow_lift1 (a : Rat) (n : Nat) : 1 + mpow 1 a n = pw (1 + a) n := by
  have h := mpow_lift 1 a n
  simp only [Rat.one_mul] at h; exact h

/-! ## 1. Proof of work: depth is a mult-add power, 1/2 is the identity

  An attacker with share q, honest share 1 − q, z blocks behind, catches up
  with probability r^z where r = q/(1−q).  Put a = r − 1.  Then r^z = 1 + a^(⊕z).
-/

/-- The step a = r − 1 written without division. -/
theorem step_formula (q : Rat) (h : q ≠ 1) : (q / (1 - q) - 1) * (1 - q) = 2 * q - 1 := by
  have h1 : (1 - q) ≠ 0 := by grind
  have : q / (1 - q) * (1 - q) = q := Rat.div_mul_cancel h1
  grind

/-- At the 51% line the step is exactly the mult-add identity. -/
theorem half_is_identity : (1/2 : Rat) / (1 - 1/2) - 1 = 0 := by decide +kernel

/-- Powers of the identity do nothing: at q = 1/2 waiting more blocks buys no security. -/
theorem depth_useless (n : Nat) : mpow 1 0 n = 0 := by
  have h := mpow_lift1 0 n
  rw [show (1 : Rat) + 0 = 1 by decide +kernel, pw_one] at h
  grind

/-- Catch-up probability as a mult-add power, checked exactly for q = 1/4, z = 6:
    r = 1/3, a = −2/3, and 1 + a^(⊕6) = (1/3)^6. -/
theorem quarter_six : 1 + mpow 1 (-2/3) 6 = 1 / 729 := by decide +kernel

/-! ## 2. BFT voting: where 1/3 comes from

  Two quorums of size Q among n validators overlap in at least 2Q − n.
  Safety: the overlap must hold an honest validator, 2Q − n > f.
  Liveness: a quorum must form without the f faulty ones, Q = n − f.
-/

theorem bft_safe (n f : Nat) (h : 3 * f < n) : 2 * (n - f) > n + f := by omega
theorem bft_unsafe (n f : Nat) (h : n ≤ 3 * f) : 2 * (n - f) ≤ n + f := by omega
/-- Three quorums sharing an honest validator need f < n/4. -/
theorem three_quorums (n f : Nat) (h : 4 * f < n) : 3 * (n - f) > 2 * n + f := by omega
/-- The longest-chain rule needs only a majority: one "quorum" against the rest. -/
theorem one_quorum (n f : Nat) (h : 2 * f < n) : n - f > f := by omega

/-! ## 3. Committees and delegation: the layers of the cube

  Three independent members, each faulty with chance p: the cube's layers hold
  (1−p)^3, 3p(1−p)^2, 3p^2(1−p), p^3 (Pascal 1,3,3,1).  A majority of three
  goes bad on the top two layers.
-/

def maj3 (p : Rat) : Rat := p^3 + 3 * p^2 * (1 - p)

theorem layers_sum (p : Rat) :
    (1 - p)^3 + 3 * p * (1 - p)^2 + 3 * p^2 * (1 - p) + p^3 = 1 := by grind
theorem maj3_closed (p : Rat) : maj3 p = 3 * p^2 - 2 * p^3 := by unfold maj3; grind
theorem maj3_half : maj3 (1/2) = 1/2 := by decide +kernel
theorem maj3_mirror (p : Rat) : maj3 (1 - p) = 1 - maj3 p := by unfold maj3; grind
/-- The pull toward 0 below 1/2 and toward 1 above it. -/
theorem maj3_drift (p : Rat) : maj3 p - p = -(p * (1 - p) * (1 - 2 * p)) := by unfold maj3; grind
theorem maj3_shrinks (p : Rat) (h0 : 0 < p) (h1 : p < 1/2) : maj3 p < p := by
  have hp : 0 < p * (1 - p) := Rat.mul_pos h0 (by grind)
  have hq : 0 < p * (1 - p) * (1 - 2 * p) := Rat.mul_pos hp (by grind)
  have := maj3_drift p
  grind

/-! ## 4. Candidate B fails: the Sybil corner

  Weighting stake by log(1+s) makes weight additive under mult-add.  Real stake
  adds with +.  In the multiplicative coordinate W = 1 + s, splitting s1 + s2
  into two accounts turns 1 + s1 + s2 into (1 + s1)(1 + s2): the gain is the
  corner square s1·s2, never negative.
-/

theorem sybil_corner (s1 s2 : Rat) : (1 + s1) * (1 + s2) - (1 + (s1 + s2)) = s1 * s2 := by grind
theorem sybil_gain (s1 s2 : Rat) (h1 : 0 ≤ s1) (h2 : 0 ≤ s2) :
    1 + (s1 + s2) ≤ (1 + s1) * (1 + s2) := by
  have := Rat.mul_nonneg h1 h2; grind

/-- An additive weight is split-proof: n equal accounts weigh what one does. -/
theorem additive_split (f : Rat → Rat) (hf : ∀ a b, f (a + b) = f a + f b) (s : Rat) :
    ∀ n : Nat, f ((n : Rat) * s) = (n : Rat) * f s := by
  have h0 : f 0 = 0 := by have := hf 0 0; simp at this; grind
  intro n
  induction n with
  | zero => simp [h0]
  | succ n ih =>
    have : ((n + 1 : Nat) : Rat) * s = (n : Rat) * s + s := by push_cast; grind
    rw [this, hf, ih]; push_cast; grind

/-! ## 5. Candidate A fails: fusing stake with age

  Two parties.  The attacker holds 1/5 of the stake and has been silent for 10
  epochs; the honest side holds 4/5 and last produced 1 epoch ago.
-/

/-- Coin age s·(1+t): the 20% attacker outweighs the 80% honest side. -/
theorem coin_age_flip : (4/5 : Rat) * (1 + 1) < (1/5) * (1 + 10) := by decide +kernel
/-- Fused exponent (1+s)^(1+t) − 1: the same flip, more strongly. -/
theorem fused_age_flip : pw (1 + 4/5) (1 + 1) - 1 < pw (1 + 1/5) (1 + 10) - 1 := by
  decide +kernel
/-- Plain stake keeps the honest majority. -/
theorem linear_keeps : (1/5 : Rat) < 4/5 := by decide +kernel

/-! ## 6. Slashing and long-range attacks: k = −1 powers -/

/-- m independent watchers each catching a double vote with chance d:
    caught = d ⊕₋₁ d ⊕₋₁ … = 1 − (1 − d)^m.  Same for the share of an old
    validator set that has left after t epochs at churn c. -/
theorem caught (d : Rat) (m : Nat) : mpow (-1) d m = 1 - pw (1 - d) m := by
  have h := mpow_lift (-1) d m
  have e : (1 : Rat) + -1 * d = 1 - d := by grind
  rw [e] at h; grind

/-! ## 7. Verifiable delay functions: doubling up the cube

  A VDF computes x^(2^T) by T squarings.  In mult-add coordinates squaring is
  a ⊕ a.  Pietrzak's proof checks the two halves: (x^(2^m))^(2^n) = x^(2^(m+n)).
-/

def dbl : Nat → Rat → Rat
  | 0, a => a
  | n + 1, a => madd (dbl n a) (dbl n a)

theorem dbl_lift (a : Rat) (n : Nat) : 1 + dbl n a = pw (1 + a) (2 ^ n) := by
  induction n with
  | zero => simp [dbl, pw]
  | succ n ih =>
    rw [dbl, madd_lift, ih, Nat.pow_succ, pw_mul]
    simp [pw]

theorem pietrzak_halves (x : Rat) (m n : Nat) : pw (pw x (2 ^ m)) (2 ^ n) = pw x (2 ^ (m + n)) := by
  rw [Nat.pow_add, pw_mul]

/-! ## 8. Energy, from published figures (TWh per year) -/

/-- Ethereum before the Merge ~23 TWh, after ~0.0026 TWh (CCRI 2022): over 99.98% less. -/
theorem merge_cut : (23 - 26/10000 : Rat) / 23 > 9998/10000 := by decide +kernel

end ConsensusFused
