/-
  OsgoodInverse: when the Osgood test can be turned around
  =========================================================

  The exact statements behind the page "The inverse of the Osgood test"
  (Nov 25, 2026). Written for core Lean 4.37.0-rc1 in the style of the shared
  library (FusedAlgebra), no mathlib, no `sorry`.

  LIBRARY NOTE: copied from osgood-inverse/lean/OsgoodInverse.lean (Nov 25).
  It compiled against the library on its first build with no changes beyond
  this header note: no sorry, standard axioms only.

  The Osgood integral is a clock: in the coordinate G(y) = ∫ dy/f(y) the
  solution of y' = f(y) moves at unit speed, so the integral is the time to
  reach the end. Real analysis (integrals, limits) is not in core Lean, so the
  statements here are the exact algebraic and discrete forms of that picture:

  1. The Riccati climb y' = (1+y)²: its solution t/(1−t) obeys the exact
     difference law Δy = h·(1+y(t))(1+y(t+h)) = h·(1 + y(t) ⊞ y(t+h)), and it
     passes every height before t = 1.
  2. The discrete climb with that law hits the wall: its clock 1/(1+y) falls by
     exactly h per step, so it reaches 0 (y = ∞) after 1/h steps. The rung-1
     climb y ↦ y ⊞ h never does. Through the mult-add inverse −u/(1+u) the
     blow-up is a walk at constant speed into the wall −1: blow-up at infinity
     and leaving the wall are one question.
  3. Non-uniqueness: the trapezoid scheme of x' = 2√x has a whole fan of
     solutions out of 0 (one for every start time c); the same scheme for the
     Lipschitz x' = x has only the zero solution.
  4. The test moves down one rung exactly: Cauchy condensation. The blocks
     [2ʲ, 2ʲ⁺¹) start at 1 + 1^(⊞j), have length 2ʲ (the lift factor), and
     the dyadic partial sums are squeezed by the condensed sums both ways, so
     the test for a passes iff the test for j ↦ 2ʲ·a(2ʲ) passes.
     Instances: 1/k (rung 0) grows at least like n/2; 1/k² stays below 2.
  5. Direction decides: the rotation field |x|²·Jx and the radial field
     |x|²·x have the same speed r³, but the first is perpendicular to x (the
     size never changes) and the second pushes straight out (r² obeys the
     Riccati law). A size bound alone cannot force blow-up.

  Not here (needs real analysis): the integrals themselves, the blow-up time
  formula as a limit, the continuous non-uniqueness proof, the rung-height
  growth of fractional iterates, the Navier–Stokes statements.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Pi
import FusedAlgebra.Topics.NSContinuum

namespace MultAdd.OsgoodInverse

open MultAdd Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1–2. The Riccati climb y' = (1+y)² and the wall it hits -/

section Field
variable {K : Type} [Field K]

/-- The clock of the Riccati climb: reciprocals of the lift. -/
theorem recip_diff (u v : K) (hu : u ≠ 0) (hv : v ≠ 0) :
    1 / v - 1 / u = (u - v) / (u * v) := by
  grind

/-- The solution of y' = (1+y)², y(0) = 0. -/
def ric (t : K) : K := t / (1 - t)

/-- Under the lift it is a reciprocal: 1 + y = 1/(1 − t). The time left, 1 − t,
    is the clock. -/
theorem ric_lift (t : K) (h : 1 - t ≠ 0) : 1 + ric t = 1 / (1 - t) := by
  unfold ric; grind

/-- **The exact difference law.** y(t+h) − y(t) = h·(1+y(t))(1+y(t+h)):
    the trapezoid-free discrete form of y' = (1+y)², with no error term. -/
theorem ric_step (t h : K) (h1 : 1 - t ≠ 0) (h2 : 1 - (t + h) ≠ 0) :
    ric (t + h) - ric t = h * (1 + ric t) * (1 + ric (t + h)) := by
  unfold ric; grind

/-- The same law in mult-add: the step is h times the lift of y(t) ⊞ y(t+h). -/
theorem ric_step_multAdd (t h : K) (h1 : 1 - t ≠ 0) (h2 : 1 - (t + h) ≠ 0) :
    ric (t + h) - ric t = h * (1 + multAdd (ric t) (ric (t + h))) := by
  rw [ric_step t h h1 h2, multAdd_lift]; grind

/-- Read in the clock coordinate, the climb is uniform: one step of the law
    lowers 1/(1+y) by exactly h. -/
theorem riccati_clock_step (a b h : K) (ha : 1 + a ≠ 0) (hb : 1 + b ≠ 0)
    (hs : b - a = h * (1 + a) * (1 + b)) : 1 / (1 + b) = 1 / (1 + a) - h := by
  grind

/-- n·h without casts. -/
def steps (h : K) : Nat → K
  | 0 => 0
  | n + 1 => steps h n + h

theorem steps_add (h : K) (m n : Nat) : steps h (m + n) = steps h m + steps h n := by
  induction n with
  | zero => show steps h m = steps h m + 0; grind
  | succ n ih => show steps h (m + n) + h = steps h m + (steps h n + h); rw [ih]; grind

/-- **The discrete Riccati climb runs on a straight clock.** As long as the lift
    stays off the wall, 1/(1+yₙ) = 1/(1+y₀) − n·h. -/
theorem riccati_clock (y : Nat → K) (h : K)
    (hs : ∀ n, y (n + 1) - y n = h * (1 + y n) * (1 + y (n + 1))) :
    ∀ n, (∀ k, k ≤ n → 1 + y k ≠ 0) → 1 / (1 + y n) = 1 / (1 + y 0) - steps h n := by
  intro n
  induction n with
  | zero => intro _; show 1 / (1 + y 0) = 1 / (1 + y 0) - 0; grind
  | succ n ih =>
    intro hne
    have h1 := ih (fun k hk => hne k (by omega))
    have h2 := riccati_clock_step (y n) (y (n + 1)) h (hne n (by omega)) (hne (n + 1) (by omega)) (hs n)
    show 1 / (1 + y (n + 1)) = 1 / (1 + y 0) - (steps h n + h)
    rw [h2, h1]; grind

/-- **Blow-up.** If the clock starts at m·h, the climb cannot take m steps with
    finite values: at step m the clock reads 0, and no finite y has
    1/(1+y) = 0. (Continuous form: the blow-up time is ∫₀^∞ ds/(1+s)² = 1.) -/
theorem riccati_hits_wall (y : Nat → K) (h : K) (m : Nat)
    (hs : ∀ n, y (n + 1) - y n = h * (1 + y n) * (1 + y (n + 1)))
    (hstart : 1 / (1 + y 0) = steps h m) :
    ¬ (∀ k, k ≤ m → 1 + y k ≠ 0) := by
  intro hne
  have hc := riccati_clock y h hs m hne
  rw [hstart] at hc
  have hm := hne m (Nat.le_refl m)
  have hz : 1 / (1 + y m) = 0 := by rw [hc]; grind
  have hone : (1 + y m) * (1 / (1 + y m)) = 1 := by grind
  rw [hz] at hone
  grind

/-- **The mirror.** Through the mult-add inverse v = −u/(1+u), a Riccati step
    is a step of exactly −h: the climb to infinity is a walk at constant speed
    into the wall −1. Blow-up at infinity and leaving the wall are one question,
    swapped by the mult-add inverse. -/
theorem mirror_step (a b h : K) (ha : 1 + a ≠ 0) (hb : 1 + b ≠ 0)
    (hs : b - a = h * (1 + a) * (1 + b)) : multAddInv b - multAddInv a = -h := by
  have hc := riccati_clock_step a b h ha hb hs
  unfold multAddInv
  grind

/-- The mirrored walk reaches the wall exactly when the climb blows up:
    1 + (−u/(1+u)) = 1/(1+u). -/
theorem mirror_lift (u : K) (hu : 1 + u ≠ 0) : 1 + multAddInv u = 1 / (1 + u) := by
  unfold multAddInv; grind

/-- **The rung-1 climb never reaches the wall.** If yₙ₊₁ = yₙ ⊞ h then
    1 + yₙ = (1 + y₀)(1 + h)ⁿ: the clock log₁₊ y advances by log(1+h) per step
    forever (continuous form: ∫ ds/(1+s) = log(1+s) → ∞). -/
theorem rung1_lift (y : Nat → K) (h : K) (hs : ∀ n, y (n + 1) = multAdd (y n) h) :
    ∀ n, 1 + y n = (1 + y 0) * (1 + h) ^ n := by
  intro n
  induction n with
  | zero => rw [Semiring.pow_zero]; grind
  | succ n ih => rw [Semiring.pow_succ, hs n, multAdd_lift, ih]; grind

end Field

/-- On the rationals the rung-1 climb stays on the safe side of the wall forever. -/
theorem rung1_never_wall (y : Nat → Rat) (h : Rat) (hh : 0 < 1 + h) (h0 : 0 < 1 + y 0)
    (hs : ∀ n, y (n + 1) = multAdd (y n) h) : ∀ n, 0 < 1 + y n := by
  intro n
  induction n with
  | zero => exact h0
  | succ n ih =>
    have e : 1 + y (n + 1) = (1 + y n) * (1 + h) := by rw [hs n, multAdd_lift]
    rw [e]
    have := Rat.mul_lt_mul_of_pos_left hh ih
    grind

/-- With h = 1/4 the Riccati clock started at 1 reads 0 after 4 steps:
    the climb from y₀ = 0 blows up at step 4 = 1/h. -/
theorem steps_quarter : steps (1 / 4 : Rat) 4 = 1 := by
  decide +kernel

theorem riccati_quarter_blows_up (y : Nat → Rat) (h0 : y 0 = 0)
    (hs : ∀ n, y (n + 1) - y n = (1 / 4) * (1 + y n) * (1 + y (n + 1))) :
    ¬ (∀ k, k ≤ 4 → 1 + y k ≠ 0) :=
  riccati_hits_wall y (1 / 4) 4 hs (by rw [h0, steps_quarter] <;> decide +kernel)

/-- The continuous climb passes every height before t = 1:
    y(1 − 1/(M+2)) = M + 1. -/
theorem ric_near_one (M : Nat) : ric (1 - 1 / ((M : Rat) + 2)) = (M : Rat) + 1 := by
  have h0 : (0 : Rat) ≤ (M : Rat) := by exact_mod_cast Nat.zero_le M
  have hM : (M : Rat) + 2 ≠ 0 := by grind
  have hwq : 1 / ((M : Rat) + 2) * ((M : Rat) + 2) = 1 := by grind
  have hw0 : 1 / ((M : Rat) + 2) ≠ 0 := by intro h; rw [h] at hwq; grind
  unfold ric
  have e : (1 : Rat) - (1 - 1 / ((M : Rat) + 2)) = 1 / ((M : Rat) + 2) := by grind
  rw [e]
  grind

/-- … and that time is still before 1: the time left 1/(M+2) is positive. -/
theorem time_left_pos (M : Nat) : 0 < 1 / ((M : Rat) + 2) := by
  have h0 : (0 : Rat) ≤ (M : Rat) := by exact_mod_cast Nat.zero_le M
  have hp : (0 : Rat) ≤ (M : Rat) + 2 := by grind
  have hM : (M : Rat) + 2 ≠ 0 := by grind
  apply Classical.byContradiction
  intro hn
  have hle : 1 / ((M : Rat) + 2) ≤ 0 := by grind
  have hmul := Rat.mul_le_mul_of_nonneg_left hle hp
  have e : ((M : Rat) + 2) * (1 / ((M : Rat) + 2)) = 1 := by grind
  grind

/-! ## 3. Non-uniqueness out of 0: x' = 2√x versus x' = x

  The trapezoid scheme of x' = 2√x on the integer grid (step 1) is
  x(n+1) = x(n) + √x(n) + √x(n+1). For every start time c the curve
  x_c(n) = (n − c)² (zero until c, Nat subtraction) solves it with √x_c = n − c:
  the solution may sit at 0 as long as it likes and then leave. This is the
  discrete face of ∫₀ dr/(2√r) = √r being finite. -/

/-- The fan member that leaves 0 at time c. -/
def fan (c n : Nat) : Nat := (n - c) ^ 2

/-- Its square root. -/
def root (c n : Nat) : Nat := n - c

theorem root_sq (c n : Nat) : root c n ^ 2 = fan c n := rfl

/-- **Every fan member solves the trapezoid scheme of x' = 2√x.** -/
theorem fan_step (c n : Nat) : fan c (n + 1) = fan c n + root c n + root c (n + 1) := by
  unfold fan root
  by_cases h : c ≤ n
  · obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le h
    have e1 : c + m + 1 - c = m + 1 := by omega
    have e2 : c + m - c = m := by omega
    rw [e1, e2] <;> grind
  · have e1 : n + 1 - c = 0 := by omega
    have e2 : n - c = 0 := by omega
    rw [e1, e2] <;> simp

/-- Every fan member starts at 0. -/
theorem fan_zero (c : Nat) : fan c 0 = 0 := by
  simp [fan]

/-- **Different start times give different solutions**, so the scheme has
    infinitely many solutions from x(0) = 0. -/
theorem fan_distinct (c d : Nat) (h : c < d) : fan c ≠ fan d := by
  intro heq
  have hd := congrFun heq d
  unfold fan at hd
  rw [Nat.sub_self] at hd
  obtain ⟨k, hk⟩ : ∃ k, d - c = k + 1 := ⟨d - c - 1, by omega⟩
  rw [hk] at hd
  grind

/-- **Lipschitz means unique.** The trapezoid scheme of x' = x,
    x(n+1) − x(n) = (h/2)(x(n) + x(n+1)), started at 0 stays at 0 (for h ≠ 2). -/
theorem lipschitz_unique (x : Nat → Rat) (h : Rat) (hh : 2 - h ≠ 0) (h0 : x 0 = 0)
    (hs : ∀ n, x (n + 1) - x n = h / 2 * (x n + x (n + 1))) : ∀ n, x n = 0 := by
  intro n
  induction n with
  | zero => exact h0
  | succ n ih =>
    have e := hs n
    rw [ih] at e
    have hp : x (n + 1) * (2 - h) = 0 := by grind
    have hq : x (n + 1) = x (n + 1) * (2 - h) / (2 - h) := by grind
    rw [hq, hp]; grind

/-! ## 4. The test moves down one rung exactly (Cauchy condensation)

  The substitution s = E(x) (1 + s = eˣ) turns ∫ ds/((1+s)g(s)) into
  ∫ dx/g(E(x)): the test for f = (1+s)·g is the test for g∘E, one rung down.
  The discrete form is Cauchy condensation: sample at the points
  1 + 1^(⊞j) = 2ʲ, weight each sample by the block length 2ʲ (the lift
  factor). -/

/-- Σ_{i<len} a(s + i). -/
def fsum (a : Nat → Rat) (s : Nat) : Nat → Rat
  | 0 => 0
  | len + 1 => fsum a s len + a (s + len)

theorem fsum_append (a : Nat → Rat) (s m : Nat) :
    ∀ n, fsum a s (m + n) = fsum a s m + fsum a (s + m) n
  | 0 => by show fsum a s m = fsum a s m + 0; grind
  | n + 1 => by
    show fsum a s (m + n) + a (s + (m + n)) = fsum a s m + (fsum a (s + m) n + a (s + m + n))
    rw [fsum_append a s m n]
    have e : s + m + n = s + (m + n) := by omega
    rw [e]; grind

theorem fsum_mono (a b : Nat → Rat) (s : Nat) :
    ∀ len, (∀ i, i < len → a (s + i) ≤ b (s + i)) → fsum a s len ≤ fsum b s len
  | 0 => fun _ => by show (0 : Rat) ≤ 0; grind
  | len + 1 => fun h => by
    have ih := fsum_mono a b s len (fun i hi => h i (by omega))
    have hl := h len (by omega)
    show fsum a s len + a (s + len) ≤ fsum b s len + b (s + len)
    grind

theorem fsum_nonneg (a : Nat → Rat) (ha : ∀ k, 0 ≤ a k) (s : Nat) : ∀ len, 0 ≤ fsum a s len
  | 0 => by show (0 : Rat) ≤ 0; grind
  | len + 1 => by
    have ih := fsum_nonneg a ha s len
    have := ha (s + len)
    show 0 ≤ fsum a s len + a (s + len)
    grind

theorem one_le_pow2 : ∀ j : Nat, 1 ≤ 2 ^ j
  | 0 => by decide
  | j + 1 => by have := one_le_pow2 j; rw [Nat.pow_succ]; omega

theorem pow2_succ (j : Nat) : 2 ^ (j + 1) = 2 ^ j + 2 ^ j := by
  rw [Nat.pow_succ]; omega

/-- A constant over a block of length 2ʲ sums to 2ʲ·c. -/
theorem const_block (c : Rat) (j : Nat) : ∀ s, fsum (fun _ => c) s (2 ^ j) = (2 : Rat) ^ j * c := by
  induction j with
  | zero => intro s; show (0 : Rat) + c = (2 : Rat) ^ 0 * c; grind
  | succ j ih =>
    intro s
    rw [pow2_succ, fsum_append, ih s, ih (s + 2 ^ j), Rat.pow_succ]
    grind

/-- **The blocks are the mult-add lattice.** The j-th block starts at
    2ʲ = 1 + 1^(⊞j), the lift of the j-fold mult-add of 1. -/
theorem block_start (j : Nat) : 2 ^ j = mpow (1 : Nat) j + 1 :=
  (MultAdd.NSContinuum.mpow_one_nat j).symm

/-- The j-th block [2ʲ, 2ʲ⁺¹). -/
def blk (a : Nat → Rat) (j : Nat) : Rat := fsum a (2 ^ j) (2 ^ j)

/-- The terms only shrink from index 1 on. -/
def Antitone1 (a : Nat → Rat) : Prop := ∀ i k, 1 ≤ i → i ≤ k → a k ≤ a i

theorem blk_le (a : Nat → Rat) (ha : Antitone1 a) (j : Nat) :
    blk a j ≤ (2 : Rat) ^ j * a (2 ^ j) := by
  unfold blk
  rw [← const_block (a (2 ^ j)) j (2 ^ j)]
  apply fsum_mono
  intro i _
  exact ha (2 ^ j) (2 ^ j + i) (one_le_pow2 j) (by omega)

theorem le_blk (a : Nat → Rat) (ha : Antitone1 a) (j : Nat) :
    (2 : Rat) ^ j * a (2 ^ (j + 1)) ≤ blk a j := by
  unfold blk
  rw [← const_block (a (2 ^ (j + 1))) j (2 ^ j)]
  apply fsum_mono
  intro i hi
  have h1 := one_le_pow2 j
  have h2 := pow2_succ j
  exact ha (2 ^ j + i) (2 ^ (j + 1)) (by omega) (by omega)

/-- Dyadic partial sums: a(1) + … + a(2ⁿ − 1). -/
def dsum (a : Nat → Rat) (n : Nat) : Rat := fsum a 1 (2 ^ n - 1)

theorem dsum_succ (a : Nat → Rat) (n : Nat) : dsum a (n + 1) = dsum a n + blk a n := by
  unfold dsum blk
  have h1 := one_le_pow2 n
  have e : 2 ^ (n + 1) - 1 = (2 ^ n - 1) + 2 ^ n := by rw [pow2_succ]; omega
  rw [e, fsum_append]
  have e2 : 1 + (2 ^ n - 1) = 2 ^ n := by omega
  rw [e2]

/-- Condensed partial sums: Σ_{j<n} 2ʲ·a(2ʲ). -/
def csum (a : Nat → Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => csum a n + (2 : Rat) ^ n * a (2 ^ n)

/-- **Upper squeeze.** a(1) + … + a(2ⁿ − 1) ≤ Σ_{j<n} 2ʲ·a(2ʲ). -/
theorem dsum_le_csum (a : Nat → Rat) (ha : Antitone1 a) : ∀ n, dsum a n ≤ csum a n
  | 0 => by show fsum a 1 (2 ^ 0 - 1) ≤ 0; show (0 : Rat) ≤ 0; grind
  | n + 1 => by
    rw [dsum_succ]
    have ih := dsum_le_csum a ha n
    have hb := blk_le a ha n
    show dsum a n + blk a n ≤ csum a n + (2 : Rat) ^ n * a (2 ^ n)
    grind

/-- **Lower squeeze.** Σ_{j≤n} 2ʲ·a(2ʲ) ≤ a(1) + 2·(a(1) + … + a(2ⁿ − 1)). -/
theorem csum_le_dsum (a : Nat → Rat) (ha : Antitone1 a) :
    ∀ n, csum a (n + 1) ≤ a 1 + 2 * dsum a n
  | 0 => by
    show 0 + (2 : Rat) ^ 0 * a (2 ^ 0) ≤ a 1 + 2 * fsum a 1 (2 ^ 0 - 1)
    show 0 + (2 : Rat) ^ 0 * a 1 ≤ a 1 + 2 * 0
    grind
  | n + 1 => by
    have ih := csum_le_dsum a ha n
    have hb := le_blk a ha n
    rw [dsum_succ]
    show csum a (n + 1) + (2 : Rat) ^ (n + 1) * a (2 ^ (n + 1)) ≤ a 1 + 2 * (dsum a n + blk a n)
    rw [Rat.pow_succ]
    grind

/-- **The test is the same one rung down.** For shrinking nonnegative terms,
    the dyadic partial sums are bounded iff the condensed sums are. -/
theorem condensation (a : Nat → Rat) (ha : Antitone1 a) (hpos : ∀ k, 0 ≤ a k) :
    (∃ B, ∀ n, csum a n ≤ B) ↔ (∃ B, ∀ n, dsum a n ≤ B) := by
  constructor
  · intro ⟨B, hB⟩
    exact ⟨B, fun n => Rat.le_trans (dsum_le_csum a ha n) (hB n)⟩
  · intro ⟨B, hB⟩
    refine ⟨a 1 + 2 * B, fun n => ?_⟩
    have hB0 : 0 ≤ B := by
      have h := hB 0
      have : dsum a 0 = 0 := rfl
      grind
    have h1 := hpos 1
    cases n with
    | zero => show (0 : Rat) ≤ a 1 + 2 * B; grind
    | succ n =>
      have := csum_le_dsum a ha n
      have := hB n
      grind

theorem lt_pow2 : ∀ n : Nat, n < 2 ^ n
  | 0 => by decide
  | n + 1 => by have := lt_pow2 n; have := one_le_pow2 n; rw [Nat.pow_succ]; omega

/-- Every partial sum a(1) + … + a(N) sits under a dyadic one, so the dyadic
    sums decide boundedness of all partial sums. -/
theorem fsum_le_dsum (a : Nat → Rat) (hpos : ∀ k, 0 ≤ a k) (N : Nat) :
    fsum a 1 N ≤ dsum a N := by
  unfold dsum
  have h := lt_pow2 N
  have e : 2 ^ N - 1 = N + (2 ^ N - 1 - N) := by omega
  rw [e, fsum_append]
  have := fsum_nonneg a hpos (1 + N) (2 ^ N - 1 - N)
  grind

/-! ### Instances: rung 0 passes, one power past it fails -/

/-- The naturals inside the rationals, without casts. -/
def nr : Nat → Rat
  | 0 => 0
  | k + 1 => nr k + 1

theorem nr_add (m : Nat) : ∀ n, nr (m + n) = nr m + nr n
  | 0 => by show nr m = nr m + 0; grind
  | n + 1 => by show nr (m + n) + 1 = nr m + (nr n + 1); rw [nr_add m n]; grind

theorem nr_nonneg : ∀ k, 0 ≤ nr k
  | 0 => by show (0 : Rat) ≤ 0; grind
  | k + 1 => by have := nr_nonneg k; show 0 ≤ nr k + 1; grind

theorem nr_pos (k : Nat) (hk : 1 ≤ k) : 0 < nr k := by
  obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
  have := nr_nonneg m
  show 0 < nr m + 1
  grind

theorem nr_mono (i k : Nat) (h : i ≤ k) : nr i ≤ nr k := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le h
  rw [nr_add]
  have := nr_nonneg d
  grind

theorem nr_pow2 : ∀ j, nr (2 ^ j) = (2 : Rat) ^ j
  | 0 => by show nr 1 = (2 : Rat) ^ 0; show (0 : Rat) + 1 = (2 : Rat) ^ 0; grind
  | j + 1 => by rw [pow2_succ, nr_add, nr_pow2 j, Rat.pow_succ]; grind

theorem one_div_pos (p : Rat) (hp : 0 < p) : 0 < 1 / p := by
  apply Classical.byContradiction
  intro hn
  have hle : 1 / p ≤ 0 := by grind
  have hmul := Rat.mul_le_mul_of_nonneg_left hle (by grind : (0 : Rat) ≤ p)
  have e : p * (1 / p) = 1 := by grind
  grind

/-- The harmonic terms 1/k: the discrete ∫ ds/s, rung 0. -/
def harm (k : Nat) : Rat := 1 / nr k

theorem harm_nonneg (k : Nat) : 0 ≤ harm k := by
  unfold harm
  rcases Nat.eq_zero_or_pos k with h | h
  · subst h; show (0 : Rat) ≤ 1 / 0; decide +kernel
  · have := one_div_pos (nr k) (nr_pos k h)
    grind

theorem harm_anti : Antitone1 harm := by
  intro i k hi hik
  unfold harm
  exact PiMultAdd.div_le_div_of (nr_pos i hi) (nr_mono i k hik)
    (by decide +kernel : (0 : Rat) ≤ 1) (Rat.le_refl)

/-- One rung down, the harmonic terms become all ones: 2ʲ · (1/2ʲ) = 1. -/
theorem csum_harm : ∀ n, csum harm n = nr n
  | 0 => rfl
  | n + 1 => by
    show csum harm n + (2 : Rat) ^ n * (1 / nr (2 ^ n)) = nr n + 1
    rw [csum_harm n, nr_pow2]
    have hp : (2 : Rat) ^ n ≠ 0 := by
      have := nr_pos (2 ^ n) (one_le_pow2 n); rw [nr_pow2] at this; grind
    grind

/-- **Rung 0 passes.** 1 + 1/2 + … + 1/(2ⁿ − 1) ≥ n/2: the clock never stops. -/
theorem harmonic_grows (n : Nat) : nr n ≤ 2 * dsum harm n := by
  have h := csum_le_dsum harm harm_anti n
  rw [csum_harm] at h
  have h1 : harm 1 = 1 := by show 1 / nr 1 = 1; show 1 / ((0 : Rat) + 1) = 1; decide +kernel
  show nr n ≤ 2 * dsum harm n
  have : nr (n + 1) = nr n + 1 := rfl
  grind

/-- The inverse squares 1/k²: one full power past rung 0. -/
def isq (k : Nat) : Rat := 1 / (nr k * nr k)

theorem isq_anti : Antitone1 isq := by
  intro i k hi hik
  unfold isq
  have hp := nr_pos i hi
  have hm := nr_mono i k hik
  have hpk : 0 ≤ nr k := nr_nonneg k
  have s1 := Rat.mul_le_mul_of_nonneg_left hm (by grind : (0 : Rat) ≤ nr i)
  have s2 := Rat.mul_le_mul_of_nonneg_left hm hpk
  have hpp : 0 < nr i * nr i := by
    have := Rat.mul_lt_mul_of_pos_left hp hp; grind
  have hle : nr i * nr i ≤ nr k * nr k := by grind
  exact PiMultAdd.div_le_div_of hpp hle (by decide +kernel : (0 : Rat) ≤ 1) (Rat.le_refl)

/-- One rung down, the inverse squares become a geometric series:
    Σ_{j<n} 2ʲ/4ʲ = 2 − 2/2ⁿ. -/
theorem csum_isq : ∀ n, csum isq n = 2 - 2 * (1 / (2 : Rat) ^ n)
  | 0 => by show (0 : Rat) = 2 - 2 * (1 / (2 : Rat) ^ 0); grind
  | n + 1 => by
    show csum isq n + (2 : Rat) ^ n * (1 / (nr (2 ^ n) * nr (2 ^ n))) = 2 - 2 * (1 / (2 : Rat) ^ (n + 1))
    rw [csum_isq n, nr_pow2, Rat.pow_succ]
    have hp : (2 : Rat) ^ n ≠ 0 := by
      have := nr_pos (2 ^ n) (one_le_pow2 n); rw [nr_pow2] at this; grind
    grind

/-- **One power past rung 0 fails the test (the sum converges).**
    1 + 1/4 + 1/9 + … + 1/(2ⁿ − 1)² ≤ 2 for every n. -/
theorem inverse_squares_bounded (n : Nat) : dsum isq n ≤ 2 := by
  have h := dsum_le_csum isq isq_anti n
  rw [csum_isq] at h
  have hp : 0 < (2 : Rat) ^ n := by
    have := nr_pos (2 ^ n) (one_le_pow2 n); rw [nr_pow2] at this; exact this
  have := one_div_pos _ hp
  grind

/-! ## 5. Direction decides: same speed, opposite outcomes -/

/-- The rotation field |x|²·Jx in the plane. -/
def rot1 (x y : Int) : Int := (x ^ 2 + y ^ 2) * (-y)
def rot2 (x y : Int) : Int := (x ^ 2 + y ^ 2) * x

/-- The radial field |x|²·x. -/
def rad1 (x y : Int) : Int := (x ^ 2 + y ^ 2) * x
def rad2 (x y : Int) : Int := (x ^ 2 + y ^ 2) * y

/-- Both have speed r³ (speed² = r⁶): the size alone fails the Osgood test. -/
theorem rot_speed (x y : Int) : rot1 x y ^ 2 + rot2 x y ^ 2 = (x ^ 2 + y ^ 2) ^ 3 := by
  unfold rot1 rot2; grind

theorem rad_speed (x y : Int) : rad1 x y ^ 2 + rad2 x y ^ 2 = (x ^ 2 + y ^ 2) ^ 3 := by
  unfold rad1 rad2; grind

/-- **The rotation pushes sideways:** x·F = 0, so d/dt |x|² = 2x·F = 0 and the
    size never changes, though the speed fails the Osgood test. -/
theorem rot_perp (x y : Int) : x * rot1 x y + y * rot2 x y = 0 := by
  unfold rot1 rot2; grind

/-- **The radial field pushes straight out:** x·F = r⁴, so ρ = r² obeys
    ρ' = 2ρ², the Riccati law of §1, and blows up at time 1/(2ρ₀). -/
theorem rad_out (x y : Int) : x * rad1 x y + y * rad2 x y = (x ^ 2 + y ^ 2) ^ 2 := by
  unfold rad1 rad2; grind

end MultAdd.OsgoodInverse
