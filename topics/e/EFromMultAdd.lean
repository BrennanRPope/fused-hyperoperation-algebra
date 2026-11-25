/-
  e from the circle operation (called mult-add before 1.0)
  ===============

  Part of Fused Hyperoperation Algebra. Companion page: ../e/index.html.
  Pictures: figures.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`, only Lean's standard
  axioms (propext, Classical.choice, Quot.sound). Run it with
  `lean EFromMultAdd.lean`. Checked with Lean 4.37.0-rc1.

  The circle operation is a ◯ b = a + b + ab = (1+a)(1+b) − 1. Combine 1/n with itself
  n times by the circle operation and add back the unit corner: that is the staircase
  stair n = (1 + 1/n)ⁿ. On the n-cube, every switch weighs 1/n, so the
  corners with k switches on weigh C(n,k)/nᵏ in total.

  Proved:
  * mpow_eq, mpow_add: n copies of x under the circle operation make (1+x)ⁿ − 1, and the
                  circle power obeys the exponential law.
  * binomial:     (1 + x)ⁿ = Σₖ C(n,k)·xᵏ (the cube's corner pieces fill the box).
  * layer_mul_fact, layer_fact_le: C(n,k)·k! counts ordered picks of k
                  different axes, so C(n,k)·k! ≤ nᵏ.
  * collisions:   the ordered picks that repeat an axis number at most
                  C(k,2)·nᵏ⁻¹.
  * wt_le, wt_ge: layer k weighs between 1/k! − C(k,2)/(n·k!) and 1/k!.
  * stair_lt_succ: the staircase climbs: (1+1/n)ⁿ < (1+1/(n+1))ⁿ⁺¹
                  (from Bernoulli's inequality, proved here too).
  * stair_le_s, s_lt_three: (1+1/n)ⁿ ≤ 1/0! + … + 1/n! < 3.
  * u_step, stair_le_u: u N = s N + 1/(N!·N) only comes down, and every
                  step of the staircase sits under every u N.
  * s_le_stair:   for n ≥ N, the staircase is within 3·C(N,2)/n of s N.
  * kernel checks: s 10 and u 10 trap e between 2.7182818 and 2.7182819.

  What Lean does NOT do here: core Lean has no real numbers or limits.
  The last step, "so the staircase and the sums close on one number, e",
  is the completeness of the real line. Everything before it is proved.

  All of this is classical (Bernoulli, Euler, Cauchy). What is new is only
  the formal check in core Lean and the reading on the circle-operation cube.
-/

namespace CircleOp.E

/-- `layer n k`: how many corners of the n-cube have exactly k switches on
    (grade k). Built by Pascal's rule: adding one more axis to the cube
    copies it twice, once with the new switch off and once with it on. -/
def layer : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 1
  | n + 1, k + 1 => layer n k + layer n (k + 1)

theorem layer_zero_right (n : Nat) : layer n 0 = 1 := by cases n <;> rfl

theorem layer_succ_succ (n k : Nat) : layer (n + 1) (k + 1) = layer n k + layer n (k + 1) := rfl

theorem layer_eq_zero : ∀ n k, n < k → layer n k = 0
  | 0, 0, h => absurd h (by decide)
  | 0, _ + 1, _ => rfl
  | _ + 1, 0, h => absurd h (by omega)
  | n + 1, k + 1, h => by
    rw [layer_succ_succ, layer_eq_zero n k (by omega), layer_eq_zero n (k + 1) (by omega)]

theorem layer_self : ∀ n, layer n n = 1
  | 0 => rfl
  | n + 1 => by rw [layer_succ_succ, layer_self n, layer_eq_zero n (n + 1) (by omega)]

/-- Absorption: (k+1)·C(n+1,k+1) = (n+1)·C(n,k). -/
theorem absorb : ∀ n k, (k + 1) * layer (n + 1) (k + 1) = (n + 1) * layer n k := by
  intro n
  induction n with
  | zero =>
    intro k
    cases k with
    | zero => rfl
    | succ k => simp [layer]
  | succ n ih =>
    intro k
    cases k with
    | zero =>
      have h := ih 0
      rw [layer_succ_succ]
      simp only [layer_zero_right] at h ⊢
      omega
    | succ j =>
      have h1 := ih j
      have h2 := ih (j + 1)
      have hp := layer_succ_succ n j
      rw [layer_succ_succ (n + 1) (j + 1)]
      generalize layer (n + 1) (j + 1) = X at *
      generalize layer (n + 1) (j + 1 + 1) = Y at *
      generalize layer n j = A at *
      generalize layer n (j + 1) = B at *
      subst hp
      grind

/-- Factorial. -/
def fact : Nat → Nat
  | 0 => 1
  | k + 1 => (k + 1) * fact k

theorem fact_pos : ∀ k, 0 < fact k
  | 0 => by decide
  | k + 1 => Nat.mul_pos (by omega) (fact_pos k)

/-- Distinct picks: n·(n-1)·…·(n-k+1), the ordered ways to pick k different axes. -/
def fall (n : Nat) : Nat → Nat
  | 0 => 1
  | k + 1 => fall n k * (n - k)

/-- Pairs among k picks: tri k = 0 + 1 + … + (k-1) = C(k,2). -/
def tri : Nat → Nat
  | 0 => 0
  | k + 1 => tri k + k

/-- (k+1)·C(n,k+1) = (n-k)·C(n,k). -/
theorem layer_step (n k : Nat) : (k + 1) * layer n (k + 1) = (n - k) * layer n k := by
  cases n with
  | zero =>
    rw [layer_eq_zero 0 (k + 1) (by omega)]; simp
  | succ m =>
    cases Nat.lt_or_ge m k with
    | inl h =>
      rw [layer_eq_zero (m + 1) (k + 1) (by omega)]
      have h0 : m + 1 - k = 0 := by omega
      rw [h0]; simp
    | inr h =>
      have ha := absorb m k
      have hp := layer_succ_succ m k
      cases k with
      | zero =>
        simp only [layer_zero_right, Nat.sub_zero] at *
        rw [hp] at ha
        generalize layer (m + 1) 1 = X at *
        generalize layer m 1 = B at *
        grind
      | succ j =>
        -- (j+2) C(m+1, j+2) = (m - j) C(m+1, j+1)
        have ha2 := absorb m (j + 1)
        have hp2 := layer_succ_succ m (j + 1)
        rw [hp2] at ha2
        have hpj := layer_succ_succ m j
        have hsub : m + 1 - (j + 1) = m - j := by omega
        rw [hsub, hpj]
        have hab := absorb m j
        rw [hpj] at hab
        generalize layer m (j + 1 + 1) = C at *
        generalize layer m (j + 1) = B at *
        generalize layer m j = A at *
        have hmj : m = (m - j) + j := by omega
        generalize m - j = d at *
        subst hmj
        grind

/-- **Picking with order.** C(n,k)·k! = n(n-1)…(n-k+1). -/
theorem layer_mul_fact (n : Nat) : ∀ k, layer n k * fact k = fall n k
  | 0 => by simp [layer_zero_right, fact, fall]
  | k + 1 => by
    have ih := layer_mul_fact n k
    have hs := layer_step n k
    show layer n (k + 1) * ((k + 1) * fact k) = fall n k * (n - k)
    rw [← ih, ← Nat.mul_assoc, Nat.mul_comm (layer n (k+1)) (k+1), hs]
    generalize layer n k = A; generalize fact k = F; generalize n - k = d
    grind

theorem fall_le_pow (n : Nat) : ∀ k, fall n k ≤ n ^ k
  | 0 => by simp [fall]
  | k + 1 => by
    show fall n k * (n - k) ≤ n ^ k * n
    exact Nat.mul_le_mul (fall_le_pow n k) (Nat.sub_le n k)

/-- **The layer bound.** C(n,k)·k! ≤ nᵏ: distinct picks are at most all picks. -/
theorem layer_fact_le (n k : Nat) : layer n k * fact k ≤ n ^ k := by
  rw [layer_mul_fact]; exact fall_le_pow n k

/-- **Collisions are few.** n·nᵏ ≤ n·fall n k + C(k,2)·nᵏ: the ordered picks that
    repeat an axis number at most C(k,2)·nᵏ⁻¹. -/
theorem collisions (n : Nat) : ∀ k, k ≤ n → n * n ^ k ≤ n * fall n k + tri k * n ^ k
  | 0, _ => by simp [fall, tri]
  | k + 1, h => by
    have ih := collisions n k (by omega)
    have hF := fall_le_pow n k
    show n * (n ^ k * n) ≤ n * (fall n k * (n - k)) + (tri k + k) * (n ^ k * n)
    have hsplit : n = (n - k) + k := by omega
    have h1 : n * (n * n ^ k) ≤ n * (n * fall n k + tri k * n ^ k) := Nat.mul_le_mul_left n ih
    have h2 : n * (k * fall n k) ≤ n * (k * n ^ k) := Nat.mul_le_mul_left n (Nat.mul_le_mul_left k hF)
    generalize n ^ k = A at *
    generalize fall n k = F at *
    generalize n - k = d at *
    subst hsplit
    grind


theorem pow2_le_fact : ∀ N, 2 ^ N ≤ fact (N + 1)
  | 0 => by decide
  | N + 1 => by
    have ih := pow2_le_fact N
    show 2 ^ N * 2 ≤ (N + 2) * fact (N + 1)
    have : 2 ^ N * 2 ≤ fact (N + 1) * (N + 2) := Nat.mul_le_mul ih (by omega)
    rw [Nat.mul_comm (N + 2)]; exact this

theorem tri_mono {j k : Nat} (h : j ≤ k) : tri j ≤ tri k := by
  induction k with
  | zero => have : j = 0 := by omega
            subst this; exact Nat.le_refl _
  | succ k ih =>
    rcases Nat.lt_or_ge j (k + 1) with hj | hj
    · exact Nat.le_trans (ih (by omega)) (by simp [tri])
    · have : j = k + 1 := by omega
      subst this; exact Nat.le_refl _

/-- **Bernoulli, with room to spare.** (c+1)^(k+2) + (k+1)(c+1)^k ≤ c^(k+2) + (k+2)(c+1)^(k+1). -/
theorem bernoulli (c : Nat) : ∀ k, (c + 1) ^ (k + 1) ≤ c ^ (k + 1) + (k + 1) * (c + 1) ^ k
  | 0 => by simp
  | k + 1 => by
    have ih := bernoulli c k
    have h := Nat.mul_le_mul_left c ih
    show (c + 1) ^ k * (c + 1) * (c + 1) ≤ c ^ (k + 1) * c + (k + 2) * ((c + 1) ^ k * (c + 1))
    rw [Nat.pow_succ] at ih h
    generalize (c + 1) ^ k = X at *
    generalize c ^ (k + 1) = Y at *
    grind

theorem bernoulli_strict (c k : Nat) :
    (c + 1) ^ (k + 2) + (k + 1) * (c + 1) ^ k ≤ c ^ (k + 2) + (k + 2) * (c + 1) ^ (k + 1) := by
  have ih := bernoulli c k
  have h := Nat.mul_le_mul_left c ih
  show (c + 1) ^ k * (c + 1) * (c + 1) + (k + 1) * (c + 1) ^ k ≤ c ^ (k + 1) * c + (k + 2) * ((c + 1) ^ k * (c + 1))
  rw [Nat.pow_succ] at ih h
  generalize (c + 1) ^ k = X at *
  generalize c ^ (k + 1) = Y at *
  grind

/-- **The staircase climbs** (cleared of fractions):
    (n+1)^(2n+1) < nⁿ·(n+2)^(n+1) for n ≥ 1, i.e. (1+1/n)ⁿ < (1+1/(n+1))ⁿ⁺¹. -/
theorem stair_nat_lt (m : Nat) :
    (m + 2) ^ (2 * (m + 1) + 1) < (m + 1) ^ (m + 1) * (m + 3) ^ (m + 2) := by
  have hb := bernoulli_strict ((m + 1) * (m + 3)) m
  have hc : (m + 1) * (m + 3) + 1 = (m + 2) ^ 2 := by grind
  rw [hc, ← Nat.pow_mul, ← Nat.pow_mul, ← Nat.pow_mul, Nat.mul_pow] at hb
  have e1 : (m + 2) ^ (2 * (m + 2)) = (m + 2) ^ (2 * m) * (m + 2) ^ 4 := by
    have : 2 * (m + 2) = 2 * m + 4 := by omega
    rw [this, Nat.pow_add]
  have e2 : (m + 2) ^ (2 * (m + 1)) = (m + 2) ^ (2 * m) * (m + 2) ^ 2 := by
    have : 2 * (m + 1) = 2 * m + 2 := by omega
    rw [this, Nat.pow_add]
  have e3 : (m + 2) ^ (2 * (m + 1) + 1) = (m + 2) ^ (2 * m) * (m + 2) ^ 3 := by
    have : 2 * (m + 1) + 1 = 2 * m + 3 := by omega
    rw [this, Nat.pow_add]
  have e4 : (m + 1) ^ (m + 2) = (m + 1) * (m + 1) ^ (m + 1) := by rw [Nat.pow_succ, Nat.mul_comm]
  rw [e1, e2, e4, Nat.mul_assoc (m + 1)] at hb
  rw [e3]
  have hQ : 0 < (m + 2) ^ (2 * m) := Nat.pow_pos (by omega)
  generalize (m + 2) ^ (2 * m) = Q at *
  generalize (m + 1) ^ (m + 1) * (m + 3) ^ (m + 2) = R at *
  have hb' : (m + 1) * (Q * (m + 2) ^ 3 + Q) ≤ (m + 1) * R := by
    grind
  have := Nat.le_of_mul_le_mul_left hb' (by omega)
  omega


/-! ## The circle operation on the rationals -/

/-- **The circle operation.** a ◯ b = a + b + ab = (1+a)(1+b) − 1. -/
def circ (a b : Rat) : Rat := a + b + a * b

/-- x combined with itself by the circle operation n times (0 is the circle-operation identity). -/
def mpow (x : Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => circ (mpow x n) x

/-- n copies of x under the circle operation make the box (1+x)ⁿ, minus the unit corner. -/
theorem mpow_eq (x : Rat) : ∀ n, mpow x n = (1 + x) ^ n - 1
  | 0 => by simp [mpow] <;> grind
  | n + 1 => by
    simp only [mpow, circ]; rw [mpow_eq x n, Rat.pow_succ]; grind

theorem pow_add' (a : Rat) (m : Nat) : ∀ n, a ^ (m + n) = a ^ m * a ^ n
  | 0 => by simp
  | n + 1 => by rw [← Nat.add_assoc, Rat.pow_succ, Rat.pow_succ, pow_add' a m n]; grind

/-- The exponential law: m steps then n steps is m + n steps. -/
theorem mpow_add (x : Rat) (m n : Nat) : mpow x (m + n) = circ (mpow x m) (mpow x n) := by
  rw [mpow_eq, mpow_eq, mpow_eq, pow_add']; simp only [circ]; grind

/-! ## Sums -/

def rsum (f : Nat → Rat) : Nat → Rat
  | 0 => 0
  | m + 1 => rsum f m + f m

theorem rsum_add (f g : Nat → Rat) : ∀ m, rsum (fun k => f k + g k) m = rsum f m + rsum g m
  | 0 => by simp [rsum] <;> grind
  | m + 1 => by simp only [rsum]; rw [rsum_add f g m]; grind

theorem rsum_mul (c : Rat) (f : Nat → Rat) : ∀ m, rsum (fun k => c * f k) m = c * rsum f m
  | 0 => by simp [rsum] <;> grind
  | m + 1 => by simp only [rsum]; rw [rsum_mul c f m]; grind

theorem rsum_shift (f : Nat → Rat) : ∀ m, rsum f (m + 1) = f 0 + rsum (fun k => f (k + 1)) m
  | 0 => by simp [rsum] <;> grind
  | m + 1 => by rw [rsum, rsum_shift f m]; simp only [rsum]; grind

theorem rsum_le (f g : Nat → Rat) : ∀ m, (∀ k, k < m → f k ≤ g k) → rsum f m ≤ rsum g m
  | 0, _ => Rat.le_refl
  | m + 1, h => by
    simp only [rsum]
    have := rsum_le f g m (fun k hk => h k (by omega))
    have := h m (by omega)
    grind

theorem rsum_nonneg (f : Nat → Rat) : ∀ m, (∀ k, k < m → 0 ≤ f k) → 0 ≤ rsum f m
  | 0, _ => Rat.le_refl
  | m + 1, h => by
    simp only [rsum]
    have := rsum_nonneg f m (fun k hk => h k (by omega))
    have := h m (by omega)
    grind

/-- More nonnegative terms, bigger sum. -/
theorem rsum_mono (f : Nat → Rat) (hf : ∀ k, 0 ≤ f k) {m m' : Nat} (h : m ≤ m') :
    rsum f m ≤ rsum f m' := by
  induction m' with
  | zero => have : m = 0 := by omega
            subst this; exact Rat.le_refl
  | succ m' ih =>
    rcases Nat.lt_or_ge m (m' + 1) with hm | hm
    · have := ih (by omega); have := hf m'; simp only [rsum]; grind
    · have : m = m' + 1 := by omega
      subst this; exact Rat.le_refl

/-! ## The cube's corners, weighted -/

/-- **The cube's corner count.** (1 + x)ⁿ = Σₖ C(n,k)·xᵏ. -/
theorem binomial (n : Nat) (x : Rat) :
    (1 + x) ^ n = rsum (fun k => (layer n k : Rat) * x ^ k) (n + 1) := by
  induction n with
  | zero => simp [rsum, layer] <;> grind
  | succ n ih =>
    have hT : rsum (fun k => (layer n (k + 1) : Rat) * x ^ (k + 1)) (n + 1)
        = rsum (fun k => (layer n k : Rat) * x ^ k) (n + 1) - 1 := by
      rw [rsum_shift (fun k => (layer n k : Rat) * x ^ k), rsum, layer_eq_zero n (n + 1) (by omega)]
      simp [layer_zero_right] <;> grind
    have step : rsum (fun k => (layer (n + 1) k : Rat) * x ^ k) (n + 1 + 1)
        = x * rsum (fun k => (layer n k : Rat) * x ^ k) (n + 1)
          + rsum (fun k => (layer n k : Rat) * x ^ k) (n + 1) := by
      rw [rsum_shift]
      have e : (fun k => ((layer (n + 1) (k + 1) : Nat) : Rat) * x ^ (k + 1)) =
          (fun k => x * ((layer n k : Rat) * x ^ k) + (layer n (k + 1) : Rat) * x ^ (k + 1)) := by
        funext k; rw [layer_succ_succ, Rat.natCast_add, Rat.pow_succ]; grind
      rw [e, rsum_add, rsum_mul, hT]; simp [layer_zero_right]; grind
    rw [Rat.pow_succ, ih, step]; grind


/-! ## Helpers on the rationals -/

theorem le_div_of_mul_le {a b c : Rat} (hc : 0 < c) (h : a * c ≤ b) : a ≤ b / c :=
  Rat.not_lt.mp (fun h' => absurd ((Rat.div_lt_iff hc).mp h') (Rat.not_lt.mpr h))

theorem one_div_nonneg {a : Rat} (ha : 0 ≤ a) : 0 ≤ 1 / a := by
  rcases Rat.le_iff_lt_or_eq.mp ha with h | h
  · exact le_div_of_mul_le h (by grind)
  · rw [← h]; decide +kernel

theorem one_div_le_one_div {a b : Rat} (hb : 0 < b) (h : b ≤ a) : 1 / a ≤ 1 / b := by
  have ha : 0 < a := by grind
  apply le_div_of_mul_le hb
  have h1 : 1 / a * b ≤ 1 / a * a := Rat.mul_le_mul_of_nonneg_left h (one_div_nonneg (Rat.le_of_lt ha))
  have h2 : 1 / a * a = 1 := Rat.div_mul_cancel (Rat.ne_of_gt ha)
  grind

theorem natcast_ne {a : Nat} (h : a ≠ 0) : (a : Rat) ≠ 0 :=
  fun h' => h (Rat.natCast_eq_zero_iff.mp h')

theorem fact_cast_pos (k : Nat) : (0 : Rat) < (fact k : Rat) := Rat.natCast_pos.mpr (fact_pos k)

theorem inv_pow_mul (N : Rat) (hN : N ≠ 0) : ∀ k : Nat, (1 / N) ^ k * N ^ k = 1
  | 0 => by simp
  | k + 1 => by rw [Rat.pow_succ, Rat.pow_succ]; have := inv_pow_mul N hN k; grind

theorem mul_pow' (a b : Rat) : ∀ n : Nat, (a * b) ^ n = a ^ n * b ^ n
  | 0 => by simp
  | n + 1 => by rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_succ, mul_pow' a b n]; grind

/-! ## e from the circle operation -/

/-- **The staircase.** 1/n combined with itself by the circle operation n times, plus the unit corner. -/
def stair (n : Nat) : Rat := (1 + 1 / (n : Rat)) ^ n

theorem stair_eq_mpow (n : Nat) : stair n = mpow (1 / (n : Rat)) n + 1 := by
  rw [mpow_eq]; unfold stair; grind

/-- The weight of layer k of the n-cube when every switch weighs 1/n. -/
def wt (n k : Nat) : Rat := (layer n k : Rat) * (1 / (n : Rat)) ^ k

theorem stair_eq_cube (n : Nat) : stair n = rsum (wt n) (n + 1) := binomial n _

theorem wt_nonneg (n k : Nat) : 0 ≤ wt n k := by
  unfold wt
  have h1 : (0 : Rat) ≤ (layer n k : Rat) := Rat.natCast_nonneg
  have h2 : (0 : Rat) ≤ (1 / (n : Rat)) ^ k := Rat.pow_nonneg (one_div_nonneg Rat.natCast_nonneg)
  have := Rat.mul_le_mul_of_nonneg_left h2 h1
  grind

/-- **Each layer is at most 1/k!.** C(n,k)/nᵏ ≤ 1/k!. -/
theorem wt_le (n k : Nat) (hn : 1 ≤ n) : wt n k ≤ 1 / (fact k : Rat) := by
  have hN : (n : Rat) ≠ 0 := natcast_ne (by omega)
  have hu := inv_pow_mul (n : Rat) hN k
  have hb : (layer n k : Rat) * (fact k : Rat) ≤ (n : Rat) ^ k := by
    have := Rat.natCast_le_natCast.mpr (layer_fact_le n k)
    rw [Rat.natCast_mul, Rat.natCast_pow] at this; exact this
  have hu0 : 0 ≤ (1 / (n : Rat)) ^ k := Rat.pow_nonneg (one_div_nonneg Rat.natCast_nonneg)
  have h := Rat.mul_le_mul_of_nonneg_right hb hu0
  apply le_div_of_mul_le (fact_cast_pos k)
  unfold wt
  generalize (1 / (n : Rat)) ^ k = u at *
  generalize (n : Rat) ^ k = P at *
  generalize (layer n k : Rat) = C at *
  generalize (fact k : Rat) = K at *
  grind

/-- **Each layer is at least 1/k! minus its collisions.** 1/k! ≤ C(n,k)/nᵏ + C(k,2)/(n·k!). -/
theorem wt_ge (n k : Nat) (hn : 1 ≤ n) (hk : k ≤ n) :
    1 / (fact k : Rat) ≤ wt n k + (tri k : Rat) / ((n : Rat) * (fact k : Rat)) := by
  have hN0 : (0 : Rat) < (n : Rat) := Rat.natCast_pos.mpr (by omega)
  have hN : (n : Rat) ≠ 0 := Rat.ne_of_gt hN0
  have hK := fact_cast_pos k
  have hu := inv_pow_mul (n : Rat) hN k
  have hc : (n : Rat) * (n : Rat) ^ k ≤ (n : Rat) * ((layer n k : Rat) * (fact k : Rat)) + (tri k : Rat) * (n : Rat) ^ k := by
    have := Rat.natCast_le_natCast.mpr (collisions n k hk)
    rw [← layer_mul_fact] at this
    simp only [Rat.natCast_mul, Rat.natCast_add, Rat.natCast_pow] at this; exact this
  have hu0 : 0 ≤ (1 / (n : Rat)) ^ k := Rat.pow_nonneg (one_div_nonneg Rat.natCast_nonneg)
  have h := Rat.mul_le_mul_of_nonneg_right hc hu0
  have hK1 : 1 / (fact k : Rat) * (fact k : Rat) = 1 := Rat.div_mul_cancel (Rat.ne_of_gt hK)
  have hNK : 0 < (n : Rat) * (fact k : Rat) := Rat.mul_pos hN0 hK
  have key : (1 / (fact k : Rat) - wt n k) * ((n : Rat) * (fact k : Rat)) ≤ (tri k : Rat) := by
    unfold wt
    generalize (1 / (n : Rat)) ^ k = u at *
    generalize (n : Rat) ^ k = P at *
    generalize 1 / (fact k : Rat) = r at *
    generalize (layer n k : Rat) = C at *
    generalize (fact k : Rat) = K at *
    generalize (tri k : Rat) = T at *
    generalize (n : Rat) = N at *
    have e1 : N * P * u = N := by
      have : N * P * u = N * (u * P) := by grind
      rw [this, hu]; grind
    have e2 : (N * (C * K) + T * P) * u = N * C * K * u + T := by
      have : (N * (C * K) + T * P) * u = N * C * K * u + T * (u * P) := by grind
      rw [this, hu]; grind
    have e3 : (r - C * u) * (N * K) = N * (r * K) - N * C * K * u := by grind
    rw [e3, hK1]; rw [e1, e2] at h; grind
  have := le_div_of_mul_le hNK key
  grind

/-- The partial sums of e: s N = 1/0! + 1/1! + … + 1/N!. -/
def s (N : Nat) : Rat := rsum (fun k => 1 / (fact k : Rat)) (N + 1)

/-- The upper partner: u N = s N + 1/(N!·N). -/
def u (N : Nat) : Rat := s N + 1 / ((fact N : Rat) * (N : Rat))

theorem stair_le_s (n : Nat) (hn : 1 ≤ n) : stair n ≤ s n := by
  rw [stair_eq_cube]
  exact rsum_le _ _ _ (fun k _ => wt_le n k hn)

theorem s_bound : ∀ N, s N ≤ 3 - 2 / (2 : Rat) ^ N
  | 0 => by decide +kernel
  | N + 1 => by
    have ih := s_bound N
    have hP : (0 : Rat) < (2 : Rat) ^ N := Rat.pow_pos (by decide)
    have hf : 1 / (fact (N + 1) : Rat) ≤ 1 / (2 : Rat) ^ N := by
      apply one_div_le_one_div hP
      have := Rat.natCast_le_natCast.mpr (pow2_le_fact N)
      rw [Rat.natCast_pow] at this; exact this
    show rsum (fun k => 1 / (fact k : Rat)) (N + 1) + 1 / (fact (N + 1) : Rat) ≤ 3 - 2 / (2 : Rat) ^ (N + 1)
    have hs : s N = rsum (fun k => 1 / (fact k : Rat)) (N + 1) := rfl
    rw [← hs, Rat.pow_succ]
    generalize (2 : Rat) ^ N = P at *
    generalize 1 / (fact (N + 1) : Rat) = r at *
    have e : 2 / (P * 2) = 1 / P := by grind
    rw [e]; grind

/-- **The partial sums never reach 3.** -/
theorem s_lt_three (N : Nat) : s N < 3 := by
  have h := s_bound N
  have hP : (0 : Rat) < (2 : Rat) ^ N := Rat.pow_pos (by decide)
  have : (0 : Rat) < 2 / (2 : Rat) ^ N := (Rat.lt_div_iff hP).mpr (by grind)
  grind

theorem s_mono {N M : Nat} (h : N ≤ M) : s N ≤ s M :=
  rsum_mono _ (fun k => one_div_nonneg (Rat.le_of_lt (fact_cast_pos k))) (by omega)

theorem s_le_u (N : Nat) : s N ≤ u N := by
  unfold u
  have : 0 ≤ 1 / ((fact N : Rat) * (N : Rat)) :=
    one_div_nonneg (by have := Rat.mul_le_mul_of_nonneg_left (Rat.natCast_nonneg (a := N)) (Rat.le_of_lt (fact_cast_pos N)); grind)
  grind

/-- **The upper partner only comes down.** u (N+1) ≤ u N for N ≥ 1. -/
theorem u_step (N : Nat) (hN : 1 ≤ N) : u (N + 1) ≤ u N := by
  have hF := fact_cast_pos N
  have hN1 : (1 : Rat) ≤ (N : Rat) := by
    have := Rat.natCast_le_natCast.mpr hN; simpa using this
  show rsum (fun k => 1 / (fact k : Rat)) (N + 1) + 1 / (fact (N + 1) : Rat)
      + 1 / ((fact (N + 1) : Rat) * ((N + 1 : Nat) : Rat))
      ≤ rsum (fun k => 1 / (fact k : Rat)) (N + 1) + 1 / ((fact N : Rat) * (N : Rat))
  have e : (fact (N + 1) : Rat) = ((N : Rat) + 1) * (fact N : Rat) := by
    show (((N + 1) * fact N : Nat) : Rat) = _
    rw [Rat.natCast_mul, Rat.natCast_add]; simp
  have e2 : ((N + 1 : Nat) : Rat) = (N : Rat) + 1 := by rw [Rat.natCast_add]; simp
  rw [e, e2]
  generalize (fact N : Rat) = F at *
  generalize (N : Rat) = x at *
  have hgoal : (1 / ((x + 1) * F) + 1 / ((x + 1) * F * (x + 1))) * (F * x) ≤ 1 := by
    have : (1 / ((x + 1) * F) + 1 / ((x + 1) * F * (x + 1))) * (F * x) = x * (x + 2) / ((x + 1) * (x + 1)) := by grind
    rw [this]
    have hpos : (0 : Rat) < (x + 1) * (x + 1) := Rat.mul_pos (by grind) (by grind)
    exact Rat.not_lt.mp (fun h => by have := (Rat.lt_div_iff hpos).mp h; grind)
  have hFx : 0 < F * x := Rat.mul_pos hF (by grind)
  have := le_div_of_mul_le hFx hgoal
  grind

theorem u_anti {N M : Nat} (hN : 1 ≤ N) (h : N ≤ M) : u M ≤ u N := by
  induction M with
  | zero => omega
  | succ M ih =>
    rcases Nat.lt_or_ge N (M + 1) with hm | hm
    · exact Rat.le_trans (u_step M (by omega)) (ih (by omega))
    · have : N = M + 1 := by omega
      subst this; exact Rat.le_refl

/-- **The ceiling.** Every step of the staircase sits under every u N. -/
theorem stair_le_u (n N : Nat) (hn : 1 ≤ n) (hN : 1 ≤ N) : stair n ≤ u N := by
  rcases Nat.le_total n N with h | h
  · exact Rat.le_trans (stair_le_s n hn) (Rat.le_trans (s_mono h) (s_le_u N))
  · exact Rat.le_trans (stair_le_s n hn) (Rat.le_trans (s_le_u n) (u_anti hN h))

/-- **The floor.** For n ≥ N, the staircase is within 3·C(N,2)/n of s N. -/
theorem s_le_stair (n N : Nat) (hn : 1 ≤ n) (hNn : N ≤ n) :
    s N - 3 * (tri N : Rat) / (n : Rat) ≤ stair n := by
  have hN0 : (0 : Rat) < (n : Rat) := Rat.natCast_pos.mpr (by omega)
  obtain ⟨c, hc⟩ : ∃ c : Rat, c = (tri N : Rat) / (n : Rat) := ⟨_, rfl⟩
  have hc0 : 0 ≤ c := by
    rw [hc]; exact le_div_of_mul_le hN0 (by have : (0:Rat) ≤ (tri N : Rat) := Rat.natCast_nonneg; grind)
  -- drop the layers above N
  have h1 : rsum (wt n) (N + 1) ≤ stair n := by
    rw [stair_eq_cube]; exact rsum_mono _ (wt_nonneg n) (by omega)
  -- each kept layer is at least 1/k! - c/k!
  have h2 : rsum (fun k => 1 / (fact k : Rat) + (-c) * (1 / (fact k : Rat))) (N + 1) ≤ rsum (wt n) (N + 1) := by
    apply rsum_le
    intro k hk
    have hg := wt_ge n k hn (by omega)
    have hK := fact_cast_pos k
    have ht : (tri k : Rat) / ((n : Rat) * (fact k : Rat)) ≤ c * (1 / (fact k : Rat)) := by
      have htk : (tri k : Rat) ≤ (tri N : Rat) := Rat.natCast_le_natCast.mpr (tri_mono (by omega))
      apply Rat.not_lt.mp
      intro hlt
      have hNK : 0 < (n : Rat) * (fact k : Rat) := Rat.mul_pos hN0 hK
      have := (Rat.lt_div_iff hNK).mp hlt
      have e : c * (1 / (fact k : Rat)) * ((n : Rat) * (fact k : Rat)) = (tri N : Rat) := by
        rw [hc]; grind
      grind
    grind
  rw [rsum_add, rsum_mul] at h2
  have hs : rsum (fun k => 1 / (fact k : Rat)) (N + 1) = s N := rfl
  rw [hs] at h2
  have h3 := s_lt_three N
  have hs0 : 0 ≤ s N := rsum_nonneg _ _ (fun k _ => one_div_nonneg (Rat.le_of_lt (fact_cast_pos k)))
  have h4 : c * s N ≤ c * 3 := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt h3) hc0
  have e : 3 * (tri N : Rat) / (n : Rat) = c * 3 := by rw [hc]; grind
  rw [e]
  grind

/-- **The staircase climbs.** (1 + 1/n)ⁿ < (1 + 1/(n+1))ⁿ⁺¹ for n ≥ 1. -/
theorem stair_lt_succ (n : Nat) (hn : 1 ≤ n) : stair n < stair (n + 1) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have hA : stair (m + 1) * ((m + 1 : Nat) : Rat) ^ (m + 1) = ((m + 2 : Nat) : Rat) ^ (m + 1) := by
    unfold stair; rw [← mul_pow']; congr 1
    have : ((m + 1 : Nat) : Rat) ≠ 0 := natcast_ne (by omega)
    rw [Rat.natCast_add, Rat.natCast_add] at *; grind
  have hB : stair (m + 2) * ((m + 2 : Nat) : Rat) ^ (m + 2) = ((m + 3 : Nat) : Rat) ^ (m + 2) := by
    unfold stair; rw [← mul_pow']; congr 1
    have : ((m + 2 : Nat) : Rat) ≠ 0 := natcast_ne (by omega)
    rw [Rat.natCast_add, Rat.natCast_add] at *; grind
  have hnat := Rat.natCast_lt_natCast.mpr (stair_nat_lt m)
  simp only [Rat.natCast_mul, Rat.natCast_pow] at hnat
  have e : ((m + 2 : Nat) : Rat) ^ (2 * (m + 1) + 1) = ((m + 2 : Nat) : Rat) ^ (m + 1) * ((m + 2 : Nat) : Rat) ^ (m + 2) := by
    rw [← pow_add']; congr 1; omega
  rw [e, ← hA, ← hB] at hnat
  have hP : (0 : Rat) < ((m + 1 : Nat) : Rat) ^ (m + 1) := Rat.pow_pos (Rat.natCast_pos.mpr (by omega))
  have hS : (0 : Rat) < ((m + 2 : Nat) : Rat) ^ (m + 2) := Rat.pow_pos (Rat.natCast_pos.mpr (by omega))
  have hPS := Rat.mul_pos hP hS
  apply Rat.lt_of_mul_lt_mul_right (c := ((m + 1 : Nat) : Rat) ^ (m + 1) * ((m + 2 : Nat) : Rat) ^ (m + 2)) _ (Rat.le_of_lt hPS)
  generalize ((m + 1 : Nat) : Rat) ^ (m + 1) = P at *
  generalize ((m + 2 : Nat) : Rat) ^ (m + 2) = S at *
  generalize ((m + 3 : Nat) : Rat) ^ (m + 2) = T at *
  grind

/-! ## Small checks, computed exactly by the kernel -/

/-- Layer k of the 100-cube against 1/k!: k = 3 gives 0.9702/6 versus 1/6. -/
example : wt 100 3 = 4851 / 30000 := by decide +kernel
/-- e is trapped: s 10 = 2.71828180… and u 10 = 2.71828183…. -/
example : (27182818 / 10000000 : Rat) < s 10 ∧ u 10 < 27182819 / 10000000 := by decide +kernel
/-- The staircase at n = 100 is about 2.7048. -/
example : (2704 / 1000 : Rat) < stair 100 ∧ stair 100 < 2705 / 1000 := by decide +kernel

end CircleOp.E
