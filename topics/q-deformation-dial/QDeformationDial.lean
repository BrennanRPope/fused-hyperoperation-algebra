/-
  Is the plane dial's d Hoffman–Ihara's q-deformation? (DRAFT, version 1)
  =======================================================================

  Part of Fused Hyperoperation Algebra. Page and checks: ../ (q-deformation-dial).
  Written for the Mathlib library (Lean 4.37.0-rc1). Built only by CI's report-only Mathlib step; the
  machine that wrote it has no Lean.

  Proved here:
  * `qnum_add`: q-integers add by our line law, [m+n]_q = [m] + [n] + (q − 1)[m][n]. So q = 1 + k.
  * `q_letter_product`: the q-MZV letter values fᵢ(n) = q^((i−1)n) / [n]^i obey Hoffman–Ihara's
    zᵢ ◇ zⱼ = z_(i+j) + (1 − q) z_(i+j−1); the correction is q^n = 1 + (q − 1)[n], the line dial again.
  * `jet_ratio`: jets a + bN with N² = kN − d multiply, and the ratio b/a combines by the plane law.
  * `orbit_add`: the orbit xₙ of 1 + x₀N satisfies x_(m+n) = F(x_m, x_n).
  * `orbit_d_zero`: with d = 0 the orbit is x₀ [n]_q with q = 1 + k x₀ (one base).
  * `orbit_two_bases`: with k = α + β, d = αβ, (α − β) Bₙ = Pⁿ − Qⁿ and (α − β) Aₙ = αQⁿ − βPⁿ,
    P = 1 + x₀α, Q = 1 + x₀β: two bases, the (p,q)-numbers. d = 0 exactly when one base is 1.

  Known: q-integers and their addition; Hoffman–Ihara (2017) §6.4; (p,q)-numbers (Chakrabarti and
  Jagannathan 1991). Ours: reading q as the line dial's k and d as the second base.
-/
import Mathlib

namespace QDeformationDial

/-- The q-integer [n]_q = 1 + q + … + q^(n−1). -/
noncomputable def qnum (q : ℝ) (n : ℕ) : ℝ := ∑ j ∈ Finset.range n, q ^ j

theorem qnum_succ (q : ℝ) (n : ℕ) : qnum q (n + 1) = qnum q n + q ^ n := by
  unfold qnum; rw [Finset.sum_range_succ]

/-- q^n = 1 + (q − 1)[n]_q. -/
theorem pow_eq (q : ℝ) (n : ℕ) : q ^ n = 1 + (q - 1) * qnum q n := by
  induction n with
  | zero => simp [qnum]
  | succ n ih => rw [qnum_succ, pow_succ, ih]; ring

theorem qnum_add' (q : ℝ) (m n : ℕ) : qnum q (m + n) = qnum q m + q ^ m * qnum q n := by
  induction n with
  | zero => simp [qnum]
  | succ n ih => rw [← add_assoc, qnum_succ, ih, qnum_succ, pow_add]; ring

/-- q-integers add by the line law with k = q − 1. -/
theorem qnum_add (q : ℝ) (m n : ℕ) :
    qnum q (m + n) = qnum q m + qnum q n + (q - 1) * qnum q m * qnum q n := by
  rw [qnum_add', pow_eq q m]; ring

/-- The q-MZV letter value fᵢ(n) = q^((i−1)n) / [n]^i, written with i = a + 1. -/
noncomputable def f (q : ℝ) (n a : ℕ) : ℝ := q ^ (a * n) / qnum q n ^ (a + 1)

/-- Hoffman–Ihara's q letter product on the letter values. -/
theorem q_letter_product (q : ℝ) (n a b : ℕ) (h : qnum q n ≠ 0) :
    f q n a * f q n b = f q n (a + b + 1) + (1 - q) * f q n (a + b) := by
  unfold f
  have e1 : q ^ ((a + b + 1) * n) = q ^ (a * n) * q ^ (b * n) * q ^ n := by
    rw [← pow_add, ← pow_add]; congr 1; ring
  have e2 : q ^ ((a + b) * n) = q ^ (a * n) * q ^ (b * n) := by
    rw [← pow_add]; congr 1; ring
  rw [e1, e2, pow_eq q n]
  first | (field_simp; ring) | field_simp

/-! ### Jets a + bN with N² = kN − d -/

/-- Product of jets (a + bN)(c + eN) with N² = kN − d. -/
def jmul (k d : ℝ) (u v : ℝ × ℝ) : ℝ × ℝ :=
  (u.1 * v.1 - d * u.2 * v.2, u.1 * v.2 + u.2 * v.1 + k * u.2 * v.2)

/-- The plane law. -/
noncomputable def F (k d x y : ℝ) : ℝ := (x + y + k * x * y) / (1 - d * x * y)

theorem jet_ratio (k d a b c e : ℝ) (ha : a ≠ 0) (hc : c ≠ 0)
    (h : a * c - d * b * e ≠ 0) :
    (jmul k d (a, b) (c, e)).2 / (jmul k d (a, b) (c, e)).1 = F k d (b / a) (e / c) := by
  have h' : 1 - d * (b / a) * (e / c) ≠ 0 := by
    intro h0; apply h
    have : a * c * (1 - d * (b / a) * (e / c)) = a * c - d * b * e := by
      first | (field_simp; ring) | field_simp
    rw [← this, h0, mul_zero]
  show (a * e + b * c + k * b * e) / (a * c - d * b * e)
      = (b / a + e / c + k * (b / a) * (e / c)) / (1 - d * (b / a) * (e / c))
  rw [div_eq_div_iff h h']
  first | (field_simp; ring) | field_simp

theorem jmul_comm (k d : ℝ) (u v : ℝ × ℝ) : jmul k d u v = jmul k d v u := by
  ext <;> dsimp only [jmul] <;> ring

theorem jmul_assoc (k d : ℝ) (u v w : ℝ × ℝ) :
    jmul k d (jmul k d u v) w = jmul k d u (jmul k d v w) := by
  ext <;> dsimp only [jmul] <;> ring

/-- Powers of the jet 1 + x₀N. -/
def jpow (k d x0 : ℝ) : ℕ → ℝ × ℝ
  | 0 => (1, 0)
  | n + 1 => jmul k d (jpow k d x0 n) (1, x0)

theorem jpow_add (k d x0 : ℝ) (m n : ℕ) :
    jpow k d x0 (m + n) = jmul k d (jpow k d x0 m) (jpow k d x0 n) := by
  induction n with
  | zero => simp [jpow, jmul]
  | succ n ih =>
    rw [← add_assoc]
    show jmul k d (jpow k d x0 (m + n)) (1, x0)
      = jmul k d (jpow k d x0 m) (jmul k d (jpow k d x0 n) (1, x0))
    rw [ih, jmul_assoc]

/-- The orbit xₙ = Bₙ / Aₙ adds by the plane law. -/
theorem orbit_add (k d x0 : ℝ) (m n : ℕ)
    (hm : (jpow k d x0 m).1 ≠ 0) (hn : (jpow k d x0 n).1 ≠ 0)
    (h : (jpow k d x0 m).1 * (jpow k d x0 n).1 - d * (jpow k d x0 m).2 * (jpow k d x0 n).2 ≠ 0) :
    (jpow k d x0 (m + n)).2 / (jpow k d x0 (m + n)).1
      = F k d ((jpow k d x0 m).2 / (jpow k d x0 m).1) ((jpow k d x0 n).2 / (jpow k d x0 n).1) := by
  rw [jpow_add]
  exact jet_ratio k d _ _ _ _ hm hn h

/-- With d = 0 the orbit is one base: Aₙ = 1 and Bₙ = x₀ [n]_q with q = 1 + k x₀. -/
theorem orbit_d_zero (k x0 : ℝ) (n : ℕ) :
    jpow k 0 x0 n = (1, x0 * qnum (1 + k * x0) n) := by
  induction n with
  | zero => simp [jpow, qnum]
  | succ n ih =>
    show jmul k 0 (jpow k 0 x0 n) (1, x0) = _
    rw [ih, qnum_succ, pow_eq (1 + k * x0) n]
    ext <;> dsimp only [jmul] <;> ring

/-- Two bases: with k = α + β and d = αβ, the orbit is built from Pⁿ and Qⁿ. -/
theorem orbit_two_bases (α β x0 : ℝ) (n : ℕ) :
    (α - β) * (jpow (α + β) (α * β) x0 n).2 = (1 + x0 * α) ^ n - (1 + x0 * β) ^ n ∧
    (α - β) * (jpow (α + β) (α * β) x0 n).1 = α * (1 + x0 * β) ^ n - β * (1 + x0 * α) ^ n := by
  induction n with
  | zero => simp [jpow]
  | succ n ih =>
    obtain ⟨hB, hA⟩ := ih
    rw [show jpow (α + β) (α * β) x0 (n + 1)
        = jmul (α + β) (α * β) (jpow (α + β) (α * β) x0 n) (1, x0) from rfl]
    dsimp only [jmul]
    constructor
    · rw [pow_succ, pow_succ]
      linear_combination (1 + (α + β) * x0) * hB + x0 * hA
    · rw [pow_succ, pow_succ]
      linear_combination hA - α * β * x0 * hB

/-- d = αβ is zero exactly when one base 1 + x₀α or 1 + x₀β equals 1 (for x₀ ≠ 0). -/
theorem d_zero_iff_base_one (α β x0 : ℝ) (hx : x0 ≠ 0) :
    α * β = 0 ↔ (1 + x0 * α = 1 ∨ 1 + x0 * β = 1) := by
  constructor
  · intro h
    rcases mul_eq_zero.1 h with h | h
    · left; rw [h]; ring
    · right; rw [h]; ring
  · rintro (h | h)
    · have : x0 * α = 0 := by linarith
      rcases mul_eq_zero.1 this with h1 | h1
      · exact absurd h1 hx
      · rw [h1]; ring
    · have : x0 * β = 0 := by linarith
      rcases mul_eq_zero.1 this with h1 | h1
      · exact absurd h1 hx
      · rw [h1]; ring

end QDeformationDial
