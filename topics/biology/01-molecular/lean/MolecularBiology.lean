/-
  Fused Hyperoperation Algebra: molecular biology (DNA, RNA, proteins).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base and the
  Lee–Yang module (for its 2×2 matrices and the wall theorem for rings).
  Intended library path: FusedAlgebra/Topics/Biology/MolecularBiology.lean

  1. Helix–coil (Zimm–Bragg) is the dial gas: with the two states relabelled, its transfer
     matrix is the Lee–Yang ring matrix with dial k = 1/σ and fugacity z = σ·s. So every odd
     ring of a helix–coil chain vanishes at s = −1, the wall, for every σ.
  2. Two binding sites: Z = 1 + a + b + ω·ab = 1 + (a ◯_ω b). The cooperativity ω is the dial.
  3. MWC allostery is a sum of two mult-add powers, dial 1 (R state) and dial c (T state).
  4. Copying fidelity: catches combine by the "or" law (dial −1); Eigen's error threshold is
     where the fitness advantage is the mult-add inverse of the copy loss.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.LeeYang

namespace MolecularBiology
open CircleOp LeeYang LeeYang.Mat2

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

theorem star_def (a b : R) : a ◯ b = a + b + a * b := rfl

/-! ## 1. Helix–coil is the dial gas -/

/-- Zimm–Bragg transfer matrix, states ordered (helix, coil): a helix after a helix weighs s,
    a helix after a coil weighs σ·s (nucleation), a coil weighs 1. -/
def zb (σ s : R) : Mat2 R := ⟨s, 1, σ * s, 1⟩

/-- Relabel the two states (conjugate by the swap). -/
def swap (x : Mat2 R) : Mat2 R := ⟨x.d, x.c, x.b, x.a⟩

theorem swap_mul (x y : Mat2 R) : swap (mul x y) = mul (swap x) (swap y) := by
  apply Mat2.ext' <;> simp only [swap, mul] <;> grind

theorem swap_one : swap (one : Mat2 R) = one := by
  apply Mat2.ext' <;> simp only [swap, one]

theorem swap_pow (x : Mat2 R) (n : Nat) : swap (pow x n) = pow (swap x) n := by
  induction n with
  | zero => simp only [pow]; exact swap_one
  | succ n ih => simp only [pow]; rw [swap_mul, ih]

theorem tr_swap (x : Mat2 R) : tr (swap x) = tr x := by
  simp only [tr, swap]; grind

/-- Relabelling the states does not change any ring's partition function. -/
theorem zb_ring_relabel (σ s : R) (n : Nat) : tr (pow (swap (zb σ s)) n) = tr (pow (zb σ s) n) := by
  rw [← swap_pow, tr_swap]

/-- **Helix–coil = dial gas.** If k·σ = 1, the relabelled Zimm–Bragg matrix is the Lee–Yang
    ring matrix with dial k and fugacity σ·s. -/
theorem zb_is_dial_gas (k σ s : R) (h : k * σ = 1) : swap (zb σ s) = ringT k (σ * s) := by
  apply Mat2.ext' <;> simp only [swap, zb, ringT]
  have : k * (σ * s) = (k * σ) * s := by grind
  rw [this, h]; grind

/-- So the helix–coil ring has the same partition function as the dial-gas ring. -/
theorem zb_ring_eq_dial_ring (k σ s : R) (h : k * σ = 1) (n : Nat) :
    tr (pow (zb σ s) n) = tr (pow (ringT k (σ * s)) n) := by
  rw [← zb_ring_relabel, zb_is_dial_gas k σ s h]

theorem zb_tr (σ s : R) : tr (zb σ s) = s + 1 := by simp only [zb, tr]

theorem zb_det (σ s : R) : det (zb σ s) = s * (1 - σ) := by simp only [zb, det]; grind

/-- **Wall theorem for helix–coil.** At s = −1 every odd ring vanishes, for every σ. -/
theorem zb_odd_ring_wall (σ : R) (m : Nat) : tr (pow (zb σ (-1)) (2 * m + 1)) = 0 :=
  tr_odd_pow_zero _ (by rw [zb_tr]; grind) m

/-- σ = 1 (no nucleation penalty): the matrix has determinant 0, so the chain is a free
    gas whose ring partition function is (1 + s)ⁿ. -/
theorem zb_free (s : R) : det (zb 1 s) = 0 := by rw [zb_det]; grind

/-! ## 2. Two binding sites: cooperativity is the dial -/

/-- The two-site binding polynomial. -/
def bind2 (ω a b : R) : R := 1 + a + b + ω * a * b

theorem bind2_is_circK (ω a b : R) : bind2 ω a b = 1 + circK ω a b := by
  unfold bind2 circK; grind

/-- Independent sites (ω = 1): the polynomial factors, plain mult-add. -/
theorem bind2_independent (a b : R) : bind2 1 a b = (1 + a) * (1 + b) := by
  unfold bind2; grind

theorem bind2_independent' (a b : R) : bind2 1 a b = 1 + (a ◯ b) := by
  unfold bind2 circ; grind

/-- Exclusive sites (ω = 0): only one can be filled, plain addition. -/
theorem bind2_exclusive (a b : R) : bind2 0 a b = 1 + a + b := by
  unfold bind2; grind

/-- The deviation from independence is exactly (ω − 1) times the corner cell. -/
theorem bind2_excess (ω a b : R) : bind2 ω a b - (1 + a) * (1 + b) = (ω - 1) * (a * b) := by
  unfold bind2; grind

/-- Half saturation for identical sites: the bound fraction (x + ωx²)/(1 + 2x + ωx²) is 1/2
    exactly when ω·x² = 1. Stated without division. -/
theorem bind2_half (ω x : R) : 2 * (x + ω * x * x) - bind2 ω x x = ω * x * x - 1 := by
  unfold bind2; grind

/-! ## 3. MWC allostery: two mult-add powers -/

/-- R state: the plain mult-add power, lifted. -/
theorem mwc_R (a : R) (n : Nat) : (1 + a) ^ n = 1 + mpow a n := (mpow_lift a n).symm

/-- T state: the weighted power at dial c, lifted. -/
theorem mwc_T (c a : R) (n : Nat) : (1 + c * a) ^ n = 1 + c * wpow c a n := (wpow_lift c a n).symm

/-- The MWC binding polynomial (1 + a)ⁿ + L(1 + c·a)ⁿ written as two mult-add powers. -/
theorem mwc_poly (L c a : R) (n : Nat) :
    (1 + a) ^ n + L * (1 + c * a) ^ n = (1 + mpow a n) + L * (1 + c * wpow c a n) := by
  rw [mwc_R, mwc_T]

/-- With c = 1 both states bind alike and MWC collapses to independent sites. -/
theorem mwc_c_one (L a : R) (n : Nat) :
    (1 + a) ^ n + L * (1 + 1 * a) ^ n = (1 + L) * (1 + a) ^ n := by
  have : (1 : R) * a = a := by grind
  rw [this]; grind

/-! ## 4. Copying fidelity -/

/-- A mistake survives two checks only if both miss: catches combine by the "or" law. -/
theorem catch_two (c₁ c₂ : R) : 1 - circK (-1) c₁ c₂ = (1 - c₁) * (1 - c₂) := by
  unfold circK; grind

theorem catch_three (c₁ c₂ c₃ : R) :
    1 - circK (-1) (circK (-1) c₁ c₂) c₃ = (1 - c₁) * (1 - c₂) * (1 - c₃) := by
  unfold circK; grind

/-- Copying a genome of L bases with error μ per base keeps the master with chance
    (1 − μ)^L = 1 + mpow(−μ, L): the copy loss is a mult-add power. -/
theorem copy_loss (μ : R) (L : Nat) : (1 - μ) ^ L = 1 + mpow (-μ) L := by
  rw [mpow_lift]
  have : (1 : R) - μ = 1 + -μ := by grind
  rw [this]

/-- Eigen: the master's growth factor σ·(1 + l) minus 1 is the mult-add of the advantage
    s = σ − 1 and the copy loss l. -/
theorem eigen_excess (s l : R) : (1 + s) * (1 + l) - 1 = s ◯ l := by
  unfold circ; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- **Eigen's error threshold is the mult-add inverse.** Away from the wall (1 + l ≠ 0), the
    master just holds its own (s ◯ l = 0) exactly when s is the mult-add inverse of l. -/
theorem eigen_threshold (s l : K) (h : 1 + l ≠ 0) : s ◯ l = 0 ↔ s = circInv l := by
  unfold circ circInv
  constructor
  · intro e; grind
  · intro e; rw [e]; grind

/-- A site sitting at the wall −1/ω makes the binding polynomial forget the other site. -/
theorem bind2_wall (ω a : K) (hω : ω ≠ 0) : bind2 ω a (-1 / ω) = 1 - 1 / ω := by
  rw [bind2_is_circK, circK_wall ω a hω]; grind

end Field

end MolecularBiology
