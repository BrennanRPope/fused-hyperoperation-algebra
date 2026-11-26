/-
  Face-power dial (DRAFT): one roundness dial aimed at any polytope's faces.

  A polytope with m faces is given by its scaled face normals N i (facets at N i · x = 1, centre at 0)
  and its offset a (a = 0 for centrally symmetric shapes; a = (d - 1)/2 for the d-simplex).
  The shape at whole-number dial n is
      shape N a n = { x | Σ_i |a + N i · x|^n ≤ m (a + 1)^n }
  i.e. the power mean of the shifted face values is at most a + 1.

  Proved here, with no `sorry`:
    * `poly_subset_shape`: the polytope sits inside every shape of the dial.
    * `face_le_of_mem_shape`: inside the shape at n, every face value is at most m^(1/n)(a + 1)
      (stated without roots: |a + N i · x|^n ≤ m (a + 1)^n). This is the 1/n sharpening rate.
    * `iInter_shape`: the polytope is exactly what all the shapes share (the dial's sharp end).
    * `shape_two`: at n = 2 the shape is a ball, when the normals sum to zero and form a tight frame
      (true for every regular polytope at its centre, and for any polytope at its analytic centre after
      a linear change of coordinates).
  Numerical companions: shape-dials/checks.py (29 checks).
-/
import Mathlib

namespace MathlibMap.FaceDial

open Finset

variable {m d : ℕ}

/-- The value of face `i` at `x`: `N i · x`. -/
def face (N : Fin m → Fin d → ℝ) (i : Fin m) (x : Fin d → ℝ) : ℝ := ∑ j, N i j * x j

/-- The polytope, written with the offset: every `|a + N i · x| ≤ a + 1`. -/
def poly (N : Fin m → Fin d → ℝ) (a : ℝ) : Set (Fin d → ℝ) :=
  {x | ∀ i, |a + face N i x| ≤ a + 1}

/-- The face-power shape at whole-number dial `n`. -/
def shape (N : Fin m → Fin d → ℝ) (a : ℝ) (n : ℕ) : Set (Fin d → ℝ) :=
  {x | ∑ i, |a + face N i x| ^ n ≤ (m : ℝ) * (a + 1) ^ n}

theorem poly_subset_shape (N : Fin m → Fin d → ℝ) {a : ℝ} (ha : 0 ≤ a + 1) (n : ℕ) :
    poly N a ⊆ shape N a n := by
  intro x hx
  simp only [poly, Set.mem_setOf_eq] at hx
  simp only [shape, Set.mem_setOf_eq]
  have hterm : ∀ i ∈ (univ : Finset (Fin m)), |a + face N i x| ^ n ≤ (a + 1) ^ n := by
    intro i _
    have h0 : 0 ≤ |a + face N i x| := abs_nonneg _
    have h1 := hx i
    gcongr
  calc ∑ i, |a + face N i x| ^ n ≤ ∑ _i : Fin m, (a + 1) ^ n := sum_le_sum hterm
    _ = (m : ℝ) * (a + 1) ^ n := by simp

theorem face_le_of_mem_shape (N : Fin m → Fin d → ℝ) {a : ℝ} {n : ℕ} {x : Fin d → ℝ}
    (hx : ∑ i, |a + face N i x| ^ n ≤ (m : ℝ) * (a + 1) ^ n) (i : Fin m) :
    |a + face N i x| ^ n ≤ (m : ℝ) * (a + 1) ^ n := by
  have hle : |a + face N i x| ^ n ≤ ∑ j, |a + face N j x| ^ n :=
    single_le_sum (f := fun j => |a + face N j x| ^ n)
      (fun j _ => pow_nonneg (abs_nonneg _) n) (mem_univ i)
  exact le_trans hle hx

/-- The sharp end of the dial: the polytope is exactly the intersection of all the shapes. -/
theorem iInter_shape (N : Fin m → Fin d → ℝ) {a : ℝ} (ha : 0 < a + 1) :
    (⋂ n : ℕ, shape N a n) = poly N a := by
  ext x
  simp only [Set.mem_iInter, shape, poly, Set.mem_setOf_eq]
  constructor
  · intro hx i
    by_contra h
    push_neg at h
    have hq : 1 < |a + face N i x| / (a + 1) := by
      rw [one_lt_div ha]; exact h
    obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt (m : ℝ) hq
    have h1 := face_le_of_mem_shape N (hx n) i
    have hpos : 0 < (a + 1) ^ n := pow_pos ha n
    have hne : a + 1 ≠ 0 := ne_of_gt ha
    have hsplit : (|a + face N i x| / (a + 1)) ^ n * (a + 1) ^ n = |a + face N i x| ^ n := by
      rw [← mul_pow]
      congr 1
      field_simp
    rw [← hsplit] at h1
    have h2 : (|a + face N i x| / (a + 1)) ^ n ≤ (m : ℝ) := le_of_mul_le_mul_right h1 hpos
    linarith
  · intro hx n
    have hmem : x ∈ poly N a := by
      simp only [poly, Set.mem_setOf_eq]; exact hx
    have := poly_subset_shape N (le_of_lt ha) n hmem
    simpa only [shape, Set.mem_setOf_eq] using this

/-- The round point of the dial: at `n = 2` the shape is a ball, when the normals sum to zero and
form a tight frame (`Σ_i (N i · x)² = c ‖x‖²`). -/
theorem shape_two (N : Fin m → Fin d → ℝ) (a c : ℝ) (hsum : ∀ j, ∑ i, N i j = 0)
    (htight : ∀ x : Fin d → ℝ, ∑ i, face N i x ^ 2 = c * ∑ j, x j ^ 2) :
    shape N a 2 = {x | c * ∑ j, x j ^ 2 ≤ (m : ℝ) * (2 * a + 1)} := by
  ext x
  simp only [shape, Set.mem_setOf_eq, sq_abs]
  have hlin : ∑ i, face N i x = 0 := by
    simp only [face]
    rw [Finset.sum_comm]
    simp [← Finset.sum_mul, hsum]
  have expand : ∑ i, (a + face N i x) ^ 2
      = (m : ℝ) * a ^ 2 + 2 * a * ∑ i, face N i x + ∑ i, face N i x ^ 2 := by
    simp only [add_sq, Finset.sum_add_distrib, Finset.mul_sum, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul] <;> ring
  rw [expand, hlin, htight x]
  constructor <;> intro h <;> nlinarith [h]

end MathlibMap.FaceDial
