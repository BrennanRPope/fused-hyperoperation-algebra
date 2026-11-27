import Mathlib

/-!
# Klein bottle × circle inside the periodic box (DRAFT)

The glide g(x, y, z) = (x + π, -y, z) glues the Klein bottle × circle. Two small facts used in
`topics/klein-flow/periodic-box/README.md`:

* gliding twice is the plain shift by 2π along x, so the Klein bottle box unfolds into the periodic box;
* the glide's linear part has determinant -1, so vorticity (a pseudovector, carried by det • D) comes
  back from the twist reversed along the x and z directions and kept along y: a vortex tube's partner
  across the twist spins the other way.
-/

namespace MathlibMap.KleinBox

/-- The glide that glues the Klein bottle × circle. -/
noncomputable def glide (p : ℝ × ℝ × ℝ) : ℝ × ℝ × ℝ := (p.1 + Real.pi, -p.2.1, p.2.2)

/-- Gliding twice is the shift by 2π along x. -/
theorem glide_glide (x y z : ℝ) : glide (glide (x, y, z)) = (x + 2 * Real.pi, y, z) := by
  have h : x + Real.pi + Real.pi = x + 2 * Real.pi := by ring
  show (x + Real.pi + Real.pi, -(-y), z) = (x + 2 * Real.pi, y, z)
  rw [neg_neg, h]

/-- The glide's linear part, diag(1, -1, 1). -/
def glideLin : Matrix (Fin 3) (Fin 3) ℝ := Matrix.diagonal ![1, -1, 1]

/-- The glide reverses orientation. -/
theorem glideLin_det : glideLin.det = -1 := by
  simp [glideLin, Matrix.det_diagonal, Fin.prod_univ_three]

/-- A vortex tube pointing along z returns from the twist pointing the other way. -/
theorem tube_z_flips : glideLin.det • glideLin.mulVec ![0, 0, 1] = ![0, 0, -1] := by
  rw [glideLin_det]
  ext i
  fin_cases i <;> simp [glideLin, Matrix.mulVec_diagonal]

/-- So does a tube pointing along the glide direction x. -/
theorem tube_x_flips : glideLin.det • glideLin.mulVec ![1, 0, 0] = ![-1, 0, 0] := by
  rw [glideLin_det]
  ext i
  fin_cases i <;> simp [glideLin, Matrix.mulVec_diagonal]

/-- A tube pointing along the flipped direction y keeps its sense. -/
theorem tube_y_keeps : glideLin.det • glideLin.mulVec ![0, 1, 0] = ![0, 1, 0] := by
  rw [glideLin_det]
  ext i
  fin_cases i <;> simp [glideLin, Matrix.mulVec_diagonal]

end MathlibMap.KleinBox
