/-
  Fused Hyperoperation Algebra: Planck–Einstein from "the rate depends only on the energy"
  (topics/why-i/rotations-force-i, paper Proposition "Planck–Einstein from the energy zero").
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/PlanckEinstein.lean.

  Assume only that A(H)² = −g(H)² for one function g. With the shift rule, ħA(x·1) = Z + xJ, so
  (Z + xJ)² is a polynomial of degree 2 in x (`rate_sq_quadratic`): g² must be quadratic,
  g(E)² = α(E − E₀)² + m². Measuring energies from E₀ and writing s = √α, J' = J/s:
      J'² = −1,   BJ' + J'B = −2sH,   B² = −(s²H² + m²).
  From these, in any ring with no 2-torsion and s cancellable:
    * J' and B commute with H (energy conservation again);
    * M := B − sJ'H anticommutes with J', commutes with H, and M² = −m².
  So m ≠ 0 needs a second quarter turn M/m anticommuting with J' (quaternion-like); on a level of
  real dimension 2 that is impossible (paper), so m = 0 and A(H) = J'·sH: Planck–Einstein with ħ = 1/s.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.PlanckEinstein
  (namespace CircleOp.PlanckEinstein, the library's convention; no proof changed).
-/
namespace CircleOp.PlanckEinstein
open Lean.Grind

variable {R : Type} [Ring R]

theorem dbl_sub (a b : R) (h : -(a + a) = -(b + b)) : (a - b) + (a - b) = 0 := by grind
theorem eq_of_sub (a b : R) (h : a - b = 0) : a = b := by grind
theorem add_comm' (a b : R) : a + b = b + a := by grind

/-- If a commutes with x and y, it commutes with x*y. -/
theorem comm_mul (a x y : R) (hx : x * a = a * x) (hy : y * a = a * y) :
    x * y * a = a * (x * y) := by
  calc x * y * a = x * (y * a) := Semiring.mul_assoc _ _ _
    _ = x * (a * y) := by rw [hy]
    _ = (x * a) * y := (Semiring.mul_assoc _ _ _).symm
    _ = (a * x) * y := by rw [hx]
    _ = a * (x * y) := Semiring.mul_assoc _ _ _

/-- The squared generator along the shifts is quadratic: (Z + xJ)² = Z² + x(ZJ + JZ) + x²J². -/
theorem rate_sq_quadratic (Z J x : R) (hx : x * J = J * x) (hxZ : x * Z = Z * x) :
    (Z + x * J) * (Z + x * J) = Z * Z + x * (Z * J + J * Z) + x * x * (J * J) := by
  rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
  have h1 : Z * (x * J) = x * (Z * J) := by
    rw [← Semiring.mul_assoc, ← hxZ, Semiring.mul_assoc]
  have h2 : x * J * Z = x * (J * Z) := Semiring.mul_assoc _ _ _
  have h3 : x * J * (x * J) = x * x * (J * J) := by
    rw [Semiring.mul_assoc, ← Semiring.mul_assoc J x J, ← hx, Semiring.mul_assoc,
      ← Semiring.mul_assoc]
  rw [h1, h2, h3, Semiring.left_distrib x]; grind

theorem J_comm_anticomm (J B : R) (hJ : J * J = -1) :
    J * (B * J + J * B) = (B * J + J * B) * J := by
  have e1 : J * (B * J + J * B) = J * B * J + (J * J) * B := by grind
  have e2 : (B * J + J * B) * J = B * (J * J) + J * B * J := by grind
  rw [e1, e2, hJ]; grind

/-- [J', H] = 0 from J'² = −1 and BJ' + J'B = −2sH, with s central and cancellable. -/
theorem J_comm_H (J B H s : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (canc : ∀ x : R, s * x = 0 → x = 0) (sJ : s * J = J * s) (hJ : J * J = -1)
    (hBJ : B * J + J * B = -(s * H + s * H)) : J * H = H * J := by
  have h := J_comm_anticomm J B hJ
  rw [hBJ, Ring.mul_neg, Ring.neg_mul, Semiring.left_distrib, Semiring.right_distrib] at h
  have e1 : J * (s * H) = s * (J * H) := by rw [← Semiring.mul_assoc, ← sJ, Semiring.mul_assoc]
  have e2 : s * H * J = s * (H * J) := Semiring.mul_assoc _ _ _
  rw [e1, e2] at h
  have d := no2 _ (dbl_sub _ _ h)
  have : s * (J * H - H * J) = 0 := by
    rw [Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, ← Ring.sub_eq_add_neg]; exact d
  exact eq_of_sub _ _ (canc _ this)

theorem comm_sq (J H : R) (h : J * H = H * J) : H * H * J = J * (H * H) := by
  calc H * H * J = H * (H * J) := Semiring.mul_assoc _ _ _
    _ = H * (J * H) := by rw [h]
    _ = (H * J) * H := (Semiring.mul_assoc _ _ _).symm
    _ = (J * H) * H := by rw [h]
    _ = J * (H * H) := Semiring.mul_assoc _ _ _

/-- [B, H] = 0 (energy conservation) also with a rest term m²: B² = −(s²H² + m²). -/
theorem B_comm_H (J B H s m2 : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (canc : ∀ x : R, s * x = 0 → x = 0) (sJ : s * J = J * s) (sB : s * B = B * s)
    (mJ : m2 * J = J * m2) (hJ : J * J = -1)
    (hBJ : B * J + J * B = -(s * H + s * H)) (hB : B * B = -(s * s * (H * H) + m2)) :
    B * H = H * B := by
  have hJH := J_comm_H J B H s no2 canc sJ hJ hBJ
  have sq : (s * s * (H * H) + m2) * J = J * (s * s * (H * H) + m2) := by
    rw [Semiring.right_distrib, Semiring.left_distrib, mJ]
    have : s * s * (H * H) * J = J * (s * s * (H * H)) :=
      comm_mul J (s * s) (H * H) (comm_mul J s s sJ sJ) (comm_sq J H hJH)
    rw [this]
  have h : B * (B * J + J * B) = (B * J + J * B) * B := by
    have e1 : B * (B * J + J * B) = (B * B) * J + B * J * B := by grind
    have e2 : (B * J + J * B) * B = B * J * B + J * (B * B) := by grind
    rw [e1, e2, hB, Ring.neg_mul, Ring.mul_neg, sq, add_comm']
  rw [hBJ, Ring.mul_neg, Ring.neg_mul, Semiring.left_distrib, Semiring.right_distrib] at h
  have e1 : B * (s * H) = s * (B * H) := by rw [← Semiring.mul_assoc, ← sB, Semiring.mul_assoc]
  have e2 : s * H * B = s * (H * B) := Semiring.mul_assoc _ _ _
  rw [e1, e2] at h
  have d := no2 _ (dbl_sub _ _ h)
  have : s * (B * H - H * B) = 0 := by
    rw [Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, ← Ring.sub_eq_add_neg]; exact d
  exact eq_of_sub _ _ (canc _ this)

/-- M := B − sJ'H anticommutes with J'. -/
theorem M_anticomm_J (J B H s : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (canc : ∀ x : R, s * x = 0 → x = 0) (sJ : s * J = J * s) (hJ : J * J = -1)
    (hBJ : B * J + J * B = -(s * H + s * H)) :
    J * (B - s * (J * H)) + (B - s * (J * H)) * J = 0 := by
  have hJH := J_comm_H J B H s no2 canc sJ hJ hBJ
  have a : J * (s * (J * H)) = -(s * H) := by
    rw [← Semiring.mul_assoc, ← sJ, Semiring.mul_assoc, ← Semiring.mul_assoc J J H, hJ,
      Ring.neg_mul, Semiring.one_mul, Ring.mul_neg]
  have b : s * (J * H) * J = -(s * H) := by
    rw [Semiring.mul_assoc, Semiring.mul_assoc, ← hJH, ← Semiring.mul_assoc J J H, hJ,
      Ring.neg_mul, Semiring.one_mul, Ring.mul_neg]
  have e : J * (B - s * (J * H)) + (B - s * (J * H)) * J
      = (B * J + J * B) - (J * (s * (J * H)) + s * (J * H) * J) := by
    rw [Ring.sub_eq_add_neg, Ring.sub_eq_add_neg, Semiring.left_distrib, Semiring.right_distrib,
      Ring.mul_neg, Ring.neg_mul]; grind
  rw [e, hBJ, a, b]; grind

/-- M := B − sJ'H squares to −m²: the rest term is a second quarter turn (times m). -/
theorem M_sq (J B H s m2 : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (canc : ∀ x : R, s * x = 0 → x = 0) (sJ : s * J = J * s) (sB : s * B = B * s)
    (sH : s * H = H * s) (mJ : m2 * J = J * m2) (hJ : J * J = -1)
    (hBJ : B * J + J * B = -(s * H + s * H)) (hB : B * B = -(s * s * (H * H) + m2)) :
    (B - s * (J * H)) * (B - s * (J * H)) = -m2 := by
  have hJH := J_comm_H J B H s no2 canc sJ hJ hBJ
  have hBH := B_comm_H J B H s m2 no2 canc sJ sB mJ hJ hBJ hB
  have hMJ := M_anticomm_J J B H s no2 canc sJ hJ hBJ
  generalize hM : B - s * (J * H) = M at hMJ ⊢
  have hBdef : B = s * (J * H) + M := by rw [← hM]; grind
  have hMH : M * H = H * M := by
    rw [← hM, Ring.sub_eq_add_neg, Semiring.right_distrib, Semiring.left_distrib, Ring.neg_mul,
      Ring.mul_neg, hBH]
    have : s * (J * H) * H = H * (s * (J * H)) :=
      comm_mul H s (J * H) sH (comm_mul H J H hJH rfl)
    rw [this]
  have sM : s * M = M * s := by
    rw [← hM, Ring.sub_eq_add_neg, Semiring.left_distrib, Semiring.right_distrib, sB, Ring.mul_neg,
      Ring.neg_mul]
    have : s * (s * (J * H)) = s * (J * H) * s := by
      rw [Semiring.mul_assoc s (J * H) s, Semiring.mul_assoc J H s, ← sH, ← Semiring.mul_assoc J s H,
        ← sJ, Semiring.mul_assoc]
    rw [this]
  have K2 : s * (J * H) * (s * (J * H)) = -(s * s * (H * H)) := by
    have q : J * H * (J * H) = -(H * H) := by
      rw [Semiring.mul_assoc, ← Semiring.mul_assoc H J H, ← hJH, Semiring.mul_assoc,
        ← Semiring.mul_assoc J J, hJ, Ring.neg_mul, Semiring.one_mul]
    have c : J * H * s = s * (J * H) := comm_mul s J H sJ.symm sH.symm
    calc s * (J * H) * (s * (J * H)) = s * ((J * H) * (s * (J * H))) := Semiring.mul_assoc _ _ _
      _ = s * ((J * H * s) * (J * H)) := by rw [← Semiring.mul_assoc (J * H) s]
      _ = s * ((s * (J * H)) * (J * H)) := by rw [c]
      _ = s * (s * ((J * H) * (J * H))) := by rw [Semiring.mul_assoc s (J * H)]
      _ = -(s * s * (H * H)) := by rw [q, Ring.mul_neg, Ring.mul_neg, ← Semiring.mul_assoc]
  have EJ : M * J = -(J * M) := by grind
  have cross : s * (J * H) * M + M * (s * (J * H)) = 0 := by
    have : M * (s * (J * H)) = -(s * (J * H) * M) := by
      rw [← Semiring.mul_assoc, ← sM, Semiring.mul_assoc s M, ← Semiring.mul_assoc M J H, EJ,
        Ring.neg_mul, Ring.mul_neg, Semiring.mul_assoc J M H, hMH, ← Semiring.mul_assoc J H M,
        ← Semiring.mul_assoc]
    rw [this]; grind
  have expand : B * B = s * (J * H) * (s * (J * H)) + (s * (J * H) * M + M * (s * (J * H))) + M * M := by
    rw [hBdef, Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]; grind
  rw [hB, K2, cross] at expand
  grind

/-- Where the second quarter turn is excluded (M = 0), the rest term vanishes and
    B = s·J'H: the rate is s·(E − E₀), Planck–Einstein with ħ = 1/s. -/
theorem planck_einstein (J B H s m2 : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (canc : ∀ x : R, s * x = 0 → x = 0) (sJ : s * J = J * s) (sB : s * B = B * s)
    (sH : s * H = H * s) (mJ : m2 * J = J * m2) (hJ : J * J = -1)
    (hBJ : B * J + J * B = -(s * H + s * H)) (hB : B * B = -(s * s * (H * H) + m2))
    (h0 : B - s * (J * H) = 0) : m2 = 0 ∧ B = s * (J * H) := by
  have h := M_sq J B H s m2 no2 canc sJ sB sH mJ hJ hBJ hB
  rw [h0, Semiring.mul_zero] at h
  constructor
  · grind
  · exact eq_of_sub _ _ h0

end CircleOp.PlanckEinstein
