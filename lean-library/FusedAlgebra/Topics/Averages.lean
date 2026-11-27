/-
  Fused Hyperoperation Algebra: the shift rule (A4) from "shifting the energy zero changes no
  measured average" (topics/why-i/rotations-force-i, paper Proposition "(A4) from measured averages").
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/Averages.lean.

  R is any ring (not necessarily commutative), e.g. real n×n matrices; S is any commutative ring of
  scalars, e.g. the reals, mapped into the centre of R by ι (x ↦ x·1). `Vt` stands for Vᵀ.
    * `comm_of_invariant`: if V is invertible with inverse Vᵀ and VᵀOV = O, then OV = VO.
    * `shift_coeff`: B² = −H² and (B + βJ)² = −(H + c)² give β(BJ + JB) = β² − c² − 2cH.
    * `no_dynamics_without_phase`: if the shift adds nothing, (A3) at H, H + 1, H − 1 forces 1 = 0
      in a ring without 2-torsion: with no antisymmetric matrix commuting with every observable,
      no dynamics satisfies (A1), (A3) and (A0′).
    * `beta_ne_zero`: the unit shift adds a nonzero multiple of J when 1 and H are independent.
    * `sigma_from_shifts`: the shifts c = 1 and c = 2 alone, with 1 and H independent over S, force
      β(H, 2) = 2β(H, 1) and σ := β(H, 1) to satisfy σ² = 1, and give BJ + JB = −2σH.
    * `flipped_unit`: J′ := σJ then has J′² = −1 and BJ′ + J′B = −2H. These are the equations of
      `EnergyShift.lean`, whose theorems then give [J′, H] = [B, H] = 0 and (B − J′H)² = 0.
  Not covered here (left to the paper): the step from equal averages to VᵀOV = O (it uses that both
  sides are symmetric), differentiating V_t at t = 0, the use of (T₀) to make every shift a
  multiple of one J, and the continuity and connectedness argument that makes σ(H) the same for
  every H.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.Averages
  (namespace CircleOp.Averages, the library's convention; no proof changed).
-/
namespace CircleOp.Averages
open Lean.Grind

/-- Scalars S mapped into the centre of R, as x ↦ x·1 maps the reals into the real matrices. -/
structure Scalars (S R : Type) [CommRing S] [Ring R] where
  ι : S → R
  map_add : ∀ a b : S, ι (a + b) = ι a + ι b
  map_mul : ∀ a b : S, ι (a * b) = ι a * ι b
  map_one : ι 1 = 1
  map_neg : ∀ a : S, ι (-a) = -ι a
  central : ∀ (a : S) (x : R), ι a * x = x * ι a

variable {S R : Type} [CommRing S] [Ring R]

theorem map_zero (k : Scalars S R) : k.ι 0 = 0 := by
  have h := k.map_add 0 0
  have : (0 : S) + 0 = 0 := by grind
  rw [this] at h; grind

theorem map_sub (k : Scalars S R) (a b : S) : k.ι (a - b) = k.ι a - k.ι b := by
  rw [Ring.sub_eq_add_neg, k.map_add, k.map_neg, ← Ring.sub_eq_add_neg]

theorem smul_smul (k : Scalars S R) (a b : S) (x : R) :
    k.ι a * (k.ι b * x) = k.ι (a * b) * x := by
  rw [k.map_mul, Semiring.mul_assoc]

/-- An invertible V with inverse Vᵀ (VVᵀ = 1) and VᵀOV = O commutes with O. -/
theorem comm_of_invariant (V Vt O : R) (hr : V * Vt = 1)
    (h : Vt * O * V = O) : O * V = V * O := by
  have e : V * (Vt * O * V) = V * O := by rw [h]
  rw [← Semiring.mul_assoc, ← Semiring.mul_assoc, hr, Semiring.one_mul] at e
  exact e

/-- (A3) at H and at H + c, with the shift ħA(H + c) = B + βJ:
    β(BJ + JB) = (β² − c²) − 2cH. -/
theorem shift_coeff (k : Scalars S R) (J B H : R) (β c : S) (hJ : J * J = -1)
    (hB : B * B = -(H * H))
    (hs : (B + k.ι β * J) * (B + k.ι β * J) = -((H + k.ι c) * (H + k.ι c))) :
    k.ι β * (B * J + J * B) = k.ι (β * β - c * c) - k.ι (c + c) * H := by
  have hβ := k.central β
  have hc := k.central c
  have e1 : (B + k.ι β * J) * (B + k.ι β * J)
      = B * B + k.ι β * (B * J + J * B) + k.ι β * k.ι β * (J * J) := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
    have a1 : B * (k.ι β * J) = k.ι β * (B * J) := by
      rw [← Semiring.mul_assoc, ← hβ B, Semiring.mul_assoc]
    have a2 : k.ι β * J * B = k.ι β * (J * B) := Semiring.mul_assoc _ _ _
    have a3 : k.ι β * J * (k.ι β * J) = k.ι β * k.ι β * (J * J) := by
      rw [Semiring.mul_assoc, ← Semiring.mul_assoc J (k.ι β) J, ← hβ J, Semiring.mul_assoc,
        ← Semiring.mul_assoc]
    rw [a1, a2, a3, Semiring.left_distrib (k.ι β)]; grind
  have e2 : (H + k.ι c) * (H + k.ι c) = H * H + k.ι (c + c) * H + k.ι c * k.ι c := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib, ← hc H,
      k.map_add, Semiring.right_distrib]; grind
  rw [e1, hB, hJ, e2, Ring.mul_neg, Semiring.mul_one] at hs
  rw [map_sub, k.map_mul, k.map_mul]
  grind

/-- If the shift adds nothing, B² = −H², B² = −(H + 1)² and B² = −(H − 1)² give 1 + 1 = 0. -/
theorem shift_nonzero (B H : R) (h0 : B * B = -(H * H))
    (hp : B * B = -((H + 1) * (H + 1))) (hm : B * B = -((H - 1) * (H - 1))) :
    (1 : R) + 1 = 0 := by
  rw [h0] at hp hm
  have ep : (H + 1) * (H + 1) = H * H + (H + H) + 1 := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib,
      Semiring.mul_one, Semiring.one_mul, Semiring.mul_one]; grind
  have em : (H - 1) * (H - 1) = H * H - (H + H) + 1 := by
    rw [Ring.sub_eq_add_neg, Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib,
      Ring.mul_neg, Ring.neg_mul, Ring.mul_neg, Semiring.mul_one, Semiring.one_mul,
      Semiring.mul_one]; grind
  rw [ep] at hp; rw [em] at hm
  grind

/-- So in a ring without 2-torsion where the shift adds nothing, 1 = 0: no dynamics satisfies
    (A1), (A3) and (A0′) when no antisymmetric matrix commutes with every observable. -/
theorem no_dynamics_without_phase (B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (h0 : B * B = -(H * H))
    (hp : B * B = -((H + 1) * (H + 1))) (hm : B * B = -((H - 1) * (H - 1))) :
    (1 : R) = 0 :=
  no2 1 (shift_nonzero B H h0 hp hm)

/-- Rearranging a − cH = b − dH. -/
theorem lin (a b c d H : R) (h : a - c * H = b - d * H) : (a - b) + (d - c) * H = 0 := by
  rw [Ring.sub_eq_add_neg d, Semiring.right_distrib, Ring.neg_mul]
  rw [Ring.sub_eq_add_neg, Ring.sub_eq_add_neg b] at h
  grind

/-- 1 and H are independent over the scalars. -/
def Indep (k : Scalars S R) (H : R) : Prop :=
  ∀ x y : S, k.ι x + k.ι y * H = 0 → x = 0 ∧ y = 0

/-- The unit shift is a nonzero multiple of J: β(H, 1) = 0 would give 1 + 2H = 0. -/
theorem beta_ne_zero (k : Scalars S R) (J B H : R) (hJ : J * J = -1) (hind : Indep k H)
    (hB : B * B = -(H * H))
    (hs : (B + k.ι 0 * J) * (B + k.ι 0 * J) = -((H + k.ι 1) * (H + k.ι 1))) :
    (1 : S) = 0 := by
  have h := shift_coeff k J B H 0 1 hJ hB hs
  rw [map_zero, Semiring.zero_mul, map_sub, k.map_mul, k.map_mul, map_zero, Semiring.zero_mul,
    k.map_one, Semiring.one_mul] at h
  have h' : k.ι 1 + k.ι (1 + 1) * H = 0 := by rw [k.map_one]; grind
  exact (hind 1 (1 + 1) h').1

/-- The shifts c = 1 and c = 2 fix the sign. With b₁ = β(H, 1), b₂ = β(H, 2), 1 and H independent
    over S, no 2-torsion in S and b₁ invertible: b₂ = 2b₁, b₁² = 1, and BJ + JB = −2b₁H. -/
theorem sigma_from_shifts (k : Scalars S R) (J B H : R) (b1 b2 w : S)
    (no2 : ∀ x : S, x + x = 0 → x = 0) (hinv : w * b1 = 1)
    (hJ : J * J = -1) (hind : Indep k H)
    (hB : B * B = -(H * H))
    (h1 : (B + k.ι b1 * J) * (B + k.ι b1 * J) = -((H + k.ι 1) * (H + k.ι 1)))
    (h2 : (B + k.ι b2 * J) * (B + k.ι b2 * J) = -((H + k.ι (1 + 1)) * (H + k.ι (1 + 1)))) :
    b2 = b1 + b1 ∧ b1 * b1 = 1 ∧ B * J + J * B = -(k.ι (b1 + b1) * H) := by
  have e1 := shift_coeff k J B H b1 1 hJ hB h1
  have e2 := shift_coeff k J B H b2 (1 + 1) hJ hB h2
  generalize hSm : B * J + J * B = Sm at e1 e2
  -- b₂·(e1) and b₁·(e2) have the same left side b₁b₂S.
  have m1 : k.ι b2 * (k.ι b1 * Sm)
      = k.ι (b2 * (b1 * b1 - 1 * 1)) - k.ι (b2 * (1 + 1)) * H := by
    rw [e1, Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, smul_smul, k.map_mul,
      ← Ring.sub_eq_add_neg, ← k.map_mul]
  have m2 : k.ι b1 * (k.ι b2 * Sm)
      = k.ι (b1 * (b2 * b2 - (1 + 1) * (1 + 1))) - k.ι (b1 * ((1 + 1) + (1 + 1))) * H := by
    rw [e2, Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, smul_smul, k.map_mul,
      ← Ring.sub_eq_add_neg, ← k.map_mul]
  have swap : k.ι b2 * (k.ι b1 * Sm) = k.ι b1 * (k.ι b2 * Sm) := by
    rw [smul_smul, smul_smul]; congr 2; grind
  rw [m1, m2] at swap
  -- So x + yH = 0 with x, y the scalars below.
  have hx : k.ι (b2 * (b1 * b1 - 1 * 1) - b1 * (b2 * b2 - (1 + 1) * (1 + 1)))
      + k.ι (b1 * ((1 + 1) + (1 + 1)) - b2 * (1 + 1)) * H = 0 := by
    rw [map_sub, map_sub]; exact lin _ _ _ _ H swap
  have ⟨ix, iy⟩ := hind _ _ hx
  have hb2 : b2 = b1 + b1 := by
    have : (b1 + b1 - b2) + (b1 + b1 - b2) = 0 := by grind
    have := no2 _ this; grind
  have hsq : b1 * b1 = 1 := by
    rw [hb2] at ix
    have z : b1 * (b1 * b1 - 1) = 0 := no2 _ (by grind)
    have : w * (b1 * (b1 * b1 - 1)) = 0 := by rw [z]; grind
    have : b1 * b1 - 1 = 0 := by
      have e : w * (b1 * (b1 * b1 - 1)) = (w * b1) * (b1 * b1 - 1) := by grind
      rw [e, hinv] at this; grind
    grind
  refine ⟨hb2, hsq, ?_⟩
  -- b₁S = (b₁² − 1) − 2H = −2H; multiply by b₁.
  have e3 : k.ι b1 * (k.ι b1 * Sm) = k.ι b1 * (-(k.ι (1 + 1) * H)) := by
    rw [e1]
    have : b1 * b1 - 1 * 1 = 0 := by grind
    rw [this, map_zero]; grind
  rw [smul_smul, hsq, k.map_one, Semiring.one_mul, Ring.mul_neg, smul_smul] at e3
  rw [e3]; congr 3; grind

/-- With J′ = σJ and σ² = 1, J′² = −1 and BJ′ + J′B = −2H: Eq. (coeffs) for J′. -/
theorem flipped_unit (k : Scalars S R) (J B H : R) (s : S) (hJ : J * J = -1)
    (hss : s * s = 1) (hS : B * J + J * B = -(k.ι (s + s) * H)) :
    (k.ι s * J) * (k.ι s * J) = -1 ∧ B * (k.ι s * J) + (k.ι s * J) * B = -(H + H) := by
  have hs := k.central s
  constructor
  · rw [Semiring.mul_assoc, ← Semiring.mul_assoc J (k.ι s) J, ← hs J, Semiring.mul_assoc,
      ← Semiring.mul_assoc, ← k.map_mul, hss, k.map_one, hJ, Semiring.one_mul]
  · have : B * (k.ι s * J) + (k.ι s * J) * B = k.ι s * (B * J + J * B) := by
      rw [← Semiring.mul_assoc, ← hs B, Semiring.mul_assoc, Semiring.mul_assoc,
        Semiring.left_distrib]
    rw [this, hS, Ring.mul_neg, smul_smul]
    have : s * (s + s) = 1 + 1 := by grind
    rw [this, k.map_add, k.map_one, Semiring.right_distrib, Semiring.one_mul]

end CircleOp.Averages
