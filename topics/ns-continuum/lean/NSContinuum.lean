/-
  NSContinuum: the dimension continuum and the 3D Navier–Stokes gap
  ==================================================================

  The exact statements behind the "dimension continuum" page (Nov 25, 2026).
  Written for core Lean 4.37.0-rc1 in the style of the shared library
  (FusedAlgebra), no mathlib, no `sorry`.

  STATUS: NOT COMPILED. No Lean toolchain was available in the session that
  wrote this file (installing one needs Brennan's yes). To check it, copy it to
  lean-library/FusedAlgebra/Topics/NSContinuum.lean, add it to build.sh and
  FusedAlgebra.lean, and run ./build.sh.

  What is here:
  1. The dyadic shell is the mult-add cube: doubling a box in n dimensions adds
     2^(nj) · 1^(◯n) new cells, and those new cells fall into the cube's layers
     1..n (how many coordinates doubled).
  2. The fan: every scaling threshold is a ray from the pivot (n, α) = (−2, 0);
     lying on the same ray is an equivalence, so sliding the dimension and
     sliding the dissipation are one dial. 3D NS sits on the ray of 2D NS with
     dissipation power 4/5.
  3. The K41 edge of the dyadic model: the constant-flux margin 6α − (n+2)
     vanishes at n = 4 when α = 1.
  4. The dyadic energy identity: the nonlinearity only passes energy upward;
     the only energy leaving a truncated model is the flux out of the top.
  5. Planar flow has an empty stretching corner: ω·∇u = 0 identically.

  What is NOT here (needs real analysis): Osgood integrals, the rung height of
  the dissipation deficiency, anything about actual Navier–Stokes solutions.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube
import FusedAlgebra.Topics.Topology.CubeProofs

namespace CircleOp.NSContinuum

open CircleOp CircleOp.Cube

/-! ## 1. The dyadic shell is the mult-add cube -/

/-- Mult-adding 1 to itself n times counts the nonempty corners of the n-cube:
    1^(◯n) + 1 = 2ⁿ. -/
theorem mpow_one_nat (n : Nat) : mpow (1 : Nat) n + 1 = 2 ^ n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show mpow (1 : Nat) n + 1 + mpow (1 : Nat) n * 1 + 1 = 2 ^ (n + 1)
    rw [Nat.pow_succ, ← ih, Nat.mul_one]
    omega

/-- Doubling the side of an n-dimensional box of side 2ʲ multiplies its cells by 2ⁿ. -/
theorem box_doubles (n j : Nat) : 2 ^ (n * (j + 1)) = 2 ^ (n * j) * 2 ^ n := by
  rw [Nat.mul_succ, Nat.pow_add]

/-- **The dyadic annulus is a mult-add power.** The new cells added when a box of
    side 2ʲ doubles number 2^(nj) · 1^(◯n): the old box times the n-fold
    mult-add of 1. -/
theorem annulus_cells (n j : Nat) :
    2 ^ (n * (j + 1)) - 2 ^ (n * j) = 2 ^ (n * j) * mpow (1 : Nat) n := by
  have h := mpow_one_nat n
  have e : mpow (1 : Nat) n = 2 ^ n - 1 := by omega
  rw [box_doubles, e, Nat.mul_sub_one]

/-- **The annulus is graded by the cube's layers.** 2ⁿ − 1 = Σ_{k=1..n} C(n,k):
    the new cells of a doubled box split by how many of the n coordinates landed
    in the upper half (k = 1 strips, k = 2 corners, …, k = n the far corner). -/
theorem annulus_graded (n : Nat) :
    (2 : Int) ^ n - 1 = bsum (fun k => (layer n (k + 1) : Int) * 1 ^ (k + 1)) n := by
  have hb := binomial n 1
  have hc : cubeSum n 1 = 1 + bsum (fun k => (layer n (k + 1) : Int) * 1 ^ (k + 1)) n := by
    unfold cubeSum; rw [bsum_shift]; simp [layer_zero_right]
  have two : (1 + 1 : Int) = 2 := rfl
  rw [two] at hb
  rw [hb, hc]
  omega

/-! ## 2. The fan of scaling rays

  Dissipation is measured in twentieths: A = 20·α, so 3D Navier–Stokes is
  (n, A) = (3, 20), Lions' threshold in 3D is (3, 25), 2D Navier–Stokes is (2, 20).
  Two points lie on the same scaling ray when A·(n'+2) = A'·(n+2): the ratio
  γ = 2α/(n+2) (dissipation power over the nonlinearity's growth) is the same. -/

/-- Same scaling ray through the pivot (n, A) = (−2, 0). -/
def onRay (n A n' A' : Int) : Prop := A * (n' + 2) = A' * (n + 2)

theorem onRay_refl (n A : Int) : onRay n A n A := rfl

theorem onRay_symm {n A n' A' : Int} (h : onRay n A n' A') : onRay n' A' n A := by
  unfold onRay at *; exact h.symm

/-- Lying on one ray is transitive (away from the pivot dimension n = −2), so the
    rays split the plane: sliding the dimension and sliding the dissipation are
    one dial. -/
theorem onRay_trans {n A n' A' n'' A'' : Int} (hne : n' + 2 ≠ 0)
    (h1 : onRay n A n' A') (h2 : onRay n' A' n'' A'') : onRay n A n'' A'' := by
  unfold onRay at *
  have h : (A * (n'' + 2)) * (n' + 2) = (A'' * (n + 2)) * (n' + 2) := by
    calc (A * (n'' + 2)) * (n' + 2) = (A * (n' + 2)) * (n'' + 2) := by grind
      _ = (A' * (n + 2)) * (n'' + 2) := by rw [h1]
      _ = (A' * (n'' + 2)) * (n + 2) := by grind
      _ = (A'' * (n' + 2)) * (n + 2) := by rw [h2]
      _ = (A'' * (n + 2)) * (n' + 2) := by grind
  exact Int.eq_of_mul_eq_mul_right hne h

/-- Lions' 3D threshold α = 5/4 is on the same ray as 2D Navier–Stokes (α = 1):
    both are critical. -/
theorem lions_on_2d_ray : onRay 3 25 2 20 := by unfold onRay; decide

/-- 3D Navier–Stokes (α = 1) is on the same ray as 2D flow with dissipation power
    α = 4/5. (2D is smooth there anyway, for a reason scaling cannot see: §5.) -/
theorem ns3d_on_ray_of_2d_four_fifths : onRay 3 20 2 16 := by unfold onRay; decide

/-- 3D Navier–Stokes is off the critical ray: the scaling gap is a quarter power. -/
theorem ns3d_supercritical : ¬ onRay 3 20 3 25 := by unfold onRay; decide

/-- The three known thresholds are rays: critical α = (n+2)/4, the dyadic
    regularity line α = (n+2)/5 (Barbato–Morandin–Romito) and the dyadic blow-up
    line α = (n+2)/6 (Cheskidov). In twentieths: 5(n+2), 4(n+2), (10/3)(n+2);
    all three vanish at the pivot n = −2. -/
theorem fan_pivot : 5 * ((-2 : Int) + 2) = 0 ∧ 4 * ((-2 : Int) + 2) = 0 ∧ 10 * ((-2 : Int) + 2) = 0 := by
  decide

/-- 3D Navier–Stokes lies exactly on the dyadic regularity line A = 4(n+2). -/
theorem ns3d_on_bmr_line : (20 : Int) = 4 * (3 + 2) := by decide

/-! ## 3. The K41 edge of the dyadic model

  Constant energy flux λʲ aⱼ³ = const gives aⱼ ∝ λ^(−j/3). Dissipation beats
  transfer at small scales iff 2γ > 2/3, i.e. 6α > n + 2. -/

/-- The constant-flux margin 6α − (n+2) (α = 1 here). -/
def k41Margin (n : Int) : Int := 6 - (n + 2)

theorem k41_edge_at_four : k41Margin 4 = 0 := by decide
theorem k41_3d_positive : 0 < k41Margin 3 := by decide
theorem k41_5d_negative : k41Margin 5 < 0 := by decide

/-! ## 4. The dyadic energy identity -/

/-- A finite sum Σ_{k<m} f k. -/
def ssum (f : Nat → Int) : Nat → Int
  | 0 => 0
  | m + 1 => ssum f m + f m

/-- Energy flux from shell j to shell j+1 in the dyadic model
    aⱼ' = c(j−1)·a(j−1)² − c(j)·aⱼ·a(j+1) (times aⱼ, for the energy). -/
def flux (c a : Nat → Int) (j : Nat) : Int := c j * a j * a j * a (j + 1)

/-- Shell j's energy gain from the nonlinearity. -/
def gain (c a : Nat → Int) : Nat → Int
  | 0 => -(flux c a 0)
  | j + 1 => flux c a j - flux c a (j + 1)

/-- The gain is aⱼ times the model's nonlinear term. -/
theorem gain_is_model (c a : Nat → Int) (j : Nat) :
    gain c a (j + 1) = a (j + 1) * (c j * a j ^ 2 - c (j + 1) * a (j + 1) * a (j + 2)) := by
  simp only [gain, flux]; grind

/-- **Energy only moves.** Summed over shells 0..N, the nonlinear gains telescope
    to minus the flux out of the top shell. -/
theorem energy_telescopes (c a : Nat → Int) : ∀ N, ssum (gain c a) (N + 1) = -(flux c a N)
  | 0 => by simp [ssum, gain]
  | N + 1 => by
    show ssum (gain c a) (N + 1) + gain c a (N + 1) = -(flux c a (N + 1))
    rw [energy_telescopes c a N]
    simp only [gain]
    omega

/-! ## 5. The empty stretching corner of planar flow

  Vortex stretching is (ω·∇)u, the i-th component Σⱼ ωⱼ ∂ⱼuᵢ. For planar flow,
  u = (u₁(x,y), u₂(x,y), 0) has ω = (0, 0, w) and nothing depends on z, so the
  only term left, w·∂_z uᵢ, is zero: the corner that feeds 3D growth is empty. -/

/-- Stretching with velocity gradient G (G i j = ∂ⱼuᵢ) and vorticity ω. -/
def stretch (G : Fin 3 → Fin 3 → Int) (ω : Fin 3 → Int) (i : Fin 3) : Int :=
  G i 0 * ω 0 + G i 1 * ω 1 + G i 2 * ω 2

theorem planar_stretching_zero (G : Fin 3 → Fin 3 → Int) (w : Int)
    (hz : ∀ i, G i 2 = 0) (i : Fin 3) :
    stretch G (fun k => if k = 2 then w else 0) i = 0 := by
  simp [stretch, hz]

end CircleOp.NSContinuum
