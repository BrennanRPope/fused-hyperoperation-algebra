/-
  Fused Hyperoperation Algebra: Planck–Einstein from "shifting the energy zero changes no measured
  average" (topics/why-i/rotations-force-i, paper Proposition "Planck–Einstein from measured
  averages"). Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/AveragesPE.lean.

  R is any ring (not necessarily commutative), e.g. real n×n matrices. Scalars are central elements.
  Only (A3′), A(H)² = −h(H), is assumed of the rate, with h(H) a matrix commuting with everything
  it meets below. Algebraic steps of the proof:
    * `second_deriv`: with X = Y + D, X² − 2XY + Y² = [Y, D] + D² (the t² term of e^{tX}e^{−tY},
      times 2).
    * `trace_comm_zero`, `theta_zero`: for any trace-like map, tr([A, K]K) = 0, so [A, K] = θK
      with tr(K²) cancellable forces θ = 0.
    * `scalar_line_expand`, `scalar_line_solve`: with [Z, K] = 0, (A3′) on the scalar line gives
      2φZK + φ²K² = h(0) − h(x); two distinct nonzero values of φ make K² and ZK scalars.
    * `Z_of_ZJ`: Z = −(ZJ)J when J² = −1.
    * `energy_identity`, `J_comm_diff`: with [B, J] = 0, 2ψJB = ψ² − (h(H + c) − h(H)), so J
      commutes with h(H + c) − h(H).
    * `delta_additive`, `affine_of_shift`: u(x + y) = u(x) + Δ(y) makes Δ additive with Δ(0) = 0
      and u(x) = u(0) + Δ(x).
    * `same_sign`, `shift_linear`: two levels sharing the shift ψ, with fixed signs, have the same
      sign, and then ψ(c) = sc with X = s(λ − E₀).
  Not covered here (left to the paper): differentiation, the use of (T₀), the step from "J commutes
  with every h(tH + c) − h(tH)" to [J, H] = 0 (periods of h and eigenprojections), continuity
  (which fixes the signs and makes the additive Δ linear), and the ordering that gives tr(K²) < 0.
-/
namespace AveragesPE
open Lean.Grind

variable {R : Type} [Ring R]

/-- The second derivative of e^{tX}e^{−tY} at t = 0 is X² − 2XY + Y²; with X = Y + D it is
    [Y, D] + D². -/
theorem second_deriv (Y D : R) :
    (Y + D) * (Y + D) - ((Y + D) * Y + (Y + D) * Y) + Y * Y = (Y * D - D * Y) + D * D := by
  rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
  grind

/-- For any additive map tr with tr(xy) = tr(yx), tr([A, K]K) = 0. -/
theorem trace_comm_zero {T : Type} [Ring T] (tr : R → T)
    (hsub : ∀ x y : R, tr (x - y) = tr x - tr y) (hcyc : ∀ x y : R, tr (x * y) = tr (y * x))
    (A K : R) : tr ((A * K - K * A) * K) = 0 := by
  rw [Ring.sub_eq_add_neg, Semiring.right_distrib, Ring.neg_mul, ← Ring.sub_eq_add_neg, hsub]
  have : tr (K * A * K) = tr (A * K * K) := by
    rw [Semiring.mul_assoc, hcyc, Semiring.mul_assoc]
  rw [this]; grind

/-- If [A, K] = θK with tr(θK²) = θ'·tr(K²) and tr(K²) cancellable, then θ' = 0. -/
theorem theta_zero {T : Type} [Ring T] (tr : R → T)
    (hsub : ∀ x y : R, tr (x - y) = tr x - tr y) (hcyc : ∀ x y : R, tr (x * y) = tr (y * x))
    (A K θ : R) (θ' : T) (h : A * K - K * A = θ * K)
    (hθ : tr (θ * K * K) = θ' * tr (K * K))
    (hcancel : ∀ y : T, y * tr (K * K) = 0 → y = 0) : θ' = 0 := by
  have z := trace_comm_zero tr hsub hcyc A K
  rw [h, hθ] at z
  exact hcancel _ z

/-- (A3′) on the scalar line: with A(x·1) = Z + φK, [Z, K] = 0, φ central, A(0)² = Z² = −h₀ and
    A(x·1)² = −hₓ, we get 2φ(ZK) + φ²K² = h₀ − hₓ. -/
theorem scalar_line_expand (Z K φ h0 hx : R) (hφ : ∀ x : R, φ * x = x * φ)
    (hZK : Z * K = K * Z) (h0' : Z * Z = -h0) (hx' : (Z + φ * K) * (Z + φ * K) = -hx) :
    (φ + φ) * (Z * K) + φ * φ * (K * K) = h0 - hx := by
  have e : (Z + φ * K) * (Z + φ * K) = Z * Z + (φ + φ) * (Z * K) + φ * φ * (K * K) := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
    have a1 : Z * (φ * K) = φ * (Z * K) := by
      rw [← Semiring.mul_assoc, ← hφ Z, Semiring.mul_assoc]
    have a2 : φ * K * Z = φ * (Z * K) := by rw [Semiring.mul_assoc, hZK]
    have a3 : φ * K * (φ * K) = φ * φ * (K * K) := by
      rw [Semiring.mul_assoc, ← Semiring.mul_assoc K φ K, ← hφ K, Semiring.mul_assoc,
        ← Semiring.mul_assoc]
    rw [a1, a2, a3, Semiring.right_distrib]; grind
  rw [e, h0'] at hx'
  grind

/-- Two values p, q of φ: from 2p(ZK) + p²K² = a and 2q(ZK) + q²K² = b (p, q central),
    pq(p − q)·K² = qa − pb and 2pq(q − p)·ZK = q²a − p²b. So when pq(p − q) is an invertible
    scalar and a, b are scalars, K² and ZK are scalars. -/
theorem scalar_line_solve (P K p q a b : R) (hp : ∀ x : R, p * x = x * p)
    (ha : (p + p) * P + p * p * K = a) (hb : (q + q) * P + q * q * K = b) :
    p * q * (p - q) * K = q * a - p * b ∧
      (p * q * (q - p) + p * q * (q - p)) * P = q * q * a - p * p * b := by
  have hp2 : ∀ x : R, p * p * x = x * (p * p) := by
    intro x; rw [Semiring.mul_assoc, hp x, ← Semiring.mul_assoc, hp x, Semiring.mul_assoc]
  constructor
  · have f1 : q * ((p + p) * P) = p * q * P + p * q * P := by
      rw [← Semiring.mul_assoc, Semiring.left_distrib, ← hp q, Semiring.right_distrib]
    have f2 : q * (p * p * K) = p * q * p * K := by
      rw [← Semiring.mul_assoc, ← Semiring.mul_assoc q p p, ← hp q]
    have f3 : p * ((q + q) * P) = p * q * P + p * q * P := by
      rw [← Semiring.mul_assoc, Semiring.left_distrib, Semiring.right_distrib]
    have f4 : p * (q * q * K) = p * q * q * K := by
      rw [← Semiring.mul_assoc, ← Semiring.mul_assoc]
    have f5 : p * q * (p - q) * K = p * q * p * K - p * q * q * K := by
      rw [Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, Semiring.right_distrib,
        Ring.neg_mul, ← Ring.sub_eq_add_neg]
    rw [← ha, ← hb, Semiring.left_distrib q, Semiring.left_distrib p, f1, f2, f3, f4, f5]
    grind
  · have g1 : q * q * ((p + p) * P) = p * q * q * P + p * q * q * P := by
      rw [← Semiring.mul_assoc (q * q), Semiring.left_distrib, ← hp (q * q), ← Semiring.mul_assoc,
        Semiring.right_distrib]
    have g2 : q * q * (p * p * K) = p * p * (q * q) * K := by
      rw [← Semiring.mul_assoc, ← hp2 (q * q)]
    have g3 : p * p * (q * q * K) = p * p * (q * q) * K := by
      rw [← Semiring.mul_assoc]
    have g4 : p * p * ((q + q) * P) = p * q * p * P + p * q * p * P := by
      have : p * p * q = p * q * p := by rw [Semiring.mul_assoc, Semiring.mul_assoc, hp q]
      rw [← Semiring.mul_assoc, Semiring.left_distrib, this, Semiring.right_distrib]
    have g5 : (p * q * (q - p) + p * q * (q - p)) * P
        = (p * q * q * P - p * q * p * P) + (p * q * q * P - p * q * p * P) := by
      have : p * q * (q - p) * P = p * q * q * P - p * q * p * P := by
        rw [Ring.sub_eq_add_neg, Semiring.left_distrib, Ring.mul_neg, Semiring.right_distrib,
          Ring.neg_mul, ← Ring.sub_eq_add_neg]
      rw [Semiring.right_distrib, this]
    rw [← ha, ← hb, Semiring.left_distrib (q * q), Semiring.left_distrib (p * p), g1, g2, g3, g4,
      g5]
    grind

/-- Z = −(ZJ)J when J² = −1: if ZJ is a scalar, Z is a multiple of J. -/
theorem Z_of_ZJ (Z J : R) (hJ : J * J = -1) : Z = -(Z * J * J) := by
  rw [Semiring.mul_assoc, hJ, Ring.mul_neg, Semiring.mul_one]; grind

/-- Energy conservation, algebraic part: with [B, J] = 0, J² = −1, ψ central,
    B² = −h(H) and (B + ψJ)² = −h(H + c): 2ψ·JB = ψ² − (h(H + c) − h(H)). -/
theorem energy_identity (B J ψ hH hc : R) (hψ : ∀ x : R, ψ * x = x * ψ)
    (hBJ : B * J = J * B) (hJ : J * J = -1) (hB : B * B = -hH)
    (hs : (B + ψ * J) * (B + ψ * J) = -hc) :
    (ψ + ψ) * (J * B) = ψ * ψ - (hc - hH) := by
  have e : (B + ψ * J) * (B + ψ * J) = B * B + (ψ + ψ) * (J * B) + ψ * ψ * (J * J) := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
    have a1 : B * (ψ * J) = ψ * (J * B) := by
      rw [← Semiring.mul_assoc, ← hψ B, Semiring.mul_assoc, hBJ]
    have a2 : ψ * J * B = ψ * (J * B) := Semiring.mul_assoc _ _ _
    have a3 : ψ * J * (ψ * J) = ψ * ψ * (J * J) := by
      rw [Semiring.mul_assoc, ← Semiring.mul_assoc J ψ J, ← hψ J, Semiring.mul_assoc,
        ← Semiring.mul_assoc]
    rw [a1, a2, a3, Semiring.right_distrib]; grind
  rw [e, hB, hJ, Ring.mul_neg, Semiring.mul_one] at hs
  grind

/-- So J commutes with h(H + c) − h(H): it commutes with ψ² and with JB. -/
theorem J_comm_diff (B J ψ hH hc : R) (hψ : ∀ x : R, ψ * x = x * ψ)
    (hBJ : B * J = J * B) (hJ : J * J = -1) (hB : B * B = -hH)
    (hs : (B + ψ * J) * (B + ψ * J) = -hc) :
    J * (hc - hH) = (hc - hH) * J := by
  have e := energy_identity B J ψ hH hc hψ hBJ hJ hB hs
  have d : hc - hH = ψ * ψ - (ψ + ψ) * (J * B) := by rw [e]; grind
  have c1 : J * (ψ * ψ) = ψ * ψ * J := by
    rw [← Semiring.mul_assoc, ← hψ J, Semiring.mul_assoc, ← hψ J, ← Semiring.mul_assoc]
  have hc2 : ∀ x : R, (ψ + ψ) * x = x * (ψ + ψ) := by
    intro x; rw [Semiring.right_distrib, Semiring.left_distrib, hψ x]
  have c2 : J * ((ψ + ψ) * (J * B)) = (ψ + ψ) * (J * B) * J := by
    rw [← Semiring.mul_assoc J, ← hc2 J, Semiring.mul_assoc, Semiring.mul_assoc (ψ + ψ) (J * B) J,
      Semiring.mul_assoc J B J, hBJ]
  rw [d, Ring.sub_eq_add_neg, Semiring.left_distrib J, Semiring.right_distrib _ _ J, Ring.mul_neg,
    Ring.neg_mul, c1, c2]

/-- One plane: u(x + y) = u(x) + Δ(y) for all x, y makes Δ additive. -/
theorem delta_additive {S : Type} [Ring S] (u Δ : S → S) (h : ∀ x y : S, u (x + y) = u x + Δ y)
    (y y' : S) : Δ (y + y') = Δ y + Δ y' := by
  have a := h 0 (y + y')
  have b := h 0 y
  have c := h y y'
  have e : (0 : S) + (y + y') = 0 + y + y' := by grind
  have e' : (0 : S) + y = y := by grind
  rw [e, h (0 + y) y', e'] at a
  rw [e'] at b
  rw [b] at a
  grind

/-- …with Δ(0) = 0 and u(x) = u(0) + Δ(x): u is affine once Δ is linear (continuity, paper). -/
theorem affine_of_shift {S : Type} [Ring S] (u Δ : S → S) (h : ∀ x y : S, u (x + y) = u x + Δ y)
    (x : S) : Δ 0 = 0 ∧ u x = u 0 + Δ x := by
  constructor
  · have a := h 0 0
    have : (0 : S) + 0 = 0 := by grind
    rw [this] at a; grind
  · have a := h 0 x
    have : (0 : S) + x = x := by grind
    rw [this] at a; exact a

/-- Every Hamiltonian: two levels λ ≠ λ′ share the shift ψ. If X + ψ(c) = s(λ + c − E₀) and
    X′ + ψ(c) = s′(λ′ + c − E₀) at c = 0 and c = 1 with fixed signs s, s′ (continuity fixes them,
    paper), then s = s′. -/
theorem same_sign {S : Type} [CommRing S] (X X' ψ0 ψ1 s s' l l' E0 : S)
    (a0 : X + ψ0 = s * (l + 0 - E0)) (a1 : X + ψ1 = s * (l + 1 - E0))
    (b0 : X' + ψ0 = s' * (l' + 0 - E0)) (b1 : X' + ψ1 = s' * (l' + 1 - E0)) : s = s' := by
  grind

/-- …and then, with ψ(0) = 0, X = s(λ − E₀) and ψ(c) = sc. -/
theorem shift_linear {S : Type} [CommRing S] (X ψ0 ψc c s l E0 : S) (hψ0 : ψ0 = 0)
    (a0 : X + ψ0 = s * (l + 0 - E0)) (ac : X + ψc = s * (l + c - E0)) :
    X = s * (l - E0) ∧ ψc = s * c := by
  constructor <;> grind

end AveragesPE
