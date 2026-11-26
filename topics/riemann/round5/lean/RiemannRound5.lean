/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, fifth look: the untried dial pairs (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib). Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/RiemannRound5.lean

  First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound5
  (namespace CircleOp.RiemannRound5, the library's convention; no proof changed).

  Divisions are cleared, so every statement is a ring identity. Nothing here touches the zeros of zeta.

  1. Height s with homothety k: the homothety log l_k(X) = (1 − X^(−k))/k turns × into the circle
     operation at −k and undoes each prime factor of the prime cube exactly; the hypernatural log
     at height s = 1 + k has the same bend, and the third-order gap is k(k + 1)/12, so the two bends
     meet only at k = 0 and k = −1.
  2. Ratio r with order q: q-multinomials as (q;q)-ratios (two small cases), and the necklace step
     behind the order-dial identity.
  3. Light: a boost keeps x·p; turning the zero speed by i keeps its size.
  4. Height s across the roundness kernels: a non-real pair makes Laguerre's expression negative,
     a real pair keeps it positive; two real zeros give h* = 1 − 1/N_eff.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound5

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. Height s with homothety k -/

/-- With a = X^(−k), b = Y^(−k), k·l_k(XY) = 1 − ab and k·l_k(X) = 1 − a. The product law
    k·l(XY) = k·l(X) + k·l(Y) − (k·l(X))(k·l(Y)) is the circle operation at −k, cleared of k. -/
theorem lk_product (a b : R) :
    1 - a * b = (1 - a) + (1 - b) - (1 - a) * (1 - b) := by grind

/-- One prime factor of the prime cube: X = (1 − k u)^(−1/k), so X^(−k) = 1 − k u and
    k·l_k(X) = 1 − (1 − k u) = k u. The homothety log returns u = p^(−s) exactly. -/
theorem lk_factor (k u : R) : 1 - (1 - k * u) = k * u := by grind

/-- Third-order coefficients: the hypernatural log at height s has s²/4 + s/12, the homothety log
    has (k² + 3k + 2)/6. At s = 1 + k (same bend) their gap, times 12, is k(k + 1). -/
theorem bends_meet_twice (k : R) :
    3 * ((1 + k) * (1 + k)) + (1 + k) - 2 * (k * k + 3 * k + 2) = k * (k + 1) := by grind

/-! ## 2. Ratio r with order q -/

/-- Two letters, once each: the two orderings weigh 1 + q, and the q-multinomial
    (q;q)₂/((q;q)₁(q;q)₁) is the same: (1 − q)(1 − q²) = ((1 − q)(1 − q))(1 + q). -/
theorem qmulti_11 (q : R) : (1 - q) * (1 - q * q) = ((1 - q) * (1 - q)) * (1 + q) := by grind

/-- Exponents (2, 1): the three orderings weigh 1 + q + q², and
    (q;q)₃ = (q;q)₂ · (q;q)₁ · (1 + q + q²). -/
theorem qmulti_21 (q : R) :
    (1 - q) * (1 - q * q) * (1 - q * q * q)
      = ((1 - q) * (1 - q * q)) * (1 - q) * (1 + q + q * q) := by grind

/-- Necklace step at N = 2: log of a homothety cube at dial c has p^(−2s)-coefficient c²/2, and
    zeta(s)^c · zeta(2s)^M(c,2) supplies c/2 + M with 2M = c² − c. -/
theorem necklace_two (c : R) : c + (c * c - c) = c * c := by grind

/-! ## 3. Light -/

/-- A boost multiplies the light-cone coordinates by λ and μ with λμ = 1, so x·p is kept:
    the Berry–Keating orbits xp = E are boost orbits. -/
theorem boost_keeps_xp (lam mu a b : R) (h : lam * mu = 1) : (lam * a) * (mu * b) = a * b := by
  grind

/-- Wick rotation turns the zero speed v into i·v. With i² = −1, (i v)(−i v) = v²:
    the same speed, at a right angle, so the zero leaves the line. -/
theorem light_same_speed (i v : R) (h : i * i = -1) : (i * v) * ((-i) * v) = v * v := by
  grind

/-! ## 4. Height s across the roundness kernels -/

/-- A non-real pair: f = (x − a)² + b². Laguerre's expression f′² − f·f″ = 2(x − a)² − 2b², negative
    at x = a when b ≠ 0, so the critical height there is above 1. (Compare round 4's real pair: + 2c.) -/
theorem pair_offline_detected (x a b : R) :
    (2 * (x - a)) * (2 * (x - a)) - ((x - a) * (x - a) + b * b) * 2
      = 2 * ((x - a) * (x - a)) - 2 * (b * b) := by grind

/-- Two real zeros, a_i = 1/(x − r_i): f′/f = a₁ + a₂, f″/f = 2a₁a₂, so
    (1 − h*)(a₁ + a₂)² = a₁² + a₂², i.e. h* = 1 − 1/N_eff with N_eff = (a₁ + a₂)²/(a₁² + a₂²). -/
theorem neff_two (a1 a2 : R) :
    (a1 + a2) * (a1 + a2) - 2 * a1 * a2 = a1 * a1 + a2 * a2 := by grind

end Ring
end CircleOp.RiemannRound5
