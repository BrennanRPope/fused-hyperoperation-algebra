/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, fourth look: what has to stay positive along the rewind,
  and the barrier chain as composed maps (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib). Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/RiemannRound4.lean

  First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound4
  (namespace CircleOp.RiemannRound4, the library's convention; no proof changed).

  Jets: h0 .. h4 stand for H, H', H'', H''', H'''' at one point; the de Bruijn–Newman flow
  dH/dt = −H'' sends h_k to −h_(k+2). Divisions are cleared (multiplied through by a power
  of h2), so every statement is a ring identity. Nothing here touches the zeros of zeta.

  1. A lone pair: the local fused discriminant −8H/H'' at the extremum is s² − 4p, and
     Laguerre's expression H'² − H·H'' is at least half of it.
  2. The Laguerre ladder: the flow moves L = H'² − H·H'' by −L'' + 2·L(H').
  3. The exact rate of the local fused discriminant at an extremum (h1 = 0):
     dD/dt = 8 + D·(log H'')'', and sin is a steady state.
  4. Letters of the barrier chain: z ↦ u·z + c; the order gap is c·(u − v); a letter with
     c = 1 has the Euler factor as its centre.
  5. The link dial: 1 + k·u + d·u² = (1 − α u)(1 − β u) with α + β = −k, α·β = d.
  6. The hypernatural log at height s: −(ln_s g)″ = (s·f′² − f·f″)/f² (g = f/f(x0)); at an extremum
     the height drops out, and s = 1 − 1/n is Laguerre's degree-n inequality.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound4

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. A lone pair -/

/-- For H = (x − a)² − c (roots a ± √c, fused point s = 2a, p = a² − c), the local fused
    discriminant −8·H(a)/H''(a) = −8·(−c)/2 equals s² − 4p = 4c. Stated with H'' = 2 cleared. -/
theorem pair_local_disc (a c : R) :
    (-8) * (-c) = 2 * ((2 * a) * (2 * a) - 4 * (a * a - c)) := by grind

/-- Laguerre's expression of the lone pair: H'² − H·H'' = 2(x − a)² + 2c, so its least value
    is 2c, half the fused discriminant 4c. -/
theorem pair_laguerre (x a c : R) :
    (2 * (x - a)) * (2 * (x - a)) - ((x - a) * (x - a) - c) * 2
      = 2 * ((x - a) * (x - a)) + 2 * c := by grind

/-- The flow on a lone pair: H = (x − a)² − c − 2t solves the flow (H_t = −2 = −H''), and the
    local fused discriminant −8H(a)/H''(a) is 4c + 8t: a straight line through the merger. -/
theorem pair_flow_disc (c t : R) : (-8) * (-c - 2 * t) = 2 * (4 * c + 8 * t) := by grind

/-! ## 2. The Laguerre ladder -/

/-- With L = h1² − h0·h2, the flow gives L_t = −2h1h3 + h2² + h0h4 (from h_k ↦ −h_(k+2)),
    L'' = h2² − h0h4, and L(H') = h2² − h1h3. The ladder identity: L_t = −L'' + 2·L(H'). -/
theorem laguerre_ladder (h0 h1 h2 h3 h4 : R) :
    (2 * h1 * (-h3) - (-h2) * h2 - h0 * (-h4))
      = -(h2 * h2 - h0 * h4) + 2 * (h2 * h2 - h1 * h3) := by grind

/-- For the lone pair the ladder reads 4 − 8 = −4 under the rewind: the diffusion term L'' = 4
    is beaten 2 : 1 by the sink 2·L(H') = 8 (h0 = u² − c, h1 = 2u, h2 = 2, h3 = h4 = 0). -/
theorem pair_ladder_rates (u c : R) :
    (2 * 2 - (u * u - c) * 0 = 4) ∧ (2 * 2 - (2 * u) * 0 = 4) := by
  constructor <;> grind

/-! ## 3. The exact rate of the local fused discriminant -/

/-- At an extremum (h1 = 0) with D = −8h0/h2 and κ = h4/h2 − (h3/h2)² = (log H'')'':
    the time derivative of D along the moving extremum, times h2³, is
    −8(−h2² + h0·h4)·h2 + 8·h0·h3², and 8 + D·κ times h2³ is 8h2³ − 8h0(h4h2 − h3²).
    They agree: dD/dt = 8 + D·κ. -/
theorem local_disc_rate (h0 h2 h3 h4 : R) :
    (-8) * (-(h2 * h2) + h0 * h4) * h2 + 8 * h0 * (h3 * h3)
      = 8 * (h2 * h2 * h2) - 8 * h0 * (h4 * h2 - h3 * h3) := by grind

/-- sin is a steady state: its jets at a crest satisfy h2 = −h0, h3 = 0, h4 = h0, and then
    8·h2³ − 8·h0·(h4·h2 − h3²) = 0, so dD/dt = 0 (the zeros of sin do not move). -/
theorem sin_steady (h0 : R) :
    8 * ((-h0) * (-h0) * (-h0)) - 8 * h0 * (h0 * (-h0) - 0 * 0) = 0 := by grind

/-! ## 4. Letters of the barrier chain -/

/-- Two letters z ↦ u z + c and z ↦ v z + c: the two orders differ by c·(u − v). With c = 0
    (shared centre 0, the Euler product) they commute. -/
theorem letter_order_gap (u v c z : R) :
    (u * (v * z + c) + c) - (v * (u * z + c) + c) = c * (u - v) := by grind

/-- A letter with c = 1 fixes z exactly when (1 − u)·z = 1: its centre is the Euler factor
    1/(1 − u), u = p^(−s). -/
theorem letter_centre (u z : R) (h : (1 - u) * z = 1) : u * z + 1 = z := by grind

/-! ## 5. The link dial -/

/-- A degree-2 link 1 + k u + d u² splits as (1 − α u)(1 − β u) when α + β = −k and αβ = d. -/
theorem link_split (α β u : R) :
    (1 - α * u) * (1 - β * u) = 1 + (-(α + β)) * u + (α * β) * (u * u) := by grind

/-! ## 6. The hypernatural log as a Laguerre dial -/

/-- Laguerre at height s, at an extremum (f′ = 0): s·0² − f·f″ = −f·f″ for every s, so the local
    fused discriminant does not depend on the height of the hypernatural log. -/
theorem height_free_at_extremum (s f0 f2 : R) : s * (0 * 0) - f0 * f2 = -(f0 * f2) := by grind

/-- Height s = 1 − 1/n (stated as s·n = n − 1): n times the height-s Laguerre expression is
    Laguerre's degree-n expression (n − 1)·f′² − n·f·f″. -/
theorem height_is_degree (s n f0 f1 f2 : R) (h : s * n = n - 1) :
    n * (s * (f1 * f1) - f0 * f2) = (n - 1) * (f1 * f1) - n * (f0 * f2) := by
  have : n * (s * (f1 * f1)) = (s * n) * (f1 * f1) := by grind
  grind

end Ring
end CircleOp.RiemannRound4
