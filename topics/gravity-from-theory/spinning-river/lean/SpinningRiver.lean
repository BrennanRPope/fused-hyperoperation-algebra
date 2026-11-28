/-
DRAFT (Nov 26, 2026). First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.SpinningRiver.
The spinning river: the cheap algebra behind spinning-river/README.md.

Legs of a cell, in the order (time, radial, azimuth). The map from the flat lattice's legs to the
cell's legs adds v·(time leg) and −v·α·(azimuth leg) to the radial leg: two shears, the dial's
flat middle. Builds on FusedAlgebra.Topics.JetFields (M3, M3.mul).
-/
import FusedAlgebra.Topics.JetFields

namespace CircleOp.SpinningRiver

open CircleOp
open CircleOp.JetFields
open Lean.Grind

section
variable {R : Type} [CommRing R]

/-- The cell map minus the identity: the radial row picks up v (from time) and −v·α (from azimuth). -/
def shear (v α : R) : M3 R := ⟨0, 0, 0, v, 0, -(v * α), 0, 0, 0⟩

/-- The cell map itself, I + shear. -/
def cellMap (v α : R) : M3 R := ⟨1, 0, 0, v, 1, -(v * α), 0, 0, 1⟩

/-- **The cell map is a dual number:** its shear part squares to zero. -/
theorem shear_sq (v α : R) : (shear v α).mul (shear v α) = ⟨0, 0, 0, 0, 0, 0, 0, 0, 0⟩ := by
  simp only [shear, M3.mul, M3.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- **Two cell maps compose by adding their shears** (boosts and turns at the flat middle add):
    (I + S₁)(I + S₂) = I + S₁ + S₂. -/
theorem cellMap_mul (v α w β : R) :
    (cellMap v α).mul (cellMap w β) = ⟨1, 0, 0, v + w, 1, -(v * α) - w * β, 0, 0, 1⟩ := by
  simp only [cellMap, M3.mul, M3.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- The cell map keeps volume: no layer but the 1 survives (P = Q = R = 0 for the shear). -/
theorem shear_layers (v α : R) :
    (shear v α).tr = 0 ∧ (shear v α).e2 = 0 ∧ (shear v α).det = 0 := by
  simp only [shear, M3.tr, M3.e2, M3.det, det9]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- **Kerr's circular-orbit period, cleared of roots.** With r = s², M = m²:
    m·(s⁶ − m²a²) = (m s³ − m² a)(s³ + m a), so T = 2π (s⁶ − m²a²)/(m s³ − m² a) = 2π (s³/m + a). -/
theorem kerr_period (s m a : R) :
    m * (s^6 - m^2 * a^2) = (m * s^3 - m^2 * a) * (s^3 + m * a) := by grind

/-! ### Round 7: deriving the turn (checks D1–D6 of spinning-river/checks.py) -/

/-- **Area law (D1a).** A shell turning at rate Ω has turn f = Ω (r²+b²) sin²θ; its clock gap 2π f is
    2 Ω times the area π (r²+b²) sin²θ of the latitude circle (π kept as a symbol). -/
theorem area_law (pi Ω r b s2 : R) :
    2 * pi * (Ω * (r^2 + b^2) * s2) = 2 * Ω * (pi * (r^2 + b^2) * s2) := by grind

/-- **Rim at light speed (D4).** The innermost shell turns at Ω₀ with Ω₀ b² = κ; if its rim speed Ω₀ b
    is 1, then b = κ. -/
theorem rim_at_light (κ b Ω0 : R) (hΩ : Ω0 * b^2 = κ) (hrim : Ω0 * b = 1) : b = κ := by grind

/-- **The t–φ block (D6a).** With V = v², V ρ² = 2 M r, ρ² = r² + a² c², c² + s² = 1:
    g_tφ² − g_tt g_φφ = s² (r² − 2 M r + a²), where g_tt = V − 1, g_tφ = −V a s², g_φφ = (r²+a²) s² + V a² s⁴. -/
theorem tphi_block (V ρ2 M r a c2 s2 : R)
    (hv : V * ρ2 = 2 * M * r) (hρ : ρ2 = r^2 + a^2 * c2) (hc : c2 + s2 = 1) :
    (V * a * s2)^2 - (V - 1) * ((r^2 + a^2) * s2 + V * a^2 * s2^2) = s2 * (r^2 - 2 * M * r + a^2) := by
  grind

/-- **The horizon turns at the shell rate (D6b).** Where r² − 2 M r + a² = 0, the dragging rate
    −g_tφ / g_φφ equals a / (r² + a²), written without division:
    (−g_tφ) (r² + a²) = a g_φφ. -/
theorem horizon_rate (V ρ2 M r a c2 s2 : R)
    (hv : V * ρ2 = 2 * M * r) (hρ : ρ2 = r^2 + a^2 * c2) (hc : c2 + s2 = 1)
    (hor : r^2 - 2 * M * r + a^2 = 0) :
    (V * a * s2) * (r^2 + a^2) = a * ((r^2 + a^2) * s2 + V * a^2 * s2^2) := by
  grind

/-! ### Round 8: deriving the speed law (checks S1–S3) -/

/-- **The complex distance (S1a).** With an element i, i·i = −1, and a lattice point
    x² + y² = (r² + a²) s², z = r c, c² + s² = 1: the squared distance to the source moved up the
    axis by i·a is (r − i a c)². -/
theorem complex_distance (i r a c s2 : R) (hi : i * i = -1) (hc : c^2 + s2 = 1) :
    (r^2 + a^2) * s2 + (r * c - i * a)^2 = (r - i * a * c)^2 := by grind

/-- **Its real part gives the speed law (S2a).** (r − i a c)(r + i a c) = r² + a² c² = ρ², so
    M/(r − i a c) = M (r + i a c)/ρ², whose real part is M r/ρ². -/
theorem conj_product (i r a c : R) (hi : i * i = -1) :
    (r - i * a * c) * (r + i * a * c) = r^2 + a^2 * c^2 := by grind

/-! ### Round 9: the shifted source's rays (check T1) -/

/-- **The ray has unit length (T1).** For real p, q with k = p + q × p (before dividing by |p|²):
    |p + q × p|² = (p·p)(1 + q·q) − (p·q)². With p·p − q·q = 1 and p·q = 0 (E = p + iq a complex unit
    vector) this is (p·p)², so k/|p|² is a unit vector. -/
theorem ray_norm (p0 p1 p2 q0 q1 q2 : R) :
    (p0 + (q1 * p2 - q2 * p1))^2 + (p1 + (q2 * p0 - q0 * p2))^2 + (p2 + (q0 * p1 - q1 * p0))^2
      = (p0^2 + p1^2 + p2^2) * (1 + q0^2 + q1^2 + q2^2) - (p0 * q0 + p1 * q1 + p2 * q2)^2 := by grind

/-! ### Round 11: two bodies at rest (checks B3, B4) -/

/-- **Two bodies combine by velocity addition (B3).** Cross-multiplied form of
    (1 − u)/(1 + u) = ((1 − u₁)/(1 + u₁))·((1 − u₂)/(1 + u₂)) for u = (u₁ + u₂)/(1 + u₁u₂):
    on the axis each body's clock factor (1 − u)/(1 + u) multiplies, and the u's add like speeds. -/
theorem velocity_addition (u1 u2 : R) :
    ((1 + u1 * u2) - (u1 + u2)) * ((1 + u1) * (1 + u2))
      = ((1 + u1 * u2) + (u1 + u2)) * ((1 - u1) * (1 - u2)) := by grind

/-- **The strut's push (B4).** With e^γ = (ℓ² − (M₁+M₂)²)/(ℓ² − (M₁−M₂)²) on the strut,
    (e^−γ − 1)/4 = M₁M₂/(ℓ² − (M₁+M₂)²): the numerators differ by exactly 4 M₁M₂. -/
theorem strut_push (l M1 M2 : R) :
    (l^2 - (M1 - M2)^2) - (l^2 - (M1 + M2)^2) = 4 * (M1 * M2) := by grind


/-! ### Round 12: two moving bodies (check O1) -/

/-- **Clock rates multiply (O1, rule B).** Each body alone slows clocks by 1 − 2U + 2U² at slow-orbit order.
    Multiplying the two bodies' factors gives the one-body form with U = U₁ + U₂, cross term 4U₁U₂ included,
    plus only terms of higher order (U₁U₂ times another U). -/
theorem clock_multiply (a b : R) :
    (1 - 2 * a + 2 * a^2) * (1 - 2 * b + 2 * b^2)
      = 1 - 2 * (a + b) + 2 * (a + b)^2 - 4 * (a * b) * (a + b) + 4 * (a^2 * b^2) := by grind

/-- **Round 13: the 2PN miss, two ways of writing it** (times 6, to stay in a ring). With total
mass 1, body 1 of mass a and body 2 of mass 1 − a, the three rules' redshift coefficient of x³ is
off by a(1 − a)(23a² − 50a + 36)/3. In the usual variables ν = a(1 − a) and Δ = (1 − a) − a that is
ν(15/2 + 9/2·Δ − 23ν/3). -/
theorem two_pn_miss (a : R) :
    2 * a * (1 - a) * (23 * a^2 - 50 * a + 36)
      = a * (1 - a) * (45 + 27 * ((1 - a) - a) - 46 * (a * (1 - a))) := by grind

/-- **Round 14: adding gravity's own field energy (check O10)**, times 3 to stay in a ring. With body 1 of
mass a and body 2 of mass 1 − a, the closed-form field-energy pieces leave the miss a(1 − a)(5 − 10a + 3a²).
The four static field-on-field numbers α = −1/2, β = 1, γ = 0, δ = −2 change body 1's x³ coefficient by the
Newtonian-order map in the second bracket. What is left is a(1 − a)²(3 − a)/3. -/
theorem fourth_rule_residual (a : R) :
    3 * (a * (1 - a) * (5 - 10 * a + 3 * a^2))
      + (8 * (1 - a) * (a * (1 - a)) + 16 * (1 - a)^2 * (a * (1 - a)) - 24 * a * (1 - a)^3
         + 6 * a^2 * (1 - a) - 12 * a * (1 - a)^2)
      = a * (1 - a)^2 * (3 - a) := by grind

/-- **The leftover is one point-like term (check O10).** Raising β by 1/4 and lowering δ by 1/4 (an extra
X̂ ⊃ m_a m_b²/(4 r₁₂² r_b) at each body) cancels the residual exactly, for every mass ratio:
the map at Δβ = 1/4, Δδ = −1/4, times 3, is minus the residual above. -/
theorem fourth_rule_point_term (a : R) :
    ((1 - a) * (a * (1 - a)) - (1 - a)^2 * (a * (1 - a))
      - 3 * a * (1 - a)^2) = -(a * (1 - a)^2 * (3 - a)) := by grind

/-- **Round 15: the leftover derived (check O11)**, everything times 4 to stay in a ring. A body's potential
is M/r plus, from the delay, −(M/2) s with s = a·n its own acceleration along the direction to the field point;
write u = 1/r, so 2 × V = 2Mu − Ms. The cross term of g₀₀'s −2V² is then −4 times the acceleration piece already
in X̂, Xacc = −(M²/2) s u (4 Xacc = −2M² s u). So g₀₀'s self part is 8 Xacc − 4 Xacc, which is the
"quarter point term" of round 14. -/
theorem own_field_squared (M s u : R) :
    -2 * ((2 * M * u - M * s)^2 - (2 * M * u)^2 - (M * s)^2) = -4 * (-2 * M^2 * s * u) := by grind

/-- **Round 16: Kerr's orbit factor, split (checks KS2–KS4).** In s = √(M/r), Kerr's
T = 1 − 6s² + 8χs³ − 3χ²s⁴ is the sum of its pieces: −2 from P·P and −4 from 2E·P (−4 = 2² × the single
copy's −1), +8 from 2E·S (= 2 × the single copy's 4), and −3 from E·E (Newton's pull of the shifted source). -/
theorem kerr_factor_pieces (s χ : R) :
    1 - 6 * s^2 + 8 * χ * s^3 - 3 * χ^2 * s^4
      = 1 + ((-2) + (-1) * 2^2) * s^2 + 2 * 4 * χ * s^3 + (-3) * χ^2 * s^4 := by grind

/-- **Rule A's push (check KS4).** Boosting the rest field h₀₀ = 2U, hₓₓ = 2U to speed v gives
h₀ₓ = −(h₀₀ + hₓₓ) v: half from the time part, half from the space part. -/
theorem boost_push (U v : R) : -((2 * U) + (2 * U)) * v = -4 * U * v := by grind

end

/-! ### Rounds 16 to 31 over a field (the topic file's `section Field`; proved by grind in place of
linear_combination, mul_eq_zero and rw with sub_div) -/
section Field
variable {K : Type} [Field K]

/-- **The ring with Kerr's mass, spin and shape (check KS5).** A ring of radius R with R² = 2a² (Kerr's
quadrupole) and spin per unit mass R·v = a turns at v² = 1/2, i.e. c/√2. -/
theorem ring_speed (a Rr v : K) (ha : a ≠ 0) (hR : Rr^2 = 2 * a^2) (hv : Rr * v = a) :
    2 * v^2 = 1 := by
  have h : a^2 * (2 * v^2 - 1) = 0 := by grind
  have ha2 : a^2 ≠ 0 := by intro h0; apply ha; grind
  grind

/-- **Round 17: the tower of shapes sums to one term (check TW1).** Newton's orbit factor around the shifted
source is (R² − 4a²)/(R² − a²) in flat radius R. In the shift's own radius r, r² = R² − a², it is exactly
1 − 3a²/r²: Kerr's whole tower of shapes leaves a single −3. -/
theorem shifted_orbit_factor (R a r : K) (hr : r^2 = R^2 - a^2) (h : r ≠ 0) :
    (R^2 - 4 * a^2) / (R^2 - a^2) = 1 - 3 * a^2 / r^2 := by
  have h1 : R^2 - a^2 = r^2 := hr.symm
  have hr2 : r^2 ≠ 0 := by intro h0; apply h; grind
  rw [h1]
  grind

/-- **Round 18: Laplace's rule around the shifted source (check TL1).** In the equator Newton's radial factor
(R² − 4a²)/(R² − a²) and up-and-down factor (R² + 2a²)/(R² − a²) add to exactly 2: empty space off the ring. -/
theorem shifted_laplace (R a : K) (h : R^2 - a^2 ≠ 0) :
    (R^2 - 4 * a^2) / (R^2 - a^2) + (R^2 + 2 * a^2) / (R^2 - a^2) = 2 := by
  grind

/-- **Round 18: Kerr's radial plus vertical factor (check TL5).** With s² = M/r, Kerr's equatorial radial factor
1 − 6s² + 8χs³ − 3χ²s⁴ and up-and-down factor 1 − 4χs³ + 3χ²s⁴ add to 2 − 6s² + 4χs³: the shape term cancels. -/
theorem kerr_laplace_sum (s χ : K) :
    (1 - 6 * s^2 + 8 * χ * s^3 - 3 * χ^2 * s^4) + (1 - 4 * χ * s^3 + 3 * χ^2 * s^4)
      = 2 - 6 * s^2 + 4 * χ * s^3 := by
  grind

/-- **Round 19: Kerr's polar orbits at Newtonian strength (check PL1).** With M = 0, the radial-to-polar prefactor of
Kerr's polar orbit, (r² + a²)³/(r²(r² + a²)²), is the shifted source's 1 + a²/r². -/
theorem kerr_polar_prefactor_newton (r a : K) (hr : r ≠ 0) (h : r^2 + a^2 ≠ 0) :
    (r^2 + a^2)^3 / (r^2 * (r^2 + a^2)^2) = 1 + a^2 / r^2 := by
  have hr2 : r^2 ≠ 0 := by intro h0; apply hr; grind
  grind

/-- **Round 19: the same prefactor without spin (check PL1)** is Schwarzschild's 1 − 6M/r. -/
theorem kerr_polar_prefactor_schwarzschild (r M : K) (hr : r ≠ 0) :
    (r^6 - 6 * M * r^5) / (r^2 * r^4) = 1 - 6 * M / r := by
  grind

/-- **Round 20: the river's time part alone, over the poles (check PX2).** With only the E·E piece, a polar
orbit has radial rate² `M(r²+a²)/(r+M)` and Carter constant `Mr²/(r+M)`; their ratio is the shifted
source's prefactor `1 + a²/r²`, at every strength `M`. -/
theorem ee_alone_polar_prefactor (r a M : K) (hr : r ≠ 0) (hM : M ≠ 0) (hrM : M + r ≠ 0) :
    (M * (a^2 + r^2) / (M + r)) / (M * r^2 / (M + r)) = 1 + a^2 / r^2 := by
  have hr2 : r^2 ≠ 0 := by intro h0; apply hr; grind
  grind

/-- **Round 20: and its elliptic modulus (check PX2)** `m = (1 − E²) a²/Kc` with `E² = r/(r+M)` is the
shifted source's `a²/r²`. -/
theorem ee_alone_polar_modulus (r a M : K) (hr : r ≠ 0) (hM : M ≠ 0) (hrM : M + r ≠ 0) :
    (1 - r / (M + r)) * a^2 / (M * r^2 / (M + r)) = a^2 / r^2 := by
  have hr2 : r^2 ≠ 0 := by intro h0; apply hr; grind
  grind

/-- **Round 20: the equator's cancellation (check PX4).** The shape-meeting-mass term the push brings in
(E·P +20, P·P +8) is taken out exactly by the twist (S·S −16, E·S −12). -/
theorem equator_push_twist_cancel : (20 : Int) + 8 + (-16) + (-12) = 0 := by decide

/-- **Round 20: over the poles the twist is absent (check PX3)**, so the push's E·P +16 and P·P −1 stand:
Kerr's +15. -/
theorem polar_push_total : (0 : Int) + 16 + (-1) = 15 := by decide

/-- **Round 21: what the jet rule adds (check JB1).** Push and twist as one jet give the field `(E + P + S)²`;
added as two velocities they give `(E + P)² + (E + S)² − E²`. The difference is the cross piece `2PS`. -/
theorem jet_minus_two_velocities (E P S : K) :
    (E + P + S)^2 - ((E + P)^2 + (E + S)^2 - E^2) = 2 * P * S := by
  grind

/-- **Round 21: the push on an orbit answers the turning shell (check JB4).** With one jet the push on a circular
or spherical orbit is `2Mr(E + S)/Δ`, so `k·p = E + S + 2Mr(E + S)/Δ = (r² + a²)(E + S)/Δ`, with
`Δ = r² + a² − 2Mr` and `S = −aL/(r² + a²)`. -/
theorem jet_push_on_orbit (E S r a M : K) (hD : r^2 + a^2 - 2 * M * r ≠ 0) :
    E + S + 2 * M * r * (E + S) / (r^2 + a^2 - 2 * M * r)
      = (r^2 + a^2) * (E + S) / (r^2 + a^2 - 2 * M * r) := by
  grind

/-- **Round 21: two velocities' spin tower (check JB2).** The spin term times `(1 − 2u)(1 − 3u)` adds
`−40 χu^{5/2} + 48 χu^{7/2}` to Kerr's `8 χu^{3/2}`. -/
theorem two_velocity_spin_factor (u : K) : 8 * ((1 - 2 * u) * (1 - 3 * u)) = 8 - 40 * u + 48 * u^2 := by
  grind

/-- **Round 21: over the poles, two velocities turn the plane slower by `Δ/(r² + a²)` (check JB4).** -/
theorem polar_turning_ratio (r a M : K) (h : r^2 + a^2 ≠ 0) :
    (r^2 + a^2 - 2 * M * r) / (r^2 + a^2) = 1 - 2 * M * r / (r^2 + a^2) := by
  grind

/-! Round 22: one law for both wobbles of a circular equatorial orbit (wobble_law.py). v = √(M/r), q = a/r,
D = 1 − 3v² + 2qv, x = L − aE. -/

/-- In-and-out wobble: D − 3(v − q)² is Kerr's 1 − 6v² + 8qv − 3q². -/
theorem kerr_radial_wobble (v q : K) :
    (1 - 3 * v^2 + 2 * q * v) - 3 * (v - q)^2 = 1 - 6 * v^2 + 8 * q * v - 3 * q^2 := by
  grind

/-- Up-and-down wobble: D + 3(v − q)² is Kerr's 1 − 4qv + 3q². -/
theorem kerr_vertical_wobble (v q : K) :
    (1 - 3 * v^2 + 2 * q * v) + 3 * (v - q)^2 = 1 - 4 * q * v + 3 * q^2 := by
  grind

/-- The circular orbit's binding: D − (E√D)² = v²(1 − (2v − q)²), with E√D = 1 − 2v² + qv. -/
theorem circular_binding (v q : K) :
    (1 - 3 * v^2 + 2 * q * v) - (1 - 2 * v^2 + q * v)^2 = v^2 * (1 - (2 * v - q)^2) := by
  grind

/-- Mino-time radial stiffness on a circular orbit. (The library adds the hypothesis 2 ≠ 0: in characteristic 2
both orbit conditions lose their M-terms and the identity fails; the topic file's linear_combination uses 1/2.) With e = E² − 1 and K0 = L² + a²(1 − E²), the orbit sits where
R/r = r³e + 2Mr² − K0 r + 2Mx² and R' = 4r³e + 6Mr² − 2K0 r + 2Mx² vanish; then r·(−R''/2) = Mr² − 3Mx². -/
theorem mino_radial_stiffness (r e M K0 x : K)
    (h0 : r^3 * e + 2 * M * r^2 - K0 * r + 2 * M * x^2 = 0)
    (h1 : 4 * r^3 * e + 6 * M * r^2 - 2 * K0 * r + 2 * M * x^2 = 0)
    (h2 : (2 : K) ≠ 0) :
    r * K0 - 6 * r^3 * e - 6 * M * r^2 = M * r^2 - 3 * M * x^2 := by
  have h : 2 * ((r * K0 - 6 * r^3 * e - 6 * M * r^2) - (M * r^2 - 3 * M * x^2)) = 0 := by grind
  grind

/-- Mino-time up-and-down stiffness on the same orbit: r·K0 = Mr² + 3Mx² (with 2 ≠ 0, as above). -/
theorem mino_vertical_stiffness (r e M K0 x : K)
    (h0 : r^3 * e + 2 * M * r^2 - K0 * r + 2 * M * x^2 = 0)
    (h1 : 4 * r^3 * e + 6 * M * r^2 - 2 * K0 * r + 2 * M * x^2 = 0)
    (h2 : (2 : K) ≠ 0) :
    r * K0 = M * r^2 + 3 * M * x^2 := by
  have h : 2 * (r * K0 - (M * r^2 + 3 * M * x^2)) = 0 := by grind
  grind

/-! Round 23: Kerr's orbits are Newton's orbits around the ring-shifted mass plus one mass term
(polar_coupling.py). P = E(r² + a²) − aL, Δ = r² − 2Mr + a², Carter's K = Q + (L − aE)², s = sin θ, c = cos θ. -/

/-- **Radial: Kerr = Newton's ring + 2Mr(K − a²).** Carter's R = P² − Δ(r² + K) equals Newton's radial function
for the ring-shifted mass in Mino time, (r² + a²)((E² − 1)r² + 2Mr − (Q + L²)) + a²L², plus one term. -/
theorem kerr_radial_newton_plus_one (r a M E L Q : K) :
    (E * (r^2 + a^2) - a * L)^2 - (r^2 - 2 * M * r + a^2) * (r^2 + (Q + (L - a * E)^2)) =
      ((r^2 + a^2) * ((E^2 - 1) * r^2 + 2 * M * r - (Q + L^2)) + a^2 * L^2)
        + 2 * M * r * ((Q + (L - a * E)^2) - a^2) := by
  grind

/-- **Polar: Kerr's Θ is Newton's**, multiplied through by s²: Q − c²(a²(1 − E²) + L²/s²) = (Q + L²) +
(E² − 1)a²c² − L²/s². -/
theorem kerr_polar_is_newton (Q L a E s c : K) (hsc : s^2 + c^2 = 1) :
    Q * s^2 - c^2 * (a^2 * (1 - E^2) * s^2 + L^2) = (Q + L^2) * s^2 + (E^2 - 1) * a^2 * c^2 * s^2 - L^2 := by
  grind

/-- **Clock: Kerr = E × Newton's + 2MrP/Δ**, multiplied through by Δ. -/
theorem kerr_clock_newton_plus_one (r a M E L s c : K) (hsc : s^2 + c^2 = 1) :
    (r^2 + a^2) * (E * (r^2 + a^2) - a * L) - (r^2 - 2 * M * r + a^2) * (a * (a * E * s^2 - L)) =
      (r^2 - 2 * M * r + a^2) * (E * (r^2 + a^2 * c^2)) + 2 * M * r * (E * (r^2 + a^2) - a * L) := by
  grind

/-- **Turning: Kerr = Newton's + (a/(r² + a²))·2MrP/Δ**, multiplied through by (r² + a²)Δ: the L/s² parts cancel. -/
theorem kerr_turning_newton_plus_one (r a M E L : K) :
    a * (E * (r^2 + a^2) - a * L) * (r^2 + a^2) - a * E * (r^2 - 2 * M * r + a^2) * (r^2 + a^2)
        + a^2 * L * (r^2 - 2 * M * r + a^2) = 2 * M * r * a * (E * (r^2 + a^2) - a * L) := by
  grind

/-- Carter's R in Mino form for every orbit: r⁴(E² − 1) + 2Mr³ − W r² + 2MrK − a²Q, W = Q + L² − a²(E² − 1). -/
theorem kerr_radial_mino_form (r a M E L Q : K) :
    (E * (r^2 + a^2) - a * L)^2 - (r^2 - 2 * M * r + a^2) * (r^2 + (Q + (L - a * E)^2)) =
      r^4 * (E^2 - 1) + 2 * M * r^3 - (Q + L^2 - a^2 * (E^2 - 1)) * r^2 + 2 * M * r * (Q + (L - a * E)^2)
        - a^2 * Q := by
  grind

/-- **The stiffness law for every spherical orbit** (round 22's, now with Q): if R and R' vanish at r, then
r²·(−R''/2) = Mr³ − 3MKr + 4a²Q. Needs 2 ≠ 0. -/
theorem spherical_radial_stiffness (r e M W Kc a Q : K) (h2 : (2 : K) ≠ 0)
    (h0 : r^4 * e + 2 * M * r^3 - W * r^2 + 2 * M * r * Kc - a^2 * Q = 0)
    (h1 : 4 * r^3 * e + 6 * M * r^2 - 2 * W * r + 2 * M * Kc = 0) :
    r^2 * (W - 6 * r^2 * e - 6 * M * r) = M * r^3 - 3 * M * Kc * r + 4 * a^2 * Q := by
  have h : 2 * (r^2 * (W - 6 * r^2 * e - 6 * M * r) - (M * r^3 - 3 * M * Kc * r + 4 * a^2 * Q)) = 0 := by grind
  grind

/-- ... and r²W = Mr³ + 3MKr − 2a²Q (the r² coefficient; over the poles W = Q(1 + m)). Needs 2 ≠ 0. -/
theorem spherical_vertical_stiffness (r e M W Kc a Q : K) (h2 : (2 : K) ≠ 0)
    (h0 : r^4 * e + 2 * M * r^3 - W * r^2 + 2 * M * r * Kc - a^2 * Q = 0)
    (h1 : 4 * r^3 * e + 6 * M * r^2 - 2 * W * r + 2 * M * Kc = 0) :
    r^2 * W = M * r^3 + 3 * M * Kc * r - 2 * a^2 * Q := by
  have h : 2 * (r^2 * W - (M * r^3 + 3 * M * Kc * r - 2 * a^2 * Q)) = 0 := by grind
  grind

/-- Over the poles at order (a/r)², doubled: 3(1 − 4u)(1 − 6u) + 12u(5 − 8u) = 3(1 + 10u − 8u²), i.e.
(3/2)(1 − 4u)(1 − 6u) + 6u(5 − 8u) = (3/2)(1 + 10u − 8u²) wherever 2 ≠ 0. -/
theorem polar_k2_law (u : K) :
    3 * (1 - 4 * u) * (1 - 6 * u) + 12 * u * (5 - 8 * u) = 3 * (1 + 10 * u - 8 * u^2) := by
  grind

/-! ### Round 25: the one mass term, piece by piece

The Kerr–Schild shell with the field entering through (k·p)², k·p = E + P + S (time part, push, twist). -/

/-- **The river's time part alone is Newton's ring with mass M E²**: −(r² + a²)·C_r, with C_r multiplied through by
r² + a², is Newton's radial function with M replaced by M E². -/
theorem time_part_is_heavier_newton (r a M E L C : K) :
    -(-(a^2 * L^2) - (r^2 + a^2) * ((E^2 - 1) * r^2) + (r^2 + a^2) * C - 2 * M * r * E^2 * (r^2 + a^2)) =
      (r^2 + a^2) * ((E^2 - 1) * r^2 + 2 * (M * E^2) * r - C) + a^2 * L^2 := by
  grind

/-- **The one term is shared three ways** (multiplied through by r² + a²): time part 2Mr(E² − 1)(r² + a²), push
2Mr[C − (E² − 1)r² − a²L²/(r² + a²)], twist 2Mr[a²L²/(r² + a²) − 2aEL], with K = C − L² + (L − aE)². -/
theorem one_term_split (r a M E L C : K) :
    (r^2 + a^2) * (2 * M * r * ((C - L^2 + (L - a * E)^2) - a^2)) =
      (r^2 + a^2) * (2 * M * r * (E^2 - 1) * (r^2 + a^2))
        + 2 * M * r * (C * (r^2 + a^2) - (E^2 - 1) * r^2 * (r^2 + a^2) - a^2 * L^2)
        + 2 * M * r * (a^2 * L^2 - 2 * a * E * L * (r^2 + a^2)) := by
  grind

/-- **The push's second-order parts cancel when the pieces form one square**: with weights from (αE + βP + γS)²,
the B² term (from E·P and P·S) equals P·P times the rest. -/
theorem river_square_cancels (α β γ E S : K) :
    (α * β * E + β * γ * S)^2 = β^2 * (α^2 * E^2 + γ^2 * S^2 + 2 * (α * γ) * E * S) := by
  grind

/-- **The clock term split** (multiplied through by (r² + a²)Δ, Δ = r² + a² − 2Mr): time 2MrE, twist 2MrS,
push 4M²r²E/Δ, push–twist 4M²r²S/Δ, with S = −aL/(r² + a²). -/
theorem clock_split (r a M E L : K) :
    2 * M * r * (E * (r^2 + a^2) - a * L) * (r^2 + a^2) =
      2 * M * r * E * ((r^2 + a^2) * (r^2 + a^2 - 2 * M * r)) + 2 * M * r * (-(a * L)) * (r^2 + a^2 - 2 * M * r)
        + 4 * M^2 * r^2 * E * (r^2 + a^2) + 4 * M^2 * r^2 * (-(a * L)) := by
  grind

/-- **Per unit energy, Kerr is Newton's ring + push + twist**: with C̃ = C/E², L̃ = L/E and Newton's energy
1 − 1/E², E² times [Newton (mass M) + 2Mr(C̃ − (1 − 1/E²)r² − 2aL̃)] is Carter's radial function. -/
theorem kerr_per_unit_energy (r a M E L C : K) (hE : E ≠ 0) :
    E^2 * ((r^2 + a^2) * ((1 - 1 / E^2) * r^2 + 2 * M * r - C / E^2) + a^2 * (L / E)^2
        + 2 * M * r * (C / E^2 - (1 - 1 / E^2) * r^2 - 2 * a * (L / E))) =
      (E * (r^2 + a^2) - a * L)^2 - (r^2 - 2 * M * r + a^2) * (r^2 + (C - L^2 + (L - a * E)^2)) := by
  grind

/-- **Round 26, per unit energy at the equator**: Newton's ring plus push plus twist adds up to
(r² + a²)(e r⁴ + 2M(1 − e) r³ − (L² − e a²) r² + 2M r (L − a)²). -/
theorem equator_split_per_unit_energy (r a M e L : K) :
    (r^2 + a^2) * (e * r^4 + 2 * M * r^3 - (L^2 - e * a^2) * r^2 + 2 * M * r * a^2)
        + ((r^2 + a^2) * (-(2 * M * e * r^3) + 2 * M * r * L^2) - 2 * M * r * a^2 * L^2)
        + ((r^2 + a^2) * (-(4 * M * r * a * L)) + 2 * M * r * a^2 * L^2) =
      (r^2 + a^2) * (e * r^4 + 2 * M * (1 - e) * r^3 - (L^2 - e * a^2) * r^2 + 2 * M * r * (L - a)^2) := by
  grind

/-- The square's three sources: the ring's a², the push's L², the twist's −2aL. -/
theorem square_from_three_sources (a L : K) : a^2 + L^2 - 2 * a * L = (L - a)^2 := by
  grind

/-- **The stiffness law for any radial function** e r⁴ + 2αr³ − W r² + 2βr on a circular orbit (R/r = 0 and
R' = 0 at r): r·(−R''/2) = αr² − 3β. Needs 2 ≠ 0. Kerr per unit energy: α = M(1 − e), β = M(L − a)². -/
theorem generic_radial_stiffness (r e α W β : K) (h2 : (2 : K) ≠ 0)
    (h0 : e * r^3 + 2 * α * r^2 - W * r + 2 * β = 0)
    (h1 : 4 * e * r^3 + 6 * α * r^2 - 2 * W * r + 2 * β = 0) :
    r * (W - 6 * e * r^2 - 6 * α * r) = α * r^2 - 3 * β := by
  have h : 2 * (r * (W - 6 * e * r^2 - 6 * α * r) - (α * r^2 - 3 * β)) = 0 := by grind
  grind

/-- ... and rW = αr² + 3β: the same two pieces with the opposite sign of 3β. -/
theorem generic_vertical_stiffness (r e α W β : K) (h2 : (2 : K) ≠ 0)
    (h0 : e * r^3 + 2 * α * r^2 - W * r + 2 * β = 0)
    (h1 : 4 * e * r^3 + 6 * α * r^2 - 2 * W * r + 2 * β = 0) :
    r * W = α * r^2 + 3 * β := by
  have h : 2 * (r * W - (α * r^2 + 3 * β)) = 0 := by grind
  grind

/-- **Round 22's law, term by term**: Newton's ring's 1 − 3q², the push's −6v² and the twist's +8qv make
D − 3(v − q)², with D = 1 − 3v² + 2qv. -/
theorem equator_ratio_three_sources (v q : K) :
    (1 - 3 * q^2) + (-(6 * v^2)) + 8 * q * v = (1 - 3 * v^2 + 2 * q * v) - 3 * (v - q)^2 := by
  grind

/-- **Kerr–Newman is Newton's ring around the mass function plus one term**: Carter's Kerr–Newman radial function
equals Newton's spinning ring with 2mr = 2Mr − Q², plus (2Mr − Q²)(K − a²), K = C − L² + (L − aE)². -/
theorem kn_ring_plus_one_term (r a M Q E L C : K) :
    (E * (r^2 + a^2) - a * L)^2 - (r^2 - 2 * M * r + a^2 + Q^2) * (r^2 + (C - L^2 + (L - a * E)^2)) =
      (r^2 + a^2) * ((E^2 - 1) * r^2 + (2 * M * r - Q^2) - C) + a^2 * L^2
        + (2 * M * r - Q^2) * ((C - L^2 + (L - a * E)^2) - a^2) := by
  grind

/-- **Round 26's equator split with charge** (multiplied through by D = r² + a², f = 2Mr − Q²): Newton's ring
e r⁴ + f r² − (L² − e a²) r² + f a², the push −f e r² + f L² − f a²L²/D and the twist −2faL + f a²L²/D add up to
e r⁴ + f(1 − e) r² − (L² − e a²) r² + f(L − a)². -/
theorem kn_equator_split (r a M Q e L : K) :
    (r^2 + a^2) * (e * r^4 + (2 * M * r - Q^2) * r^2 - (L^2 - e * a^2) * r^2 + (2 * M * r - Q^2) * a^2)
        + ((r^2 + a^2) * (-((2 * M * r - Q^2) * e * r^2) + (2 * M * r - Q^2) * L^2)
            - (2 * M * r - Q^2) * a^2 * L^2)
        + ((r^2 + a^2) * (-(2 * (2 * M * r - Q^2) * a * L)) + (2 * M * r - Q^2) * a^2 * L^2) =
      (r^2 + a^2) * (e * r^4 + (2 * M * r - Q^2) * (1 - e) * r^2 - (L^2 - e * a^2) * r^2
        + (2 * M * r - Q^2) * (L - a)^2) := by
  grind

/-- **The stiffness law with a constant term** e r⁴ + 2αr³ − W r² + 2βr − γ (R = R' = 0 at r):
r²·(−R''/2) = αr³ − 3βr + 4γ. Kerr–Newman per unit energy: α = M(1 − e), β = M(L − a)², γ = Q²(L − a)².
(Round 23's spherical law with a²K → γ.) Needs 2 ≠ 0. -/
theorem charged_radial_stiffness (r e α W β γ : K) (h2 : (2 : K) ≠ 0)
    (h0 : e * r^4 + 2 * α * r^3 - W * r^2 + 2 * β * r - γ = 0)
    (h1 : 4 * e * r^3 + 6 * α * r^2 - 2 * W * r + 2 * β = 0) :
    r^2 * (W - 6 * e * r^2 - 6 * α * r) = α * r^3 - 3 * β * r + 4 * γ := by
  have h : 2 * (r^2 * (W - 6 * e * r^2 - 6 * α * r) - (α * r^3 - 3 * β * r + 4 * γ)) = 0 := by grind
  grind

/-- **Kerr–Newman's equatorial law, three sources**: Newton's ring 1 − 3q² + s²(1 + q²)/V², the push
−6v² + 3s² − s⁴/V² and the twist 8qV (v² = V² + s²) make Kerr's law in V plus one charge term. -/
theorem kn_ratio_three_sources (v V q s : K) (hV : V ≠ 0) (hv : v^2 = V^2 + s^2) :
    (1 - 3 * q^2 + s^2 * (1 + q^2) / V^2) + (-(6 * v^2) + 3 * s^2 - s^4 / V^2) + 8 * q * V =
      (1 - 3 * V^2 + 2 * q * V) - 3 * (V - q)^2 + s^2 / V^2 * (1 + q^2 - s^2 - 3 * V^2) := by
  rw [hv]
  grind

/-- **Light past the hole, time part alone** (round 28, a = 0, null ray, E = 1): with A = r² and
C = −r² + b² − 2Mr the radial function −A·C is r²(r² + 2Mr − b²), Newton's pull on light. Its only turning
points are where r² + 2Mr = b², so it has no photon sphere and no shadow; the bending is 2 arctan(M/b). -/
theorem time_part_light (r M b : K) :
    -(r^2) * (-(r^2) + b^2 - 2 * M * r) = r^2 * (r^2 + 2 * M * r - b^2) := by
  grind

/-- **The push makes the photon sphere** (round 28, a = 0): time part + push give Schwarzschild's
R = r³ − b²(r − 2M) (divided by r). A circular light orbit (R = R' = 0, r ≠ 0) sits at r = 3M.
Needs 2 ≠ 0. -/
theorem push_photon_sphere (r b M : K) (h2 : (2 : K) ≠ 0) (hr : r ≠ 0)
    (h0 : r^3 - b^2 * (r - 2 * M) = 0) (h1 : 3 * r^2 - b^2 = 0) : r = 3 * M := by
  have h : r^2 * (2 * (r - 3 * M)) = 0 := by grind
  grind

/-- **The push sets the shadow's size** (round 28, a = 0): at the photon sphere r = 3M the critical
impact parameter satisfies b² = 3r² = 27M², so the shadow's radius on the sky is √27 M. -/
theorem push_shadow_radius (r b M : K) (h1 : 3 * r^2 - b^2 = 0) (hr3 : r - 3 * M = 0) :
    b^2 = 27 * M^2 := by
  grind

/-- **What the push adds to light** (round 29, no spin): Schwarzschild's null radial function r⁴ − b²r² + 2Mb²r is
the time part's r⁴ + 2Mr³ − b²r² (Newton's light) plus the push's −2Mr(r² − b²). -/
theorem push_null_term (r M b : K) :
    r^4 - b^2 * r^2 + 2 * M * b^2 * r = (r^4 + 2 * M * r^3 - b^2 * r^2) - 2 * M * r * (r^2 - b^2) := by
  grind

/-- **One dial, two tests** (round 29): with the push at strength λ (k·p = E + λP), the perihelion factor
(2 + 2γ − β)/3 is λ² and the delay gives γ = λ², so β + γ = 2 whatever λ is. -/
theorem one_dial_beta_gamma (β γ l : K) (hp : 2 + 2 * γ - β = 3 * l^2) (hg : γ = l^2) : β + γ = 2 := by
  grind

/-- **The square ties perihelion to the push** (round 29): the perihelion coefficient 2(2μ² + ν) from the EP weight μ
and the PP weight ν becomes 6λ² when the field is one square (μ = λ, ν = λ²). -/
theorem square_dial_perihelion (l : K) : 2 * (2 * l^2 + l^2) = 6 * l^2 := by
  grind

/-- **The third-level phase weight** (round 30): a change c₆ in the v⁶ term of E′/F moves the 3PN phase coefficient by
40 c₆/((6 − 8)(6 − 5)) = −20 c₆, and a change δe₃ in the energy moves c₆ by 4 δe₃, so φ₆ moves by −80 δe₃. -/
theorem third_level_phase_weight (de3 : K) :
    (6 - 8) * (6 - 5) * (-20 : K) = 40 ∧ -20 * (4 * de3) = -80 * de3 := by
  constructor <;> grind

/-- **The π² floor** (round 30): missing only relativity's −(205/96)π²ν in the energy moves φ₆ by −(1025/6)π²ν
(stated with the 96 and the 6 cleared, so it holds in any field). -/
theorem pi_floor_phase (p2 n : K) : 6 * (-80 * (205 * p2 * n)) = 96 * (-1025 * p2 * n) := by
  grind

/-- **The boost-safe shape** (round 31): a π² term whose velocity weights satisfy 2α + γ = 0 = 2β + γ depends on
the bodies' velocities only through v₁ − v₂, so on a circle it carries the relative speed alone. -/
theorem boost_safe_shape (α v1 v2 : K) : α * v1^2 + α * v2^2 + (-2 * α) * v1 * v2 = α * (v1 - v2)^2 := by
  grind

/-- **The time part alone on a circle** (round 31): field energy's time part has weights −1/8, −1/8, 1; on a circle
(v₁ = m₂v, v₂ = −m₁v, m₁ + m₂ = 1) eight times its bracket is −1 − 6ν with ν = m₁m₂, so it would put π² into ν²
as well, where relativity has none. -/
theorem time_part_circle (m1 : K) :
    -((1 - m1)^2 + m1^2) - 8 * (m1 * (1 - m1)) = -1 - 6 * (m1 * (1 - m1)) := by
  grind

/-- **The lag cancels the shape** (round 32): with D₁, D₂ the two bodies' velocity derivatives (written x, y),
eight times the time part, −(x² + y²) − 8xy, plus four times the full time derivative (x + y)², is 3(x² + y²):
the time part is (3/8)(D₁² + D₂²) up to a total time derivative, the piece the lag then removes. -/
theorem lag_cancels_shape (x y : K) : -(x^2 + y^2) - 8 * x * y + 4 * (x + y)^2 = 3 * (x^2 + y^2) := by
  grind

/-- **What is left is a time derivative** (round 32): eight times (time part + lag), −(x² + y²) − 8xy − 3(x² + y²),
is −4(x + y)², i.e. the sum is −½(D₁ + D₂)², the velocity part of −½ d²F/dt². -/
theorem lagged_sum_is_time_derivative (x y : K) :
    (-(x^2 + y^2) - 8 * x * y) + (-3 * (x^2 + y^2)) = -4 * (x + y)^2 := by
  grind

end Field

end CircleOp.SpinningRiver
