/-
  Laws from the hypersphere: the exact algebra behind the hypersphere page and its orbits section
  ===============================================================================================

  DRAFT. Companion checks: ../checks.py (57 checks) and ../orbits/checks_orbits.py (21 checks).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports Base only.
  First compiled in the shared library on Nov 26. Every proof is `grind`, `decide` or a short
  induction; the facts that evaluate Rat definitions or Rat numerals (sections 7 to 9) use
  `unfold …; grind` or `grind`, because `decide` can't evaluate Rat arithmetic in core Lean.
  Names updated to 1.0 (CircleOp.circ, called MultAdd.multAdd before 1.0).

  What is here, and what is known vs ours:
    1. Quaternions (known: Hamilton 1843; Euler's four-square identity 1748).
       A pure-imaginary unit squares to −1 (every tangent direction at 1 on S³ is an "i"),
       ij = k = −ji, and the norm multiplies, which is what lets u = r·ω split.
    2. Fock's stereographic map (known: Fock 1935) puts each momentum on the unit 3-sphere.
    3. Kepler orbits as tilted great circles (known): with p0² = −2E, (h·p0)² + e² = 1.
    4. Level n holds n² orbital states (known): 1 + 3 + … + (2n − 1) = n².
    5. Electric orbits (ours as an exact equality; implied in qualitative form by Biedenharn 1983,
       not stated as numbers in anything we read): the circular orbit with the same L moves at speed k/L, so the disc
       of speed-radius k/L is that orbit's velocity, and Sommerfeld's advance per lap equals
       its Thomas turn per lap.
    6. Hydrogen (exact rewrite of Dirac, ours as a rewrite): the quantum defect δ = κx/(1+x),
       with x the disc area per lap over 2π, gives (κ − δ)² = κ² − α², Dirac's defect.
    7. Gravity (known PPN formulas, our reading): in units of hodograph area / c² the perihelion
       advance is 2(2 + 2γ − β) and the geodetic turn is 1 + 2γ. GR: 6 and 3, so the orbit turns
       twice its gyroscope. Nordström: −1 and −1, equal, like the electric case.
       First order only; ../orbits/checks_orbits.py P5b shows GR's 6 drifts in strong fields.
    8. Gravity's 6, split exactly (ours as a split; the frequencies are textbook; ../orbits/six/checks_six.py):
       for a near-circular Schwarzschild orbit, (radial ÷ angular frequency)² = 1 − 6M/r =
       (1 − 2M/r)(1 − 4v²), v the local orbital speed, beside the gyroscope's 1 − 3M/r = (1 − 2M/r)(1 − v²).
       Same space factor; the orbit's velocity factor runs at 2v. Both fuse by the circle operation.
       The general rule T = (N − 4v²)/B for any static metric is checked symbolically in sympy, not here.
       The 4v² is v² (special relativity) + 3v² (gravity pulls on energy): ../orbits/twov/checks_twov.py.
-/
import FusedAlgebra.Base

namespace Hypersphere
open CircleOp

local infixl:65 " ✳ " => circ

/-! ## 1. Quaternions over the integers -/

structure Q where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq, Repr

def mul (p q : Q) : Q :=
  { a := p.a * q.a - p.b * q.b - p.c * q.c - p.d * q.d
    b := p.a * q.b + p.b * q.a + p.c * q.d - p.d * q.c
    c := p.a * q.c - p.b * q.d + p.c * q.a + p.d * q.b
    d := p.a * q.d + p.b * q.c - p.c * q.b + p.d * q.a }

def norm (p : Q) : Int := p.a * p.a + p.b * p.b + p.c * p.c + p.d * p.d

def qi : Q := ⟨0, 1, 0, 0⟩
def qj : Q := ⟨0, 0, 1, 0⟩
def qk : Q := ⟨0, 0, 0, 1⟩
def neg (p : Q) : Q := ⟨-p.a, -p.b, -p.c, -p.d⟩

theorem i_mul_j : mul qi qj = qk := by decide
theorem j_mul_i : mul qj qi = neg qk := by decide
theorem i_sq : mul qi qi = ⟨-1, 0, 0, 0⟩ := by decide

/-- **Every tangent direction at 1 squares to minus its length squared**: on the unit
    3-sphere each unit tangent at 1 is an "i" (a whole 2-sphere of them). -/
theorem tangent_sq (x y z : Int) :
    mul ⟨0, x, y, z⟩ ⟨0, x, y, z⟩ = ⟨-(x * x + y * y + z * z), 0, 0, 0⟩ := by
  simp only [mul, Q.mk.injEq]; grind

/-- **The norm multiplies** (Euler's four-square identity): what makes the polar split
    u = r·ω work on S³. -/
theorem norm_mul (p q : Q) : norm (mul p q) = norm p * norm q := by
  simp only [norm, mul]; grind

/-! ## 2. Fock's map: a momentum p with |p|² = s lands on the unit 3-sphere -/

/-- The point (2 p0 p, s − p0²)/(s + p0²) has length 1: |2 p0 p|² + (s − p0²)² = (s + p0²)². -/
theorem fock_on_sphere (s p0 : Rat) :
    4 * p0 * p0 * s + (s - p0 * p0) * (s - p0 * p0) = (s + p0 * p0) * (s + p0 * p0) := by
  grind

/-! ## 3. Kepler orbits are great circles, tilted by the eccentricity -/

/-- With k = m = 1, e² = 1 + 2 E h². Set p0² = −2E (the Fock sphere's radius squared): then
    (h p0)² + e² = 1, so (h p0, e) is the cosine and sine of the great circle's tilt. -/
theorem kepler_tilt (E h2 e2 : Rat) (h : e2 = 1 + 2 * E * h2) :
    h2 * (-2 * E) + e2 = 1 := by
  grind

/-! ## 4. Level n holds n² orbital states -/

/-- Number of orbital states up to level n: ℓ = 0 … n − 1, each with 2ℓ + 1. -/
def deg : Nat → Nat
  | 0 => 0
  | n + 1 => deg n + (2 * n + 1)

theorem deg_eq (n : Nat) : deg n = n * n := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [deg, ih]; grind

/-! ## 5. Electric orbits: the disc radius is the circular orbit's speed -/

/-- Circular orbit under k/r² (k = m = c = 1): γ v² r = 1 and L = γ v r give L·v = 1,
    so the disc's speed-radius k/L is exactly the speed v of the circular orbit with that L. -/
theorem circular_L_v (g v r : Rat) (h : g * v * v * r = 1) : (g * v * r) * v = 1 := by
  grind

/-- Hence Sommerfeld's factor 1/(1 − (k/L)²) equals the circular orbit's γ² = 1/(1 − v²):
    the orbit's advance per lap equals that orbit's Thomas turn per lap. -/
theorem sommerfeld_eq_thomas (L v : Rat) (h : L * v = 1) :
    1 - (1 / L) * (1 / L) = 1 - v * v := by
  have hL : L ≠ 0 := by intro h0; rw [h0] at h; grind
  have e : 1 / L = v := by grind
  rw [e]

/-! ## 6. Hydrogen: Dirac's defect from the disc area -/

/-- Let f = 1/√(1 − α²/κ²) (so f²(κ² − α²) = κ²), x = f − 1 the disc area per lap over 2π,
    and δ = κx/(1 + x), i.e. δ f = κ(f − 1). Then (κ − δ)² = κ² − α²: Dirac's quantum defect,
    so E = mc²/√(1 + α²/(n − δ)²) is Dirac's fine structure. -/
theorem dirac_defect (κ α δ f : Rat) (hf0 : f ≠ 0)
    (hδ : δ * f = κ * (f - 1)) (hf : f * f * (κ * κ - α * α) = κ * κ) :
    (κ - δ) * (κ - δ) = κ * κ - α * α := by
  have h1 : (κ - δ) * f = κ := by grind
  have h2 : f * f * ((κ - δ) * (κ - δ)) = f * f * (κ * κ - α * α) := by
    have : f * f * ((κ - δ) * (κ - δ)) = ((κ - δ) * f) * ((κ - δ) * f) := by grind
    rw [this, h1, hf]
  have hff : f * f ≠ 0 := by
    intro h0
    have : f = 0 := by grind
    exact hf0 this
  grind

/-! ## 7. Gravity: the coefficients in units of hodograph area / c² (first order) -/

/-- PPN perihelion advance per orbit, 2π(2 + 2γ − β)M/p, in units of area/c² = πM/p. -/
def perihelionCoef (γ β : Rat) : Rat := 2 * (2 + 2 * γ - β)

/-- PPN geodetic turn per circular orbit, 2π(γ + ½)M/r, in units of area/c² = πM/r. -/
def geodeticCoef (γ : Rat) : Rat := 1 + 2 * γ

theorem gr_perihelion : perihelionCoef 1 1 = 6 := by unfold perihelionCoef; grind
theorem gr_geodetic : geodeticCoef 1 = 3 := by unfold geodeticCoef; grind
theorem nordstrom_perihelion : perihelionCoef (-1) (1 / 2) = -1 := by unfold perihelionCoef; grind
theorem nordstrom_geodetic : geodeticCoef (-1) = -1 := by unfold geodeticCoef; grind

/-- **GR's orbit turns twice as much as its gyroscope** (first order). -/
theorem gr_orbit_twice_gyro : perihelionCoef 1 1 = 2 * geodeticCoef 1 := by unfold perihelionCoef geodeticCoef; grind

/-- In Nordström's scalar gravity the orbit turns exactly as much as its gyroscope
    (first order), as the electric orbit does exactly (§5). -/
theorem nordstrom_orbit_eq_gyro : perihelionCoef (-1) (1 / 2) = geodeticCoef (-1) := by unfold perihelionCoef geodeticCoef; grind

/-- In general: perihelion − geodetic = 3 + 2γ − 2β. It vanishes for Nordström and is 3 for GR. -/
theorem orbit_minus_gyro (γ β : Rat) :
    perihelionCoef γ β - geodeticCoef γ = 3 + 2 * γ - 2 * β := by
  unfold perihelionCoef geodeticCoef; grind

/-! ## 8. Gravity's 6, split exactly -/

/-- Local orbital speed squared of a circular Schwarzschild orbit, x = M/r: v² = x/(1 − 2x). -/
def vsq (x : Rat) : Rat := x / (1 - 2 * x)

/-- **The orbit's split**: (1 − 2x)(1 − 4v²) = 1 − 6x. -/
theorem orbit_split (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    (1 - 2 * x) * (1 - 4 * vsq x) = 1 - 6 * x := by
  unfold vsq; grind

/-- **The gyroscope's split** (as on the Gravity Probe B page): (1 − 2x)(1 − v²) = 1 − 3x. -/
theorem gyro_split (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    (1 - 2 * x) * (1 - vsq x) = 1 - 3 * x := by
  unfold vsq; grind

/-- The orbit's fractions fuse by the circle operation: −6x = (−2x) ✳ (−4v²). -/
theorem orbit_fractions_fuse (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    -(6 * x) = (-(2 * x)) ✳ (-(4 * vsq x)) := by
  show -(6 * x) = -(2 * x) + -(4 * vsq x) + -(2 * x) * -(4 * vsq x)
  unfold vsq; grind

/-- At the innermost stable orbit r = 6M the orbit's disc radius 2v reaches light speed. -/
theorem isco_two_v : 4 * vsq (1 / 6) = 1 := by unfold vsq; grind

/-- **Where the 2v comes from** (../orbits/twov): the energy potential E(1 − 1/√(1 − 2x)) has Newton's closure
    number N = (1 − 5x)/(1 − 2x), which is 1 − 3v². -/
theorem energy_pull_closure (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    (1 - 5 * x) / (1 - 2 * x) = 1 - 3 * vsq x := by
  unfold vsq; grind

/-- So 1 − 6x = (1 − 2x)·[(1 − 3v²) − v²]: space, times the energy pull (3v²) and special relativity (v²),
    which add to 4v² = (2v)². -/
theorem six_is_two_one_three (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    (1 - 2 * x) * ((1 - 3 * vsq x) - vsq x) = 1 - 6 * x := by
  unfold vsq; grind

/-- At the photon sphere r = 3M the gyroscope's disc radius v reaches light speed. -/
theorem photon_sphere_v : vsq (1 / 3) = 1 := by unfold vsq; grind

/-! ## 9. Gravity's 6 through the double copy (Kerr-Schild form, ../orbits/dcopy)

A probe in Schwarzschild, in Kerr-Schild form, has mass shell η(p,p) − (2M/r)(k·p)² = −1, with k·p = E + p_r.
The single copy (Coulomb) has η(p,p) − (2Q/r)(k·p) = −1: one power of k·p, not two. -/

/-- The square has three pieces: Newton's pull E², the cross term 2E·p, and space p². -/
theorem kp_square (e p : Rat) : (e + p) * (e + p) = e * e + 2 * (e * p) + p * p := by grind

/-- First-order turn per lap, in units of π M²/L², when the three pieces carry weights a, b, c
    (checked by simulated orbits, dcopy K7–K8). -/
def ksCoef (a b c : Rat) : Rat := 4 * (b * b) + 2 * (a * c)

theorem ks_all : ksCoef 1 1 1 = 6 := by unfold ksCoef; grind
theorem ks_newton_only : ksCoef 1 0 0 = 0 := by unfold ksCoef; grind
theorem ks_newton_space : ksCoef 1 0 1 = 2 := by unfold ksCoef; grind
theorem ks_newton_cross : ksCoef 1 1 0 = 4 := by unfold ksCoef; grind

/-- The single copy's cross term has strength 1 (Sommerfeld's 1); gravity's has strength 2, from the 2 in 2E·p,
    so it gives 2² = 4. With space's 2 that is the 6. -/
theorem six_is_space_plus_square : (2 : Rat) + 2 * 2 * 1 = 6 := by grind

/-- The clock-warp reading's special relativity (1) and energy pull (3) are the double copy's one square. -/
theorem one_plus_three_is_square : (1 : Rat) + 3 = 2 * 2 := by grind

/-- Completing the square in p_r (a canonical shift) turns the Kerr-Schild shell into a static one with
    1/A = 1 + φ + φ²/(1 − φ), which is 1/(1 − φ): Schwarzschild's clock factor, φ = 2M/r. -/
theorem ks_clock (φ : Rat) (h : 1 - φ ≠ 0) : 1 + φ + φ * φ / (1 - φ) = 1 / (1 - φ) := by grind

end Hypersphere
