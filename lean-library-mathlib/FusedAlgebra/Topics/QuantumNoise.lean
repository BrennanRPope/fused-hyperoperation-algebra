module

public import FusedAlgebra.Basic

/-!
# Quantum noise through the mult-add cube

The noise channels of a qubit compose by weighted mult-add `wAdd k a b = a + b + k * a * b`:

* `k = -1`: amplitude damping (T1 decay), "any fault counts" (`amplitude_damping`);
* `k = -2`: bit-flip and phase-flip, "only an odd number of faults counts" (`bit_flip`,
  `stacked_bit_flips`);
* `k = -4 / 3`: depolarizing noise on one qubit (`depolarizing`).

## Main definitions and results

* `QuantumNoise.Pauli`: a Pauli channel, the four chances of doing nothing, `X`, `Y` or `Z`, with
  composition `Pauli.comp` and the three coordinates `Pauli.lx`, `Pauli.ly`, `Pauli.lz` of its
  point in the cube `[-1, 1]³`. `comp_lx`, `comp_ly`, `comp_lz`: composition multiplies the
  coordinates. `faces`: the four faces of the tetrahedron are the four chances.
* `maj3`, `vote`: the 3-qubit repetition code's failure chance and its form in the polarization
  `l = 1 - 2 * p`; `vote_helps`: between the wall and the clean corner the vote always helps.
* `gx_below`, `gx_above`, `gz_below`, `gz_above`: the thresholds of the two-level concatenated
  Shor code, bracketed between explicit rationals.
* `macwilliams_rep3`: the MacWilliams identity for the 3-repetition code.

`Pauli.comp` is multiplication in the group algebra of the Klein four-group, and the coordinates
`lx`, `ly`, `lz` are its three non-trivial characters; this file does not make that connection to
mathlib's `MonoidAlgebra`. That these formulas describe the quantum channels (density matrices)
is not checked here. The core-Lean original worked over `ℚ`; here the identities hold in any
commutative ring, or any field of characteristic zero where a division appears.
-/

@[expose] public section

namespace QuantumNoise

open MultAdd

variable {R : Type*} {K : Type*}

/-! ### The three channel rules -/

section CommRing

variable [CommRing R]

/-- Amplitude damping: survival `1 - g` multiplies, the `k = -1` rule. -/
theorem amplitude_damping (g h : R) : 1 - wAdd (-1) g h = (1 - g) * (1 - h) := by
  unfold wAdd; ring

/-- Bit-flip (and phase-flip): the polarization `1 - 2 * p` multiplies, the `k = -2` rule. -/
theorem bit_flip (p q : R) : 1 - 2 * wAdd (-2) p q = (1 - 2 * p) * (1 - 2 * q) := by
  unfold wAdd; ring

/-- The rule `p + q - p * q` is the `k = -1` rule. It is right for amplitude damping. For
bit-flips it is off by exactly `p * q`. -/
theorem call_rule_gap (p q : R) : wAdd (-1) p q - wAdd (-2) p q = p * q := by
  unfold wAdd; ring

/-- `n` stacked bit-flips: the chance `P` of an odd number of flips satisfies
`1 - 2 * P = (1 - 2 * p) ^ n`. -/
theorem stacked_bit_flips (p : R) (n : ℕ) : 1 - 2 * wpow (-2) p n = (1 - 2 * p) ^ n := by
  have h := one_add_mul_wpow (-2) p n
  have e : (1 : R) + -2 * p = 1 - 2 * p := by ring
  rw [e] at h
  rw [← h]; ring

end CommRing

/-- Depolarizing noise on one qubit: `1 - 4 / 3 * p` multiplies, the `k = -4 / 3` rule. -/
theorem depolarizing [Field K] (p q : K) :
    1 - 4 / 3 * wAdd (-4 / 3) p q = (1 - 4 / 3 * p) * (1 - 4 / 3 * q) := by
  unfold wAdd; ring

/-! ### Pauli channels: the tetrahedron in the cube -/

/-- A Pauli channel: the chances of doing nothing, `X`, `Y` or `Z`. -/
@[ext]
structure Pauli (R : Type*) where
  /-- The chance of doing nothing. -/
  i : R
  /-- The chance of an `X` error. -/
  x : R
  /-- The chance of a `Y` error. -/
  y : R
  /-- The chance of a `Z` error. -/
  z : R

namespace Pauli

variable [CommRing R]

/-- The `X` coordinate of a channel's point in the cube `[-1, 1]³`: how much of the `X` Bloch
direction survives. -/
def lx (p : Pauli R) : R := p.i + p.x - p.y - p.z

/-- The `Y` coordinate of a channel's point in the cube `[-1, 1]³`. -/
def ly (p : Pauli R) : R := p.i - p.x + p.y - p.z

/-- The `Z` coordinate of a channel's point in the cube `[-1, 1]³`. -/
def lz (p : Pauli R) : R := p.i - p.x - p.y + p.z

/-- Doing one channel then another: `X * Y = Z` (up to phase) and so on. -/
def comp (p q : Pauli R) : Pauli R where
  i := p.i * q.i + p.x * q.x + p.y * q.y + p.z * q.z
  x := p.i * q.x + p.x * q.i + p.y * q.z + p.z * q.y
  y := p.i * q.y + p.y * q.i + p.x * q.z + p.z * q.x
  z := p.i * q.z + p.z * q.i + p.x * q.y + p.y * q.x

end Pauli

section CommRing

variable [CommRing R]

theorem comp_lx (p q : Pauli R) : (p.comp q).lx = p.lx * q.lx := by
  simp only [Pauli.comp, Pauli.lx]; ring

theorem comp_ly (p q : Pauli R) : (p.comp q).ly = p.ly * q.ly := by
  simp only [Pauli.comp, Pauli.ly]; ring

theorem comp_lz (p q : Pauli R) : (p.comp q).lz = p.lz * q.lz := by
  simp only [Pauli.comp, Pauli.lz]; ring

/-- The total chance of a composed channel is the product of the totals. -/
theorem comp_total (p q : Pauli R) :
    (p.comp q).i + (p.comp q).x + (p.comp q).y + (p.comp q).z =
      (p.i + p.x + p.y + p.z) * (q.i + q.x + q.y + q.z) := by
  simp only [Pauli.comp]; ring

/-- The four faces of the tetrahedron are the four chances (times `4`). So "inside the
tetrahedron" is exactly "all chances are nonnegative". -/
theorem faces (p : Pauli R) (h : p.i + p.x + p.y + p.z = 1) :
    1 + p.lx + p.ly + p.lz = 4 * p.i ∧ 1 + p.lx - p.ly - p.lz = 4 * p.x ∧
    1 - p.lx + p.ly - p.lz = 4 * p.y ∧ 1 - p.lx - p.ly + p.lz = 4 * p.z := by
  simp only [Pauli.lx, Pauli.ly, Pauli.lz]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> linear_combination -h

/-! ### Codes: the 3-qubit repetition code -/

/-- Majority vote of `3` fails when `2` or `3` copies flip (the top two layers of the
`3`-cube). -/
def maj3 (p : R) : R := 3 * p ^ 2 * (1 - p) + p ^ 3

/-- Phase-flips go straight through the bit-flip code: an odd number of them is a logical
phase-flip, so the polarization is cubed. -/
theorem odd3_polar (q : R) : 1 - 2 * (3 * q * (1 - q) ^ 2 + q ^ 3) = (1 - 2 * q) ^ 3 := by ring

end CommRing

section Field

variable [Field K]

/-- The majority vote of `3` in the polarization coordinate `l = 1 - 2 * p`. -/
def vote (l : K) : K := (3 * l - l ^ 3) / 2

/-- The vote only has fixed points at `-1`, `0` and `1` (the mirror, the wall, the clean
corner). -/
theorem vote_sub_self [NeZero (2 : K)] (l : K) : vote l - l = l * (1 - l) * (1 + l) / 2 := by
  unfold vote; field_simp; ring

theorem maj3_polar [NeZero (2 : K)] (p : K) : 1 - 2 * maj3 p = vote (1 - 2 * p) := by
  unfold maj3 vote; field_simp; ring

/-- Between the wall and the clean corner the vote always helps. -/
theorem vote_helps [LinearOrder K] [IsStrictOrderedRing K] (l : K) (h0 : 0 < l) (h1 : l < 1) :
    l < vote l := by
  have a : 0 < 1 - l := sub_pos.2 h1
  rw [← sub_pos, vote_sub_self]
  positivity

/-! ### The Shor code's thresholds, bracketed exactly

Two levels of the concatenated Shor code send the bit-flip polarization `x` to
`vote (vote x ^ 9)` and the phase-flip polarization `z` to `vote (vote (z ^ 3)) ^ 3`. A point
that moves up is inside the safe box; one that moves down is outside. In flip chances: bit-flips
cross between 9.65% and 9.70%, phase-flips between 7.45% and 7.50%. -/

/-- One level of the Shor code on the bit-flip polarization. -/
def gx (x : K) : K := vote (vote x ^ 9)

/-- One level of the Shor code on the phase-flip polarization. -/
def gz (z : K) : K := vote (vote (z ^ 3)) ^ 3

theorem gx_below : gx (806 / 1000 : ℚ) < 806 / 1000 := by norm_num [gx, vote]

theorem gx_above : gx (807 / 1000 : ℚ) > 807 / 1000 := by norm_num [gx, vote]

theorem gz_below : gz (850 / 1000 : ℚ) < 850 / 1000 := by norm_num [gz, vote]

theorem gz_above : gz (851 / 1000 : ℚ) > 851 / 1000 := by norm_num [gz, vote]

/-! ### MacWilliams: switching to the dual code is flipping the cube's signs -/

/-- The MacWilliams identity for the 3-repetition code. -/
theorem macwilliams_rep3 [NeZero (2 : K)] (x y : K) :
    ((x + y) ^ 3 + (x - y) ^ 3) / 2 = x ^ 3 + 3 * x * y ^ 2 := by
  field_simp; ring

/-- The chance a 3-repetition code misses an error (all three flip), computed two ways: on the
code's own cube, and from the dual code's layers in the polarization `1 - 2 * p`. -/
theorem undetected_two_ways [CharZero K] (p : K) :
    ((1 - p) ^ 3 + p ^ 3) - (1 - p) ^ 3 = (2 / 8) * (1 + 3 * (1 - 2 * p) ^ 2) - (1 - p) ^ 3 := by
  ring

/-! ### The echo: composing with the mult-add inverse -/

/-- Composing with the mult-add inverse undoes the channel. This is
`MultAdd.multAdd_multAddInv`. -/
theorem madd_inverse (a : K) (h : 1 + a ≠ 0) : a ✳ (-a / (1 + a)) = 0 :=
  multAdd_multAddInv a h

end Field

end QuantumNoise
