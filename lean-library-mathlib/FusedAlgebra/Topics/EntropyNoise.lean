module

public import FusedAlgebra.Topics.QuantumNoise
public import FusedAlgebra.Topics.InformationTheory
public import FusedAlgebra.Topics.Finance

/-!
# Entropy and quantum noise on the mult-add cube

This file joins the information-theory kit (`InfoFused`) to the Pauli channels of
`QuantumNoise`.

## Main results

* The entropy weight. `certainty_absorbs_iff`: certainty absorbs under `wAdd k` exactly when
  `k = -1`, so Tsallis' line `k = 1 - q` meets the "or" line only at `q = 2`.
  `escort_escort`: escorts of escorts multiply their orders.
* The Pauli tetrahedron. `collision_parseval`: the collision sum `∑ p ^ 2` is
  `(1 + |λ| ^ 2) / 4`, so the collision (Rényi-2) entropy has spheres for level sets;
  `midsphere_iff`: `∑ p ^ 2 ≤ 1 / 2` exactly inside the midsphere `|λ| = 1`.
  `octahedron_in_ball`, `ball_in_cube`: the octahedron `|λ|₁ ≤ 1` sits inside that ball, which
  sits inside the cube. `compose_stays_in_ball`, `compose_stays_in_octahedron`,
  `compose_stays_in_cube`: composing with any channel keeps a channel inside each of the three.
  `collision_series`, `collision_parallel`: series and parallel composition of collision sums.
* Stacking. `pattern_linear_entropy`: the `n`-fold fault pattern's linear entropy is the `k = -1`
  climb. `cancellation_hides`: the net channel is purer than its fault pattern,
  `((1 + x) / 2) ^ n ≤ (1 + x ^ n) / 2`.
* The 3-bit code. `synd_diagonal`: syndromes are the body diagonals of the 3-cube;
  `alarm_one_coin`: an alarm leaves exactly one fresh coin; `chain_rule_pair`: the entropy chain
  rule on one diagonal, for any `L` with `L (x * y) = L x + L y`.
* Nishimori. `nishimori_tanh`, `kw_lens_polar`: on the Nishimori line `tanh K = 1 - 2 * p`, and
  the Kramers–Wannier lens sends it to the error odds `p / (1 - p)`.

The inequalities are elementary (`cancellation_hides` is the two-point power mean inequality,
`FinanceMultAdd.arith_mean_pow_le`). Nothing here uses real numbers or logarithms: Shannon and
von Neumann entropies as numbers, capacities and thresholds are not part of this file. The
core-Lean original worked over `ℚ`; here the statements hold in any commutative ring, field of
characteristic zero or linearly ordered field, as each one needs.
-/

@[expose] public section

namespace EntropyNoise

open MultAdd QuantumNoise InfoFused

variable {R : Type*} {K : Type*}

/-! ### The entropy weight -/

/-- Certainty swallows everything, `wAdd k a 1 = 1` for every `a`, exactly when `k = -1`. The
Tsallis laws (`k = 1 - q`) meet the "or" laws only there, at `q = 2`. -/
theorem certainty_absorbs_iff [CommRing R] (k : R) : (∀ a : R, wAdd k a 1 = 1) ↔ k = -1 := by
  constructor
  · intro h
    have h1 := h 1
    unfold wAdd at h1
    linear_combination h1
  · rintro rfl a
    unfold wAdd; ring

/-- Raising every entry to `a` then to `b` is raising to `a * b`: escorts of escorts. -/
theorem escort_powers [Monoid R] (l : List R) (a b : ℕ) :
    (l.map (· ^ a)).map (· ^ b) = l.map (· ^ (a * b)) := by
  simp [pow_mul]

/-- The escort of order `n` of a list of weights: `xᵢ ^ n / ∑ xⱼ ^ n`. -/
def escort [DivisionRing K] (n : ℕ) (l : List K) : List K :=
  l.map fun x ↦ x ^ n / (l.map (· ^ n)).sum

/-- Escort of an escort: `P_b (P_a p) = P_{a * b} p`, when the normalizer is nonzero. -/
theorem escort_escort [Field K] (l : List K) (a b : ℕ) (ha : (l.map (· ^ a)).sum ≠ 0) :
    escort b (escort a l) = escort (a * b) l := by
  unfold escort
  generalize (l.map (· ^ a)).sum = Z at ha ⊢
  have hpow : ∀ x : K, (x ^ a / Z) ^ b = x ^ (a * b) / Z ^ b := fun x ↦ by
    rw [div_pow, pow_mul]
  have hsum : ((l.map fun x ↦ x ^ a / Z).map (· ^ b)).sum =
      (l.map (· ^ (a * b))).sum / Z ^ b := by
    simp only [List.map_map, Function.comp_def, hpow]
    simp only [div_eq_mul_inv, List.sum_map_mul_right]
  rw [hsum, List.map_map]
  refine List.map_congr_left fun x _ ↦ ?_
  simp only [Function.comp_apply, hpow]
  exact div_div_div_cancel_right₀ (pow_ne_zero b ha) _ _

/-! ### The Pauli tetrahedron: three norms of `λ` -/

section Tetrahedron

/-- The squared length of a channel's point `λ` in the cube. -/
def len2 [CommRing R] (p : Pauli R) : R := p.lx ^ 2 + p.ly ^ 2 + p.lz ^ 2

/-- The collision sum `∑ p ^ 2`, whose `-log` is the collision (Rényi-2) entropy. -/
def coll [Semiring R] (p : Pauli R) : R := p.i ^ 2 + p.x ^ 2 + p.y ^ 2 + p.z ^ 2

/-- The clean corner of the tetrahedron: the identity channel. -/
def corner [Zero R] [One R] : Pauli R := ⟨1, 0, 0, 0⟩

theorem corner_point [CommRing R] :
    (corner : Pauli R).lx = 1 ∧ (corner : Pauli R).ly = 1 ∧ (corner : Pauli R).lz = 1 ∧
      coll (corner : Pauli R) = 1 := by
  simp [corner, Pauli.lx, Pauli.ly, Pauli.lz, coll]

section Field

variable [Field K] [CharZero K]

/-- Parseval on the tetrahedron: `∑ p ^ 2 = (1 + |λ| ^ 2) / 4`. -/
theorem collision_parseval (p : Pauli K) (h : p.i + p.x + p.y + p.z = 1) :
    coll p = (1 + len2 p) / 4 := by
  rw [eq_div_iff (by norm_num)]
  simp only [coll, len2, Pauli.lx, Pauli.ly, Pauli.lz]
  linear_combination (p.i + p.x + p.y + p.z + 1) * h

/-- The Choi state's linear entropy `1 - ∑ p ^ 2` is `(3 - |λ| ^ 2) / 4`: zero at the corners,
largest (`3 / 4`) at the centre, constant on spheres about the centre. -/
theorem linear_entropy_choi (p : Pauli K) (h : p.i + p.x + p.y + p.z = 1) :
    1 - coll p = (3 - len2 p) / 4 := by
  rw [collision_parseval p h]; ring

/-- The centre of the tetrahedron: the completely depolarizing channel. -/
def centre : Pauli K := ⟨1 / 4, 1 / 4, 1 / 4, 1 / 4⟩

/-- The midpoint of an edge of the tetrahedron: `X` with chance `1 / 2`. -/
def edgeMid : Pauli K := ⟨1 / 2, 1 / 2, 0, 0⟩

theorem centre_point :
    (centre : Pauli K).lx = 0 ∧ (centre : Pauli K).ly = 0 ∧ (centre : Pauli K).lz = 0 ∧
      coll (centre : Pauli K) = 1 / 4 := by
  simp only [centre, Pauli.lx, Pauli.ly, Pauli.lz, coll]; norm_num

/-- The edge midpoint is `(1, 0, 0)`: on the cube's face, the sphere and the octahedron at once,
with collision sum exactly `1 / 2`. -/
theorem edge_midpoint :
    (edgeMid : Pauli K).lx = 1 ∧ (edgeMid : Pauli K).ly = 0 ∧ (edgeMid : Pauli K).lz = 0 ∧
      coll (edgeMid : Pauli K) = 1 / 2 ∧ len2 (edgeMid : Pauli K) = 1 := by
  simp only [edgeMid, Pauli.lx, Pauli.ly, Pauli.lz, coll, len2]; norm_num

/-- Series: the collision sum of a composed channel, through Parseval. -/
theorem collision_series (p q : Pauli K) (hp : p.i + p.x + p.y + p.z = 1)
    (hq : q.i + q.x + q.y + q.z = 1) :
    coll (p.comp q) =
      (1 + (p.lx * q.lx) ^ 2 + (p.ly * q.ly) ^ 2 + (p.lz * q.lz) ^ 2) / 4 := by
  have ht : (p.comp q).i + (p.comp q).x + (p.comp q).y + (p.comp q).z = 1 := by
    rw [comp_total, hp, hq, one_mul]
  rw [collision_parseval _ ht, len2, comp_lx, comp_ly, comp_lz]; ring

/-- The depolarizing diagonal: total fault chance `s` spread evenly over `X`, `Y`, `Z`. -/
def depol (s : K) : Pauli K := ⟨1 - s, s / 3, s / 3, s / 3⟩

theorem depol_point (s : K) :
    (depol s).lx = 1 - 4 / 3 * s ∧ (depol s).ly = 1 - 4 / 3 * s ∧
      (depol s).lz = 1 - 4 / 3 * s := by
  simp only [depol, Pauli.lx, Pauli.ly, Pauli.lz]
  refine ⟨?_, ?_, ?_⟩ <;> ring

end Field

/-- Squaring the cells of the outer-product grid is the grid of the squares. -/
theorem map_sq_outer [CommSemiring R] (xs ys : List R) :
    (outer xs ys).map (· ^ 2) = outer (xs.map (· ^ 2)) (ys.map (· ^ 2)) := by
  induction xs with
  | nil => simp
  | cons x xs ih => simp [outer_cons, ih, mul_pow, Function.comp_def]

/-- Parallel (two qubits, independent faults): the collision sum of the 16 joint fault chances
is the product, so linear entropy composes by the union rule (`InfoFused.gini_simpson`). -/
theorem collision_parallel [CommSemiring R] (p q : Pauli R) :
    ((outer [p.i, p.x, p.y, p.z] [q.i, q.x, q.y, q.z]).map (· ^ 2)).sum = coll p * coll q := by
  rw [map_sq_outer, outer_sum]
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, coll]; ring

section Order

variable [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]

/-- `|a * b| ≤ u` whenever `|a| ≤ u` and `|b| ≤ 1`. -/
theorem abs_mul_le_of_abs_le_one {a b u : R} (ha : |a| ≤ u) (hb : |b| ≤ 1) : |a * b| ≤ u := by
  rw [abs_mul]
  exact (mul_le_of_le_one_right (abs_nonneg a) hb).trans ha

/-- Octahedron inside the ball: `|x| + |y| + |z| ≤ 1` gives `x ^ 2 + y ^ 2 + z ^ 2 ≤ 1`. -/
theorem octahedron_in_ball (x y z : R) (h : |x| + |y| + |z| ≤ 1) :
    x ^ 2 + y ^ 2 + z ^ 2 ≤ 1 := by
  have hx := abs_nonneg x
  have hy := abs_nonneg y
  have hz := abs_nonneg z
  rw [← sq_abs x, ← sq_abs y, ← sq_abs z]
  nlinarith

/-- Ball inside the cube. -/
theorem ball_in_cube (x y z : R) (h : x ^ 2 + y ^ 2 + z ^ 2 ≤ 1) : x ^ 2 ≤ 1 := by
  linarith [sq_nonneg y, sq_nonneg z]

/-- Composition multiplies coordinates (`QuantumNoise.comp_lx` and its companions), and a
coordinate of a channel lies in `[-1, 1]`. So composing never leaves the ball. -/
theorem compose_stays_in_ball (ax ay az bx b_y bz : R) (hx : bx ^ 2 ≤ 1) (hy : b_y ^ 2 ≤ 1)
    (hz : bz ^ 2 ≤ 1) :
    (ax * bx) ^ 2 + (ay * b_y) ^ 2 + (az * bz) ^ 2 ≤ ax ^ 2 + ay ^ 2 + az ^ 2 := by
  simp only [mul_pow]
  linarith [mul_le_of_le_one_right (sq_nonneg ax) hx, mul_le_of_le_one_right (sq_nonneg ay) hy,
    mul_le_of_le_one_right (sq_nonneg az) hz]

/-- Composing never leaves the octahedron. -/
theorem compose_stays_in_octahedron {ax ay az bx b_y bz u v w : R} (hx : |ax| ≤ u)
    (hy : |ay| ≤ v) (hz : |az| ≤ w) (hbx : |bx| ≤ 1) (hby : |b_y| ≤ 1) (hbz : |bz| ≤ 1) :
    |ax * bx| ≤ u ∧ |ay * b_y| ≤ v ∧ |az * bz| ≤ w :=
  ⟨abs_mul_le_of_abs_le_one hx hbx, abs_mul_le_of_abs_le_one hy hby,
    abs_mul_le_of_abs_le_one hz hbz⟩

/-- Composing never leaves any cube `max |λᵢ| ≤ c` about the centre (classical capacity cannot
rise). -/
theorem compose_stays_in_cube {a b c : R} (ha : |a| ≤ c) (hb : |b| ≤ 1) : |a * b| ≤ c :=
  abs_mul_le_of_abs_le_one ha hb

end Order

section OrderedField

variable [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- `∑ p ^ 2 ≤ 1 / 2` (collision entropy at least one bit) exactly inside the midsphere. -/
theorem midsphere_iff (p : Pauli K) (h : p.i + p.x + p.y + p.z = 1) :
    coll p ≤ 1 / 2 ↔ len2 p ≤ 1 := by
  rw [collision_parseval p h]; constructor <;> intro h' <;> linarith

/-- Along the diagonal `(λ, λ, λ)` the octahedron's face is `λ = 1 / 3` and the midsphere is
`λ ^ 2 = 1 / 3`, so the octahedron is strictly inside there. -/
theorem diagonal_crossings : (1 / 3 : K) + 1 / 3 + 1 / 3 = 1 ∧ 3 * (1 / 3 : K) ^ 2 < 1 := by
  norm_num

end OrderedField

end Tetrahedron

/-! ### Stacking: what cancellation hides -/

section Stacking

/-- One bit-flip channel's collision sum in the polarization `λ = 1 - 2 * p`. -/
theorem bit_collision [Field K] [CharZero K] (p : K) :
    (1 - p) ^ 2 + p ^ 2 = (1 + (1 - 2 * p) ^ 2) / 2 := by
  ring

/-- The `n`-fold fault pattern (a corner of the `n`-cube) has collision sum `c ^ n`, so its
linear entropy climbs at weight `-1`. -/
theorem pattern_linear_entropy [CommRing R] (c : R) (n : ℕ) :
    1 - c ^ n = wpow (-1) (1 - c) n := by
  rw [wpow_neg_one, sub_sub_cancel]

/-- The net channel after `n` stacked flips is purer than its fault pattern:
`((1 + x) / 2) ^ n ≤ (1 + x ^ n) / 2` for `x = λ ^ 2`. This is the two-point power mean
inequality `FinanceMultAdd.arith_mean_pow_le`. -/
theorem cancellation_hides [Field K] [LinearOrder K] [IsStrictOrderedRing K] (x : K)
    (h0 : 0 ≤ x) (n : ℕ) : ((1 + x) / 2) ^ n ≤ (1 + x ^ n) / 2 := by
  simpa using FinanceMultAdd.arith_mean_pow_le zero_le_one h0 n

/-- At `n = 2` the purity gap is a perfect square. -/
theorem gap_two [Field K] [CharZero K] (l : K) :
    (1 + l ^ 4) / 2 - ((1 + l ^ 2) / 2) ^ 2 = (1 - l ^ 2) ^ 2 / 4 := by
  ring

end Stacking

/-! ### The 3-bit code: syndromes are body diagonals -/

/-- The two parity checks of the 3-bit code on a corner `(a, b, c)` of the cube. -/
def synd (a b c : Bool) : Bool × Bool := (a != b, b != c)

/-- Opposite corners give the same syndrome. -/
theorem synd_antipodal : ∀ a b c : Bool, synd (!a) (!b) (!c) = synd a b c := by decide

/-- Corners with the same syndrome are equal or opposite: each syndrome is one body diagonal of
the 3-cube. -/
theorem synd_diagonal : ∀ a b c d e f : Bool, synd a b c = synd d e f →
    (a = d ∧ b = e ∧ c = f) ∨ (a = !d ∧ b = !e ∧ c = !f) := by decide

section Code

variable [CommRing R]

/-- No alarm: both ends of the main diagonal. Alarm: each of the three other diagonals. -/
theorem syndromes_total (p : R) :
    ((1 - p) ^ 3 + p ^ 3) + 3 * (p * (1 - p) ^ 2 + p ^ 2 * (1 - p)) = 1 := by ring

/-- The vote picks the lighter end of each diagonal; the failures add up to the top two layers of
the cube (`QuantumNoise.maj3`). -/
theorem decoder_failure (p : R) : p ^ 3 + 3 * (p ^ 2 * (1 - p)) = maj3 p := by
  unfold maj3; ring

end Code

section CodeField

variable [Field K]

theorem no_alarm_polar [CharZero K] (p : K) :
    (1 - p) ^ 3 + p ^ 3 = (1 + 3 * (1 - 2 * p) ^ 2) / 4 := by ring

theorem alarm_polar [CharZero K] (p : K) :
    p * (1 - p) ^ 2 + p ^ 2 * (1 - p) = (1 - (1 - 2 * p) ^ 2) / 4 := by ring

/-- After an alarm, the leftover doubt is exactly one fresh coin of chance `p`. -/
theorem alarm_one_coin (p q : K) (hpq : p + q = 1) (hp : p ≠ 0) (hq : q ≠ 0) :
    p * q * q / (p * q * q + p * p * q) = q := by
  have e : p * q * q + p * p * q = p * q := by linear_combination p * q * hpq
  rw [e, mul_div_cancel_left₀ q (mul_ne_zero hp hq)]

/-- The chain rule on one diagonal, for any `L` with `L (x * y) = L x + L y` (`log` is one):
`a * L a + b * L b = s * L s + (a * L (a / s) + b * L (b / s))` with `s = a + b`. Entropy of the
faults is entropy of the syndrome plus entropy left on the diagonal. -/
theorem chain_rule_pair (L : K → K) (hL : ∀ x y, L (x * y) = L x + L y) (a b : K)
    (hs : a + b ≠ 0) :
    a * L a + b * L b =
      (a + b) * L (a + b) + (a * L (a / (a + b)) + b * L (b / (a + b))) := by
  have la := hL (a + b) (a / (a + b))
  have lb := hL (a + b) (b / (a + b))
  rw [mul_div_cancel₀ _ hs] at la lb
  linear_combination a * la + b * lb

end CodeField

/-! ### Nishimori: the polarization is `tanh` of the coupling -/

section Nishimori

variable [Field K]

/-- On the Nishimori line `exp (-2 * K) = w = p / (1 - p)`. Then
`tanh K = (1 - w) / (1 + w) = 1 - 2 * p`. -/
theorem nishimori_tanh (p : K) (hp : 1 - p ≠ 0) :
    (1 - p / (1 - p)) / (1 + p / (1 - p)) = 1 - 2 * p := by
  have e₁ : 1 - p / (1 - p) = (1 - 2 * p) / (1 - p) := by field_simp; ring
  have e₂ : 1 + p / (1 - p) = 1 / (1 - p) := by field_simp; ring
  rw [e₁, e₂]; field_simp

/-- The Kramers–Wannier lens `M t = (1 - t) / (1 + t)` sends `λ = 1 - 2 * p` to
`p / (1 - p)`. -/
theorem kw_lens_polar [NeZero (2 : K)] (p : K) :
    (1 - (1 - 2 * p)) / (1 + (1 - 2 * p)) = p / (1 - p) := by
  rw [show 1 - (1 - 2 * p) = 2 * p by ring, show 1 + (1 - 2 * p) = 2 * (1 - p) by ring,
    mul_div_mul_left _ _ two_ne_zero]

end Nishimori

end EntropyNoise
