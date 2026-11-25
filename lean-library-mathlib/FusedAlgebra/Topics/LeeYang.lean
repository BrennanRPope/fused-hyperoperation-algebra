module

public import FusedAlgebra.Basic

/-!
# Lee–Yang zeros and the dial

Two ways to build a partition function from the dial `k` (with `q = k + 1`).

* **The dial gas.** One edge, `1 + wAdd k x y = 1 + x + y + k * x * y`, is the two-site lattice
  gas with pair weight `k` (`k = 1` free, `k = 0` hard core). On a ring of `n` sites the transfer
  matrix is `ringT k z = !![1, z; 1, k * z]`, whose trace is the weighted lift `1 + k * z`. When
  the trace is `0` (`z` on the wall `-1 / k`), `T ^ 2 = -det T • 1`, so every odd power of `T`
  has trace `0`: every odd ring vanishes on the wall, for every `k`.
* **The Potts ring.** Its partition function is
  `(Q + v) ^ n + (Q - 1) * v ^ n = (Q + v) ^ n * (1 + k * t ^ n)` with `t = v / (Q + v)` and
  `k = Q - 1`: the weighted lift of `t ^ n`. Proper colourings are `v = -1`. Rings vanish at
  `Q = 0` and `Q = 1`, at `Q = 2` exactly when odd, and never at `Q = -1`.

## Main results

* `sq_eq_trace_smul_sub_det_smul_one`: Cayley–Hamilton for `2 × 2` matrices.
* `pow_two_mul_add_one_of_trace_eq_zero`, `pow_two_mul_of_trace_eq_zero`: powers of a traceless
  `2 × 2` matrix.
* `odd_ring_vanishes_on_wall`, `even_ring_on_wall`: the wall theorem for rings.
* `triangle_count`, `square_count`, `pentagon_count`: the transfer matrix counts correctly.
* `ising_ring_transfer`, `pottsRing_lift`, `cycle_chromatic`, `parity_never_zero`.

The mathematics is known: the transfer-matrix method, the Cayley–Hamilton theorem (in mathlib for
every size as `Matrix.aeval_self_charpoly`), the chromatic polynomial of a cycle. The core-Lean
original used its own record `Mat2`; here matrices are mathlib's `Matrix (Fin 2) (Fin 2) R`.
-/

@[expose] public section

namespace LeeYang

open MultAdd Matrix

variable {R : Type*} [CommRing R]

/-! ### A. The dial gas -/

/-- One edge of the dial gas: the two-site lattice gas with pair weight `k` is the lift of the
weighted mult-add of the two fugacities. -/
theorem edge_gas (k x y : R) : 1 + wAdd k x y = 1 + x + y + k * x * y := by
  unfold wAdd; ring

/-- At `k = 1` (`q = 2`) the two sites are independent: the gas factors. -/
theorem edge_gas_free (x y : R) : 1 + wAdd 1 x y = (1 + x) * (1 + y) := by
  rw [wAdd_one, one_add_multAdd]

/-- At `k = 0` (`q = 1`) the two sites exclude each other: hard core. -/
theorem edge_gas_hardcore (x y : R) : 1 + wAdd 0 x y = 1 + x + y := by
  rw [wAdd_weight_zero, add_assoc]

/-- The ring transfer matrix: from an empty site the next site weighs `1` (empty) or `z`
(occupied); from an occupied site, `1` or `k * z`. -/
def ringT (k z : R) : Matrix (Fin 2) (Fin 2) R := !![1, z; 1, k * z]

/-- Its trace is the weighted lift of `z`. -/
theorem trace_ringT (k z : R) : (ringT k z).trace = 1 + k * z := by
  simp [ringT, trace_fin_two]

theorem det_ringT (k z : R) : (ringT k z).det = z * (k - 1) := by
  simp [ringT, det_fin_two]; ring

/-- Cayley–Hamilton for `2 × 2` matrices: `x ^ 2 = tr x • x - det x • 1`. The general theorem is
`Matrix.aeval_self_charpoly`. -/
theorem sq_eq_trace_smul_sub_det_smul_one (x : Matrix (Fin 2) (Fin 2) R) :
    x ^ 2 = x.trace • x - x.det • (1 : Matrix (Fin 2) (Fin 2) R) := by
  rw [eta_fin_two x]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pow_two, trace_fin_two_of, det_fin_two_of] <;> ring

/-- A traceless `2 × 2` matrix squares to a scalar. -/
theorem sq_of_trace_eq_zero (x : Matrix (Fin 2) (Fin 2) R) (h : x.trace = 0) :
    x ^ 2 = (-x.det) • (1 : Matrix (Fin 2) (Fin 2) R) := by
  rw [sq_eq_trace_smul_sub_det_smul_one, h, zero_smul, zero_sub, neg_smul]

/-- Even powers of a traceless `2 × 2` matrix are scalar: `x ^ (2 * m) = (-det x) ^ m • 1`. -/
theorem pow_two_mul_of_trace_eq_zero (x : Matrix (Fin 2) (Fin 2) R) (h : x.trace = 0) (m : ℕ) :
    x ^ (2 * m) = (-x.det) ^ m • (1 : Matrix (Fin 2) (Fin 2) R) := by
  rw [pow_mul, sq_of_trace_eq_zero x h, smul_pow, one_pow]

/-- With trace `0`, every odd power is a multiple of `x`:
`x ^ (2 * m + 1) = (-det x) ^ m • x`. -/
theorem pow_two_mul_add_one_of_trace_eq_zero (x : Matrix (Fin 2) (Fin 2) R) (h : x.trace = 0)
    (m : ℕ) : x ^ (2 * m + 1) = (-x.det) ^ m • x := by
  rw [pow_succ, pow_two_mul_of_trace_eq_zero x h, smul_mul_assoc, one_mul]

/-- Odd powers of a traceless `2 × 2` matrix are traceless. -/
theorem trace_pow_two_mul_add_one_of_trace_eq_zero (x : Matrix (Fin 2) (Fin 2) R)
    (h : x.trace = 0) (m : ℕ) : (x ^ (2 * m + 1)).trace = 0 := by
  rw [pow_two_mul_add_one_of_trace_eq_zero x h, trace_smul, h, smul_zero]

/-- **The wall theorem for rings.** The dial gas on a ring of `2 * m + 1` sites, `tr (T ^ n)`,
vanishes whenever `1 + k * z = 0`, that is at `z = -1 / k`, for every weight `k`. -/
theorem odd_ring_vanishes_on_wall (k z : R) (hw : 1 + k * z = 0) (m : ℕ) :
    (ringT k z ^ (2 * m + 1)).trace = 0 :=
  trace_pow_two_mul_add_one_of_trace_eq_zero _ (by rw [trace_ringT]; exact hw) m

/-- On the wall an even ring is `2 * (z * (1 - k)) ^ m`, not zero (unless `k = 1`). -/
theorem even_ring_on_wall (k z : R) (hw : 1 + k * z = 0) (m : ℕ) :
    (ringT k z ^ (2 * m)).trace = 2 * (z * (1 - k)) ^ m := by
  rw [pow_two_mul_of_trace_eq_zero _ (by rw [trace_ringT]; exact hw) m, det_ringT, trace_smul,
    trace_one, Fintype.card_fin, smul_eq_mul, show -(z * (k - 1)) = z * (1 - k) by ring]
  push_cast; ring

/-- The transfer matrix counts correctly: the 3-ring (triangle) is
`1 + 3 * z + 3 * k * z ^ 2 + k ^ 3 * z ^ 3` (each pair of occupied sites shares an edge). -/
theorem triangle_count (k z : R) :
    (ringT k z ^ 3).trace = 1 + 3 * z + 3 * k * z ^ 2 + k ^ 3 * z ^ 3 := by
  simp [ringT, pow_succ, trace_fin_two]; ring

/-- The 4-ring: `1 + 4 * z + (2 + 4 * k) * z ^ 2 + 4 * k ^ 2 * z ^ 3 + k ^ 4 * z ^ 4`. -/
theorem square_count (k z : R) :
    (ringT k z ^ 4).trace
      = 1 + 4 * z + (2 + 4 * k) * z ^ 2 + 4 * k ^ 2 * z ^ 3 + k ^ 4 * z ^ 4 := by
  simp [ringT, pow_succ, trace_fin_two]; ring

/-- The 5-ring:
`1 + 5 * z + (5 + 5 * k) * z ^ 2 + (5 * k + 5 * k ^ 2) * z ^ 3 + 5 * k ^ 3 * z ^ 4 + k ^ 5 * z ^ 5`.
-/
theorem pentagon_count (k z : R) :
    (ringT k z ^ 5).trace =
      1 + 5 * z + (5 + 5 * k) * z ^ 2 + (5 * k + 5 * k ^ 2) * z ^ 3 + 5 * k ^ 3 * z ^ 4
        + k ^ 5 * z ^ 5 := by
  simp [ringT, pow_succ, trace_fin_two]; ring

/-- The triangle splits off the wall: `1 + k * z` times a quadratic. -/
theorem triangle_factor (k z : R) :
    1 + 3 * z + 3 * k * z ^ 2 + k ^ 3 * z ^ 3
      = (1 + k * z) * (1 + (3 - k) * z + k ^ 2 * z ^ 2) := by
  ring

/-- The quadratic's discriminant: `(3 - k) ^ 2 - 4 * k ^ 2 = 3 * (1 - k) * (3 + k)`. It is `≤ 0`
for `k ≥ 1` (the two other zeros are a conjugate pair with product `1 / k ^ 2`, on the circle
`|z| = 1 / k`: Lee–Yang), and `> 0` for `-3 < k < 1` (real zeros). In mathlib's notation the left
side is `discrim (k ^ 2) (3 - k) 1`. -/
theorem triangle_disc (k : R) : (3 - k) ^ 2 - 4 * k ^ 2 = 3 * (1 - k) * (3 + k) := by
  ring

theorem discrim_triangle (k : R) : discrim (k ^ 2) (3 - k) 1 = 3 * (1 - k) * (3 + k) := by
  rw [discrim]; ring

theorem triangle_disc_nonpos {S : Type*} [CommRing S] [LinearOrder S] [IsStrictOrderedRing S]
    (k : S) (hk : 1 ≤ k) : (3 - k) ^ 2 - 4 * k ^ 2 ≤ 0 := by
  rw [triangle_disc]
  nlinarith [mul_nonneg (sub_nonneg.2 hk) (add_nonneg (by norm_num : (0 : S) ≤ 3)
    (zero_le_one.trans hk))]

/-- The ring's eigenvalue discriminant:
`(1 + k * z) ^ 2 - 4 * (z * (k - 1)) = k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1`. Its zeros are the
ends of the zero arc (Yang–Lee edges). At `k = 1` it is `(1 + z) ^ 2`. -/
theorem ring_disc (k z : R) :
    (1 + k * z) ^ 2 - 4 * (z * (k - 1)) = k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 := by
  ring

theorem ring_disc_free (z : R) : (1 + 1 * z) ^ 2 - 4 * (z * (1 - 1)) = (1 + z) ^ 2 := by
  ring

/-- The ring's zero family. A ring of `n` sites vanishes exactly when the two eigenvalues have
ratio `e ^ (i φ)` with `φ = π (2 j + 1) / n`, which is `(1 + k z) ^ 2 = c * z * (k - 1)` with
`c = 4 cos ^ 2 (φ / 2)` in `[0, 4]`. Each `c` gives a quadratic in `z`. -/
theorem ring_zero_family (k z c : R) :
    (1 + k * z) ^ 2 - c * (z * (k - 1)) = k ^ 2 * z ^ 2 + (2 * k - c * (k - 1)) * z + 1 := by
  ring

/-- `c = 0` is the wall, doubled: `(1 + k * z) ^ 2`. -/
theorem ring_family_wall (k z : R) : (1 + k * z) ^ 2 - 0 * (z * (k - 1)) = (1 + k * z) ^ 2 := by
  ring

/-- `c = 4` is the edge of the zero arc (the eigenvalue discriminant). -/
theorem ring_family_edge (k z : R) :
    (1 + k * z) ^ 2 - 4 * (z * (k - 1)) = k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 := ring_disc k z

/-- In `ζ = k * z` the family reads `k * (ζ ^ 2 + 1) + (2 * k - c * (k - 1)) * ζ = 0` after
multiplying by `k`: so `ζ + 1 / ζ = c * (1 - 1 / k) - 2` runs from `-2` (the wall, `c = 0`) to
`2 - 4 / k` (`c = 4`). -/
theorem ring_family_zeta (k z c : R) :
    k * (k ^ 2 * z ^ 2 + (2 * k - c * (k - 1)) * z + 1)
      = k * ((k * z) ^ 2 + 1) + (2 * k - c * (k - 1)) * (k * z) := by
  ring

/-- The two edge points `z₁`, `z₂` of the zero arc satisfy `k ^ 2 * (z₁ * z₂) = 1` (Vieta), so in
`ζ = k * z` they are a reciprocal pair, symmetric under `ζ ↦ 1 / ζ` (the up–down flip of the
spins). -/
theorem ring_edges_vieta (k z z1 z2 : R) (hp : k ^ 2 * (z1 * z2) = 1)
    (hs : k ^ 2 * (z1 + z2) = 2 * k - 4) :
    k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 = k ^ 2 * (z - z1) * (z - z2) := by
  linear_combination -hp + z * hs

/-! ### B. The Potts ring: the dial is `k = Q - 1` -/

/-- Partition function of the `Q`-state Potts model on a ring of `n` sites with edge weight
`1 + v * δ` (eigenvalues of the edge matrix: `Q + v` once, `v` with multiplicity `Q - 1`). -/
def pottsRing (Q v : R) (n : ℕ) : R := (Q + v) ^ n + (Q - 1) * v ^ n

/-- For `Q = 2` (Ising) the formula is the trace of the transfer matrix
`!![1 + v, 1; 1, 1 + v]`. -/
theorem ising_ring_transfer (v : R) (n : ℕ) :
    ((!![1 + v, 1; 1, 1 + v] : Matrix (Fin 2) (Fin 2) R) ^ n).trace = pottsRing 2 v n := by
  have key : ∀ n : ℕ, ∃ a b : R,
      (!![1 + v, 1; 1, 1 + v] : Matrix (Fin 2) (Fin 2) R) ^ n = !![a, b; b, a] ∧
        a + b = (2 + v) ^ n ∧ a - b = v ^ n := by
    intro n
    induction n with
    | zero => exact ⟨1, 0, by simp [one_fin_two], by simp, by simp⟩
    | succ n ih =>
      obtain ⟨a, b, h, h1, h2⟩ := ih
      refine ⟨a * (1 + v) + b, a + b * (1 + v), ?_, ?_, ?_⟩
      · rw [pow_succ, h, mul_fin_two]
        simp [add_comm]
      · rw [pow_succ, ← h1]; ring
      · rw [pow_succ, ← h2]; ring
  obtain ⟨a, b, h, h1, h2⟩ := key n
  rw [h, pottsRing, ← h1, ← h2]
  simp [trace_fin_two]
  ring

/-- **The Potts ring is the weighted lift of `t ^ n`.** With `t = v / (Q + v)` and `k = Q - 1`:
`Z = (Q + v) ^ n * (1 + k * t ^ n)`. -/
theorem pottsRing_lift (Q v t : R) (n : ℕ) (ht : t * (Q + v) = v) :
    pottsRing Q v n = (Q + v) ^ n * (1 + (Q - 1) * t ^ n) := by
  rw [pottsRing]
  nth_rw 2 [← ht]
  rw [mul_pow]; ring

/-- Two separate rings together: their weights multiply, which in the lift is the weighted
mult-add of the two `t ^ n`, at the dial `k = Q - 1`. -/
theorem two_rings (Q t A B : R) (n m : ℕ) :
    (A * (1 + (Q - 1) * t ^ n)) * (B * (1 + (Q - 1) * t ^ m))
      = (A * B) * (1 + (Q - 1) * wAdd (Q - 1) (t ^ n) (t ^ m)) := by
  rw [one_add_mul_wAdd]; ring

/-- Proper colourings: `v = -1`. The chromatic polynomial of the `n`-cycle. -/
theorem cycle_chromatic (Q : R) (n : ℕ) :
    pottsRing Q (-1) n = (Q - 1) ^ n + (-1) ^ n * (Q - 1) := by
  rw [pottsRing, ← sub_eq_add_neg, mul_comm]

/-- Colourings sit on the wall: `v = -1` gives `t * (Q + v) = v` with `t = -1 / (Q - 1)`, that
is `1 + k * t = 0` for `k = Q - 1`. (Stated without division.) -/
theorem colourings_on_wall (Q t : R) (h : (Q - 1) * t = -1) : t * (Q + (-1)) = -1 := by
  linear_combination h

/-- `q = 0` (the "or" point, `k = -1`) is a zero of every ring. -/
theorem ring_zero_at_q0 (n : ℕ) : pottsRing (0 : R) (-1) n = 0 := by
  rw [pottsRing]; ring

/-- `q = 1` (plain addition, `k = 0`) is a zero of every ring with at least one site. -/
theorem ring_zero_at_q1 (n : ℕ) : pottsRing (1 : R) (-1) (n + 1) = 0 := by
  simp [pottsRing]

/-- `q = 2` (mult-add, `k = 1`) is a zero exactly for odd rings. -/
theorem ring_zero_at_q2_odd (m : ℕ) : pottsRing (2 : R) (-1) (2 * m + 1) = 0 := by
  rw [cycle_chromatic]; simp [pow_succ, pow_mul]; norm_num

theorem ring_at_q2_even (m : ℕ) : pottsRing (2 : R) (-1) (2 * m) = 2 := by
  rw [cycle_chromatic, pow_mul]; norm_num

/-- `q = -1` (parity, `k = -2`): the ring gives `(-1) ^ n * (2 ^ n - 2)`, the signed number of
acyclic orientations of the cycle. -/
theorem ring_at_parity (n : ℕ) : pottsRing (-1 : R) (-1) n = (-1) ^ n * (2 ^ n - 2) := by
  rw [pottsRing, show (-1 : R) + -1 = -1 * 2 by ring, mul_pow]; ring

/-- The parity point is never a zero of a ring with two or more sites (over `ℤ`). -/
theorem parity_never_zero (n : ℕ) (hn : 2 ≤ n) : pottsRing (-1 : ℤ) (-1) n ≠ 0 := by
  rw [ring_at_parity]
  have h4 : (2 : ℤ) ^ 2 ≤ 2 ^ n := pow_le_pow_right₀ (by norm_num) hn
  exact mul_ne_zero (pow_ne_zero _ (by norm_num)) (by omega)

end LeeYang
