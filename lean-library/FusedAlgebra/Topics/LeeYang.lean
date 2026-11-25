/-
  Fused Hyperoperation Algebra: Lee–Yang zeros and the dial.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/LeeYang.lean

  LIBRARY NOTE: copied from lee-yang/lean/LeeYang.lean (Nov 25), unchanged apart from this
  note. No theorem states exactly what a library theorem states; the closest are
  edge_gas_hardcore and edge_gas_free, which add 1 to MultAdd.wAdd_weight_zero and
  MultAdd.wAdd_one.

  Two ways to build a partition function from the dial k (with q = k + 1):

  A. The dial gas. One edge: 1 + (x ⊞ₖ y) = 1 + x + y + k·xy is the two-site lattice gas
     with pair weight k (k = 1 free, k = 0 hard core). On a ring of n sites the transfer
     matrix is T = [[1, z], [1, k·z]], whose trace is the weighted lift 1 + k·z.
     When the trace is 0 (z on the wall −1/k), T² = −det(T)·I, so every odd power of T has
     trace 0: every odd ring vanishes exactly on the wall, for every k.
  B. The Potts ring. Its partition function is (Q + v)ⁿ + (Q − 1)·vⁿ = (Q + v)ⁿ·(1 + k·tⁿ)
     with t = v/(Q + v) and k = Q − 1: the weighted lift of tⁿ. Proper colourings are v = −1,
     t = −1/k, the wall. Rings vanish at Q = 0 and Q = 1, at Q = 2 exactly when odd, and
     never at Q = −1 (the parity point): there the count is ±(2ⁿ − 2).
-/
import FusedAlgebra.Base

namespace LeeYang
open MultAdd

set_option linter.unusedSectionVars false

/-! ## A 2×2 matrix record (self-contained) -/

structure Mat2 (R : Type) where
  a : R
  b : R
  c : R
  d : R

theorem Mat2.ext' {R : Type} {x y : Mat2 R} (h1 : x.a = y.a) (h2 : x.b = y.b)
    (h3 : x.c = y.c) (h4 : x.d = y.d) : x = y := by
  cases x; cases y
  dsimp only at h1 h2 h3 h4
  subst h1; subst h2; subst h3; subst h4; rfl

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

namespace Mat2
def one : Mat2 R := ⟨1, 0, 0, 1⟩
def mul (x y : Mat2 R) : Mat2 R :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def smul (s : R) (x : Mat2 R) : Mat2 R := ⟨s * x.a, s * x.b, s * x.c, s * x.d⟩
def tr (x : Mat2 R) : R := x.a + x.d
def det (x : Mat2 R) : R := x.a * x.d - x.b * x.c
def pow (x : Mat2 R) : Nat → Mat2 R
  | 0 => one
  | n + 1 => mul x (pow x n)
end Mat2

open Mat2

/-! ## A. The dial gas -/

/-- One edge of the dial gas: the two-site lattice gas with pair weight k is the lift of
    the weighted mult-add of the two fugacities. -/
theorem edge_gas (k x y : R) : 1 + wAdd k x y = 1 + x + y + k * x * y := by
  unfold wAdd; grind

/-- At k = 1 (q = 2) the two sites are independent: the gas factors. -/
theorem edge_gas_free (x y : R) : 1 + wAdd 1 x y = (1 + x) * (1 + y) := by
  unfold wAdd; grind

/-- At k = 0 (q = 1) the two sites exclude each other: hard core. -/
theorem edge_gas_hardcore (x y : R) : 1 + wAdd 0 x y = 1 + x + y := by
  unfold wAdd; grind

/-- The ring transfer matrix: from an empty site the next site weighs 1 (empty) or z
    (occupied); from an occupied site, 1 or k·z. -/
def ringT (k z : R) : Mat2 R := ⟨1, z, 1, k * z⟩

/-- Its trace is the weighted lift of z. -/
theorem ringT_tr (k z : R) : tr (ringT k z) = 1 + k * z := by
  unfold ringT tr; rfl

theorem ringT_det (k z : R) : det (ringT k z) = z * (k - 1) := by
  unfold ringT det; grind

/-- Cayley–Hamilton for 2×2 matrices: x² = tr(x)·x − det(x)·I. -/
theorem cayley_hamilton (x : Mat2 R) :
    mul x x = ⟨tr x * x.a - det x, tr x * x.b, tr x * x.c, tr x * x.d - det x⟩ := by
  apply Mat2.ext' <;> simp only [mul, tr, det] <;> grind

/-- With trace 0, every odd power is a multiple of x: x^(2m+1) = (−det x)^m · x. -/
theorem odd_pow_of_tr_zero (x : Mat2 R) (h : tr x = 0) (m : Nat) :
    pow x (2 * m + 1) = smul ((-det x) ^ m) x := by
  induction m with
  | zero =>
    apply Mat2.ext' <;> simp only [pow, mul, one, smul] <;> grind
  | succ m ih =>
    have e : 2 * (m + 1) + 1 = (2 * m + 1) + 1 + 1 := by omega
    rw [e]
    simp only [pow]
    simp only [pow] at ih
    rw [ih]
    have hd : x.d = -x.a := by simp only [tr] at h; grind
    apply Mat2.ext' <;> simp only [mul, smul, det] <;> rw [hd] <;> grind

/-- Odd powers of a traceless 2×2 matrix are traceless. -/
theorem tr_odd_pow_zero (x : Mat2 R) (h : tr x = 0) (m : Nat) : tr (pow x (2 * m + 1)) = 0 := by
  rw [odd_pow_of_tr_zero x h m]
  simp only [tr, smul] at h ⊢
  grind

/-- **The wall theorem for rings.** The dial gas on a ring of 2m + 1 sites, tr(Tⁿ), vanishes
    whenever 1 + k·z = 0, that is at z = −1/k, for every weight k. -/
theorem odd_ring_vanishes_on_wall (k z : R) (hw : 1 + k * z = 0) (m : Nat) :
    tr (pow (ringT k z) (2 * m + 1)) = 0 :=
  tr_odd_pow_zero _ (by rw [ringT_tr]; exact hw) m

/-- Even powers of a traceless 2×2 matrix are scalar: x^(2m) = (−det x)^m · I. -/
theorem even_pow_of_tr_zero (x : Mat2 R) (h : tr x = 0) (m : Nat) :
    pow x (2 * m) = smul ((-det x) ^ m) one := by
  induction m with
  | zero => apply Mat2.ext' <;> simp only [pow, one, smul] <;> grind
  | succ m ih =>
    have e : 2 * (m + 1) = 2 * m + 1 + 1 := by omega
    rw [e]
    simp only [pow]
    rw [ih]
    have hd : x.d = -x.a := by simp only [tr] at h; grind
    apply Mat2.ext' <;> simp only [mul, smul, one, det] <;> rw [hd] <;> grind

/-- On the wall an even ring is 2·(z(1 − k))^m, not zero (unless k = 1). -/
theorem even_ring_on_wall (k z : R) (hw : 1 + k * z = 0) (m : Nat) :
    tr (pow (ringT k z) (2 * m)) = 2 * (z * (1 - k)) ^ m := by
  rw [even_pow_of_tr_zero _ (by rw [ringT_tr]; exact hw) m, ringT_det]
  simp only [tr, smul, one]
  have : -(z * (k - 1)) = z * (1 - k) := by grind
  rw [this]; grind

/-- The transfer matrix counts correctly: the 3-ring (triangle) is
    1 + 3z + 3k·z² + k³·z³ (each pair of occupied sites shares an edge). -/
theorem triangle_count (k z : R) :
    tr (pow (ringT k z) 3) = 1 + 3 * z + 3 * k * z ^ 2 + k ^ 3 * z ^ 3 := by
  simp only [pow, mul, one, ringT, tr]; grind

/-- The 4-ring: 1 + 4z + (2 + 4k)z² + 4k²z³ + k⁴z⁴. -/
theorem square_count (k z : R) :
    tr (pow (ringT k z) 4) = 1 + 4 * z + (2 + 4 * k) * z ^ 2 + 4 * k ^ 2 * z ^ 3 + k ^ 4 * z ^ 4 := by
  simp only [pow, mul, one, ringT, tr]; grind

/-- The 5-ring: 1 + 5z + (5 + 5k)z² + (5k + 5k²)z³ + 5k³z⁴ + k⁵z⁵. -/
theorem pentagon_count (k z : R) :
    tr (pow (ringT k z) 5) =
      1 + 5 * z + (5 + 5 * k) * z ^ 2 + (5 * k + 5 * k ^ 2) * z ^ 3 + 5 * k ^ 3 * z ^ 4
        + k ^ 5 * z ^ 5 := by
  simp only [pow, mul, one, ringT, tr]; grind

/-- The triangle splits off the wall: (1 + kz) times a quadratic. -/
theorem triangle_factor (k z : R) :
    1 + 3 * z + 3 * k * z ^ 2 + k ^ 3 * z ^ 3 = (1 + k * z) * (1 + (3 - k) * z + k ^ 2 * z ^ 2) := by
  grind

/-- The quadratic's discriminant: (3 − k)² − 4k² = 3(1 − k)(3 + k). It is ≤ 0 for k ≥ 1
    (the two other zeros are a conjugate pair with product 1/k², on the circle |z| = 1/k:
    Lee–Yang), and > 0 for −3 < k < 1 (real zeros). -/
theorem triangle_disc (k : R) : (3 - k) ^ 2 - 4 * k ^ 2 = 3 * (1 - k) * (3 + k) := by
  grind

theorem triangle_disc_nonpos (k : Int) (hk : 1 ≤ k) : (3 - k) ^ 2 - 4 * k ^ 2 ≤ 0 := by
  rw [triangle_disc]
  have h1 : 1 - k ≤ 0 := by omega
  have h2 : 0 ≤ 3 + k := by omega
  have : 3 * (1 - k) * (3 + k) = -(3 * (k - 1) * (3 + k)) := by grind
  rw [this]
  have : 0 ≤ 3 * (k - 1) * (3 + k) :=
    Int.mul_nonneg (Int.mul_nonneg (by omega) (by omega)) h2
  omega

/-- The ring's eigenvalue discriminant: (1 + kz)² − 4z(k − 1) = k²z² + (4 − 2k)z + 1.
    Its zeros are the ends of the zero arc (Yang–Lee edges). At k = 1 it is (1 + z)². -/
theorem ring_disc (k z : R) :
    (1 + k * z) ^ 2 - 4 * (z * (k - 1)) = k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 := by
  grind

theorem ring_disc_free (z : R) : (1 + 1 * z) ^ 2 - 4 * (z * (1 - 1)) = (1 + z) ^ 2 := by
  grind

/-- The ring's zero family. A ring of n sites vanishes exactly when the two eigenvalues have
    ratio e^{iφ} with φ = π(2j+1)/n, which is (1 + kz)² = c·z(k − 1) with c = 4cos²(φ/2) in
    [0, 4]. Each c gives a quadratic in z: -/
theorem ring_zero_family (k z c : R) :
    (1 + k * z) ^ 2 - c * (z * (k - 1)) = k ^ 2 * z ^ 2 + (2 * k - c * (k - 1)) * z + 1 := by
  grind

/-- c = 0 is the wall, doubled: (1 + kz)². -/
theorem ring_family_wall (k z : R) : (1 + k * z) ^ 2 - 0 * (z * (k - 1)) = (1 + k * z) ^ 2 := by
  grind

/-- c = 4 is the edge of the zero arc (the eigenvalue discriminant). -/
theorem ring_family_edge (k z : R) :
    (1 + k * z) ^ 2 - 4 * (z * (k - 1)) = k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 := ring_disc k z

/-- In ζ = k·z the family reads k·(ζ² + 1) + (2k − c(k − 1))·ζ = 0 after multiplying by k: so
    ζ + 1/ζ = c(1 − 1/k) − 2 runs from −2 (the wall, c = 0) to 2 − 4/k (c = 4). -/
theorem ring_family_zeta (k z c : R) :
    k * (k ^ 2 * z ^ 2 + (2 * k - c * (k - 1)) * z + 1)
      = k * ((k * z) ^ 2 + 1) + (2 * k - c * (k - 1)) * (k * z) := by
  grind

/-- The two edge points z₁, z₂ of the zero arc satisfy k²·z₁z₂ = 1 (Vieta), so in ζ = k·z they
    are a reciprocal pair, symmetric under ζ ↦ 1/ζ (the up–down flip of the spins). -/
theorem ring_edges_vieta (k z z1 z2 : R) (hp : k ^ 2 * (z1 * z2) = 1)
    (hs : k ^ 2 * (z1 + z2) = 2 * k - 4) :
    k ^ 2 * z ^ 2 + (4 - 2 * k) * z + 1 = k ^ 2 * (z - z1) * (z - z2) := by
  grind

/-! ## B. The Potts ring: the dial is k = Q − 1 -/

/-- Partition function of the Q-state Potts model on a ring of n sites with edge weight
    1 + v·δ (eigenvalues of the edge matrix: Q + v once, v with multiplicity Q − 1). -/
def pottsRing (Q v : R) (n : Nat) : R := (Q + v) ^ n + (Q - 1) * v ^ n

/-- For Q = 2 (Ising) the formula is the trace of the transfer matrix [[1+v, 1], [1, 1+v]]. -/
theorem ising_ring_transfer (v : R) (n : Nat) :
    tr (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n) = pottsRing 2 v n := by
  have key : ∀ n, (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).a + (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).b
        = (2 + v) ^ n ∧
      (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).a - (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).b = v ^ n ∧
      (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).d = (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).a ∧
      (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).c = (pow (⟨1 + v, 1, 1, 1 + v⟩ : Mat2 R) n).b := by
    intro n
    induction n with
    | zero => simp only [pow, one]; grind
    | succ n ih =>
      obtain ⟨h1, h2, h3, h4⟩ := ih
      simp only [pow, mul]
      refine ⟨?_, ?_, ?_, ?_⟩
      · rw [h4]; have : (2 + v) ^ (n + 1) = (2 + v) * (2 + v) ^ n := by grind
        rw [this, ← h1]; grind
      · rw [h4]; have : v ^ (n + 1) = v * v ^ n := by grind
        rw [this, ← h2]; grind
      · rw [h3, h4]; grind
      · rw [h3, h4]; grind
  obtain ⟨h1, h2, h3, _⟩ := key n
  simp only [tr, pottsRing]
  rw [h3]
  have : (2 : R) - 1 = 1 := by grind
  rw [this]; grind

/-- **The Potts ring is the weighted lift of tⁿ.** With t = v/(Q + v) and k = Q − 1:
    Z = (Q + v)ⁿ·(1 + k·tⁿ). -/
theorem pottsRing_lift {K : Type} [Lean.Grind.Field K] (Q v t : K) (n : Nat)
    (ht : t * (Q + v) = v) :
    pottsRing Q v n = (Q + v) ^ n * (1 + (Q - 1) * t ^ n) := by
  unfold pottsRing
  have : v ^ n = (t * (Q + v)) ^ n := by rw [ht]
  rw [this]
  have hp : ∀ m : Nat, (t * (Q + v)) ^ m = t ^ m * (Q + v) ^ m := by
    intro m; induction m with
    | zero => grind
    | succ m ih => grind
  rw [hp]; grind

/-- Two separate rings together: their weights multiply, which in the lift is the weighted
    mult-add of the two tⁿ, at the dial k = Q − 1. -/
theorem two_rings (Q t A B : R) (n m : Nat) :
    (A * (1 + (Q - 1) * t ^ n)) * (B * (1 + (Q - 1) * t ^ m))
      = (A * B) * (1 + (Q - 1) * wAdd (Q - 1) (t ^ n) (t ^ m)) := by
  have h := wAdd_lift (Q - 1) (t ^ n) (t ^ m)
  rw [h]; grind

/-- Proper colourings: v = −1. The chromatic polynomial of the n-cycle. -/
theorem cycle_chromatic (Q : R) (n : Nat) :
    pottsRing Q (-1) n = (Q - 1) ^ n + (-1) ^ n * (Q - 1) := by
  unfold pottsRing
  have h : Q + -1 = Q - 1 := by grind
  rw [h]; grind

/-- Colourings sit on the wall: v = −1 gives t·(Q + v) = v with t = −1/(Q − 1), i.e.
    1 + k·t = 0 for k = Q − 1. (Stated without division: (Q − 1)·t = −1 and t(Q − 1) = −1.) -/
theorem colourings_on_wall (Q t : R) (h : (Q - 1) * t = -1) : t * (Q + (-1)) = -1 := by
  grind

/-- q = 0 (the "or" point, k = −1) is a zero of every ring. -/
theorem ring_zero_at_q0 (n : Nat) : pottsRing (0 : R) (-1) n = 0 := by
  unfold pottsRing
  have h1 : (0 : R) + -1 = -1 := by grind
  have h2 : (0 : R) - 1 = -1 := by grind
  rw [h1, h2]; grind

/-- q = 1 (plain addition, k = 0) is a zero of every ring with at least one site. -/
theorem ring_zero_at_q1 (n : Nat) : pottsRing (1 : R) (-1) (n + 1) = 0 := by
  unfold pottsRing
  have h1 : (1 : R) + -1 = 0 := by grind
  have h2 : (1 : R) - 1 = 0 := by grind
  rw [h1, h2]
  have : (0 : R) ^ (n + 1) = 0 := by grind
  rw [this]; grind

theorem one_pow' (n : Nat) : (1 : R) ^ n = 1 := by
  induction n with
  | zero => grind
  | succ n ih =>
    have : (1 : R) ^ (n + 1) = 1 ^ n * 1 := by grind
    rw [this, ih]; grind

/-- q = 2 (mult-add, k = 1) is a zero exactly for odd rings. -/
theorem ring_zero_at_q2_odd (m : Nat) : pottsRing (2 : R) (-1) (2 * m + 1) = 0 := by
  unfold pottsRing
  have h1 : ((2 : R) + -1) = 1 := by grind
  have h2 : ((2 : R) - 1) = 1 := by grind
  rw [h1, h2]
  have : ((-1 : R)) ^ (2 * m + 1) = -1 := by
    induction m with
    | zero => grind
    | succ m ih =>
      have e : 2 * (m + 1) + 1 = (2 * m + 1) + 2 := by omega
      rw [e]
      have : ((-1 : R)) ^ (2 * m + 1 + 2) = (-1) ^ (2 * m + 1) * ((-1) * (-1)) := by grind
      rw [this, ih]; grind
  have h1' : (1 : R) ^ (2 * m + 1) = 1 := one_pow' _
  rw [this, h1']; grind

theorem ring_at_q2_even (m : Nat) : pottsRing (2 : R) (-1) (2 * m) = 2 := by
  unfold pottsRing
  have h1 : ((2 : R) + -1) = 1 := by grind
  have h2 : ((2 : R) - 1) = 1 := by grind
  rw [h1, h2]
  have : ((-1 : R)) ^ (2 * m) = 1 := by
    induction m with
    | zero => grind
    | succ m ih =>
      have e : 2 * (m + 1) = 2 * m + 2 := by omega
      rw [e]
      have : ((-1 : R)) ^ (2 * m + 2) = (-1) ^ (2 * m) * ((-1) * (-1)) := by grind
      rw [this, ih]; grind
  have h1' : (1 : R) ^ (2 * m) = 1 := one_pow' _
  rw [this, h1']; grind

/-- q = −1 (parity, k = −2): the ring gives (−1)ⁿ(2ⁿ − 2), the signed number of acyclic
    orientations of the cycle. -/
theorem ring_at_parity (n : Nat) : pottsRing (-1 : R) (-1) n = (-1) ^ n * (2 ^ n - 2) := by
  unfold pottsRing
  have h1 : ((-1 : R) + -1) = (-1) * 2 := by grind
  rw [h1]
  have hp : ∀ m : Nat, ((-1 : R) * 2) ^ m = (-1) ^ m * 2 ^ m := by
    intro m; induction m with
    | zero => grind
    | succ m ih => grind
  rw [hp]; grind

end Ring

/-- The parity point is never a zero of a ring with two or more sites (over Int). -/
theorem two_pow_ge_four (n : Nat) (hn : 2 ≤ n) : (4 : Int) ≤ 2 ^ n := by
  induction n with
  | zero => omega
  | succ m ih =>
    rcases Nat.lt_or_ge m 2 with hm | hm
    · rcases (show m = 1 by omega) with h; subst h; decide
    · have := ih hm; rw [Int.pow_succ]; omega

theorem neg_one_pow_cases (n : Nat) : ((-1 : Int)) ^ n = 1 ∨ ((-1 : Int)) ^ n = -1 := by
  induction n with
  | zero => left; rfl
  | succ m ih =>
    rw [Int.pow_succ]
    rcases ih with h | h <;> rw [h] <;> decide

theorem parity_never_zero (n : Nat) (hn : 2 ≤ n) : pottsRing (-1 : Int) (-1) n ≠ 0 := by
  rw [ring_at_parity]
  have h2 := two_pow_ge_four n hn
  intro h
  rcases neg_one_pow_cases n with hs | hs <;> rw [hs] at h <;> omega

end LeeYang
