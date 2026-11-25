module

public import FusedAlgebra.Topics.QAcross.Knots

/-!
# Partitions, theta and primes (the counting `q`)

Power series in `q` are truncated here: a list of `N` integer coefficients (`QKnots.Poly`) means
"mod `q ^ N`". So every theorem of this file is a **finite check** up to the stated power, done
by the kernel (`decide +kernel`). The full theorems are standard (Euler 1750, Jacobi 1829,
Ramanujan 1919) and are not proved here.

| statement | checked |
|---|---|
| Euler's pentagonal theorem, `∏ (1 - qⁿ) = ∑ (-1) ^ k * q ^ (k * (3 * k - 1) / 2)` | mod `q ^ 30` |
| Euler's product, `∏ (1 - q ^ n) * ∑ p(n) * q ^ n = 1` | mod `q ^ 30` |
| the partition numbers `p(0), …, p(50)` | computed exactly |
| Ramanujan's congruences mod `5`, `7`, `11` | up to `p(50)` |
| Jacobi's triple product, `∏ (1 - q²ⁿ) * (1 + q²ⁿ⁻¹) ^ 2 = ∑ q ^ (k ^ 2)` | mod `q ^ 30` |
| Jacobi's four-square theorem | mod `q ^ 30` |

In mathlib the partition numbers are `Fintype.card (Nat.Partition n)`, with generating function
`Nat.Partition.hasProd_genFun`; Lagrange's four-square theorem is `Nat.sum_four_squares`. The
pentagonal theorem, the triple product, Jacobi's four-square count and Ramanujan's congruences
are not in mathlib. This file does not connect its computed table `partList` to
`Nat.Partition`.

The mult-add link: `∏ (1 - q ^ n)` is `1 +` the mult-add of the list `-q, -q ^ 2, -q ^ 3, …`
(`one_add_multAddList_neg_pow`, an instance of `MultAdd.one_add_multAddList`, valid in any
commutative ring; `euler_is_multAdd_30` is the same statement computed on truncated series). Its
corner expansion is a signed sum over sets of distinct parts, and the pentagonal theorem says
those corners cancel down to the pentagonal numbers.
-/

@[expose] public section

namespace QPartitions

open MultAdd
open QKnots (Poly padd pscale pmul pzeros trim)

/-! ### Truncated series -/

/-- The series `1` (mod `q ^ N`). -/
def one (N : ℕ) : Poly := 1 :: List.replicate (N - 1) 0

/-- The product of two series, mod `q ^ N`. -/
def tmul (N : ℕ) (p q : Poly) : Poly := (pmul p q).take N

/-- `p * (1 - q ^ k)` mod `q ^ N`. -/
def oneMinus (N k : ℕ) (p : Poly) : Poly := (padd p (pscale (-1) (pzeros k p))).take N

/-- `p * (1 + q ^ k)` mod `q ^ N`. -/
def onePlus (N k : ℕ) (p : Poly) : Poly := (padd p (pzeros k p)).take N

/-- The loop of `divOneMinus`: each new coefficient is the old one plus the coefficient `k`
places back; `buf` holds the last `k` results. -/
def divAux : Poly → Poly → Poly
  | [], _ => []
  | a :: as, buf => (a + buf.headD 0) :: divAux as (buf.tail ++ [a + buf.headD 0])

/-- `p / (1 - q ^ k)`, to the length of `p`. -/
def divOneMinus (k : ℕ) (p : Poly) : Poly := divAux p (List.replicate k 0)

/-! ### Euler's product and the pentagonal numbers -/

/-- `∏ k ∈ [1, N), (1 - q ^ k)` mod `q ^ N`. -/
def eulerProd (N : ℕ) : Poly :=
  (List.range' 1 (N - 1)).foldl (fun acc k ↦ oneMinus N k acc) (one N)

/-- The coefficient of `q ^ m` in `∑ k : ℤ, (-1) ^ k * q ^ (k * (3 * k - 1) / 2)`: `1` or `-1` at
the pentagonal numbers `0, 1, 2, 5, 7, 12, 15, 22, 26, …`, else `0`. -/
def pentCoeff (m : ℕ) : ℤ :=
  (List.range (m + 2)).foldl (fun acc k ↦
    acc + (if k * (3 * k - 1) / 2 = m then (-1) ^ k else 0)
        + (if 1 ≤ k ∧ k * (3 * k + 1) / 2 = m then (-1) ^ k else 0)) 0

/-- Euler's pentagonal theorem, checked mod `q ^ 30`. -/
theorem pentagonal_30 : eulerProd 30 = (List.range 30).map pentCoeff := by decide +kernel

theorem eulerProd_30_explicit :
    eulerProd 30 = [1, -1, -1, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0] := by decide +kernel

/-- Mult-add on truncated series: `x + y + x * y` mod `q ^ N`. -/
def tma (N : ℕ) (x y : Poly) : Poly := (padd (padd x y) (pmul x y)).take N

/-- The mult-add of `-q, -q ^ 2, …, -q ^ (N - 1)`, mod `q ^ N`. -/
def eulerMA (N : ℕ) : Poly :=
  (List.range' 1 (N - 1)).foldl (fun acc k ↦ tma N acc (pzeros k [-1])) []

/-- `1 + (-q ✳ -q ^ 2 ✳ ⋯)` is Euler's product, computed mod `q ^ 30`. -/
theorem euler_is_multAdd_30 : trim (padd [1] (eulerMA 30)) = trim (eulerProd 30) := by
  decide +kernel

/-- Euler's product is `1 +` the mult-add of `-q ^ k`: the list lift `one_add_multAddList`, in
any ring and for any list of exponents. -/
theorem one_add_multAddList_neg_pow {R : Type*} [Ring R] (q : R) (l : List ℕ) :
    1 + multAddList (l.map fun k ↦ -q ^ k) = (l.map fun k ↦ 1 - q ^ k).prod := by
  rw [one_add_multAddList, List.map_map]
  simp only [Function.comp_def, sub_eq_add_neg]

/-! ### Partition numbers -/

/-- `∑ p(n) * q ^ n = ∏ 1 / (1 - q ^ k)`, mod `q ^ N`. -/
def partGen (N : ℕ) : Poly :=
  (List.range' 1 (N - 1)).foldl (fun acc k ↦ divOneMinus k acc) (one N)

/-- The partition numbers `p(0), …, p(50)`. -/
def partList : Poly :=
  [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792,
   1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977,
   21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273,
   173525, 204226]

theorem partitions_to_50 : partGen 51 = partList := by decide +kernel

/-- Euler's product times the partition series is `1` (mod `q ^ 30`): the two are inverse. -/
theorem euler_inverse_30 : tmul 30 (eulerProd 30) (partGen 30) = one 30 := by decide +kernel

/-- Ramanujan: `p(5 * n + 4) ≡ 0 [ZMOD 5]`, for `p(4), p(9), …, p(49)`. -/
theorem ramanujan_5 : ∀ n : Fin 10, partList.getD (5 * n.val + 4) 0 % 5 = 0 := by decide

/-- Ramanujan: `p(7 * n + 5) ≡ 0 [ZMOD 7]`, for `p(5), p(12), …, p(47)`. -/
theorem ramanujan_7 : ∀ n : Fin 7, partList.getD (7 * n.val + 5) 0 % 7 = 0 := by decide

/-- Ramanujan: `p(11 * n + 6) ≡ 0 [ZMOD 11]`, for `p(6), p(17), p(28), p(39), p(50)`. -/
theorem ramanujan_11 : ∀ n : Fin 5, partList.getD (11 * n.val + 6) 0 % 11 = 0 := by decide

/-- The congruences are special to `5`, `7`, `11`: for the prime `13` the same pattern
`p(13 * n + r) ≡ 0` already fails for every residue `r` by `n ≤ 2`. -/
theorem no_ramanujan_13 :
    ∀ r : Fin 13, ∃ n : Fin 3, partList.getD (13 * n.val + r.val) 0 % 13 ≠ 0 := by decide

/-! ### Theta, the triple product and four squares -/

/-- `θ(q) = ∑ k : ℤ, q ^ (k ^ 2) = 1 + 2 * q + 2 * q ^ 4 + 2 * q ^ 9 + ⋯`, mod `q ^ N`. -/
def theta (N : ℕ) : Poly :=
  (List.range N).map fun i ↦ if i = 0 then 1 else if IsSquare i then 2 else 0

/-- Jacobi's triple product `∏ (1 - q ^ (2 * n)) * (1 + q ^ (2 * n - 1)) ^ 2`, mod `q ^ N`. -/
def tripleProd (N : ℕ) : Poly :=
  (List.range' 1 (N - 1)).foldl
    (fun acc k ↦ oneMinus N (2 * k) (onePlus N (2 * k - 1) (onePlus N (2 * k - 1) acc))) (one N)

/-- The product side equals `θ`, mod `q ^ 30`. -/
theorem triple_product_30 : tripleProd 30 = theta 30 := by decide +kernel

/-- The sum of the divisors `d` of `n` with `¬ 4 ∣ d`. -/
def sigmaStar (n : ℕ) : ℕ := ∑ d ∈ n.divisors with ¬ 4 ∣ d, d

/-- Jacobi's count: `r₄(n) = 8 * σ*(n)` for `n ≥ 1`, and `r₄(0) = 1`. -/
def r4 (n : ℕ) : ℤ := if n = 0 then 1 else 8 * (sigmaStar n : ℤ)

/-- `θ ^ 4` counts the ways to write `n` as a sum of four squares (signs and order counted);
Jacobi's theorem says that count is `8 *` (the sum of the divisors not divisible by `4`).
Checked mod `q ^ 30`. -/
theorem jacobi_four_squares_30 :
    tmul 30 (tmul 30 (theta 30) (theta 30)) (tmul 30 (theta 30) (theta 30)) =
      (List.range 30).map r4 := by
  decide +kernel

/-- For a prime `p`, `r₄(p) = 8 * (p + 1)`: the only divisors of a prime are `1` and itself.
(The original checked the primes below `30`; this is every prime.) -/
theorem r4_prime {p : ℕ} (hp : p.Prime) : r4 p = 8 * (p + 1) := by
  have h4 : ¬ 4 ∣ p := fun h ↦ by
    have := (Nat.prime_dvd_prime_iff_eq Nat.prime_two hp).1 (dvd_trans (by decide) h)
    subst this; revert h; decide
  have h1 : (1 : ℕ) ≠ p := hp.one_lt.ne
  rw [r4, ite_eq_right hp.ne_zero, sigmaStar, hp.divisors, Finset.sum_filter,
    Finset.sum_pair h1, ite_eq_left (by decide), ite_eq_left h4]
  push_cast; ring

theorem four_squares_primes :
    r4 2 = 24 ∧ r4 3 = 32 ∧ r4 5 = 48 ∧ r4 7 = 64 ∧ r4 11 = 96 ∧ r4 13 = 112 ∧
    r4 17 = 144 ∧ r4 19 = 160 ∧ r4 23 = 192 ∧ r4 29 = 240 := by decide +kernel

/-- Jacobi's count is positive, for every `n`. With Jacobi's theorem this is Lagrange's
four-square theorem, which mathlib proves directly as `Nat.sum_four_squares`. (The original
checked `n < 30`.) -/
theorem r4_pos (n : ℕ) : 0 < r4 n := by
  rw [r4]
  split_ifs with h
  · exact one_pos
  · have h1 : 1 ∈ n.divisors.filter (fun d ↦ ¬ 4 ∣ d) := by
      simp [h]
    have : 0 < sigmaStar n := Finset.sum_pos' (fun _ _ ↦ Nat.zero_le _) ⟨1, h1, one_pos⟩
    positivity

/-- Every `n < 30` is a sum of four squares (Lagrange's theorem, read off Jacobi's count). -/
theorem lagrange_below_30 : ∀ n : Fin 30, 0 < r4 n.val := fun n ↦ r4_pos n

end QPartitions
