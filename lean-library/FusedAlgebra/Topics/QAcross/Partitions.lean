/-
  Fused Hyperoperation Algebra: partitions, theta and primes (the counting q)
  ===========================================================================
  Core Lean 4 only (no mathlib), no `sorry`.

  LIBRARY NOTE: copied from q-across/lean/Partitions.lean (Nov 25) and compiled here for the
  first time: no sorry, standard axioms only; no changes were needed. No statement changed.

  Power series in q are truncated: a list of N coefficients means "mod q^N". So every theorem
  here is a FINITE CHECK up to the stated power. The full theorems are standard and proved
  elsewhere (Euler 1750, Jacobi 1829, Ramanujan 1919); they are not proved here.

    Euler's pentagonal theorem   Π(1 − qⁿ) = Σ (−1)ᵏ q^(k(3k−1)/2)          checked mod q³⁰
    Euler's product              Π(1 − qⁿ) · Σ p(n) qⁿ = 1                   checked mod q³⁰
    Partition numbers            p(0..50)                                    computed exactly
    Ramanujan's congruences      p(5n+4) ≡ 0 mod 5, p(7n+5) ≡ 0 mod 7,
                                 p(11n+6) ≡ 0 mod 11                         checked to p(50)
    Jacobi's triple product      Π(1 − q²ⁿ)(1 + q²ⁿ⁻¹)² = θ(q) = Σ q^(k²)   checked mod q³⁰
    Jacobi's four-square theorem θ(q)⁴ = 1 + 8 Σₙ (Σ_{d|n, 4∤d} d) qⁿ       checked mod q³⁰

  The mult-add link: Π(1 − qⁿ) is 1 + (the mult-add of the list −q, −q², −q³, …), by the
  library's foldCirc_lift. Its corner expansion is a signed sum over sets of distinct parts, and
  the pentagonal theorem says those corners cancel down to the pentagonal numbers (Franklin's
  picture: slide the top row of a Ferrers diagram onto its staircase edge). The q here is the
  Tate gluing q: η(τ) = q^(1/24) Π(1 − qⁿ), θ and the partition series all live on the Tate
  family, q = e^(2πiτ).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.QAcross.Knots

namespace QPartitions
open CircleOp
open QKnots (Poly padd pscale pmul pzeros trim negOnePow)

/-! ## Truncated series -/

/-- The series 1 (mod q^N). -/
def one (N : Nat) : Poly := 1 :: List.replicate (N - 1) 0

def tmul (N : Nat) (p q : Poly) : Poly := (pmul p q).take N

/-- p·(1 − q^k) mod q^N. -/
def oneMinus (N k : Nat) (p : Poly) : Poly := (padd p (pscale (-1) (pzeros k p))).take N

/-- p·(1 + q^k) mod q^N. -/
def onePlus (N k : Nat) (p : Poly) : Poly := (padd p (pzeros k p)).take N

/-- p / (1 − q^k): each new coefficient is the old one plus the coefficient k places back.
    `buf` holds the last k results. -/
def divAux : Poly → Poly → Poly
  | [], _ => []
  | a :: as, buf => (a + buf.headD 0) :: divAux as (buf.tail ++ [a + buf.headD 0])

def divOneMinus (k : Nat) (p : Poly) : Poly := divAux p (List.replicate k 0)

/-! ## Euler's product and the pentagonal numbers -/

/-- Π_{k=1}^{N−1} (1 − q^k) mod q^N. -/
def eulerProd (N : Nat) : Poly := (List.range' 1 (N - 1)).foldl (fun acc k => oneMinus N k acc) (one N)

/-- Coefficient of q^m in Σ_{k∈ℤ} (−1)ᵏ q^(k(3k−1)/2): +1 or −1 at the pentagonal numbers
    0, 1, 2, 5, 7, 12, 15, 22, 26, …, else 0. -/
def pentCoeff (m : Nat) : Int :=
  (List.range (m + 2)).foldl (fun acc k =>
    acc + (if k * (3 * k - 1) / 2 == m then negOnePow k else 0)
        + (if decide (1 ≤ k) && k * (3 * k + 1) / 2 == m then negOnePow k else 0)) 0

/-- Euler's pentagonal theorem, checked mod q³⁰. -/
theorem pentagonal_30 : eulerProd 30 = (List.range 30).map pentCoeff := by decide +kernel

theorem eulerProd_30_explicit :
    eulerProd 30 = [1, -1, -1, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0] := by decide +kernel

/-- The same product as a mult-add: x ◯ y = x + y + xy on truncated series. -/
def tma (N : Nat) (x y : Poly) : Poly := (padd (padd x y) (pmul x y)).take N

/-- −q ◯ −q² ◯ … ◯ −q^(N−1) mod q^N. -/
def eulerCirc (N : Nat) : Poly := (List.range' 1 (N - 1)).foldl (fun acc k => tma N acc (pzeros k [-1])) []

/-- 1 + (−q ◯ −q² ◯ …) is Euler's product (the list lift foldCirc_lift, seen on series). -/
theorem euler_is_circ_30 : trim (padd [1] (eulerCirc 30)) = trim (eulerProd 30) := by
  decide +kernel

/-! ## Partition numbers -/

/-- Σ p(n) qⁿ = Π 1/(1 − q^k), mod q^N. -/
def partGen (N : Nat) : Poly := (List.range' 1 (N - 1)).foldl (fun acc k => divOneMinus k acc) (one N)

/-- p(0), …, p(50). -/
def partList : Poly :=
  [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792,
   1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977,
   21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273,
   173525, 204226]

theorem partitions_to_50 : partGen 51 = partList := by decide +kernel

/-- Euler's product times the partition series is 1 (mod q³⁰): the two are inverse. -/
theorem euler_inverse_30 : tmul 30 (eulerProd 30) (partGen 30) = one 30 := by decide +kernel

/-- Ramanujan: p(5n+4) ≡ 0 mod 5, for p(4), p(9), …, p(49). -/
theorem ramanujan_5 : ∀ n : Fin 10, partList.getD (5 * n.val + 4) 0 % 5 = 0 := by decide

/-- Ramanujan: p(7n+5) ≡ 0 mod 7, for p(5), p(12), …, p(47). -/
theorem ramanujan_7 : ∀ n : Fin 7, partList.getD (7 * n.val + 5) 0 % 7 = 0 := by decide

/-- Ramanujan: p(11n+6) ≡ 0 mod 11, for p(6), p(17), p(28), p(39), p(50). -/
theorem ramanujan_11 : ∀ n : Fin 5, partList.getD (11 * n.val + 6) 0 % 11 = 0 := by decide

/-- The congruences are special to 5, 7, 11: for the prime 13 the same pattern p(13n+r) ≡ 0
    already fails for every residue r by n ≤ 2 (13n + r ≤ 38). -/
theorem no_ramanujan_13 :
    ∀ r : Fin 13, ∃ n : Fin 3, partList.getD (13 * n.val + r.val) 0 % 13 ≠ 0 := by decide

/-! ## Theta, the triple product and four squares -/

def isSq (i : Nat) : Bool := (List.range (i + 1)).any (fun j => j * j == i)

/-- θ(q) = Σ_{k∈ℤ} q^(k²) = 1 + 2q + 2q⁴ + 2q⁹ + …, mod q^N. -/
def theta (N : Nat) : Poly := (List.range N).map (fun i => if i == 0 then 1 else if isSq i then 2 else 0)

/-- Jacobi's triple product: Π (1 − q²ⁿ)(1 + q²ⁿ⁻¹)², mod q^N. -/
def tripleProd (N : Nat) : Poly :=
  (List.range' 1 (N - 1)).foldl
    (fun acc k => oneMinus N (2 * k) (onePlus N (2 * k - 1) (onePlus N (2 * k - 1) acc))) (one N)

/-- The product side equals θ, mod q³⁰. (θ is the function the Tate gluing makes periodic.) -/
theorem triple_product_30 : tripleProd 30 = theta 30 := by decide +kernel

/-- Σ of the divisors d of n with 4 ∤ d. -/
def sigmaStar (n : Nat) : Nat :=
  ((List.range' 1 n).filter (fun d => n % d == 0 && d % 4 != 0)).foldl (fun x y => x + y) 0

/-- Jacobi's count: r₄(n) = 8·σ*(n) for n ≥ 1, r₄(0) = 1. -/
def r4 (n : Nat) : Int := if n == 0 then 1 else 8 * (sigmaStar n : Int)

/-- θ⁴ counts the ways to write n as a sum of four squares (signs and order counted); Jacobi's
    theorem says that count is 8 × (sum of divisors not divisible by 4). Checked mod q³⁰. -/
theorem jacobi_four_squares_30 :
    tmul 30 (tmul 30 (theta 30) (theta 30)) (tmul 30 (theta 30) (theta 30)) =
      (List.range 30).map r4 := by
  decide +kernel

/-- So for a prime p, r₄(p) = 8(p + 1): the prime's only divisors are 1 and p. -/
theorem four_squares_primes :
    r4 2 = 24 ∧ r4 3 = 32 ∧ r4 5 = 48 ∧ r4 7 = 64 ∧ r4 11 = 96 ∧ r4 13 = 112 ∧
    r4 17 = 144 ∧ r4 19 = 160 ∧ r4 23 = 192 ∧ r4 29 = 240 := by decide +kernel

/-- Every n < 30 is a sum of four squares (Lagrange's theorem, checked in this range). -/
theorem lagrange_below_30 : ∀ n : Fin 30, 0 < r4 n.val := by decide +kernel

end QPartitions
