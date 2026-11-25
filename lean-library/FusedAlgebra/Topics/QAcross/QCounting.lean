/-
  Fused Hyperoperation Algebra: q-counting, the quantum group and the Ising ring
  =============================================================================
  Core Lean 4 only (no mathlib), no `sorry`.

  LIBRARY NOTE: copied from q-across/lean/QCounting.lean (Nov 25) and compiled here for the
  first time: no sorry, standard axioms only; one fix: the section's `[CommRing R]` became `[Lean.Grind.CommRing R]`
  (core Lean's ring class; the bare name doesn't resolve without `open Lean.Grind`). No statement changed.

  1. Gaussian binomials [n, k]_q (standard): q-Pascal, and the GENERAL theorem that at q = 1
     they are the layers of the mult-add cube (gauss_at_one, every n and k). At q = 2 they count
     subspaces over the two-element field, at q = −1 they fold Pascal in half.
  2. The q-binomial theorem (standard, Cauchy/Gauss): (1 + t)(1 + qt)⋯(1 + q^(n−1) t), i.e.
     1 + (t ⊞ qt ⊞ … ⊞ q^(n−1) t), has layer k equal to q^(k(k−1)/2)·[n, k]_q. Checked n ≤ 6.
     Geometrically: the corners of the cube with k switches on, weighted by q^(sum of their
     positions), are the lattice paths in a k × (n − k) box weighted by the area under them.
  3. U_q(sl₂) (standard, Drinfeld–Jimbo): the relation [E, F] = (K − K⁻¹)/(q − q⁻¹) holds on
     every irreducible representation because of one ring identity, proved here for any
     commutative ring.
  4. The Ising ring (standard transfer matrix): with T = [[1, u], [u, 1]], the trace of Tⁿ is
     (1 + u)ⁿ + (1 − u)ⁿ = 2 + u^(⊞n) + (−u)^(⊞n), every n, any commutative ring.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube
import FusedAlgebra.Topics.QAcross.Knots

namespace QCounting
open MultAdd
open QKnots (Poly padd pzeros trim)

/-! ## 1. Gaussian binomials -/

/-- [n, k]_q as a polynomial in q, by q-Pascal: [n+1, k+1] = [n, k] + q^(k+1)·[n, k+1]. -/
def gauss : Nat → Nat → Poly
  | _, 0 => [1]
  | 0, _ + 1 => []
  | n + 1, k + 1 => padd (gauss n k) (pzeros (k + 1) (gauss n (k + 1)))

/-- Evaluate a polynomial at x. -/
def peval (x : Int) : Poly → Int
  | [] => 0
  | c :: p => c + x * peval x p

theorem peval_padd (x : Int) : ∀ p q : Poly, peval x (padd p q) = peval x p + peval x q
  | [], q => by
    show peval x q = 0 + peval x q
    omega
  | a :: p, [] => by
    show peval x (a :: p) = peval x (a :: p) + 0
    omega
  | a :: p, b :: q => by
    show (a + b) + x * peval x (padd p q) = (a + x * peval x p) + (b + x * peval x q)
    rw [peval_padd x p q]
    grind

theorem peval_pzeros (x : Int) : ∀ (k : Nat) (p : Poly), peval x (pzeros k p) = pw x k * peval x p
  | 0, p => by
    show peval x p = 1 * peval x p
    omega
  | k + 1, p => by
    show 0 + x * peval x (pzeros k p) = pw x k * x * peval x p
    rw [peval_pzeros x k p]
    grind

/-- At q = 1 the Gaussian binomial is the binomial coefficient: the layer of the mult-add cube.
    (Every n and k, not a finite check.) -/
theorem gauss_at_one : ∀ n k : Nat, peval 1 (gauss n k) = (MultAdd.Cube.layer n k : Int)
  | n, 0 => by
    rw [MultAdd.Cube.layer_zero_right]
    simp only [gauss, peval] <;> decide
  | 0, k + 1 => by
    simp only [gauss, peval, MultAdd.Cube.layer] <;> decide
  | n + 1, k + 1 => by
    simp only [gauss]
    rw [peval_padd, peval_pzeros, gauss_at_one n k, gauss_at_one n (k + 1),
      MultAdd.Cube.layer_succ_succ, pw_one]
    omega

/-- The four points of the dial on [4, 2]_q = 1 + q + 2q² + q³ + q⁴:
    q = 1 counts 6 = C(4,2); q = 2 counts the 35 planes in a 4-dimensional space over
    the two-element field; q = 0 gives 1; q = −1 gives 2 = C(2,1) (Pascal folded in half). -/
theorem gauss_4_2 :
    trim (gauss 4 2) = [1, 1, 2, 1, 1] ∧ peval 1 (gauss 4 2) = 6 ∧ peval 2 (gauss 4 2) = 35 ∧
    peval 0 (gauss 4 2) = 1 ∧ peval (-1) (gauss 4 2) = 2 := by
  decide +kernel

/-- At q = −1, [n, k] is C(⌊n/2⌋, ⌊k/2⌋), or 0 when n is even and k odd. Checked n ≤ 8. -/
theorem gauss_at_neg_one :
    ∀ n : Fin 9, ∀ k : Fin 9, k.val ≤ n.val →
      peval (-1) (gauss n.val k.val) =
        (if n.val % 2 == 0 && k.val % 2 == 1 then 0
         else (MultAdd.Cube.layer (n.val / 2) (k.val / 2) : Int)) := by
  decide +kernel

/-! ## 2. The q-binomial theorem is a mult-add fold -/

/-- Polynomials in t whose coefficients are polynomials in q. -/
def padd2 : List Poly → List Poly → List Poly
  | [], Q => Q
  | p :: P, [] => p :: P
  | p :: P, q :: Q => padd p q :: padd2 P Q

/-- Multiply by (1 + qⁱ·t). -/
def mulLin (i : Nat) (P : List Poly) : List Poly := padd2 P ([] :: P.map (pzeros i))

/-- (1 + t)(1 + qt)⋯(1 + q^(n−1) t) = 1 + (t ⊞ qt ⊞ … ⊞ q^(n−1) t). -/
def qbinLHS (n : Nat) : List Poly := (List.range n).foldl (fun P i => mulLin i P) [[1]]

/-- Σₖ q^(k(k−1)/2)·[n, k]_q·tᵏ. -/
def qbinRHS (n : Nat) : List Poly := (List.range (n + 1)).map (fun k => pzeros (k * (k - 1) / 2) (gauss n k))

theorem q_binomial_theorem : ∀ n : Fin 7, (qbinLHS n.val).map trim = (qbinRHS n.val).map trim := by
  decide +kernel

/-! ## 3. The quantum group U_q(sl₂) -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- On the irreducible representation of dimension n + 1, with basis v₀ … vₙ, E and F move
    between neighbours with factors [k]_q and [n − k]_q (symmetric q-integers
    [m] = (qᵐ − q⁻ᵐ)/(q − q⁻¹)), and K acts by q^(n − 2k). The relation
    [E, F] = (K − K⁻¹)/(q − q⁻¹) on v_k is [k+1][n−k] − [k][n−k+1] = [n − 2k]. Clearing the
    denominators (q − q⁻¹)² and writing x = qᵏ, y = q^(n−k) (xi, yi, qi their inverses), it is
    this ring identity, so the relation holds on every irreducible representation. -/
theorem uq_sl2_commutator (q qi x xi y yi : R) :
    (q * x - qi * xi) * (y - yi) - (x - xi) * (q * y - qi * yi) = (q - qi) * (xi * y - x * yi) := by
  grind

/-- K E K⁻¹ = q² E: E moves v_k to v_{k−1}, whose weight q^(n−2k+2) is q² times v_k's. -/
theorem uq_sl2_weight (q x : R) : q * q * x = q * (q * x) := by grind

/-- The classical limit q → 1: [E, F] = H, i.e. (a + 1)·b − a·(b + 1) = b − a. -/
theorem sl2_classical (a b : Int) : (a + 1) * b - a * (b + 1) = b - a := by grind

/-- The quantum dimension of the 2-dimensional representation is [2] = q + q⁻¹; with q = A²
    it is minus the knot loop value d = −A² − A⁻². -/
theorem quantum_dim_is_loop (A Ai : R) : A * A + Ai * Ai = -(-(A * A) - Ai * Ai) := by grind

/-! ## 4. The Ising ring -/

/-- One more site: multiply by the transfer matrix T = [[1, u], [u, 1]]
    (u = e^(−2βJ), the cost of a domain wall). Tⁿ = [[a, b], [b, a]]; we track (a, b). -/
def isingStep (u : R) (p : R × R) : R × R := (p.1 + u * p.2, u * p.1 + p.2)

def ising (u : R) : Nat → R × R
  | 0 => (1, 0)
  | n + 1 => isingStep u (ising u n)

/-- The two eigen-directions: a + b = (1 + u)ⁿ and a − b = (1 − u)ⁿ. -/
theorem ising_eigen (u : R) (n : Nat) :
    (ising u n).1 + (ising u n).2 = pw (1 + u) n ∧ (ising u n).1 - (ising u n).2 = pw (1 - u) n := by
  induction n with
  | zero =>
    exact ⟨by show (1 : R) + 0 = 1; grind, by show (1 : R) - 0 = 1; grind⟩
  | succ n ih =>
    obtain ⟨h1, h2⟩ := ih
    constructor
    · show (ising u n).1 + u * (ising u n).2 + (u * (ising u n).1 + (ising u n).2) = pw (1 + u) n * (1 + u)
      rw [← h1]; grind
    · show (ising u n).1 + u * (ising u n).2 - (u * (ising u n).1 + (ising u n).2) = pw (1 - u) n * (1 - u)
      rw [← h2]; grind

/-- The partition function of n spins on a ring (up to the factor e^(nβJ)) is the trace:
    Z = (1 + u)ⁿ + (1 − u)ⁿ. -/
theorem ising_trace (u : R) (n : Nat) : 2 * (ising u n).1 = pw (1 + u) n + pw (1 - u) n := by
  obtain ⟨h1, h2⟩ := ising_eigen u n
  grind

/-- ... and in mult-add terms: Z = 2 + u^(⊞n) + (−u)^(⊞n). Only the even layers of the cube
    survive (a ring has an even number of domain walls). -/
theorem ising_multAdd (u : R) (n : Nat) : 2 * (ising u n).1 = 2 + mpow u n + mpow (-u) n := by
  have h := ising_trace u n
  rw [pw_eq_pow, pw_eq_pow] at h
  have a := mpow_lift u n
  have b := mpow_lift (-u) n
  have e : (1 : R) + -u = 1 - u := by grind
  rw [e] at b
  grind

end Ring

end QCounting
