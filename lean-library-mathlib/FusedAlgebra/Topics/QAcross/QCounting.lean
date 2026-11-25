module

public import FusedAlgebra.Cube
public import FusedAlgebra.Topics.QAcross.Knots

/-!
# `q`-counting, the quantum group and the Ising ring

1. **Gaussian binomials** `[n, k]_q` (standard; not in mathlib). They are built by `q`-Pascal as
   coefficient lists (`gauss`), so that finite checks can be done by the kernel. The general
   theorem `gauss_at_one` says that at `q = 1` they are the binomial coefficients `n.choose k`,
   the layers of the mult-add cube, for every `n` and `k`. At `q = 2` they count subspaces over
   the field with two elements, at `q = -1` they fold Pascal's triangle in half.
2. **The `q`-binomial theorem** (standard, Cauchy and Gauss): the product
   `(1 + t) * (1 + q * t) * ⋯ * (1 + q ^ (n - 1) * t)`, which is
   `1 + (t ✳ q * t ✳ ⋯ ✳ q ^ (n - 1) * t)`, has layer `k` equal to
   `q ^ (k * (k - 1) / 2) * [n, k]_q`. This is a finite check for `n ≤ 6`
   (`q_binomial_theorem`).
3. **`U_q(sl₂)`** (standard, Drinfeld and Jimbo): the relation `[E, F] = (K - K⁻¹) / (q - q⁻¹)`
   holds on every irreducible representation because of one ring identity
   (`uq_sl2_commutator`), valid in any commutative ring.
4. **The Ising ring** (the standard transfer matrix): with `T = !![1, u; u, 1]`, the trace of
   `T ^ n` is `(1 + u) ^ n + (1 - u) ^ n = 2 + mpow u n + mpow (-u) n`, for every `n` and in any
   commutative ring (`trace_transfer_pow`, `ising_multAdd`).
-/

@[expose] public section

namespace QCounting

open MultAdd Polynomial
open QKnots (Poly padd pzeros trim toPolynomial)

/-! ### 1. Gaussian binomials -/

/-- The Gaussian binomial `[n, k]_q` as a polynomial in `q`, by `q`-Pascal:
`[n + 1, k + 1] = [n, k] + q ^ (k + 1) * [n, k + 1]`. -/
def gauss : ℕ → ℕ → Poly
  | _, 0 => [1]
  | 0, _ + 1 => []
  | n + 1, k + 1 => padd (gauss n k) (pzeros (k + 1) (gauss n (k + 1)))

@[simp]
theorem gauss_zero_right (n : ℕ) : gauss n 0 = [1] := by cases n <;> rfl

@[simp]
theorem gauss_zero_succ (k : ℕ) : gauss 0 (k + 1) = [] := rfl

/-- `q`-Pascal. -/
theorem gauss_succ_succ (n k : ℕ) :
    gauss (n + 1) (k + 1) = padd (gauss n k) (pzeros (k + 1) (gauss n (k + 1))) := rfl

/-- Evaluate a coefficient list at `x` (Horner's rule). This is the computable form of
`Polynomial.eval`; see `peval_eq_eval`. -/
def peval (x : ℤ) : Poly → ℤ
  | [] => 0
  | c :: p => c + x * peval x p

theorem peval_eq_eval (x : ℤ) (p : Poly) : peval x p = (toPolynomial p).eval x := by
  induction p with
  | nil => simp [peval]
  | cons c p ih => simp [peval, ih]

theorem peval_padd (x : ℤ) (p q : Poly) : peval x (padd p q) = peval x p + peval x q := by
  simp [peval_eq_eval, QKnots.toPolynomial_padd]

theorem peval_pzeros (x : ℤ) (k : ℕ) (p : Poly) : peval x (pzeros k p) = x ^ k * peval x p := by
  simp [peval_eq_eval, QKnots.toPolynomial_pzeros]

/-- `q`-Pascal for the values of the Gaussian binomials. -/
theorem peval_gauss_succ_succ (x : ℤ) (n k : ℕ) :
    peval x (gauss (n + 1) (k + 1)) =
      peval x (gauss n k) + x ^ (k + 1) * peval x (gauss n (k + 1)) := by
  rw [gauss_succ_succ, peval_padd, peval_pzeros]

/-- At `q = 1` the Gaussian binomial is the binomial coefficient: the layer of the mult-add
cube. (Every `n` and `k`, not a finite check.) -/
theorem gauss_at_one : ∀ n k : ℕ, peval 1 (gauss n k) = n.choose k
  | n, 0 => by simp [peval]
  | 0, k + 1 => by simp [peval]
  | n + 1, k + 1 => by
    rw [peval_gauss_succ_succ, gauss_at_one n k, gauss_at_one n (k + 1), Nat.choose_succ_succ]
    simp

/-- The four points of the dial on `[4, 2]_q = 1 + q + 2 * q ^ 2 + q ^ 3 + q ^ 4`: `q = 1` counts
`6 = Nat.choose 4 2`; `q = 2` counts the `35` planes in a `4`-dimensional space over the field
with two elements; `q = 0` gives `1`; `q = -1` gives `2 = Nat.choose 2 1` (Pascal folded in
half). -/
theorem gauss_4_2 :
    trim (gauss 4 2) = [1, 1, 2, 1, 1] ∧ peval 1 (gauss 4 2) = 6 ∧ peval 2 (gauss 4 2) = 35 ∧
    peval 0 (gauss 4 2) = 1 ∧ peval (-1) (gauss 4 2) = 2 := by
  decide +kernel

/-- At `q = -1`, `[n, k]` is `(n / 2).choose (k / 2)`, or `0` when `n` is even and `k` is odd.
Checked for `n ≤ 8`. -/
theorem gauss_at_neg_one :
    ∀ n : Fin 9, ∀ k : Fin 9, k.val ≤ n.val →
      peval (-1) (gauss n.val k.val) =
        (if n.val % 2 = 0 ∧ k.val % 2 = 1 then 0 else ((n.val / 2).choose (k.val / 2) : ℤ)) := by
  decide +kernel

/-! ### 2. The `q`-binomial theorem is a mult-add of a list -/

/-- Sum of two polynomials in `t` whose coefficients are polynomials in `q`. -/
def padd2 : List Poly → List Poly → List Poly
  | [], Q => Q
  | p :: P, [] => p :: P
  | p :: P, q :: Q => padd p q :: padd2 P Q

/-- Multiply by `1 + q ^ i * t`. -/
def mulLin (i : ℕ) (P : List Poly) : List Poly := padd2 P ([] :: P.map (pzeros i))

/-- `(1 + t) * (1 + q * t) * ⋯ * (1 + q ^ (n - 1) * t)`, which is
`1 + (t ✳ q * t ✳ ⋯ ✳ q ^ (n - 1) * t)`. -/
def qbinLHS (n : ℕ) : List Poly := (List.range n).foldl (fun P i ↦ mulLin i P) [[1]]

/-- `∑ k, q ^ (k * (k - 1) / 2) * [n, k]_q * t ^ k`. -/
def qbinRHS (n : ℕ) : List Poly :=
  (List.range (n + 1)).map fun k ↦ pzeros (k * (k - 1) / 2) (gauss n k)

/-- The `q`-binomial theorem, checked for `n ≤ 6`. -/
theorem q_binomial_theorem :
    ∀ n : Fin 7, (qbinLHS n.val).map trim = (qbinRHS n.val).map trim := by
  decide +kernel

/-! ### 3. The quantum group `U_q(sl₂)` -/

section Ring

variable {R : Type*} [CommRing R]

/-- On the irreducible representation of dimension `n + 1`, with basis `v₀, …, vₙ`, `E` and `F`
move between neighbours with factors `[k]_q` and `[n - k]_q` (symmetric `q`-integers
`[m] = (q ^ m - q⁻¹ ^ m) / (q - q⁻¹)`), and `K` acts by `q ^ (n - 2 * k)`. The relation
`[E, F] = (K - K⁻¹) / (q - q⁻¹)` on `v_k` is `[k + 1] * [n - k] - [k] * [n - k + 1] = [n - 2 * k]`.
Clearing the denominators `(q - q⁻¹) ^ 2` and writing `x = q ^ k`, `y = q ^ (n - k)` (`xi`, `yi`,
`qi` their inverses), it is this ring identity, so the relation holds on every irreducible
representation. -/
theorem uq_sl2_commutator (q qi x xi y yi : R) :
    (q * x - qi * xi) * (y - yi) - (x - xi) * (q * y - qi * yi) =
      (q - qi) * (xi * y - x * yi) := by
  ring

/-- The classical limit `q → 1`: `[E, F] = H`, that is `(a + 1) * b - a * (b + 1) = b - a`. -/
theorem sl2_classical (a b : R) : (a + 1) * b - a * (b + 1) = b - a := by ring

/-- The quantum dimension of the `2`-dimensional representation is `[2] = q + q⁻¹`; with
`q = A ^ 2` it is minus the knot loop value `d = -A ^ 2 - A⁻¹ ^ 2`. -/
theorem quantum_dim_is_loop (A Ai : R) : A ^ 2 + Ai ^ 2 = -(-A ^ 2 - Ai ^ 2) := by ring

/-! ### 4. The Ising ring -/

/-- One more site: multiply by the transfer matrix `T = !![1, u; u, 1]` (`u = exp (-2 * β * J)`,
the cost of a domain wall). `T ^ n = !![a, b; b, a]`, and we track `(a, b)`. -/
def isingStep (u : R) (p : R × R) : R × R := (p.1 + u * p.2, u * p.1 + p.2)

/-- The first row `(a, b)` of the `n`-th power of the transfer matrix. -/
def ising (u : R) (n : ℕ) : R × R := (isingStep u)^[n] (1, 0)

@[simp]
theorem ising_zero (u : R) : ising u 0 = (1, 0) := rfl

theorem ising_succ (u : R) (n : ℕ) : ising u (n + 1) = isingStep u (ising u n) :=
  Function.iterate_succ_apply' _ _ _

/-- `ising u n` is the first row of the `n`-th power of the transfer matrix. -/
theorem transfer_pow (u : R) (n : ℕ) :
    !![1, u; u, 1] ^ n =
      !![(ising u n).1, (ising u n).2; (ising u n).2, (ising u n).1] := by
  induction n with
  | zero => simp [Matrix.one_fin_two]
  | succ n ih =>
    rw [pow_succ, ih, ising_succ, isingStep, Matrix.mul_fin_two]
    congr 1; ext i j
    fin_cases i <;> fin_cases j <;> simp <;> ring

/-- The two eigen-directions: `a + b = (1 + u) ^ n` and `a - b = (1 - u) ^ n`. -/
theorem ising_eigen (u : R) (n : ℕ) :
    (ising u n).1 + (ising u n).2 = (1 + u) ^ n ∧ (ising u n).1 - (ising u n).2 = (1 - u) ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
    obtain ⟨h₁, h₂⟩ := ih
    rw [ising_succ, isingStep, pow_succ, pow_succ, ← h₁, ← h₂]
    constructor <;> ring

/-- The partition function of `n` spins on a ring (up to the factor `exp (n * β * J)`) is the
trace: `Z = (1 + u) ^ n + (1 - u) ^ n`. -/
theorem ising_trace (u : R) (n : ℕ) : 2 * (ising u n).1 = (1 + u) ^ n + (1 - u) ^ n := by
  obtain ⟨h₁, h₂⟩ := ising_eigen u n
  rw [← h₁, ← h₂]; ring

/-- The trace of the `n`-th power of the transfer matrix. -/
theorem trace_transfer_pow (u : R) (n : ℕ) :
    (!![1, u; u, 1] ^ n).trace = (1 + u) ^ n + (1 - u) ^ n := by
  rw [transfer_pow, Matrix.trace_fin_two_of, ← two_mul, ising_trace]

/-- In mult-add terms: `Z = 2 + mpow u n + mpow (-u) n`. Only the even layers of the cube
survive (a ring has an even number of domain walls). -/
theorem ising_multAdd (u : R) (n : ℕ) : 2 * (ising u n).1 = 2 + mpow u n + mpow (-u) n := by
  rw [ising_trace, mpow_eq, mpow_eq, ← sub_eq_add_neg]; ring

end Ring

end QCounting
