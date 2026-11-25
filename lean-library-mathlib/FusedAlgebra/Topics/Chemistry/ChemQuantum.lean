module

public import FusedAlgebra.Basic

/-!
# Quantum chemistry: independent fragments compose by mult-add

Write a molecule's wavefunction in intermediate normalisation, `Ψ = Φ * (1 + c)`, where `c` is
everything beyond the reference determinant. Two molecules that do not interact have
`Ψ_AB = Ψ_A * Ψ_B`, so `1 + c_AB = (1 + c_A) * (1 + c_B)`, that is `c_AB = c_A ✳ c_B`
(`fragments_multAdd`).

The corner `c_A * c_B` is "both molecules excited at once". Truncated CI (CISD) keeps only
`c_A + c_B`: it is the `k = 0` law (`MultAdd.wAdd_weight_zero`), and its size-consistency error
is exactly the dropped corner (`cisd_drops_corner`). Coupled cluster writes `1 + c = exp T`, the
log of the lift, so `T` adds and the corner is kept automatically. `N` equal fragments are
`MultAdd.one_add_mpow` and `MultAdd.mpow_add`.

The `N`-molecule model: each molecule has a reference and a doubly excited state `2 * Δ` above
it, coupled by `g`. CISD's energy `E` solves `E * (E - 2 * Δ) = N * g ^ 2` (`Secular`); the exact
energy is `N * e₁`. `cisd_gap` and `cisd_not_extensive` say CISD is not size-extensive.

Spin coupling counts states by mult-add: `spin_multiplicity`, and the Clebsch–Gordan series of
two equal spins, `clebsch_gordan_square` (the sum of the first odd numbers is a square).

This is textbook quantum chemistry (size consistency of CI and coupled cluster) read as
statements about mult-add. The core-Lean original stated the model over `ℚ`; here it holds in any
commutative ring without zero divisors.
-/

@[expose] public section

namespace ChemQuantum

open MultAdd Finset

variable {R : Type*}

section CommRing

variable [CommRing R]

/-- Two independent fragments: the excess amplitude mult-adds. This is `one_add_multAdd`. -/
theorem fragments_multAdd (cA cB : R) : (1 + cA) * (1 + cB) = 1 + (cA ✳ cB) :=
  (one_add_multAdd cA cB).symm

/-- What truncated CI drops when two fragments are joined: exactly the corner. -/
theorem cisd_drops_corner (cA cB : R) : (cA ✳ cB) - (cA + cB) = cA * cB := by
  unfold multAdd; ring

/-- Three fragments: the corner expansion has pairs and one triple, nothing else. -/
theorem three_fragments (a b c : R) :
    a ✳ b ✳ c = a + b + c + (a * b + a * c + b * c) + a * b * c := by
  unfold multAdd; ring

/-- What CISD drops from three fragments: every pair and the triple. -/
theorem cisd_drops_three (a b c : R) :
    (a ✳ b ✳ c) - (a + b + c) = a * b + a * c + b * c + a * b * c := by
  unfold multAdd; ring

/-! ### The `N`-molecule CISD model -/

/-- The CISD secular equation of `N` molecules: `E * (E - 2 * Δ) = N * g ^ 2`. -/
def Secular (D g E N : R) : Prop := E * (E - 2 * D) = N * g ^ 2

/-- If `e₁` solves the one-molecule problem, `N * e₁` misses the `N`-molecule CISD equation by
exactly `N * (N - 1) * e₁ ^ 2`, the pairs of molecules both excited. -/
theorem cisd_gap (D g e N : R) (h1 : Secular D g e 1) :
    (N * e) * (N * e - 2 * D) - N * g ^ 2 = N * (N - 1) * e ^ 2 := by
  unfold Secular at h1
  linear_combination N * h1

/-- So `N * e₁` solves `N`-molecule CISD only when `N = 0`, `N = 1` or `e₁ = 0`: CISD is not
size-extensive. -/
theorem cisd_not_extensive [NoZeroDivisors R] (D g e N : R) (hN0 : N ≠ 0) (hN1 : N ≠ 1)
    (he : e ≠ 0) (h1 : Secular D g e 1) : ¬ Secular D g (N * e) N := by
  intro hN
  have hz : N * (N - 1) * e ^ 2 = 0 := by
    rw [← cisd_gap D g e N h1, sub_eq_zero]
    exact hN
  exact mul_ne_zero (mul_ne_zero hN0 (sub_ne_zero.2 hN1)) (pow_ne_zero 2 he) hz

end CommRing

/-! ### Spin coupling counts states by mult-add -/

/-- Multiplicities multiply, so `2 * s` (multiplicity minus one) mult-adds. This is
`multAdd_add_one`. -/
theorem spin_multiplicity [Semiring R] (s1 s2 : R) : (s1 + 1) * (s2 + 1) = (s1 ✳ s2) + 1 :=
  (multAdd_add_one s1 s2).symm

/-- The Clebsch–Gordan series for two equal spins `j ⊗ j` (`2 * j = n`): the multiplicities
`2 * S + 1` for `S = 0, …, n` sum to `(n + 1) ^ 2`. The sum of the first odd numbers is a
square. -/
theorem clebsch_gordan_square (n : ℕ) : ∑ S ∈ range (n + 1), (2 * S + 1) = (n + 1) ^ 2 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [sum_range_succ, ih]; ring

end ChemQuantum
