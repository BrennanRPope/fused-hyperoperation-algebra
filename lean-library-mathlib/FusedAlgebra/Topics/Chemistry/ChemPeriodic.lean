module

public import FusedAlgebra.Basic

/-!
# The periodic law: shells are squares, periods are diagonals

Hydrogen's shell `n` has subshells `l = 0, …, n - 1` with `2 * l + 1` orbitals each. They are the
L-shaped gnomons of an `n × n` square, so the shell holds `n ^ 2` orbitals (`gnomons_fill`), and
`n ^ 2 = 1 + (n - 1) ✳ (n - 1)` (`shell_square`).

Madelung fills subshells along the diagonals `n + l = s`. Diagonal `s` holds `2 * ⌈s / 2⌉ ^ 2`
electrons (`diagCap`, `diagCap_closed`): 2, 2, 8, 8, 18, 18, 32, 32, which are Janet's left-step
periods (`diagonals`); the usual table merges the first two (`periodLen`, `periods`), giving
noble gases at 2, 10, 18, 36, 54, 86, 118 (`nobleZ`, `noble_gases`).

The dial: the `q`-integer `[n]_q` is `n` copies of `1` under weighted mult-add with `k = q - 1`
(`q_integer`, `q_integer_classical`).

The sum of the first odd numbers being a square is classical; the rest is counting.
-/

@[expose] public section

namespace ChemPeriodic

open MultAdd Finset

/-- The shell square: `(m + 1) ^ 2 = 1 + m ✳ m`, with `m = n - 1`. -/
theorem shell_square {R : Type*} [CommSemiring R] (m : R) : (m + 1) ^ 2 = 1 + (m ✳ m) := by
  rw [one_add_multAdd]; ring

/-- Orbitals in subshells `l = 0, …, m`: the sum of `2 * l + 1`. -/
def shellOrbitals (m : ℕ) : ℕ := ∑ l ∈ range (m + 1), (2 * l + 1)

/-- The gnomons fill the square: shell `n = m + 1` holds `n ^ 2` orbitals. -/
theorem gnomons_fill (m : ℕ) : shellOrbitals m = (m + 1) ^ 2 := by
  unfold shellOrbitals
  induction m with
  | zero => rfl
  | succ m ih => rw [sum_range_succ, ih]; ring

/-- Each gnomon is corner plus two edges of the next square: `2 * l + 1 = 1 + l + l`. -/
theorem gnomon (l : ℕ) : (l + 1) ^ 2 - l ^ 2 = 1 + l + l :=
  Nat.sub_eq_of_eq_add (by ring)

/-- With spin, shell `n` holds `2 * n ^ 2` electrons. -/
theorem shell_capacity (m : ℕ) : 2 * shellOrbitals m = 2 * (m + 1) ^ 2 := by
  rw [gnomons_fill]

/-- SO(4): two independent spins with multiplicity `n` each give `n * n` states, and `n - 1`
mult-adds with itself. This is `one_add_multAdd`. -/
theorem two_spins {R : Type*} [CommSemiring R] (a b : R) : (a + 1) * (b + 1) = 1 + (a ✳ b) := by
  rw [one_add_multAdd]; ring

/-! ### Madelung diagonals and the periods -/

/-- Electrons on the Madelung diagonal `n + l = s`: subshells `l = 0, …, ⌊(s - 1) / 2⌋`, each
holding `2 * (2 * l + 1)`. -/
def diagCap (s : ℕ) : ℕ := 2 * shellOrbitals ((s - 1) / 2)

/-- Diagonal `s` holds `2 * ⌈s / 2⌉ ^ 2` electrons, in closed form. -/
theorem diagCap_closed (s : ℕ) : diagCap s = 2 * ((s - 1) / 2 + 1) ^ 2 := by
  unfold diagCap; rw [gnomons_fill]

theorem diagonals :
    (List.range 8).map (fun i ↦ diagCap (i + 1)) = [2, 2, 8, 8, 18, 18, 32, 32] := by
  simp only [diagCap_closed]; decide

/-- The usual table merges the first two diagonals into period 1. -/
def periodLen (p : ℕ) : ℕ := 2 * (p / 2 + 1) ^ 2

theorem periods : (List.range 7).map (fun i ↦ periodLen (i + 1)) = [2, 8, 8, 18, 18, 32, 32] := by
  decide

/-- Noble gases: running totals of the period lengths. -/
def nobleZ (p : ℕ) : ℕ := ∑ i ∈ range p, periodLen (i + 1)

theorem noble_gases :
    (List.range 7).map (fun i ↦ nobleZ (i + 1)) = [2, 10, 18, 36, 54, 86, 118] := by
  decide

/-- Period `p` (`≥ 2`) and the diagonal pair it covers agree: periods `2 * k` and `2 * k + 1`
have the length of diagonals `2 * k + 1` and `2 * k + 2`. -/
theorem period_is_diagonal (k : ℕ) : periodLen (2 * k + 1) = diagCap (2 * k + 2) := by
  rw [diagCap_closed, periodLen]
  have h1 : (2 * k + 1) / 2 = k := by omega
  have h2 : (2 * k + 2 - 1) / 2 = k := by omega
  rw [h1, h2]

/-! ### The dial: `q`-integers -/

/-- `n` copies of `1` under `wAdd k` lift to `(1 + k) ^ n`: with `q = k + 1`, that is `q ^ n`, so
the weighted power is the `q`-integer `(q ^ n - 1) / (q - 1)`. -/
theorem q_integer {R : Type*} [CommSemiring R] (k : R) (n : ℕ) :
    1 + k * wpow k 1 n = (1 + k) ^ n := by
  simpa using one_add_mul_wpow k 1 n

/-- At `k = 0` (`q = 1`) the weighted power of `1` is `n` itself: the ordinary shell. -/
theorem q_integer_classical {R : Type*} [CommSemiring R] (n : ℕ) : wpow (0 : R) 1 n = n := by
  induction n with
  | zero => simp
  | succ n ih => simp [wpow_succ, ih]

end ChemPeriodic
