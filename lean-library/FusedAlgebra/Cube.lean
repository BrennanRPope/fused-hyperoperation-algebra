/-
  FusedAlgebra.Cube: the layers of the mult-add cube
  ==================================================

  1 + (a ◯ a ◯ … ◯ a) = (1 + a)^n, and the n-cube behind it has
  `layer n k` corners with k switches on (Pascal's triangle). Shared by the
  e, topology (Lucas, Fermat, AKS, Kummer, Wilson) modules.

  From topology/lean/CubeProofs.lean and e/EFromMultAdd.lean, which each
  built the same layers. Core Lean only, no `sorry`.
-/
import FusedAlgebra.Base

namespace CircleOp.Cube

/-! ## The cube's layers -/

/-- `layer n k`: how many corners of the n-cube have exactly k switches on
    (grade k). Built by Pascal's rule: adding one more axis to the cube
    copies it twice, once with the new switch off and once with it on. -/
def layer : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 1
  | n + 1, k + 1 => layer n k + layer n (k + 1)

theorem layer_zero_right (n : Nat) : layer n 0 = 1 := by cases n <;> rfl

theorem layer_succ_succ (n k : Nat) : layer (n + 1) (k + 1) = layer n k + layer n (k + 1) := rfl

theorem layer_eq_zero : ∀ n k, n < k → layer n k = 0
  | 0, 0, h => absurd h (by decide)
  | 0, _ + 1, _ => rfl
  | _ + 1, 0, h => absurd h (by omega)
  | n + 1, k + 1, h => by
    rw [layer_succ_succ, layer_eq_zero n k (by omega), layer_eq_zero n (k + 1) (by omega)]

theorem layer_self : ∀ n, layer n n = 1
  | 0 => rfl
  | n + 1 => by rw [layer_succ_succ, layer_self n, layer_eq_zero n (n + 1) (by omega)]

/-- Absorption: (k+1)·C(n+1,k+1) = (n+1)·C(n,k). -/
theorem absorb : ∀ n k, (k + 1) * layer (n + 1) (k + 1) = (n + 1) * layer n k := by
  intro n
  induction n with
  | zero =>
    intro k
    cases k with
    | zero => rfl
    | succ k => simp [layer]
  | succ n ih =>
    intro k
    cases k with
    | zero =>
      have h := ih 0
      rw [layer_succ_succ]
      simp only [layer_zero_right] at h ⊢
      omega
    | succ j =>
      have h1 := ih j
      have h2 := ih (j + 1)
      have hp := layer_succ_succ n j
      rw [layer_succ_succ (n + 1) (j + 1)]
      generalize layer (n + 1) (j + 1) = X at *
      generalize layer (n + 1) (j + 1 + 1) = Y at *
      generalize layer n j = A at *
      generalize layer n (j + 1) = B at *
      subst hp
      grind

/-! ## Primes -/

/-- Factorial. -/
def fact : Nat → Nat
  | 0 => 1
  | k + 1 => (k + 1) * fact k

theorem fact_pos : ∀ k, 0 < fact k
  | 0 => by decide
  | k + 1 => Nat.mul_pos (by omega) (fact_pos k)


/-- (k+1)·C(n,k+1) = (n-k)·C(n,k). -/
theorem layer_step (n k : Nat) : (k + 1) * layer n (k + 1) = (n - k) * layer n k := by
  cases n with
  | zero =>
    rw [layer_eq_zero 0 (k + 1) (by omega)]; simp
  | succ m =>
    cases Nat.lt_or_ge m k with
    | inl h =>
      rw [layer_eq_zero (m + 1) (k + 1) (by omega)]
      have h0 : m + 1 - k = 0 := by omega
      rw [h0]; simp
    | inr h =>
      have ha := absorb m k
      have hp := layer_succ_succ m k
      cases k with
      | zero =>
        simp only [layer_zero_right, Nat.sub_zero] at *
        rw [hp] at ha
        generalize layer (m + 1) 1 = X at *
        generalize layer m 1 = B at *
        grind
      | succ j =>
        -- (j+2) C(m+1, j+2) = (m - j) C(m+1, j+1)
        have ha2 := absorb m (j + 1)
        have hp2 := layer_succ_succ m (j + 1)
        rw [hp2] at ha2
        have hpj := layer_succ_succ m j
        have hsub : m + 1 - (j + 1) = m - j := by omega
        rw [hsub, hpj]
        have hab := absorb m j
        rw [hpj] at hab
        generalize layer m (j + 1 + 1) = C at *
        generalize layer m (j + 1) = B at *
        generalize layer m j = A at *
        have hmj : m = (m - j) + j := by omega
        generalize m - j = d at *
        subst hmj
        grind

end CircleOp.Cube
