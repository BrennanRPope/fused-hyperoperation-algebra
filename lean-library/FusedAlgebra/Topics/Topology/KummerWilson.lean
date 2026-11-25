/-
  Kummer and Wilson on the mult-add cube
  ======================================

  Part of Fused Mult-Add Algebra. Builds on CubeProofs.lean (the cube's
  layer counts, primes, Fermat). Pictures: ../fig_kummer_wilson.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`, only Lean's standard
  axioms. Nothing is taken as given. To run it:
      lean -o CubeProofs.olean CubeProofs.lean
      LEAN_PATH=.:$(lean --print-prefix)/lib/lean lean KummerWilson.lean

  Proved:
  * carries_dsum: adding digit by digit in base p, every carry turns p
                  units of digit sum into 1:
                  (p−1)·carries + s(a+b) = s(a) + s(b).
  * dsum_succ:    adding 1 trades each trailing digit p−1 for a carry.
  * layer_val_dsum: (p−1)·v_p(C(j+k,k)) = s(j) + s(k) − s(j+k), from the
                  absorption rule of the cube's layers.
  * kummer:       the exact power of p dividing C(j+k,k), the number of
                  corners of the (j+k)-cube with k switches on, is the
                  number of carries when you add j and k in base p
                  (Kummer 1852).
  * S_pow:        inclusion–exclusion over the cube's corners:
                  Σₖ (−1)ᵏ C(n,k)(c−k)ᵐ is n! for m = n and 0 for m < n.
  * wilson:       (p−1)! ≡ −1 (mod p) for every prime p (Lagrange 1771).
  * wilson_iff:   for n ≥ 2, (n−1)! ≡ −1 (mod n) exactly when n is prime.

  Both theorems are classical. What is new here is only the formal check
  in core Lean and the reading on the cube and the mult-add circle.

  LIBRARY NOTE: imports the topology CubeProofs module; fact now comes from FusedAlgebra.Cube.
-/
import FusedAlgebra.Topics.Topology.CubeProofs

namespace CircleOp.Cube

/-! ## Kummer: carries -/

/-- `val p x`: how many times p divides x (0 when x = 0). -/

def val (p x : Nat) : Nat :=
  if h : 2 ≤ p ∧ 0 < x ∧ x % p = 0 then 1 + val p (x / p) else 0
termination_by x
decreasing_by exact Nat.div_lt_self h.2.1 (by omega)

set_option linter.unusedVariables false in
/-- `dsum p x`: the sum of x's digits in base p. -/
def dsum (p x : Nat) : Nat :=
  if h : 2 ≤ p ∧ 0 < x then x % p + dsum p (x / p) else 0
termination_by x
decreasing_by exact Nat.div_lt_self h.2 (by omega)

set_option linter.unusedVariables false in
/-- `carries p a b c`: how many carries happen when you add a and b in
    base p, starting with carry c. -/
def carries (p a b c : Nat) : Nat :=
  if h : 2 ≤ p ∧ 0 < a + b then
    (a % p + b % p + c) / p + carries p (a / p) (b / p) ((a % p + b % p + c) / p)
  else 0
termination_by a + b
decreasing_by
  have h1 := Nat.div_le_self a p
  have h2 := Nat.div_le_self b p
  rcases Nat.eq_zero_or_pos a with ha | ha
  · have := Nat.div_lt_self (n := b) (by omega) (by omega : 1 < p); omega
  · have := Nat.div_lt_self ha (by omega : 1 < p); omega

theorem dsum_eq {p : Nat} (hp : 2 ≤ p) (x : Nat) : dsum p x = x % p + dsum p (x / p) := by
  rcases Nat.eq_zero_or_pos x with rfl | hx
  · rw [dsum]; simp [dsum]
  · rw [dsum]; simp [hp, hx]

theorem dsum_zero (p : Nat) : dsum p 0 = 0 := by rw [dsum]; simp

theorem val_of_mod_ne {p x : Nat} (h : x % p ≠ 0) : val p x = 0 := by
  rw [val]; simp [h]

theorem val_of_mod_eq {p x : Nat} (hp : 2 ≤ p) (hx : 0 < x) (h : x % p = 0) :
    val p x = 1 + val p (x / p) := by
  rw [val]; simp [hp, hx, h]


theorem carries_eq {p : Nat} (hp : 2 ≤ p) (a b c : Nat) (hab : 0 < a + b) :
    carries p a b c = (a % p + b % p + c) / p +
      carries p (a / p) (b / p) ((a % p + b % p + c) / p) := by
  rw [carries]; simp [hp, hab]

theorem carries_zero {p : Nat} (c : Nat) : carries p 0 0 c = 0 := by
  rw [carries]; simp

/-- Adding digit by digit: each carry turns p units of digit sum into 1. -/
theorem carries_dsum {p : Nat} (hp : 2 ≤ p) : ∀ s a b c, a + b = s → c ≤ 1 →
    (p - 1) * carries p a b c + dsum p (a + b + c) = dsum p a + dsum p b + c := by
  intro s
  induction s using Nat.strongRecOn with
  | _ s ih =>
  intro a b c hs hc
  rcases Nat.eq_zero_or_pos (a + b) with h0 | hpos
  · have ha : a = 0 := by omega
    have hb : b = 0 := by omega
    subst ha; subst hb
    simp only [Nat.zero_add]
    rw [carries_zero, dsum_eq hp c, dsum_zero]
    have hc' : c % p = c := Nat.mod_eq_of_lt (by omega)
    have hc'' : c / p = 0 := Nat.div_eq_of_lt (by omega)
    rw [hc', hc'', dsum_zero]; omega
  · rw [carries_eq hp a b c hpos]
    have ra := Nat.mod_add_div a p
    have rb := Nat.mod_add_div b p
    have ma := Nat.mod_lt a (by omega : p > 0)
    have mb := Nat.mod_lt b (by omega : p > 0)
    generalize hr : a % p + b % p + c = r
    have rr := Nat.mod_add_div r p
    have mr := Nat.mod_lt r (by omega : p > 0)
    have hc1 : r / p ≤ 1 := by
      have : r < p * 2 := by omega
      exact Nat.le_of_lt_succ ((Nat.div_lt_iff_lt_mul (by omega)).2 (by grind))
    have hsum : a + b + c = p * (a / p + b / p + r / p) + r % p := by
      rw [Nat.mul_add, Nat.mul_add]; omega
    have dm := divmod_of_eq (q := a / p + b / p + r / p) mr
    rw [← hsum] at dm
    have hlt : a / p + b / p < s := by
      have h1 := Nat.div_le_self a p
      have h2 := Nat.div_le_self b p
      rcases Nat.eq_zero_or_pos a with ha | ha
      · have := Nat.div_lt_self (n := b) (by omega) (by omega : 1 < p); omega
      · have := Nat.div_lt_self ha (by omega : 1 < p); omega
    have IH := ih _ hlt (a / p) (b / p) (r / p) rfl hc1
    rw [dsum_eq hp (a + b + c), dm.1, dm.2, dsum_eq hp a, dsum_eq hp b]
    have e1 : (p - 1) * (r / p + carries p (a / p) (b / p) (r / p)) =
        (p - 1) * (r / p) + (p - 1) * carries p (a / p) (b / p) (r / p) := Nat.mul_add _ _ _
    have e2 : (p - 1) * (r / p) + r / p = p * (r / p) := by
      rw [Nat.sub_one_mul]; have : r / p ≤ p * (r / p) := Nat.le_mul_of_pos_left _ (by omega); omega
    omega

/-- Adding 1: the digit sum goes up by 1, minus p−1 for each trailing digit p−1. -/
theorem dsum_succ {p : Nat} (hp : 2 ≤ p) : ∀ x, dsum p (x + 1) + (p - 1) * val p (x + 1) = dsum p x + 1 := by
  intro x
  induction x using Nat.strongRecOn with
  | _ x ih =>
  have rx := Nat.mod_add_div x p
  have mx := Nat.mod_lt x (by omega : p > 0)
  rcases Nat.lt_or_ge (x % p + 1) p with h | h
  · have e : x + 1 = p * (x / p) + (x % p + 1) := by omega
    have dm := divmod_of_eq (q := x / p) h
    rw [← e] at dm
    rw [val_of_mod_ne (by omega), dsum_eq hp (x + 1), dm.1, dm.2, dsum_eq hp x]
    omega
  · have e : x + 1 = p * (x / p + 1) + 0 := by rw [Nat.mul_add]; omega
    have dm := divmod_of_eq (q := x / p + 1) (by omega : 0 < p)
    rw [← e] at dm
    have hx : 0 < x := by
      rcases Nat.eq_zero_or_pos x with h0 | h0
      · subst h0; simp at h; omega
      · exact h0
    have IH := ih (x / p) (Nat.div_lt_self hx (by omega))
    rw [val_of_mod_eq hp (by omega) dm.2, dm.1, dsum_eq hp (x + 1), dm.1, dm.2, dsum_eq hp x]
    have e1 : (p - 1) * (1 + val p (x / p + 1)) = (p - 1) + (p - 1) * val p (x / p + 1) := by
      rw [Nat.mul_add, Nat.mul_one]
    omega

/-! ### Exact powers of p -/

/-- `Val p x e`: p^e divides x, and p^(e+1) does not. -/
def Val (p x e : Nat) : Prop := p ^ e ∣ x ∧ ¬ p ^ (e + 1) ∣ x

theorem val_spec {p : Nat} (hp : 2 ≤ p) : ∀ x, 0 < x → Val p x (val p x) := by
  intro x
  induction x using Nat.strongRecOn with
  | _ x ih =>
  intro hx
  by_cases h : x % p = 0
  · rw [val_of_mod_eq hp hx h]
    have rx := Nat.mod_add_div x p
    have hq : 0 < x / p := by
      rcases Nat.eq_zero_or_pos (x / p) with h0 | h0
      · rw [h0] at rx; omega
      · exact h0
    have IH := ih (x / p) (Nat.div_lt_self hx (by omega)) hq
    have ex : x = p * (x / p) := by omega
    have pp : 0 < p := by omega
    generalize x / p = q at IH ex ⊢
    subst ex
    constructor
    · rw [ Nat.add_comm, Nat.pow_succ, Nat.mul_comm (p ^ _) p]
      exact Nat.mul_dvd_mul_left p IH.1
    · rw [Nat.add_comm 1, Nat.pow_succ, Nat.mul_comm (p ^ _) p, Nat.mul_dvd_mul_iff_left pp]
      exact IH.2
  · rw [val_of_mod_ne h]
    refine ⟨by simp, ?_⟩
    simp only [Nat.zero_add, Nat.pow_one]
    intro hd; exact h (Nat.mod_eq_zero_of_dvd hd)

theorem Val.unique {p x e f : Nat} (he : Val p x e) (hf : Val p x f) : e = f := by
  rcases Nat.lt_trichotomy e f with h | h | h
  · exact absurd (Nat.dvd_trans (Nat.pow_dvd_pow p (by omega)) hf.1) he.2
  · exact h
  · exact absurd (Nat.dvd_trans (Nat.pow_dvd_pow p (by omega)) he.1) hf.2

theorem euclid {p a b : Nat} (hp : IsPrime p) (ha : ¬ p ∣ a) (h : p ∣ a * b) : p ∣ b := by
  have hc : Nat.Coprime p a := by
    rcases hp.2 _ (Nat.gcd_dvd_left p a) with g | g
    · exact g
    · exact absurd (g ▸ Nat.gcd_dvd_right p a) ha
  exact hc.dvd_of_dvd_mul_left h

theorem Val.mul {p a b e f : Nat} (hp : IsPrime p) (ha : Val p a e) (hb : Val p b f) :
    Val p (a * b) (e + f) := by
  have pp : 0 < p := by have := hp.1; omega
  obtain ⟨a', rfl⟩ := ha.1
  obtain ⟨b', rfl⟩ := hb.1
  have na : ¬ p ∣ a' := by
    intro hd; apply ha.2; rw [Nat.pow_succ]; exact Nat.mul_dvd_mul_left _ hd
  have nb : ¬ p ∣ b' := by
    intro hd; apply hb.2; rw [Nat.pow_succ]; exact Nat.mul_dvd_mul_left _ hd
  have e1 : p ^ e * a' * (p ^ f * b') = p ^ (e + f) * (a' * b') := by
    rw [Nat.pow_add]; grind
  rw [e1]
  refine ⟨Nat.dvd_mul_right _ _, ?_⟩
  rw [Nat.pow_succ, Nat.mul_dvd_mul_iff_left (Nat.pow_pos pp)]
  intro hd
  exact nb (euclid hp na hd)

theorem layer_pos : ∀ n k, k ≤ n → 0 < layer n k
  | n, 0, _ => by rw [layer_zero_right]; decide
  | 0, _ + 1, h => absurd h (by omega)
  | n + 1, k + 1, h => by
    rw [layer_succ_succ]; have := layer_pos n k (by omega); omega

/-- The cube count and the digit sums: (p−1)·v_p(C(j+k,k)) = s(j) + s(k) − s(j+k). -/
theorem layer_val_dsum {p : Nat} (hp : IsPrime p) : ∀ k j,
    (p - 1) * val p (layer (j + k) k) + dsum p (j + k) = dsum p j + dsum p k := by
  have p2 := hp.1
  intro k
  induction k with
  | zero =>
    intro j
    rw [Nat.add_zero, layer_zero_right, val_of_mod_ne (by rw [Nat.mod_eq_of_lt (by omega)]; decide),
      dsum_zero]
    simp
  | succ k ih =>
    intro j
    have IH := ih j
    have ab := absorb (j + k) k
    have L := val_spec p2 (layer (j + k + 1) (k + 1)) (layer_pos _ _ (by omega))
    have R := val_spec p2 (layer (j + k) k) (layer_pos _ _ (by omega))
    have v1 := val_spec p2 (j + k + 1) (by omega)
    have v2 := val_spec p2 (k + 1) (by omega)
    have hL := Val.mul hp v2 L
    have hR := Val.mul hp v1 R
    rw [ab] at hL
    have hu := Val.unique hL hR
    have hu' := congrArg (fun t => (p - 1) * t) hu
    simp only [Nat.mul_add] at hu'
    have s1 := dsum_succ p2 (j + k)
    have s2 := dsum_succ p2 k
    show (p - 1) * val p (layer (j + k + 1) (k + 1)) + dsum p (j + k + 1) = dsum p j + dsum p (k + 1)
    omega

/-- **Kummer's theorem (1852).** The exact power of p that divides
    C(j+k, k), the number of corners of the (j+k)-cube with k switches on,
    is the number of carries when you add j and k in base p. -/
theorem kummer {p : Nat} (hp : IsPrime p) (j k : Nat) :
    p ^ carries p j k 0 ∣ layer (j + k) k ∧ ¬ p ^ (carries p j k 0 + 1) ∣ layer (j + k) k := by
  have p2 := hp.1
  have h1 := layer_val_dsum hp k j
  have h2 := carries_dsum p2 (j + k) j k 0 rfl (by omega)
  rw [Nat.add_zero] at h2
  have e : (p - 1) * val p (layer (j + k) k) = (p - 1) * carries p j k 0 := by omega
  have e' := Nat.eq_of_mul_eq_mul_left (by omega) e
  have := val_spec p2 (layer (j + k) k) (layer_pos _ _ (by omega))
  rw [e'] at this
  exact this


/-! ## Wilson: inclusion–exclusion on the cube -/

/-- Signed corner sum over the n-cube: a corner with k switches on counts
    with sign (−1)ᵏ and weight f k. -/
def S (n : Nat) (f : Nat → Int) : Int :=
  bsum (fun k => (-1 : Int) ^ k * (layer n k : Int) * f k) (n + 1)

theorem S_add (n : Nat) (f g : Nat → Int) : S n (fun k => f k + g k) = S n f + S n g := by
  unfold S; rw [← bsum_add]; congr 1; funext k; grind

theorem S_mul (n : Nat) (c : Int) (f : Nat → Int) : S n (fun k => c * f k) = c * S n f := by
  unfold S; rw [← bsum_mul]; congr 1; funext k; grind

theorem S_zero (n : Nat) : S n (fun _ => 0) = 0 := by
  have := S_mul n 0 (fun _ => 0); simp at this ⊢; exact this

theorem S_sub (n : Nat) (f g : Nat → Int) : S n (fun k => f k - g k) = S n f - S n g := by
  have h1 := S_add n f (fun k => (-1) * g k)
  rw [S_mul] at h1
  have e : (fun k => f k - g k) = (fun k => f k + (-1) * g k) := by funext k; grind
  rw [e, h1]; grind

theorem S_bsum (n : Nat) (g : Nat → Nat → Int) : ∀ m,
    S n (fun k => bsum (fun i => g i k) m) = bsum (fun i => S n (g i)) m
  | 0 => by simp only [bsum]; exact S_zero n
  | m + 1 => by simp only [bsum]; rw [S_add, S_bsum n g m]

/-- Pascal's rule for the signed sum: adding an axis takes a difference. -/
theorem S_succ (n : Nat) (f : Nat → Int) : S (n + 1) f = S n f - S n (fun k => f (k + 1)) := by
  unfold S
  rw [bsum_shift]
  have e : (fun k => (-1 : Int) ^ (k + 1) * ((layer (n + 1) (k + 1) : Nat) : Int) * f (k + 1)) =
      (fun k => (-1 : Int) ^ (k + 1) * (layer n (k + 1) : Int) * f (k + 1) +
        (-1) * ((-1 : Int) ^ k * (layer n k : Int) * f (k + 1))) := by
    funext k; rw [layer_succ_succ, Int.natCast_add, Int.pow_succ]; grind
  rw [e, bsum_add, bsum_mul]
  have h2 : bsum (fun k => (-1 : Int) ^ k * (layer n k : Int) * f k) (n + 1) =
      (-1 : Int) ^ 0 * (layer n 0 : Int) * f 0 +
        bsum (fun k => (-1 : Int) ^ (k + 1) * (layer n (k + 1) : Int) * f (k + 1)) (n + 1) := by
    rw [bsum_shift, bsum, layer_eq_zero n (n + 1) (by omega)]; simp
  rw [h2]; simp only [layer_zero_right]; omega

theorem layer_succ_self (n : Nat) : layer (n + 1) n = n + 1 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [layer_succ_succ, ih, layer_self]

/-- **Inclusion–exclusion on the cube.** Σₖ (−1)ᵏ C(n,k) (c − k)ᵐ is n!
    when m = n and 0 when m < n. With c = n it counts the orderings of n
    things: corner K of the n-cube is the set of values a map is made to
    miss, and the signed count of maps that miss them leaves exactly the
    maps that miss nothing. -/
theorem S_pow : ∀ n m, m ≤ n → ∀ c : Int,
    S n (fun k => (c - k) ^ m) = if m = n then (fact n : Int) else 0 := by
  intro n
  induction n with
  | zero =>
    intro m hm c
    have : m = 0 := by omega
    subst this; simp [S, bsum, layer, fact]
  | succ n ih =>
    intro m hm c
    rw [S_succ]
    have hd : ∀ k : Nat, (c - (k : Int)) ^ m - (c - ((k + 1 : Nat) : Int)) ^ m =
        bsum (fun i => (layer m i : Int) * ((c - 1) - k) ^ i) m := by
      intro k
      have hb := binomial m ((c - 1) - k)
      unfold cubeSum at hb
      rw [bsum, layer_self] at hb
      have e1 : 1 + ((c - 1) - k) = c - (k : Int) := by grind
      have e2 : c - ((k + 1 : Nat) : Int) = (c - 1) - k := by push_cast; grind
      rw [e1] at hb; rw [e2, hb]; simp
    have e3 : S n (fun k => (c - k) ^ m) - S n (fun k => (c - ((k + 1 : Nat) : Int)) ^ m) =
        S n (fun k => bsum (fun i => (layer m i : Int) * ((c - 1) - k) ^ i) m) := by
      rw [← S_sub]; congr 1; funext k; exact hd k
    rw [e3]
    have e4 := S_bsum n (fun i (k : Nat) => (layer m i : Int) * ((c - 1) - (k : Int)) ^ i) m
    rw [e4]
    have e5 : ∀ i, i < m → S n (fun k => (layer m i : Int) * ((c - 1) - (k : Int)) ^ i) =
        (layer m i : Int) * (if i = n then (fact n : Int) else 0) := by
      intro i hi
      rw [S_mul, ih i (by omega) (c - 1)]
    rw [bsum_congr _ _ m e5]
    rcases Nat.lt_or_ge m (n + 1) with hlt | hge
    · have : ∀ i, i < m → (layer m i : Int) * (if i = n then (fact n : Int) else 0) = 0 := by
        intro i hi; have : i ≠ n := by omega
        simp [this]
      rw [bsum_congr _ (fun _ => 0) m this]
      have hz : ∀ m, bsum (fun _ => (0 : Int)) m = 0 := by
        intro m; induction m with
        | zero => rfl
        | succ m ih => simp [bsum, ih]
      rw [hz]; have : m ≠ n + 1 := by omega
      simp [this]
    · have hm' : m = n + 1 := by omega
      subst hm'
      rw [bsum]
      have : ∀ i, i < n → (layer (n + 1) i : Int) * (if i = n then (fact n : Int) else 0) = 0 := by
        intro i hi; have : i ≠ n := by omega
        simp [this]
      rw [bsum_congr _ (fun _ => 0) n this]
      have hz : ∀ m, bsum (fun _ => (0 : Int)) m = 0 := by
        intro m; induction m with
        | zero => rfl
        | succ m ih => simp [bsum, ih]
      rw [hz, layer_succ_self]
      simp [fact]


/-! ### Wilson -/

theorem bsum_meq {p : Nat} (f g : Nat → Int) : ∀ m, (∀ k, k < m → MEq p (f k) (g k)) →
    MEq p (bsum f m) (bsum g m)
  | 0, _ => MEq.refl _ _
  | m + 1, h => by
    simp only [bsum]
    exact MEq.add (bsum_meq f g m (fun k hk => h k (by omega))) (h m (by omega))

theorem bsum_one : ∀ m, bsum (fun _ => (1 : Int)) m = m
  | 0 => rfl
  | m + 1 => by simp only [bsum]; rw [bsum_one m]; push_cast; rfl

/-- Along the (p−1)-cube, the layer counts alternate ±1 mod p. -/
theorem alt_layer {p : Nat} (hp : IsPrime p) : ∀ k, k ≤ p - 1 →
    MEq p ((-1 : Int) ^ k * (layer (p - 1) k : Int)) 1 := by
  obtain ⟨q, hq⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.1; omega⟩
  intro k
  induction k with
  | zero => intro _; simp [layer_zero_right]; exact MEq.refl _ _
  | succ k ih =>
    intro hk
    have IH := ih (by omega)
    have hd := prime_dvd_middle hp (k := k + 1) (by omega) (by omega)
    rw [hq, layer_succ_succ] at hd
    rw [hq, Nat.add_sub_cancel] at IH ⊢
    obtain ⟨c, hc⟩ := hd
    obtain ⟨d, hd⟩ := IH
    refine ⟨(-1 : Int) ^ (k + 1) * c + d, ?_⟩
    have e : (layer q (k + 1) : Int) = ((q + 1 : Nat) : Int) * c - layer q k := by
      have := congrArg (fun t : Nat => (t : Int)) hc; push_cast at this; omega
    rw [e, Int.pow_succ]
    grind

/-- Fermat for a unit: if p does not divide a, then a^(p−1) ≡ 1. -/
theorem fermat_unit {p : Nat} (hp : IsPrime p) (a : Nat) (ha : ¬ p ∣ a) :
    MEq p ((a : Int) ^ (p - 1)) 1 := by
  obtain ⟨q, hq⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.1; omega⟩
  have h := fermat hp a
  unfold MEq at h ⊢
  have e : (a : Int) ^ p - a = (a : Int) * ((a : Int) ^ (p - 1) - 1) := by
    rw [hq, Nat.add_sub_cancel, Int.pow_succ]; grind
  rw [e, Int.ofNat_dvd_left, Int.natAbs_mul, Int.natAbs_natCast] at h
  exact Int.ofNat_dvd_left.mpr (euclid hp ha h)

/-- **Wilson's theorem (Lagrange 1771).** (p − 1)! ≡ −1 (mod p) for every
    prime p. Lean's route: inclusion–exclusion on the (p−1)-cube counts
    the (p−1)! orderings; mod p every layer but the last adds exactly 1,
    so the count is p − 1 ≡ −1. -/
theorem wilson {p : Nat} (hp : IsPrime p) : MEq p (fact (p - 1)) (-1) := by
  obtain ⟨q, hq⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.1; omega⟩
  have hq1 : 1 ≤ q := by have := hp.1; omega
  subst hq
  have hS := S_pow q q (Nat.le_refl q) (q : Int)
  simp only [ite_true] at hS
  rw [Nat.add_sub_cancel, ← hS]
  unfold S
  rw [bsum]
  have hlast : (-1 : Int) ^ q * (layer q q : Int) * ((q : Int) - (q : Nat)) ^ q = 0 := by
    rw [Int.sub_self]
    obtain ⟨r, rfl⟩ : ∃ r, q = r + 1 := ⟨q - 1, by omega⟩
    simp [Int.pow_succ]
  rw [hlast, Int.add_zero]
  have hterms : ∀ k, k < q → MEq (q + 1) ((-1 : Int) ^ k * (layer q k : Int) * ((q : Int) - k) ^ q)
      ((fun _ => (1 : Int)) k) := by
    intro k hk
    have h1 := alt_layer hp k (by omega)
    rw [Nat.add_sub_cancel] at h1
    have e : ((q : Int) - (k : Int)) = ((q - k : Nat) : Int) := by
      rw [Int.ofNat_sub (by omega)]
    have h2 := fermat_unit hp (q - k) (by
      intro hd; have := Nat.le_of_dvd (by omega) hd; omega)
    rw [Nat.add_sub_cancel] at h2
    rw [e]
    have := MEq.mul h1 h2
    simpa using this
  have hb := bsum_meq _ _ q hterms
  rw [bsum_one] at hb
  refine MEq.trans hb ⟨1, ?_⟩
  push_cast; omega

/-- n! divides into every number from 1 to n. -/
theorem dvd_fact : ∀ n d, 1 ≤ d → d ≤ n → d ∣ fact n
  | 0, d, h1, h2 => absurd h2 (by omega)
  | n + 1, d, h1, h2 => by
    rcases Nat.lt_or_ge d (n + 1) with h | h
    · exact Nat.dvd_trans (dvd_fact n d h1 (by omega)) ⟨n + 1, by simp [fact, Nat.mul_comm]⟩
    · have : d = n + 1 := by omega
      subst this; exact ⟨fact n, rfl⟩

/-- **Wilson's test.** For n ≥ 2, (n − 1)! ≡ −1 (mod n) exactly when n is prime. -/
theorem wilson_iff (n : Nat) (hn : 2 ≤ n) : MEq n (fact (n - 1)) (-1) ↔ IsPrime n := by
  constructor
  · intro h
    refine ⟨hn, fun d hd => ?_⟩
    rcases Nat.eq_zero_or_pos d with h0 | h0
    · subst h0; rw [Nat.zero_dvd] at hd; omega
    have hdn := Nat.le_of_dvd (by omega) hd
    rcases Nat.lt_or_ge d n with hlt | hge
    · have h1 : d ∣ fact (n - 1) := dvd_fact (n - 1) d h0 (by omega)
      have h2 : n ∣ fact (n - 1) + 1 := by
        unfold MEq at h
        have : ((fact (n - 1) : Int) - -1) = ((fact (n - 1) + 1 : Nat) : Int) := by push_cast; omega
        rw [this] at h
        exact Int.natCast_dvd_natCast.mp h
      have h3 := (Nat.dvd_add_right h1).mp (Nat.dvd_trans hd h2)
      left; exact Nat.dvd_one.mp h3
    · right; omega
  · intro hp; exact wilson hp

example : fact 6 = 720 ∧ 720 % 7 = 6 := by decide
example : carries 3 5 4 0 = 2 ∧ layer 9 4 = 126 ∧ val 3 126 = 2 := by decide +kernel

end CircleOp.Cube
