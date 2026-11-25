/-
  Three classical proofs on the mult-add cube
  ==========================================

  Part of Fused Mult-Add Algebra. Companion to HopfInvariantOne.lean and
  the topology page. Pictures: ../fig_cube_proofs.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`, only Lean's standard
  axioms. Run it with `lean CubeProofs.lean`. Nothing is taken as given:
  everything is built here from the cube's layer counts.

  `layer n k` is the number of corners of the n-cube with k switches on.
  It is defined by Pascal's rule, which is one mult-add step: multiplying
  by (1 + x) adds an axis, and every old corner appears twice, once with
  the new switch off and once with it on.

  Proved:
  * binomial:     (1 + x)ⁿ = Σₖ layer n k · xᵏ. The cube's corner pieces
                  add up to the box (the geometric underpinning).
  * prime_dvd_middle: for a prime p, p divides every middle layer of the
                  p-cube.
  * aks_middle:   for n ≥ 2, n divides every middle layer of the n-cube
                  EXACTLY when n is prime. Equivalently (1 + x)ⁿ ≡ 1 + xⁿ
                  mod n, the identity behind the AKS primality test.
  * lucas_step, lucas: C(n,k) mod p is the product of the C's of the
                  base-p digits (Lucas 1878). Split the axes into blocks
                  of p; mod p, only corners whose blocks are all-off or
                  all-on count.
  * frobenius:    (1 + x)ᵖ ≡ 1 + xᵖ mod p: the mult-add p-th power is the
                  plain p-th power mod p.
  * fermat:       aᵖ ≡ a mod p for every natural number a (Fermat 1640).

  All three theorems are classical. What is new here is only the
  formal check in core Lean and the reading on the cube.

  LIBRARY NOTE: the cube layers (layer, absorb) now come from FusedAlgebra.Cube.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube

namespace MultAdd.Cube

/-! ## Primes -/

/-- A prime: at least 2, and its only divisors are 1 and itself. -/
def IsPrime (p : Nat) : Prop := 2 ≤ p ∧ ∀ d, d ∣ p → d = 1 ∨ d = p

theorem coprime_of_lt {p k : Nat} (hp : IsPrime p) (hk0 : 0 < k) (hkp : k < p) : Nat.Coprime p k := by
  rcases hp.2 _ (Nat.gcd_dvd_left p k) with h | h
  · exact h
  · have : Nat.gcd p k ≤ k := Nat.le_of_dvd hk0 (Nat.gcd_dvd_right p k)
    omega

/-- **The middle of the prime cube.** For a prime p, every middle layer
    of the p-cube has a multiple of p corners. -/
theorem prime_dvd_middle {p k : Nat} (hp : IsPrime p) (hk0 : 0 < k) (hkp : k < p) :
    p ∣ layer p k := by
  obtain ⟨q, rfl⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.1; omega⟩
  obtain ⟨j, rfl⟩ : ∃ j, k = j + 1 := ⟨k - 1, by omega⟩
  have h := absorb q j
  have hd : q + 1 ∣ (j + 1) * layer (q + 1) (j + 1) := ⟨layer q j, by rw [h]⟩
  exact (coprime_of_lt hp hk0 hkp).dvd_of_dvd_mul_left hd

/-! ## Congruence mod p, on the integers -/

def MEq (p : Nat) (x y : Int) : Prop := (p : Int) ∣ x - y

theorem MEq.refl (p : Nat) (x : Int) : MEq p x x := ⟨0, by simp⟩
theorem MEq.symm {p x y} (h : MEq p x y) : MEq p y x := by
  obtain ⟨c, hc⟩ := h; exact ⟨-c, by grind⟩
theorem MEq.trans {p x y z} (h1 : MEq p x y) (h2 : MEq p y z) : MEq p x z := by
  obtain ⟨c, hc⟩ := h1; obtain ⟨d, hd⟩ := h2; exact ⟨c + d, by grind⟩
theorem MEq.add {p x y u v} (h1 : MEq p x y) (h2 : MEq p u v) : MEq p (x + u) (y + v) := by
  obtain ⟨c, hc⟩ := h1; obtain ⟨d, hd⟩ := h2; exact ⟨c + d, by grind⟩
theorem MEq.mul {p x y u v} (h1 : MEq p x y) (h2 : MEq p u v) : MEq p (x * u) (y * v) := by
  obtain ⟨c, hc⟩ := h1; obtain ⟨d, hd⟩ := h2
  exact ⟨c * u + y * d, by
    have : x = y + p * c := by omega
    have : u = v + p * d := by omega
    subst_vars; grind⟩
theorem MEq.zero_of_dvd {p : Nat} {x : Int} (h : (p : Int) ∣ x) : MEq p x 0 := by
  obtain ⟨c, hc⟩ := h; exact ⟨c, by simp [hc]⟩

/-- Reading n in base p: n = p·q + r with r < p pins down q and r. -/
theorem divmod_of_eq {p q r : Nat} (hr : r < p) : (p * q + r) / p = q ∧ (p * q + r) % p = r := by
  have hp : 0 < p := by omega
  constructor
  · rw [Nat.add_comm, Nat.add_mul_div_left r q hp, Nat.div_eq_of_lt hr]; simp
  · rw [Nat.add_comm, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hr]

/-! ## Lucas's theorem -/

/-- The block count: split the n axes into blocks of p plus a remainder,
    and count corners whose blocks are each all-off or all-on. -/
def blockCount (p n k : Nat) : Nat := layer (n / p) (k / p) * layer (n % p) (k % p)

theorem blockCount_eq {p a b c d : Nat} (hb : b < p) (hd : d < p) :
    blockCount p (p * a + b) (p * c + d) = layer a c * layer b d := by
  unfold blockCount
  rw [(divmod_of_eq (q := a) hb).1, (divmod_of_eq (q := a) hb).2,
      (divmod_of_eq (q := c) hd).1, (divmod_of_eq (q := c) hd).2]

theorem layer_zero_left (k : Nat) : layer 0 k = if k = 0 then 1 else 0 := by
  cases k <;> rfl

/-- **Lucas, one step.** C(n,k) ≡ C(n div p, k div p) · C(n mod p, k mod p). -/
theorem lucas_step {p : Nat} (hp : IsPrime p) :
    ∀ n k, MEq p (layer n k) (blockCount p n k) := by
  have hp0 : 0 < p := by have := hp.1; omega
  intro n
  induction n with
  | zero =>
    intro k
    cases k with
    | zero =>
      unfold blockCount
      rw [Nat.zero_div, Nat.zero_mod]; exact MEq.refl _ _
    | succ k =>
      unfold blockCount
      rw [Nat.zero_div, Nat.zero_mod]
      have hk := Nat.div_add_mod (k + 1) p
      by_cases h : (k + 1) % p = 0
      · have hq : (k + 1) / p ≠ 0 := by
          intro hq; rw [hq, h] at hk; simp at hk
        simp only [layer_zero_left ((k + 1) / p), hq, ↓reduceIte]
        simp [layer]; exact MEq.refl _ _
      · simp only [layer_zero_left ((k + 1) % p), h, ↓reduceIte]
        simp [layer]; exact MEq.refl _ _
  | succ n ih =>
    intro k
    obtain ⟨a, b, hb, hn⟩ : ∃ a b, b < p ∧ n = p * a + b :=
      ⟨n / p, n % p, Nat.mod_lt n hp0, (Nat.div_add_mod n p).symm⟩
    cases k with
    | zero =>
      unfold blockCount
      rw [Nat.zero_div, Nat.zero_mod, layer_zero_right, layer_zero_right, layer_zero_right]
      exact MEq.refl _ _
    | succ k =>
      obtain ⟨c, d, hd, hk⟩ : ∃ c d, d < p ∧ k = p * c + d :=
        ⟨k / p, k % p, Nat.mod_lt k hp0, (Nat.div_add_mod k p).symm⟩
      have IH := MEq.add (ih k) (ih (k + 1))
      rw [layer_succ_succ]
      push_cast
      refine MEq.trans IH ?_
      -- now compare block counts
      rw [hn, hk]
      by_cases hb1 : b + 1 < p <;> by_cases hd1 : d + 1 < p
      · -- no carry anywhere
        rw [show p * a + b + 1 = p * a + (b + 1) by omega, show p * c + d + 1 = p * c + (d + 1) by omega,
          blockCount_eq hb hd, blockCount_eq hb hd1, blockCount_eq hb1 hd1, layer_succ_succ]
        push_cast; rw [Int.mul_add]; exact MEq.refl _ _
      · -- k carries, n doesn't
        have hdp : d + 1 = p := by omega
        rw [show p * a + b + 1 = p * a + (b + 1) by omega, show p * c + d + 1 = p * (c + 1) + 0 by grind,
          blockCount_eq hb hd, blockCount_eq hb hp0, blockCount_eq hb1 hp0,
          layer_eq_zero b d (by omega), layer_zero_right, layer_zero_right]
        simp; exact MEq.refl _ _
      · -- n carries, k doesn't: the prime middle appears
        have hbp : b + 1 = p := by omega
        rw [show p * a + b + 1 = p * (a + 1) + 0 by grind, show p * c + d + 1 = p * c + (d + 1) by omega,
          blockCount_eq hb hd, blockCount_eq hb hd1, blockCount_eq hp0 hd1,
          layer_zero_left]
        simp only [show d + 1 ≠ 0 by omega, ↓reduceIte]
        have hmid : p ∣ layer b d + layer b (d + 1) := by
          rw [← layer_succ_succ, hbp]; exact prime_dvd_middle hp (by omega) hd1
        obtain ⟨t, ht⟩ := hmid
        refine ⟨layer a c * t, ?_⟩
        have : ((layer b d : Nat) : Int) + (layer b (d + 1) : Nat) = (p : Int) * t := by
          rw [← Int.natCast_add, ht]; simp
        simp only [Int.natCast_mul, Nat.mul_zero, Int.natCast_zero, Int.sub_zero]
        rw [← Int.mul_add, this]; grind
      · -- both carry
        have hbp : b + 1 = p := by omega
        have hdp : d + 1 = p := by omega
        rw [show p * a + b + 1 = p * (a + 1) + 0 by grind, show p * c + d + 1 = p * (c + 1) + 0 by grind,
          blockCount_eq hb hd, blockCount_eq hb hp0, blockCount_eq hp0 hp0,
          show d = b by omega, layer_self, layer_zero_right, layer_zero_right, layer_succ_succ]
        simp; exact MEq.refl _ _

/-- Digit by digit in base p: the product of the layers of each digit pair. -/
def digitProd (p : Nat) : Nat → Nat → Nat → Nat
  | 0, _, _ => 1
  | f + 1, n, k => layer (n % p) (k % p) * digitProd p f (n / p) (k / p)

/-- **Lucas's theorem.** For n, k with at most f base-p digits,
    C(n,k) ≡ product over digits of C(nᵢ, kᵢ)  (mod p). -/
theorem lucas {p : Nat} (hp : IsPrime p) :
    ∀ f n k, n < p ^ f → k < p ^ f → MEq p (layer n k) (digitProd p f n k) := by
  intro f
  induction f with
  | zero =>
    intro n k hn hk
    have h1 : n = 0 := by simp at hn; omega
    have h2 : k = 0 := by simp at hk; omega
    subst h1; subst h2; exact MEq.refl _ _
  | succ f ih =>
    intro n k hn hk
    rw [Nat.pow_succ, Nat.mul_comm] at hn hk
    have h2 := ih (n / p) (k / p) (Nat.div_lt_of_lt_mul hn) (Nat.div_lt_of_lt_mul hk)
    refine MEq.trans (lucas_step hp n k) ?_
    unfold blockCount
    simp only [digitProd, Int.natCast_mul]
    rw [Int.mul_comm (layer (n / p) (k / p) : Int)]
    exact MEq.mul (MEq.refl _ _) h2

/-! ## Primes are exactly where the middle vanishes (the AKS identity) -/

theorem exists_prime_factor : ∀ n, 2 ≤ n → ∃ p, IsPrime p ∧ p ∣ n := by
  intro n
  induction n using Nat.strongRecOn with
  | _ n ih =>
    intro hn
    rcases Classical.em (IsPrime n) with hpn | hpn
    · exact ⟨n, hpn, Nat.dvd_refl n⟩
    · have hex : ∃ d, d ∣ n ∧ d ≠ 1 ∧ d ≠ n := by
        apply Classical.byContradiction
        intro h
        apply hpn
        refine ⟨hn, fun d hd => ?_⟩
        apply Classical.byContradiction
        intro h'
        exact h ⟨d, hd, by omega, by omega⟩
      obtain ⟨d, hd, hd1, hdn⟩ := hex
      have hd0 : d ≠ 0 := by
        intro h0; subst h0; have := Nat.eq_zero_of_zero_dvd hd; omega
      have hdlt : d < n := by have := Nat.le_of_dvd (by omega) hd; omega
      obtain ⟨p, hp, hpd⟩ := ih d hdlt (by omega)
      exact ⟨p, hp, Nat.dvd_trans hpd hd⟩

/-- **The AKS identity on the cube.** For n ≥ 2: every middle layer of the
    n-cube has a multiple of n corners exactly when n is prime.
    (Equivalently (1 + x)ⁿ ≡ 1 + xⁿ mod n: the mult-add n-th power equals
    the plain n-th power, which is the identity behind the AKS test.) -/
theorem aks_middle (n : Nat) (hn : 2 ≤ n) :
    (∀ k, 0 < k → k < n → n ∣ layer n k) ↔ IsPrime n := by
  constructor
  · intro hall
    apply Classical.byContradiction
    intro hnp
    obtain ⟨p, hp, hpn⟩ := exists_prime_factor n hn
    have hpne : p ≠ n := fun h => hnp (h ▸ hp)
    have hp2 := hp.1
    obtain ⟨m, hm⟩ := hpn
    obtain ⟨m', rfl⟩ : ∃ m', m = m' + 1 := ⟨m - 1, by
      cases m with
      | zero => simp at hm; omega
      | succ m => omega⟩
    have hm' : n = p * m' + p := by rw [hm, Nat.mul_succ]
    have hm1 : 1 ≤ m' := by
      cases m' with
      | zero => simp at hm'; omega
      | succ _ => omega
    have hpm : p ≤ p * m' := Nat.le_mul_of_pos_right p hm1
    obtain ⟨t, ht⟩ := hall p (by omega) (by omega)
    obtain ⟨q, hq⟩ : ∃ q, n = q + 1 := ⟨n - 1, by omega⟩
    obtain ⟨j, hj⟩ : ∃ j, p = j + 1 := ⟨p - 1, by omega⟩
    have hab := absorb q j
    rw [← hq, ← hj, ht] at hab
    have hcancel : layer q j = p * t := by
      apply Nat.eq_of_mul_eq_mul_left (show 0 < n by omega)
      rw [← hab]; grind
    have hjp : j < p := by omega
    have hqform : q = p * m' + j := by omega
    have hl := lucas_step hp q j
    have hb : blockCount p q j = 1 := by
      unfold blockCount
      rw [hqform, (divmod_of_eq (q := m') hjp).1, (divmod_of_eq (q := m') hjp).2,
        Nat.div_eq_of_lt hjp, Nat.mod_eq_of_lt hjp, layer_zero_right, layer_self]
    rw [hb, hcancel] at hl
    obtain ⟨c, hc⟩ := hl
    have h1 : (p : Int) ∣ 1 := ⟨t - c, by push_cast at hc; grind⟩
    have := Int.le_of_dvd (by decide) h1
    omega
  · intro hp k hk0 hkn
    exact prime_dvd_middle hp hk0 hkn

/-! ## The binomial theorem as the cube's corner count -/

def bsum (f : Nat → Int) : Nat → Int
  | 0 => 0
  | m + 1 => bsum f m + f m

theorem bsum_add (f g : Nat → Int) : ∀ m, bsum (fun k => f k + g k) m = bsum f m + bsum g m
  | 0 => by simp [bsum]
  | m + 1 => by simp only [bsum]; rw [bsum_add f g m]; grind

theorem bsum_mul (c : Int) (f : Nat → Int) : ∀ m, bsum (fun k => c * f k) m = c * bsum f m
  | 0 => by simp [bsum]
  | m + 1 => by simp only [bsum]; rw [bsum_mul c f m]; grind

theorem bsum_shift (f : Nat → Int) : ∀ m, bsum f (m + 1) = f 0 + bsum (fun k => f (k + 1)) m
  | 0 => by simp [bsum]
  | m + 1 => by
    rw [bsum, bsum_shift f m]; simp only [bsum]; grind

theorem bsum_congr (f g : Nat → Int) : ∀ m, (∀ k, k < m → f k = g k) → bsum f m = bsum g m
  | 0, _ => rfl
  | m + 1, h => by
    simp only [bsum]; rw [bsum_congr f g m (fun k hk => h k (by omega)), h m (by omega)]

theorem bsum_dvd (d : Int) (f : Nat → Int) : ∀ m, (∀ k, k < m → d ∣ f k) → d ∣ bsum f m
  | 0, _ => ⟨0, by simp [bsum]⟩
  | m + 1, h => by
    simp only [bsum]
    exact Int.dvd_add (bsum_dvd d f m (fun k hk => h k (by omega))) (h m (by omega))

/-- The cube sum: every corner of the n-cube, weighted x per switch that is on. -/
def cubeSum (n : Nat) (x : Int) : Int := bsum (fun k => (layer n k : Int) * x ^ k) (n + 1)

/-- **The cube's corner count.** (1 + x)ⁿ = Σₖ C(n,k) xᵏ: a box of side
    1 + x cut at 1 on every axis falls into the n-cube's corner pieces. -/
theorem binomial : ∀ n x, (1 + x) ^ n = cubeSum n x := by
  intro n x
  induction n with
  | zero => simp [cubeSum, bsum, layer]
  | succ n ih =>
    have hT : bsum (fun k => (layer n (k + 1) : Int) * x ^ (k + 1)) (n + 1) = cubeSum n x - 1 := by
      have hc : cubeSum n x = 1 + bsum (fun k => (layer n (k + 1) : Int) * x ^ (k + 1)) n := by
        unfold cubeSum; rw [bsum_shift]; simp [layer_zero_right]
      rw [hc, bsum, layer_eq_zero n (n + 1) (by omega)]; simp; omega
    have step : cubeSum (n + 1) x = x * cubeSum n x + cubeSum n x := by
      unfold cubeSum
      rw [bsum_shift]
      have e : (fun k => ((layer (n + 1) (k + 1) : Nat) : Int) * x ^ (k + 1)) =
          (fun k => x * ((layer n k : Int) * x ^ k) + (layer n (k + 1) : Int) * x ^ (k + 1)) := by
        funext k; rw [layer_succ_succ, Int.natCast_add, Int.pow_succ]; grind
      rw [e, bsum_add, bsum_mul]
      have hT' := hT; unfold cubeSum at hT'
      rw [hT']; simp [layer_zero_right]; grind
    rw [Int.pow_succ, ih, step]; grind

/-! ## Fermat's little theorem -/

/-- **Mult-add powers are plain powers mod p.** (1 + x)ᵖ ≡ 1 + xᵖ: in the
    p-cube only the two end corners survive mod p. -/
theorem frobenius {p : Nat} (hp : IsPrime p) (x : Int) : MEq p ((1 + x) ^ p) (1 + x ^ p) := by
  obtain ⟨q, hq⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.1; omega⟩
  rw [binomial]
  unfold cubeSum
  rw [hq, bsum, bsum_shift]
  have hmid : (p : Int) ∣ bsum (fun k => ((layer (q + 1) (k + 1) : Nat) : Int) * x ^ (k + 1)) q := by
    apply bsum_dvd
    intro k hk
    have := prime_dvd_middle hp (k := k + 1) (by omega) (by omega)
    rw [hq] at this
    obtain ⟨c, hc⟩ := this
    exact ⟨c * x ^ (k + 1), by rw [hc, hq]; push_cast; grind⟩
  obtain ⟨c, hc⟩ := hmid
  refine ⟨c, ?_⟩
  rw [layer_self, layer_zero_right, ← hq]
  simp only [Int.natCast_one, Int.pow_zero, Int.one_mul, Int.mul_one]
  rw [hq] at hc ⊢
  push_cast at hc ⊢
  omega

/-- **Fermat's little theorem.** aᵖ ≡ a (mod p) for every natural number a.
    Each step a → a + 1 adds one mult-add unit, and frobenius says the p-th
    power passes straight through it mod p. -/
theorem fermat {p : Nat} (hp : IsPrime p) : ∀ a : Nat, MEq p ((a : Int) ^ p) a := by
  have hp1 : 1 ≤ p := by have := hp.1; omega
  intro a
  induction a with
  | zero =>
    obtain ⟨q, hq⟩ : ∃ q, p = q + 1 := ⟨p - 1, by omega⟩
    rw [hq]; simp [Int.pow_succ]; exact MEq.refl _ _
  | succ a ih =>
    have h := frobenius hp (a : Int)
    have e : ((a + 1 : Nat) : Int) = 1 + (a : Int) := by push_cast; grind
    rw [e]
    exact MEq.trans h (MEq.add (MEq.refl _ 1) ih)

/-! ## Small checks -/

/-- 9 is not prime, and the 9-cube shows it: layer 3 has 84 corners, not a multiple of 9. -/
example : layer 9 3 = 84 ∧ ¬ (9 ∣ 84) := by decide
/-- Lucas in base 7: C(10,4) = 210, and 10 = 13₇, 4 = 04₇ gives C(1,0)·C(3,4) = 0 ≡ 210. -/
example : layer 10 4 = 210 ∧ digitProd 7 2 10 4 = 0 ∧ 210 % 7 = 0 := by decide

end MultAdd.Cube
