/-
  FreyPieces: the elementary pieces of the Wiles route to Fermat's Last Theorem
  ============================================================================

  Core Lean 4 (4.37.0-rc1), no mathlib, no `sorry`, no new axioms. Built on the shared
  library (FusedAlgebra.Base for mult-add, CubeProofs for `IsPrime`, Elliptic for the node).

  FLT n says: no x, y ≠ 0 in Nat with x^n + y^n = z^n.

  What is proved here (the deep steps of the route, modularity of semistable curves and
  Ribet's level lowering, are NOT proved here; these are the elementary pieces around them):

  1. Exponent reduction. FLT for d gives FLT for every multiple of d; every n ≥ 3 has 4 or an
     odd prime as a divisor; so FLT 4 plus FLT for odd primes gives all n ≥ 3, and FLT 4,
     FLT 3 plus FLT for primes p ≥ 5 (the Frey range) gives all n ≥ 3.
  2. Normalization. Euclid's lemma (Nat and Int); any solution shrinks to a coprime one;
     in a coprime solution exactly one of x, y, z is even; over Int, a coprime solution for
     odd p ≥ 5 can be re-signed and reordered so that a ≡ 3 mod 4, b even, 32 ∣ b^p.
  3. The Frey curve y² = x(x − A)(x + B) with A + B = C: its 2-torsion, its discriminant, the
     change of variables to the minimal model Y² + XY = X³ + a2·X² + a4·X, the invariants
     c4 = A² + AB + B² and 256·Δ = (ABC)², semistability (no prime divides both Δ and c4,
     including the node at 2), the reduction mod q | A to the node family, whose chord law
     is mult-add, and 256·Δ = ((abc)²)^p, so every odd prime's exact exponent in Δ is a
     multiple of p.
  4. Mult-add restatement: x^n + y^n = z^n ⇔ (y/x)^n = mpow (z/x − 1) n over Rat.
  5. Level 2: the counts μ = 3, ν2 = 1, ν3 = 0, ν∞ = 2 for Γ0(2), which put the genus of
     X0(2) at 0 through the (cited, not proved) genus formula.

  LIBRARY NOTE: copied from fermat/lean/FreyPieces.lean (Nov 25) as Topics/Fermat/FreyPieces.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Elliptic
import FusedAlgebra.Topics.Topology.CubeProofs

namespace FreyPieces
open MultAdd
open MultAdd.Cube (IsPrime exists_prime_factor MEq)

/-! ## 1. Exponent reduction -/

/-- FLT n: no x, y ≠ 0 and z in Nat with x^n + y^n = z^n. -/
def FLT (n : Nat) : Prop := ∀ x y z : Nat, x ≠ 0 → y ≠ 0 → x ^ n + y ^ n ≠ z ^ n

/-- 1a. FLT for d gives FLT for every multiple n of d (x^n = (x^(n/d))^d). -/
theorem flt_of_dvd {d n : Nat} (hdn : d ∣ n) (hd : FLT d) : FLT n := by
  obtain ⟨k, rfl⟩ := hdn
  intro x y z hx hy h
  apply hd (x ^ k) (y ^ k) (z ^ k) ((Nat.pos_iff_ne_zero.mp (Nat.pow_pos (Nat.pos_of_ne_zero hx)))) ((Nat.pos_iff_ne_zero.mp (Nat.pow_pos (Nat.pos_of_ne_zero hy))))
  rw [← Nat.pow_mul, ← Nat.pow_mul, ← Nat.pow_mul, Nat.mul_comm k d]
  exact h

/-- A prime other than 2 is odd. -/
theorem prime_ne_two_odd {p : Nat} (hp : IsPrime p) (h2 : p ≠ 2) : p % 2 = 1 := by
  rcases Nat.mod_two_eq_zero_or_one p with h | h
  · exfalso
    rcases hp.2 2 (Nat.dvd_of_mod_eq_zero h) with h' | h'
    · omega
    · exact h2 h'.symm
  · exact h

/-- 1b. Every n ≥ 3 is divisible by 4 or by an odd prime. -/
theorem four_or_odd_prime (n : Nat) (hn : 3 ≤ n) :
    4 ∣ n ∨ ∃ p, IsPrime p ∧ p ≠ 2 ∧ p ∣ n := by
  have oddcase : ∀ m, 2 ≤ m → m % 2 = 1 → ∃ p, IsPrime p ∧ p ≠ 2 ∧ p ∣ m := by
    intro m hm hodd
    obtain ⟨p, hp, hpm⟩ := exists_prime_factor m hm
    refine ⟨p, hp, ?_, hpm⟩
    rintro rfl
    have := Nat.mod_eq_zero_of_dvd hpm
    omega
  rcases Nat.mod_two_eq_zero_or_one n with h | h
  · obtain ⟨m, rfl⟩ := Nat.dvd_of_mod_eq_zero h
    rcases Nat.mod_two_eq_zero_or_one m with h' | h'
    · left
      obtain ⟨k, rfl⟩ := Nat.dvd_of_mod_eq_zero h'
      exact ⟨k, by omega⟩
    · right
      obtain ⟨p, hp, hp2, hpm⟩ := oddcase m (by omega) h'
      exact ⟨p, hp, hp2, Nat.dvd_trans hpm (Nat.dvd_mul_left m 2)⟩
  · right; exact oddcase n (by omega) h

/-- 1c. FLT 4 and FLT for every odd prime give FLT for every n ≥ 3. -/
theorem flt_reduction (h4 : FLT 4) (hp : ∀ p, IsPrime p → p ≠ 2 → FLT p) :
    ∀ n, 3 ≤ n → FLT n := by
  intro n hn
  rcases four_or_odd_prime n hn with h | ⟨p, hpp, hp2, hpn⟩
  · exact flt_of_dvd h h4
  · exact flt_of_dvd hpn (hp p hpp hp2)

/-- 1d. FLT 4, FLT 3 and FLT for every prime p ≥ 5 (the Frey range) give FLT for every n ≥ 3. -/
theorem flt_reduction_frey (h4 : FLT 4) (h3 : FLT 3) (hp : ∀ p, IsPrime p → 5 ≤ p → FLT p) :
    ∀ n, 3 ≤ n → FLT n := by
  apply flt_reduction h4
  intro p hpp hp2
  have h1 := hpp.1
  rcases (show p = 3 ∨ p = 4 ∨ 5 ≤ p by omega) with h | h | h
  · subst h; exact h3
  · subst h
    rcases hpp.2 2 ⟨2, rfl⟩ with h' | h' <;> omega
  · exact hp p hpp h

/-! ## 2. Normalization -/

/-! ### 2a. Euclid's lemma -/

/-- 2a. Euclid's lemma on Nat: a prime dividing a product divides a factor. -/
theorem euclid {p a b : Nat} (hp : IsPrime p) (h : p ∣ a * b) : p ∣ a ∨ p ∣ b := by
  rcases hp.2 _ (Nat.gcd_dvd_left p a) with h1 | h1
  · right; exact Nat.Coprime.dvd_of_dvd_mul_left h1 h
  · left; rw [← h1]; exact Nat.gcd_dvd_right p a

/-- 2a. Euclid's lemma on Int (through natAbs). -/
theorem euclid_int {p : Nat} {a b : Int} (hp : IsPrime p) (h : (p : Int) ∣ a * b) :
    (p : Int) ∣ a ∨ (p : Int) ∣ b := by
  rw [← Int.natAbs_dvd_natAbs, Int.natAbs_mul] at h
  simp only [Int.natAbs_natCast] at h
  rcases euclid hp h with h' | h'
  · left; rw [← Int.natAbs_dvd_natAbs]; simpa using h'
  · right; rw [← Int.natAbs_dvd_natAbs]; simpa using h'

/-- An odd prime is coprime to 2 (so to every power of 2). -/
theorem odd_prime_coprime_two {q : Nat} (hq : IsPrime q) (h2 : q ≠ 2) : Nat.Coprime q 2 := by
  have := prime_ne_two_odd hq h2
  rcases hq.2 _ (Nat.gcd_dvd_left q 2) with h | h
  · exact h
  · exfalso
    have hd : q ∣ 2 := h ▸ Nat.gcd_dvd_right q 2
    have h3 := Nat.le_of_dvd (by omega) hd
    have h4 := hq.1
    omega

/-! ### 2b. Coprime reduction -/

/-- 2b. A solution with x, y ≠ 0 shrinks (divide by gcd x y) to one with gcd x' y' = 1. -/
theorem coprime_reduction {n x y z : Nat} (hn : 1 ≤ n) (hx : x ≠ 0) (hy : y ≠ 0)
    (h : x ^ n + y ^ n = z ^ n) :
    ∃ x' y' z', x' ≠ 0 ∧ y' ≠ 0 ∧ Nat.gcd x' y' = 1 ∧ x' ^ n + y' ^ n = z' ^ n
      ∧ x' ∣ x ∧ y' ∣ y := by
  have hg : 0 < Nat.gcd x y := Nat.gcd_pos_of_pos_left y (Nat.pos_of_ne_zero hx)
  generalize hgd : Nat.gcd x y = g at hg
  obtain ⟨x', hx'⟩ : g ∣ x := hgd ▸ Nat.gcd_dvd_left x y
  obtain ⟨y', hy'⟩ : g ∣ y := hgd ▸ Nat.gcd_dvd_right x y
  have hcop : Nat.gcd x' y' = 1 := by
    have := Nat.gcd_mul_left g x' y'
    rw [← hx', ← hy', hgd] at this
    have h2 : g * 1 = g * Nat.gcd x' y' := by rw [Nat.mul_one]; exact this
    exact (Nat.eq_of_mul_eq_mul_left hg h2).symm
  have hsum : z ^ n = g ^ n * (x' ^ n + y' ^ n) := by
    rw [← h, hx', hy', Nat.mul_pow, Nat.mul_pow, Nat.mul_add]
  have hgz : g ∣ z := by
    rw [← Nat.pow_dvd_pow_iff (show n ≠ 0 by omega), hsum]
    exact Nat.dvd_mul_right _ _
  obtain ⟨z', hz'⟩ := hgz
  refine ⟨x', y', z', ?_, ?_, hcop, ?_, ⟨g, by rw [hx', Nat.mul_comm]⟩, ⟨g, by rw [hy', Nat.mul_comm]⟩⟩
  · rintro rfl; simp at hx'; exact hx hx'
  · rintro rfl; simp at hy'; exact hy hy'
  · rw [hz', Nat.mul_pow] at hsum
    exact (Nat.eq_of_mul_eq_mul_left ((Nat.pow_pos hg)) hsum).symm


/-- 2b. So FLT n holds once coprime solutions are ruled out. -/
theorem flt_of_coprime {n : Nat} (hn : 1 ≤ n)
    (h : ∀ x y z : Nat, x ≠ 0 → y ≠ 0 → Nat.gcd x y = 1 → x ^ n + y ^ n ≠ z ^ n) : FLT n := by
  intro x y z hx hy hxyz
  obtain ⟨x', y', z', hx', hy', hc, he, -, -⟩ := coprime_reduction hn hx hy hxyz
  exact h x' y' z' hx' hy' hc he

/-! ### 2c. Parity -/

/-- x^n has the parity of x for n ≥ 1. -/
theorem pow_mod_two {x n : Nat} (hn : 1 ≤ n) : x ^ n % 2 = x % 2 := by
  rw [Nat.pow_mod]
  rcases Nat.mod_two_eq_zero_or_one x with h | h <;> rw [h]
  · rw [Nat.zero_pow (by omega)]
  · rw [Nat.one_pow]

/-- 2c. In a coprime solution exactly one of x, y, z is even. -/
theorem coprime_parity {n x y z : Nat} (hn : 1 ≤ n) (hc : Nat.gcd x y = 1)
    (h : x ^ n + y ^ n = z ^ n) :
    (x % 2 = 0 ∧ y % 2 = 1 ∧ z % 2 = 1) ∨ (x % 2 = 1 ∧ y % 2 = 0 ∧ z % 2 = 1)
      ∨ (x % 2 = 1 ∧ y % 2 = 1 ∧ z % 2 = 0) := by
  have hx := pow_mod_two (x := x) hn
  have hy := pow_mod_two (x := y) hn
  have hz := pow_mod_two (x := z) hn
  have hm : (x ^ n + y ^ n) % 2 = z ^ n % 2 := by rw [h]
  rw [Nat.add_mod, hx, hy, hz] at hm
  have hnot : ¬ (x % 2 = 0 ∧ y % 2 = 0) := by
    rintro ⟨h1, h2⟩
    have : 2 ∣ Nat.gcd x y := Nat.dvd_gcd (Nat.dvd_of_mod_eq_zero h1) (Nat.dvd_of_mod_eq_zero h2)
    rw [hc] at this; omega
  omega

/-! ### 2d. Frey normalization over Int -/

/-- (−x)^(2k+1) = −x^(2k+1). -/
theorem neg_pow_odd (x : Int) (k : Nat) : (-x) ^ (2 * k + 1) = -(x ^ (2 * k + 1)) := by
  induction k with
  | zero => simp [Int.pow_succ]
  | succ k ih =>
    rw [show 2 * (k + 1) + 1 = (2 * k + 1) + 1 + 1 by omega, Int.pow_succ, Int.pow_succ, ih,
      Int.pow_succ x (2 * k + 1 + 1), Int.pow_succ x (2 * k + 1)]
    grind

/-- (−x)^p = −x^p for odd p. -/
theorem neg_pow_of_odd {p : Nat} (hp : p % 2 = 1) (x : Int) : (-x) ^ p = -(x ^ p) := by
  obtain ⟨k, rfl⟩ : ∃ k, p = 2 * k + 1 := ⟨p / 2, by omega⟩
  exact neg_pow_odd x k

/-- Over Int, x^n % 2 = x % 2 for n ≥ 1. -/
theorem int_pow_emod_two (x : Int) {n : Nat} (hn : 1 ≤ n) : x ^ n % 2 = x % 2 := by
  induction n with
  | zero => omega
  | succ n ih =>
    rcases Nat.eq_zero_or_pos n with h0 | h0
    · subst h0; simp [Int.pow_succ]
    · rw [Int.pow_succ, Int.mul_emod, ih h0]
      rcases (show x % 2 = 0 ∨ x % 2 = 1 by omega) with h | h <;> rw [h] <;> decide

/-- From o1^p + e^p + o2^p = 0 with o1 odd, e even: re-sign to a' ≡ 3 mod 4, b' even. -/
theorem normalize_core {p : Nat} (hp : p % 2 = 1) {o1 e o2 : Int}
    (h1 : o1 % 2 = 1) (he : e % 2 = 0) (h : o1 ^ p + e ^ p + o2 ^ p = 0) :
    ∃ a' b' c' : Int, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0
      ∧ (a' * b' * c').natAbs = (o1 * e * o2).natAbs := by
  rcases (show o1 % 4 = 3 ∨ o1 % 4 = 1 by omega) with h4 | h4
  · refine ⟨o1, e, -o2, ?_, h4, he, ?_⟩
    · rw [neg_pow_of_odd hp]; omega
    · simp [Int.natAbs_mul, Int.natAbs_neg]
  · refine ⟨-o1, -e, o2, ?_, by omega, by omega, ?_⟩
    · rw [neg_pow_of_odd hp, neg_pow_of_odd hp]; omega
    · simp [Int.natAbs_mul, Int.natAbs_neg]

/-- b even and p ≥ 5 give 32 ∣ b^p. -/
theorem two_pow_dvd {p : Nat} (h5 : 5 ≤ p) {b : Int} (hb : b % 2 = 0) : (32 : Int) ∣ b ^ p := by
  obtain ⟨k, rfl⟩ : ∃ k, b = 2 * k := ⟨b / 2, by omega⟩
  obtain ⟨m, rfl⟩ : ∃ m, p = 5 + m := ⟨p - 5, by omega⟩
  rw [mul_pow', pow_add']
  exact ⟨2 ^ m * k ^ (5 + m), by rw [show (2:Int) ^ 5 = 32 by decide]; grind⟩

/-- 2d. Frey normalization: a coprime Int solution for odd p ≥ 5 becomes one with a' ≡ 3 mod 4, b' even, 32 ∣ b'^p, same |abc|. -/
theorem frey_normalize {p : Nat} (hp : p % 2 = 1) (h5 : 5 ≤ p) {a b c : Int}
    (hcop : Int.gcd a b = 1) (habc : a * b * c ≠ 0) (h : a ^ p + b ^ p = c ^ p) :
    ∃ a' b' c' : Int, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0
      ∧ (a' * b' * c').natAbs = (a * b * c).natAbs ∧ a' * b' * c' ≠ 0 ∧ (32 : Int) ∣ b' ^ p := by
  have hp1 : 1 ≤ p := by omega
  have ha := int_pow_emod_two a hp1
  have hb := int_pow_emod_two b hp1
  have hc := int_pow_emod_two c hp1
  have hnot : ¬ (a % 2 = 0 ∧ b % 2 = 0) := by
    rintro ⟨h1, h2⟩
    have d1 : 2 ∣ a.natAbs := Int.natAbs_dvd_natAbs.mpr (Int.dvd_of_emod_eq_zero h1)
    have d2 : 2 ∣ b.natAbs := Int.natAbs_dvd_natAbs.mpr (Int.dvd_of_emod_eq_zero h2)
    have : 2 ∣ Int.gcd a b := Nat.dvd_gcd d1 d2
    rw [hcop] at this; omega
  have hm : (a ^ p + b ^ p) % 2 = c ^ p % 2 := by rw [h]
  rw [Int.add_emod, ha, hb, hc] at hm
  have hnc := neg_pow_of_odd hp c
  suffices H : ∃ a' b' c' : Int, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0
      ∧ (a' * b' * c').natAbs = (a * b * c).natAbs by
    obtain ⟨a', b', c', h1, h2, h3, h4⟩ := H
    refine ⟨a', b', c', h1, h2, h3, h4, ?_, two_pow_dvd h5 h3⟩
    intro h0
    rw [h0] at h4
    exact habc (Int.natAbs_eq_zero.mp h4.symm)
  rcases (show (a % 2 = 0 ∧ b % 2 = 1 ∧ c % 2 = 1) ∨ (a % 2 = 1 ∧ b % 2 = 0)
      ∨ (a % 2 = 1 ∧ b % 2 = 1 ∧ c % 2 = 0) by omega) with ⟨h1, h2, h3⟩ | ⟨h1, h2⟩ | ⟨h1, h2, h3⟩
  · obtain ⟨a', b', c', e1, e2, e3, e4⟩ := normalize_core hp h2 h1 (o2 := -c) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_⟩
    rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]; grind
  · obtain ⟨a', b', c', e1, e2, e3, e4⟩ := normalize_core hp h1 h2 (o2 := -c) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_⟩
    rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]
  · obtain ⟨a', b', c', e1, e2, e3, e4⟩ :=
      normalize_core hp h1 (show (-c) % 2 = 0 by omega) (o2 := b) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_⟩
    rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]; grind

/-! ## 3. The Frey curve y² = x(x − A)(x + B) -/

/-- The Frey cubic x(x − A)(x + B). -/
def freyRHS (A B x : Int) : Int := x * (x - A) * (x + B)

/-! ### 3a. Two-torsion -/

/-- 3a. The 2-torsion points (0,0), (A,0), (−B,0) lie on E. -/
theorem frey_two_torsion (A B : Int) :
    (0 : Int) * 0 = freyRHS A B 0 ∧ (0 : Int) * 0 = freyRHS A B A ∧ (0 : Int) * 0 = freyRHS A B (-B) := by
  unfold freyRHS; refine ⟨?_, ?_, ?_⟩ <;> grind

/-! ### 3b. Discriminant -/

/-- 3b. 16·(product of root differences)² = 16(ABC)² when A + B = C. -/
theorem frey_disc_roots (A B C : Int) (hC : A + B = C) :
    16 * ((0 - A) * (A - (-B)) * (0 - (-B))) ^ 2 = 16 * (A * B * C) ^ 2 := by
  subst hC; grind

/-! ### 3c. Minimal model -/

/-- 3c. The identity behind x = 4X, y = 8Y + 4X. -/
theorem frey_change_of_vars (A B X Y : Int) :
    (8 * Y + 4 * X) ^ 2 - 4 * X * (4 * X - A) * (4 * X + B)
      = 64 * (Y ^ 2 + X * Y - X ^ 3) - 16 * (B - A - 1) * X ^ 2 + 4 * (A * B) * X := by grind

/-- 3c. A ≡ 3 mod 4 and 16 ∣ B give 4 ∣ B − A − 1 and 16 ∣ AB. -/
theorem frey_minimal_divisibility {A B : Int} (hA : A % 4 = 3) (hB : 16 ∣ B) :
    4 ∣ B - A - 1 ∧ 16 ∣ A * B := by
  exact ⟨by omega, Int.dvd_trans hB (Int.dvd_mul_left A B)⟩

/-- a2 = (B − A − 1)/4 of the minimal model. -/
def a2of (A B : Int) : Int := (B - A - 1) / 4
/-- a4 = −AB/16 of the minimal model. -/
def a4of (A B : Int) : Int := -((A * B) / 16)

/-- 3c. With these integers, (4X, 8Y + 4X) is on E exactly when Y² + XY = X³ + a2X² + a4X. -/
theorem frey_minimal_model {A B : Int} (hA : A % 4 = 3) (hB : 16 ∣ B) :
    4 * a2of A B = B - A - 1 ∧ 16 * a4of A B = -(A * B) ∧
    ∀ X Y : Int, (8 * Y + 4 * X) ^ 2 = freyRHS A B (4 * X) ↔
      Y ^ 2 + X * Y = X ^ 3 + a2of A B * X ^ 2 + a4of A B * X := by
  have h1 : 4 * a2of A B = B - A - 1 := by unfold a2of; omega
  have h2 : 16 * a4of A B = -(A * B) := by
    unfold a4of
    obtain ⟨k, hk⟩ := hB
    have : A * B = 16 * (A * k) := by rw [hk]; grind
    rw [this]; generalize A * k = m; omega
  refine ⟨h1, h2, fun X Y => ?_⟩
  have key := frey_change_of_vars A B X Y
  unfold freyRHS
  generalize a2of A B = a2 at h1
  generalize a4of A B = a4 at h2
  constructor
  · intro h; grind
  · intro h; grind

/-! ### 3d. Invariants of the minimal model -/

/-- Weierstrass coefficients a1, a2, a3, a4, a6. -/
structure Weier where
  a1 : Int
  a2 : Int
  a3 : Int
  a4 : Int
  a6 : Int

namespace Weier
/-- b2 = a1² + 4a2. -/
def b2 (E : Weier) : Int := E.a1 ^ 2 + 4 * E.a2
/-- b4 = 2a4 + a1a3. -/
def b4 (E : Weier) : Int := 2 * E.a4 + E.a1 * E.a3
/-- b6 = a3² + 4a6. -/
def b6 (E : Weier) : Int := E.a3 ^ 2 + 4 * E.a6
/-- b8 = a1²a6 + 4a2a6 − a1a3a4 + a2a3² − a4². -/
def b8 (E : Weier) : Int :=
  E.a1 ^ 2 * E.a6 + 4 * E.a2 * E.a6 - E.a1 * E.a3 * E.a4 + E.a2 * E.a3 ^ 2 - E.a4 ^ 2
/-- c4 = b2² − 24b4. -/
def c4 (E : Weier) : Int := E.b2 ^ 2 - 24 * E.b4
/-- Δ = −b2²b8 − 8b4³ − 27b6² + 9b2b4b6. -/
def disc (E : Weier) : Int :=
  -(E.b2 ^ 2 * E.b8) - 8 * E.b4 ^ 3 - 27 * E.b6 ^ 2 + 9 * E.b2 * E.b4 * E.b6
end Weier

/-- The model Y² + XY = X³ + a2X² + a4X (a1 = 1, a3 = a6 = 0). -/
def freyModel (a2 a4 : Int) : Weier := ⟨1, a2, 0, a4, 0⟩

/-- 3d. c4 of the Frey model is A² + AB + B². -/
theorem frey_c4 {A B a2 a4 : Int} (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) :
    (freyModel a2 a4).c4 = A ^ 2 + A * B + B ^ 2 := by
  simp only [Weier.c4, Weier.b2, Weier.b4, freyModel]
  grind

/-- 3d. 256·Δ of the Frey model is (ABC)². -/
theorem frey_disc {A B C a2 a4 : Int} (hC : A + B = C) (h2 : 4 * a2 = B - A - 1)
    (h4 : 16 * a4 = -(A * B)) :
    256 * (freyModel a2 a4).disc = (A * B * C) ^ 2 := by
  subst hC
  simp only [Weier.disc, Weier.b2, Weier.b4, Weier.b6, Weier.b8, freyModel]
  have hb2 : 1 ^ 2 + 4 * a2 = B - A := by grind
  rw [hb2]
  grind

/-! ### 3e. Semistability: only nodes -/

/-- A prime cannot divide two coprime integers. -/
theorem prime_dvd_gcd_false {q : Nat} (hq : IsPrime q) {A B : Int} (hcop : Int.gcd A B = 1)
    (hA : (q : Int) ∣ A) (hB : (q : Int) ∣ B) : False := by
  have d1 : q ∣ A.natAbs := by simpa using Int.natAbs_dvd_natAbs.mpr hA
  have d2 : q ∣ B.natAbs := by simpa using Int.natAbs_dvd_natAbs.mpr hB
  have : q ∣ Int.gcd A B := Nat.dvd_gcd d1 d2
  rw [hcop] at this
  have := Nat.le_of_dvd (by omega) this
  have := hq.1
  omega

/-- 3e. If gcd(A, B) = 1, a prime dividing ABC does not divide A² + AB + B². -/
theorem frey_c4_not_dvd {A B C : Int} (hC : A + B = C) (hcop : Int.gcd A B = 1) {q : Nat}
    (hq : IsPrime q) (hd : (q : Int) ∣ A * B * C) : ¬ (q : Int) ∣ A ^ 2 + A * B + B ^ 2 := by
  intro hc
  rcases euclid_int hq hd with hab | hC'
  · rcases euclid_int hq hab with hA | hB
    · have : (q : Int) ∣ B * B := by
        have e : B * B = (A ^ 2 + A * B + B ^ 2) - A * (A + B) := by grind
        rw [e]; exact Int.dvd_sub hc (Int.dvd_trans hA (Int.dvd_mul_right A _))
      rcases euclid_int hq this with h | h <;> exact prime_dvd_gcd_false hq hcop hA h
    · have : (q : Int) ∣ A * A := by
        have e : A * A = (A ^ 2 + A * B + B ^ 2) - B * (A + B) := by grind
        rw [e]; exact Int.dvd_sub hc (Int.dvd_trans hB (Int.dvd_mul_right B _))
      rcases euclid_int hq this with h | h <;> exact prime_dvd_gcd_false hq hcop h hB
  · have : (q : Int) ∣ A * A := by
      have e : A * A = (A ^ 2 + A * B + B ^ 2) - B * C := by subst hC; grind
      rw [e]; exact Int.dvd_sub hc (Int.dvd_trans hC' (Int.dvd_mul_left B C))
    rcases euclid_int hq this with h | h <;>
    · have hB : (q : Int) ∣ B := by
        have e : B = C - A := by omega
        rw [e]; exact Int.dvd_sub hC' h
      exact prime_dvd_gcd_false hq hcop h hB

/-- 3e. Semistability: a prime dividing Δ never divides c4 (node, not cusp). -/
theorem frey_semistable {A B C a2 a4 : Int} (hC : A + B = C) (hcop : Int.gcd A B = 1)
    (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) {q : Nat} (hq : IsPrime q)
    (hd : (q : Int) ∣ (freyModel a2 a4).disc) : ¬ (q : Int) ∣ (freyModel a2 a4).c4 := by
  rw [frey_c4 h2 h4]
  apply frey_c4_not_dvd hC hcop hq
  have h256 : (q : Int) ∣ (A * B * C) * (A * B * C) := by
    rw [show (A * B * C) * (A * B * C) = (A * B * C) ^ 2 by grind, ← frey_disc hC h2 h4]
    exact Int.dvd_trans hd (Int.dvd_mul_left 256 _)
  rcases euclid_int hq h256 with h | h <;> exact h

/-- 3e. At 2: A odd and 32 ∣ B make c4 odd and 4 ∣ Δ, so 2 is a bad prime with a node. -/
theorem frey_at_two {A B C a2 a4 : Int} (hC : A + B = C) (hA : A % 2 = 1) (hB : 32 ∣ B)
    (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) :
    (freyModel a2 a4).c4 % 2 = 1 ∧ 4 ∣ (freyModel a2 a4).disc := by
  obtain ⟨m, rfl⟩ := hB
  constructor
  · rw [frey_c4 h2 h4]
    obtain ⟨k, rfl⟩ : ∃ k, A = 2 * k + 1 := ⟨A / 2, by omega⟩
    have : (2 * k + 1) ^ 2 + (2 * k + 1) * (32 * m) + (32 * m) ^ 2
        = 2 * (2 * k * k + 2 * k + 32 * k * m + 16 * m + 512 * m * m) + 1 := by grind
    rw [this]; omega
  · have hD := frey_disc hC h2 h4
    have : (A * (32 * m) * C) ^ 2 = 1024 * (A * m * C) ^ 2 := by grind
    rw [this] at hD
    exact ⟨(A * m * C) ^ 2, by omega⟩

/-! ### 3f. Node shape and the mult-add chord law -/

/-- 3f. For an odd prime q dividing A, q divides a4. -/
theorem frey_a4_mod_q {A B a4 : Int} (h4 : 16 * a4 = -(A * B)) {q : Nat} (hq : IsPrime q)
    (h2 : q ≠ 2) (hA : (q : Int) ∣ A) : (q : Int) ∣ a4 := by
  have h16 : (q : Int) ∣ 16 * a4 := by
    rw [h4]; exact Int.dvd_neg.mpr (Int.dvd_trans hA (Int.dvd_mul_right A B))
  rw [← Int.natAbs_dvd_natAbs, Int.natAbs_mul] at h16
  simp only [Int.natAbs_natCast] at h16
  have hc : Nat.Coprime q 16 := Nat.Coprime.pow_right 4 (odd_prime_coprime_two hq h2)
  have := Nat.Coprime.dvd_of_dvd_mul_left hc h16
  rw [← Int.natAbs_dvd_natAbs]; simpa using this

/-- 3f. Then mod q the model is the node family Y² + 1·XY = X³ + a2X². -/
theorem frey_reduces_to_node {a2 a4 : Int} {q : Nat} (h : (q : Int) ∣ a4) (X Y : Int) :
    MEq q (Y * Y + X * Y - (X * X * X + a2 * X * X + a4 * X))
          (Y * Y + 1 * X * Y - (X * X * X + a2 * X * X)) := by
  unfold MEq
  rw [show Y * Y + X * Y - (X * X * X + a2 * X * X + a4 * X) - (Y * Y + 1 * X * Y - (X * X * X + a2 * X * X))
      = a4 * (-X) by grind]
  exact Int.dvd_trans h (Int.dvd_mul_right a4 _)

section Chord
variable {R : Type} [Lean.Grind.CommRing R]

/-- The x-coordinate of the second meeting point of the line y = t·x with the node of slopes m1, m2. -/
def nodeX (m1 m2 t : R) : R := (t - m1) * (t - m2)

/-- The collinearity determinant of three points. -/
def chordDet (x1 y1 x2 y2 x3 y3 : R) : R := x1 * (y2 - y3) - y1 * (x2 - x3) + (x2 * y3 - x3 * y2)

/-- 3f. (m2 − m1)·det = Vandermonde · (∏(tᵢ − m1) − ∏(tᵢ − m2)). -/
theorem node_chord_identity (m1 m2 t1 t2 t3 : R) :
    (m2 - m1) * chordDet (nodeX m1 m2 t1) (t1 * nodeX m1 m2 t1) (nodeX m1 m2 t2)
        (t2 * nodeX m1 m2 t2) (nodeX m1 m2 t3) (t3 * nodeX m1 m2 t3)
      = (t1 - t2) * (t2 - t3) * (t3 - t1) *
        ((t1 - m1) * (t2 - m1) * (t3 - m1) - (t1 - m2) * (t2 - m2) * (t3 - m2)) := by
  unfold chordDet nodeX; grind
end Chord

/-- 3f. The point (x, t·x), x = nodeX, lies on y² + a1xy = x³ + a2x² with a1 = −(m1+m2), a2 = −m1m2 (library node_param). -/
theorem node_point (m1 m2 t : Int) :
    let x := nodeX m1 m2 t
    (t * x) * (t * x) + (-(m1 + m2)) * x * (t * x) = x * x * x + (-(m1 * m2)) * x * x := by
  intro x
  have h := EllipticLead.node_param (-(m1 + m2)) (-(m1 * m2)) t
  have e : x = t * t + (-(m1 + m2)) * t - (-(m1 * m2)) := by
    simp only [x, nodeX]; rw [EllipticLead.slopes_factor]
  rw [e]; exact h

/-- Nonzero times nonzero is nonzero in Rat. -/
theorem rat_mul_ne_zero {a b : Rat} (ha : a ≠ 0) (hb : b ≠ 0) : a * b ≠ 0 := by
  intro h; rcases Rat.mul_eq_zero.mp h with h | h <;> contradiction

/-- (1 + d/b)·b = b + d for b ≠ 0. -/
theorem one_add_div_mul {b d : Rat} (h : b ≠ 0) : (1 + d / b) * b = b + d := by
  rw [Rat.add_mul, Rat.div_mul_cancel h, Rat.one_mul]

/-- 3f. Link: three points of the node are collinear exactly when s1 ⊞ s2 ⊞ s3 = 0, sᵢ = (m2 − m1)/(tᵢ − m2). -/
theorem node_chord_multAdd (m1 m2 t1 t2 t3 : Rat) (hm : m1 ≠ m2)
    (h12 : t1 ≠ t2) (h23 : t2 ≠ t3) (h31 : t3 ≠ t1)
    (h1 : t1 ≠ m2) (h2 : t2 ≠ m2) (h3 : t3 ≠ m2) :
    chordDet (nodeX m1 m2 t1) (t1 * nodeX m1 m2 t1) (nodeX m1 m2 t2)
        (t2 * nodeX m1 m2 t2) (nodeX m1 m2 t3) (t3 * nodeX m1 m2 t3) = 0 ↔
      multAdd (multAdd ((m2 - m1) / (t1 - m2)) ((m2 - m1) / (t2 - m2))) ((m2 - m1) / (t3 - m2)) = 0 := by
  have key := node_chord_identity m1 m2 t1 t2 t3
  have sub_ne : ∀ u v : Rat, u ≠ v → u - v ≠ 0 := fun u v h e => h (by grind)
  have e1 := sub_ne _ _ h1
  have e2 := sub_ne _ _ h2
  have e3 := sub_ne _ _ h3
  have hlift : (1 + multAdd (multAdd ((m2 - m1) / (t1 - m2)) ((m2 - m1) / (t2 - m2))) ((m2 - m1) / (t3 - m2)))
      * ((t1 - m2) * (t2 - m2) * (t3 - m2)) = (t1 - m1) * (t2 - m1) * (t3 - m1) := by
    rw [multAdd_lift, multAdd_lift]
    have f1 := one_add_div_mul (d := m2 - m1) e1
    have f2 := one_add_div_mul (d := m2 - m1) e2
    have f3 := one_add_div_mul (d := m2 - m1) e3
    generalize 1 + (m2 - m1) / (t1 - m2) = u1 at f1 ⊢
    generalize 1 + (m2 - m1) / (t2 - m2) = u2 at f2 ⊢
    generalize 1 + (m2 - m1) / (t3 - m2) = u3 at f3 ⊢
    grind
  have hV : (t1 - t2) * (t2 - t3) * (t3 - t1) ≠ 0 :=
    rat_mul_ne_zero (rat_mul_ne_zero (sub_ne _ _ h12) (sub_ne _ _ h23)) (sub_ne _ _ h31)
  have hB : (t1 - m2) * (t2 - m2) * (t3 - m2) ≠ 0 := rat_mul_ne_zero (rat_mul_ne_zero e1 e2) e3
  have hd : m2 - m1 ≠ 0 := sub_ne _ _ (fun h => hm h.symm)
  generalize chordDet (nodeX m1 m2 t1) (t1 * nodeX m1 m2 t1) (nodeX m1 m2 t2)
        (t2 * nodeX m1 m2 t2) (nodeX m1 m2 t3) (t3 * nodeX m1 m2 t3) = D at key ⊢
  generalize multAdd (multAdd ((m2 - m1) / (t1 - m2)) ((m2 - m1) / (t2 - m2))) ((m2 - m1) / (t3 - m2)) = S at hlift ⊢
  generalize (t1 - t2) * (t2 - t3) * (t3 - t1) = V at key hV
  generalize (t1 - m2) * (t2 - m2) * (t3 - m2) = Pb at key hlift hB
  generalize (t1 - m1) * (t2 - m1) * (t3 - m1) = Pa at key hlift
  constructor
  · intro hD
    rw [hD] at key
    have : Pa - Pb = 0 := by
      rcases Rat.mul_eq_zero.mp (show V * (Pa - Pb) = 0 by grind) with h | h
      · exact absurd h hV
      · exact h
    have : S * Pb = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h | h
    · exact h
    · exact absurd h hB
  · intro hS
    rw [hS] at hlift
    have : (m2 - m1) * D = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h | h
    · exact absurd h hd
    · exact h

/-! ### 3g. The Tate parameter is a p-th power -/

/-- 3g. With A = a^p, B = b^p, C = c^p: 256·Δ = ((abc)²)^p. -/
theorem frey_disc_pth_power {a b c Δ : Int} {p : Nat} (hΔ : 256 * Δ = (a ^ p * b ^ p * c ^ p) ^ 2) :
    256 * Δ = ((a * b * c) ^ 2) ^ p := by
  rw [hΔ, ← mul_pow', ← mul_pow', ← Int.pow_mul, ← Int.pow_mul, Nat.mul_comm]

/-- Pull every factor q out of a nonzero M: M = q^k·m with q ∤ m. -/
theorem factor_out {q : Nat} (hq : 2 ≤ q) : ∀ M, M ≠ 0 → ∃ k m, M = q ^ k * m ∧ ¬ q ∣ m := by
  intro M
  induction M using Nat.strongRecOn with
  | _ M ih =>
    intro hM
    by_cases hd : q ∣ M
    · obtain ⟨M', rfl⟩ := hd
      have hM' : M' ≠ 0 := by rintro rfl; simp at hM
      have hlt : M' < q * M' := by
        have : 1 * M' < q * M' := Nat.mul_lt_mul_of_pos_right (by omega) (by omega)
        omega
      obtain ⟨k, m, hk, hm⟩ := ih M' hlt hM'
      exact ⟨k + 1, m, by rw [hk, Nat.pow_succ]; grind, hm⟩
    · exact ⟨0, M, by simp, hd⟩

/-- A prime not dividing m does not divide m^p. -/
theorem not_dvd_pow {q m : Nat} (hq : IsPrime q) (hm : ¬ q ∣ m) : ∀ p, ¬ q ∣ m ^ p := by
  intro p
  induction p with
  | zero =>
    intro h
    rw [Nat.pow_zero] at h
    have h3 := Nat.le_of_dvd (by decide) h
    have h4 := hq.1
    omega
  | succ p ih =>
    intro h
    rw [Nat.pow_succ] at h
    rcases euclid hq h with h | h
    · exact ih h
    · exact hm h

/-- q^a ∣ q^b·m with q ∤ m forces a ≤ b. -/
theorem pow_dvd_bound {q a b m : Nat} (hq : 2 ≤ q) (hm : ¬ q ∣ m) (h : q ^ a ∣ q ^ b * m) : a ≤ b := by
  apply Classical.byContradiction
  intro hab
  have h1 : q ^ (b + 1) ∣ q ^ b * m := Nat.dvd_trans (Nat.pow_dvd_pow q (by omega)) h
  rw [Nat.pow_succ] at h1
  exact hm (Nat.dvd_of_mul_dvd_mul_left (Nat.pow_pos (by omega)) h1)

/-- The exact exponent of a prime in a p-th power is a multiple of p. -/
theorem exact_exponent_pth_power {q M e : Nat} (p : Nat) (hq : IsPrime q) (hM : M ≠ 0)
    (he : q ^ e ∣ M ^ p) (he1 : ¬ q ^ (e + 1) ∣ M ^ p) : p ∣ e := by
  obtain ⟨k, m, rfl, hm⟩ := factor_out hq.1 M hM
  rw [Nat.mul_pow, ← Nat.pow_mul] at he he1
  have hmp := not_dvd_pow hq hm p
  have h1 := pow_dvd_bound hq.1 hmp he
  have h2 : ¬ e + 1 ≤ k * p := fun h =>
    he1 (Nat.dvd_trans (Nat.pow_dvd_pow q h) (Nat.dvd_mul_right _ _))
  exact ⟨k, by rw [Nat.mul_comm]; omega⟩

/-- 3g bonus: for every odd prime q, the exact exponent of q in Δ is a multiple of p. -/
theorem frey_disc_valuation {a b c Δ : Int} {p : Nat}
    (hΔ : 256 * Δ = (a ^ p * b ^ p * c ^ p) ^ 2) (habc : a * b * c ≠ 0)
    {q : Nat} (hq : IsPrime q) (h2 : q ≠ 2) {e : Nat}
    (he : (q : Int) ^ e ∣ Δ) (he1 : ¬ (q : Int) ^ (e + 1) ∣ Δ) : p ∣ e := by
  have hP := frey_disc_pth_power hΔ
  have hN : 256 * Δ.natAbs = ((a * b * c).natAbs ^ 2) ^ p := by
    have := congrArg Int.natAbs hP
    simpa [Int.natAbs_mul, Int.natAbs_pow] using this
  have hcop : ∀ f, Nat.Coprime (q ^ f) 256 := fun f =>
    Nat.Coprime.pow f 8 (odd_prime_coprime_two hq h2)
  have tr : ∀ f, (q : Int) ^ f ∣ Δ ↔ q ^ f ∣ ((a * b * c).natAbs ^ 2) ^ p := by
    intro f
    rw [← Int.natAbs_dvd_natAbs, Int.natAbs_pow, Int.natAbs_natCast, ← hN]
    constructor
    · intro h; exact Nat.dvd_trans h (Nat.dvd_mul_left _ _)
    · intro h; exact Nat.Coprime.dvd_of_dvd_mul_left (hcop f) h
  have hM : (a * b * c).natAbs ^ 2 ≠ 0 := by
    have : (a * b * c).natAbs ≠ 0 := by simpa using habc
    exact Nat.pos_iff_ne_zero.mp (Nat.pow_pos (Nat.pos_of_ne_zero this))
  exact exact_exponent_pth_power p hq hM ((tr e).mp he) (fun h => he1 ((tr (e + 1)).mpr h))

/-! ## 4. Mult-add restatement -/

/-- (y/x)^n · x^n = y^n for x ≠ 0. -/
theorem rat_div_pow {x : Rat} (y : Rat) (hx : x ≠ 0) (n : Nat) : (y / x) ^ n * x ^ n = y ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_succ]
    have h := Rat.div_mul_cancel (a := y) hx
    have e : ∀ W w X : Rat, W * w * (X * x) = (W * X) * (w * x) := by intros; grind
    rw [e, ih, h]

/-- 4. x^n + y^n = z^n exactly when (y/x)^n = mpow (z/x − 1) n: an n-th power equals a mult-add n-th power. -/
theorem flt_multAdd (x y z : Rat) (hx : x ≠ 0) (n : Nat) :
    x ^ n + y ^ n = z ^ n ↔ (y / x) ^ n = mpow (z / x - 1) n := by
  rw [mpow_eq]
  have hxn : x ^ n ≠ 0 := by
    induction n with
    | zero => simp
    | succ n ih => rw [Rat.pow_succ]; intro h; rcases Rat.mul_eq_zero.mp h with h | h <;> contradiction
  have e1 : (1 + (z / x - 1)) = z / x := by grind
  rw [e1]
  have hy := rat_div_pow y hx n
  have hz := rat_div_pow z hx n
  generalize (y / x) ^ n = u at hy ⊢
  generalize (z / x) ^ n = v at hz ⊢
  constructor
  · intro h
    have : (u - (v - 1)) * x ^ n = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h | h
    · grind
    · contradiction
  · intro h; subst h; grind

/-! ## 5. Level 2 target: X0(2) has genus 0 -/

/-- μ = [SL2(Z) : Γ0(2)] = #P¹(Z/2) = 3, counted as the nonzero vectors of (Z/2)² (units of Z/2 are {1}). -/
def mu2 : Nat := ((List.range 4).filter (fun v => v ≠ 0)).length
/-- ν2 = #{x mod 2 : x² + 1 ≡ 0}. -/
def nu2 : Nat := ((List.range 2).filter (fun x => (x * x + 1) % 2 = 0)).length
/-- ν3 = #{x mod 2 : x² + x + 1 ≡ 0}. -/
def nu3 : Nat := ((List.range 2).filter (fun x => (x * x + x + 1) % 2 = 0)).length
/-- ν∞ = Σ_{d ∣ 2} φ(gcd(d, 2/d)) = number of divisors of 2 (each gcd is 1). -/
def nuInf : Nat := ((List.range 3).filter (fun d => d ≠ 0 ∧ 2 % d = 0)).length

/-- 5. μ = 3, ν2 = 1, ν3 = 0, ν∞ = 2 for Γ0(2), by computation. -/
theorem level2_counts : mu2 = 3 ∧ nu2 = 1 ∧ nu3 = 0 ∧ nuInf = 2 := by decide

/-- 5. The genus formula 12g = 12 + μ − 3ν2 − 4ν3 − 6ν∞ (cited, not proved here) gives g(X0(2)) = 0. -/
theorem genus_X0_2 (g : Int)
    (hformula : 12 * g = 12 + (mu2 : Int) - 3 * nu2 - 4 * nu3 - 6 * nuInf) : g = 0 := by
  obtain ⟨h1, h2, h3, h4⟩ := level2_counts
  rw [h1, h2, h3, h4] at hformula
  omega


end FreyPieces

/-! ## Axiom check: only propext, Classical.choice, Quot.sound may appear. -/
