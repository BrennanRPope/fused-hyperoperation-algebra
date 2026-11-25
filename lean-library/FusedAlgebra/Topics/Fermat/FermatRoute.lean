/-
  FermatRoute: the whole of Fermat's Last Theorem, minus exactly one deep step.

  This file joins the three proof files:
    Fermat4.lean     exponent 4 (Fermat's descent)            fully proved
    Fermat3.lean     exponent 3 (Euler's descent in ℤ[ω])     fully proved
    FreyPieces.lean  reductions, normalization, the Frey curve fully proved

  What it proves:
    • `flt_of_prime_case`: if FLT holds for every prime p ≥ 5, it holds for every n ≥ 3.
    • `flt_of_no_frey_curve`: FLT for every n ≥ 3 follows from ONE statement, `NoFreyCurve`:
      no coprime solution a^p + b^p = c^p (p ≥ 5 prime, abc ≠ 0) exists in the
      normalized shape a ≡ 3 (mod 4), b even, 32 ∣ b^p. That is exactly the shape in which
      the Frey curve y² = x(x − a^p)(x + b^p) is semistable (FreyPieces 3e) with discriminant
      a p-th power up to 2^8 (3g).

  What is NOT proved here: `NoFreyCurve` itself. It is the content of Wiles and Taylor–Wiles
  (semistable curves are modular), Mazur (the mod-p picture is irreducible) and Ribet (level
  lowering to N = 2), finished by the genus-0 count of X0(2) (FreyPieces 5). It stays an
  explicit hypothesis, never an axiom and never `sorry`.

  LIBRARY NOTE: copied from fermat/lean/FermatRoute.lean (Nov 25) as Topics/Fermat/FermatRoute.
-/
import FusedAlgebra.Topics.Fermat.Fermat3
import FusedAlgebra.Topics.Fermat.Fermat4
import FusedAlgebra.Topics.Fermat.FreyPieces

open CircleOp
open CircleOp.Cube (IsPrime)
open FreyPieces

namespace FermatRoute

/-- FLT for every n ≥ 3 follows from the prime case p ≥ 5 alone (n = 3 and n = 4 are proved). -/
theorem flt_of_prime_case (hp : ∀ p, IsPrime p → 5 ≤ p → FLT p) : ∀ n, 3 ≤ n → FLT n :=
  flt_reduction_frey fermat4 Fermat3.flt3_nat hp

/-- A prime dividing x^n divides x (Int). -/
theorem prime_dvd_of_dvd_pow {q : Nat} (hq : IsPrime q) {x : Int} :
    ∀ {n : Nat}, (q : Int) ∣ x ^ n → (q : Int) ∣ x
  | 0, h => by
      exfalso
      rw [Int.pow_zero, Int.ofNat_dvd_left] at h
      have := Nat.le_of_dvd (by decide) h
      have := hq.1; simp at *; omega
  | n + 1, h => by
      rw [Int.pow_succ] at h
      rcases euclid_int hq h with h1 | h1
      · exact prime_dvd_of_dvd_pow hq h1
      · exact h1

/-- Two nonzero Ints with no common prime factor are coprime. -/
theorem gcd_one_of_no_prime {a c : Int} (ha : a ≠ 0)
    (H : ∀ q, IsPrime q → (q : Int) ∣ a → (q : Int) ∣ c → False) : Int.gcd a c = 1 := by
  have hpos : 0 < Int.gcd a c := Int.gcd_pos_of_ne_zero_left c ha
  by_cases h1 : Int.gcd a c = 1
  · exact h1
  · obtain ⟨q, hq, hqd⟩ := CircleOp.Cube.exists_prime_factor (Int.gcd a c) (by omega)
    exfalso
    apply H q hq
    · exact Int.ofNat_dvd_left.mpr (Nat.dvd_trans hqd (Nat.gcd_dvd_left _ _))
    · exact Int.ofNat_dvd_left.mpr (Nat.dvd_trans hqd (Nat.gcd_dvd_right _ _))

/-- x divides x^p for p ≥ 1. -/
theorem dvd_self_pow (x : Int) {p : Nat} (hp : 1 ≤ p) : x ∣ x ^ p := by
  obtain ⟨m, rfl⟩ : ∃ m, p = m + 1 := ⟨p - 1, by omega⟩
  rw [Int.pow_succ]; exact Int.dvd_mul_left _ _

/-- In a solution a^p + b^p = c^p with gcd(a, b) = 1, also gcd(a, c) = 1. -/
theorem cop_left {p : Nat} {a b c : Int} (ha : a ≠ 0) (hcop : Int.gcd a b = 1)
    (h : a ^ p + b ^ p = c ^ p) (hp : 1 ≤ p) : Int.gcd a c = 1 := by
  apply gcd_one_of_no_prime ha
  intro q hq hqa hqc
  have hb : (q : Int) ∣ b ^ p := by
    have : b ^ p = c ^ p - a ^ p := by omega
    rw [this]
    exact Int.dvd_sub (Int.dvd_trans hqc (dvd_self_pow c hp))
      (Int.dvd_trans hqa (dvd_self_pow a hp))
  have hqb := prime_dvd_of_dvd_pow hq hb
  have : q ∣ Int.gcd a b := Nat.dvd_gcd (Int.ofNat_dvd_left.mp hqa) (Int.ofNat_dvd_left.mp hqb)
  rw [hcop] at this
  have := Nat.le_of_dvd (by decide) this
  have := hq.1; omega

/-- Int.gcd only sees absolute values. -/
theorem gcd_natAbs_congr {a b a' b' : Int} (ha : a'.natAbs = a.natAbs) (hb : b'.natAbs = b.natAbs) :
    Int.gcd a' b' = Int.gcd a b := by
  show Nat.gcd a'.natAbs b'.natAbs = Nat.gcd a.natAbs b.natAbs
  rw [ha, hb]

/-- 2d, with coprimality kept: re-signing keeps |a'| = |o1|, |b'| = |e|. -/
theorem normalize_core' {p : Nat} (hp : p % 2 = 1) {o1 e o2 : Int}
    (h1 : o1 % 2 = 1) (he : e % 2 = 0) (h : o1 ^ p + e ^ p + o2 ^ p = 0) :
    ∃ a' b' c' : Int, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0
      ∧ (a' * b' * c').natAbs = (o1 * e * o2).natAbs
      ∧ a'.natAbs = o1.natAbs ∧ b'.natAbs = e.natAbs := by
  rcases (show o1 % 4 = 3 ∨ o1 % 4 = 1 by omega) with h4 | h4
  · refine ⟨o1, e, -o2, ?_, h4, he, ?_, rfl, rfl⟩
    · rw [neg_pow_of_odd hp]; omega
    · simp [Int.natAbs_mul, Int.natAbs_neg]
  · refine ⟨-o1, -e, o2, ?_, by omega, by omega, ?_, Int.natAbs_neg o1, Int.natAbs_neg e⟩
    · rw [neg_pow_of_odd hp, neg_pow_of_odd hp]; omega
    · simp [Int.natAbs_mul, Int.natAbs_neg]

/-- Frey normalization keeping gcd(a', b') = 1: the input the Frey curve needs. -/
theorem frey_normalize_cop {p : Nat} (hp : p % 2 = 1) (h5 : 5 ≤ p) {a b c : Int}
    (hcop : Int.gcd a b = 1) (habc : a * b * c ≠ 0) (h : a ^ p + b ^ p = c ^ p) :
    ∃ a' b' c' : Int, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0
      ∧ Int.gcd a' b' = 1 ∧ a' * b' * c' ≠ 0 ∧ (32 : Int) ∣ b' ^ p := by
  have hp1 : 1 ≤ p := by omega
  have ha0 : a ≠ 0 := fun e => habc (by rw [e]; simp)
  have hb0 : b ≠ 0 := fun e => habc (by rw [e]; simp)
  have hac : Int.gcd a c = 1 := cop_left ha0 hcop h hp1
  have hbc : Int.gcd b c = 1 :=
    cop_left hb0 (by rw [Int.gcd_comm]; exact hcop) (by omega) hp1
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
      ∧ (a' * b' * c').natAbs = (a * b * c).natAbs ∧ Int.gcd a' b' = 1 by
    obtain ⟨a', b', c', h1, h2, h3, h4, hg⟩ := H
    refine ⟨a', b', c', h1, h2, h3, hg, ?_, two_pow_dvd h5 h3⟩
    intro h0
    rw [h0] at h4
    exact habc (Int.natAbs_eq_zero.mp h4.symm)
  rcases (show (a % 2 = 0 ∧ b % 2 = 1 ∧ c % 2 = 1) ∨ (a % 2 = 1 ∧ b % 2 = 0)
      ∨ (a % 2 = 1 ∧ b % 2 = 1 ∧ c % 2 = 0) by omega) with ⟨h1, h2, h3⟩ | ⟨h1, h2⟩ | ⟨h1, h2, h3⟩
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core' hp h2 h1 (o2 := -c) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]; grind
    · rw [gcd_natAbs_congr e5 e6, Int.gcd_comm]; exact hcop
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core' hp h1 h2 (o2 := -c) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]
    · rw [gcd_natAbs_congr e5 e6]; exact hcop
  · obtain ⟨a', b', c', e1, e2, e3, e4, e5, e6⟩ :=
      normalize_core' hp h1 (show (-c) % 2 = 0 by omega) (o2 := b) (by rw [hnc]; omega)
    refine ⟨a', b', c', e1, e2, e3, ?_, ?_⟩
    · rw [e4]; simp [Int.natAbs_mul, Int.natAbs_neg]; grind
    · rw [gcd_natAbs_congr e5 (e6.trans (Int.natAbs_neg c))]; exact hac

/-- The one deep statement: no normalized coprime Fermat solution for a prime p ≥ 5
    (equivalently, the semistable Frey curve with these numbers cannot exist). -/
def NoFreyCurve : Prop :=
  ∀ (p : Nat), IsPrime p → 5 ≤ p → ∀ a b c : Int,
    Int.gcd a b = 1 → a * b * c ≠ 0 → a % 4 = 3 → b % 2 = 0 → (32 : Int) ∣ b ^ p →
    a ^ p + b ^ p ≠ c ^ p

/-- Fermat's Last Theorem for every n ≥ 3, given only `NoFreyCurve`. -/
theorem flt_of_no_frey_curve (H : NoFreyCurve) : ∀ n, 3 ≤ n → FLT n := by
  apply flt_of_prime_case
  intro p hp h5
  apply flt_of_coprime (by omega)
  intro x y z hx hy hcop hxyz
  have hodd : p % 2 = 1 := prime_ne_two_odd hp (by omega)
  have hI : (x : Int) ^ p + (y : Int) ^ p = (z : Int) ^ p := by exact_mod_cast hxyz
  have hg : Int.gcd (x : Int) (y : Int) = 1 := by
    rw [Int.gcd_natCast_natCast]; exact hcop
  have hz : (z : Int) ≠ 0 := by
    intro e
    have hz0 : z = 0 := by exact_mod_cast e
    subst hz0
    rw [Nat.zero_pow (by omega)] at hxyz
    have : x ^ p = 0 := by omega
    exact hx (Nat.pow_eq_zero.mp this).1
  have habc : (x : Int) * (y : Int) * (z : Int) ≠ 0 := by
    have hx' : (x : Int) ≠ 0 := by exact_mod_cast hx
    have hy' : (y : Int) ≠ 0 := by exact_mod_cast hy
    exact Int.mul_ne_zero (Int.mul_ne_zero hx' hy') hz
  obtain ⟨a', b', c', he, ha4, hb2, hg', hne, h32⟩ := frey_normalize_cop hodd h5 hg habc hI
  exact H p hp h5 a' b' c' hg' hne ha4 hb2 h32 he

end FermatRoute

