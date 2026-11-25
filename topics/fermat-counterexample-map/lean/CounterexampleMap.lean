/-
  CounterexampleMap: what a Fermat counterexample would force, where the forced facts clash,
  and which methods cannot work
  ==========================================================================================

  Core Lean 4 (4.37.0-rc1), no mathlib, no `sorry`, no new axioms. Builds on FreyPieces and
  FermatRoute (same folder in the library). Deep theorems appear only as fields of the
  structure `DeepArrows` (named hypotheses), never as axioms.

  STATUS: written in this thread WITHOUT a Lean compiler (none installed here). Every
  identity was checked in Python/sympy and every numeric claim in PARI (see checks.py and
  checks_output.txt next to this folder). The library thread should compile it as
  `FusedAlgebra.Topics.Fermat.CounterexampleMap` and fix any step Lean rejects.

  A hypothetical counterexample a^p + b^p = c^p, p ≥ 5 prime, normalised by FermatRoute so
  that gcd(a, b) = 1, a ≡ 3 mod 4, b even, 32 ∣ b^p. Its Frey curve y² = x(x − a^p)(x + b^p).

  What is new here (FreyPieces already has the 2-torsion, the minimal model, c4, 256·Δ = (ABC)²,
  semistability, the node and its mult-add chord law, and v_q(Δ) ≡ 0 mod p for odd q):

  A. c6 of the minimal model, the identity c4³ − c6² = 1728·Δ, and gcd(c4³, Δ) = 1, so the
     j-invariant j = c4³/Δ has exactly Δ as its denominator (every bad prime is a pole of j).
  B. The bad primes are exactly the primes dividing ABC, so the conductor is rad(abc).
  C. The exponent of 2: Δ = 2^(2pv − 8)·(odd) where v = v₂(b), and p never divides 2pv − 8.
     So 2 is the one prime that stays in the level of the mod-p representation (level 2).
  D. The genus of X0(N) by the standard formula, computed for N = 1..13: zero at every level
     up to 10, at 12 and at 13, one at 11.
  E. The squeeze: Wiles (modularity), Mazur (irreducibility) and Ribet (level lowering to 2),
     with dim S2(Γ0(N)) = genus (Shimura) as the fourth cited arrow, give a contradiction.
     Combined with FermatRoute, the four arrows give FLT for every n ≥ 3.
  F. Methods that cannot work, each with a checkable witness:
     F1. No congruence test: for every modulus n a fake triple passes every normalisation.
     F2. Sophie Germain's case I: mod 11, fifth powers force 11 ∣ xyz (pushes divisibility only).
     F3. Size: a Fermat triple has abc < c³, so an explicit abc bound c^p < (abc)² gives p < 6.
     F4. Lamé's gap: in Z[√−5] no element has norm 2 or 3, yet 2·3 = (1+√−5)(1−√−5) = 6.
     F5. Kummer's gap: 37 divides the numerator of B_32 (37 is irregular); the Wieferich
         primes 1093 and 3511 pass the 2^(p−1) ≡ 1 mod p² test.
     F6. Near misses that look like solutions, and the small modulus that kills each.
-/
import FusedAlgebra.Topics.Fermat.FermatRoute

namespace CounterexampleMap
open CircleOp
open CircleOp.Cube (IsPrime)
open FreyPieces

/-! ## A. c6, the c4³ − c6² = 1728Δ identity, and the j-invariant's denominator -/

/-- c6 = −b2³ + 36·b2·b4 − 216·b6. -/
def c6 (E : Weier) : Int := -(E.b2 ^ 3) + 36 * E.b2 * E.b4 - 216 * E.b6

/-- A1. 2·c6 of the Frey model is (A − B)(A + 2B)(2A + B). -/
theorem frey_c6 {A B a2 a4 : Int} (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) :
    2 * c6 (freyModel a2 a4) = (A - B) * (A + 2 * B) * (2 * A + B) := by
  simp only [c6, Weier.b2, Weier.b4, Weier.b6, freyModel]
  have hb2 : (1 : Int) ^ 2 + 4 * a2 = B - A := by grind
  rw [hb2]
  grind

/-- A2. The polynomial identity behind j − 1728: 4c4³ − (2c6)² = 27(ABC)². -/
theorem frey_j_identity (A B : Int) :
    4 * (A ^ 2 + A * B + B ^ 2) ^ 3 - ((A - B) * (A + 2 * B) * (2 * A + B)) ^ 2
      = 27 * (A * B * (A + B)) ^ 2 := by
  grind

/-- A3. c4³ − c6² = 1728·Δ for the Frey model. -/
theorem frey_c4_c6_disc {A B C a2 a4 : Int} (hC : A + B = C) (h2 : 4 * a2 = B - A - 1)
    (h4 : 16 * a4 = -(A * B)) :
    (freyModel a2 a4).c4 ^ 3 - c6 (freyModel a2 a4) ^ 2 = 1728 * (freyModel a2 a4).disc := by
  have hc4 := frey_c4 h2 h4
  have hc6 := frey_c6 h2 h4
  have hD := frey_disc hC h2 h4
  have hj := frey_j_identity A B
  subst hC
  generalize (freyModel a2 a4).c4 = K at *
  generalize c6 (freyModel a2 a4) = c at *
  generalize (freyModel a2 a4).disc = D at *
  subst hc4
  have hsq : 4 * c ^ 2 = ((A - B) * (A + 2 * B) * (2 * A + B)) ^ 2 := by
    rw [← hc6]; grind
  grind

/-- A4. No prime dividing Δ divides c4³, so j = c4³/Δ is in lowest terms with denominator Δ:
    every bad prime is a pole of j (potentially multiplicative reduction, Tate's curve). -/
theorem frey_j_denominator {A B C a2 a4 : Int} (hC : A + B = C) (hcop : Int.gcd A B = 1)
    (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) {q : Nat} (hq : IsPrime q)
    (hd : (q : Int) ∣ (freyModel a2 a4).disc) : ¬ (q : Int) ∣ (freyModel a2 a4).c4 ^ 3 := by
  have hs := frey_semistable hC hcop h2 h4 hq hd
  generalize (freyModel a2 a4).c4 = K at hs ⊢
  intro h
  have e : K ^ 3 = K * (K * K) := by grind
  rw [e] at h
  rcases euclid_int hq h with h1 | h1
  · exact hs h1
  · rcases euclid_int hq h1 with h3 | h3 <;> exact hs h3

/-! ## B. The bad primes are exactly the primes of ABC: conductor = rad(abc) -/

/-- An odd prime does not divide 256 (so it can be cancelled from 256·Δ). -/
theorem odd_prime_coprime_256 {q : Nat} (hq : IsPrime q) (h2 : q ≠ 2) : Nat.Coprime q 256 :=
  Nat.Coprime.pow_right 8 (odd_prime_coprime_two hq h2)

/-- B1. With A odd and 32 ∣ B, a prime divides Δ exactly when it divides ABC. Because the
    curve is semistable (FreyPieces 3e), its conductor is the product of these primes,
    rad(ABC) = rad(abc). -/
theorem frey_bad_primes {A B C a2 a4 : Int} (hC : A + B = C) (hA : A % 2 = 1) (hB : 32 ∣ B)
    (h2 : 4 * a2 = B - A - 1) (h4 : 16 * a4 = -(A * B)) {q : Nat} (hq : IsPrime q) :
    (q : Int) ∣ (freyModel a2 a4).disc ↔ (q : Int) ∣ A * B * C := by
  have hD := frey_disc hC h2 h4
  by_cases hq2 : q = 2
  · subst hq2
    constructor
    · intro _
      obtain ⟨m, rfl⟩ := hB
      exact ⟨A * (16 * m) * C, by grind⟩
    · intro _
      obtain ⟨_, h4d⟩ := frey_at_two hC hA hB h2 h4
      obtain ⟨k, hk⟩ := h4d
      exact ⟨2 * k, by rw [hk]; grind⟩
  · have hc := odd_prime_coprime_256 hq hq2
    constructor
    · intro h
      have h' : (q : Int) ∣ (A * B * C) * (A * B * C) := by
        rw [show (A * B * C) * (A * B * C) = (A * B * C) ^ 2 by grind, ← hD]
        exact Int.dvd_trans h ⟨256, by grind⟩
      rcases euclid_int hq h' with h1 | h1 <;> exact h1
    · intro h
      have h' : (q : Int) ∣ 256 * (freyModel a2 a4).disc := by
        rw [hD]; exact Int.dvd_trans h ⟨A * B * C, by grind⟩
      rw [← Int.natAbs_dvd_natAbs, Int.natAbs_mul] at h'
      simp only [Int.natAbs_natCast] at h'
      have := Nat.Coprime.dvd_of_dvd_mul_left hc h'
      rw [← Int.natAbs_dvd_natAbs]; simpa using this

/-! ## C. The exponent of 2, and why 2 is the prime that stays in the level -/

/-- C1. If b = 2^v·m then ((a·b·c)²)^p = 256 · 2^(2pv − 8) · ((a·m·c)²)^p, once 2pv ≥ 8.
    With a, m, c odd this says the minimal discriminant (256·Δ = ((abc)²)^p) has exactly
    2^(2pv − 8) as its 2-part. -/
theorem disc_two_part (a m c : Int) {v p : Nat} (h8 : 8 ≤ 2 * p * v) :
    ((a * (2 ^ v * m) * c) ^ 2) ^ p = 256 * (2 ^ (2 * p * v - 8) * ((a * m * c) ^ 2) ^ p) := by
  have e1 : (a * (2 ^ v * m) * c) ^ 2 = (2 ^ v) ^ 2 * (a * m * c) ^ 2 := by grind
  have e2 : (((2 : Int) ^ v) ^ 2) ^ p = 2 ^ (2 * p * v) := by
    rw [← Int.pow_mul, ← Int.pow_mul]; congr 1; grind
  rw [e1, mul_pow', e2]
  generalize 2 * p * v = M at h8 ⊢
  obtain ⟨k, rfl⟩ : ∃ k, M = 8 + k := ⟨M - 8, by omega⟩
  rw [Nat.add_sub_cancel_left, pow_add']
  grind

/-- C2. For a prime p ≥ 5 and v ≥ 1, p never divides 2pv − 8: the exponent of 2 in the
    minimal discriminant is never a multiple of p, so the mod-p representation stays ramified
    at 2 and the level cannot drop below 2. (The odd primes all drop: FreyPieces 3g.) -/
theorem two_exponent_not_multiple {p v : Nat} (hp : IsPrime p) (h5 : 5 ≤ p) (hv : 1 ≤ v) :
    ¬ (p : Int) ∣ 2 * (p : Int) * v - 8 := by
  intro h
  have hpv : (p : Int) ∣ 2 * (p : Int) * v := ⟨2 * v, by grind⟩
  have h8 : (p : Int) ∣ 8 := by
    have := Int.dvd_sub hpv h
    rwa [show 2 * (p : Int) * v - (2 * (p : Int) * v - 8) = 8 by omega] at this
  have h8n : p ∣ 2 * (2 * 2) := by
    have : ((8 : Nat) : Int) = 8 := rfl
    rw [← this] at h8
    exact Int.natCast_dvd_natCast.mp h8
  rcases euclid hp h8n with h | h
  · have := Nat.le_of_dvd (by decide) h; omega
  · rcases euclid hp h with h | h <;> (have := Nat.le_of_dvd (by decide) h; omega)

/-! ## D. The genus of X0(N) by the standard formula -/

/-- φ(N), counted. -/
def phi (N : Nat) : Nat := (List.range N).countP (fun i => Nat.gcd i N = 1)
/-- Pairs (c, d) mod N with gcd(c, d, N) = 1. -/
def primPairs (N : Nat) : Nat :=
  (List.range (N * N)).countP (fun i => Nat.gcd (Nat.gcd (i / N) (i % N)) N = 1)
/-- μ(N) = [SL2(Z) : Γ0(N)] = #P¹(Z/N) = primPairs N / φ(N). -/
def mu (N : Nat) : Nat := primPairs N / phi N
/-- ν2(N) = #{x mod N : x² + 1 ≡ 0}. -/
def nu2 (N : Nat) : Nat := (List.range N).countP (fun x => (x * x + 1) % N = 0)
/-- ν3(N) = #{x mod N : x² + x + 1 ≡ 0}. -/
def nu3 (N : Nat) : Nat := (List.range N).countP (fun x => (x * x + x + 1) % N = 0)
/-- ν∞(N) = Σ_{d ∣ N} φ(gcd(d, N/d)), the number of cusps. -/
def nuInf (N : Nat) : Nat :=
  ((List.range (N + 1)).filter (fun d => d ≠ 0 ∧ N % d = 0)).foldl
    (fun s d => s + phi (Nat.gcd d (N / d))) 0
/-- 12·genus by the standard formula 12g = 12 + μ − 3ν2 − 4ν3 − 6ν∞ (the formula is cited). -/
def g12 (N : Nat) : Int := 12 + (mu N : Int) - 3 * nu2 N - 4 * nu3 N - 6 * nuInf N

/-- D1. Level 2: μ = 3, ν2 = 1, ν3 = 0, ν∞ = 2, so 12g = 0. -/
theorem g12_two : mu 2 = 3 ∧ nu2 2 = 1 ∧ nu3 2 = 0 ∧ nuInf 2 = 2 ∧ g12 2 = 0 := by decide +kernel

/-- D2. 12·genus for N = 1, …, 13: zero except at 11 (genus 1, the first weight-2 cusp form). -/
theorem g12_table :
    (List.range' 1 13).map g12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0] := by decide +kernel

/-! ## E. The squeeze -/

/-- A normalised counterexample for a prime p ≥ 5 (the shape FermatRoute reduces FLT to). -/
def Counterexample (p : Nat) (a b c : Int) : Prop :=
  IsPrime p ∧ 5 ≤ p ∧ Int.gcd a b = 1 ∧ a * b * c ≠ 0 ∧ a % 4 = 3 ∧ b % 2 = 0 ∧
    (32 : Int) ∣ b ^ p ∧ a ^ p + b ^ p = c ^ p

/-- The four deep arrows, as named hypotheses. Everything else the squeeze needs is proved:
    semistability and conductor rad(abc) (FreyPieces 3e, B1), the odd exponents being multiples
    of p (FreyPieces 3g, which is what Ribet's theorem consumes), the exponent at 2 never being
    one (C2, which is why the level lands on 2 and not 1), and 12·g(X0(2)) = 0 (D1). -/
structure DeepArrows where
  /-- dim S2(Γ0(N)), the weight-2 cusp forms of level N. -/
  dimS2 : Nat → Nat
  /-- Shimura: dim S2(Γ0(N)) equals the genus of X0(N). -/
  dim_eq_genus : ∀ N, 1 ≤ N → 12 * (dimS2 N : Int) = g12 N
  /-- "The Frey curve of (p, a, b, c) is modular." -/
  Modular : Nat → Int → Int → Int → Prop
  /-- "Its mod-p Galois representation is irreducible." -/
  Irreducible : Nat → Int → Int → Int → Prop
  /-- Wiles, Taylor–Wiles (1995): every semistable elliptic curve over Q is modular. -/
  wiles : ∀ p a b c, Counterexample p a b c → Modular p a b c
  /-- Mazur (1978), as used by Serre: a semistable curve with all 2-torsion rational has
      irreducible mod-p representation for p ≥ 5. -/
  mazur : ∀ p a b c, Counterexample p a b c → Irreducible p a b c
  /-- Ribet (1990): irreducible + modular + unramified at every odd bad prime (exponents ≡ 0
      mod p) + finite at p ⇒ it comes from a weight-2 newform of level 2. -/
  ribet : ∀ p a b c, Counterexample p a b c → Modular p a b c → Irreducible p a b c →
    1 ≤ dimS2 2

/-- E1. The squeeze: the four arrows leave no room for a counterexample. -/
theorem squeeze (H : DeepArrows) : ∀ p a b c, ¬ Counterexample p a b c := by
  intro p a b c hc
  have h1 := H.ribet p a b c hc (H.wiles p a b c hc) (H.mazur p a b c hc)
  have h2 := H.dim_eq_genus 2 (by decide)
  have h3 : g12 2 = 0 := g12_two.2.2.2.2
  rw [h3] at h2
  omega

/-- E2. With FermatRoute, the four arrows give Fermat's Last Theorem for every n ≥ 3. -/
theorem flt_of_deep_arrows (H : DeepArrows) : ∀ n, 3 ≤ n → FLT n :=
  FermatRoute.flt_of_no_frey_curve
    (fun p hp h5 a b c hg hne ha hb h32 he => squeeze H p a b c ⟨hp, h5, hg, hne, ha, hb, h32, he⟩)

/-! ## F. Methods that cannot work -/

/-- F1. No congruence can rule a counterexample out: for every modulus n ≠ 0 the fake triple
    (−1, 32n, −1) has a ≡ 3 mod 4, 32 ∣ b, gcd(a, b) = 1, abc ≠ 0, and a^p + b^p ≡ c^p mod n. -/
theorem fake_passes_mod (n : Int) (hn : n ≠ 0) {p : Nat} (hp : 1 ≤ p) :
    (-1 : Int) % 4 = 3 ∧ (32 : Int) ∣ 32 * n ∧ Int.gcd (-1) (32 * n) = 1 ∧
      (-1) * (32 * n) * (-1) ≠ 0 ∧ n ∣ (-1) ^ p + (32 * n) ^ p - (-1) ^ p := by
  refine ⟨by decide, ⟨n, rfl⟩, by simp [Int.gcd], by omega, ?_⟩
  have e : (-1 : Int) ^ p + (32 * n) ^ p - (-1) ^ p = (32 * n) ^ p := by omega
  rw [e]
  exact Int.dvd_trans ⟨32, by grind⟩ (FermatRoute.dvd_self_pow (32 * n) hp)

/-- F2a. Fifth powers mod 11 are only 0, 1, 10. -/
theorem germain_residues :
    ∀ r, r < 11 → ∀ s, s < 11 → ∀ t, t < 11 →
      (r ^ 5 + s ^ 5) % 11 = t ^ 5 % 11 → r = 0 ∨ s = 0 ∨ t = 0 := by decide +kernel

/-- F2b. Sophie Germain's auxiliary prime for p = 5: x⁵ + y⁵ = z⁵ forces 11 ∣ x, y or z.
    This pushes divisibility around (it is how case I is proved) but can never contradict:
    case II, where p divides abc, escapes every such prime. -/
theorem germain_mod_11 (x y z : Nat) (h : x ^ 5 + y ^ 5 = z ^ 5) : 11 ∣ x ∨ 11 ∣ y ∨ 11 ∣ z := by
  have ex : x ^ 5 % 11 = (x % 11) ^ 5 % 11 := Nat.pow_mod x 5 11
  have ey : y ^ 5 % 11 = (y % 11) ^ 5 % 11 := Nat.pow_mod y 5 11
  have ez : z ^ 5 % 11 = (z % 11) ^ 5 % 11 := Nat.pow_mod z 5 11
  have e : ((x % 11) ^ 5 + (y % 11) ^ 5) % 11 = (z % 11) ^ 5 % 11 := by
    rw [Nat.add_mod, ← ex, ← ey, ← ez, ← Nat.add_mod, h]
  rcases germain_residues (x % 11) (Nat.mod_lt _ (by decide)) (y % 11) (Nat.mod_lt _ (by decide))
      (z % 11) (Nat.mod_lt _ (by decide)) e with h1 | h1 | h1
  · exact Or.inl (Nat.dvd_of_mod_eq_zero h1)
  · exact Or.inr (Or.inl (Nat.dvd_of_mod_eq_zero h1))
  · exact Or.inr (Or.inr (Nat.dvd_of_mod_eq_zero h1))

/-- F3a. For a + b = c with a, b ≥ 1: abc < c³ (so rad of a Fermat triple is below c³). -/
theorem abc_lt_cube {a b c : Nat} (ha : 1 ≤ a) (hb : 1 ≤ b) (hc : a + b = c) :
    a * b * c < c * c * c := by
  have h1 : a * b < c * c := Nat.mul_lt_mul_of_lt_of_le (by omega) (by omega) (by omega)
  exact Nat.mul_lt_mul_of_lt_of_le h1 (Nat.le_refl c) (by omega)

/-- F3b. An explicit abc bound of the shape c^p < (abc)² applied to (a^p, b^p, c^p) leaves only
    p < 6. (abc itself is not proved, so this is a conditional route, not a proof.) -/
theorem explicit_abc_kills {a b c p : Nat} (ha : 1 ≤ a) (hb : 1 ≤ b) (hc : a + b = c)
    (habc : c ^ p < (a * b * c) ^ 2) : p < 6 := by
  have h3 := abc_lt_cube ha hb hc
  have hsq : (a * b * c) ^ 2 < c ^ 6 := by
    have e1 : (a * b * c) ^ 2 = (a * b * c) * (a * b * c) := by grind
    have e2 : c ^ 6 = (c * c * c) * (c * c * c) := by grind
    rw [e1, e2]
    exact Nat.mul_lt_mul_of_lt_of_le h3 (Nat.le_of_lt h3) (by
      have : 1 ≤ c * c * c := Nat.mul_pos (Nat.mul_pos (by omega) (by omega)) (by omega)
      omega)
  apply Classical.byContradiction
  intro hp
  have : c ^ 6 ≤ c ^ p := Nat.pow_le_pow_right (by omega) (by omega)
  omega

/-- F4a. Squares mod 5 are 0, 1, 4. -/
theorem sq_mod_five (x : Int) : (x * x) % 5 = 0 ∨ (x * x) % 5 = 1 ∨ (x * x) % 5 = 4 := by
  rw [Int.mul_emod]
  have h0 : 0 ≤ x % 5 := Int.emod_nonneg x (by decide)
  have h1 : x % 5 < 5 := Int.emod_lt_of_pos x (by decide)
  generalize x % 5 = r at h0 h1 ⊢
  have : r = 0 ∨ r = 1 ∨ r = 2 ∨ r = 3 ∨ r = 4 := by omega
  rcases this with rfl | rfl | rfl | rfl | rfl <;> omega

/-- F4b. Lamé's gap (1847): in Z[√−5] nothing has norm x² + 5y² equal to 2 or 3, so 2 and 3
    are irreducible there, yet 6 = 2·3 has the second factorisation (1 + √−5)(1 − √−5).
    Unique factorisation, which elementary descent needs, fails; for Fermat it fails in
    Z[ζ_p] from p = 23 on (class number of Q(ζ23) is 3, cited). -/
theorem no_norm_two_three (x y : Int) : x ^ 2 + 5 * y ^ 2 ≠ 2 ∧ x ^ 2 + 5 * y ^ 2 ≠ 3 := by
  have key := sq_mod_five x
  have e : x ^ 2 + 5 * y ^ 2 = x * x + 5 * (y * y) := by grind
  rw [e]
  constructor <;> intro h <;> omega

/-- Multiplication in Z[√−5] on pairs (a, b) = a + b√−5. -/
def zmul (x y : Int × Int) : Int × Int := (x.1 * y.1 - 5 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)
/-- The norm a² + 5b². -/
def znorm (x : Int × Int) : Int := x.1 ^ 2 + 5 * x.2 ^ 2

/-- F4c. The norm is multiplicative. -/
theorem znorm_mul (x y : Int × Int) : znorm (zmul x y) = znorm x * znorm y := by
  simp only [znorm, zmul]; grind

/-- F4d. The two factorisations of 6. -/
theorem six_two_ways : zmul (2, 0) (3, 0) = (6, 0) ∧ zmul (1, 1) (1, -1) = (6, 0) := by decide +kernel

/-- F5a. Kummer's gap: 37 divides the numerator 7709321041217 of B_32 (B_32 = −7709321041217/510),
    so 37 is irregular and Kummer's 1850 proof does not cover it. -/
theorem irregular_37 : 7709321041217 % 37 = 0 := by decide +kernel

/-- F5b. The Wieferich primes: 2^(p−1) ≡ 1 mod p² for p = 1093 and 3511. Wieferich (1909)
    showed case I needs this; it holds for these two, so the criterion alone cannot finish. -/
theorem wieferich_1093_3511 : 2 ^ 1092 % (1093 * 1093) = 1 ∧ 2 ^ 3510 % (3511 * 3511) = 1 := by decide +kernel

/-- F6. Two famous near misses, each killed by a tiny modulus: 3 for the first (3 divides 3987
    and 4365 but not 4472), 2 for the second (even + odd cannot be even). -/
theorem near_misses :
    3987 ^ 12 + 4365 ^ 12 ≠ 4472 ^ 12 ∧ (3987 ^ 12 + 4365 ^ 12) % 3 = 0 ∧ 4472 ^ 12 % 3 = 1 ∧
    1782 ^ 12 + 1841 ^ 12 ≠ 1922 ^ 12 ∧ (1782 ^ 12 + 1841 ^ 12) % 2 = 1 ∧ 1922 ^ 12 % 2 = 0 := by decide +kernel

end CounterexampleMap
