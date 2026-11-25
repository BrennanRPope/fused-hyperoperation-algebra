# Fermat counterexample map

What a hypothetical counterexample a^p + b^p = c^p (p ≥ 5 prime) would force, where the
forced facts collide, and which methods are known not to work. Asked for Nov 25 as a
counterexample hunt "for method": build up what can't work.

## Files

- `lean/CounterexampleMap.lean`: 24 theorems in core Lean 4 (no mathlib, no `sorry`; the four
  deep theorems are fields of the structure `DeepArrows`, never axioms). Meant to live at
  `FusedAlgebra/Topics/Fermat/CounterexampleMap.lean`, importing `FermatRoute`.
  **Not yet compiled**: this thread had no Lean toolchain. Every identity was checked in
  sympy and every number in PARI/Python (below). The library thread should add it to
  `build.sh` after `FermatRoute`, compile, and fix any step Lean rejects.
- `checks.py` / `checks_output.txt`: 11 checks (sympy + PARI via cypari2), about 20 s.

## The forced properties, elementary to deep

| # | Property of the Frey curve y² = x(x − a^p)(x + b^p) | Status |
|---|---|---|
| 1 | Normalisation: gcd 1, a ≡ 3 mod 4, b even, 32 ∣ b^p | proved (FreyPieces 2, FermatRoute) |
| 2 | Full rational 2-torsion | proved (FreyPieces 3a) |
| 3 | Δ = 16(abc)^{2p}; minimal model, Δ_min = (abc)^{2p}/2^8 | proved (FreyPieces 3b–3d) |
| 4 | c4 = A²+AB+B², 2c6 = (A−B)(A+2B)(2A+B), c4³ − c6² = 1728Δ | c4 proved earlier; c6 and identity new (A1–A3) |
| 5 | Semistable (nodes only, including at 2) | proved (FreyPieces 3e) |
| 6 | Bad primes = primes of abc, so conductor = rad(abc) | new (B1) |
| 7 | j = c4³/Δ in lowest terms: every bad prime is a pole of j | new (A4) |
| 8 | Exponent of every odd prime (p included) in Δ is a multiple of p | proved (FreyPieces 3g) |
| 9 | Exponent of 2 is 2p·v₂(b) − 8, never a multiple of p | new (C1, C2) |
| 10 | So the mod-p representation is unramified outside 2p, finite at p, ramified at 2: Serre level 2, weight 2 | consequence of 8–9 plus Tate curves (cited) |
| 11 | Irreducible mod p | **Mazur** (cited, `DeepArrows.mazur`) |
| 12 | Modular, of level rad(abc) | **Wiles, Taylor–Wiles** (cited, `DeepArrows.wiles`) |
| 13 | Comes from a weight-2 newform of level 2 | **Ribet** (cited, `DeepArrows.ribet`) |
| 14 | dim S2(Γ0(N)) = genus X0(N); genus formula | Shimura (cited); formula evaluated for N ≤ 13 (D1, D2), genus(X0(2)) = 0 |

**The squeeze (E1):** 12 + 13 say dim S2(Γ0(2)) ≥ 1; 14 says it is 0. `squeeze` proves that
the four cited arrows leave no counterexample; `flt_of_deep_arrows` chains it to FLT for every
n ≥ 3 via FermatRoute. Ribet's arrow is the one that consumes property 8 (and 9 says where it
stops); Wiles's arrow is what makes the curve modular at all; Mazur's arrow is the one real
curves escape through (checks.py section 10: the curve from −125 + 128 = 3 has every exponent
a multiple of 3 except at 3 itself, so its mod-3 level collapses, and there is no contradiction
only because it has a rational 3-isogeny).

## Methods that cannot work (section F of the Lean file, sections 4–8 of checks.py)

- **Congruences / local tests.** A fake triple (−1, 32n, −1) passes every normalisation and
  satisfies the equation mod any n (F1). Mod-q obstructions exist (q = 11, 41, 71, 101 for p = 5)
  but they only force q ∣ abc, and q-adic solutions with abc ≠ 0 always exist.
- **Sophie Germain / Legendre auxiliary primes.** Prove case I (p ∤ abc) prime by prime
  (F2: mod 11 for p = 5; checks list the first auxiliary prime for every p < 100); case II escapes.
- **Elementary descent / Lamé's factorisation.** Needs unique factorisation in Z[ζ_p], which fails
  from p = 23 (class number 3, cited). Toy witness: Z[√−5], no norm 2 or 3, 6 = 2·3 = (1+√−5)(1−√−5) (F4).
- **Kummer's regular primes.** Covers regular p only; 37 is the first irregular prime
  (37 ∣ numerator of B_32, F5a); irregular primes are infinite (Jensen 1915); it is not known whether
  regular ones are.
- **Wieferich-type criteria.** Case I needs 2^{p−1} ≡ 1 mod p²; 1093 and 3511 satisfy it (F5b), so the
  criterion alone cannot finish.
- **Size / abc.** A Fermat triple has abc < c³ (F3a), so its abc quality exceeds p/3; for p = 5 that
  is above 5/3, more than the record 1.6299 (Reyssat). An explicit bound c < rad(abc)² would kill
  every p ≥ 6 (F3b), but abc is unproved.
- **Near misses** (F6): 3987¹² + 4365¹² vs 4472¹² dies mod 3; 1782¹² + 1841¹² vs 1922¹² dies mod 2.
  Each fake dies to some modulus, no modulus kills them all.

## Mult-add

At each bad prime the curve's smooth points form the multiplicative group, whose chord law is
mult-add (FreyPieces 3f). The Frey condition "exponent ≡ 0 mod p" is the statement that Tate's
parameter is a p-th power up to a unit, i.e. a mult-add p-th power: this is Tate's uniformisation
restated, clean but not new. No other step gets a cleaner view from the rungs.
