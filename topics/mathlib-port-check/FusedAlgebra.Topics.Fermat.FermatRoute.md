# FusedAlgebra.Topics.Fermat.FermatRoute — review

Theorems (9): EQUIV 3 (flt_of_prime_case, frey_normalize_cop, flt_of_no_frey_curve), MATHLIB 4 (prime_dvd_of_dvd_pow, gcd_one_of_no_prime, dvd_self_pow, gcd_natAbs_congr), GENERAL 1 (cop_left → FreyPieces.gcd_eq_one_of_add_pow_eq, hypothesis a ≠ 0 dropped), DROPPED-OK 1 (normalize_core' → strengthened FreyPieces.normalize_core). Definition NoFreyCurve: EQUIV (IsPrime → Nat.Prime only).

No WEAKER / DIFFERENT / MISSING / VACUOUS items.

Checked: FLT n (original: x, y ≠ 0, z may be 0) vs mathlib FermatLastTheoremFor n (also c ≠ 0). The exponent is unchanged and FreyPieces.fermatLastTheoremFor_iff proves the two equal for every n, so the switch is exact. IsPrime (2 ≤ p ∧ only divisors 1, p) is Nat.Prime. NoFreyCurve still uses Int.gcd a b = 1, unchanged.

Stays an assumption (in both original and port): NoFreyCurve (Wiles / Taylor–Wiles / Mazur / Ribet). Nothing deep is proved. Worth knowing: the port's new lemma noFreyCurve_of_flt shows NoFreyCurve is just a special case of FLT. So flt_of_no_frey_curve is an elementary reduction, not progress on the deep step. That was equally true of the original.

Route changes:
- flt_of_prime_case now gets exponents 3 and 4 from fermat4 / Fermat3.flt3_nat, which in the port are corollaries of mathlib (fermatLastTheoremFour / fermatLastTheoremThree), not the user's descents.
- cop_left is proved with IsCoprime algebra instead of the original's prime-factor argument (prime_dvd_of_dvd_pow + gcd_one_of_no_prime).

New port theorems with no original: 1 (noFreyCurve_of_flt).
