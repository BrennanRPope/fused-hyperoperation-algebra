# FusedAlgebra.Foundations.FusedPeano

Theorems (65): SAME 32, MATHLIB 19, EQUIV 7, DROPPED-OK 6, MISSING 1.

## MISSING
- `add_right_cancel'`. Original: `{a b c : N} (h : add a c = add b c) : a = b` on any model. The table names mathlib's `add_right_cancel`, but that lemma needs `IsRightCancelAdd N`. A CommSemiring obtained by `Equiv.commSemiring` does not provide it, and the port declares no such instance. The table's own note admits this ("no IsCancelAdd N instance is declared"). No port statement gives right cancellation of + on N. It is true and takes two lines from `ringEquivNat.injective` and `Nat.add_right_cancel`, but as ported it is not stated. To settle in Lean, check whether `example [FusedPeano N] (a b c : N) (h : a + c = b + c) : a = b := add_right_cancel h` elaborates (I expect it fails).

## Medium-confidence EQUIV (route change that matters)
The port builds + and * on N by transporting ℕ's structure along `equivNat`. The original defined them by choice from sadd and madd. The two are the same operations, because `succ_add_eq_sadd` (with succ_inj) and `madd_eq` (with cancellation in N) pin the transported + and * to the original definitions. Both are proved from `toN_sadd` / `toN_madd`, and those come from the axioms by induction.
- `toN_add`: toN (a+b) = toN a + toN b. In the port this is true by construction (map_add of ringEquivNat.symm), not a theorem about the sadd-derived addition. The content now lives in sadd_eq / succ_add_eq_sadd.
- `toN_mul`: the same situation, with madd_eq.
- The Peano recursions add_zero', add_succ', succ_add', add_assoc', mul_succ', mul_comm' and the others are now mathlib lemmas of the transported semiring. The original derived them from the fused axioms. That derivation is no longer in the port.

## Definitions
- `lad` → `MultAdd.rung`, checked for all n, a, b. Both give rung 0 = b+1, rung 1 a 0 = a+1 and rung (n+2) a 0 = a (via hyperoperation_ge_two_eq_self), and rung_succ_succ is lad's iteration law. So lad = rung, and toN_lad → toN_rung and lad_one/two/three → rung_one/two/three keep their meaning.
- `hyp`/`hypBase` → mathlib `hyperoperation`. The bases a, 0, 1 at n = 1, 2, ≥3 and the successor at n = 0 match exactly.
- `lad_lift` (DROPPED-OK, medium). In the port, rung is defined through the lift, so rung_add_one is only positivity. The original theorem "the schema ladder is the conjugated hyperoperation, for every n" is now split across rung_zero, rung_one_zero_right, rung_add_two_zero_right and rung_succ_succ. Taken together they say the same thing, but no single theorem states it.
- fromN → equivNat (the unique inverse of toN). recN and ladN are unchanged.

Other route notes: the axiom-level rung-1 and rung-2 proofs (sadd_comm, madd_sadd, madd_assoc, …) are kept verbatim. instCommSemiring is now a mathlib CommSemiring plus N ≃+* ℕ, which is stronger, but it is built by transport, not from the derived laws.

New port theorems with no original: 12, plus the 2 rung base lemmas referenced by the ladBase row.
