# LogicSets.lean

- 35 theorems, namespace `LogicSets`, core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Imports `FusedAlgebra.Base` and `FusedAlgebra.Rungs` (for `fexp`, `sadd`, `fexp_tower`, `fexp_sadd`).
- Compiles clean: `LEAN_PATH=/tmp/lib/build lean lean/LogicSets.lean` (Nov 25, 2026). A `#print axioms`
  run on idem_defect, hamS_odds, hamT_mirror, multAdd_definable, fold_xor_trues, one_time_pad,
  bits_closed_iff, zero_definable showed only propext, Classical.choice, Quot.sound.
- Intended module path: `FusedAlgebra/Topics/LogicSets.lean` (add to build.sh and FusedAlgebra.lean).

Contents: dial table; idempotent defect (k+1)(k+2)ab; OR/XOR/AND keep idempotents; distributivity and
absorption defects and their vanishing on idempotents; OR from XOR+AND and back; lattice order; OR wall;
bits closed only at k = -1, -2; n trues fold to [n]_0 = 1 and [n]_{-1} = n mod 2; one-time pad (XOR wall 1/2);
Hamacher t-conorm as landscape law on p + q = -1; landmarks gamma = 0, 1, 2; miss formula; Hamacher = dial on
odds; De Morgan mirror to the t-norm; Hamacher on bits is OR; fexp at base 1 injective; mult-add, shifted add
and 0 definable from the fused exponent; lift iso for mult-add alone.
Not repeated (already in Topics/Statistics): truth tables, De Morgan for OR, inclusion–exclusion.
