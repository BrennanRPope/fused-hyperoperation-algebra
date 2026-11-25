# FusedAlgebra.Topics.QAcross.QDial

Theorems (16): SAME 6, GENERAL 10.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items. Definition `qnum` changed from the recursion to a geometric Finset sum; same values (qnum_succ recovers the recursion), so nothing downstream changes meaning.

Route changes:
- `qnum_lift` (now `one_add_mul_qnum`): the original proved it through the mult-add route (`qnum_eq_wpow` + `wpow_lift`); the port cites mathlib's `Finset.mul_geom_sum`.
- `qnum_add_std` (now `qnum_add_eq_add_pow_mul`): the original derived the standard rule from the mult-add rule `qnum_add`; the port proves the standard rule directly and derives `qnum_add` from it (direction reversed).
- Dial values / truth tables: `decide` over ℤ became `ring`/`simp` over any ring.

New port theorems with no original: 0.
