# FusedAlgebra.Topics.Relativity

Theorems (30): GENERAL 30. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

ℚ is generalised to a field / comm ring / linearly ordered field. Five theorems take `[NeZero (2 : K)]` (`dop_back`, `dop_wall`, `rung_natural_q`, `rung_parabola`, `rung_raw_walls`). That is automatic over ℚ and genuinely needed in characteristic 2.

The three removed hypotheses were each checked to be unneeded (Lean's x/0 = 0 makes both sides agree):
- `velC_units` without c ≠ 0: at c = 0 both sides are 0.
- `landscape_walls` without c ≠ 0: at c = 0 both sides are 1.
- `kick_step` without kD + u·kN ≠ 0: both sides are (N+uD)/(D+uN).

Definition changes checked:
- `B2` → `!![1,u;u,1]`.
- `kN`/`kD`: the mutual recursion → entries (0,1)/(0,0) of `boost u ^ n`. By hand: kN_succ = kN + u·kD and kD_succ = kD + u·kN, with bases 0 and 1, which matches.
- `kicks` → `Function.iterate`, with the same recursion. No original theorem uses kicks.

Route changes: none of note. The proofs are field_simp/ring instead of grind.

New port theorems with no original: 8 (these include the kN/kD/kicks equation lemmas).
