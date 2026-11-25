# Brief: check that a mathlib port says the same thing as its core-Lean original

You are one of several parallel reviewers. Work read-only. Do NOT run git commands that change anything, do NOT push, do NOT install anything (no elan, lake, pip), do NOT call any mcp__hearthbot__ tool, and do NOT message anyone. There is no Lean on this machine: judge by careful reading and reasoning. Only write files under /home/user/portcheck/results/.

## Where things are
Repo checkout (branch mathlib-port): /home/claude/fused-hyperoperation-algebra
- Originals (core Lean 4, no mathlib): lean-library/FusedAlgebra/...   (Base.lean is the original of the port's Basic.lean)
- Port (mathlib): lean-library-mathlib/FusedAlgebra/...
- Porter's name tables, one per original module: lean-library-mathlib/port-notes/<Module>.tsv with columns original / status / new_name / note (status: ported, mathlib, dropped). Most modules also have port-notes/<Module>.md with the porter's remarks.
- Full list of the 1,750 original theorems: lean-library/theorems.tsv (columns: module, full theorem name, axioms). Every original THEOREM of your modules must get a verdict row (definitions too, if the tsv table lists them, but focus on theorems).
- Overview of the port: lean-library-mathlib/PORT-NOTES.md. Core definitions of the port are in lean-library-mathlib/FusedAlgebra/Basic.lean (mult-add ✳ on PreQuasiregular, wAdd, wpow, mpow, multAddList, msub, rung on mathlib hyperoperation). Originals' core is lean-library/FusedAlgebra/Base.lean and Rungs.lean.
The porter's notes are claims to check, not truth. The point of this review is that nobody has compared the restated statements line by line.

## What to do for each of your modules
For every original theorem (and every definition the theorems depend on), find its port counterpart via the table, then read BOTH statements in the source, including the `variable` lines and typeclass assumptions in scope, and the definitions they use. Decide one verdict:

- SAME: the same statement up to renaming/notation (⊞→✳, x*x→x^2, foldMA→multAddList, etc.).
- EQUIV: an equivalent restatement (say why in a few words: e.g. "u⁻¹ instead of hypothesis x*u=1, equivalent since u is forced to be x⁻¹", "conjunction split into three lemmas, all present").
- GENERAL: the port is strictly more general and the original follows from it by specialising (ℚ→any field, ℕ→semiring) with no extra hypothesis that is not automatic in the original's type. Also use GENERAL for a strengthened statement (hypothesis removed). Say what got more general.
- MATHLIB: replaced by a mathlib declaration. Check (from your knowledge of mathlib) that the named declaration plausibly states the original; say "unsure" if you can't tell.
- DROPPED-OK: dropped, and the named replacement really covers it.
- WEAKER: the port proves less (a hypothesis added that is not automatic in the original's setting, a conclusion part missing, a narrower type, a bound loosened, numbers changed).
- DIFFERENT: the port proves a different claim.
- MISSING: no counterpart found, or the named counterpart does not exist in the port.
- VACUOUS: the port's statement is trivial or emptied out (hypotheses can't hold, conclusion is True, a definition changed so the statement no longer means what the original meant, proof just restates an assumed hypothesis that carries the original's content).

Be strict about definitions: if a definition changed (e.g. list recursion became a Finset sum, a record became a mathlib Matrix, pairs of rationals became ℂ, indices shifted by one, a natural subtraction appeared), check that theorems using it still mean the same; a definition that changed meaning is a DIFFERENT/VACUOUS finding on every theorem that relies on it. For numeric `decide`/`norm_num` checks, check that the numbers are the same. For generalisations that added a hypothesis "automatic in the original's type" (e.g. nonnegativity, which is automatic over ℕ but not over ℚ), check that it really is automatic in the original; if the original was over ℚ and the hypothesis is new, that is WEAKER.

Also note, separately from the verdict, when the PROOF ROUTE changed in a way the user would care about (e.g. the original proved something by hand, the port cites mathlib for it), because the user cares that his own proofs are kept. Put that in the reason as "route: ...".

Give a confidence for each verdict: high / medium / low. Mark low when you'd want a Lean check, and say in one line what lemma would settle it.

## Output (write exactly these files, one pair per module)
1. /home/user/portcheck/results/<Module>.tsv, where <Module> is the module name as in the port-notes file name (e.g. FusedAlgebra.Topics.Finance). Tab-separated, header line:
   original<TAB>new_name<TAB>verdict<TAB>confidence<TAB>reason
   One row per original theorem (use the short or full original name as in theorems.tsv) and one row per definition whose meaning you checked. Keep reasons under ~200 characters, no tabs or newlines inside fields.
2. /home/user/portcheck/results/<Module>.md: a few lines: counts per verdict (theorems only), then every WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence item with the original statement and the port statement quoted (short), and why. Also list any route changes worth knowing. Also list port theorems that have no original (new bridge lemmas) only by count.

When done with all your modules, reply with a short summary: per module, counts per verdict, and the real findings (WEAKER/DIFFERENT/MISSING/VACUOUS) worst first, each in one line. Be honest; if a table row's claim is wrong, say so. Don't pad.
