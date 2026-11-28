# DRAFT — Round 15: the coin-rule chain as a Lean file

**Result:** `CoinRule.lean` states and proves, for Lean 4 with Mathlib: the law's lift, that the law is commutative and associative, that Tsallis entropy composes independent parts by the law with k = 1 − q, and that every Tsallis entropy with q ≠ 1 breaks the coin rule. **It has not been compiled.** Installing Lean was refused in the workspace where it was written, so no machine has checked the Lean. The same five steps pass a symbolic check in sympy (`checks.py`, 5 of 5).

![chain](lp1_chain.png)

## To check it on your own computer
1. Install Lean with elan (see the Lean community's install page).
2. `lake new coinrule math` makes a project with Mathlib, and `lake exe cache get` fetches Mathlib's prebuilt files.
3. Copy `CoinRule.lean` into the project and run `lake env lean CoinRule.lean`. No output means every proof checked.

## Ours vs known
- Known: Tsallis pseudo-additivity (Tsallis 1988), and Faddeev 1956 (honest-weight grouping gives Shannon). Faddeev is cited, not proved here.
- Ours: the packaging as one short chain from the law to "only k = 0 passes the coin rule".

## Lead?
No new lead. It's a tool: once it compiles, the core claim of rounds 7 and 8 is machine-checked.
