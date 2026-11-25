# Porting guide: core-Lean `lean-library/` → mathlib form `lean-library-mathlib/`

The original library (`../lean-library/FusedAlgebra/…`) is core Lean 4 with no mathlib. This
project restates it in mathlib style, on Lean 4.37.0-rc1 with mathlib v4.37.0-rc1.

## Where things go

- Original `lean-library/FusedAlgebra/X/Y.lean` → port `lean-library-mathlib/FusedAlgebra/X/Y.lean`
  (same path; the one exception is `Base.lean` → `Basic.lean`).
- One name-mapping table per original module: `port-notes/FusedAlgebra.X.Y.tsv` (format below).
- Free-form remarks for a module (things that didn't fit, judgment calls):
  `port-notes/FusedAlgebra.X.Y.md` (optional, short).

## File shape

```lean
module

public import FusedAlgebra.Basic        -- and any other ported modules you need

/-!
# Title

What the file is about, the main definitions and results, in mathlib's module-docstring style.
Say plainly which results are already known mathematics (and where they are in mathlib).
-/

@[expose] public section

namespace OriginalNamespace

open MultAdd   -- brings the scoped `✳` notation into scope

…

end OriginalNamespace
```

`FusedAlgebra.Basic` already imports all of `Mathlib`, so no other mathlib imports are needed.
Keep the original namespace of each module so names map one to one.

## The core API (read `FusedAlgebra/Basic.lean`, `Rungs.lean`, `Cube.lean` first)

Mult-add is `a ✳ b = a + b + a * b` (`MultAdd.multAdd`, scoped notation, precedence 65, the same
as `+`, so write `1 + (a ✳ b)` with parentheses). It is mathlib's circle operation: multiplication
in `PreQuasiregular R`. The old names changed; the full tables are
`port-notes/FusedAlgebra.Base.tsv`, `…Rungs.tsv`, `…Cube.tsv`. The ones you will hit most:

| core Lean | mathlib port |
|---|---|
| `a ⊞ b` | `a ✳ b` |
| `multAdd_lift`, `multAdd_lift'` | `one_add_multAdd`, `multAdd_add_one` |
| `mpow_lift`, `wAdd_lift`, `wpow_lift` | `one_add_mpow`, `one_add_mul_wAdd`, `one_add_mul_wpow` |
| `foldMA`, `foldMA_lift` | `multAddList`, `one_add_multAddList` |
| `pw x n`, `prodL`, `sumL`, `groups` | `x ^ n`, `List.prod`, `List.sum`, `List.sublists'` |
| `Cube.layer n k`, `Cube.fact n` | `Nat.choose n k`, `Nat.factorial n` |
| `tet`, `pent` | still `tet`, `pent`, now `hyperoperation 4`, `hyperoperation 5` |
| `Nat`, `Int`, `Rat`, `Type` | `ℕ`, `ℤ`, `ℚ`, `Type*` |

## What "mathlib form" means here

1. **Build on mathlib's structures.** Replace hand-rolled definitions by mathlib's when one fits:
   `Finset.sum`/`∑`, `List.sum`, `Nat.choose`, `Nat.factorial`, `Nat.Prime`, `Nat.gcd`, `ZMod`,
   `Matrix`, `Polynomial`, `Finset.range`, `Function.iterate`, `deriv`, group/ring/field classes…
   Keep a custom definition only when mathlib has nothing equivalent, and give it a docstring.
2. **Use mathlib's theorems and tactics instead of re-proving basics.** `ring`, `field_simp`,
   `linarith`, `nlinarith`, `positivity`, `norm_num`, `omega`/`lia`, `simp`, `gcongr`, `decide`.
   Do not use `grind` with `Lean.Grind` classes (`open Lean.Grind`, `Semiring.pow_succ`, …).
3. **Drop helper lemmas that are one-liners from mathlib**, and record each in the mapping table
   with the mathlib name that replaces it.
4. **If a headline theorem is already in mathlib**, state it here as a short corollary of the
   mathlib theorem (so the name still exists), and say so in the docstring and the table.
5. **Generalise when it is free.** If a statement about `ℚ` holds with the same proof in any
   field, commutative ring or linear ordered field, state it there (`{K : Type*} [Field K]`,
   `[LinearOrder K] [IsStrictOrderedRing K]`). If generalising costs real work, keep the
   original type.
6. **Keep the mathematical content.** Every original theorem is either ported (possibly renamed
   or generalised), replaced by a named mathlib theorem, or listed as not ported with the reason.
   Never weaken a statement silently, never use `sorry`, never add axioms.

## Style (the build must be warning-free)

- Mathlib naming: `snake_case` theorems named after their statement (`one_add_multAdd`, not
  `multAdd_lift`), `lowerCamelCase` defs, `UpperCamelCase` structures and types. Rename where the
  original breaks this, and record the rename.
- Every `def`, `structure`, `class`, `instance` with a name, and every main theorem gets a
  `/-- docstring -/`. Section headers use `/-! ### … -/`.
- Lines at most 100 characters. `fun x ↦ …` (not `=>`). No `$`; use `<|`. Two-space indentation,
  `by` at the end of the line. Unicode `ℕ ℤ ℚ ℝ`.
- Use `@[simp]` sparingly, and only on lemmas whose left side is in simp normal form.
- No unused variables or arguments, no unused `have`s.

## Checking your work

```bash
cd lean-library-mathlib
export PATH="$HOME/.elan/bin:/usr/bin:$PATH"
lake env lean FusedAlgebra/Topics/X.lean      # must print nothing: no errors, no warnings
lake build FusedAlgebra.Topics.X              # only once the line above is clean
```

Each `lake env lean` run takes about 15 to 30 seconds (it loads mathlib). Never run
`lake update`, `lake clean`, `lake exe cache get` or a bare `lake build`.

## Mapping table format

`port-notes/FusedAlgebra.X.Y.tsv`, tab-separated, with this header line and one row for every
original `def`, `structure`, `instance` and `theorem` of the module, including helpers. Use the
original name without its namespace in the first column.

```
original	status	new_name	note
```

`status` is one of:

- `ported`: restated here; `new_name` is the full new name (same or renamed).
- `mathlib`: replaced by an existing mathlib declaration; `new_name` is that declaration.
- `dropped`: removed as redundant; `new_name` says what covers it, `note` says why.
- `not-ported`: did not fit; `note` says why.

The original theorem names of every module are listed in `../lean-library/theorems.tsv`
(module, name, axioms); your table must cover every name listed there for your module.
