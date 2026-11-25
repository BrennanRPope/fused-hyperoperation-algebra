# FusedAlgebra.Topics.FractionalRungs

Theorems: 8. SAME 5, GENERAL 2, DROPPED-OK 1. No findings.

- `dial_or` → `MultAdd.wAdd_neg_one`, an identical statement.
- Rung index check: the port adds `fuse_hyperoperation : fuse (hyperoperation n) = rung n`, proved by `rfl`. With mathlib's indexing (0 = succ, 1 = add, 2 = mul, 3 = pow, 4 = tetration), `rung 1` = `sadd`, `rung 2` = mult-add, `rung 3` = `fexp`, `rung 4` = `ftet`. That matches the original Rungs.lean numbering, and `fuse_succ`/`fuse_add`/`fuse_mul` agree with it.

Route changes: `fuse_succ` is now `rfl`; otherwise none.

New port theorems with no original: 1 (`fuse_hyperoperation`).
