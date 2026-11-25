# CLAUDE.md

A fresh session should read [README.md](README.md) and then [INDEX.md](INDEX.md) before doing anything else, then [IN-PROGRESS.md](IN-PROGRESS.md).

Conventions:
- Brennan's rule: an outcome must be geometric, with a geometric underpinning we're sure of, before it counts.
- Be honest about what is already known math and what might be new.
- Lean is core Lean 4 only (no mathlib) unless a file says otherwise.
- One folder per topic under `topics/`. Keep scripts and their saved output together, and keep large generated data out of git (see DATA.md).
- Keep INDEX.md up to date when a topic or page is added.
- When a new result could help mathematicians (a formal check, a sharp open question, a measurable number, a useful negative result), add it to the "What's here that could help mathematicians" catalog in README.md, with its status and an honest novelty label. Leave out pure re-derivations.
- Keep the "Promising leads" section of README.md current as work reports: one or two lines per lead (what it is, status, link), known math kept apart from our claims.
- Design principle (Brennan, Nov 25, 2026): build on existing published work first and follow the literature's names, notation and conventions; label anything beyond it as new. See CONTRIBUTING.md.
- Every write-up has an "Assumptions" section (the hypotheses its theorems carry, any modeling assumptions, which claims are only checked numerically) and, next to it, a "Known vs new" section with sources. `python3 tools/assumption_sections.py` lists the write-ups that still lack an Assumptions section.
- Before adding a name or symbol, check NAMING-COLLISIONS.md and add a row for any new clash.
- CI (`.github/workflows/ci.yml`) rebuilds the Lean library, fails on any `sorry`, `axiom`, `opaque` or non-standard axiom, fails if ASSUMPTIONS.md is stale, and runs every `topics/**/checks*.py`. After changing the library, copy `lean-library/build/theorems.tsv` over `lean-library/theorems.tsv`, run `python3 tools/theorems_md.py` and `python3 tools/assumptions.py`, and commit the results. A check script whose expected results change needs its saved `*_output.txt` updated too.
