# What CI checks, and what it doesn't

Audit of Nov 25, 2026, done when everything was put on main. [CI](.github/workflows/ci.yml) runs on every push and pull request to main.

## Checked on every push

- **Core Lean library** (`lean-library/`): every module builds with core Lean 4; every theorem uses only the three standard axioms and no `sorry`; ASSUMPTIONS.md and theorems.tsv are current.
- **Every topic's own Lean file** (`tools/topic_lean.py`): each `topics/**/*.lean` is compiled on its own against the library, with sibling modules compiled first. 149 of 161 build. The 12 that don't are old first drafts whose compiled versions are in the library. They are listed with reasons in [tools/topic_lean_known.tsv](tools/topic_lean_known.tsv), and any other failure turns CI red.
- **Mathlib version** (`lean-library-mathlib/`, through [mathlib-build.yml](.github/workflows/mathlib-build.yml)): the port of the original 62 modules and the proof-complexity study build with Mathlib's prebuilt cache; standard axioms only, no `sorry`; every original theorem has a row in the port's name tables (`tools/port_notes.py`); `tools/CheckNames.lean` reports any "ported" row whose name doesn't resolve.
- **Mathlib map** (`topics/mathlib-map/`, same workflow, same pinned Mathlib tag `v4.37.0-rc1`): its fifteen modules build (`lake build`) and `Check.lean` prints the axioms of all 66 theorems; the job fails if the build reports a `sorry` or `sorryAx` appears. `LeeYangBand.lean`, `LeeYangStep.lean` and `LeeYangP3.lean` hold definitions only, so there is nothing in them to prove; `LeeYangP3Constants.lean` and `LeeYangCore.lean` are real proofs. `LeeYangNearZone.lean` is a real proof too, with the far-zone and near-zone computer checks as named hypotheses. `FaceDial.lean` (shape-dials, DRAFT) was written without a Lean toolchain, staged, and compiled clean in CI on Nov 26 with no proof changed; it now gates like the rest. Section 5 of `LeeYangNearZone.lean` (out to p = 8) was staged the same way, built clean, and now gates. `DepthDerived.lean` (shape-dials section 6) was staged the same way, built clean on Nov 26 after two fixes (a Mathlib lemma name and a closing `ring`, no statement changed), and now gates. `TopologyShapes.lean` (shape-dials section 7, topological shapes) was staged the same way, built clean on Nov 26 with no proof changed, and now gates. `KleinFlow.lean` (Klein-bottle flow) was staged the same way, built clean on Nov 26 with no proof changed, and now gates. No staged drafts remain.
- **Numeric checks** (`tools/run_checks.py`): every `topics/**/checks*.py`, plus the 65 check scripts with other names listed in [tools/extra_checks.txt](tools/extra_checks.txt), each compared with its saved output.

## Not checked yet

- **`topics/tesseract/TesseractCoords.lean`** (imports all of Mathlib) is built by the Mathlib job, and a failure turns CI red. It was report-only until it built cleanly after the 1.0 rename (CI, Nov 25).
- **Proof complexity's `Measure.lean` and `Tree.lean`** aren't built, because compiling them writes `yard.jsonl` and `tree.tsv`. They only measure; the proofs they measure are built.
- **The port still uses the pre-1.0 names** (`MultAdd`, `multAdd`, `wAdd`) where the core library now says `CircleOp`, `circ` and `circK`. The port builds and its name tables are checked against the library as it was when ported (`lean-library-mathlib/tools/original-theorems.tsv`). Renaming the port is a separate job.
- **Mathlib's declaration linter** (`lake exe runLinter FusedAlgebra`) passed when the port was written but isn't run in CI.
- **Scripts left out on purpose:** `topics/proof-engine/film.py` (needs a trace file from a Lean run) and every `figures.py`, `page.py` and `build.py` (they redraw pictures and pages, so there is nothing to compare).
- **Folders with results but no check script:**
  - `topics/dependency-map/`, `topics/dial-potential/` (the dark-energy search), `topics/proof-gallery/`: analysis and data scripts that print results but don't test them.
  - `topics/three-body-word-formula/`: long orbit integrations whose results are saved as JSON, with no check script yet.
  - `topics/jet-fields/`, `topics/shape-dependencies/`: their theorems are in the Lean library, which CI builds, but their numbers have no check script.
  - `topics/flt-tracking/`, `topics/mathlib-port-check/`, `topics/mathlib-plan/`, `topics/names-review/`, `topics/p-vs-np/`, `topics/paper-candidates/`: write-ups and reviews with nothing to compute.
- **Pages and figures** are not rebuilt in CI.
