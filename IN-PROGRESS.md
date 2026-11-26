# Still in progress (last synced Nov 25, 2026, 08:38 UTC)

This repo is a snapshot. Work carries on in other sessions, so any folder may get newer results after the last sync. These were written to in the 20 minutes before it:

- `topics/biology/`
- `topics/category-theory/`
- `topics/chemistry/`
- `topics/combinatorics/`
- `topics/complex-analysis/`
- `topics/differential-geometry/`
- `topics/dynamical-systems/`
- `topics/functional-analysis/`
- `topics/game-theory/`
- `topics/graph-theory/`
- `topics/groups-lie/`
- `topics/information-theory-full/`
- `topics/knot-theory/`
- `lean-library/` (grows as topics are proved; the copy here is a clean build of 3,323 theorems in 141 modules, Nov 26 ~03:50 UTC)
- `topics/lee-yang/`
- `topics/logic-sets/`
- `topics/measure-probability/`
- `topics/numerical-analysis/`
- `topics/patterns/`
- `topics/shape-directions/`

Full audit of the working folder against main, Nov 25 ~22:00 UTC: everything with content is now in the repo. No drafts are left: the one folder still being written during the audit, `topics/fractional-rungs/e-every-height/`, was finished before the push. README.md "Promising leads" tracks where each lead stands.

Nov 25 ~22:25 UTC: two folders saved mid-way as drafts at Brennan's ask, both merged into other streams: `topics/height-running-coupling/` and `topics/lucas-weights/`. Two more drafts at ~22:40: `topics/cube-weights/` and `topics/holographic-entropy-cube/` (working folder `exploration/`). Each README opens with a DRAFT banner, and INDEX.md lists them under Drafts.

Nov 25 ~23:00 UTC: `topics/height-flow/` added (renamed from height-entropy; merges the entropy family, the running coupling and the horn invariants). Its invariants part is DRAFT. `topics/height-running-coupling/` stays where it was and is reused by it. `topics/probability-entanglement/` added at ~23:20 and `topics/quantum-to-classical/` at ~23:30, both as DRAFTs. At ~23:55: `topics/precession/` (its Lean file compiled into the library), `topics/semiring-audit/` and `topics/number-theory-geometries/`, all DRAFTs. `topics/gravity-from-theory/` at ~00:10 Nov 26, DRAFT. At ~01:00 Nov 26: `topics/curved-space/` (with `radial-plane/`), `topics/newton-g/`, `topics/ns-same-shape/`, `topics/warp-bubble/` (its Lean file compiled into the library), `topics/beyond-the-cube/`, `topics/polar-coordinates/`, `topics/fusion-magnets/` and `topics/closed-leads-revisit/`, all DRAFTs; gravity-from-theory gained `velocity-lines/` and `river-rebuild/`; the library gained the unified rung-one semiring proof and its lineage file. At ~01:20: `topics/fractional-steps/`, `topics/parabolic-fingerprints/` and `topics/polar-recast/` (with `hec5/` and `literature/`), DRAFTs, and `gravity-from-theory/dimension-dial/`. `topics/hypersphere/` (with `orbits/` and `orbits/six/`) followed at ~01:40, DRAFT, its Lean file compiled into the library. At ~01:55: `topics/three-body/` (libtb.so is built from tb.c on first use, not kept), `topics/quantum-action/` (working folder `action/`; Lean compiled), `topics/qsharp-functors/`, `topics/simplex/` (with section 7) and `topics/lucas-temperley-lieb/`, all DRAFTs. At ~02:05: hypersphere updated (orbits/dcopy and orbits/twov; its Lean file now 33 theorems, compiled), `topics/shape-dial-survey/` and `topics/exceptional-shapes/` (working folder for the missing shapes) and `topics/conics/`, DRAFTs. `topics/bridge-survey/` followed at ~02:20, DRAFT. At ~03:20: the polishing step (a cheap refinement after every number search, Brennan's ask): `topics/polish-audit/` (DRAFT) lists what changed, with `polish.py` and its output in each folder it touched. Three-body's 66 orbits now close to 26–32 digits and 6 orbits, of 6 loop types, are not in Li and Liao's list (was 9 orbits of 8 types); the dial-potential search had a bracket bug, and with it fixed the sombrero (hilltop) ties overall but fails the held-out test, so it is ruled out, while the dial on a circle (cosine) stays alive. At ~03:25: `topics/three-body-word-formula/` (DRAFT): the run-based period rule, the word-to-orbit recipe, and the hyperspherical and glued-sphere tests. Its `data/external/` downloads were still being written and are left out; the two compiled libraries are built from their .c files on first use. At ~03:40: the dial names from affine geometry: `topics/affine-names/` (DRAFT) and the renamed dial terms in the notes of 21 folders (the 1.0 names and the Lean names are kept). At ~04:00: the holographic third pass in `topics/holographic-entropy-cube/` (checks3, figures3, and three Lean files now compiled into the library as HoloGluing, HoloGluingField and FuzzyFrame, 19 theorems), and the dial design rule across the rename folders: base dials are settings (odds t, similarity ratio r, growth degree β), readings are outputs (holonomy ω).

The Mathlib version of the Lean library (`lean-library-mathlib/`) was merged into main on Nov 25 at Brennan's ask ("I'd like to make sure that everything's on main"), and main's CI now builds it with the proof-complexity study on every push. No pull request to Mathlib is planned until the theory is more complete.
The Mathlib map (`topics/mathlib-map/`) was added on Nov 26, ported to the 1.0 names, and wired into the Mathlib CI job on the pinned tag. It was first compiled against Mathlib master (commit `348837d`) because the cache server was out of reach; this machine can't reach the cache either, so the first CI run on the tag is its real build check.

When syncing again, copy the current state of every folder and update INDEX.md.
