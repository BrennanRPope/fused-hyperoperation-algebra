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
- `lean-library/` (grows as topics are proved; the copy here is a clean build of 3,272 theorems in 137 modules, synced Nov 26 ~01:40 UTC)
- `topics/lee-yang/`
- `topics/logic-sets/`
- `topics/measure-probability/`
- `topics/numerical-analysis/`
- `topics/patterns/`
- `topics/shape-directions/`

Full audit of the working folder against main, Nov 25 ~22:00 UTC: everything with content is now in the repo. No drafts are left: the one folder still being written during the audit, `topics/fractional-rungs/e-every-height/`, was finished before the push. README.md "Promising leads" tracks where each lead stands.

Nov 25 ~22:25 UTC: two folders saved mid-way as drafts at Brennan's ask, both merged into other streams: `topics/height-running-coupling/` and `topics/lucas-weights/`. Two more drafts at ~22:40: `topics/cube-weights/` and `topics/holographic-entropy-cube/` (working folder `exploration/`). Each README opens with a DRAFT banner, and INDEX.md lists them under Drafts.

Nov 25 ~23:00 UTC: `topics/height-flow/` added (renamed from height-entropy; merges the entropy family, the running coupling and the horn invariants). Its invariants part is DRAFT. `topics/height-running-coupling/` stays where it was and is reused by it. `topics/probability-entanglement/` added at ~23:20 and `topics/quantum-to-classical/` at ~23:30, both as DRAFTs. At ~23:55: `topics/precession/` (its Lean file compiled into the library), `topics/semiring-audit/` and `topics/number-theory-geometries/`, all DRAFTs. `topics/gravity-from-theory/` at ~00:10 Nov 26, DRAFT. At ~01:00 Nov 26: `topics/curved-space/` (with `radial-plane/`), `topics/newton-g/`, `topics/ns-same-shape/`, `topics/warp-bubble/` (its Lean file compiled into the library), `topics/beyond-the-cube/`, `topics/polar-coordinates/`, `topics/fusion-magnets/` and `topics/closed-leads-revisit/`, all DRAFTs; gravity-from-theory gained `velocity-lines/` and `river-rebuild/`; the library gained the unified rung-one semiring proof and its lineage file. At ~01:20: `topics/fractional-steps/`, `topics/parabolic-fingerprints/` and `topics/polar-recast/` (with `hec5/` and `literature/`), DRAFTs, and `gravity-from-theory/dimension-dial/`. `topics/hypersphere/` (with `orbits/` and `orbits/six/`) followed at ~01:40, DRAFT, its Lean file compiled into the library.

The Mathlib version of the Lean library (`lean-library-mathlib/`) was merged into main on Nov 25 at Brennan's ask ("I'd like to make sure that everything's on main"), and main's CI now builds it with the proof-complexity study on every push. No pull request to Mathlib is planned until the theory is more complete.

When syncing again, copy the current state of every folder and update INDEX.md.
