# Proof shape gallery

> **Names:** some words or symbols on this page (rung) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Every theorem in ../lean-library (1,485 on Nov 25) drawn as a circle-operation shape and
grouped by shape alone. Page: https://claude.ai/artifact/NWREn3Fbw1i5reUMrYaAXZ (private).

Run in order (Python 3, no packages):
  bash ../proof-engine/lean/build.sh   # builds the engine; Extract.lean writes
                                       # ../proof-engine/lean/build/extract.jsonl
  python3 -I parse.py        # theorems.json: statement, proof, line (plus any theorem
                             # the regex misses, taken from the Lean extract)
  python3 -I features.py     # features.json: body (degree.py), frame, moves from the
                             # proof text, and moves read from the compiled proof term
  python3 -I cluster.py      # gallery.json: families, rung height, twins, ambiguity
  python3 -I build_data.py   # site/data.json for site/index.html

Shape = size of the statement (polynomial degree with library defs expanded) +
frame (=, <=, !=, divides, exists, iff) + the proof engine moves. Moves come first
from the proof's own text (lemma names and patterns, features.py), then from the
compiled proof term: Extract.lean lists the constants each proof uses (opening
Lean's auxiliary constants), minus those in its statement, and TERM_RULES maps them
to moves. Family = size + leading move (priority order in cluster.py).
Placed on the rung sheet at the highest rung each statement uses.

Nov 25 run: 1,485 theorems, 53 families, 748 exact shapes. 315 proofs name no move
in their own text; the compiled term gives a move to all but 88 of them (Lean-made
structure fields and ext lemmas, rfl-like facts, and proofs that only chain earlier
theorems or split cases). 440 tiles took their lead move from the compiled term.
Earlier run (Lean not run): 440 theorems, 154 with no move.
