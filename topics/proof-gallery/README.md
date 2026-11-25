# Proof shape gallery

Every theorem in ../lean-library (440) drawn as a mult-add shape and grouped by
shape alone. Page: https://claude.ai/artifact/NWREn3Fbw1i5reUMrYaAXZ (private).

Run in order (Python 3, no packages; the Lean files are only read):
  python3 parse.py        # theorems.json: statement, proof, line for each theorem
  python3 features.py     # features.json: body (degree.py), frame, moves
  python3 cluster.py      # gallery.json: families, rung height, twins, ambiguity
  python3 build_data.py   # site/data.json for site/index.html

Shape = size of the statement (polynomial degree with library defs expanded) +
frame (=, <=, !=, divides, exists, iff) + the proof engine moves read from the
proof text, unfolded defs and cited lemmas by a fixed rule table (features.py).
Family = size + leading move (priority order in cluster.py). Fields never used.
Placed on the rung sheet (rung-shape pick, Nov 25) at the highest rung each
statement uses. Lean was not run: move detection is by lemma names and patterns.
Known weak spot: 154 proofs closed by grind/omega/simp show no move in their own
text; they are dashed on the page.
