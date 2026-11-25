# Do shapes depend on shapes?

Question (Brennan): laws depend on other laws, and we group laws by shape. Do certain shapes depend on other shapes?

## Answer in three levels

1. **Proof level (objective).** In the compiled library (3,003 theorems, uses read from Lean's own proof terms), the proof bodies form a real order: stick → square → box / hyperbox → cube of any size. Only 27% of the 637 uses between different bodies run backwards; shuffling labels inside each Lean module gives 44% (z ≈ −11, 0 of 300 shuffles as low). Box and hyperbox are a near tie. The point body is too sparse to place. The lead moves order too: read-off facts at the bottom, dial and squash near the top (27% backwards against 42%). Caveat: much of this is expected, because big statements get proved from small lemmas. It says the library is built bottom-up; it is not a fact about nature.
2. **Families.** At the full 58-family level the graph tangles into cycles (mostly same-body families used both ways). Keep only links with 12 or more uses and it is a clean DAG of 20 edges. Main exception: cube of any size → stick, 43 uses, where a law proved for every n is used at one n. That is specialising, not building.
3. **Whole science laws (39 laws, 57 textbook links).** No order. 13 links go to a bigger body, 12 to a smaller one, 32 stay level; chance would give as many down as up. Most laws are square-dominated by weak majorities. The chain electromagnetism → special relativity → general relativity → Gravity Probe B goes hyperbox → square → square → point. The library's import order even reverses history (the EM file imports the relativity files).
4. **Exact level (composition shapes, the landscape's discriminant classes).** One clean line, proved in Lean (21 theorems, standard axioms only): addition is the root; the circle operation (called mult-add before 1.0) is addition plus the corner ab; the no-wall law is two walls turned by i; merging the two walls gives the border (addition); pinching the torus once gives the circle operation, twice gives addition. This is known math (formal group laws; singular cubics, Silverman III.2.5), now checked in this library's own terms.

## Lead verdict

The proof-shape grouping is not a lead about laws: its order is a fact about how proofs are written. The composition-shape grouping is the only notion here with an exact dependency order, so it is the one worth following: sort physical laws by their composition shape and check whether their textbook dependencies follow addition → walls → turned / torus. First look (my reading, not tested): Galilean velocity (border) → special relativity (two walls) → Thomas rotation (no wall); Shannon (border) → Tsallis (two walls). Expect much of it to reproduce known group contractions, so the test is whether anything shows up that the contraction story does not already predict.

## Tags on law links

Each of the 57 links in `laws.py` is tagged:
- `textbook`: standard textbook derivation.
- `reading`: my own reading (thermodynamics → information, GR → dark sector, EM → charge conservation).
- `lean`: the library's Lean files show the use directly (SR → Thomas, quantum noise → entropy, information → entropy).

## Files

- `data/gallery_all.json`: proof gallery rerun on the current library (3,026 entries).
- `data/extract.jsonl`: theorem → constants used, from the proof engine's Extract (3,003 theorems).
- `proof_uses.py` → `proof_uses.json`, `proof_uses_output.txt`: body flow and null, body matrix, family graph at thresholds, strongest edges, move matrix.
- `orders.py` → `orders.json`, `orders_output.txt`: best orders (exact search over subsets, fewest back uses), 300-draw nulls, pair counts, example lists.
- `laws.py` → `laws.json`, `laws_output.txt`: 39 laws, 57 tagged links, shape flow along links.
- `lean/ShapeDependencies.lean`: the 21 exact theorems. Compiled in the shared library as `FusedAlgebra/Topics/ShapeDependencies.lean` (imports `FusedAlgebra.Base` only; standard axioms only, no `sorry`; renamed for 1.0: `circ`, `circK`, dials k and d).
- `figures.py` → `figures/*.svg`: bodies, actions, families, laws, chain_physics, chain_exact, exception.
- `index.html`: the page.

## How to rerun

```
# Lean check (shared library build)
cd lean && LEAN_PATH=/mnt/project-files/math-theory/lean-library/build lean ShapeDependencies.lean

# Analysis
python3 -I proof_uses.py data/gallery_all.json data/extract.jsonl
python3 -I orders.py data/gallery_all.json data/extract.jsonl
python3 -I laws.py data/gallery_all.json
python3 -I figures.py
```

To regenerate the two data files: build a copy of the library with every module compiled (Electromagnetism and SizeDialResidual had no .olean when this was run), build the proof engine against it, run its Extract to get `extract.jsonl`, and run the gallery scripts. The gallery's `degree.py` needs a memo on its definition expansion (key `(id(defs), short, tuple(ds))`), or it runs for many minutes; `cluster.py` needs a fallback field for modules it does not list.

Nothing here was pushed to the repository.
