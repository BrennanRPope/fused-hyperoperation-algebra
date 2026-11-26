"""The results on the dial map: one row per result, with the dials it uses, how it uses
the fused operation, its status, how new it is, and the folder it lives in.
Every row was read off that folder's README on Nov 26, 2026 (see README.md here)."""

DIALS = [  # letter, name, colour
    ("k", "homothety k", "#1f77b4"),
    ("d", "curvature dial d\n(reading: curvature −Δ/4)", "#d62728"),
    ("p", "ℓp exponent p", "#2ca02c"),
    ("C", "cubic form C", "#9467bd"),
    ("t", "odds t", "#ff7f0e"),
    ("r", "similarity ratio r", "#8c564b"),
    ("q", "order dial q (or ħ)", "#e377c2"),
]

# bin: "op"    = uses a ∘ₖ b = a + b + k·ab itself, or its cube of corners
#      "plane" = uses the plane law (a + b + k·ab)/(1 − d·ab) or the wider affine/Möbius composition
#      "pic"   = uses only the dial picture: the family of shapes or a reading such as Δ sets the
#                question; the proof is other mathematics
# status: "lean" proved in Lean; "checked" exact or numeric checks (plus prose proof where stated);
#         "partial" partly proved, rest open or numeric
# new: "known" known result; "reading" known result, the placement is ours; "new" not found in print
#      in the write-up's own search (not a guarantee)
RESULTS = [
    # --- uses the fused operation itself ---
    ("Fused rungs: every rung's laws, all rungs at once", "k", "op", "lean", "reading", "hyperoperation-ladder"),
    ("∘₁ is the categorical product of partial functions", "k", "op", "lean", "known", "category-theory"),
    ("Every dial k ≠ 0 is one monoid; walls = quasispectrum; Jacobson radical", "k", "op", "lean", "reading", "mathlib-map"),
    ("Sugeno λ-measure and α-determinantal 'no point' law", "kt", "op", "lean", "reading", "dial-point-processes"),
    ("Power set with odds; the unique half (odds √2 − 1)", "tk", "op", "lean", "reading", "fractional-powerset-dials"),
    ("Fractional power set needs the absence of choice", "t", "op", "partial", "new", "fractional-powerset"),
    ("Prime products from fermions to bosons (primon gas)", "kt", "op", "checked", "reading", "dial-survey"),
    ("Velocity addition is ∘₁; rapidities add", "k", "op", "lean", "reading", "relativity"),
    ("Schwarzschild clock and ruler deficits are circle inverses", "k", "op", "checked", "reading", "general-relativity"),
    ("Holographic entropy cone = corner sum at k = −1; fuzzy cuts keep monogamy", "kt", "op", "checked", "new", "holographic-entropy-cube"),
    ("Tunnelling: averaged odds compose by ∘ₖ with k = −2K", "kd", "op", "checked", "reading", "tunnelling"),
    ("Rulecraft: same-target rules fuse at k = −1; sort by target", "kq", "op", "checked", "reading", "rulecraft-history"),
    ("Pap's generated semiring, one unified proof", "k", "op", "lean", "known", "semiring-audit"),
    ("Lee–Yang zeros of the dial gas on rings", "k", "op", "lean", "reading", "lee-yang"),
    ("Decrement paths: Shapley split, and the gap when causes depend", "t", "op", "checked", "new", "decrement-paths"),
    ("Fractal half boxes: dimension as a fractional circle power", "tr", "op", "lean", "reading", "fractal-dials"),
    # --- uses the plane law or wider composition ---
    ("Map of one-parameter composition laws (the k, d plane)", "kd", "plane", "lean", "known", "landscape"),
    ("Hirzebruch's χ_y genera on one straight line", "kd", "plane", "checked", "known", "landscape"),
    ("Light sets the dial: d = −1/c²", "d", "plane", "lean", "reading", "electromagnetism"),
    ("Frame dragging from boosting Newton's pull", "d", "plane", "checked", "reading", "gravity-from-theory"),
    ("River rebuilt: Schwarzschild exact, Kerr by one turn", "d", "plane", "lean", "reading", "gravity-from-theory"),
    ("Swapping two maps moves the result by an oriented area", "q", "plane", "lean", "known", "order-dial"),
    # --- dial picture only ---
    ("Elliptic modulus on the curvature dial (pendulum)", "d", "pic", "lean", "reading", "elliptic-modulus"),
    ("Force reach: carrier mass² = Δ/4; flux tube = dimension", "d", "pic", "checked", "reading", "force-reach"),
    ("Hidden circle: α = 4G|Δ|, a short-range force at 1/3 gravity", "d", "pic", "partial", "known", "force-reach/new-force"),
    ("Rulecraft: Δ counts anchors; loop product b·c; clash count", "dq", "pic", "checked", "reading", "rulecraft-history"),
    ("Riemann round 3: barrier chain = Euler product; no road", "d", "pic", "checked", "known", "riemann/round3"),
    ("Lee–Yang on the ℓp ball: real zeros iff n ≤ 1 + p²", "p", "pic", "partial", "new", "shape-dial-survey"),
    ("Bell's bound slides along p: 4 → 2√2 → 4", "p", "pic", "checked", "reading", "shape-dial-survey"),
    ("Cubic form: round to triangle; conics; affine spheres", "Cd", "pic", "checked", "known", "affine-sphere"),
]

if __name__ == "__main__":
    from collections import Counter
    print(len(RESULTS), Counter(r[2] for r in RESULTS), Counter(r[3] for r in RESULTS), Counter(r[4] for r in RESULTS))
