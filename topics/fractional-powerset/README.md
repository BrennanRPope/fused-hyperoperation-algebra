# Fractional power set: a dial between X and P(X) (Nov 25, 2026)

> **Names:** some words or symbols on this page (mult-add, rung, wall, twist) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Task: "Find a fractional power set for any in-between rung" (Brennan's ask: half as a first test, really a family P^t for all t, P^s∘P^t = P^{s+t}, P^1 = P, by size). Also maps grades of the axiom of choice against candidate levers (rung, q, twist, flow position, symmetry group).

Page: published as an artifact (see PAGES-INDEX.md, "fractional power set"); source index.html here.

## Files
- checks.py / checks_output.txt: 18 of 18 pass. Pure Python (fractions), no libraries.
- lean/FractionalPowerset.lean: 11 theorems, core Lean 4, no imports, no sorry. UNCOMPILED here (no toolchain); for the library thread to compile and merge as Topics/FractionalPowerset.
- index.html: the page (inline SVG plates, same style as infinities/).

## Findings
1. Finite wall: no non-decreasing f: N->N with f(f(n)) = 2^n (fails at n=0,1). From n=2 up a strictly increasing one exists (3,4,8,9,10,11,16,32,64,128,129,...).
2. With choice, a monotone half power set needs k < F(k) < 2^k, so GCH kills it; in Merimovich's gap-2 world (2^l = l^{++} for all l, from a large cardinal) successor is a half power set. Square root of "jump d rungs" exists iff d even (Lean: half_shift_even).
3. With choice, no dyadic dial: halvings P^{1/2^k}(k) strictly decrease forever; impossible for well-ordered sizes (Lean: no_dyadic_dial on N).
4. Without choice: obstruction gone (Jech–Sochor; Karagila 2014 with DC_k). No result found either way: open.
5. Natural constructions blocked: half iterate of e^x - 1 = x + x^2/4 + x^3/48 + 0 + x^5/3840 - 7x^6/92160 (Baker 1958, OEIS A052104); n!*coef = 1/2 at n=2, negative at n=6.
6. Mostowski model: subsets of atoms with least support of size n number 2*3^n; power set turns orbit growth o(n) into ~2^{o(n)} per support (pairs: 2^{2n^2}). Heuristic: a half step needs a symmetry group with orbit growth between polynomial and exponential (e.g. S_inf wr S_inf, partition numbers ~ e^{c sqrt n}); Macpherson/Simon gaps constrain.
7. Choice grades by "which rung is flat on infinite sizes": successor flat (m+1=m) <=> every infinite set D-infinite; addition flat (2m=m) strictly between; multiplication flat (m^2=m) <=> AC (Tarski 1924), the circle operation ∘₁ (called mult-add before 1.0) flat <=> AC; exponent never flat (Cantor). Landscape is a partial order (BPI vs DC incomparable).
8. Levers: rung (real, 3 steps); symmetry group (real, landscape); DC height (real ladder, doesn't touch the dial); k (two-way switch: k=0 add rung, k>=1 mult rung); twist (no work, size arithmetic commutes); flow t (lives on growth rates).

## Verdicts
- Fractional power set: lead worth following (open in ZF; next step: permutation model from a group with e^{c sqrt n} orbit growth).
- Choice as position: real picture, modest lead (check whether aleph_0*m = m is its own grade, maybe in Howard–Rubin).

## Sources
Merimovich arXiv math/0005179; Karagila arXiv 1212.4396; Halbeisen–Shelah arXiv math/9308220 and math/0010268; OEIS A052104; Simon arXiv 1810.06531; Gitik arXiv math/0212405.
