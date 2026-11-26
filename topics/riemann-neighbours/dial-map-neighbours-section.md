<!-- To insert in topics/dial-map/README.md as a new section "8. Neighbours", after section 7 and before "Assumptions". DRAFT, Nov 26 2026. -->

## 8. Neighbours: finding associated work on the map (DRAFT, Nov 26)

Brennan: "this is how i feel like the map will be most useful is to be able to find associated work". Every problem placed on the map gets a short **Neighbours** table: the other problems on the same layer with the same shape of statement. Proved neighbours come first, because those are the ones we can learn a method from; then false ones, because they rule methods out; then open ones.

| Column | What goes in it |
|---|---|
| Neighbour | the problem's name |
| Where | same layer, one layer up, or one layer down |
| Status | proved (by whom, year), open, or false |
| Key ingredient | the one idea that makes the proof work, or why it fails |
| Can our dials say it? | yes / partly / no, with the dial |
| Mathlib | what it could build on, and the main gap |

### First example: the Riemann hypothesis

RH sits on Layer 1 (k = −1 folded over the prime cube gives 1/ζ) and Layer 2 (each local factor 1 + k·u + d·u² should be on the elliptic side). Full write-up: [riemann-neighbours](../riemann-neighbours/README.md).

| Neighbour | Where | Status | Key ingredient | Our dials? | Mathlib |
|---|---|---|---|---|---|
| Ihara zeta, Ramanujan graphs | Layer 1 → 2 | proved for LPS/Margulis graphs (1988); every degree, bipartite (Marcus–Spielman–Srivastava 2015) | fold = det(I − uB), a polynomial; Ramanujan bound | yes: Ramanujan is \|k̂\| ≤ 2 | walk counts, Hermitian adjacency; no Ihara zeta |
| Curves over F_p (Hasse, Weil) | Layer 2 | proved (Hasse 1933, Weil 1948, Deligne 1974) | Frobenius; degree form positive definite | yes: the degree form is our norm a² + kab + db² | `WeierstrassCurve.localPolynomial` is our dial point; no Hasse bound |
| Selberg zeta | Layer 2 | proved (Selberg 1956); λ₁ ≥ 1/4 open for congruence surfaces | self-adjoint Laplacian | yes: factor s² − s + λ is (−1, λ), elliptic iff λ ≥ 1/4 | nothing |
| Lee–Yang | Layer 1 with edges | proved (Lee–Yang 1952, Asano 1970) | ferromagnetic coupling, contraction | yes: two spins are ∘ₖ with k = 1/a; more spins are k on each edge (our dial gas) | nothing |
| Davenport–Heilbronn | Layer 2 only | **false**: zero at 0.8085 + 85.6993 i | functional equation without a fold | shows Layer 1 is needed | Dirichlet L-functions are there |
| Dirichlet, Dedekind, Hecke, Selberg class | Layers 1 and 2 | open (GRH) | none yet | same as ζ | Euler products, functional equations, non-vanishing on Re s = 1 |
| Bost–Connes | Layer 3 | partition function ζ(β), phase transition at β = 1 (Bost–Connes 1995) | sees the pole, not the zeros | yes: the q dial | nothing |

**What the neighbours say together.** Every proved row has an operator whose determinant is the prime fold and a positivity that puts its factors on the elliptic side. Our dials say the positivity. The operator for ζ is unknown (Hilbert–Pólya), which is why RH stays open on the map.
