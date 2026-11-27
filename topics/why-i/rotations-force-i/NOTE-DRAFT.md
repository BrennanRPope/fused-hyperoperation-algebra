**DRAFT. Not sent. No recipient chosen.**

Subject: Does "only energy differences matter" force one shared i? Is this already known?

Hi [name],

I have a short argument that quantum theory over the reals needs one quarter turn J (J² = −1) shared by every Hamiltonian. I'd like to know whether it's already in print, and whether any of the assumptions is doing more work than it looks like.

**Setting.** States are real unit vectors in ℝⁿ. Observables are symmetric matrices, and any observable can be the energy of some arrangement. Each Hamiltonian H generates dynamics ψ' = A(H)ψ, with:

- A1 (rotations): A(H) is antisymmetric.
- A2 (conservation): [A(H), H] = 0.
- A3 (Planck–Einstein): A(H)² = −H²/ħ².
- A4 (only energy differences matter): A(H + c·1) = A(H) + c·A(1) for every real c.

**Claim.** A(H) = JH/ħ for one fixed J, the same for every H. J is antisymmetric, J² = −1, and J commutes with every observable.

**Proof.**

1. A3 at H = 1 gives J := ħA(1) with J² = −1. By A1, J is antisymmetric.
2. Let B = ħA(H). A2 for H + c gives [B, H] + c[J, H] = 0 for all c, so [J, H] = 0.
3. A3 for H + c gives (B + cJ)² = −(H + c)². Matching powers of c gives B² = −H² and BJ + JB = −2H.
4. Write B = JH + E. Then E is antisymmetric, anticommutes with J and commutes with H. Squaring gives B² = −H² + E², so E² = 0.
5. For antisymmetric E, tr(EᵀE) = −tr(E²) = 0, so E = 0.

**Consequences.**

- Real QM, where every symmetric matrix is an observable, can't satisfy A4: only scalars commute with all symmetric matrices.
- With J shared, the observables lie inside the complex Hermitian matrices. Adding "every symmetric operator commuting with J is an observable" rules out the quaternionic case.
- For two systems, I add a sixth assumption: a shift of the total energy has no owner, so J⊗1 = 1⊗J on joint states. That cuts the real tensor product down to the complex (balanced) one, where J⊗J = −1 and local tomography holds.

**What I know is already out there.**

- Each single dynamics carrying its own positive-frequency J: Segal; Ashtekar–Magnon; Aste, arXiv:1905.12894; Vaughan, arXiv:2609.23041.
- A universal J: postulated by Stueckelberg (1960); derived from Poincaré symmetry by Moretti–Oppio, arXiv:1611.09029.
- The hidden J and balanced tensor product in recent real formulations: Bang–Cho–Baek, arXiv:2607.05865.

**My questions.**

- Has the energy-zero step (A4 plus Planck–Einstein forcing one J for all Hamiltonians) been published, perhaps in the older real-Hilbert-space literature?
- Is A4, read as "a shift adds a fixed turn and nothing else", a fair statement of "only energy differences matter", or does it smuggle in the answer?

I've checked each step symbolically and numerically (sympy, 15 checks), and I can send the full write-up with pictures.

Thanks,
Brennan
