**DRAFT. Not sent.** You send it yourself, from your own mail, with the paper's PDF attached (compile it in Overleaf first).

- **To:** [name and address]
- **Where the address came from:** choose a recipient; tailored versions are kept in project files
- **Before sending:** replace [one line about you] with a sentence on who you are, and check the attachment is the current PDF.

---

Subject: Does "only energy differences matter" force a single i in real quantum theory?

Hi [name],

I have a short argument that quantum theory over the reals needs one complex structure shared by every Hamiltonian, and I'd like to know whether it is already in print before I post it. [one line about you]

The setting is a real Hilbert space where each Hamiltonian H generates dynamics dψ/dt = A(H)ψ. Assume A(H) is antisymmetric, the Planck-Einstein relation A(H)² = −H²/ħ², and that shifting the energy zero only adds a fixed rotation, A(H + c·1) = A(H) + c·A(1). Then J = ħA(1) satisfies J² = −1, and A(H) = JH/ħ with [J, H] = 0 for every H. The proof is five lines: expand A(H + c·1)² in powers of c and match terms, then a trace argument removes what is left. Energy conservation comes out rather than going in.

Each input can be weakened. The shift rule follows from asking that two non-interacting systems evolve by a rule fixed by their total Hamiltonian, since a constant in the total energy has no owner. Antisymmetry follows from stability, given observables that fix states up to phase. Planck-Einstein follows from asking only that an energy eigenstate's rotation rate be some fixed function of its energy, and ħ comes out as a constant. For two systems the same argument gives J⊗1 = 1⊗J, which cuts the real tensor product down to the complex one. The algebraic steps are checked in Lean 4, and every claim has a numerical check.

The closest work I found is Alfsen and Shultz's dynamical correspondence and Barnum, Müller and Ududec's energy observability postulate, which both assume a linear map from all observables to generators. Here linearity is assumed only for shifts of the energy zero. Have you seen this argument before? The draft is attached.

Thanks,
Brennan
