# The k dial on the cusp: the chain becomes q-integers  (DRAFT)

**Plain result.** Put the homothety dial on the chain's step, so each move adds with our law
a ∘ₖ b = a + b + kab instead of a + b. Because 1 + k(a ∘ₖ b) = (1 + ka)(1 + kb), the chain's sums
become exactly the q-integers [m]_q = (q^m − 1)/(q − 1), with q = 1 + k and m running over the
same fractions as before (checked on the tree, to 1 part in 10¹²). So **the k dial on the chain is
the q dial**, with q = 1 + k. The chain sum turns into

  F_k(β) = Σ φ(m)·[m]_q^(−β),

a q-version of Knauf's ζ(β−1)/ζ(β). For any k > 0 it converges for every Re β > 0, with no cusp
and no wobble: the k dial does make the slow tail fast.

**Dialling back to k = 0** works, but only in step with k. Written out exactly, F_k is the plain
totient series with a smooth Bose-shaped weight (x/(eˣ − 1))^β on each term, x = m·log(1 + k)
(picture A2). The chain itself cuts the fractions off raggedly by depth, which is where the cusp
came from. A smooth cut has a clean expansion, but its next terms carry ζ's own zeros at the same
height. So the miss shrinks only like k (picture A1): about 1% at k = 0.0002. Adding the known next
shapes (k and k log k) doesn't help. That is less accurate than fitting the cusp out of chains up to
length 33 (cusp-size, about 0.1%).

## Grades

| | worth following? | could reach a proof? |
|---|---|---|
| k dial on the cusp | **Modestly.** It's a clean framework fact: on the chain, the k dial and the q dial are one move (q = 1 + k), and the k dial turns the cusp's sharp cut into a smooth Bose-weighted one. As a way to the zeros it is slower than the cusp fit. | No. Dialling back to k = 0 is the classical smoothed Dirichlet series, whose error terms are governed by ζ's zeros themselves, so it can only restate RH. |

## Pictures

- `figures/figA_k_dial_rewind.png`: A1, distance from ζ(β−1)/ζ(β) as k is dialled down, at four heights, against the cusp fit; A2, the weight each denominator gets: ragged for the chain of length 33, smooth for the k dial.

## Known vs ours

Known: q-integers and q-analogues of ζ that tend to ζ as q → 1 (for example Kaneko–Kurokawa–Wakayama);
the circle operation's (1 + ka)(1 + kb) form (Perlis 1942, Jacobson 1945); smoothed Dirichlet sums and
their explicit-formula error terms (classical).
Ours (not found in a quick search): the k dial on Knauf's chain giving exactly Σ φ(m)[m]_q^(−β) with
q = 1 + k; the comparison of the smooth and ragged cuts on the line Re β = 3/2.

## Files

`identity.py` → `identity.txt` (tree sums = q-integers), `rewind.py` → `rewind.json` (dial back to 0;
the known pole term H(2−β)·L^(β−2)/ζ(2) is tiny on this line), `extrap.py` → `extrap.json` (next
shapes), `figs.py`, `checks.py` (9 checks). Panel A2 reads `/tmp/knauf_levels/levels_33.npz` from
`../length-smoothing/stream.py`.
