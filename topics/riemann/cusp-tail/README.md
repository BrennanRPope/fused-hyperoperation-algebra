# Cusp tail: the chain's wobble is mostly its cusp, and can be fitted out  (DRAFT)

**Plain result.** Finite Knauf chains have zeros next to ζ's zeros (on the line Re β = 3/2), but
each misses by about ±0.01, swinging from side to side as the chain grows (length-smoothing). This
piece finds where most of that swing comes from: the chain's **cusp**, the edge states near the
fractions 1/n that the Farey map's neutral (parabolic) point produces. At each zero ρ the chain's
error behaves like A·n^(1−ρ), the same shape as the tail of a partial ζ sum. Multiply the error by
n^(ρ−1) and it stops spinning and sits in one place (picture A).

Knowing the shape, the cusp tail can be fitted out: model Z_n(β) as F(β) + A(β)·n^(1−β) over a
window of lengths and solve for F. The zeros of F miss the line by about 0.0005 from length 29 to 33,
2 to 3 times better than blind averaging over the same window on lengths 29 to 32 (at 33 the two
meet). With only two lengths the fit already beats a single chain at every length from 16 to 33.
Wrong shapes in the same fit do 1.4 to 20 times worse (picture C).

**Correction to length-smoothing.** That piece reported a "swing law" of about 0.05 × height
radians per length. The turn is really about −γ·log((n+1)/n) ≈ −γ/n per length, slowing as chains
grow; 0.05 × height was its average over lengths 16 to 30 (where 1/n averages about 0.045). The
law's upshot (higher zeros swing faster) stands; the "steady angle" does not.

## Grades

| | worth following? | could reach a proof? |
|---|---|---|
| Cusp tail | **Yes.** It names the source of the wobble with a checkable shape, and a fit built on that shape beats blind averaging at the longest lengths. Next: work out A(β) from the cusp itself (fractions with one large continued-fraction digit) instead of fitting it, so that only the remainder is left. | Not as it stands. It splits the wobble into an explicit cusp part and a smaller remainder, which is the kind of split a proof would need, but the remainder still has to be bounded at every height, and that is as hard as RH. |

Limits: lengths up to 33, first 5 to 8 zeros, numerical. The fitted miss drops in a step near
length 29 and then sits near 0.0005; whether it keeps falling is not settled at these lengths. For
zeros above height 32 the cusp explains less of the swing (picture A3).

## Pictures

- `figures/figA_swing_is_cusp.png`: A1, the chain's error at ζ's 4th zero spins round 0 as the
  chain grows; A2, multiplied by the cusp's spin n^(ρ−1) it sits still; A3, how still it sits, for
  6 zeros.
- `figures/figB_cusp_removed.png`: sideways miss against length, for one chain, the averaged chain,
  and the chain with its cusp tail fitted out.
- `figures/figC_null_test.png`: the same fit with wrong shapes.

## Known vs ours

Known: Knauf's chain and its limit ζ(β−1)/ζ(β) for β > 2; its states are Stern–Brocot fractions,
so length n holds every fraction with denominator up to n+1 exactly φ(q) times (checks.py confirms
this); the chain's phase transition at β = 2 comes from the parabolic point of the Farey map
(Kleban–Özlük, Prellberg, Fiala–Kleban–Özlük, who also study finite-size asymptotics at the
critical point); tails of partial Dirichlet sums go like n^(1−s); fitting known tail shapes is
ordinary Richardson-style extrapolation.

Ours (not found in a quick literature search): that the zero wobble of finite chains on
Re β = 3/2 is mostly the cusp's n^(1−ρ) spin; fitting that tail out to locate the zeros; the
null tests; the correction to the swing law.

## Files

- `probe.py`: per-length turn against −γ/n.
- `demod.py` → `demod.json`: raw and cusp-demodulated errors, and how still each sits.
- `correct.py` → `correct.json`: zeros with the cusp fitted out over 2 to 8 lengths.
- `combo.py` → `combo.json`: cusp fit over half the chain.
- `null.py` → `null.json`: wrong shapes.
- `figs.py`: pictures. `checks.py`: 13 checks.
- Lengths above 20 use the level counts streamed by `../length-smoothing/stream.py` (in
  /tmp/knauf_levels; rerun it to 33 before rerunning these scripts).
