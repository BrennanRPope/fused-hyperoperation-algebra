# DRAFT · From an orbit's word to its start

*Three equal bodies, zero angular momentum, all starting in a line. Orange in the pictures is known math, blue is ours, red does not work, grey is neither.*

## The answer first

**There is no formula, published or found here, that turns an orbit's word into its starting push.** Simple word counts miss the start by about 15 times the gap between neighbouring orbits.

**What does work is a recipe, and a sharper period law.** The word's letters, read like digits, plus a new rule that gives the period from the word's runs, rebuild 32 of 68 held-out orbits from the word alone. The period rule is the lead worth following: it is 12 times tighter than the known one.

**Brennan's two extra ideas do not add a lead.** In hyperspherical coordinates the period, word length and stability come out no simpler, and gluing three whole hyperspheres loses the words; the one gluing that keeps them is the shape sphere, already known.

## The question

Every periodic orbit in this family is named by a word: each time the three bodies pass through a straight line, write the letter of the body in the middle (a or b), in capitals when it passes the other way. The figure-eight is aBAb.

Each orbit also has a starting point: the bodies start in a line at (-1,0), (1,0), (0,0), the two outer bodies get the same push (v₁, v₂) and the middle one gets the opposite push, twice as hard. So the question is whether some function takes the word to (v₁, v₂) and to the period.

We used the 695 orbits that Li and Liao published in 2017, with their starts, periods and words.

## 1. The period: known, and a sharper version

![Period from the word: the known straight line, the run prices, and how much the misses shrink](figures/fig1-period-from-the-word.png)

**Known (Dmitrašinović and Šuvakov, 2015).** Put every orbit at the same energy, and its period grows in step with the word's length: about 2.46 time units per letter. Tested on orbits left out of the fit, that rule misses by 3.3% on a typical orbit, and by up to 31%.

**Ours: each run has a price.** Split the word (read round in a circle) into runs, a run being letters that all turn the same way before a reversal. aBAb has two runs of 2. Then the period is the sum of a fixed price per run:

| run length | price | price per letter | runs seen in 695 words |
|---|---|---|---|
| 1 | 1.37 | 1.37 | 24 |
| 2 | 4.62 | 2.31 | 8,052 |
| 3 | 7.43 | 2.48 | 10,196 |
| 4 | 9.97 | 2.49 | 7,830 |
| 5 | 12.38 | 2.48 | 56 |

Five numbers, fitted with each test orbit left out, predict the period to **0.28%** on a typical orbit and never worse than 1.5%. That is 12 times tighter than counting letters. On 54 orbits of our own, found by a different search, it lands within 1% for 52 of them.

![Each run gets one colour on the shape sphere; the run prices add up to the period](figures/fig2-runs-on-the-shape-sphere.png)

The picture of why: on the shape sphere (each point is a triangle shape, the three crosses are collisions), a run of 2 is a figure-eight-like swing between two collision points, and a run of 4 is a butterfly-like loop round one. The time each kind of piece takes is nearly fixed, so the times add. Runs of 2 are cheaper per letter than longer runs, which is why counting letters alone misses.

We did not find this rule in what we could read. The straight-line law is known; the run prices are, as far as we found, new.

## 2. A formula for the start from word counts: does not work

![Where the orbits start, and how far word-count formulas miss](figures/fig3-counting-letters-misses.png)

Starts crowd onto a few curves in the quarter disc, with a typical gap of 0.004 between an orbit and its nearest neighbour. We tried straight-line fits on letter counts, run counts and symmetry class, and nearest look-alike words. On held-out orbits they miss by about 0.06, fifteen times the gap. **No word-count formula places an orbit.** Nothing in the literature claims one either: the 2015 paper predicts periods only.

## 3. Letters as digits: the word does pin down a region

![Butterfly I: starts kept while their letters keep agreeing with the word](figures/fig4-letters-as-digits.png)

The right way to read a word is like the digits of a number. Run every start in the quarter disc (502,041 of them on an 800 × 800 grid), write down its first 32 letters, and keep the starts that agree with the word letter by letter. Each letter keeps about 85% of the survivors.

![Each letter narrows the start like a digit](figures/fig5-each-letter-narrows.png)

But the survivors shrink to a thin strip, not a point. Every start on that strip shadows the orbit for a while. The word alone cannot say where along the strip the orbit is; something else has to close the loop.

## 4. The recipe: letters, then the run-law period, then a polish

![Word to start for orbit I.A9: the letter strip, the closest return at the predicted period, and the polished orbit](figures/fig6-word-to-start.png)

That something is the period. The recipe uses only the word:

1. Keep the grid starts whose letters match the word.
2. Predict the period from the run prices (refitted with this orbit left out) and, along the strip, find the start that comes back closest to itself at that period.
3. Polish that guess with a standard shooting method until the orbit closes.

We then compared the result with the catalog start, which the recipe never saw.

![From the word alone, 68 held-out orbits of every length](figures/fig7-word-alone-results.png)

Out of 68 held-out orbits, from 4 to 136 letters:

- **32 land on the catalog orbit** (start within a millionth): 27 on the first pass, 5 more on a second, finer pass.
- **6 land on another closed orbit with exactly the same word** and the same period to 5 or 6 digits, but a different start. For 3 of them the path on the shape sphere is clearly different from the catalog orbit's; for the other 3 we could not rule out that it is the catalog orbit with the bodies relabelled.
- **5 close up into an orbit with a different word**, and **25 find nothing** that closes.

Short words work best: 8 of 11 up to 40 letters, 13 of 22 from 41 to 80, and 11 of 35 above 80. Long words fail because the 800 × 800 letter map is too coarse to hold their thin strips; the best guesses for the misses are typically between a ten-thousandth and a few thousandths off, close but outside the polish's reach.

![The figure-eight and a second closed orbit with the same word](figures/fig8-same-word-two-orbits.png)

**A word does not always name one orbit.** For the figure-eight's word aBAb, the polish also found a second closed orbit (in an earlier run): same word, same scale-free period to 6 digits (9.237681 against 9.237683), a different start (0.3394, 0.5362 against 0.3471, 0.5327), and each body on its own curve instead of all three on one. Polishing from points between the two lands on one or the other, so they are two separate orbits, not a family. We have not checked whether this one is already known.

## 5. Invented words: does not work

We also made 30 new words, by splicing halves of catalog words in the same symmetric pattern and keeping only those not in Li and Liao's list, and asked the recipe for orbits with exactly those words. **None gave a new orbit.** Two "hits" were the figure-eight's word written twice and three times over, which is just the figure-eight going round again. For most of the others, no start anywhere on the map follows the invented word, repeated, for more than 5 to 26 letters. Words are not free to choose: splicing breaks whatever rule real orbits obey. That rule is a good question in itself.

## 6. Hyperspherical coordinates: the words' home, but nothing gets simpler

Brennan asked to try the problem on a hypersphere, in polar coordinates there. Once the centre of mass is taken out, three bodies in a plane are one point in four dimensions. Polar coordinates there split it into one size, the hyperradius R, and a direction on a 3-sphere. Forgetting the overall turn leaves Montgomery's shape sphere, where the three collisions are the punctures and the words come from (known). So this is the natural home for the words, and we rewrote our 66 polished orbits and Li and Liao's 695 in it.

![Six of our polished orbits: direction on the shape sphere and size over one period](figures/fig9-hyperspherical-orbits.png)

Then we asked whether the period, the word and the stability come out simpler in these coordinates. Three new clocks were tried against ordinary time, with the same held-out test for each:

![Which clock the word predicts best, the run law on our orbits, and stability per clock](figures/fig10-hyperspherical-scores.png)

- **Hyper-angle time** (time divided by R², the clock the direction on the sphere naturally keeps): the word predicts it to 4.7% at best, against 0.3% for ordinary time.
- **Length of the path on the shape sphere**: 3.3%.
- **Breaths of the size R** (how many times R swells per period): they do not match letters or runs, and the period per breath misses by 28%.
- **Stability**: how fast nearby orbits peel away, per unit of each clock, spreads by 80% to 115% of its average whatever the clock. No clock turns it into a constant rate.

One useful thing did come out. On our 66 polished orbits the run law works on every line-up the orbit makes, not on the shortened word: 62 of 66 within 1% (the other 4 within 2.3%), against 52 of 66 when letters that cancel are struck out first. In the catalog the two are the same, which is why this did not show earlier.

**Verdict:** the hypersphere is where the words live, which was already known. In its own coordinates the period, word length and stability come out no simpler, so it does not open a new lead. The sharp law stays in ordinary time.

## 7. Gluing three hyperspheres: only one way works, and it is the shape sphere

Brennan also asked whether gluing three hyperspheres together might help. The words are loops round the collision points, written in two letters, so a glued shape can only help if it has two independent loops. That decides it:

![Ways of gluing three spheres, and which ones can hold the words](figures/fig11-gluing-spheres.png)

1. One sphere per pair of bodies, each touching the next in a ring: one loop only. One letter cannot write these words.
2. Three spheres in a chain, or joined along discs: no loops at all.
3. One cap per collision (bodies 1 and 2 close, 2 and 3, 1 and 3), glued along their rims, with the collision points themselves removed: two independent loops, exactly the two-letter words. This is the shape sphere with its three punctures (known, Montgomery).
4. In four dimensions, the same thing is the 3-sphere with three linked circles removed, one circle per collision. Its loops are the words times one extra overall turn (known).

**Verdict:** gluing whole spheres loses the words. The gluing that keeps them gives back the shape sphere, so it adds nothing new, and it says nothing about the run prices. It does not help.

## Should this be followed up?

**The run law, yes.** It is simple, it is 12 times tighter than the known law, it holds on orbits from a different search, and it has a picture behind it: an orbit is built from a few kinds of pieces whose times add. The natural next steps are to test it on the larger 2024 catalogs and with unequal masses or nonzero angular momentum, and to ask why a run of 2 costs 4.62 and a run of 4 costs 9.97.

**A closed formula for the start, no.** The letters behave like digits of a chaotic map, so the start is a limit, not an expression. The recipe above is the constructive answer: word in, orbit out, for about half of the real words we tried (32 of 68, and most of the short ones). Making it reliable for long words needs a finer letter map, not a formula.

## Sorting board

| | what | verdict |
|---|---|---|
| Known | period ≈ 2.46 × letters (2015), 3.3% typical miss | holds |
| Known | letters as symbolic dynamics (Tanikawa and Mikkola, and others) | holds, used here |
| Ours | period = sum of five run prices, 0.28% typical miss | holds on all 695 and on 54 of ours |
| Ours | word → start recipe: letters, run-law period, polish | 32 of 68 held-out orbits land on the catalog orbit |
| Ours | a second orbit with the figure-eight's word | found, not checked against newer catalogs |
| Does not work | word counts → start formula | misses by 15× the gap |
| Does not work | spliced, invented words → new orbits | 0 of 30 |
| Known | hypersphere and shape sphere as the home of the words (Montgomery) | holds |
| Does not work | hyper-angle time, shape-sphere length or breaths of R as a simpler period | 3.3% to 28% misses, against 0.3% |
| Does not work | a clock that makes stability a constant rate | spread 80% to 115% |
| Ours | run law counts every line-up, not the shortened word | 62 of 66 of our orbits within 1% |
| Does not work | gluing three whole spheres | one loop or none, the words need two |

## What we could not check

The larger 2024 catalog of Hristov, Hristova and coauthors and the authors' own data pages were out of reach from here, so the run law and the twin orbit were not tested against them.

## Sources

- Dmitrašinović and Šuvakov, "Topological dependence of Kepler's third law for collisionless periodic three-body orbits", Phys. Lett. A (2015), arXiv 1507.08096.
- Li and Liao, "More than six hundred new families of Newtonian periodic planar collisionless three-body orbits", Sci. China (2017), arXiv 1705.00527; data from the sjtu-liao/three-body repository on GitHub.
- Hristov, Hristova and others, 2024 catalog of periodic three-body orbits, arXiv 2404.16526 (not read; data blocked here).
- Šuvakov and Dmitrašinović, "Three classes of Newtonian three-body planar periodic orbits", PRL (2013).
- Tanikawa and Mikkola, symbolic dynamics of the free-fall and planar three-body problems.
- Montgomery, "The N-body problem, the braid group, and action-minimizing periodic solutions", Nonlinearity (1998), and "Infinitely many syzygies" (2002): the shape sphere, its punctures and the words.

*Code and data: math-theory/three-body-word-formula (hyper.py, hyper_fit.py and hyper_figures.py for the hypersphere and gluing sections, period_runs.py for the run law, predict2.py and retry.py for the recipe, predict_new.py for invented words, figures.py for the pictures).*
