SPEC=dict(
title="Ecology and the circle operation",
h1="A predator's time budget puts Beddington–DeAngelis and Crowley–Martin on one dial, and Tsallis diversity has its ceiling at the wall",
lede="Model a predator that searches, handles prey, and gets drawn into fights with other predators. If a handling predator is drawn in at k times the rate of a searching one, the feeding rate is αN/(1 + x ∘ₖ y), with x the handling load and y the interference load. k = 0 is the Beddington–DeAngelis functional response, k = 1 is Crowley–Martin, and k above 1 means a predator holding prey attracts trouble. Separately, Tsallis diversity is the q-logarithm of Hill's effective number, so independent communities combine by the dial with k = 1 − q. For q above 1 the wall 1/(q − 1) is exactly the ceiling that diversity approaches as species are added.",
stats=[("11 / 11","numeric checks"),("10","Lean theorems, no sorry"),("1","modest lead"),],
hero_fig="figures/f1_functional_response.png",
hero_cap="Left: feeding rate per predator as predator density rises, with prey plentiful, for four values of the dial. Right: the 95% spread of fitted k across 200 simulated feeding experiments with true k = 0.5. Only the large, low-noise design separates the two classic forms.",
known="Beddington (1975), DeAngelis et al. (1975) and Crowley and Martin (1989) are the two classic predator-interference responses, and Skalski and Gilliam (2001) found that neither wins across datasets. Tsallis pseudo-additivity, S(A·B) = S(A) + S(B) + (1 − q)S(A)S(B), is well known in physics, and Jost (2006) brought the Hill-number view to ecology. The 'additive' and 'multiplicative' nulls for multiple stressors are standard (Folt et al. 1999; Crain et al. 2008). In a short search I did not find the one-parameter family joining Beddington–DeAngelis and Crowley–Martin, or its four-state derivation. It is simple enough that it may exist under another name.",
lead="Modest, noted rather than followed. The dial turns 'which functional response?' into 'what is k?', a single number with a behavioural meaning that a standard feeding experiment can estimate. A fitted k above 1 would point to kleptoparasitism, a handling predator drawing others in. But the simulations show the cost: separating k = 0 from k = 1 takes about 120 trials at 10% noise, and at a more typical 25% noise even 120 trials do not. That may be why published datasets disagree. Taking this further needs those datasets, which I did not have here.",
kicker="Biology series · 10 of 11",
sections=[
("fr","Functional response","1. Predator interference is a dial",
"""<p>A predator is in one of four states: searching, handling a prey item (entered at rate αN, left at rate 1/h), interfering after meeting another predator while searching (rate γP, left at 1/τ), or interfering while handling, entered at k·γP. The steady-state share of time spent searching is one over</p>
<div class="eq">1 + x + y + k·xy = 1 + (x ∘<sub>k</sub> y), &nbsp; x = αNh, &nbsp; y = γPτ</div>
<p>so the feeding rate is αN/(1 + x ∘<sub>k</sub> y). The steady state of the four-state chain matches the formula to 1 × 10⁻¹⁴, and a stochastic simulation of one predator over 200,000 time units ate within 0.1% of it.</p>
<ul><li><b>k = 0:</b> handlers never interfere. This is Beddington–DeAngelis, 1 + x + y.</li>
<li><b>k = 1:</b> handlers interfere as often as searchers. This is Crowley–Martin, (1 + x)(1 + y).</li>
<li><b>k above 1:</b> holding prey draws other predators in, as with food theft.</li></ul>
<p>In a simulated 30-trial experiment crossing six prey densities with five predator densities (true k = 0.5, 4% noise), the dial fit recovered k = 0.56 and beat both classic forms on AIC. With more realistic noise the picture changes: the spread of fitted k is 0.23–0.98 for 30 trials at 10% noise and 0.02–3.1 at 25%; 120 trials narrow it to 0.34–0.72 and 0.15–1.6.</p>"""),
("div","Diversity","2. Diversity: the wall is the ceiling",
"""<p>Hill's effective number of species D<sub>q</sub> multiplies when two independent communities are combined, and Tsallis (HCDT) entropy is its q-logarithm, S<sub>q</sub> = ln<sub>q</sub>(D<sub>q</sub>). So the entropies combine by</p>
<div class="eq">S(A·B) = S(A) ∘<sub>1−q</sub> S(B)</div>
<p>q = 1 is dial 0, Shannon entropy adding. For q above 1 the wall −1/k is 1/(q − 1), and that is exactly the ceiling of S<sub>q</sub>: at q = 2, communities of 10, 100 and 10,000 equally common species give 0.9, 0.99 and 0.9999, never reaching 1. A community at the ceiling absorbs any other, as a wall should. This is one of the few places in the series where the wall is a natural bound rather than an unphysical point.</p>"""),
("other","Stressors","3. Multiple stressors and colonisation",
"""<p>Ecologists judge two stressors synergistic or antagonistic against a null. Adding proportional losses is dial 0, and independent action, the usual multiplicative null, is the 'or' law, dial −1, the same choice drug combinations face on the cell biology page. Island colonisation from n sources is wpow(−1, c, n).</p>"""),
],
table_h2="Where the circle operation (called mult-add before 1.0) shows up in ecology",
table=[
["Predator interference","Feeding rate αN/(1 + x ∘ₖ y); k = 0 BD, 1 CM, above 1 theft","Beddington 1975; Crowley and Martin 1989","Exact derivation; one-parameter family possibly new"],
["Diversity","Tsallis = ln_q(Hill); combines by dial 1 − q; wall = ceiling 1/(q − 1)","Tsallis 1988; Jost 2006","Exact, known"],
["Multiple stressors","Additive null dial 0, multiplicative null 'or'","Folt 1999; Crain 2008","Exact, known"],
["Colonisation","wpow(−1, c, n)","Textbook","Relabelling"],
],
lean_file="biology/10-ecology/lean/Ecology.lean",
checks="11 of 11 pass",
sources=[
("Skalski and Gilliam 2001, Functional responses with predator interference: viable alternatives to the Holling type II model, Ecology 82:3083",None),
("Crowley and Martin 1989, Functional responses and interference within and between year classes of a dragonfly population, J. N. Am. Benthol. Soc. 8:211",None),
("Jost 2006, Entropy and diversity, Oikos 113:363",None),
("Cell biology page in this series (drug-combination nulls)","https://claude.ai/artifact/MD9HyXQFA5M5ai5okJdWEX"),
],
footnote="All feeding data on this page are simulated; the noise levels (10% and 25%) are my assumptions, not taken from experiments. The descriptions of Beddington–DeAngelis and Crowley–Martin, and Skalski and Gilliam's finding, come from a web search; the other citations are from memory.",
)
