SPEC=dict(
title="Neuroscience and the circle operation",
h1="On a passive dendrite, two synapses combine by an exact lens law, and the measured shunting rule is its dial",
lede="For any passive linear dendritic tree, the soma's response to two conductance synapses is exactly V₁₂ = (a + b + k·ab)/(1 − d·ab), where a and b are the responses to each alone and k and d depend only on transfer impedances and reversal potentials. The rule Hao and colleagues measured in hippocampal neurons in 2009, V ≈ EPSP + IPSP + k·EPSP·IPSP, is the dial part of that law. In one compartment it is (a + b − 2ab)/(1 − ab), and the reversal potential is the wall.",
stats=[("15 / 15","numeric checks"),("17","Lean theorems, no sorry"),("0","leads"),],
hero_fig="figures/f1_dendrite_dial.png",
hero_cap="Model dendrite (60-node trunk, 25-node oblique branch), inhibition 10 nodes from the soma. Left: for four inhibition strengths and a range of excitations, the amount the combined response falls short of the linear sum lies on one line, slope k, whatever the strengths: the dial. Right: k for excitation at each node. Beyond the inhibition k is nearly flat at the geometric value Z_II/(|E_I|·Z_sI); between the inhibition and the soma it falls. Inhibition on the oblique branch shunts excitation in that branch about three times harder than excitation on the trunk.",
known="The physics is classical steady-state cable theory. Koch, Poggio and Torre worked out shunting interactions on dendrites in the early 1980s, Gidon and Segev (2012) gave exact transfer-resistance formulas for how much inhibition shunts at each site, and Hao, Wang, Dan, Poo and Zhang (2009) measured the bilinear rule and derived it from a two-port analysis. I could not read Hao's derivation (their equation is an image I could not open), so they may already have the exact form. What this page adds is the packaging: the exact law as a lens law with closed-form k and d, the factorisation that turns summation into multiplication, and Lean proofs from the circuit equations.",
lead="None worth following. The connection is real and exact, and the dial appears in the experimental literature in exactly Brennan's form a + b + k·ab, which no other page in this series has had. But the content is known circuit theory, and the one new prediction, the d term beyond the dial, changes the combined response by at most 0.04 mV at realistic amplitudes in my model. That is below what a patch-clamp recording can resolve, so it gives no experiment to run.",
kicker="Biology series · 7 of 11",
sections=[
("one","One compartment","1. One compartment: the lens law with k = −2, d = 1",
"""<p>Measure voltage in units of the synaptic reversal potential and conductance in units of the leak. A synapse of conductance g holds the cell at v = g/(1 + g). Two synapses together give (g₁ + g₂)/(1 + g₁ + g₂), and in terms of the separate responses that is</p>
<div class="eq">v₁₂ = (v₁ + v₂ − 2v₁v₂)/(1 − v₁v₂)</div>
<p>a lens law with k = −2 and d = 1. In odds v/(1 − v) the conductances come back and simply add, dial 0. The reversal potential v = 1 is the wall: a synapse that clamps the cell at its reversal makes every other input irrelevant. All proved in Lean.</p>"""),
("tree","Any tree","2. Any passive tree: an exact lens law with geometric constants",
"""<p>Take any passive linear network, two synapses at sites i and j with reversal potentials E<sub>i</sub>, E<sub>j</sub>, and the soma at s. Solving the two-port equations gives, exactly,</p>
<div class="eq">V₁₂ = (a + b + k·ab)/(1 − d·ab)<br>k = −Z<sub>ij</sub>(1/(E<sub>j</sub>Z<sub>si</sub>) + 1/(E<sub>i</sub>Z<sub>sj</sub>)), &nbsp; d = Z<sub>ij</sub>²/(E<sub>i</sub>E<sub>j</sub>Z<sub>si</sub>Z<sub>sj</sub>)</div>
<p>where Z are transfer impedances. The constants depend on where the synapses are and on their reversal potentials, never on how strong they are. I checked this on 200 random trees with random sites, reversals and strengths; the largest relative error was 1 × 10⁻¹³. The Lean file proves it from the circuit equations.</p>
<p>Writing k = −(A + B) and d = A·B, the law factorises: 1 − A·V₁₂ = (1 − Aa)(1 − Ab)/(1 − ABab), and the same for B. So in the coordinate T(v) = (1 − Av)/(1 − Bv) the two responses simply multiply, and T − 1 combines by plain ∘₁. With one synapse at the soma, B is 1 in reversal units, the law sits on the line k + d = −1, and that synapse's reversal is again an absorbing wall.</p>"""),
("hao","Experiment","3. The measured rule is the dial",
"""<p>Hao and colleagues (2009) paired excitatory and inhibitory inputs on CA1 pyramidal cells and found the soma response is well described by EPSP + IPSP + k·EPSP·IPSP, with k about 0.14 per mV, the same across amplitudes, roughly constant when excitation was farther out than inhibition and smaller when it was nearer the soma. That is the dial a + b + k·ab, and it is the lens law with the small d term dropped.</p>
<p>On a model dendrite with EPSPs of 0.2–5 mV and IPSPs of 0.3–2.6 mV, the best bilinear fit gives k = 0.2172 per mV and the exact formula gives 0.2172 (I did not fit the model to their cells, so the value differs from theirs). The location pattern follows from the formula: for excitation anywhere beyond the inhibitory site, the dominant part of k is Z<sub>II</sub>/(|E<sub>I</sub>|·Z<sub>sI</sub>), which does not depend on where the excitation is. In the model k varied 4% over the whole distal trunk and fell to 65% of that next to the soma. Inhibition on an oblique branch gave k of 0.59 per mV for excitation in that branch and at most 0.18 on the trunk, the branch-local shunting Hao reported.</p>
<p>What the dial leaves out is the 1 − d·ab denominator. In the model that term is 1.7% of the shunting term and at most 0.04 mV of the response.</p>"""),
("other","Smaller","4. Two small placements",
"""<p>The chance that at least one of n release sites at a synapse fires, each with probability p, is 1 − (1 − p)ⁿ, the "or" power wpow(−1). Weber's law makes n just-noticeable steps multiply a stimulus by (1 + w)ⁿ, a bare circle-operation power. Both are checked and proved and neither says anything new.</p>"""),
],
table_h2="Where the circle operation (called mult-add before 1.0) shows up in neuroscience",
table=[
["One-compartment summation","Lens law k = −2, d = 1; odds add; reversal is the wall","Textbook","Exact, known"],
["Two synapses on any passive tree","Exact lens law, k and d from transfer impedances; T-coordinate multiplies","Classical cable theory","Exact; closed form packaged here"],
["Hao et al. shunting rule","The dial a + b + k·ab, first-order part of the lens law","Hao et al. 2009","Exact correspondence; d term too small to measure"],
["Location dependence of k","Distal k = Z_II/(|E_I|·Z_sI), flat; falls toward the soma","Hao et al. 2009; Gidon and Segev 2012","Follows from the formula"],
["Release sites, Weber","'Or' power; bare circle-operation power","Textbook","Relabelling only"],
],
lean_file="biology/07-neuroscience/lean/NeuroScience.lean",
checks="15 of 15 pass",
sources=[
("Hao, Wang, Dan, Poo and Zhang 2009, An arithmetic rule for spatial summation of excitatory and inhibitory inputs in pyramidal neurons, PNAS 106:21906","https://pmc.ncbi.nlm.nih.gov/articles/PMC2799885/"),
("Koch, Poggio and Torre 1983, Nonlinear interactions in a dendritic tree, PNAS 80:2799",None),
("Gidon and Segev 2012, Principles governing the operation of synaptic inhibition in dendrites, Neuron 75:330",None),
("Cell biology page in this series (drug combinations, the same dial)","https://claude.ai/artifact/MD9HyXQFA5M5ai5okJdWEX"),
],
footnote="The model dendrite is a passive compartmental cable I chose for Hao-like amplitudes (axial 30, leak 0.02, soma 0.3, in arbitrary conductance units; excitatory reversal 65 mV and inhibitory −5 mV above rest). It is not fitted to any cell, so its k is not theirs. Hao's k ≈ 0.14 per mV and the location pattern are from their paper, which I read; the Koch and Gidon references are from memory and I did not open them.",
)
