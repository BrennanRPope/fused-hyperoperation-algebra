SPEC=dict(
title="Cell biology and the circle operation",
h1="Drug combinations: the standard interaction number is the dial, with the Highest Single Agent model at its wall",
lede="Chou and Talalay's combination index, v₁ + v₂ + α·v₁v₂, is exactly v₁ ∘<sub>α</sub> v₂. On that dial α = 0 is Loewe additivity and α = −1 is the Highest Single Agent model, whose rule that one drug makes the other irrelevant is the library's wall theorem. Bliss independence is plain ∘₁ on effect odds. I fitted α to three real mixture experiments on algae and duckweed: one sits at Loewe, two sit near the HSA end.",
stats=[("17 / 17","numeric checks"),("13","Lean theorems, no sorry"),("3","real mixture data sets placed on the dial")],
hero_fig="figures/f1_combination_dial.png",
hero_cap="Left: the iso-effect lines of the combination index at four dial settings. α = −1 is the square corner, the Highest Single Agent model. Right: α fitted to three published binary-mixture experiments, with 95% ranges from the likelihood profile.",
known="Everything named here is known pharmacology and toxicology. Loewe additivity (1926), Bliss independence (1939), the Highest Single Agent model, and Chou and Talalay's combination index with its interaction term (1984) are standard; so is the multi-target model of cell kill. Other one-parameter families also run from Loewe to HSA, for instance Hewlett's model, which the data authors used for one of these sets. What this page adds is the placement: the combination index is the weighted circle operation, α = −1 is HSA, and HSA's defining property is the wall.",
lead="Not a strong one, so I noted it here rather than following it further. The placement is exact and gives one number that puts HSA, Loewe and Chou's 'nonexclusive' case on a single line, and real data land at sensible places on it. But drug-synergy scoring already has several one-parameter families that do the same job, so the dial is a tidy common language, not a new prediction. Bliss being plain ∘₁ on effect odds is exact and pleasing, but it is a rewriting of 'survivals multiply'.",
kicker="Biology series · 3 of 11",
sections=[
("ci","The dial","1. The combination index is the weighted circle operation",
"""<p>Mix two drugs. For the effect the mixture reaches, write v₁ for the dose of drug 1 in the mix divided by the dose of drug 1 alone that gives that same effect, and v₂ likewise. Chou and Talalay's combination index is</p>
<div class="eq">CI = v₁ + v₂ + α·v₁v₂ = v₁ ∘<sub>α</sub> v₂</div>
<p>and the mixture is called additive when CI = 1. The interaction number α is the dial. Three settings carry names:</p>
<ul><li><b>α = 0:</b> Loewe additivity, the straight iso-effect line.</li>
<li><b>α = −1:</b> the index is 1 − (1 − v₁)(1 − v₂), the "or" law, so CI = 1 exactly when one drug alone is already at its own effective dose. That is the Highest Single Agent model: the mixture does as well as its better component and no better. I checked this on 300 random cases with unequal curve shapes; it holds exactly.</li>
<li><b>α = 1:</b> Chou's "mutually nonexclusive" case, plain ∘₁.</li></ul>
<p>The wall carries a meaning here. At v₂ = −1/α the index is −1/α whatever v₁ is. For α = −1 the wall is v₂ = 1: once drug 2 alone gives the effect, drug 1 is irrelevant. That sentence is the definition of HSA, and Lean proves it from the library's wall theorem.</p>
<p>One correction to a common shorthand: Chou's α = 1 is not Bliss independence, even when every curve has Hill slope 1. On random cases the two differ by up to 91% in effect odds.</p>"""),
("bliss","Bliss","2. Bliss independence is plain ∘₁ on effect odds",
"""<p>Bliss independence says the surviving fractions multiply: 1 − E<sub>AB</sub> = (1 − E<sub>A</sub>)(1 − E<sub>B</sub>). In effects that is the "or" law, dial −1. In effect odds o = E/(1 − E), since 1 + o = 1/(1 − E), it is</p>
<div class="eq">o<sub>AB</sub> = o<sub>A</sub> ∘₁ o<sub>B</sub></div>
<p>exactly, for any shape of dose–response curve (checked to 2×10⁻¹⁴, proved in Lean). With Hill slope 1, Loewe additivity is odds simply adding, dial 0 on the same coordinate. So on effect odds, Loewe and Bliss are dial 0 and dial 1, the same pair as the two two-nutrient growth laws on the metabolism page.</p>"""),
("data","Real data","3. Three real mixture experiments",
"""<p>The drc package ships three binary-mixture toxicity experiments by Cedergreen and colleagues, analysed by Sørensen et al. (2007): growth of a green alga under glyphosate and metsulfuron, and growth of duckweed under acifluorfen and diquat, and under mecoprop and terbuthylazine. Each uses fixed mixing ratios at a series of dilutions. I fitted log-logistic curves for the single compounds and one α for all the mixtures, by least squares, and got 95% ranges from the likelihood profile.</p>
<ul><li><b>Glyphosate + metsulfuron:</b> α = −0.32 (−0.75 to 0.45). Loewe is not rejected (p = 0.34), matching the data authors' note that concentration addition fits.</li>
<li><b>Acifluorfen + diquat:</b> α = −0.94 (−0.95 to −0.85). Loewe is rejected (p ≈ 10⁻¹¹) and so is HSA itself (p = 0.005): just above the HSA end. The mixing ratios agree with each other, α from −1.06 to −0.80.</li>
<li><b>Mecoprop + terbuthylazine:</b> α = −0.98 (−1.00 to −0.80). Loewe is rejected (p ≈ 2×10⁻⁵); HSA is not (p = 0.49).</li></ul>
<p>Bliss independence is rejected in all three. Adding a lower plateau to every curve moves each α by 0.01 at most. For the alga the separate mixing ratios scatter widely (−1.0 to 0.4), so a single α there is a summary, not a law.</p>
<figure><img src="figures/f2_duckweed_fit.png" alt="Duckweed growth rate against total dose for two mixing ratios, with the Loewe fit and the fitted dial."><figcaption>Acifluorfen + diquat at two mixing ratios. The fitted dial (α = −0.94) sits slightly right of Loewe: the mixture is weaker than dose addition predicts.</figcaption></figure>"""),
("kill","Cell kill","4. Cell kill and radiation",
"""<p>In the multi-target model a cell dies only if all n targets are hit, so it survives if any target survives: survival 1 − (1 − e<sup>−D/D₀</sup>)ⁿ is the dial −1 power of the per-target survival. That is exact and textbook.</p>
<p>The linear-quadratic radiation model is a near miss worth recording. Splitting a dose d₁ + d₂ into two fractions saves 2β·d₁d₂ of log-kill, a corner term, and incomplete repair scales it between 0 and 1. But no fixed dial k makes the split effect a weighted circle operation of the two single-fraction effects: the k needed runs from 0.55 to 0.17 across ordinary doses. It is a quadratic, not the circle operation.</p>"""),
],
table_h2="Where the circle operation shows up in cell biology",
table=[
["Combination index","CI = v₁ ∘<sub>α</sub> v₂; α is the dial","Chou and Talalay 1984","Exact, known"],
["Highest Single Agent","α = −1; its 'other drug irrelevant' rule is the wall","Long-standing reference model","Exact placement"],
["Bliss independence","Plain ∘₁ on effect odds, any curve shape","Bliss 1939","Exact, known"],
["Real mixtures","α = −0.32 (alga), −0.94 and −0.98 (duckweed)","Data: Cedergreen; Sørensen et al. 2007","Measured dial values"],
["Multi-target kill","Survival is wpow(−1)","Textbook radiobiology","Exact, known"],
["LQ split dose","Corner term but no fixed dial","Textbook radiobiology","Negative result"],
],
lean_file="biology/03-cell/lean/CellBiology.lean",
checks="17 of 17 pass",
sources=[
("Chou and Talalay 1984, Quantitative analysis of dose-effect relationships, Adv. Enzyme Regul. 22:27",None),
("Sørensen, Cedergreen, Skovgaard and Streibig 2007, An isobole-based statistical model and test for synergism/antagonism in binary mixture toxicity experiments, Environ. Ecol. Stat. 14:383",None),
("Data: drc package data sets glymet, acidiq and mecter (via Rdatasets)","https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/drc/acidiq.csv"),
("Metabolism page in this series","https://claude.ai/artifact/8sJ7P1mBE6PeqqNEDg7JhQ"),
],
footnote="Which compound the data sets' percentage column refers to is not stated in the documentation I read, so the figure labels give the mixing ratio only. α is symmetric in the two compounds, so the fitted values do not depend on that labelling.",
)
