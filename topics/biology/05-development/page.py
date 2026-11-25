SPEC=dict(
title="Developmental biology and the circle operation",
h1="Growth curves are one circle-operation power of a fading excess, and logistic growth is the inverse again",
lede="The Richards family of growth curves, which holds the logistic, Gompertz and von Bertalanffy curves as special cases, is exactly size/asymptote = 1 + mpow(x, p), where the excess x fades exponentially and p = −1/ν. Logistic growth is p = −1, the circle-operation inverse, the same saturation as an enzyme's Michaelis–Menten curve. Gompertz growth in discrete time iterates a single circle-operation power. Three real growth data sets put the exponent in three different places.",
stats=[("12 / 12","numeric checks"),("8","Lean theorems, no sorry"),("0","leads"),],
hero_fig="figures/f1_growth_family.png",
hero_cap="Left: four members of the Richards family, all of the form 1 + mpow(x, −1/ν). Right: three real growth data sets with the best Richards curve and the logistic and Gompertz fits. Soybean and loblolly reject the named shapes; the orange trees cannot tell them apart.",
known="All of this is known. The Richards curve (1959) unifies the logistic, Gompertz and von Bertalanffy growth laws through one shape parameter ν, and allometry is Huxley's power law. Writing the family as 1 + mpow(x, −1/ν) is a relabelling of Richards' own formula. What the page adds is the link to the rest of the series: the logistic case is the same circle-operation inverse that gives enzyme saturation, and the Gompertz case is an iterated circle-operation power.",
lead="None. The placement is exact but it is a relabelling: the shape parameter ν already plays the role of the exponent, and the data show no preferred value. It is noted here as glue for the series, not followed.",
kicker="Biology series · 5 of 11",
sections=[
("family","Growth curves","1. One power, many growth laws",
"""<p>The Richards growth curve is y = K·(1 + Q·e<sup>−rt</sup>)<sup>−1/ν</sup>. Put x = Q·e<sup>−rt</sup>, an excess that fades exponentially. Then</p>
<div class="eq">y / K = 1 + mpow(x, −1/ν)</div>
<p>a single circle-operation power. The named growth laws are settings of the power:</p>
<ul><li><b>Logistic (ν = 1, power −1):</b> y/K = 1/(1 + x) = 1 + inv(x). This is the circle-operation inverse, the same saturation as Michaelis–Menten kinetics on the metabolism page.</li>
<li><b>Von Bertalanffy (ν = −⅓, power 3):</b> y/K = (1 + x)³ with x between −1 and 0.</li>
<li><b>Gompertz (ν → 0):</b> the limit, checked numerically to 10⁻⁴.</li></ul>
<p>In discrete time Gompertz growth is y<sub>t+1</sub>/K = (y<sub>t</sub>/K)<sup>λ</sup>, so each step is the circle-operation power x ↦ (1 + x)<sup>λ</sup> − 1. Powers compose by multiplying exponents, so n steps are one power λⁿ. Lean proves both facts for whole-number exponents.</p>
<p>Allometry, Huxley's y ∝ x<sup>a</sup>, says the growth of y is the circle-operation power of the growth of x. Kleiber's ¾ law on the metabolism page is the same statement.</p>"""),
("data","Real data","2. Where real growth sits",
"""<p>I fitted the Richards curve to three classic growth data sets, pooling all individuals in each, and tested the three named shapes against the free fit:</p>
<ul><li><b>Soybean leaf weight:</b> ν = 2.54 (power −0.39). Logistic, Gompertz and von Bertalanffy are all rejected (p ≤ 0.005).</li>
<li><b>Loblolly pine height:</b> ν = −0.81 (power 1.24), past von Bertalanffy. All three named shapes are rejected.</li>
<li><b>Orange tree trunks:</b> ν = 0.93, but the data cannot tell the three shapes apart (all p > 0.4).</li></ul>
<p>So the exponent varies from data set to data set, as growth biologists already know, and nothing singles out a special value. Pooling different plots, years and trees mixes individuals with different curves, so these ν values describe the average curve, not each organism.</p>
<p>I also considered phyllotaxis and morphogen gradients and did not find a place where the circle operation (called mult-add before 1.0) applies: phyllotaxis is about continued fractions and the golden angle, and gradients from two sources simply add. I ran no computation on either.</p>"""),
],
table_h2="Where the circle operation shows up in development",
table=[
["Richards growth family","y/K = 1 + mpow(x, −1/ν), x fading exponentially","Richards 1959","Exact relabelling"],
["Logistic growth","Power −1: the circle-operation inverse, as in enzyme saturation","Verhulst 1838","Exact, known"],
["Gompertz growth","Iterated circle-operation power; steps compose to λⁿ","Gompertz 1825","Exact, known"],
["Allometry","Growth of y = mpow(growth of x, a)","Huxley 1932","Relabelling only"],
["Real growth data","ν = 2.54, −0.81, 0.93: no shared exponent","Standard data sets","No lead"],
],
lean_file="biology/05-development/lean/Development.lean",
checks="12 of 12 pass",
sources=[
("Richards 1959, A flexible growth function for empirical use, J. Exp. Bot. 10:290",None),
("Data: Soybean (nlme), Loblolly and Orange (R datasets), via Rdatasets","https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/datasets/Loblolly.csv"),
("Metabolism page in this series (saturation as the inverse)","https://claude.ai/artifact/8sJ7P1mBE6PeqqNEDg7JhQ"),
],
footnote="Fits are ordinary least squares on pooled data with three free parameters plus ν; p-values compare a fixed-ν fit against the free one with an F test. A mixed-effects fit per individual would be the careful version.",
)
