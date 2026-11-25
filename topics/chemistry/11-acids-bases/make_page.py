import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
curve = r["phos_curve"]
picture = "plot(document.getElementById('pic'),{x:[0,3],y:[0,14],ml:52,mr:20,xticks:[0,0.5,1,1.5,2,2.5,3],yticks:[0,2,4,6,8,10,12,14],xlab:'protons removed per H₃PO₄',ylab:'pH',series:%s,extra:({X,Y})=>`%s`});" % (
 json.dumps([dict(pts=curve, c='var(--blue)', w=3, name='', short=''),
             dict(pts=[[0.5, 2.14], [1.5, 7.20], [2.5, 12.34]], c='var(--red)', line=False, dots=True, r=6, name='', short='')]),
 "".join(f'<text x="${{X({x}+0.08)}}" y="${{Y({p})+5}}">pKa {p}</text>' for x, p in ((0.5, 2.14), (1.5, 7.20), (2.5, 12.34))))
g = r["gaps"]
pts = [[g[n][0], g[n][1]] for n in ("citric", "phosphoric", "arsenic")]
labs = [dict(x=g["citric"][0], y=g["citric"][1], t="citric", dx=10, dy=4), dict(x=g["phosphoric"][0], y=g["phosphoric"][1], t="phosphoric", dx=-10, dy=-10, a="end"),
        dict(x=g["arsenic"][0], y=g["arsenic"][1], t="arsenic", dx=10, dy=16)]
gap_js = """
plot(document.getElementById('gaps-fig'),{x:[0,6],y:[0,6],ml:44,mr:20,xticks:[0,1,2,3,4,5,6],yticks:[0,1,2,3,4,5,6],xlab:'lower gap pKa₂ − pKa₁',ylab:'upper gap pKa₃ − pKa₂',
 series:[{pts:[[0,0],[6,6]],c:'var(--green)',w:2.5},{pts:[[0.477,0],[0.477,6]],c:'var(--red)',w:2.5,dash:'6 4'},
 {pts:%P%,c:'var(--blue)',line:false,dots:true,r:7,labels:%L%}],
 extra:({X,Y})=>`<text x="${X(2.6)}" y="${Y(2.0)}" fill="var(--green)">pairwise: equal gaps</text><text x="${X(0.6)}" y="${Y(5.6)}" fill="var(--red)">dial fold: lower gap = log 3</text>`});
""".replace("%P%", json.dumps(pts)).replace("%L%", json.dumps(labs))
kp = r["k_pair"]; ek = r["pair_kJ"]; pr = r["pKa3_pred"]
sections = [
 dict(id="logs", nav="pH and logs", h="pH turns the lift into addition, and independent groups fold at k = −1",
  body=f"""<p>Chaining dissociation steps multiplies their constants, so pKa values add and K − 1 composes by the circle operation (called mult-add before 1.0), as in subject 7. The fraction of an acid deprotonated is a logistic curve in pH (Henderson–Hasselbalch). For two independent groups, the chance that at least one has lost its proton is α₁ ⊞₋₁ α₂, the failure-union law again.</p>
<p>The isoelectric point is an average, not a circle operation: glycine (pKa 2.34 and 9.60) is neutral at pH {r['glycine_pI']:.2f}, exactly their mean, because h² = K₁K₂ there.</p>"""),
 dict(id="buffers", nav="Buffers", h="Buffers: capacities add, and mixing combines ratios by the mediant",
  body="""<p>Two buffers in one solution have a total capacity equal to the sum of their separate terms at the shared pH, the plain addition law (k = 0). Mixing two buffers made from one acid combines their base-to-acid ratios as the mediant (b₁ + b₂)/(a₁ + a₂), which always lies between the two. The full equilibrium calculation agrees with the mediant formula to 0.001 pH units in the test case. Neither of these is the circle operation.</p>"""),
 dict(id="gaps", nav="Triprotic test", h="A real test: triprotic acids follow pairwise interactions, not the dial folded over three sites",
  body=f"""<p>Subject 8 showed that folding the dial over n sites gives the MWC model with exclusive binding. For three sites that fold makes the protonation polynomial 1 + 3x + 3kx² + k²x³, and a short calculation shows that this pins the lower gap pKa₂ − pKa₁ at log 3 = 0.48, whatever k is. A uniform interaction between every pair of protons instead gives 1 + 3x + 3kx² + k³x³, which makes the two gaps equal.</p>
<p>Real acids pick the second. From one standard table: citric acid has gaps 1.63 and 1.64, phosphoric 5.06 and 5.14, arsenic 4.72 and 4.59 (arsenic values vary by about 0.2 between sources). Each pair is close to equal and every gap is far from 0.48. Fitting the pairwise factor gives about {ek['citric']:.1f} kJ/mol of repulsion per proton pair for citric acid and {ek['phosphoric']:.0f} kJ/mol for phosphoric. Using the first two pKas to predict the third gives {pr['citric']:.2f} for citric (table 6.40), {pr['phosphoric']:.2f} for phosphoric (12.34) and {pr['arsenic']:.2f} for arsenic (11.52).</p>
<p>So the dial's many-site fold, which matched MWC exactly, does not describe polyprotic acids. That is a useful negative: it marks where the fold stops applying. The pairwise picture that does fit is the classical Bjerrum and Kirkwood–Westheimer electrostatic view.</p>""",
  fig=True, figcls="", cap="Lower and upper pKa gaps for three triprotic acids (blue). The green diagonal is the uniform pairwise prediction, equal gaps; the dashed red line is the three-site dial fold, which fixes the lower gap at log 3.", js=gap_js),
]
claims = [
 ("Chained Ka multiply; independent groups combine at k = −1", "tidy reading", "Python; Lean <code>chained_K</code>, <code>indep_groups</code>"),
 ("Isoelectric point is where h² = K₁K₂ (glycine 5.97)", "known", "Python; Lean <code>isoelectric</code>"),
 ("Buffer capacities add; mixing buffers combines ratios by the mediant", "known", "full equilibrium solve; Lean <code>mediant_above</code>, <code>mediant_below</code>"),
 ("Three-site dial fold pins pKa₂ − pKa₁ at log 3 for every k", "new here", "Python; Lean <code>fold_three</code>, <code>fold_gap_log3</code>, <code>fold_unequal_gaps</code>"),
 ("Uniform pairwise interaction gives equal pKa gaps", "known", "Python; Lean <code>pairwise_equal_gaps</code>"),
 ("Citric, phosphoric and arsenic acids have near-equal gaps, far from log 3", "checked against a table", "tabulated pKas"),
]
lean = [("indep_groups, chained_K","k = −1 for independent groups, the circle operation for chained constants"),("isoelectric","h² = K₁K₂ at the isoelectric point"),
 ("mediant_above, mediant_below","the buffer-mixing mediant lies between the ratios"),("fold_three","the three-site dial fold's coefficients 3, 3k, k²"),
 ("fold_gap_log3, fold_unequal_gaps","the fold pins one gap at log 3 and breaks gap equality"),("pairwise_equal_gaps","pairwise interaction gives equal gaps")]
n = build("index.html", title="Acids and bases", tile=("11","pH","Acids"),
 h1="Acids and bases add in logs and fold at k = −1, and triprotic acids reject the dial's many-site fold in favour of pairwise interactions",
 lead="pH turns the lift into addition, independent groups combine by failure-union, and buffer mixing is a mediant. The real test is triprotic acids: folding the dial over three sites would pin one pKa gap at log 3, while pairwise interactions make the two gaps equal. Citric, phosphoric and arsenic acids all have near-equal gaps.",
 figs=[("14 / 14","checks pass"),("9","Lean theorems, compiled"),("1.63 / 1.64","citric acid pKa gaps"),("0.48","what the dial fold would force")],
 picture=picture, caption="Titration of phosphoric acid: pH against protons removed. The three half-equivalence points sit at the three pKas, spaced by nearly equal gaps of about 5.",
 verdict_h="Lead check: not a lead",
 verdict="The triprotic test is clean and decisive, but it goes against the dial: its many-site fold does not describe polyprotic acids, and the picture that fits is classical pairwise electrostatics. The rest restates known log and mixing rules. I noted the negative result and didn't follow it further.",
 sections=sections, claims=claims, lean_file="chemistry/11-acids-bases/lean/ChemAcids.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/11-acids-bases/: <code>checks.py</code> (14 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemAcids.lean</code>, <code>make_page.py</code>.",
 sources='Sources: pKa values from a general-chemistry ionization-constant appendix (citric 3.13, 4.76, 6.40; phosphoric 2.14, 7.20, 12.34; arsenic 2.21, 6.93, 11.52; arsenic varies between tables); Bjerrum (1923) and Kirkwood and Westheimer (1938) for the electrostatic pairwise picture; Monod, Wyman and Changeux (1965).',
 series_nav="Chemistry pass, part 11 of 13: acids and bases. Before: electrochemistry. Next: organic chemistry.")
print("page bytes", n)
