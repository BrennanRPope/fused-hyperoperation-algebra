import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
def e_q(x, q):
    if abs(q-1) < 1e-12: return math.exp(x)
    b = 1 + (1-q)*x
    return b**(1/(1-q)) if b > 0 else 0.0
ts = [0.02*i for i in range(0, 201)]
cols = {0: 'var(--red)', 0.5: 'var(--yel)', 1: 'var(--blue)', 2: 'var(--green)', 3: 'var(--ink)'}
ser = [dict(pts=[[t, e_q(-t, n)] for t in ts], c=c, name=f"order {n}", short=f"{n}", dy={0:20,0.5:4,1:-12,2:4,3:14}[n]) for n, c in cols.items()]
picture = "plot(document.getElementById('pic'),{x:[0,4],y:[0,1],ml:52,mr:110,xlab:'time × rate constant (A₀ = 1)',ylab:'reactant left = Tsallis q-exponential, q = order',series:%s,extra:({X,Y})=>`<text x=\"${X(1.06)}\" y=\"${Y(0)-8}\">order 0: gone at 1</text><text x=\"${X(2)}\" y=\"${Y(0)-28}\" text-anchor=\"middle\">order ½: gone at 2</text>`});" % json.dumps(ser)

ks = [-3.5 + 0.05*i for i in range(0, 71)]
osc = [[-n, n*(n-1)/2] for n in [0.05*i for i in range(0, 71)]]
land_js = """
plot(document.getElementById('landscape-fig'),{x:[-3.6,0.6],y:[-0.5,4.2],ml:44,mr:20,xticks:[-3,-2,-1,0],xlab:'k',ylab:'d   (landscape law (a + b + k·ab)/(1 − d·ab))',
 series:[{pts:%B%,c:'var(--yel)',w:4},{pts:%O%,c:'var(--blue)',dash:'6 4',w:2},
 {pts:[[0,0],[-1,0],[-2,1]],c:'var(--blue)',line:false,dots:true,r:7,labels:[{x:0,y:0,t:'order 0',dx:-8,dy:-12,a:'end'},{x:-1,y:0,t:'order 1',dx:0,dy:22,a:'middle'},{x:-2,y:1,t:'order 2 (exact, border)',dx:10,dy:16}]},
 {pts:[[-3,3]],c:'var(--red)',line:false,dots:true,r:7,labels:[{x:-3,y:3,t:'order 3: osculates only',dx:10,dy:4}]}],
 extra:({X,Y})=>`<text x="${X(-3.3)}" y="${Y(3.9)}">turning (d &gt; k²/4)</text><text x="${X(-2.9)}" y="${Y(0.35)}">two walls</text><text x="${X(-3.55)}" y="${Y(1.25)}" fill="var(--yel)">border d = k²/4</text>`});
""".replace("%B%", json.dumps([[k, k*k/4] for k in ks if k*k/4 <= 4.2])).replace("%O%", json.dumps(osc))
cat = r["catalyst_700K"]; hot = r["hot_vs_463"]
sections = [
 dict(id="first", nav="First order", h="First order: consecutive times and parallel channels both compose at k = −1",
  body="""<p>For a first-order reaction the fraction left after time t is e<sup>−kt</sup>, so the fractions left over consecutive times multiply and the conversions combine by failure-union, X(t₁ + t₂) = X(t₁) ⊞₋₁ X(t₂). Two parallel channels add their rate constants (the k = 0 law), so again the surviving fractions multiply and the conversions combine at k = −1. The branching ratio between the channels stays fixed at every time.</p>"""),
 dict(id="second", nav="Second order", h="Second order is an exact landscape law on the border",
  body="""<p>For 2A → products, 1/[A] grows linearly in time, so the conversion's odds X/(1 − X) add over consecutive times. Composing two conversions therefore gives</p>
<div class="eq">X(t₁ + t₂) = (X₁ + X₂ − 2X₁X₂)/(1 − X₁X₂)</div>
<p>which is the landscape law at (k, d) = (−2, 1). That point is on the border k² = 4d, the same country as the flat triangle and Lavoisier's mass law. Zeroth order composes by plain addition, (0, 0), and first order by (−1, 0). All three are checked against numerical integration, and Lean proves the second-order identity.</p>"""),
 dict(id="landscape", nav="Order n", h="Order n is a Tsallis q-exponential, and only orders 0, 1 and 2 are landscape laws",
  body=f"""<p>For a reaction of order n the fraction left is the Tsallis q-exponential e<sub>q</sub>(−kt) with q = n, and consecutive times compose by Tsallis's q-product. Behind this sits a clean identity: the q-logarithm turns ordinary products into the dial, ln<sub>q</sub>(xy) = ln<sub>q</sub>x ⊞<sub>1−q</sub> ln<sub>q</sub>y.</p>
<p>Is every order a landscape law? No. Fitting the landscape form to third-order kinetics needs a q that drifts from {r['q_order3'][0]:.2f} to {r['q_order3'][1]:.2f} across conversions, so it is not one law. Near zero conversion, though, order n always touches the landscape at (−n, n(n − 1)/2), which gives exactly (0, 0), (−1, 0) and (−2, 1) for the first three orders. That curve's discriminant is n(2 − n): two walls for orders between 0 and 2, the border at 0 and 2, and the turning country above 2. This is a small curiosity, not a lead: away from zero conversion it describes nothing.</p>
<p>The wall (all reactant gone) is reached in finite time exactly when n &lt; 1: order ½ empties at kt = 2, order 0 at kt = 1, while orders 1 and up only approach it (picture at the top).</p>""",
  fig=True, figcls="", cap="The landscape plane. Orders 0, 1 and 2 are exact landscape laws (blue dots), with order 2 on the border. Other orders only touch the dashed curve (−n, n(n−1)/2) near zero conversion; order 3 (red) is shown where it touches.", js=land_js),
 dict(id="arrhenius", nav="Arrhenius and N₂", h="Arrhenius: barriers add in the exponent, so rates multiply, and nitrogen fixation needs both heat and a catalyst",
  body=f"""<p>The rate constant k = A·e<sup>−Ea/RT</sup> puts barriers and inverse temperatures in the exponent, so rate ratios over successive changes multiply. At 700 K a catalyst that lowers a barrier by 10, 50 or 100 kJ/mol multiplies the rate by {cat['10']:.1f}, {cat['50']:.1e} and {cat['100']:.1e}.</p>
<p>This closes the nitrogen-fixation thread. Thermodynamics (subject 7) wants ammonia made below about 463 K, but the N≡N bond's hard first step (subject 4) is far too slow there. Raising the temperature to 700 K speeds a 100 kJ/mol step by about {hot['100']:.0e}, at the cost of equilibrium yield, which the high pressure (subject 8) wins back. The iron catalyst lowers the barrier so that this compromise works at all.</p>"""),
]
claims = [
 ("First order and parallel channels: conversions compose at k = −1", "tidy reading", "numerical integration; Lean <code>first_order_compose</code>, <code>parallel_channels</code>"),
 ("Second order: conversions compose by the landscape law at (−2, 1), on the border", "new here", "numerical integration; Lean <code>second_order_lens</code>, <code>second_order_border</code>"),
 ("Order n survival is the Tsallis q-exponential with q = n; times compose by the q-product", "known", "orders ½ to 3"),
 ("ln_q(xy) = ln_q x ⊞_(1−q) ln_q y", "known", "Python; Lean <code>qlog_two</code>, <code>qlog_zero</code>"),
 ("Order 3 is not a landscape law; order n osculates (−n, n(n−1)/2) with discriminant n(2 − n)", "new here", "Python; Lean <code>osculating_disc</code>, <code>osculating_small</code>"),
 ("Order n < 1 empties in finite time", "known", "Python; Lean <code>half_order_wall</code>"),
 ("Catalyst factors e^(ΔE/RT) at 700 K; 700 K vs 463 K for a 100 kJ/mol step", "known", "Python; Lean <code>arrhenius_steps</code>"),
]
lean = [("first_order_compose, parallel_channels, first_order_point","first order at k = −1"),("landNum, landDen, second_order_lens","second order is the landscape law (−2, 1)"),
 ("second_order_border","k² = 4d there"),("osculating_disc, osculating_small","the osculating curve and its discriminant"),
 ("qlog_two, qlog_zero","q-logarithms turn products into the dial"),("half_order_wall","order ½ empties at t = 2/k"),("arrhenius_steps","rate ratios multiply")]
n = build("index.html", title="Rates and Arrhenius", tile=("9","k","Rates"),
 h1="First order composes at k = −1, second order is an exact landscape law on the border, and nitrogen fixation needs heat plus a catalyst",
 lead="Consecutive times and parallel channels combine first-order conversions by failure-union. Second-order conversions compose by (a + b − 2ab)/(1 − ab), the landscape point (−2, 1) on the border. Order n is a Tsallis q-exponential; only orders 0, 1 and 2 are landscape laws. Arrhenius puts barriers in the exponent, which is why catalysts and heat matter for ammonia.",
 figs=[("19 / 19","checks pass"),("11","Lean theorems, compiled"),("(−2, 1)","second order: border"),("×6.6e3","700 K vs 463 K, 100 kJ/mol step")],
 picture=picture, caption="Fraction of reactant left against time for orders 0 to 3, each a Tsallis q-exponential with q equal to the order. Orders below 1 hit zero in finite time; order 1 and above only approach it.",
 verdict_h="Lead check: not a lead",
 verdict="Second-order kinetics composing by an exact border law is a clean reading I haven't seen stated, but it is the familiar fact that 1/[A] grows linearly, so it predicts nothing new. The osculating curve for order n is a curiosity that describes nothing away from zero conversion. Neither is worth following.",
 sections=sections, claims=claims, lean_file="chemistry/09-rates/lean/ChemRates.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/09-rates/: <code>checks.py</code> (19 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemRates.lean</code>, <code>make_page.py</code>.",
 sources='Sources: Tsallis, J. Stat. Phys. 52, 479 (1988) and later work on q-exponentials; Arrhenius (1889); standard kinetics texts for integrated rate laws.',
 series_nav="Chemistry pass, part 9 of 13: rates. Before: equilibrium. Next: electrochemistry, with photosynthesis.")
print("page bytes", n)
