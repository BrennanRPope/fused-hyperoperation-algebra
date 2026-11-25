import sys, json, math
import numpy as np
sys.path.insert(0, "..")
from page_builder import build
from scipy.optimize import brentq
r = json.load(open("results.json"))
names = ["oxalic","malonic","succinic","glutaric","adipic","pimelic","suberic","azelaic","sebacic"]
om = r["diacid_omega"]
pts = [[i, v] for i, (n, v) in enumerate(om.items())]
picture = """
plot(document.getElementById('pic'),{x:[-0.5,8.8],y:[0.002,2],ylog:true,ml:56,mr:%MR%,xticks:[0,1,2,3,4,5,6,7,8],yticks:[0.002,0.01,0.1,1],yfmt:v=>String(v),xlab:'CH₂ groups between the two acid groups',ylab:'dial ω = 4·10^(−ΔpKa)',
 hlines:[{y:1,c:'var(--yel)',label:'k = 1: independent',short:'k = 1'}],
 series:[{pts:%P%,c:'var(--blue)',dots:true,r:6,labels:%L%}]});
""".replace("%P%", json.dumps(pts)).replace("%L%", json.dumps([{"x": p[0], "y": p[1], "t": names[i], "dx": 8, "dy": 16 if i % 2 else -8} for i, p in enumerate(pts) if i in (0, 1, 2, 4, 8)])).replace("%MR%", "230")

# zeros of 1 + 2x + w x^2 as w varies
zr = []
for i in range(0, 121):
    w = 10**(-1 + 0.02*i)
    for rt in np.roots([w, 2, 1]):
        zr.append([float(rt.real), float(rt.imag), w])
zeros_js = """
{const el=document.getElementById('zeros-fig');const Z=%Z%;const draw=()=>{const W=el.clientWidth,H=el.clientHeight;const m={l:40,r:16,t:20,b:40};
const x0=-3.2,x1=0.4;const sc=(W-m.l-m.r)/(x1-x0);const cy=(m.t+H-m.b)/2;const X=v=>m.l+(v-x0)*sc,Y=v=>cy-v*sc;
let s=`<line x1="${m.l}" x2="${W-m.r}" y1="${cy}" y2="${cy}" stroke="var(--mute)"/><line x1="${X(0)}" x2="${X(0)}" y1="${m.t}" y2="${H-m.b}" stroke="var(--mute)"/>`;
s+=`<circle cx="${X(-0.5)}" cy="${cy}" r="${0.5*sc}" fill="none" stroke="var(--red)" stroke-dasharray="4 4" opacity=".6"/>`;
[-3,-2,-1,0].forEach(v=>{s+=`<text x="${X(v)}" y="${H-m.b+18}" text-anchor="middle" style="font-size:12px">${v}</text>`});
Z.forEach(p=>{if(p[0]<x0)return;const c=p[2]<=1?'var(--blue)':'var(--red)';s+=`<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="2.6" fill="${c}"/>`});
s+=`<circle cx="${X(-1)}" cy="${cy}" r="7" fill="none" stroke="var(--yel)" stroke-width="3"/>`;
s+=`<text x="${X(-1)-10}" y="${cy-14}" text-anchor="end">ω = 1: double zero at −1</text>`;
s+=`<text x="${X(-2.2)}" y="${cy+22}" text-anchor="middle" fill="var(--blue)">ω &lt; 1: real zeros</text>`;
s+=`<text x="${X(-0.5)}" y="${Y(0.62)}" text-anchor="middle" fill="var(--red)">ω &gt; 1: on this circle</text>`;
s+=`<text x="${(m.l+W-m.r)/2}" y="${H-6}" text-anchor="middle" style="font-size:12px">complex plane of the zero x</text>`;
el.innerHTML=s};const ro=new ResizeObserver(draw);ro.observe(el);draw();}
""".replace("%Z%", json.dumps(zr))

def Yf(x, k, n, h=1e-7):
    P = lambda t: 1 + ((1 + k*t)**n - 1)/k
    return (P(x*(1+h)) - P(x*(1-h)))/(2*h)/(n*P(x))
def nH(k, n):
    x0 = brentq(lambda x: Yf(x, k, n) - 0.5, 1e-12, 1e6); h = 1e-4
    lg = lambda x: math.log(Yf(x, k, n)/(1 - Yf(x, k, n)))
    return (lg(x0*(1+h)) - lg(x0*(1-h)))/(math.log(1+h) - math.log(1-h))
kk = [10**(0.001 + 0.05*i) for i in range(0, 101)]
hill_js = """
plot(document.getElementById('mwc-fig'),{x:[1,1e5],y:[1,4],xlog:true,ml:44,mr:110,xticks:[1,10,100,1000,10000,100000],xlab:'dial k (MWC allosteric constant L = k − 1)',ylab:'Hill coefficient of the n-site fold',
 hlines:[{y:2,c:'var(--mute)'},{y:3,c:'var(--mute)'}],
 series:[{pts:%A%,c:'var(--blue)',name:'2 sites',short:'2'},{pts:%B%,c:'var(--green)',name:'3 sites',short:'3'},{pts:%C%,c:'var(--red)',name:'4 sites',short:'4'}]});
""".replace("%A%", json.dumps([[k, nH(k, 2)] for k in kk])).replace("%B%", json.dumps([[k, nH(k, 3)] for k in kk])).replace("%C%", json.dumps([[k, nH(k, 4)] for k in kk]))

hi = r["haber_ideal"]
haber_js = """
plot(document.getElementById('haber-fig'),{x:[0,320],y:[0,0.55],ml:56,mr:40,xticks:[0,50,100,150,200,250,300],yticks:[0,0.1,0.2,0.3,0.4,0.5],yfmt:v=>Math.round(v*100)+'%',xlab:'pressure, bar (700 K, 1 N₂ : 3 H₂ feed)',ylab:'equilibrium ammonia mole fraction (constant-ΔH model)',
 series:[{pts:%P%,c:'var(--blue)',dots:true},{pts:[[200,%V%]],c:'var(--red)',line:false,dots:true,r:7,labels:[{x:200,y:%V%,t:'with vdW fugacities',dx:-10,dy:-12,a:'end'}]}]});
""".replace("%P%", json.dumps([[float(P), v] for P, v in hi.items()])).replace("%V%", str(r["haber_200_vdw"]))
hill = r["hill4"]
sections = [
 dict(id="diacids", nav="Diacids", h="In a diacid the dial measures how much the two charges feel each other",
  body=f"""<p>A binding polynomial lists a molecule's states weighted by concentration. For two sites with interaction factor ω it is 1 + K₁h + K₂h + ωK₁K₂h², so P − 1 = (K₁h) ⊞<sub>ω</sub> (K₂h): the dial of subject 7 again.</p>
<p>Dicarboxylic acids measure ω directly. With two equal sites the gap between the two pKa values is log 4 − log ω, where log 4 = 0.60 is the gap two independent sites would show from counting alone. So ω = 4·10<sup>−ΔpKa</sup>. From Wikipedia's table it climbs from {om['oxalic']:.3f} for oxalic acid, where the two acid groups touch, to {om['sebacic']:.2f} for sebacic acid, eight CH₂ groups apart: from near-exclusion toward independence. It never passes 1, since like charges only repel. This is Bjerrum's 1923 electrostatic picture in dial form.</p>"""),
 dict(id="zeros", nav="Zeros and cooperativity", h="For two sites, cooperativity begins exactly where the zeros leave the real line",
  body="""<p>For two equal sites the binding polynomial is 1 + 2x + ωx², and three things change together at ω = 1, the circle operation (called mult-add before 1.0) setting:</p>
<ul><li><b>Hill coefficient:</b> it is 2√ω/(1 + √ω), above 1 exactly when ω &gt; 1.</li>
<li><b>Zeros:</b> they are real for ω ≤ 1 and leave the real line for ω &gt; 1, at distance 1/√ω from 0. Every such pair lies on one circle, the one through 0 and −1.</li>
<li><b>Dial:</b> k = ω passes the circle operation.</li></ul>
<p>Lean proves the no-real-zero half by completing the square. The link between real zeros and Hill coefficients at most 1 is known for binding polynomials in general.</p>""",
  fig=True, figcls="", cap="Zeros of 1 + 2x + ωx² as ω runs from 0.1 to 10: two real zeros (blue) meet at x = −1 when ω = 1, then leave the real line along the circle through 0 and −1 (red).", js=zeros_js),
 dict(id="mwc", nav="Many sites", h="Folding the dial over n sites gives the MWC model with exclusive binding",
  body=f"""<p>Folding the dial over n identical sites gives 1 + ((1 + kx)ⁿ − 1)/k. Expanded, a state with j ligands carries C(n, j)·k<sup>j−1</sup>: only the first ligand pays full price, and every later one is k times more welcome. Multiplied through by k this is (k − 1) + (1 + kx)ⁿ, which is the Monod–Wyman–Changeux model of a protein with a tense state that binds nothing and a relaxed state with n independent sites, with allosteric constant L = k − 1. The saturation curves match exactly (checked for n = 2, 3, 4), and Lean proves the identity.</p>
<p>So the weighted fold is not a new model of cooperativity: it is the exclusive-binding limit of MWC. For four sites its Hill coefficient climbs from 1 at k = 1 to {hill['100']:.2f} at k = 100 and {hill['1000']:.2f} at k = 1000, approaching 4 slowly.</p>""",
  fig=True, figcls="", cap="Hill coefficient of the dial folded over 2, 3 and 4 sites, against k. Each rises from 1 at the circle operation toward the number of sites.", js=hill_js),
 dict(id="haber", nav="Haber equilibrium", h="Pressure pushes nitrogen fixation forward, and non-ideality helps",
  body=f"""<p>Mass action multiplies partial pressures into the reaction quotient, so N₂ + 3H₂ ⇌ 2NH₃ goes as pressure to the minus two. At 700 K with a 1 : 3 feed, the constant-ΔH model from subject 7 gives an ammonia fraction of {100*hi['1']:.1f}% at 1 bar, {100*hi['100']:.0f}% at 100 bar and {100*hi['200']:.0f}% at 200 bar. The van der Waals fugacities from subject 6 raise the 200-bar value to {100*r['haber_200_vdw']:.0f}%. Holding ΔH constant overstates K at 700 K, so these are rough model values, not reactor data.</p>""",
  fig=True, figcls="", cap="Equilibrium ammonia fraction against pressure at 700 K in the constant-ΔH ideal-gas model, with the van der Waals correction at 200 bar.", js=haber_js),
]
claims = [
 ("Two-site binding polynomial: P − 1 = a ⊞_ω b", "tidy reading", "Python; Lean <code>binding_dial</code>"),
 ("Diacids: ω = 4·10^(−ΔpKa), from 0.004 (oxalic) to 0.75 (sebacic), always < 1", "known", "Wikipedia pKa table; Lean <code>diacid_gap</code>"),
 ("Two sites: n_H > 1 ⇔ ω > 1 ⇔ complex zeros on |x| = 1/√ω", "known", "Python; Lean <code>two_site_square</code>, <code>cooperative_no_real_zero</code>"),
 ("The n-site dial fold is MWC with exclusive binding and L = k − 1", "tidy reading", "saturation curves n = 2–4; Lean <code>fold_is_mwc</code>, <code>four_site_coeffs</code>"),
 ("Haber at 700 K: 1 bar 0.6%, 200 bar 41% (45% with vdW), constant-ΔH model", "known", "Python"),
]
lean = [("binding_dial","binding polynomial = dial"),("diacid_gap","β₂·4 = ω·β₁²"),("two_site_disc, two_site_square","discriminant 4(1 − ω)"),
 ("cooperative_no_real_zero","ω > 1 leaves no real zero"),("independent_double_zero","ω = 1 is (1 + x)²"),
 ("fold_lift, fold_is_mwc, four_site_coeffs","the n-site fold is MWC with L = k − 1"),("haber_pressure","Q scales as p⁻²"),("sq_nonneg_rat","helper")]
n = build("index.html", title="Equilibrium and mass action", tile=("8","K","Equilibrium"),
 h1="Binding polynomials run on the dial: diacids read it, cooperativity starts at the circle operation, and the n-site fold is the MWC model",
 lead="A two-site binding polynomial minus one is the dial with k = the interaction factor. Diacids measure it from 0.004 to 0.75 as the charges separate. For two sites cooperativity, zeros off the real line and k > 1 begin together, and folding the dial over n sites gives exactly the MWC model with exclusive binding, L = k − 1.",
 figs=[("12 / 12","checks pass"),("11","Lean theorems, compiled"),("L = k − 1","n-site fold = MWC"),("0.004 → 0.75","diacid dial, oxalic to sebacic")],
 picture=picture, caption="The interaction dial of nine dicarboxylic acids, from their two pKa values. It rises toward k = 1 (independent sites) as the two acid groups move apart, but stays below it because like charges repel.",
 verdict_h="Lead check: not a lead",
 verdict="Every result here is exact, and the MWC identity gives the weighted fold a clean biochemical meaning, but MWC, Bjerrum's diacid picture and the real-zeros criterion are all established. Nothing here predicts something those models don't, so I didn't follow it further.",
 sections=sections, claims=claims, lean_file="chemistry/08-equilibrium/lean/ChemEquilibrium.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/08-equilibrium/: <code>checks.py</code> (12 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemEquilibrium.lean</code>, <code>make_page.py</code>.",
 sources='Sources: <a href="https://en.wikipedia.org/wiki/Dicarboxylic_acid">Dicarboxylic acid (Wikipedia)</a> for pKa values; Bjerrum (1923); Monod, Wyman and Changeux (1965); Wyman and Gill, Binding and Linkage (1990); Lee and Yang (1952).',
 series_nav="Chemistry pass, part 8 of 13: equilibrium. Before: chemical thermodynamics. Next: rates and Arrhenius.")
print("page bytes", n)
