import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json")); pat = json.load(open("patterns.json"))
steps = list(range(1, 21))
picture = """
plot(document.getElementById('pic'),{x:[0,20.5],y:[0,1],ml:62,mr:%MR%,xticks:[1,5,10,15,20],yticks:[0,0.2,0.4,0.6,0.8,1],yfmt:v=>Math.round(v*100)+'%',xlab:'steps in a row',ylab:'overall yield = 1 − (losses folded at k = −1)',
 series:[{pts:%A%,c:'var(--green)',name:'95% per step',short:'95%'},{pts:%B%,c:'var(--blue)',name:'90% per step',short:'90%'},{pts:%C%,c:'var(--red)',name:'80% per step',short:'80%'},
 {pts:[[3,0.512]],c:'var(--red)',line:false,dots:true,r:7,labels:[{x:3,y:0.512,t:'8 fragments, convergent',dx:12,dy:-8}]},{pts:[[7,0.2097]],c:'var(--red)',line:false,dots:true,r:7,labels:[{x:7,y:0.2097,t:'8 fragments, linear',dx:12,dy:-10}]}]});
""".replace("%A%", json.dumps([[n, 0.95**n] for n in steps])).replace("%B%", json.dumps([[n, 0.9**n] for n in steps])).replace("%C%", json.dumps([[n, 0.8**n] for n in steps])).replace("%MR%", "130")

iso_js = """
{const el=document.getElementById('isotopes-fig');const P=%P%;const draw=()=>{const W=el.clientWidth,H=el.clientHeight,small=W<520;
const names=[['Cl2','Cl₂: 9 : 6 : 1'],['Br2','Br₂: 1 : 2 : 1'],['CHCl3','CHCl₃'],['C60','C₆₀']];const cw=W/4;let s='';
names.forEach(([k,t],i)=>{const p=P[k];const mx=Math.max(...p);const x0=cw*i+cw*0.12,bw=cw*0.76/p.length,base=H-46,top=34;
 p.forEach((v,j)=>{const h=(base-top)*v/mx;s+=`<rect x="${x0+j*bw+2}" y="${base-h}" width="${bw-4}" height="${h}" fill="var(--blue)"/>`;
   if(v/mx>0.02&&(bw>30||j%2==0))s+=`<text x="${x0+(j+0.5)*bw}" y="${base+16}" text-anchor="middle" style="font-size:11px">M+${j}</text>`});
 s+=`<text x="${cw*(i+0.5)}" y="20" text-anchor="middle" style="font-size:${small?12:14}px;font-weight:700">${t}</text>`});
el.innerHTML=s};const ro=new ResizeObserver(draw);ro.observe(el);draw();}
""".replace("%P%", json.dumps(pat))

med_js = """
{const el=document.getElementById('proportions-fig');const draw=()=>{const W=el.clientWidth,H=el.clientHeight;const lo=0.66,hi=0.80;const ml=30,mr=30;
const X=v=>ml+(v-lo)/(hi-lo)*(W-ml-mr);let s=`<line x1="${ml}" x2="${W-mr}" y1="${H/2}" y2="${H/2}" stroke="var(--mute)" stroke-width="2"/>`;
[0.68,0.70,0.72,0.74,0.76,0.78].forEach(v=>{s+=`<line x1="${X(v)}" x2="${X(v)}" y1="${H/2-5}" y2="${H/2+5}" stroke="var(--mute)"/><text x="${X(v)}" y="${H/2+24}" text-anchor="middle" style="font-size:12px">${Math.round(v*100)}%</text>`});
const pts=[[0.6994,'Fe₂O₃','Fe : O = 2 : 3','var(--red)',-1],[0.7773,'FeO','1 : 1','var(--red)',-1],[0.7236,'Fe₃O₄','3 : 4 = (1+2) : (1+3)','var(--blue)',1]];
pts.forEach(([v,n,t,c,up])=>{const y=H/2+(up<0?-34:52);s+=`<circle cx="${X(v)}" cy="${H/2}" r="8" fill="${c}"/><text x="${X(v)}" y="${y}" text-anchor="middle" style="font-weight:700">${n}</text><text x="${X(v)}" y="${y+17}" text-anchor="middle" style="font-size:12px">${t}</text>`});
s+=`<text x="${ml}" y="16" style="font-size:13px">iron, percent by mass</text>`;el.innerHTML=s};const ro=new ResizeObserver(draw);ro.observe(el);draw();}
"""
sections = [
 dict(id="yields", nav="Yields", h="Yields multiply, so losses the circle operation at k = −1",
  body=f"""<p>A step that keeps a fraction y of the material loses ℓ = 1 − y. Two steps keep y₁y₂, so their losses combine by the weighted circle operation (called mult-add before 1.0) at k = −1, the same failure-union law as two independent risks:</p>
<div class="eq">ℓ₁ ⊞₋₁ ℓ₂ = ℓ₁ + ℓ₂ − ℓ₁ℓ₂ = 1 − (1 − ℓ₁)(1 − ℓ₂)</div>
<p>The fold is associative and order-free, so n equal steps lose the k = −1 power 1 − (1 − ℓ)ⁿ. Ten steps at 90% each keep only 34.9%.</p>
<p>What the law says about planning is that only the number of steps a molecule passes through matters, not the total number of steps. Joining eight fragments one at a time puts the first fragment through seven couplings, 0.8⁷ = {100*r['linear']:.0f}% at 80% each. Joining them in pairs, then pairs of pairs, puts every fragment through three, 0.8³ = {100*r['convergent']:.0f}%. That is the standard case for convergent synthesis, and it is just the depth of the fold.</p>"""),
 dict(id="isotopes", nav="Isotope patterns", h="A molecule's isotope pattern is the product of its atoms' patterns",
  body=f"""<p>Write each atom's isotopes as a polynomial in a mass-step x: chlorine is 0.758 + 0.242x². A molecule's pattern is the product of its atoms' polynomials, so patterns compose by the lift. Two chlorines give (3 + x)² ≈ 9 : 6 : 1 and two bromines (1 + x)² ≈ 1 : 2 : 1, the peaks every mass spectrometrist knows.</p>
<p>The k = −1 law shows up again in the chance that a molecule carries at least one heavy atom, 1 − (1 − p)ⁿ. For C₆₀ that is {100*r['C60']['any13']:.1f}%, and the M+1 peak is {100*r['C60']['M1']:.1f}% of the M peak.</p>""",
  fig=True, figcls="", cap="Isotope patterns from multiplying atom polynomials with IUPAC abundances, each scaled to its tallest peak.", js=iso_js),
 dict(id="proportions", nav="Proportions", h="Whole-number ratios, and mixed oxides are mediants",
  body="""<p>Dalton's law of multiple proportions holds to the precision of the atomic masses: the grams of oxygen per gram of nitrogen in N₂O, NO, N₂O₃, NO₂ and N₂O₅ stand exactly 1 : 2 : 3 : 4 : 5.</p>
<p>The definite proportions of a mixed compound follow the mediant rule, because combining two compounds adds their atom counts. Magnetite Fe₃O₄ is FeO plus Fe₂O₃, so its iron-to-oxygen ratio 3 : 4 is (1 + 2) : (1 + 3), which lies between 1 : 1 and 2 : 3. The same holds for red lead Pb₃O₄ = 2PbO + PbO₂ and hausmannite Mn₃O₄. Lean proves the mediant always lies between its parents. These are counting facts, not compose by the circle operation.</p>""",
  fig=True, figcls="", cap="Magnetite's iron content sits between its two parent oxides, as the mediant of their Fe : O ratios must.", js=med_js),
 dict(id="balancing", nav="Balancing", h="Balancing an equation means finding the one direction the element counts leave free",
  body="""<p>Write each species' element counts as a column, reactants positive and products negative. A balanced equation is a vector that sums the columns to zero. For the Haber process, respiration and photosynthesis there is exactly one such direction, giving N₂ + 3H₂ → 2NH₃ and C₆H₁₂O₆ + 6O₂ ⇌ 6CO₂ + 6H₂O.</p>
<p>Nitrogenase is the instructive case. N₂ + H → NH₃ + H₂ leaves a two-way family of balanced equations, so element counts alone cannot say how much hydrogen the enzyme wastes. The measured one H₂ per N₂ is a fact about its mechanism, and it comes back in subject 13.</p>"""),
]
claims = [
 ("Losses of steps in a row combine as circK(−1); n steps lose 1 − (1 − ℓ)ⁿ", "tidy reading", "Python; Lean <code>loss_two</code>, <code>loss_n</code>, <code>yield_lift</code>"),
 ("Yield depends on the depth of the coupling tree: 21% linear vs 51% convergent for 8 fragments at 80%", "known", "Python; Lean <code>eight_fragments</code>"),
 ("Isotope patterns are products of atom polynomials (the lift); Cl₂ 9 : 6 : 1, Br₂ 1 : 2 : 1", "known", "convolution with IUPAC abundances; Lean <code>chlorine_961</code>, <code>bromine_121</code>"),
 ("Chance of at least one heavy atom is the k = −1 power", "tidy reading", "C₆, C₂₀, C₆₀; Lean <code>any_heavy</code>"),
 ("Nitrogen oxides 1 : 2 : 3 : 4 : 5; mixed oxides are mediants", "known", "atomic masses; Lean <code>multiple_proportions</code>, <code>mediant_between</code>"),
 ("Balancing = the null direction of the element matrix; nitrogenase's H₂ is not fixed by stoichiometry", "known", "sympy null spaces; Lean <code>nitrogenase_family</code>"),
]
lean = [("loss_two, loss_n, yield_lift, loss_assoc","losses fold at k = −1, yields multiply"),("ten_steps","ten steps at 90%"),("eight_fragments","tree depth 3 vs chain depth 7"),
 ("two_atoms, chlorine_961, bromine_121","isotope polynomials multiply"),("any_heavy","at least one heavy atom"),("mediant_between, magnetite, magnetite_between","mixed compounds are mediants"),
 ("multiple_proportions","N₂O … N₂O₅"),("haber_balanced, respiration_balanced, nitrogenase_family","balanced equations")]
n = build("index.html", title="Stoichiometry and proportions", tile=("5","NH₃","Stoichiometry"),
 h1="Yields multiply, so losses fold at k = −1; isotope patterns multiply too; proportions are counting",
 lead="Steps in a row keep the product of their yields, so their losses combine by the weighted circle operation at k = −1 and only a molecule's depth in the synthesis tree matters. Isotope patterns are products of atom polynomials. Multiple and definite proportions are whole-number counting, with mixed oxides as mediants.",
 figs=[("14 / 14","checks pass"),("17","Lean theorems, compiled"),("51% vs 21%","convergent vs linear, 8 fragments"),("1 : 2 : 3 : 4 : 5","nitrogen oxides")],
 picture=picture, caption="Overall yield after n steps at a fixed yield per step. The red dots compare joining eight fragments one at a time (seven couplings deep) with joining them in pairs (three deep).",
 verdict_h="Lead check: not a lead",
 verdict="The k = −1 law is exactly right for yields and for heavy-isotope odds, but chemists already work with products of yields and with isotope generating functions, so writing them as the circle operation changes no calculation. Proportions and balancing are counting and linear algebra with no the circle operation in them.",
 sections=sections, claims=claims, lean_file="chemistry/05-stoichiometry/lean/ChemStoichiometry.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/05-stoichiometry/: <code>checks.py</code> (14 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>patterns.json</code>, <code>lean/ChemStoichiometry.lean</code>, <code>make_page.py</code>.",
 sources='Sources: IUPAC representative isotopic abundances; Dalton (1808) on multiple proportions; Velluz, Valls and Mathieu (1967) and Hendrickson (1977) on convergent synthesis; Seefeldt et al. on nitrogenase stoichiometry.',
 series_nav="Chemistry pass, part 5 of 13: stoichiometry. Before: bonding and valence. Next: gas laws.")
print("page bytes", n)
