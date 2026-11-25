import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
lv = [[N, 2*math.cos(j*math.pi/(N+1))] for N in range(1, 15) for j in range(1, N+1)]
picture = """
plot(document.getElementById('pic'),{x:[0,15],y:[-3.2,3.2],ml:44,mr:%MR%,xticks:[1,2,4,6,8,10,12,14],yticks:[-3,-2,-1,0,1,2,3],xlab:'carbons in the chain, N',ylab:'π level, units of |β| (bonding up)',
 hlines:[{y:2,c:'var(--yel)',label:'border: band edge x = 2',short:'edge'},{y:-2,c:'var(--yel)',label:'border: x = −2',short:'edge'}],
 series:[{pts:%LV%,c:'var(--blue)',line:false,dots:true,r:4}],
 extra:({X,Y,W,H,m,small})=>`<rect x="${m.l}" y="${Y(3.2)}" width="${W-m.l-m.r}" height="${Y(2)-Y(3.2)}" fill="var(--red)" opacity=".10"/><rect x="${m.l}" y="${Y(-2)}" width="${W-m.l-m.r}" height="${Y(-3.2)-Y(-2)}" fill="var(--red)" opacity=".10"/>`+
  (small?'':`<text x="${W-m.r+6}" y="${Y(2.7)}">two walls: no levels</text><text x="${W-m.r+6}" y="${Y(0)+5}">turning: the band</text><text x="${W-m.r+6}" y="${Y(-2.7)}">two walls</text>`)});
""".replace("%LV%", json.dumps(lv)).replace("%MR%", "190")

xs = [i/50 for i in range(-200, 201)]
def rate(x):
    return math.acosh(abs(x)/2) if abs(x) > 2 else 0.0
def rate_alt(x, t1=1.2, t2=0.8):
    tr = (x*x - t1*t1 - t2*t2)/(t1*t2)
    return math.acosh(abs(tr)/2)/2 if abs(tr) > 2 else 0.0
bands_js = """
plot(document.getElementById('bands-fig'),{x:[-4,4],y:[0,1.4],ml:48,xticks:[-4,-3,-2,-1,0,1,2,3,4],xlab:'energy x, units of β',ylab:'decay per carbon, arccosh(|x|/2)',
 bands:[{x0:-4,x1:-2,fill:'var(--red)',label:'two walls'},{x0:-2,x1:2,fill:'var(--blue)',label:'turning: band'},{x0:2,x1:4,fill:'var(--red)',label:'two walls'}],
 mr:16,series:[{pts:%A%,c:'var(--ink)'},{pts:%B%,c:'var(--red)'}],extra:({X,Y})=>`<text x="${X(0)}" y="${Y(0.75)}" text-anchor="middle">black: even bonds</text><text x="${X(0)}" y="${Y(0.75)+18}" text-anchor="middle" fill="var(--red)">red: alternating 1.2 / 0.8</text><text x="${X(0)}" y="${Y(0.32)}" text-anchor="middle" fill="var(--red)">gap</text>`});
""".replace("%A%", json.dumps([[x, rate(x)] for x in xs])).replace("%B%", json.dumps([[x, rate_alt(x)] for x in xs]))

rings_js = """
{const el=document.getElementById('rings-fig');const draw=()=>{const W=el.clientWidth,H=el.clientHeight,small=W<560;
const rings=[{N:6,mob:false,ne:6,t:'benzene, 6 e',t2:'closed shell'},{N:4,mob:false,ne:4,t:'cyclobutadiene, 4 e',t2:'two unpaired'},{N:4,mob:true,ne:4,t:'Möbius 4-ring, 4 e',t2:'closed shell'}];
const cw=W/3,R=Math.min(cw*0.30,H*0.32);let s='';
rings.forEach((g,i)=>{const cx=cw*(i+0.5),cy=H*0.45;
 s+=`<circle cx="${cx}" cy="${cy}" r="${R}" fill="none" stroke="var(--mute)" stroke-dasharray="4 4"/>`;
 s+=`<line x1="${cx-R-20}" x2="${cx+R+20}" y1="${cy}" y2="${cy}" stroke="var(--yel)" stroke-width="2"/>`;
 const pts=[];for(let j=0;j<g.N;j++){const th=Math.PI*(2*j+(g.mob?1:0))/g.N;pts.push([cx+R*Math.sin(th),cy+R*Math.cos(th)])}
 s+=`<polygon points="${pts.map(p=>p.join(',')).join(' ')}" fill="none" stroke="var(--blue)" stroke-width="2"/>`;
 const ord=pts.map((p,j)=>j).sort((a,b)=>pts[b][1]-pts[a][1]);
 const occ=new Array(g.N).fill(0);let left=g.ne;
 // fill lowest first; degenerate pairs get one each before pairing (Hund)
 let k=0;while(left>0&&k<ord.length){const y=pts[ord[k]][1];const grp=ord.filter(j=>Math.abs(pts[j][1]-y)<1e-6);
  grp.forEach(j=>{if(left>0){occ[j]++;left--}});const pair=left>=grp.length;grp.forEach(j=>{if(pair&&left>0){occ[j]++;left--}});k+=grp.length}
 pts.forEach((p,j)=>{s+=`<line x1="${p[0]-16}" x2="${p[0]+16}" y1="${p[1]}" y2="${p[1]}" stroke="var(--ink)" stroke-width="4"/>`;
  for(let e=0;e<occ[j];e++)s+=`<text x="${p[0]-12+e*13}" y="${p[1]-6}" style="font-size:15px;font-weight:700;fill:var(--red);stroke:none">${e?'↓':'↑'}</text>`});
 s+=`<text x="${cx}" y="${H-30}" text-anchor="middle" style="font-size:${small?12:14}px">${g.t}</text><text x="${cx}" y="${H-12}" text-anchor="middle" style="font-size:${small?12:14}px">${g.t2}</text>`;});
if(!small)s+=`<text x="8" y="18" style="font-size:13px">yellow line: zero (nonbonding); levels below it bond</text>`;
el.innerHTML=s};const ro=new ResizeObserver(draw);ro.observe(el);draw();}
"""

lad = r["ladders"]
cols = {"C–C": "var(--blue)", "N–N": "var(--red)", "C–N": "var(--green)", "C–O": "var(--ink)"}
ser = []
for name, f in lad.items():
    e1, i2, i3 = f["inc"]; meas = [[1, e1], [2, e1+i2], [3, e1+i2+i3]]
    ser.append(dict(pts=meas, c=cols[name], dots=True))
    k = f["k"]; pred = [[n/10, ((1+k*e1)**(n/10)-1)/k] for n in range(10, 31)]
    ser.append(dict(pts=pred, c=cols[name], dash="6 5", w=2))
leg = json.dumps(list(cols.items()))
ladder_js = "plot(document.getElementById('ladder-fig'),{x:[0.8,3.2],y:[0,1400],ml:56,mr:16,xticks:[1,2,3],xfmt:v=>['','single','double','triple'][v],ylab:'bond enthalpy, kJ/mol',series:%s,extra:({X,Y,m})=>%s.map((c,i)=>`<circle cx=\"${m.l+18}\" cy=\"${m.t+14+i*20}\" r=\"5\" fill=\"${c[1]}\"/><text x=\"${m.l+30}\" y=\"${m.t+19+i*20}\">${c[0]}</text>`).join('')+`<text x=\"${m.l+18}\" y=\"${m.t+100}\">dashed: one dial</text>`});" % (json.dumps(ser), leg)

sections = [
 dict(id="overlap", nav="Overlap", h="Overlap is a cosine, and the bond follows the same law as a bound mass",
  body=f"""<p>Mix two atomic orbitals as a φ<sub>a</sub> + b φ<sub>b</sub>. Their overlap S works as the cosine between them, and the size of the mixture is</p>
<div class="eq">|a φ<sub>a</sub> + b φ<sub>b</sub>|² = a² + b² + 2ab·S</div>
<p>That is the same law of cosines as the invariant mass on the previous page, M² = m₁² + m₂² + 2m₁m₂c. S = 1 is the border, where the two orbitals are one and their weights simply add. |S| &lt; 1 is the turning country. For hydrogen's 1s orbitals at the H₂ bond length, S = {r['S_H2']:.3f}.</p>
<p>The cosine also explains a textbook asymmetry. Solving H c = E S c gives a bonding level (α + β)/(1 + S) and an antibonding level (α − β)/(1 − S). The antibonding one is pushed up more than the bonding one comes down, by exactly (1 + S)/(1 − S): about {r['split_ratio']:.0f} times for H₂ in the extended Hückel estimate. That is why filling both (as in He₂) costs energy overall.</p>"""),
 dict(id="bands", nav="Chains and bands", h="A carbon chain is the landscape point (x, 1): band, edge and gap are the three countries",
  body="""<p>In Hückel theory a chain of carbons obeys one recurrence: the polynomial whose roots are the levels of an N-carbon chain satisfies P(N + 2) = x·P(N + 1) − P(N), with x the energy in units of β. That is the landscape recurrence with k = x and q = 1. Lean proves its invariant P(N+1)² − P(N+2)P(N) equals q<sup>N+1</sup> in general, so a plain chain keeps it at exactly 1 at every length.</p>
<p>The three countries are the three things an electron can do along the chain:</p>
<ul><li><b>Turning, |x| &lt; 2:</b> the step multiplier is a rotation, so the wave travels. Every level of every chain lands here, at 2cos(jπ/(N+1)) (the picture at the top).</li>
<li><b>Border, |x| = 2:</b> the band edges, where the root is double.</li>
<li><b>Two walls, |x| &gt; 2:</b> a real pair λ, 1/λ, so the wave dies off by arccosh(|x|/2) per carbon. These are gaps.</li></ul>
<p>Alternate the bonds (long, short, long, as in polyacetylene) and the two-step transfer still has q = 1, but a gap of width 2|t₁ − t₂| opens in the middle of the band, inside the two-walls country. A 200-carbon chain shows it at {:.3f} against the predicted 0.400.</p>""".format(r["poly_gap_finite"]),
  fig=True, figcls="", cap="How fast a state dies off per carbon, against its energy. Inside the band (blue) it doesn't die at all. Bond alternation (red) opens a two-walls gap in the middle, the semiconductor gap of polyacetylene.", js=bands_js),
 dict(id="rings", nav="Rings", h="Rings: 4n + 2 and 4n are the two ways to close the chain",
  body=f"""<p>Close the chain into a ring and the levels become 2cos(2πj/N), the vertices of a polygon in a circle (Frost's picture). A plain ring has one bottom level and then pairs, so it closes its bonding shell at 2 + 4m electrons: Hückel's rule. Give the ring one twist (a Möbius ring) and every level comes in a pair, so it closes at 4m. Checked for every ring of 3 to 19 atoms, and Lean proves the two counts.</p>
<p>Benzene's six π electrons give {r['benzene']:.0f}β, which is 2β more than three separate double bonds: its aromatic stabilization. Cyclobutadiene's four give exactly {r['cyclobutadiene']:.0f}β, nothing extra, with two electrons left unpaired on the middle line. Butadiene, the open chain, gains {r['butadiene']-4:.3f}β.</p>""",
  fig=True, figcls="", cap="Frost diagrams: each ring's levels sit at the corners of a polygon in a circle, vertex down for a plain ring and edge down for a Möbius ring. Levels below the yellow zero line are bonding; arrows are electrons.", js=rings_js),
 dict(id="orders", nav="Bond orders", h="Bond order counts bonding minus antibonding pairs, and N₂ has three",
  body="""<p>Fill the second-row diatomics' valence orbitals and halve the difference between bonding and antibonding electrons: Li₂ 1, Be₂ 0, B₂ 1, C₂ 2, N₂ 3, O₂ 2, F₂ 1, Ne₂ 0. Those match the measured bonds, including O₂'s two unpaired electrons. This is plain counting, with no circle operation (called mult-add before 1.0) in it.</p>"""),
 dict(id="ladder", nav="Bond ladders and N₂", h="One dial does not describe single, double and triple bonds, and nitrogen shows why",
  body=f"""<p>Bond energies don't add: a double bond is not two single bonds. The natural the circle operation test is to treat each extra bond as a weighted circle operation with one dial k, fixed by the single and double bonds, and see if it predicts the triple. Using the LibreTexts table of average bond enthalpies:</p>
<ul><li><b>C–C:</b> predicts {lad['C–C']['pred']:.0f} kJ/mol against {lad['C–C']['meas']} ({100*lad['C–C']['err']:+.0f}%).</li>
<li><b>C–N:</b> predicts {lad['C–N']['pred']:.0f} against {lad['C–N']['meas']} ({100*lad['C–N']['err']:+.0f}%).</li>
<li><b>N–N:</b> predicts {lad['N–N']['pred']:.0f} against {lad['N–N']['meas']} ({100*lad['N–N']['err']:+.0f}%).</li>
<li><b>C–O:</b> predicts {lad['C–O']['pred']:.0f} against {lad['C–O']['meas']} ({100*lad['C–O']['err']:+.0f}%).</li></ul>
<p>So the dial fits carbon's ladder and fails the rest. Nitrogen is the instructive failure. Carbon's increments shrink (348, 266, 225 kJ/mol), but nitrogen's grow (163, 255, 523): its single bond is weak (lone pairs on neighbouring atoms repel) and its π bonds are strong. Breaking N≡N down to N=N costs 523 of the 941 kJ/mol, {100*r['N_first_step_share']:.0f}%. That first step is why nitrogen fixation is hard, for the Haber process and for nitrogenase alike. The thermodynamics and the rate of that step come back in subjects 7 and 9.</p>""",
  fig=True, figcls="", cap="Average bond enthalpies for single, double and triple bonds (dots), with the dial curve fixed by the first two (dashed). Carbon follows its curve; nitrogen's triple bond sits far above it.", js=ladder_js),
]
claims = [
 ("|aφa + bφb|² = a² + b² + 2ab·S, the invariant-mass law with S as cosine", "tidy reading", "Python; Lean <code>overlap_border</code>, <code>overlap_opposite</code>, <code>overlap_orthogonal</code>"),
 ("Antibonding rises more than bonding falls, by (1+S)/(1−S)", "known", "generalized eigenproblem; Lean <code>split_asymmetry</code>"),
 ("The Hückel chain recurrence is the landscape point (x, 1); its Cassini invariant is q^(N+1)", "tidy reading", "Python; Lean <code>cassini</code>, <code>huckel_unimodular</code>"),
 ("Band = turning, edge = border, gap = two walls (decay arccosh(|x|/2) per site)", "tidy reading", "Python; Lean <code>no_root_in_band</code>, <code>band_edge</code>, <code>gap_pair</code>"),
 ("Bond alternation opens a gap 2|t₁ − t₂|", "known", "Python, infinite and 200-site chains"),
 ("Hückel 4n+2 and Möbius 4n closed shells; benzene +2β, cyclobutadiene 0", "known", "rings of 3–19 atoms; Lean <code>huckel_rule</code>, <code>mobius_rule</code>"),
 ("MO bond orders Li₂ … Ne₂", "known", "Python; Lean <code>n2_triple</code>, <code>o2_double</code>, <code>f2_single</code>"),
 ("One dial fits the C–C bond ladder to 3% but misses N–N, C–N and C–O by 9–25%", "new here", "checks.py with LibreTexts bond enthalpies"),
]
lean = [("overlap_border, overlap_opposite, overlap_orthogonal, bonding_norm","the overlap law at S = 1, −1, 0, and 2 + 2S"),
 ("split_asymmetry","both splitting numerators are αS − β"),("chainP, cassini_step, cassini","the landscape recurrence and its invariant q^(n+1)"),
 ("huckel_unimodular","plain chains keep the invariant at 1"),("allyl, butadiene","the 3- and 4-carbon level polynomials"),
 ("band_square, no_root_in_band","inside the band there is no real step multiplier"),("band_edge, gap_pair","the double root at the edge; the reciprocal pair in the gap"),
 ("huckel_rule, mobius_rule, benzene_aromatic","4m + 2 and 4m"),("n2_triple, o2_double, f2_single","bond orders"),("dial_ladder","three steps of one dial, written out")]
n = build("index.html", title="Bonding and valence", tile=("4","N≡N","Bonding"),
 h1="A bond's overlap is a cosine, a carbon chain is the landscape point (x, 1), and nitrogen's triple bond breaks the dial",
 lead="Mixing orbitals follows the same law of cosines as a bound mass. Chains of carbons follow the landscape recurrence with q = 1, so the band, its edges and its gaps are the turning, border and two-walls countries. A single dial fits carbon's single-double-triple ladder but not nitrogen's, whose hard first step is the problem in nitrogen fixation.",
 figs=[("19 / 19","checks pass"),("23","Lean theorems, compiled"),("q = 1","every Hückel chain"),("56%","of N≡N's energy in its first step")],
 picture=picture, caption="Every π level of every Hückel chain from 1 to 14 carbons. They all fall in the turning band between the borders at ±2|β|; nothing lands in the two-walls region outside.",
 verdict_h="Lead check: not a lead",
 verdict="The overlap-as-cosine and band-as-country readings are exact and tie this subject to the conservation page, but they are the standard transfer-matrix picture of tight binding renamed. The bond-ladder test is a genuine check of the dial and it fails outside carbon, so it gives nothing to follow.",
 sections=sections, claims=claims, lean_file="chemistry/04-bonding/lean/ChemBonding.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/04-bonding/: <code>checks.py</code> (19 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemBonding.lean</code>, <code>make_page.py</code>.",
 sources='Sources: <a href="https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Book:_Introduction_to_Inorganic_Chemistry_(Wikibook)/12:_Resources_for_Students_and_Teachers/12.07:_Appendix_3-_Bond_Enthalpies">LibreTexts, Appendix 3: Bond Enthalpies</a>; Hückel (1931); Frost and Musulin (1953); Heilbronner (1964) on Möbius rings; Su, Schrieffer and Heeger (1979) on polyacetylene; Wolfsberg and Helmholz (1952) for β ≈ 1.75 S α.',
 series_nav="Chemistry pass, part 4 of 13: bonding and valence. Before: the periodic law. Next: stoichiometry.")
print("page bytes", n)
