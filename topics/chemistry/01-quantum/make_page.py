import sys, json, numpy as np
from math import comb
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
D, g = 1.0, 0.3
e1 = D - np.sqrt(D*D + g*g)
def ansatz(N, k):
    best = 0.0
    for c in np.linspace(-0.4, 0.05, 2251):
        a = [1.0] + [k**(m-1) * c**m for m in range(1, N+1)]
        num = den = 0.0
        for m in range(N+1):
            den += comb(N, m)*a[m]**2; num += comb(N, m)*a[m]**2*2*D*m
            if m < N: num += 2*comb(N, m+1)*(m+1)*g*a[m]*a[m+1]
        best = min(best, num/den)
    return best/(N*e1)
Ns = list(range(1, 17))
k0 = [[N, (D-np.sqrt(D*D+N*g*g))/(N*e1)] for N in Ns]
kh = [[N, ansatz(N, 0.5)] for N in Ns]
k1 = [[N, 1.0] for N in Ns]
ph = r["pyscf_h2"]; pw = r["pyscf_h2o"]
picture = f"""
plot(document.getElementById('pic'),{{mr:215,x:[1,16],y:[0.7,1.02],xlab:'number of identical molecules, far apart',ylab:'share of the true correlation energy kept',
 yfmt:v=>Math.round(v*100)+'%',
 series:[{{pts:{json.dumps(k1)},c:'var(--green)',name:'k = 1: coupled cluster (exact)',short:'k = 1',dy:-6}},
 {{pts:{json.dumps(kh)},c:'var(--blue)',name:'k = ½: half the corner',short:'k = ½'}},
 {{pts:{json.dumps(k0)},c:'var(--red)',name:'k = 0: CISD, no corner',short:'k = 0'}}],
 extra:({{X,Y,small}})=>small?'':`<text x="${{X(9)}}" y="${{Y(0.74)}}">the loss is the corner: pairs of molecules excited together</text>`}});
"""
cube_js = r"""
{const el=document.getElementById('corner-fig');const draw=()=>{const W=el.clientWidth,H=el.clientHeight,s=Math.min(W*0.42,H*0.55),ox=W*0.18,oy=H*0.82,dx=s*0.45,dy=-s*0.32;
const P=(a,b,c)=>[ox+a*s+c*dx,oy-b*s+c*dy];const V=[];for(let a=0;a<2;a++)for(let b=0;b<2;b++)for(let c=0;c<2;c++)V.push([a,b,c]);
let t='';const E=[];V.forEach(p=>V.forEach(q=>{const d=p.reduce((s,x,i)=>s+Math.abs(x-q[i]),0);if(d===1&&p<q)E.push([p,q])}));
E.forEach(([p,q])=>{const A=P(...p),B=P(...q);const kept=(p.join('')==='000');t+=`<line x1="${A[0]}" y1="${A[1]}" x2="${B[0]}" y2="${B[1]}" stroke="${kept?'var(--red)':'var(--rule)'}" stroke-width="${kept?4:2}"/>`});
const name=p=>{const n=['A','B','C'].filter((x,i)=>p[i]);return n.length?n.join('·'):'reference'};
V.forEach(p=>{const A=P(...p),o=p.reduce((a,b)=>a+b,0);const col=o===0?'var(--ink)':o===1?'var(--red)':'var(--blue)';t+=`<circle cx="${A[0]}" cy="${A[1]}" r="8" fill="${col}"/><text x="${A[0]+12}" y="${o===0?A[1]+26:A[1]-10}">${name(p)}</text>`});
el.innerHTML=t};const ro2=new ResizeObserver(draw);ro2.observe(el);draw();}
"""
sections = [
 dict(id="corner", nav="The corner", h="Two molecules that ignore each other: their wavefunctions multiply",
  body=f"""<p>Write a molecule's wavefunction as the reference determinant times 1 + c, where c holds every excited configuration (intermediate normalisation). For two molecules too far apart to interact, the exact wavefunction is the product, so</p>
<div class="eq">1 + c<sub>AB</sub> = (1 + c<sub>A</sub>)(1 + c<sub>B</sub>) &nbsp; ⟹ &nbsp; c<sub>AB</sub> = c<sub>A</sub> ✳ c<sub>B</sub> = c<sub>A</sub> + c<sub>B</sub> + c<sub>A</sub>c<sub>B</sub></div>
<p>The corner c<sub>A</sub>c<sub>B</sub> is "both molecules excited at the same time". For two two-electron molecules that is a quadruple excitation of the pair.</p>
<p>Truncated configuration interaction with singles and doubles (CISD) keeps only excitations up to two electrons, so it keeps c<sub>A</sub> + c<sub>B</sub> and throws the corner away. It is the weighted law at k = 0. That is the whole reason CISD is not size-consistent. Coupled cluster writes 1 + c = e<sup>T</sup>: T is the log of the lift, so T adds and the corner comes back by itself. In circle-operation terms, coupled cluster is the ruler that turns ✳ into +.</p>""",
  fig=True, figcls="", cap="Three fragments A, B, C on the circle-operation cube. CISD keeps the red edges out of the reference corner (one fragment excited at a time). Every blue corner is a product of fragment excitations that only the full ✳ keeps.", js=cube_js),
 dict(id="model", nav="N molecules", h="With N molecules, CISD keeps a share that falls like 1/√N",
  body=f"""<p>The checks file solves a model exactly: N identical molecules, each with a reference state and one doubly excited state 2Δ above it, coupled by g (Δ = 1, g = 0.3). The exact energy is N times one molecule. CISD's energy is Δ − √(Δ² + N g²), so its share of the correlation falls like 1/√N: 98% for two molecules, 86% for ten, under 7% for ten thousand.</p>
<p>A wavefunction that keeps the corner with weight k (every m-fold product scaled by k<sup>m−1</sup>) is exact only at k = 1. Half a corner is better than none but still loses ground as N grows. Size-consistency is not a matter of degree: only the full circle operation composes correctly.</p>
<p>Lean proves the model's key step: if e₁ solves one molecule, then N·e₁ misses the N-molecule CISD equation by exactly N(N − 1)e₁², the count of molecule pairs times the corner (<code>cisd_gap</code>, <code>cisd_not_extensive</code>).</p>"""),
 dict(id="real", nav="Real molecules", h="Real molecules show the same loss: 0.7 kcal/mol for two H₂, 5.9 for two waters",
  body=f"""<p>I ran PySCF on real molecules 100 Å apart, where they cannot interact.</p>
<ul><li><b>Two H₂ (cc-pVDZ):</b> CCSD gives exactly twice one molecule ({2*ph['cc1']:.6f} Eh). CISD gives {ph['ci2']:.6f} Eh against {2*ph['ci1']:.6f}: it loses {ph['missing_kcal']:.2f} kcal/mol, though for a single H₂ CISD and CCSD agree exactly (two electrons leave no room for a corner).</li>
<li><b>Two waters (6-31G):</b> CCSD is exactly twice one water; CISD keeps {100*pw['ci2']/(2*pw['ci1']):.1f}% and loses {abs(2*pw['ci1']-pw['ci2'])*627.5:.1f} kcal/mol, more than a hydrogen bond's worth of error from two molecules that are not even touching.</li></ul>
<p>The same counting runs through spin coupling: multiplicities multiply, (2s₁ + 1)(2s₂ + 1), so 2s combines by the circle operation, and the Clebsch–Gordan series fills the product exactly (36 cases checked, and in Lean).</p>"""),
]
claims = [
 ("Independent fragments: c<sub>AB</sub> = c<sub>A</sub> ✳ c<sub>B</sub>", "tidy reading", "200 random pairs; Lean <code>fragments_circ</code>"),
 ("CISD drops exactly the corner, so it is the k = 0 law", "tidy reading", "Lean <code>cisd_drops_corner</code>, <code>cisd_drops_three</code>"),
 ("Coupled cluster's exponential is the log of the lift, so it is size-consistent", "known", "PySCF: 2 H₂ and 2 H₂O, CCSD = 2× to 1e-6 Eh"),
 ("N-molecule CISD keeps a share falling like 1/√N", "known", "exact diagonalisation N ≤ 10; closed form to N = 160,000"),
 ("N·e₁ misses CISD by N(N−1)e₁², the molecule pairs", "tidy reading", "Lean <code>cisd_gap</code>"),
 ("Only corner weight k = 1 is size-consistent; k = ½ still drifts", "new here", "variational scan, N = 1 to 16"),
 ("Real CISD loses 0.71 kcal/mol (2 H₂) and 5.9 kcal/mol (2 H₂O)", "known", "PySCF runs in checks.py"),
 ("Spin multiplicities multiply, 2s combines by the circle operation", "tidy reading", "36 cases; Lean <code>spin_multiplicity</code>"),
]
lean = [("fragments_circ","(1 + c<sub>A</sub>)(1 + c<sub>B</sub>) = 1 + c<sub>A</sub> ✳ c<sub>B</sub>"),("cisd_drops_corner","what CISD drops is c<sub>A</sub>c<sub>B</sub>"),
 ("cisd_is_weight_zero, exact_is_weight_one","the two ends of the dial"),("n_fragments","1 + c<sub>N</sub> = (1 + c)<sup>N</sup>"),("fragments_assoc","joining N then M is joining N + M"),
 ("three_fragments, cisd_drops_three","the three-fragment corner expansion and what CISD loses"),("cisd_gap","the N(N − 1)e₁² gap"),("cisd_not_extensive","CISD is not size-extensive for N ≥ 2"),
 ("spin_multiplicity","(s₁+1)(s₂+1) = s₁ ✳ s₂ + 1"),("clebsch_gordan_square","the odd multiplicities of j ⊗ j fill a square")]
n = build("index.html", title="Quantum chemistry and the circle operation", tile=("1","Ψ","Quantum chemistry"),
 h1="Size-consistency is the circle operation: CISD fails because it drops the corner",
 lead="Two molecules that don't interact have a wavefunction equal to the product of theirs, so their excited parts compose by the circle operation. Truncated CI keeps only the edges and loses the corner; coupled cluster keeps it by working in the log of the lift.",
 figs=[("25 / 25","numeric checks pass"),("12","Lean theorems, compiled"),("5.9 kcal/mol","CISD's error on two far-apart waters"),("k = 1","the only size-consistent corner weight")],
 picture=picture, caption="Share of the true correlation energy kept as identical molecules are added, for corner weight k = 0 (CISD), ½ and 1 (coupled cluster). Model: Δ = 1, g = 0.3, solved exactly in checks.py.",
 verdict_h="Lead check: not a lead",
 verdict="Size-consistency is fully understood in quantum chemistry: it is why coupled cluster replaced truncated CI decades ago. The circle operation explains it in one line, which makes this the clearest teaching picture of the failure, but it adds no prediction. The one small new fact, that a partial corner weight still fails, follows directly from the algebra.",
 sections=sections, claims=claims, lean_file="chemistry/01-quantum/lean/ChemQuantum.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/01-quantum/: <code>checks.py</code> (25 checks, numpy and PySCF), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemQuantum.lean</code>, <code>make_page.py</code>.",
 sources="Sources: Szabo and Ostlund, <i>Modern Quantum Chemistry</i> (1982), ch. 4 on size-consistency of truncated CI; Bartlett and Musiał, Rev. Mod. Phys. 79, 291 (2007), on coupled cluster; PySCF (Sun et al., WIREs Comput. Mol. Sci. 2018).",
 series_nav="Chemistry pass, part 1 of 13: quantum chemistry. Next: conservation of mass and charge.")
print("page bytes", n)
