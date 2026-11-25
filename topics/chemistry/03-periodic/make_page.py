import sys, json
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
touch = {}
for Z, d in r["moves"].items():
    for k in d: touch[k] = touch.get(k, 0) + 1
picture = """
{const el=document.getElementById('pic');const L='spdf';const order=%s;const touch=%s;
const draw=()=>{const W=el.clientWidth,H=el.clientHeight,small=W<560;const cols=4,rows=7;
const ml=small?34:60,mr=small?10:230,mt=30,mb=26;const cw=(W-ml-mr)/cols,rh=(H-mt-mb)/rows;
const bx=(n,l)=>[ml+l*cw,mt+(n-1)*rh];let s='';
const pal=['var(--blue)','var(--yel)'];
for(let n=1;n<=7;n++)for(let l=0;l<Math.min(n,4);l++){const [x,y]=bx(n,l);const d=n+l;const idx=order.indexOf(n+L[l])+1;
 s+=`<rect x="${x+3}" y="${y+3}" width="${cw-6}" height="${rh-6}" rx="4" fill="${pal[d%%2]}" opacity="${d%%2?0.55:0.22}"/>`;
 s+=`<text x="${x+cw/2}" y="${y+rh/2+6}" text-anchor="middle" style="font-weight:700;font-size:${small?14:18}px">${n}${L[l]}</text>`;
 if(!small)s+=`<text x="${x+8}" y="${y+18}" style="font-size:13px">${idx}</text>`;
 const t=touch[n+L[l]];if(t)s+=`<circle cx="${x+cw-14}" cy="${y+16}" r="10" fill="var(--red)"/><text x="${x+cw-14}" y="${y+21}" text-anchor="middle" style="font-size:13px;fill:#fff;stroke:none">${t}</text>`;}
for(let l=0;l<4;l++)s+=`<text x="${ml+l*cw+cw/2}" y="${mt-10}" text-anchor="middle">l = ${l} (${L[l]})</text>`;
for(let n=1;n<=7;n++)s+=`<text x="${ml-10}" y="${mt+(n-0.5)*rh+5}" text-anchor="end">n = ${n}</text>`;
if(!small){const caps=[2,2,8,8,18,18,32,32];for(let d=1;d<=8;d++){const n=Math.min(d,7);const y=mt+(Math.min(d,7)-0.5)*rh;
 s+=`<text x="${W-mr+16}" y="${mt+(d-0.5)*(H-mt-mb)/8+5}">n + l = ${d}: holds ${caps[d-1]}</text>`;}}
el.innerHTML=s};const ro=new ResizeObserver(draw);ro.observe(el);draw();}
""" % (json.dumps(r["order"]), json.dumps(touch))
square_js = """
{const el=document.getElementById('square-fig');const draw=()=>{const W=el.clientWidth,H=el.clientHeight;const N=4;const gap=W<480?10:18;
const side=Math.min((W-gap*(N+1))/(1+2+3+4),H*0.8/4);const cols=['var(--blue)','var(--red)','var(--yel)','var(--green)'];let s='';let x0=gap;
for(let n=1;n<=N;n++){const y0=H*0.88-n*side;
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){const l=Math.max(i,j);s+=`<rect x="${x0+i*side+1}" y="${y0+(n-1-j)*side+1}" width="${side-2}" height="${side-2}" fill="${cols[l]}" opacity=".8"/>`}
 s+=`<text x="${x0+n*side/2}" y="${H*0.88+22}" text-anchor="middle">n = ${n}: ${n*n}</text>`;x0+=n*side+gap}
const lab=['s','p','d','f'];lab.forEach((t,i)=>{s+=`<rect x="${gap+i*70}" y="8" width="14" height="14" fill="${cols[i]}"/><text x="${gap+i*70+20}" y="20">${t}: ${2*i+1}</text>`});
el.innerHTML=s};const ro3=new ResizeObserver(draw);ro3.observe(el);draw();}
"""
sections = [
 dict(id="square", nav="Shells are squares", h="Each shell is a square, and its subshells are the square's gnomons",
  body="""<p>Hydrogen's shell n holds subshells l = 0 to n − 1 with 2l + 1 orbitals each: 1, 3, 5, 7. Those are the L-shaped gnomons that build an n × n square, corner plus two edges each time, so the shell holds n² orbitals:</p>
<div class="eq">n² = 1 + (n − 1) ∘₁ (n − 1)</div>
<p>The square is not a coincidence. Hydrogen's hidden symmetry is SO(4), which splits into two independent spins of the same size j, with 2j = n − 1. Each spin has n states, the two multiply to n · n, and decomposing j ⊗ j by the Clebsch–Gordan rule gives exactly l = 0 to n − 1. So each shell is the lift of a circle operation of two equal spins, with a geometric underpinning we can be sure of (Fock, 1935).</p>""",
  fig=True, figcls="", cap="Shells n = 1 to 4 as squares of orbitals; colour marks the subshell (s, p, d, f), each an L-shaped gnomon one wider than the last.", js=square_js),
 dict(id="periods", nav="Periods", h="The periods are the Madelung diagonals, and each holds a doubled square",
  body="""<p>Atoms fill subshells in order of n + l, then n (the Madelung rule). Each diagonal n + l = s holds 2 · ⌈s/2⌉² electrons: 2, 2, 8, 8, 18, 18, 32, 32. Those are the periods of Janet's left-step table exactly. The usual table merges the first two diagonals, giving periods 2, 8, 8, 18, 18, 32, 32 and noble gases at 2, 10, 18, 36, 54, 86, 118. Lean computes all three lists and proves that period 2k + 1 has the length of diagonal 2k + 2.</p>
<p>Why the periods come in doubled squares is the same gnomon count run along a diagonal instead of a row: a diagonal picks up one subshell of each odd size up to its length, and odd numbers sum to squares.</p>"""),
 dict(id="exceptions", nav="The exceptions", h="All 20 exceptions move one or two electrons to a neighbouring diagonal",
  body="""<p>Twenty elements have ground states that break the Madelung order (Wikipedia's Aufbau table: Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au, Ac, Th, Pa, U, Np, Cm, Lr). I compared each with the rule. Every one moves only one or two electrons, and always between subshells at most one diagonal apart. Ten of them never leave their diagonal at all: La, Ce, Gd, Ac, Th, Pa, U, Np and Cm swap between d and f on one diagonal, and lawrencium swaps 6d for 7p on diagonal 8. The other ten move an electron from an s subshell up one diagonal into d. No exception touches a p subshell except lawrencium's.</p>
<p>So the diagonals are sturdy and the exceptions are near-ties between neighbours, which matches the usual explanation (subshells on adjacent diagonals have nearly equal energies in heavy atoms). The red badges on the picture count how many exceptions touch each subshell.</p>"""),
 dict(id="dial", nav="The dial", h="The dial gives q-deformed shells, but I found no chemistry that uses them",
  body="""<p>The q-integer [n]<sub>q</sub> = 1 + q + … + q<sup>n−1</sup> is exactly n copies of 1 under the circle operation with k = q − 1, Brennan's dial q = k + 1. A q-deformed shell would hold [n]<sub>q</sub>² orbitals, and Negadi and Kibler built a q-deformed periodic table this way in the 1990s. Nothing in real atomic data asks for q ≠ 1: the shells are exact squares. This is a placement on the dial, not a finding.</p>"""),
]
claims = [
 ("Shell n holds n² = 1 + (n−1) ∘₁ (n−1) orbitals, built from gnomons", "tidy reading", "n ≤ 29; Lean <code>shell_square</code>, <code>gnomons_fill</code>"),
 ("The square is SO(4): two spins j with 2j = n − 1", "known", "Clebsch–Gordan series, n ≤ 19"),
 ("Madelung diagonals hold 2, 2, 8, 8, 18, 18, 32, 32 (Janet's periods)", "known", "Python; Lean <code>diagonals</code>, <code>diagCap_closed</code>"),
 ("Periods 2, 8, 8, 18, 18, 32, 32; noble gases at 2 … 118", "known", "Python; Lean <code>periods</code>, <code>noble_gases</code>"),
 ("Every one of the 20 exceptions moves ≤ 2 electrons across ≤ 1 diagonal", "new here", "checks.py against the 20 experimental configurations"),
 ("Ten exceptions stay on one diagonal; only Lr touches p", "new here", "checks.py"),
 ("q-integers are circle powers with k = q − 1", "tidy reading", "Python; Lean <code>q_integer</code>"),
]
lean = [("shell_square","(m+1)² = 1 + m ∘₁ m"),("gnomons_fill","the subshells fill the square"),("gnomon","each gnomon is corner + two edges"),("shell_capacity","2n² electrons"),
 ("two_spins","(a+1)(b+1) = 1 + a ∘₁ b"),("diagonals","the eight Madelung diagonals' capacities"),("diagCap_closed","diagonal s holds 2⌈s/2⌉²"),
 ("periods, noble_gases","the seven periods and noble-gas numbers"),("period_is_diagonal","period 2k+1 = diagonal 2k+2"),("q_integer, q_integer_classical","the dial gives q-integers; k = 0 gives n")]
n = build("index.html", title="The periodic law", tile=("3","Pt","Periodic law"),
 h1="Shells are squares, periods are diagonals, and every exception is a near-tie with a neighbour",
 lead="Each shell holds n² = 1 + (n − 1) ∘₁ (n − 1) orbitals because hydrogen's symmetry is two equal spins multiplied. The periods run along the Madelung diagonals, each holding a doubled square, and all twenty exceptions move one or two electrons at most one diagonal over.",
 figs=[("12 / 12","checks pass"),("12","Lean theorems, compiled"),("20 / 20","exceptions within one diagonal"),("2·⌈s/2⌉²","electrons on diagonal s")],
 picture=picture, caption="Subshells by n (rows) and l (columns), shaded by Madelung diagonal n + l; small numbers give the filling order. Red badges count the experimental exceptions that move electrons into or out of each subshell.",
 verdict_h="Lead check: not a lead",
 verdict="The square shells are SO(4), known since 1935, and the doubled-square periods are Janet's 1928 table. The exception pattern (always a neighbour on the diagonal, never more than two electrons) is a clean regularity I found no single statement of, but it is what the standard near-degeneracy explanation predicts, so it opens nothing new.",
 sections=sections, claims=claims, lean_file="chemistry/03-periodic/lean/ChemPeriodic.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/03-periodic/: <code>checks.py</code> (12 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemPeriodic.lean</code>, <code>make_page.py</code>.",
 sources='Sources: <a href="https://en.wikipedia.org/wiki/Aufbau_principle">Aufbau principle (Wikipedia)</a> for the 20 exceptions; Fock, Z. Phys. 98, 145 (1935) on SO(4); Janet’s left-step table (1928); Negadi and Kibler on the q-deformed periodic table (1990s).',
 series_nav="Chemistry pass, part 3 of 13: the periodic law. Before: conservation of mass and charge. Next: bonding and valence.")
print("page bytes", n)
