import sys, importlib.util
sys.path.insert(0, "/mnt/project-files/math-theory/prior-art-audit")
spec = importlib.util.spec_from_file_location("claims", "/mnt/project-files/math-theory/prior-art-audit/claims.py")
cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
D27, D28 = "2026-11-27", "2026-11-28"
# partly-known rows where the audit says the link itself wasn't found -> Yours (link)
LINK = {"A2b":"q>0 sandwich on S3 not found stated","A8":"sizing 600-cell Regge universes by the law",
 "B9b":"equality not found (implied by Biedenharn 1983)","B11":"probability reading not found",
 "L1":"the number sqrt6-3 not found printed","S85":"the one-line picture","M189":"the constant not found (Zolotarev 1957 unread)",
 "NT":"explicit families and plane version (line formula Barany-Verma 2025)","G2":"copula link (line is Felsenstein 1979)",
 "N11":"the identity (Owen 1972; actuarial UDD)","G5":"Lee-Yang radius test (zeros: Briggs 1985)",
 "P26":"per-run prices (Suvakov-Dmitrasinovic 2015)","A03":"cross-link","A10":"cross-link","G1":"cross-link","A124":"cross-link",
 "G3":"cross-link","G4":"cross-link","A63":"cross-link","A67":"cross-link","A71":"cross-link","A83":"cross-link","P137":"cross-link"}
SRC = {"A1":"Nesterov-Sabinin 2000-04","A7":"Bachmann-Becigneul-Ganea 2020","A2a":"Kim 2020","A3":"Ungar; Nesterov-Sabinin",
 "A4":"Nesterov-Sabinin PRD 2000","A5":"de Sitter 1917 (cosh restated)","A6":"Tegmark 1997; Dorling 1970","B1":"Cork et al. 2021; Gibbons-Steif",
 "B2":"Sussman 2013","B3":"Bartnik-McKinnon 1988","B4":"Ranada 1989","B5a":"Abbassi et al. 2013","B6":"Markovic-Shapiro 2000",
 "B7":"any GR text","B8":"circle-group commutator; Pereira et al. 2001","B9a":"standard ISCO; a rewrite","B12":"QCD","B14":"textbook diagnosis",
 "C2":"chosen, not forced (rounds 27-28)","B10":"Calcagni 2012; Adler 2011","P153":"Rahn-Doherty-Mabuchi 2002","M106":"Merimovich (one line)",
 "M25":"Spivak 2022 (Garner)","M188":"classical 1-b/a^2","M224":"Karlin-McGregor 1968","S44":"Buchstaber-Panov (Todd genus)",
 "NF8":"OEIS A008826; Comtet 1974","A76":"Cajueiro 2006","P180":"Cheskidov 2008","M264":"folklore level","S50":"folklore level",
 "M153":"folklore level","M09":"folklore level","M212":"folklore level","WS":"folklore level","A29":"cross-link, thin search","P174":"thin search"}
rows = []  # (id, item, group, note, where, date)
for cid, area, label, status, conf, pred in cl.C:
    if status == "print": g, n = "Known", SRC.get(cid, "in print")
    elif status == "new":
        if conf == "low": g, n = "Unsure", "nothing found, but the search was thin"
        else: g, n = "Yours", "nothing found after a real search"
    else:
        if cid in LINK: g, n = "Yours", "link: " + LINK[cid]
        elif conf == "low": g, n = "Unsure", SRC.get(cid, "partly known, thin search")
        else: g, n = "Known", "your framing; " + SRC.get(cid, "pieces in print")
    rows.append((cid, label, g, n, "prior-art-audit", D27))
new = [
 ("K1","Core law a+b+kab","Known","Perlis 1942, Jacobson 1945; Mathlib PreQuasiregular. The dial k is your framing","mathlib-map",D28),
 ("K2","Every k != 0 gives the same monoid","Known","one-line rescaling","mathlib-map",D28),
 ("K3","Two-dial law (a+b+kab)/(1-d ab)","Known","Coleman-McGuinness 1991 rational formal group laws; Todd genus","why-i/cayley-klein",D28),
 ("K4","Extended law with c (helix round)","Known","same family as K3: c is the repo's existing d dial, not a new dial","bridges-study/helix-from-law",D28),
 ("K5","Every (k,d) law is the tangent law, kappa = -Delta/4","Known","your framing; C-M classify by k^2-4d","why-i/cayley-klein",D28),
 ("K6","Kerr helix in imaginary time","Known","Gibbons-Hawking 1977","bridges-study/c2-one-body",D28),
 ("K7","Kerr T and S = A/4","Known","Bekenstein 1973, Hawking 1975","bridges-study",D28),
 ("K8","Page curve","Known","Page 1993","bridges-study",D28),
 ("K9","Renyi of a black hole (closed: area exponent = p dial)","Known","Renyi 1961","bridges-study",D28),
 ("K10","Local tomography = law at k = +1","Known","your framing (a restatement); Masanes-Muller 2011, Hardy 2001","spin-i-tie/mm-map",D28),
 ("K11","Smooth turns need complex steps; k = -2 pin","Known","your framing; Hardy 2001, Rabi formula","spin-i-tie/smooth-turns",D28),
 ("K12","Anticommuting flips give 3D spin, third axis needs i","Known","your framing; Pauli 1927","spin-i-tie/spin-3d",D28),
 ("K13","Max uncertainty iff anticommuting","Known","your framing; Maassen-Uffink 1988, Wootters-Fields 1989","spin-i-tie/uncertainty-sign",D28),
 ("K14","Coin rule forces k = 0, then Shannon","Known","your framing and Lean packaging; Faddeev 1956","spin-i-tie/lean-proof",D28),
 ("K15","Info causality needs more than the coin rule","Known","page's own verdict: nothing new","spin-i-tie/ic-from-coin",D28),
 ("K16","Only the exact disk passes the polygon test","Known","Barnum et al. 2009","spin-i-tie/polygon-entropy",D28),
 ("K17","Fixed cell volume: vacuum energy doesn't gravitate","Known","your framing (as why-i's A4); Unruh 1989, Henneaux-Teitelboim 1989","why-i/unimodular",D28),
 ("K18","Unimodular leak fit (BIC prefers plain)","Known","models Josset-Perez-Sudarsky 2017; the fit is our computation","why-i/unimodular-fit",D28),
 ("K19","Torsion is derived, not a dial","Known","Cartan-Schouten 1926","forces rounds",D28),
 ("K20","Chain dial k equals q = 1 + k (q-integers)","Known","q-analogs; Knauf's chain","riemann",D28),
 ("K21","Field energy as its own source","Known","Kraichnan 1955, Deser 1970","gravity-from-theory/spinning-river",D28),
 ("K22","Hyperoperation sequence","Known","Goodstein 1947; Ackermann 1928; Knuth 1976","hyperoperation-ladder",D28),
 ("K23","Commutative ladder of operations","Known","Bennett 1915; Jaramillo 2020 (arXiv 2012.04093)","hyperoperation-ladder",D28),
 ("K24","Rung n+1 fails to distribute over rung n","Known","as a fact (Venkatesh, AMS column 2022)","hyperoperation-ladder",D28),
 ("Y1","Fused whole rung n = Goodstein's H_n through x -> x+1, unique","Yours","the tie plus its uniqueness theorem; the shift itself is one line. Backs the name 'Fused Hyperoperation Algebra' for whole rungs","hyperoperation-ladder",D28),
 ("Y2","Two ladders meet at exactly one step, a+b+ab","Yours","not found in a search today (thin)","hyperoperation-ladder",D28),
 ("Y3","Law board for every rung, 80 Lean theorems, mod-7 proof on giant rungs","Yours","formal work; the mod-7 trick is standard, its use here is new","hyperoperation-ladder",D28),
 ("Y4","Lean library, 3,954 theorems","Yours","formal work","lean-library",D28),
 ("Y5","Lee-Yang holds iff a <= p on the roundness dial; p = 3 and near zone to p = 8 proved","Yours","page's search: not found in print; Kuttner/Golubov neighbour; a specialist should read","shape-dial-survey",D28),
 ("Y6","Three-slit cap on the dial, |k| < 0.07","Yours","our bound from Kauten 2017 data","spin-i-tie/three-slit-cap",D28),
 ("U1","Helix = rim speed composed with itself at d = +-1","Unsure","reading not lit-checked","bridges-study/helix-from-law",D28),
 ("U2","Unrelated coin scales leftover by 1 + kS (mixing forces k = 0)","Unsure","page: may be written somewhere","spin-i-tie/linear-mixing",D28),
 ("U3","Law at k != 0 interferes at every Sorkin order","Unsure","page: as far as found","spin-i-tie/interference-law",D28),
 ("U4","Reflection positivity on the real dial; RP gives the energy-zero step","Unsure","page: link not found, short step","why-i/reflection-positivity",D28),
 ("U5","Only the backwards-reading sum has no preferred hand; eta(0) closed form","Unsure","no lit check recorded","forces round 32",D28),
 ("U6","Extra flex on the Freudenthal gluing; graviton cone = lattice wave","Unsure","page: not found by us","forces round 34",D28),
 ("U7","SO(10) route: law on the 10-sphere, one hand = one generation","Unsure","group textbook (Georgi 1975); the route not lit-checked","physics-bridges",D28),
 ("U8","Light-speed lag cancels pi^2; 41/16 must come from field energy","Unsure","short search only","gravity-from-theory/spinning-river",D28),
 ("U9","k dial on Knauf's chain gives sum phi(m)[m]_q^-beta","Unsure","quick search only","riemann/k-dial-cusp",D28),
 ("U10","Correction to Gupta's 2PN Lagrangian","Unsure","check Eq. (7.1) before it counts","hypersphere/orbits",D28),
 ("U11","6 three-body loop types not in Li-Liao's list","Unsure","newer catalogs not checked","three-body",D28),
 ("U12","Decaying flow on a Klein bottle ends as one vortex","Unsure","numerical, no lit check","klein-flow",D28),
]
rows += new
from collections import Counter
cnt = Counter(r[2] for r in rows)
by = Counter((r[2], r[5]) for r in rows)
print(cnt, by)
import json; json.dump(rows, open("/tmp/kvy/rows.json","w"))
