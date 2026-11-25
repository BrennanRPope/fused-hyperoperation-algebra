"""figures.py: the pictures for the shape-dependency page (SVG, labelled, theme-aware classes).
Reads orders.json, proof_uses.json, laws.json. Writes figures/*.svg."""
import json, math, pathlib, collections
H = pathlib.Path(__file__).resolve().parent
O = json.load(open(H / "orders.json")); U = json.load(open(H / "proof_uses.json")); L = json.load(open(H / "laws.json"))
OUT = H / "figures"; OUT.mkdir(exist_ok=True)
STYLE = """<style>
.ink{fill:#1b1d21}.mut{fill:#5d636c}.ln{stroke:#1b1d21}.fwd{stroke:#2f6db5}.bak{stroke:#d1452f}
.fwdt{fill:#2f6db5}.bakt{fill:#d1452f}.box{fill:#ffffff;stroke:#1b1d21}
.b-point{fill:#8a8f98}.b-rod{fill:#e3a21a}.b-rect{fill:#3f8f6b}.b-box{fill:#2f6db5}.b-hyper{fill:#7652b5}.b-cube{fill:#c2483a}
.s-point{stroke:#8a8f98}.s-rod{stroke:#e3a21a}.s-rect{stroke:#3f8f6b}.s-box{stroke:#2f6db5}.s-hyper{stroke:#7652b5}.s-cube{stroke:#c2483a}
text{font-family:Inter,system-ui,sans-serif}
@media (prefers-color-scheme:dark){.ink{fill:#eceef1}.mut{fill:#a9afb8}.ln{stroke:#eceef1}.box{fill:#1d2024;stroke:#eceef1}}
</style>"""
BNAME = {"point": "Point", "rod": "Stick", "rect": "Square", "box": "Box", "hyper": "Hyperbox", "cube": "Cube of any size"}

def icon(b, x, y, s=22):
    """small picture of a body centred at x, y"""
    c = f"s-{b}"; f = f"b-{b}"
    if b == "point": return f'<circle cx="{x}" cy="{y}" r="{s*0.25}" class="{f}"/>'
    if b == "rod": return (f'<line x1="{x-s}" y1="{y}" x2="{x+s}" y2="{y}" class="{c}" stroke-width="5" stroke-linecap="round"/>')
    if b == "rect":
        return (f'<rect x="{x-s}" y="{y-s}" width="{2*s}" height="{2*s}" class="{f}" opacity=".25"/>'
                f'<rect x="{x-s}" y="{y-s}" width="{2*s}" height="{2*s}" fill="none" class="{c}" stroke-width="3"/>'
                f'<rect x="{x}" y="{y-s}" width="{s}" height="{s}" class="{f}"/>')
    def cube(cx, cy, r, op=1):
        d = r * 0.5
        top = f"{cx-r},{cy-r+d*0} {cx-r+d},{cy-r-d} {cx+r+d},{cy-r-d} {cx+r},{cy-r}"
        side = f"{cx+r},{cy-r} {cx+r+d},{cy-r-d} {cx+r+d},{cy+r-d} {cx+r},{cy+r}"
        return (f'<g opacity="{op}"><rect x="{cx-r}" y="{cy-r}" width="{2*r}" height="{2*r}" class="{f}" opacity=".35"/>'
                f'<polygon points="{top}" class="{f}" opacity=".6"/><polygon points="{side}" class="{f}"/>'
                f'<rect x="{cx-r}" y="{cy-r}" width="{2*r}" height="{2*r}" fill="none" class="{c}" stroke-width="2"/></g>')
    if b == "box": return cube(x - 3, y + 3, s * 0.85)
    if b == "hyper":
        return cube(x - 8, y + 8, s * 0.6, .55) + cube(x + 5, y - 5, s * 0.6) + \
            f'<line x1="{x-8-s*.6}" y1="{y+8-s*.6}" x2="{x+5-s*.6}" y2="{y-5-s*.6}" class="{c}" stroke-width="1.5"/>' \
            f'<line x1="{x-8+s*.6}" y1="{y+8+s*.6}" x2="{x+5+s*.6}" y2="{y-5+s*.6}" class="{c}" stroke-width="1.5"/>'
    if b == "cube":
        return cube(x - 10, y + 10, s * .45, .5) + cube(x - 2, y + 2, s * .45, .75) + cube(x + 6, y - 6, s * .45) + \
            (f'<text x="{x+s+6}" y="{y+s}" class="ink" font-size="14" font-style="italic">n</text>' if s >= 15 else '')
    return ""

def svg(w, h, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{title}">'
            f'<title>{title}</title>{STYLE}{body}</svg>')

def arcs(order, pairs, label, fname, title, sub, icons=True, minw=1, wscale=0.22, gap=190):
    n = len(order); W = 70 + gap * (n - 1) + 130; Hh = 640; y0 = 330
    X = {k: 80 + gap * i for i, k in enumerate(order)}
    out = []
    mx = max(max(p["fwd"], p["back"]) for p in pairs)
    for p in pairs:
        a, b = p["a"], p["b"]
        for cnt, up in ((p["fwd"], True), (p["back"], False)):
            if cnt < minw: continue
            x1, x2 = X[a], X[b]; r = (x2 - x1) / 2; cx = (x1 + x2) / 2
            hgt = min(r * 0.6, 185)
            sw = max(1.2, cnt * wscale)
            if up:
                d = f"M{x1},{y0-34} C{x1},{y0-34-hgt*1.33} {x2},{y0-34-hgt*1.33} {x2},{y0-34}"
                out.append(f'<path d="{d}" fill="none" class="fwd" stroke-width="{sw:.1f}" opacity=".75"/>')
                out.append(f'<polygon points="{x2-6},{y0-46} {x2+6},{y0-46} {x2},{y0-34}" class="fwdt"/>')
                out.append(f'<text x="{cx}" y="{y0-34-hgt-6}" text-anchor="middle" font-size="13" class="fwdt">{cnt}</text>')
            else:
                d = f"M{x2},{y0+34} C{x2},{y0+34+hgt*1.33} {x1},{y0+34+hgt*1.33} {x1},{y0+34}"
                out.append(f'<path d="{d}" fill="none" class="bak" stroke-width="{sw:.1f}" stroke-dasharray="6 4" opacity=".8"/>')
                out.append(f'<polygon points="{x1-6},{y0+46} {x1+6},{y0+46} {x1},{y0+34}" class="bakt"/>')
                out.append(f'<text x="{cx}" y="{y0+34+hgt+16}" text-anchor="middle" font-size="13" class="bakt">{cnt}</text>')
    for k in order:
        out.append(f'<circle cx="{X[k]}" cy="{y0}" r="30" class="box" stroke-width="1.5"/>')
        if icons: out.append(icon(k, X[k], y0, 15))
        out.append(f'<text x="{X[k]}" y="{y0+4 if not icons else y0+50}" text-anchor="middle" font-size="14" font-weight="600" class="ink">{label(k)}</text>' if icons else
                   f'<text x="{X[k]}" y="{y0+5}" text-anchor="middle" font-size="12.5" font-weight="600" class="ink">{label(k)}</text>')
    out.append(f'<text x="20" y="28" font-size="18" font-weight="700" class="ink">{title}</text>')
    out.append(f'<text x="20" y="50" font-size="13.5" class="mut">{sub}</text>')
    out.append(f'<text x="{W-20}" y="80" text-anchor="end" font-size="13" class="fwdt">above: uses that follow the order (fundamental → derivative)</text>')
    out.append(f'<text x="{W-20}" y="{Hh-14}" text-anchor="end" font-size="13" class="bakt">below, dashed: uses against it</text>')
    (OUT / fname).write_text(svg(W, Hh, "".join(out), title))

ob = O["body"]
BORD = ["rod", "rect", "box", "hyper", "cube"]
bp = []
for i, a in enumerate(BORD):
    for b in BORD[i+1:]:
        f = next((p for p in ob["pairs"] if {p["a"], p["b"]} == {a, b}), None)
        if f is None: continue
        fwd, back = (f["fwd"], f["back"]) if f["a"] == a else (f["back"], f["fwd"])
        bp.append(dict(a=a, b=b, fwd=fwd, back=back))
arcs(BORD, bp, lambda k: BNAME[k], "bodies.svg",
     "Proof bodies stand on smaller bodies",
     f"{ob['total']} uses between bodies in the compiled library; {ob['back']} ({100*ob['back']/ob['total']:.0f}%) run backwards, against {100*ob['null_mean']:.0f}% by chance")
oa = O["action"]
ACT = {"plain": "read off", "repeat": "stack", "explode": "explode", "collapse": "collapse", "measure": "measure",
       "swap": "flip", "regroup": "cut 2 ways", "dial": "dial", "flatten": "flatten", "squash": "squash", "tally": "tally"}
arcs(oa["order"], [p for p in oa["pairs"] if p["fwd"] + p["back"] >= 12], lambda k: ACT[k], "actions.svg",
     "Proof moves: read-off facts at the bottom, dial and squash near the top",
     f"{oa['total']} uses between lead moves; {oa['back']} ({100*oa['back']/oa['total']:.0f}%) backwards, against {100*oa['null_mean']:.0f}% by chance; pairs used 12+ times in all are drawn, both ways",
     icons=False, wscale=0.18, gap=108)

# ---- strong family DAG (w >= 12)
fe = [e for e in U["family_edges"] if e["w"] >= 12 and e["src"] != e["dst"]]
nodes = sorted({e["src"] for e in fe} | {e["dst"] for e in fe})
succ = collections.defaultdict(list)
for e in fe: succ[e["src"]].append(e["dst"])
depth = {}
def dep(v):
    if v in depth: return depth[v]
    preds = [e["src"] for e in fe if e["dst"] == v]
    depth[v] = 0 if not preds else 1 + max(dep(p) for p in preds)
    return depth[v]
for v in nodes: dep(v)
layers = collections.defaultdict(list)
for v in nodes: layers[depth[v]].append(v)
FAM = lambda f: f.split("-")[0]
ACTN = lambda f: ACT.get(f.split("-")[1], f.split("-")[1])
Wd, Hd = 1360, 120 + 150 * (max(layers) + 1)
pos = {}
for d, vs in layers.items():
    vs.sort(key=lambda v: (["point", "rod", "rect", "box", "hyper", "cube", "mixed"].index(FAM(v)), v))
    for i, v in enumerate(vs):
        pos[v] = (Wd * (i + 1) / (len(vs) + 1), 110 + 150 * d)
o = []
for e in fe:
    (x1, y1), (x2, y2) = pos[e["src"]], pos[e["dst"]]
    sw = max(1.5, e["w"] * 0.12)
    o.append(f'<line x1="{x1:.0f}" y1="{y1+26:.0f}" x2="{x2:.0f}" y2="{y2-30:.0f}" class="fwd" stroke-width="{sw:.1f}" opacity=".6"/>')
    o.append(f'<circle cx="{x2:.0f}" cy="{y2-30:.0f}" r="3.5" class="fwdt"/>')
    o.append(f'<text x="{x1+(x2-x1)*0.62+5:.0f}" y="{y1+26+(y2-y1-56)*0.62+4:.0f}" font-size="12" class="fwdt">{e["w"]}</text>')
for v, (x, y) in pos.items():
    b = FAM(v)
    o.append(f'<rect x="{x-70:.0f}" y="{y-26:.0f}" width="140" height="52" rx="10" class="box" stroke-width="1.3"/>')
    o.append(icon(b, x - 46, y, 11) if b != "mixed" else "")
    o.append(f'<text x="{x-22:.0f}" y="{y-3:.0f}" font-size="13" font-weight="600" class="ink">{BNAME.get(b,b).replace(" of any size"," (any n)")}</text>')
    o.append(f'<text x="{x-22:.0f}" y="{y+14:.0f}" font-size="12.5" class="mut">{ACTN(v)}</text>')
o.append('<text x="20" y="30" font-size="18" font-weight="700" class="ink">The strong links between shape families form a clean order</text>')
o.append(f'<text x="20" y="52" font-size="13.5" class="mut">Family A → family B when 12 or more proofs in B use a theorem from A ({len(fe)} links, no cycles). Below 12 uses the families tangle into cycles.</text>')
(OUT / "families.svg").write_text(svg(Wd, Hd, "".join(o), "Strong family links"))

# ---- science laws graph
nodesL = L["nodes"]; edgesL = L["edges"]
BO = {"point": 0, "rod": 1, "rect": 2, "box": 3, "hyper": 3, "cube": 4}
predsL = collections.defaultdict(list)
for e in edgesL: predsL[e["b"]].append(e["a"])
dl = {}
def dL(v):
    if v not in dl: dl[v] = 0 if not predsL[v] else 1 + max(dL(p) for p in predsL[v])
    return dl[v]
for v in nodesL: dL(v)
col = collections.defaultdict(list)
for v in nodesL: col[dl[v]].append(v)
FO = {"physics": 0, "chemistry": 1, "biology": 2}
band = collections.defaultdict(list)
for v in nodesL: band[(FO[nodesL[v]["field"]], dl[v])].append(v)
rows = {f: max(len(band[(f, d)]) for d in range(max(dl.values()) + 1)) for f in range(3)}
top = {0: 130}; top[1] = top[0] + 78 * rows[0] + 50; top[2] = top[1] + 78 * rows[1] + 50
Wl = 1560; Hl = top[2] + 78 * rows[2] + 30
posL = {}
for (f, d), vs in band.items():
    vs.sort(key=lambda v: nodesL[v]["label"])
    for i, v in enumerate(vs):
        posL[v] = (80 + d * (Wl - 170) / max(dl.values()), top[f] + 20 + i * 78)
o = []
for e in edgesL:
    (x1, y1), (x2, y2) = posL[e["a"]], posL[e["b"]]
    a, b = nodesL[e["a"]]["body"], nodesL[e["b"]]["body"]
    cls = "fwd" if BO[a] < BO[b] else ("bak" if BO[a] > BO[b] else "ln")
    dash = ' stroke-dasharray="5 4"' if cls == "bak" else ""
    op = ".25" if cls == "ln" else ".8"
    o.append(f'<path d="M{x1+10:.0f},{y1:.0f} C{(x1+x2)/2:.0f},{y1:.0f} {(x1+x2)/2:.0f},{y2:.0f} {x2-10:.0f},{y2:.0f}" fill="none" class="{cls}" stroke-width="1.6" opacity="{op}"{dash}/>')
for v, (x, y) in posL.items():
    n = nodesL[v]; b = n["body"]
    o.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="9" class="b-{b}"/>')
    lab = n["label"].split(" (")[0]
    lab = {"Conservation of mass and charge": "Mass, charge kept", "Equilibrium, mass action": "Equilibrium"}.get(lab, lab)
    if len(lab) > 22: lab = lab.split(",")[0]
    o.append(f'<text x="{x:.0f}" y="{y-14:.0f}" text-anchor="middle" font-size="11.5" class="ink">{lab}</text>')
    o.append(f'<text x="{x:.0f}" y="{y+23:.0f}" text-anchor="middle" font-size="10.5" class="mut">{BNAME[b].replace(" of any size","")} {100*n["body_share"]:.0f}%</text>')
f = L["flow"]
for fi, nm in ((0, "Physics"), (1, "Chemistry"), (2, "Biology")):
    o.append(f'<text x="20" y="{top[fi]-22}" font-size="13" font-weight="700" class="mut">{nm}</text>')
    if fi: o.append(f'<line x1="16" y1="{top[fi]-40}" x2="{Wl-16}" y2="{top[fi]-40}" class="ln" stroke-width=".6" opacity=".3"/>')
o.append('<text x="20" y="30" font-size="18" font-weight="700" class="ink">Whole science laws: their shapes do not line up with what they stand on</text>')
o.append(f'<text x="20" y="52" font-size="13.5" class="mut">{len(edgesL)} textbook links, left (fundamental) to right. Blue: to a bigger body ({f["up"]}). Red dashed: to a smaller body ({f["down"]}). Grey: same layer ({f["same"]}). Chance would give as many reds as blues.</text>')
lx = 20
for b in ["point", "rod", "rect", "box", "hyper", "cube"]:
    o.append(f'<circle cx="{lx+8}" cy="72" r="7" class="b-{b}"/><text x="{lx+20}" y="76" font-size="12.5" class="ink">{BNAME[b]}</text>')
    lx += 30 + 8 * len(BNAME[b]) + 18
(OUT / "laws.svg").write_text(svg(Wl, Hl, "".join(o), "Science laws by shape"))

# ---- zoom: Brennan's chain
chain = [("prod", 0, 0), ("em", 1, 0), ("sr", 2, 0), ("gr", 3, 0), ("gpb", 4, 0), ("thomas", 3, 1), ("pulsar", 4, 1)]
P = {k: (90 + 190 * c, 190 + 150 * r) for k, c, r in chain}
links = [("prod", "em"), ("em", "sr"), ("prod", "sr"), ("sr", "gr"), ("gr", "gpb"), ("sr", "thomas"), ("thomas", "gpb"), ("gr", "pulsar"), ("thomas", "pulsar")]
o = []
for a, b in links:
    (x1, y1), (x2, y2) = P[a], P[b]
    ba, bb = nodesL[a]["body"], nodesL[b]["body"]
    cls = "fwd" if BO[ba] < BO[bb] else ("bak" if BO[ba] > BO[bb] else "ln")
    dash = ' stroke-dasharray="6 4"' if cls == "bak" else ""
    if a == "prod" and b == "sr":
        o.append(f'<path d="M{x1},{y1-40} C{x1+80},{y1-110} {x2-80},{y2-110} {x2},{y2-40}" fill="none" class="{cls}" stroke-width="2.4"{dash} opacity=".8"/>')
    else:
        o.append(f'<line x1="{x1+44}" y1="{y1+(8 if y2>y1 else 0)}" x2="{x2-44}" y2="{y2-(8 if y2>y1 else 0)}" class="{cls}" stroke-width="2.4"{dash} opacity=".85"/>')
for k, (x, y) in P.items():
    n = nodesL[k]; b = n["body"]
    o.append(f'<rect x="{x-62}" y="{y-40}" width="124" height="80" rx="12" class="box" stroke-width="1.3"/>')
    o.append(icon(b, x, y - 8, 14))
    o.append(f'<text x="{x}" y="{y+56}" text-anchor="middle" font-size="13" font-weight="600" class="ink">{n["label"].split(" (")[0]}</text>')
    o.append(f'<text x="{x}" y="{y+72}" text-anchor="middle" font-size="12" class="mut">{BNAME[b].replace(" of any size","")}, {100*n["body_share"]:.0f}% of its proofs</text>')
o.append('<text x="20" y="30" font-size="18" font-weight="700" class="ink">Zoom: electromagnetism → special relativity → gravity</text>')
o.append('<text x="20" y="52" font-size="13.5" class="mut">The chain is real (textbook). The shapes along it go down as often as up: hyperbox → square → square → point.</text>')
o.append('<line x1="20" y1="76" x2="50" y2="76" class="fwd" stroke-width="2.4"/><text x="56" y="80" font-size="12.5" class="ink">to a bigger body</text>'
         '<line x1="190" y1="76" x2="220" y2="76" class="bak" stroke-width="2.4" stroke-dasharray="6 4"/><text x="226" y="80" font-size="12.5" class="ink">to a smaller body</text>'
         '<line x1="370" y1="76" x2="400" y2="76" class="ln" stroke-width="2.4"/><text x="406" y="80" font-size="12.5" class="ink">same body</text>')
(OUT / "chain_physics.svg").write_text(svg(1000, 520, "".join(o), "Zoom on the physics chain"))

# ---- zoom: the exact chain of composition shapes
o = []
def node(x, y, t1, t2, t3, kind):
    o.append(f'<rect x="{x-110}" y="{y-48}" width="220" height="96" rx="14" class="box" stroke-width="1.5"/>')
    # tiny wall picture: number line with walls
    o.append(f'<line x1="{x-80}" y1="{y-18}" x2="{x+80}" y2="{y-18}" class="ln" stroke-width="1.5"/>')
    if kind == "two":
        for wx in (x - 40, x + 40): o.append(f'<line x1="{wx}" y1="{y-30}" x2="{wx}" y2="{y-6}" class="bak" stroke-width="3"/>')
    elif kind == "border":
        o.append(f'<line x1="{x}" y1="{y-32}" x2="{x}" y2="{y-4}" class="bak" stroke-width="5"/>')
    elif kind == "none":
        o.append(f'<circle cx="{x}" cy="{y-18}" r="14" fill="none" class="fwd" stroke-width="2.5"/>')
    elif kind == "torus":
        o.append(f'<ellipse cx="{x}" cy="{y-18}" rx="34" ry="13" fill="none" class="fwd" stroke-width="2.5"/><ellipse cx="{x}" cy="{y-20}" rx="12" ry="4" fill="none" class="fwd" stroke-width="2"/>')
    o.append(f'<text x="{x}" y="{y+10}" text-anchor="middle" font-size="15" font-weight="700" class="ink">{t1}</text>')
    o.append(f'<text x="{x}" y="{y+27}" text-anchor="middle" font-size="12.5" class="mut">{t2}</text>')
    o.append(f'<text x="{x}" y="{y+42}" text-anchor="middle" font-size="12" class="mut">{t3}</text>')
def arrow(x1, y1, x2, y2, cls, t, t2="", dash=False, dy=-10):
    d = ' stroke-dasharray="6 4"' if dash else ""
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="{cls}" stroke-width="2.6"{d}/>')
    ang = math.atan2(y2 - y1, x2 - x1); ax, ay = x2, y2
    p1 = (ax - 12 * math.cos(ang - .4), ay - 12 * math.sin(ang - .4)); p2 = (ax - 12 * math.cos(ang + .4), ay - 12 * math.sin(ang + .4))
    o.append(f'<polygon points="{ax},{ay} {p1[0]:.1f},{p1[1]:.1f} {p2[0]:.1f},{p2[1]:.1f}" class="{cls}t"/>')
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    o.append(f'<text x="{mx}" y="{my+dy}" text-anchor="middle" font-size="12.5" class="{cls}t">{t}</text>')
    if t2:
        ty = my + dy + 15 if (dy > 0 or x1 == x2) else my + dy - 15
        if y1 == y2 and dy < 0: o[-1] = o[-1].replace(f'y="{my+dy}"', f'y="{my+dy}"'); ty = my + dy - 16
        o.append(f'<text x="{mx}" y="{ty}" text-anchor="middle" font-size="11" class="mut" font-style="italic">{t2}</text>')
node(150, 190, "Addition", "one double wall (border)", "0 periods: a + b", "border")
node(560, 190, "The circle operation", "two real walls", "1 period: a + b + ab", "two")
node(970, 190, "Torus law", "elliptic curve", "2 periods", "torus")
node(560, 470, "Turning law", "no real wall", "slope addition (a+b)/(1−ab)", "none")
arrow(262, 175, 446, 175, "fwd", "add the corner ab", "law_is_addition_plus_corner")
arrow(446, 212, 262, 212, "bak", "walls merge (limit)", "two_walls_difference", dash=True, dy=22)
arrow(672, 175, 856, 175, "fwd", "two periods", "textbook (theta functions)")
arrow(856, 212, 672, 212, "bak", "pinch once: node", "nodal_collinear, nodal_phi", dash=True, dy=22)
arrow(560, 240, 560, 420, "fwd", "", "")
o.append('<text x="572" y="372" font-size="12.5" class="fwdt">turn by i</text><text x="572" y="388" font-size="11" class="mut" font-style="italic">turn_num, turn_disc</text>')
o.append('<path d="M1010,240 C1010,330 190,330 190,240" fill="none" class="bak" stroke-width="2.2" stroke-dasharray="6 4"/><polygon points="183,252 197,252 190,240" class="bakt"/>')
o.append('<text x="760" y="334" font-size="12.5" class="bakt">pinch twice: cusp → addition</text><text x="760" y="350" font-size="11" class="mut" font-style="italic">cusp_collinear, cusp_is_addition</text>')
o.append('<text x="20" y="30" font-size="18" font-weight="700" class="ink">Zoom: the composition shapes form one clean line, proved</text>')
o.append('<text x="20" y="52" font-size="13.5" class="mut">Solid blue: what you need to build the next shape. Red dashed: the limit that takes it back. Lean theorem names in italics (21 theorems, no sorry).</text>')
o.append('<text x="20" y="74" font-size="13.5" class="mut">Zooming never changes the shape (zoom_disc), and the border is fixed by the turn (turn_keeps_border).</text>')
(OUT / "chain_exact.svg").write_text(svg(1130, 550, "".join(o), "Exact chain of composition shapes"))

# ---- zoom: the exception (any-size law used at one size)
ex = O["ex_cube_rod"]
o = []
o.append('<text x="20" y="30" font-size="18" font-weight="700" class="ink">Zoom: the main exception runs from "any size" down to a stick</text>')
o.append(f'<text x="20" y="52" font-size="13.5" class="mut">{len(ex)} proofs about one fixed size use a law stated for every size n, then set n. That is specialising, not building.</text>')
o.append(icon("cube", 120, 170, 34)); o.append(icon("rod", 640, 170, 34))
o.append('<text x="120" y="232" text-anchor="middle" font-size="14" font-weight="600" class="ink">a law for every n</text>')
o.append('<text x="640" y="232" text-anchor="middle" font-size="14" font-weight="600" class="ink">used at one n</text>')
o.append('<line x1="200" y1="170" x2="560" y2="170" class="bak" stroke-width="3" stroke-dasharray="7 5"/><polygon points="560,162 574,170 560,178" class="bakt"/>')
pick = [e for e in ex if any(s in e["user"] for s in ("three_dvd_of_cube", "alternating_corners", "nilpotent_weight_group", "frey_disc_valuation"))][:4]
y = 280
for e in pick:
    o.append(f'<text x="40" y="{y}" font-size="12.5" class="ink">{e["needed"].split(".")[-1]}</text>')
    o.append(f'<text x="330" y="{y}" font-size="12.5" class="bakt">→</text>')
    o.append(f'<text x="360" y="{y}" font-size="12.5" class="ink">{e["user"].split(".")[-1]}</text>')
    o.append(f'<text x="40" y="{y+16}" font-size="11.5" class="mut">{e["stmt"][:70].replace("<","&lt;").replace(">","&gt;").replace("&","&amp;") if False else e["stmt"][:70].replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")}</text>')
    y += 44
(OUT / "exception.svg").write_text(svg(780, y + 10, "".join(o), "The main exception"))
print("figures:", sorted(p.name for p in OUT.glob("*.svg")))
