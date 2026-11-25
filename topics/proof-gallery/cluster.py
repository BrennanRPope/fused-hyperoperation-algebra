"""cluster.py: group the proofs by shape alone, name each group by its geometry,
place each group on the rung sheet, and find cross-field twins.

Shape of a proof = (body, frame, moves):
  body   from the statement's dimension (degree.py): point, rod, rectangle, box,
         hyperbox (4 to 9 sides) or a cube of any size
  frame  how the two pictures meet: same, inside, mismatch, tiles, witness, both ways
  moves  the engine moves read off the proof's own text (features.py, full weight)
Group = body + the proof's leading action, the first move present in a fixed
order (squash, dial, flatten, regroup, swap, explode, collapse, repeat; else plain).
Groups under 4 proofs are pooled by action; the field never enters.
Height on the rung sheet (s = 0 addition, 1 mult-add, 2 fused exponent, 3 and up
the towers) comes from the highest operation the statement uses.
Output: gallery.json
"""
import json, re, pathlib
from collections import Counter, defaultdict

HERE = pathlib.Path(__file__).resolve().parent
F = json.load(open(HERE / "features.json"))

ACTIONS = ["squash", "dial", "flatten", "regroup", "swap", "explode", "collapse", "repeat"]
SHAPE_MOVES = ACTIONS + ["lift", "multadd", "extrude", "tally"]

FIELD = {
    "FusedAlgebra.Base": ("Foundation", "foundation"),
    "FusedAlgebra.Rungs": ("Rung ladder", "foundation"),
    "FusedAlgebra.Cube": ("Mult-add cube", "foundation"),
    "FusedAlgebra.Topics.ClosingRate": ("Closing rates", "economy and society"),
    "FusedAlgebra.Topics.Consensus": ("Distributed consensus", "computing and information"),
    "FusedAlgebra.Topics.DarkSector": ("Dark matter (MOND)", "physics"),
    "FusedAlgebra.Topics.DataCenter": ("Data-center energy", "computing and information"),
    "FusedAlgebra.Topics.E": ("The number e", "numbers and infinity"),
    "FusedAlgebra.Topics.Finance": ("Finance", "economy and society"),
    "FusedAlgebra.Topics.FractionalRungs": ("Fractional rungs and solitons", "physics"),
    "FusedAlgebra.Topics.FusedConstants": ("Fused constants", "foundation"),
    "FusedAlgebra.Topics.Infinities": ("Infinities", "numbers and infinity"),
    "FusedAlgebra.Topics.InformationTheory": ("Information theory", "computing and information"),
    "FusedAlgebra.Topics.NavierStokes": ("Fluids (Navier–Stokes)", "physics"),
    "FusedAlgebra.Topics.QuantumNoise": ("Quantum noise", "physics"),
    "FusedAlgebra.Topics.Rulecraft": ("Rulecraft society model", "economy and society"),
    "FusedAlgebra.Topics.Toric.ToricCharts": ("Toric geometry", "geometry and topology"),
    "FusedAlgebra.Topics.Toric.ProjectivePlane": ("Projective plane", "geometry and topology"),
    "FusedAlgebra.Topics.Topology.CubeProofs": ("Primes on the cube", "numbers and infinity"),
    "FusedAlgebra.Topics.Topology.KummerWilson": ("Kummer and Wilson", "numbers and infinity"),
    "FusedAlgebra.Topics.Topology.HopfInvariantOne": ("Hopf invariant one", "geometry and topology"),
}

def body(d):
    if d is None: return "unclear"
    if d >= 99: return "cube"
    if d >= 4: return "hyper"
    return ["point", "rod", "rect", "box"][d]

BODY_NAME = {"point": "Point", "rod": "Rod", "rect": "Rectangle", "box": "Box",
             "hyper": "Hyperbox", "cube": "Cube of any size", "unclear": "Unread shape", "mixed": "Mixed sizes"}
ACTION_NAME = {
    "plain": "read straight off", "squash": "squashed onto the hyperbola", "dial": "with a dialed corner",
    "flatten": "flattened, area to length", "regroup": "cut two ways", "swap": "flipped on its diagonal",
    "explode": "broken into corner pieces", "collapse": "with a side shrunk to nothing",
    "repeat": "stacked layer on layer",
}
ACTION_SENT = {
    "plain": "The proof only reads the picture: no re-cut, flip, squash or dial is needed.",
    "squash": "A rod is rescaled to 1/(1+a), so a rectangle of area exactly 1 sits on the hyperbola xy = 1 (inverses, division, the wall at −1).",
    "dial": "The corner piece is weighted by k, so it shrinks, grows or flips sign (the a + b + k·ab family, solitons).",
    "flatten": "Area becomes length: a product turns into a sum of exponents, logs, digit sums or carries.",
    "regroup": "The same box is cut as slabs one way and as columns the other (associativity, appending lists).",
    "swap": "The picture is reflected across its diagonal, so the same pieces come in another order.",
    "explode": "The box is pulled apart into its 2ⁿ corner pieces and they are counted by layer (Pascal, factorials, corner sums).",
    "collapse": "A side of length 0 (or a factor of 1) flattens a whole dimension away.",
    "repeat": "One layer is laid down again and again (induction): the shape grows a slice at a time.",
}
FRAME_NAME = {"same": "the two pictures are the same", "inside": "one picture fits inside the other",
              "mismatch": "the pictures cannot match", "tiles": "the picture tiles evenly",
              "witness": "a piece is found that fits", "both_ways": "each picture holds the other"}

def rung_height(f):
    """highest operation in the statement, on the sheet's scale"""
    st = f["statement"]
    if re.search(r"\bfpent\b|\bpent\b|\bH\b", st): return 4
    if re.search(r"\bftet\b|\btet\b", st): return 3
    d = f["dim"]
    if d == 0:                                   # numbers only: the operation used sets the rung
        if re.search(r"\bfexp\b|\^|\bpw\b|\bwpow\b", st): return 2
        if re.search(r"\bmultAdd\b|\bwAdd\b|\*", st): return 1
        return 0
    if d is None: return 1
    if f["moves"].get("flatten", 0) >= 1: return -1 if d <= 1 else 0
    if d >= 99: return 2
    if d >= 2: return 1
    return 0

def ambiguity(f, own):
    why = []
    if f["dim"] is None: why.append("its statement could not be read, so its size is a guess")
    if not own: why.append("no move shows in its own text; it was closed by an automatic step, so only the statement gives its shape")
    halfonly = [m for m, v in f["moves"].items() if v == 0.5 and m in ACTIONS]
    if halfonly and not own:
        why.append("it touches " + ", ".join(halfonly) + " only through a definition or lemma it cites")
    if not f["moves"] and f["dim"] is not None and f["dim"] <= 1 and f["frame"] in ("tiles", "witness", "both_ways"):
        why.append("it is about chart or divisibility logic more than about a mult-add picture")
    return why

def main():
    rows = []
    for f in F:
        own = [m for m in ACTIONS if f["moves"].get(m, 0) >= 1]
        allm = sorted(m for m in SHAPE_MOVES if f["moves"].get(m, 0) >= 1)
        lead = own[0] if own else "plain"
        fld, grp = FIELD[f["module"]]
        rows.append(dict(name=f["name"], short=f["name"].split(".")[-1], module=f["module"], field=fld, group=grp,
                         line=f["line"], doc=f["doc"], statement=f["statement"], body=body(f["dim"]),
                         dim=f["dim"], frame=f["frame"], moves=allm, lead=lead,
                         touched=sorted(m for m, v in f["moves"].items() if v == 0.5),
                         height=rung_height(f), ambiguous=ambiguity(f, own)))
    # ---- groups
    key = lambda r: (r["body"], r["lead"])
    cnt = Counter(key(r) for r in rows)
    for r in rows:
        r["cluster"] = key(r) if cnt[key(r)] >= 4 else ("mixed", r["lead"])
    cnt2 = Counter(r["cluster"] for r in rows)
    for r in rows:
        if cnt2[r["cluster"]] < 3:
            r["cluster"] = ("mixed", "odd")
    groups = defaultdict(list)
    for r in rows: groups[r["cluster"]].append(r)
    border = ["point", "rod", "rect", "box", "hyper", "cube", "mixed", "unclear"]
    aorder = ["plain"] + ACTIONS + ["odd"]
    out_groups = []
    for (b, a), members in sorted(groups.items(), key=lambda kv: (aorder.index(kv[0][1]), border.index(kv[0][0]))):
        members.sort(key=lambda r: (r["frame"], r["moves"], r["field"], r["short"]))
        hs = sorted(r["height"] for r in members)
        if a == "odd":
            name = "Odd ones out"; sent = "Shapes that occur fewer than three times in the library, kept together so they are not forced into a group they do not fit."
        else:
            name = f"{BODY_NAME[b]}, {ACTION_NAME[a]}" if b != "mixed" else f"{ACTION_NAME[a].capitalize()}, small sizes"
            sent = ACTION_SENT[a]
        fields = Counter(r["field"] for r in members)
        fr = Counter(r["frame"] for r in members)
        out_groups.append(dict(id=f"{b}-{a}", body=b, action=a, name=name, sentence=sent, count=len(members),
                               height=hs[len(hs) // 2], fields=fields.most_common(), frames=fr.most_common(),
                               ngroups=len({r["group"] for r in members}),
                               members=[r["name"] for r in members]))
    for r in rows: r["cluster"] = f"{r['cluster'][0]}-{r['cluster'][1]}"
    # ---- twins: identical full shape, different broad fields (foundation left out)
    sig = lambda r: (r["body"], r["frame"], tuple(r["moves"]))
    by = defaultdict(list)
    for r in rows: by[sig(r)].append(r)
    twins = []
    for s, members in by.items():
        if len(s[2]) < 2: continue                     # need a specific shape, not a bare one
        ext = [r for r in members if r["group"] != "foundation" and not r["ambiguous"]]
        grps = defaultdict(list)
        for r in ext: grps[r["group"]].append(r)
        if len(grps) < 2: continue
        picks = [sorted(v, key=lambda r: (len(r["statement"]), r["short"]))[0] for k, v in sorted(grps.items())]
        twins.append(dict(signature=dict(body=s[0], frame=s[1], moves=list(s[2])), size=len(members),
                          ngroups=len(grps), picks=[p["name"] for p in picks]))
    twins.sort(key=lambda t: (-t["ngroups"], -len(t["signature"]["moves"]), t["size"]))
    data = dict(proofs=rows, groups=out_groups, twins=twins[:12],
                stats=dict(total=len(rows), groups=len(out_groups),
                           ambiguous=sum(1 for r in rows if r["ambiguous"]),
                           unclear=sum(1 for r in rows if r["dim"] is None),
                           shapes=len(by)))
    json.dump(data, open(HERE / "gallery.json", "w"), indent=1, ensure_ascii=False)
    print(data["stats"])
    for g in out_groups:
        print(f'{g["count"]:4d} s={g["height"]:2d} {g["ngroups"]}grp  {g["name"]}   {g["fields"][:3]}')
    print("--- twins")
    for t in twins[:12]:
        print(t["signature"], t["size"], [n.split(".")[-1] for n in t["picks"]])

if __name__ == "__main__":
    main()
