"""features.py: read each proof's shape off its Lean source, mechanically.

The rulebook below maps what a proof touches (its statement, its proof steps,
the library definitions it unfolds and the library lemmas it cites) to the
proof engine's nine moves, plus four marks the engine page does not list:
  explode  the box cut into its 2^n corner pieces (Pascal layers, corner sums)
  tally    checked by counting cells (decide / rfl on numbers), no move drawn
  repeat   induction or recursion: one layer laid down again and again
  lists    the box has any number of sides (a List or Fin n of rods)
Nothing is hand-assigned per theorem: change a rule and every proof is redone.
Output: features.json
"""
import json, re, pathlib
from collections import Counter

HERE = pathlib.Path(__file__).resolve().parent
LIB = HERE.parent / "lean-library" / "FusedAlgebra"
T = json.load(open(HERE / "theorems.json"))

MOVES = ["lift", "multadd", "swap", "regroup", "collapse", "squash", "dial", "flatten", "extrude",
         "explode", "repeat", "tally"]

# ---------------------------------------------------------------- definitions
DEF = re.compile(r"^(?:private |noncomputable )*(?:def|abbrev)\s+(?:\{[^}]*\}\s*)?([\w.']+)", re.M)

def module_of(p):
    return ".".join(p.relative_to(LIB.parent).with_suffix("").parts)

defs = {}   # (module, short) -> body text ; short -> list of modules
for p in LIB.rglob("*.lean"):
    txt = p.read_text()
    txt = re.sub(r"/-.*?-/", " ", txt, flags=re.S)
    txt = re.sub(r"--[^\n]*", " ", txt)
    ms = list(DEF.finditer(txt))
    for i, m in enumerate(ms):
        end = len(txt)
        nxt = re.search(r"^\S", txt[m.end():], re.M)
        if nxt: end = m.end() + nxt.start()
        short = m.group(1).split(".")[-1]
        defs.setdefault(short, {})[module_of(p)] = txt[m.start():end]

FOUNDATION = ["FusedAlgebra.Base", "FusedAlgebra.Rungs", "FusedAlgebra.Cube"]

def def_body(short, module):
    d = defs.get(short)
    if not d: return None
    if module in d: return d[module]
    for f in FOUNDATION:
        if f in d: return d[f]
    if len(d) == 1: return next(iter(d.values()))
    return None

thm_by_short = {}
for t in T:
    thm_by_short.setdefault(t["short"], []).append(t)

def cited_statement(short, module):
    c = thm_by_short.get(short)
    if not c: return None
    for t in c:
        if t["module"] == module: return t["statement"]
    for t in c:
        if t["module"] in FOUNDATION: return t["statement"]
    return c[0]["statement"]

IDENT = re.compile(r"[A-Za-z_][\w'.]*")

# ---------------------------------------------------------------- the rulebook
# Each rule: (move, regex). Matched against a proof's "reach": its own text,
# plus the bodies of library definitions it names, plus the statements of
# library lemmas it cites. Order does not matter.
RULES = [
    ("lift",     r"_lift\b|\blift_inj\b|\blifts?\b|\bslide\b|\bfuse\b|\(1 \+ [^()]+\)|\([^()]+ \+ 1\)|\(1 - [^()]+\)|\(o\d - 1\)|- 1 =|1 \+ multAdd|1 \+ k \*"),
    ("multadd",  r"\bmultAdd\w*|⊞|\bmadd\b|\bfoldMA\w*|\bmpow\w*|\bmaddOnes\b|\bvariable_node\b"),
    ("swap",     r"_comm\b|\bcomm\b|_perm\b|\bperm\b|\bmul_comm\b|\badd_comm\b|\.symm\b|\bsymm\b|\bflip\b"),
    ("regroup",  r"_assoc\b|\bassoc\b|_append\b|\bappend\b|\bmul_assoc\b|\badd_assoc\b|\bouter\b"),
    ("collapse", r"_zero\b|\bzero_\w+|_one\b|\bone_\w+|\bmul_one\b|\badd_zero\b|multAdd \w+ 0\b|multAdd 0 |\bwAdd 0 |\bsub_self\b|\bsub_add_cancel\b"),
    ("squash",   r"[Ii]nv\b|⁻¹|\bmsub\b|/\s*\(?\s*[A-Za-z(]|\bdiv\w*|\bneg_one_absorbs\b|\bwall\b|\bmultAddInv\b|\btrans1?2?\b|\breciprocal\b"),
    ("dial",     r"\bwAdd\w*|\bwpow\w*|\bwOnes\b|\bweight\w*|\bclimb\b|\bscale\w*|\bk \* [a-z]|\bwkey\b|\bvadd\b|\bwt\b"),
    ("flatten",  r"\bpow_add\w*|\bpow_mul\w*|\bmpow_add\b|\bpw_add\b|\bpw_mul\b|\blog\w*|\bdsum\w*|\bcarries\w*|\bval\b|\bVal\b|\bgrowth\b|\bmall\b|\bdigitProd\b"),
    ("extrude",  r"\^|\bpw\b|\bpw_\w+|\bpow\w*|\bmpow\w*|\bwpow\w*|\bfexp\w*|\bftet\w*|\bfpent\w*|\btet\w*|\bpent\w*|\bstair\b|\bpsi\w*|\bwind\b|\bgsum\b|\biter\b|\brepeat\b"),
    ("explode",  r"\bcorner\w*|\bgroups\b|\blayers?\b|\blayer_\w+|\bfact\w*|\bchoose\b|\bbinom\w*|\bPascal\b|\bsigns\b|\bminuses\b|\bblockCount\b|\bsumL\b|\bprodL\b|\bfall\b|\btri\b|\bS_pow\b|\bS\b"),
    ("repeat",   r"\binduction\b|\brecursion\b|\bNat\.rec\b|^\s*\| (?:n|k|m|0|zero|succ|nil|cons)\b"),
    ("tally",    r"\bdecide\b|\bnative_decide\b"),
]
RULES = [(m, re.compile(r, re.M)) for m, r in RULES]

def reach(t):
    own = t["statement"] + "\n" + t["proof"]
    parts = [own]
    seen = set()
    for idt in set(IDENT.findall(own)):
        short = idt.split(".")[-1]
        if short in seen or short == t["short"]: continue
        seen.add(short)
        b = def_body(short, t["module"])
        if b: parts.append(b)
        s = cited_statement(short, t["module"])
        if s: parts.append(s)
    return "\n".join(parts), own

def frame(stmt):
    s = stmt
    if re.search(r"\bFalse\b|≠|¬|\bno_\w*", s): return "mismatch"
    if re.search(r"∣|%|\bMEq\b|≡|\bIsPrime\b|\bOneMod\b|\bIsOdd\b|\bVal\b", s): return "tiles"
    if re.search(r"↔", s): return "both_ways"
    if re.search(r"∃|\bInChart\w*", s): return "witness"
    if re.search(r"≤|<|≥|>", s): return "inside"
    if re.search(r"=", s): return "same"
    return "same"

HYP = re.compile(r"[=≤<≠∣↔∧∨¬→]|\bProp\b|\bIsPrime\b|\bOneMod\b|\bIsOdd\b|\bVal\b|\bMEq\b")

def rods(t):
    """count the rods (sides of the box) from the binders"""
    b = t["binders"]
    stmt = t["statement"]
    n = 0; lists = False; counts = 0
    for grp in re.finditer(r"[({\[]([^(){}\[\]]*?):([^(){}\[\]]*)[)}\]]", b):
        names, ty = grp.group(1).split(), grp.group(2).strip()
        if any(nm.startswith("h") and len(nm) <= 4 and (nm[1:].isdigit() or len(nm) <= 3) for nm in names) and HYP.search(ty):
            continue
        if HYP.search(ty) or ty in ("Type", "Type u") or ty.startswith("Comm") or ty.startswith("Field") or "Ring" in ty or "Semiring" in ty:
            continue
        if "List" in ty or "Fin " in ty or "→" in ty:
            lists = True; continue
        for nm in names:
            if not re.search(rf"(?<![\w']){re.escape(nm)}(?![\w'])", stmt):
                continue
            if re.search(rf"\^\s*\(?{re.escape(nm)}(?![\w'])", stmt) or \
               re.search(rf"\b(?:pw|mpow|wpow|wOnes|layers|tet|pent|iter|repeat|stair|gsum|fact|rsum|s|u|psi\w*)\b[^=≤<]*?\b{re.escape(nm)}\b\s*(?:[)=≤<]|$)", stmt) and ty == "Nat":
                counts += 1; continue
            n += 1
    if re.search(r"\bList\b|\[.*?,.*?\]", stmt): lists = True
    return n, counts, lists


# ---------------------------------------------------------------- dimension
import sys
sys.path.insert(0, str(HERE))
import degree as D

def parse_defs():
    out = {}
    for short, mods in defs.items():
        for mod, body in mods.items():
            m = re.match(r"^(?:private |noncomputable )*(?:def|abbrev)\s+(?:\{[^}]*\}\s*)?([\w.']+)(.*)", body, re.S)
            rest = m.group(2)
            params = []
            depth = 0; i = 0
            # explicit binders before the type colon
            for g in re.finditer(r"\(([^()]*?):[^()]*\)", rest.split(":=")[0]):
                params += g.group(1).split()
            rec = "termination_by" in rest or re.search(rf"(?<![\w.']){re.escape(short)}(?![\w'])", rest.split(":=", 1)[-1]) is not None
            tree = None
            patterns = bool(re.search(r"^\s*\|", rest, re.M))
            if patterns:
                # a pattern-matching definition recurses when it names itself
                rec = rec or re.search(rf"(?<![\w.']){re.escape(short)}(?![\w'])", rest.split("|", 1)[-1]) is not None
            if ":=" in rest and not rec and not patterns:
                try: tree = D.read(re.sub(r"\s+", " ", rest.split(":=", 1)[1]).strip())
                except Exception: tree = None
            out.setdefault(short, {})[mod] = (params, tree, rec)
    return out

DEFS = parse_defs()

def defs_for(module):
    d = {}
    for short, mods in DEFS.items():
        if module in mods: d[short] = mods[module]
        else:
            for f in FOUNDATION:
                if f in mods: d[short] = mods[f]; break
            else: d[short] = next(iter(mods.values()))
    return d

def dimension(t):
    try:
        e = D.read(t["statement"])
    except Exception:
        return None
    env = {}
    for grp in re.finditer(r"[({\[]([^(){}\[\]]*?):([^(){}\[\]]*)[)}\]]", t["binders"]):
        for nm in grp.group(1).split(): env[nm] = 1
    try:
        d = D.deg(e, env, defs_for(t["module"]))
        return None if d is None else min(d, D.N)
    except Exception:
        return None

# ---------------------------------------------------------------- structural moves
def sides(stmt):
    depth = 0
    for i, c in enumerate(stmt):
        if c in "([{⟨": depth += 1
        elif c in ")]}⟩": depth -= 1
        elif depth == 0 and c == "=" and stmt[i-1:i] not in ("≠", ":", "<", ">", "!") and stmt[i+1:i+2] != ">":
            return stmt[:i].strip(), stmt[i+1:].strip()
    return None

def atoms(x): return [v for k, v in D.tokenize(x) if k in ("id", "num") or v in "+-*/^⊞"]

def structural(t):
    """moves visible in the statement itself: the two sides are one shape re-cut"""
    out = set()
    sd = sides(t["statement"])
    if not sd or "∧" in t["statement"]: return out
    L, R = sd
    al, ar = atoms(L), atoms(R)
    if al and sorted(al) == sorted(ar):
        if al == ar: out.add("regroup")          # same pieces, same order, new brackets
        else: out.add("swap")                    # same pieces, new order
    for one, other in ((L, R), (R, L)):
        if re.fullmatch(r"[\w'.]+", one.strip()) and re.search(r"(?<![\w.])[01](?![\w.])", other):
            out.add("collapse")                  # f a 0 = a : a side of length 0 vanishes
    return out

def main():
    out = []
    for t in T:
        r, own = reach(t)
        hits = Counter()
        for m, rx in RULES:
            k = len(rx.findall(own))
            hits[m] += k
            if k == 0 and rx.search(r):
                hits[m] += 0.5        # touched only through a definition or a cited lemma
        n, counts, lists = rods(t)
        # a literal-only statement with no binders is a single numeric check
        f = frame(t["statement"])
        for m in structural(t): hits[m] += 1
        dim = dimension(t)
        if dim == 0 and re.search(r"\b(?:decide|rfl|norm_num|native_decide)\b", t["proof"]) and not re.search(r"[A-Za-z]", re.sub(r"\b(?:fexp|ftet|fpent|sadd|multAdd|wAdd|tet|pent|Rat|Nat|Int)\b", "", t["statement"])):
            hits["tally"] += 1
        moves = {m: hits[m] for m in MOVES if hits[m] > 0}
        cites = sorted({i.split(".")[-1] for i in IDENT.findall(t["proof"]) if i.split(".")[-1] in thm_by_short and i.split(".")[-1] != t["short"]})
        out.append(dict(name=t["name"], module=t["module"], line=t["line"], doc=t["doc"],
                        statement=t["statement"], proof=t["proof"], frame=f, rods=n, counts=counts,
                        lists=lists, dim=dim, moves=moves, cites=cites))
    json.dump(out, open(HERE / "features.json", "w"), indent=1, ensure_ascii=False)
    print("frames", Counter(o["frame"] for o in out))
    print("dim", Counter(o["dim"] for o in out))
    tot = Counter()
    for o in out:
        for m in o["moves"]: tot[m] += 1
    print("moves", tot)

if __name__ == "__main__":
    main()
