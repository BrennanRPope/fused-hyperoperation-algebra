"""degree.py: the dimension of a statement, read from its algebra.

The dimension of a circle-operation picture is the number of rods multiplied together
in its biggest piece: a + b is a rod (1), a + b + ab a rectangle (2),
(1+a)(1+b)(1+c) a box (3), (1+a)^n a cube of any size (N). So we parse each
Lean statement and take its polynomial degree in the variables, expanding the
library's own definitions (multAdd a b counts like a + b + ab). A power with a
symbolic exponent, a fold over a list, or a recursive definition run n times
gives N. Anything the parser cannot read gives None, reported as unclear.
"""
import re

N = 99        # "any number of sides"

TOK = re.compile(r"""\s*(?:
  (?P<num>\d+(?:\.\d+)?)|
  (?P<op>⁻¹|::|\+\+|:=|=>|≤|≥|≠|↔|→|∧|∨|¬|∣|⊞|[-+*/^%=<>(),⟨⟩\[\]{}:!|∘·•])|
  (?P<id>[A-Za-z_α-ωΑ-Ωℓ][\w'.₀-₉ᵢⱼₖ]*)|
  (?P<bad>\S))""", re.X)

def tokenize(s):
    out = []
    for m in TOK.finditer(s):
        if m.group("num"): out.append(("num", m.group("num")))
        elif m.group("op"): out.append(("op", m.group("op")))
        elif m.group("id"): out.append(("id", m.group("id")))
        elif m.group("bad"): out.append(("bad", m.group("bad")))
    return out

class P:
    def __init__(self, toks): self.t = toks; self.i = 0
    def peek(self, k=0): return self.t[self.i + k] if self.i + k < len(self.t) else (None, None)
    def next(self): x = self.peek(); self.i += 1; return x
    def eat(self, v):
        if self.peek()[1] == v: self.i += 1; return True
        return False

BIN = {"↔": 1, "→": 2, "∨": 3, "∧": 4, "=": 5, "≠": 5, "≤": 5, "<": 5, "≥": 5, ">": 5, "∣": 5,
       "+": 6, "-": 6, "⊞": 6, "++": 6, "::": 6, "*": 7, "/": 7, "%": 7, "∘": 7, "^": 8}
RIGHT = {"^", "→", "::", "∧", "∨"}
ATOM_STOP = set(BIN) | {"then", "else", "⁻¹", ")", "⟩", "]", ",", "}", ":", "=>", ":=", "|", None}

def parse(p, minprec=0):
    left = parse_app(p)
    while True:
        k, op = p.peek()
        if k == "op" and op in BIN and BIN[op] >= minprec:
            p.next()
            right = parse(p, BIN[op] if op in RIGHT else BIN[op] + 1)
            left = ("bin", op, left, right)
        else:
            break
    return left

def parse_atom_m(p):
    a = parse_atom(p)
    while p.peek()[1] == "⁻¹":
        p.next(); a = ("neg", a)
    while p.peek() == ("bad", ".") and p.peek(1)[0] in ("id", "num"):
        p.next(); k, m = p.next()
        a = ("app", ("id", m), [a]) if k == "id" else a
    return a

def parse_app(p):
    f = parse_atom_m(p)
    args = []
    while True:
        k, v = p.peek()
        if k is None or v in ATOM_STOP or (k == "op" and v not in ("(", "⟨", "[", "!", "¬")):
            break
        args.append(parse_atom_m(p))
    if args: return ("app", f, args)
    return f

def parse_atom(p):
    k, v = p.next()
    if k == "num": return ("num", v)
    if k == "id":
        if v in ("fun", "λ"):
            while p.peek()[1] not in ("=>", None): p.next()
            p.eat("=>"); body = parse(p); return ("fun", body)
        if v == "if":
            parse(p)
            if p.peek()[1] != "then": raise ValueError("if")
            p.next(); a = parse(p)
            if p.peek()[1] != "else": raise ValueError("if")
            p.next(); b = parse(p)
            return ("tuple", a, b)
        return ("id", v)
    if v == "(":
        if p.peek()[1] == ")": p.next(); return ("num", "0")
        e = parse(p)
        if p.peek()[1] == ":":           # type ascription (x : Rat)
            while p.peek()[1] not in (")", None): p.next()
        while p.peek()[1] == ",":        # tuple
            p.next(); e2 = parse(p); e = ("tuple", e, e2)
        p.eat(")"); return e
    if v == "⟨":
        items = [parse(p)]
        while p.eat(","): items.append(parse(p))
        p.eat("⟩"); return ("tuple",) + tuple(items)
    if v == "[":
        items = []
        if not p.eat("]"):
            items.append(parse(p))
            while p.eat(","): items.append(parse(p))
            p.eat("]")
        return ("list", items)
    if v in ("·", "•"):
        return ("id", "·")
    if v in ("-", "¬", "!"):
        return ("neg", parse_atom(p))
    if v in ("∀", "∃"):
        names = []
        while p.peek()[1] not in (",", None):
            kk, vv = p.next()
            if kk == "id": names.append(vv)
            if vv == ":":
                while p.peek()[1] not in (",", None, ")"): p.next()
        p.eat(",")
        return ("quant", v, names, parse(p))
    raise ValueError(f"cannot read {v!r}")

def read(stmt):
    p = P(tokenize(stmt))
    e = parse(p)
    if p.i < len(p.t): raise ValueError("leftover " + str(p.t[p.i:p.i+4]))
    return e

# ------------------------------------------------------------ definitions as degree rules
# The foundation's operations, as functions from argument degrees to a degree.
def _mpow(ds): return N if len(ds) > 1 else ds[0]
KNOWN = {
    "multAdd": lambda ds: ds[0] + ds[1] if len(ds) == 2 else max(ds + [0]),
    "madd": lambda ds: ds[0] + ds[1] if len(ds) == 2 else max(ds + [0]),
    "wAdd": lambda ds: ds[1] + ds[2] + ds[0] * 0 if len(ds) == 3 else max(ds + [0]),
    "sadd": lambda ds: max(ds + [0]),
    "msub": lambda ds: max(ds + [0]), "multAddInv": lambda ds: max(ds + [0]),
    "fexp": lambda ds: N if any(ds) else 0, "ftet": lambda ds: N if any(ds) else 0, "fpent": lambda ds: N if any(ds) else 0, "tet": lambda ds: N if any(ds) else 0, "pent": lambda ds: N if any(ds) else 0,
    "foldMA": lambda ds: N if any(ds) else 0, "prodL": lambda ds: N if any(ds) else 0, "sumL": lambda ds: 1, "lsum": lambda ds: 1, "groups": lambda ds: N,
}
NOTVAR = {"Rat","Nat","Int","K","R","Prop","Bool","True","False","Type","true","false","none","some","Option","List","Sum","Prod","inl","inr","c0","c1","c2","Cone"}
COUNTED = {"pw", "mpow", "wpow", "iter", "repeat", "Nat.repeat", "wOnes", "layers", "stair", "gsum", "psi", "psiT", "wind", "rsum", "dbl"}

def deg(e, env, defs):
    """env: variable -> degree; defs: name -> (params, body_tree) for non-recursive defs"""
    k = e[0]
    if k == "num": return 0
    if k == "opaque": return None
    if k == "id":
        v = e[1]
        if v in env: return env[v]
        if v in defs and not defs[v][0]:
            return deg(defs[v][1], {}, defs) if defs[v][1] is not None else 0
        base = v.split(".")[0]
        if base in env: return env[base]       # p.lx, x.1
        if v in NOTVAR or base in NOTVAR or v in KNOWN or v in COUNTED or v in defs: return 0
        if "." in v and base[:1].isupper(): return 0   # Rat.foo, Nat.bar
        return 1
    if k == "neg": return deg(e[1], env, defs)
    if k == "fun": return deg(e[1], env, defs)
    if k == "tuple": 
        ds = [deg(x, env, defs) for x in e[1:]]
        return None if None in ds else max(ds)
    if k == "list":
        ds = [deg(x, env, defs) for x in e[1]]
        return None if None in ds else (max(ds) if ds else 0)
    if k == "quant":
        env2 = dict(env); 
        for n in e[2]: env2[n] = 1
        return deg(e[3], env2, defs)
    if k == "bin":
        op, a, b = e[1], e[2], e[3]
        da, db = deg(a, env, defs), deg(b, env, defs)
        if da is None or db is None: return None
        if op in ("+", "-", "++", "::", "%", "∣", "=", "≠", "≤", "<", "≥", ">", "↔", "→", "∧", "∨"):
            return max(da, db)
        if op == "⊞": return min(N, da + db)
        if op in ("*", "∘"): return min(N, da + db)
        if op == "/": return max(da, db)
        if op == "^":
            if b[0] == "num": return min(N, da * int(float(b[1])))
            if db == 0: return da          # power by a constant name
            return N
    if k == "app":
        f, args = e[1], e[2]
        ds = [deg(x, env, defs) for x in args]
        if None in ds: return None
        name = f[1] if f[0] == "id" else None
        if name is None: return max(ds + [0])
        short = name.split(".")[-1]
        if short in KNOWN: return KNOWN[short](ds)
        if short in COUNTED or name in COUNTED:
            # last argument is the count: symbolic -> any size, numeral -> that many
            last = args[-1]
            if last[0] == "num": return min(N, max(ds[:-1] + [0]) * int(float(last[1])))
            if ds[-1] == 0 and short in ("pw", "mpow", "wpow"): return max(ds[:-1] + [0])
            return N
        if short in defs:
            params, body, rec = defs[short]
            if rec: return N if any(d > 0 for d in ds) else max(ds + [0])
            if body is None: return max(ds + [0])
            env2 = {pn: (ds[i] if i < len(ds) else 1) for i, pn in enumerate(params)}
            try:
                return deg(body, env2, defs)
            except RecursionError:
                return None
        return max(ds + [0])
    return None
