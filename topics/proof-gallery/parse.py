"""parse.py: split every Lean file in the shared library into theorems.

Reads ../lean-library/FusedAlgebra/**/*.lean (never writes there). For each
`theorem` it records: full name (with namespaces, matching theorems.tsv),
module, doc comment, binders, statement, proof text, and source line.
Output: theorems.json
"""
import json, re, pathlib, sys

LIB = pathlib.Path(__file__).resolve().parent.parent / "lean-library"
SRC = LIB / "FusedAlgebra"

TOP = re.compile(r"^(private |protected |noncomputable )*(theorem|def|abbrev|structure|inductive|instance|example|namespace|section|end|open|variable|set_option|@\[|/--|/-!|#|termination_by|decreasing_by)\b|^(/--|/-!|/-|@\[)")

def strip_comments(s):
    s = re.sub(r"/-.*?-/", " ", s, flags=re.S)
    s = re.sub(r"--[^\n]*", " ", s)
    return s

def module_name(p):
    rel = p.relative_to(LIB).with_suffix("")
    return ".".join(rel.parts)

def split_top(body):
    """split a declaration at the top-level ':=' (outside brackets) or at '|' patterns"""
    depth = 0
    i = 0
    while i < len(body):
        c = body[i]
        if c in "([{⟨": depth += 1
        elif c in ")]}⟩": depth -= 1
        elif depth == 0 and body.startswith(":=", i):
            return body[:i], body[i+2:]
        elif depth == 0 and c == "\n" and body[i+1:].lstrip(" ").startswith("|") :
            return body[:i], body[i:]
        i += 1
    return body, ""

def split_sig(sig):
    """binders vs statement: statement starts at the first top-level ':'"""
    depth = 0
    for i, c in enumerate(sig):
        if c in "([{⟨": depth += 1
        elif c in ")]}⟩": depth -= 1
        elif depth == 0 and c == ":" and not sig.startswith(":=", i):
            return sig[:i], sig[i+1:]
    return sig, ""

def parse_file(p):
    lines = p.read_text().split("\n")
    out = []
    ns = []          # stack of (kind, name)
    i = 0
    doc = None
    while i < len(lines):
        ln = lines[i]
        m = re.match(r"^namespace\s+(\S+)", ln)
        if m: ns.append(("ns", m.group(1))); i += 1; continue
        m = re.match(r"^section\b\s*(\S*)", ln)
        if m: ns.append(("sec", m.group(1))); i += 1; continue
        m = re.match(r"^end\b\s*(\S*)", ln)
        if m and ns: ns.pop(); i += 1; continue
        if ln.startswith("/--"):
            j = i; buf = []
            while True:
                buf.append(lines[j])
                if "-/" in lines[j]: break
                j += 1
            doc = re.sub(r"\s+", " ", " ".join(buf).replace("/--", "").replace("-/", "")).strip()
            i = j + 1; continue
        m = re.match(r"^(?:private |protected )?theorem\s+(\S+)", ln)
        if m:
            name = m.group(1)
            j = i + 1
            while j < len(lines) and not (lines[j] and not lines[j][0].isspace() and TOP.match(lines[j])):
                j += 1
            raw = "\n".join(lines[i:j])
            body = strip_comments(raw)
            body = re.sub(r"^(?:private |protected )?theorem\s+\S+", "", body, count=1)
            sig, proof = split_top(body)
            binders, stmt = split_sig(sig)
            prefix = ".".join(n for k, n in ns if k == "ns" and n)
            full = (prefix + "." if prefix else "") + name
            out.append(dict(name=full, short=name, module=module_name(p), line=i + 1,
                            doc=doc or "", binders=binders.strip(), statement=re.sub(r"\s+", " ", stmt).strip(),
                            proof=proof.strip()))
            doc = None
            i = j; continue
        if ln and not ln[0].isspace() and not ln.startswith("/--"):
            if not ln.startswith("@["): doc = None
        i += 1
    return out

def main():
    allt = []
    for p in sorted(SRC.rglob("*.lean")):
        allt += parse_file(p)
    tsv = {}
    for row in (LIB / "theorems.tsv").read_text().strip().split("\n"):
        mod, name, _ = row.split("\t")
        tsv[name] = mod
    names = {t["name"] for t in allt}
    miss = set(tsv) - names
    extra = names - set(tsv)
    print(f"parsed {len(allt)} theorems; tsv has {len(tsv)}; missing {len(miss)}; extra {len(extra)}")
    if miss: print("missing:", sorted(miss)[:20])
    if extra: print("extra:", sorted(extra)[:20])
    json.dump(allt, open(pathlib.Path(__file__).parent / "theorems.json", "w"), indent=1, ensure_ascii=False)

if __name__ == "__main__":
    main()
