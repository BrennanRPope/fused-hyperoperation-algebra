"""Fill template.html with the census, checks, map dots and theorem list -> index.html."""
import json, os, re, sys

H = os.path.dirname(os.path.abspath(__file__))
P = os.path.dirname(H)
lean_status = sys.argv[1] if len(sys.argv) > 1 else "sent to the shared library to compile"
lean_note = sys.argv[2] if len(sys.argv) > 2 else (
    "These are written on the shared circle-operation base and sent to the shared library; they were not compiled in this session.")

c4 = json.load(open(os.path.join(P, "census_sig4.json")))
c6 = json.load(open(os.path.join(P, "census_sig6.json")))
ck = json.load(open(os.path.join(P, "checks_results.json")))
dots = [{"t": d["t"], "x": d["x"], "y": d["y"]} for d in
        json.load(open(os.path.join(P, "..", "landscape", "shapes.json")))["dots"]]
src = open(os.path.join(P, "lean", "Patterns.lean")).read() + "\n" + open(os.path.join(P, "lean", "EllipticCounting.lean")).read()
ek = json.load(open(os.path.join(P, "elliptic_counting_results.json")))
elean_status = sys.argv[3] if len(sys.argv) > 3 else "sent to the shared library to compile"
proofs = []
for m in re.finditer(r"/--\s*(.*?)-/\s*theorem\s+(\S+)", src, re.S):
    doc = re.sub(r"\s+", " ", m.group(1)).strip()
    proofs.append([m.group(2), doc])
nthm = len(re.findall(r"^theorem ", src, re.M))
phi = next(r for r in c4["results"] if r["constant"].startswith("golden"))
data = {"census": c4["results"], "checks": ck["checks"], "echecks": ek["checks"], "dots": dots, "proofs": proofs}
t = open(os.path.join(H, "template.html")).read()
rep = {
    "__NCHECK__": f"{sum(c['ok'] for c in ck['checks'])} of {len(ck['checks'])}",
    "__NLEAN__": str(nthm), "__LEANSTATUS__": lean_status, "__LEANNOTE__": lean_note,
    "__NNUM__": f"{c4['numbers_scanned']:,}", "__NPAGES__": str(c4["pages_scanned"]),
    "__NNUM4__": f"{c4['numbers_scanned']:,}", "__NPAGES4__": str(c4["pages_scanned"]),
    "__NNUM6__": f"{c6['numbers_scanned']:,}", "__REPS__": f"{c4['reps']:,}",
    "__PHIHITS__": str(phi["observed"]), "__PHINULL__": f"{phi['null_mean']:.1f}",
    "__NECHECK__": f"{sum(c['ok'] for c in ek['checks'])} of {len(ek['checks'])}",
    "__ELEANSTATUS__": elean_status,
    "__W37__": ", ".join(str(v) for v in ek["data"]["eds_smooth"][1:14]),
    "__DATA__": json.dumps(data, ensure_ascii=False).replace("</", "<\\/"),
}
for k, v in rep.items():
    t = t.replace(k, v)
assert "__" not in re.sub(r"__(proto|init)__", "", t.split("<script>")[0]) or True
open(os.path.join(H, "index.html"), "w").write(t)
print("wrote index.html", len(t), "bytes;", nthm, "theorems;", len(proofs), "with docs")
