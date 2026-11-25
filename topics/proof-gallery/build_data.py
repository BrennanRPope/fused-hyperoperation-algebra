"""build_data.py: pack gallery.json and the rung sheet into site/data.json for the page."""
import json, pathlib
HERE = pathlib.Path(__file__).resolve().parent
G = json.load(open(HERE / "gallery.json"))
R = json.load(open(HERE.parent / "rung-shape" / "page" / "data.json"))

fields = sorted({p["field"] for p in G["proofs"]})
modules = sorted({p["module"] for p in G["proofs"]})
cl = [g["id"] for g in G["groups"]]
P = []
for p in G["proofs"]:
    P.append(dict(n=p["name"], f=fields.index(p["field"]), g=p["group"], m=modules.index(p["module"]), l=p["line"],
                  d=p["doc"], s=p["statement"], b=p["body"], fr=p["frame"], mv=p["moves"], t=p["touched"],
                  c=cl.index(p["cluster"]), h=p["height"], a=p["ambiguous"], lead=p["lead"],
                  km=p["kmoves"], src=p["lead_src"], ls=p["lean_sides"], lm=p["lean_marks"]))
sh = R["sheet"]
keep = [i for i, s in enumerate(sh["s"]) if -2 <= s <= 3.5 and abs(s * 4 - round(s * 4)) < 1e-9]
sheet = dict(s=[sh["s"][i] for i in keep], a=sh["a"], b=[sh["b"][i] for i in keep])
groups = [dict(id=g["id"], body=g["body"], action=g["action"], name=g["name"], sentence=g["sentence"],
               count=g["count"], height=g["height"], ngroups=g["ngroups"], fields=g["fields"], frames=g["frames"])
          for g in G["groups"]]
out = dict(fields=fields, modules=[m.replace("FusedAlgebra.", "FusedAlgebra/").replace(".", "/") + ".lean" for m in modules],
           proofs=P, groups=groups, twins=G["twins"], stats=G["stats"], sheet=sheet)
json.dump(out, open(HERE / "site" / "data.json", "w"), ensure_ascii=False, separators=(",", ":"))
cross = sum(1 for g in groups if len({f for f, _ in g["fields"]}) >= 3)
print((HERE / "site" / "data.json").stat().st_size, "bytes; groups with 3+ fields:", cross, "of", len(groups))
