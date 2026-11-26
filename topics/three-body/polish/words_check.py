"""DRAFT: is each loop word stable? Re-read it from the polished start at the original sampling (3000 per unit time)
and at 4x finer sampling (12000). Writes words_check.json."""
import sys, os, json
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "..")
import numpy as np, tbpy
HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, "polished.json")))
lib = open(os.path.join(HERE, "..", "data", "li-liao-2017-words.md")).read()
out = []
for r in P:
    p1, p2, T = float(r["p1"]), float(r["p2"]), float(r["T"])
    ws = {}
    for k in (3000, 12000):
        t, Xs = tbpy.traj(tbpy.suvakov(p1, p2), T, int(T * k), 1e-13)
        ws[k] = tbpy.word_from_syzygies(tbpy.syzygies(t, Xs), 3)
    out.append(dict(id=r["id"], old=r["word"], w3000=ws[3000], w12000=ws[12000],
                    same3000=tbpy.canon(ws[3000]) == tbpy.canon(r["word"]), stable=tbpy.canon(ws[3000]) == tbpy.canon(ws[12000])))
json.dump(out, open(os.path.join(HERE, "words_check.json"), "w"), indent=1)
print("same as stored at original sampling:", sum(o["same3000"] for o in out), "of", len(out))
print("unchanged at 4x finer sampling:", sum(o["stable"] for o in out), "of", len(out))
for o in out:
    if not (o["same3000"] and o["stable"]): print(o["id"], len(o["old"]), len(o["w3000"]), len(o["w12000"]), o["old"], o["w12000"])
