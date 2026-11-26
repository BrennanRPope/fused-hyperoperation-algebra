"""DRAFT: for the orbits whose loop word changed with finer sampling, read the word at 3000, 12000, 48000 and 192000
samples per unit time (in chunks), until it stops changing. Then look each one up in Li and Liao's 695-family list.
Writes words_converge.json."""
import sys, os, json
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "..")
import numpy as np, tbpy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
from checks_lib import lookup
P = {r["id"]: r for r in json.load(open(os.path.join(HERE, "polished.json")))}
W = json.load(open(os.path.join(HERE, "words_check.json")))
ids = [o["id"] for o in W if not o["stable"]] + [int(a) for a in sys.argv[2:]]
def word_at(p1, p2, T, per, chunks=64):
    X = tbpy.suvakov(p1, p2); syz = []; t0 = 0.0; dT = T / chunks
    for c in range(chunks):
        t, Xs = tbpy.traj(X, dT, max(10, int(dT * per)), 1e-13)
        syz += [(t0 + a, b, u) for a, b, u in tbpy.syzygies(t, Xs)]
        X = Xs[-1]; t0 += dT
    return tbpy.word_from_syzygies(syz, 3)
if __name__ == "__main__":
  out = []
  for i in ids:
      r = P[i]; p1, p2, T = float(r["p1"]), float(r["p2"]), float(r["T"])
      ws = {per: word_at(p1, p2, T, per) for per in (3000, 12000, 48000, 192000)}
      c = [tbpy.canon(w) for w in ws.values()]
      m, pw = lookup(ws[192000])
      out.append(dict(id=i, T=T, rmin=r["rmin"], old_word=r["word"], old_match=r["match"], words={str(k): v for k, v in ws.items()},
                      converged=c[-1] == c[-2], new_word=ws[192000], new_match=m, new_power=pw))
      print(i, f"T={T:.3f} rmin={r['rmin']:.4f}", [len(w) for w in ws.values()], "converged" if c[-1] == c[-2] else "NOT converged",
            "old:", r["match"] or "unlisted", "new:", m or "unlisted", pw, flush=True)
  json.dump(out, open(os.path.join(HERE, "words_converge.json"), "w"), indent=1)
