"""DRAFT: confirm every loop word at 48000 samples per unit time (in chunks) against the 12000 reading. Writes words_all.json."""
import sys, os, json
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "..")
import tbpy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from checks_lib import lookup
from words_converge import word_at  # noqa (runs nothing heavy on import? no: guarded below)
P = json.load(open(os.path.join(HERE, "polished.json")))
out = []
for r in P:
    p1, p2, T = float(r["p1"]), float(r["p2"]), float(r["T"])
    a, b = word_at(p1, p2, T, 12000), word_at(p1, p2, T, 48000)
    m, pw = lookup(b)
    out.append(dict(id=r["id"], T=T, word=b, same=tbpy.canon(a) == tbpy.canon(b), old_word=r["word"],
                    word_changed=tbpy.canon(b) != tbpy.canon(r["word"]), old_match=r["match"], new_match=m, power=pw))
json.dump(out, open(os.path.join(HERE, "words_all.json"), "w"), indent=1)
print("12000 and 48000 agree:", sum(o["same"] for o in out), "of", len(out))
print("word differs from the stored one:", [o["id"] for o in out if o["word_changed"]])
print("family changed:", [(o["id"], o["old_match"] or "unlisted", o["new_match"] or "unlisted") for o in out if (o["old_match"] or "") != (o["new_match"] or "")])
print("unlisted now:", sum(not o["new_match"] for o in out))
