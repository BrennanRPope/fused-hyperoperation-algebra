"""DRAFT: read the polish output (polish_output.txt), compare with results.json, re-read each loop word from the
polished start, and write polished.json. Needs ../tbpy.py with libtb.so built (numpy only)."""
import sys, os, re, json, math
from decimal import Decimal, getcontext
getcontext().prec = 50
HERE = os.path.dirname(os.path.abspath(__file__))
TB = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)   # folder holding tbpy.py + libtb.so
sys.path.insert(0, TB)
import numpy as np, tbpy
E = json.load(open(os.path.join(os.path.dirname(HERE), "results.json")))["E"]
txt = open(os.path.join(HERE, "polish_output.txt")).read()
out = []
for blk in txt.split("ORBIT ")[1:]:
    i = int(blk.split()[0]); o = E[i]; r = dict(id=i)
    r["start_halfres"] = max(abs(float(x)) for x in re.search(r"start_halfres (\S+) (\S+) (\S+)", blk).groups())
    r["double_close"] = float(re.search(r"double_close (\S+) time (\S+)", blk).group(1))
    q = re.search(r"quad0 rc (\S+) it (\d+)(.*?) time (\S+)", blk); r["quad_rc"] = int(q.group(1)); r["quad_it"] = int(q.group(2))
    steps = [float(b) for a, b in re.findall(r"\[(\S+) (\S+)\]", q.group(3))]; r["last_step"] = steps[-1] if steps else None
    z = re.search(r"quad_z (\S+) (\S+) (\S+)", blk).groups(); r["p1"], r["p2"], r["T"] = z
    r["tol_diff"] = float(re.search(r"quad_tol_diff (\S+)", blk).group(1))
    m = re.search(r"quad_close (\S+) steps (\d+) rmin (\S+) total_time (\S+)", blk)
    r["quad_close"] = float(m.group(1)); r["steps"] = int(m.group(2)); r["rmin"] = float(m.group(3)); r["time"] = float(m.group(4))
    err = max(r["tol_diff"], r["last_step"] or 0, 1e-34)
    r["digits"] = int(min(33, math.floor(-math.log10(err))))
    r["dT"] = float(abs(Decimal(z[2]) - Decimal(repr(o["T"]))))
    r["dp"] = float(max(abs(Decimal(z[0]) - Decimal(repr(o["p1"]))), abs(Decimal(z[1]) - Decimal(repr(o["p2"])))))
    p1, p2, T = float(z[0]), float(z[1]), float(z[2])
    t, Xs = tbpy.traj(tbpy.suvakov(p1, p2), T, n=40000)
    w = tbpy.canon(tbpy.word_from_syzygies(tbpy.syzygies(t, Xs)))
    r["word_same"] = (w == tbpy.canon(o["word"])); r["word"] = o["word"]; r["match"] = o["match"]; r["T_old"] = o["T"]
    out.append(r)
out.sort(key=lambda r: r["id"])
json.dump(out, open(os.path.join(HERE, "polished.json"), "w"), indent=1)
ok = [r for r in out if r["quad_rc"] == 0 and r["quad_close"] < 1e-20]
print(f"{len(out)} orbits; {len(ok)} close in quad to < 1e-20; digits min/median/max "
      f"{min(r['digits'] for r in out)}/{sorted(r['digits'] for r in out)[len(out)//2]}/{max(r['digits'] for r in out)}")
print("largest change in T:", max(r["dT"] for r in out), " in p:", max(r["dp"] for r in out))
print("loop words unchanged:", sum(r["word_same"] for r in out), "of", len(out))
print("total time (s):", round(sum(r["time"] for r in out), 1), " slowest:", max(r["time"] for r in out))
for r in out:
    if r["quad_close"] >= 1e-20 or not r["word_same"] or r["digits"] < 25:
        print("NOTE", r["id"], r["match"] or "unlisted", f"T={float(r['T']):.4f} close={r['quad_close']:.1e} digits={r['digits']} word_same={r['word_same']} rmin={r['rmin']}")
