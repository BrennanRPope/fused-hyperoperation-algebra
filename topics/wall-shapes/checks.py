"""DRAFT. Runs every check of the wall-shapes study and writes checks_output.txt."""
import json, subprocess, sys
lines = []
r = subprocess.run([sys.executable, "depth.py"], capture_output=True, text=True); assert r.returncode == 0, r.stderr
d = json.load(open("depth.json"))
names = {"D1_rayleigh_identity_all_equal": "D1 Branden/Rayleigh identity beta*gamma - alpha*delta = k^2(1+kc)",
 "D2_pass": "D2 line test: n=2 real stable iff W>=0; n>=3 only at W=0",
 "D3_pass": "D3 exact Gaussian-integer witnesses in the upper half-plane for n=3,4 and W>0, W<0",
 "D4_pass": "D4 sector law: W>0 missed iff n*alpha<=2pi, W<0 missed iff n*alpha<=pi",
 "D5_pass": "D5 Cayley dictionary: two-spin Lee-Yang = x1 x2 = tanh J",
 "D6_pass": "D6 pure n-body parity: closed-form zeros; on circle only n=2, J>0",
 "D7_pass": "D7 4-body J4 -> infinity: |z| = sqrt2 -+ 1",
 "D8_pass": "D8 pairs + 4-body: line test agrees with the diagonal map at sample points",
 "D9_pass": "D9 pair floor J2 = ln3/8 as J4 -> infinity (N=4)",
 "D10_pass": "D10 all 4-body couplings, no pairs (N=5..10): Lee-Yang fails hot (small J4), holds cold (J4=2)"}
for k, t in names.items(): lines.append(("PASS " if d[k] else "FAIL ") + t)
r = subprocess.run([sys.executable, "walls.py"], capture_output=True, text=True)
w = json.load(open("walls.json"))
ok = len(w) == 25 and all(x["break_type"] in ("positivity", "growth", "number theory", "topology") and x["shapes"] for x in w)
lines.append(("PASS " if ok else "FAIL ") + "W1 catalogue: 25 walls, each with a kind of break and at least one shape")
rk = json.load(open("rank.json"))
r1 = all(len({row["neighbour_maps"] for row in rk[n]["rows"]}) == 1 for n in rk if n.startswith("rank 1"))
r2 = all([row["neighbour_maps"] for row in rk[n]["rows"]] == sorted(row["neighbour_maps"] for row in rk[n]["rows"])
         and rk[n]["rows"][-1]["neighbour_maps"] > 10 * rk[n]["rows"][0]["neighbour_maps"] for n in rk if not n.startswith("rank 1"))
lines.append(("PASS " if r1 else "FAIL ") + "R1 rank 1 systems: neighbour maps constant across levels (finite type)")
lines.append(("PASS " if r2 else "FAIL ") + "R2 rank 2 and 3 systems: neighbour maps grow more than tenfold")
n_pass = sum(l.startswith("PASS") for l in lines)
lines.append(f"{n_pass} of {len(lines)} pass")
open("checks_output.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
