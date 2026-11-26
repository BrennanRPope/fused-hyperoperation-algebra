"""Summarise interlace_check.json into checks (I1-I4) and print them."""
import json, time
import numpy as np
import mpmath as mp
exec(open("interlace_check.py").read().split("if __name__")[0])
R = json.load(open("interlace_check.json")); RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
il = R["interlace"]
both = {k: v for k, v in il.items() if v["n_a"] > 0 and v["n_b"] > 0}
band_both = {k: v for k, v in both.items() if v["in_band"] and v["next_in_band"]}
check("I1 zeros of G_a and G_(a+1) alternate whenever both are on the line in |s| <= 30",
      all(v["ok"] for v in both.values()), f"{len(both)} pairs, {len(band_both)} with both a and a+1 in the band")
none = [k for k, v in il.items() if v["n_b"] == 0]
check("I2 for p = 1.5 the partner a + 1 is already past p and has no zeros on the line, so the pairing has nothing to pair with",
      set(none) == {"p=1.5, a=1", "p=1.5, a=1.25", "p=1.5, a=1.5"}, ", ".join(none))
# I3 collisions: count drops not explained by a zero leaving through s = 30
ev_ = {}
for p, T in R["tracks"].items():
    drops = []; gaps = []
    for i in range(1, len(T["a"])):
        z0, z1 = T["zeros"][i - 1], T["zeros"][i]
        if z0: gaps.append((T["a"][i - 1], float(min(np.diff([0.0] + z0)))))
        if len(z1) < len(z0):
            edge = (len(z0) - len(z1) == 1 and z0[-1] > S - 1.5)
            drops.append(dict(a=T["a"][i], frm=len(z0), to=len(z1), through_edge=bool(edge)))
    first = min([d["a"] for d in drops if not d["through_edge"]], default=None)
    band_gap = min(g for a, g in gaps if a <= float(p) + 1e-9)
    ev_[p] = dict(drops=drops, first_collision_a=first, min_gap_up_to_p=band_gap)
OUT["tracks"] = ev_
print("    ", ev_)
check("I3 no two zeros meet while a <= p; the first meeting comes after a passes p (p = 1.5 at a = 1.7; p = 3 none up to a = 4.5)",
      all((v["first_collision_a"] is None) or v["first_collision_a"] > float(p) for p, v in ev_.items()),
      "; ".join(f"p={p}: smallest gap {v['min_gap_up_to_p']:.2f} up to a = p" for p, v in ev_.items()))
# I4 polish: refine the zeros of three band cases with findroot, compare with the grid, recheck alternation
t0 = time.perf_counter(); worst = 0; ok = True
for p, a in [(3.0, 2.0), (4.0, 2.5), (5.0, 3.5)]:
    zg, _ = zeros(p, a); zp, c = zeros(p, a, polish=True)
    zg1, _ = zeros(p, a + 1); zp1, c1 = zeros(p, a + 1, polish=True)
    worst = max(worst, max(abs(np.array(zg) - np.array(zp))), max(abs(np.array(zg1) - np.array(zp1))))
    res = max(abs(float(ev(c, z))) for z in zp)
    ok &= interlace(zp, zp1) and res < 1e-12
OUT["polish_seconds"] = round(time.perf_counter() - t0, 2)
check("I4 polish: findroot moves each zero by under half a grid step and alternation still holds", ok and worst <= DS/2 + 1e-12,
      f"largest move {worst:.3f}, {OUT['polish_seconds']} s")
json.dump(OUT, open("interlace_summary.json", "w"), indent=1, default=float)
print(f"\n{sum(RES)} of {len(RES)} checks pass")
