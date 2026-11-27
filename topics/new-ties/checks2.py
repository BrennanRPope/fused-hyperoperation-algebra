"""Round-2 checks. Run after round2.py."""
import json, numpy as np
R=json.load(open("round2.json")); R1=json.load(open("union_ifs.json")); out=[]
ok=lambda n,c,i="": out.append((n,bool(c),i))
for name,r in R["R1"].items():
    d=[x[1] for x in r["sdim_S"]]
    ok(f"L1 lower bounds increase and stay below D*, {name}", np.all(np.diff(d)>0) and d[-1]<=r["D_trace"]+1e-9, f"s-dim(S_10)={d[-1]:.5f}, D*={r['D_trace']:.5f}")
    ok(f"L2 lower bound within 0.001 of D* at L=10, {name}", r["D_trace"]-d[-1]<1e-3)
    f=r["free_check"]; ok(f"L3 generator words are all distinct maps (free), {name}", f["words"]==f["distinct"], f"{f['distinct']} of {f['words']}")
    ok(f"L4 squeeze: lower <= D* <= upper D_11, {name}", d[-1]<=r["D_trace"]<=R1[name]["Dn"][-1])
r=R["R2"]; ok("G1 flips (negative slope) keep the count", r["match"]); ok("G2 flips: D_n above D*", r["Dn"][-1]>=r["D_trace"])
r=R["R3"]; ok("G3 plane: count = shared-centre count", r["match"]); ok("G4 plane: box count nearer D* than Moran", abs(r["box"]-r["D_trace"])<abs(r["box"]-r["D_moran"]), f"box {r['box']:.4f} D* {r['D_trace']:.4f} Moran {r['D_moran']:.4f}")
for lab,v in R["R4"].items():
    if lab.startswith("collision"): ok("G5 Möbius collision example is detected (count falls short)", (not v["match"]) and v["first_collision_length"]==4)
    else: ok(f"G6 Möbius plane laws: count = shared-centre count, {lab}", v["match"])
r=R["R5"]
ok("S1 learning: distinct strengths = every other Fibonacci number", r["counts"][:6]==[3,8,21,55,144,377])
ok("S2 learning: bits kept on average <= most it can keep <= log2 3", r["hrw_bits_increment"]<=r["bits_per_trial_max"]<=r["bits_if_all_differ"])
ok("S3 learning: growth per trial = golden ratio squared", abs(r["growth"]-((1+5**0.5)/2)**2)<1e-9)
p=sum(c for _,c,_ in out)
with open("checks2_output.txt","w") as fh:
    for n,c,i in out:
        l=f"{'PASS' if c else 'FAIL'}  {n}  {i}"; print(l); fh.write(l+"\n")
    fh.write(f"{p}/{len(out)} pass\n")
print(f"{p}/{len(out)} pass")
