"""Round 13 (DRAFT): compact summaries of dist13.json (15 MB) and control13.json (4 MB), which are not committed.
Writes dist13_summary.json and control13_summary.json. Regenerate the full files with dist13.py / control13.py."""
import json, math
d = json.load(open("dist13.json")); c = json.load(open("control13.json"))
P = sorted(d["per_pair"], key=lambda p: p[0]); rec = []; best = 0
for p in P:
    if -p[4] > best: best = -p[4]; rec.append(p)
ks = list(range(4, 261))
tail = lambda V: {k: sum(1 for v in V if v >= k) for k in ks}
json.dump(dict(amax=d["amax"], pairs=d["pairs"], not_done=d["not_done"], seconds=d["seconds"],
               row_format="a, b, rho, log2 lower bound, log2 upper bound, nodes, finished, x-word, y-word",
               closest50=sorted(d["per_pair"], key=lambda p: p[4])[:50], record_staircase=rec,
               tail_counts_dist_below_2_to_minus_k=tail([-p[4] for p in d["per_pair"]]),
               max_t_for_dist_ge_a_to_minus_t_from_lower_bounds=max(-p[3] / math.log2(p[0]) for p in d["per_pair"]),
               max_c_for_dist_ge_2_to_minus_ca_from_lower_bounds=max(-p[3] / p[0] for p in d["per_pair"]),
               max_nodes=max(p[5] for p in d["per_pair"])), open("dist13_summary.json", "w"), indent=0)
json.dump(dict(n=c["n"], not_done=c["not_done"], seconds=c["seconds"],
               row_format="seed, window, rho, log2 lower bound, log2 upper bound, nodes, finished",
               closest20=sorted(c["per_point"], key=lambda p: p[4])[:20],
               tail_counts_dist_below_2_to_minus_k=tail([-p[4] for p in c["per_point"]])), open("control13_summary.json", "w"), indent=0)
print("ok", len(rec), "records")
