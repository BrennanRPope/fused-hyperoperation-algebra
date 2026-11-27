"""Which set U of spins should 'flip' define, to put the most energy in the odd (safe) class?
odd share(U) = sum of J_S^2 over links S meeting U an odd number of times / sum over all links. One Hadamard
transform of J^2 gives it for every U at once. This is weighted MAX-CUT on the link (hyper)graph. DRAFT."""
import numpy as np, sys, json
sys.path.insert(0, "../knauf-field"); from knauf import chain
from frames import walsh
out = {}
for n in (8, 10, 12, 14, 16):
    H = np.log(chain(n).astype(float)); J = walsh(H); J[0] = 0; P = J**2; tot = P.sum()
    sig = walsh(P)*len(P)                       # sum_S P_S (-1)^{|S cap U|}
    odd = 0.5*(tot - sig)/tot
    U_all = 2**n - 1; alt = sum(1 << i for i in range(n-1, -1, -2)); alt2 = sum(1 << i for i in range(n-2, -1, -2))
    best = int(np.argmax(odd))
    pairs = {k: P[k] for k in range(len(P)) if bin(k).count("1") == 2}
    out[n] = dict(all=float(odd[U_all]), alternating=float(odd[alt]), alternating_other=float(odd[alt2]), best=float(odd[best]),
                  best_U=format(best, "0%db" % n), pair_share=float(sum(pairs.values())/tot))
    print(n, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[n].items()}, "| gray-frame flip set =", format(alt, "0%db" % n))
json.dump(out, open("bestflip.json", "w"), indent=1)
