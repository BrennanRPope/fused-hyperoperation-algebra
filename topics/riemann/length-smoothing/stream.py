"""Level counts of Knauf's chain for long chains without storing 2^n states. DRAFT.
State values at depth n are the first entries of Stern-Brocot pairs: (u, v) -> (u, u+v), (u+v, v), from (1, 1).
Writes levels_n.npz (q, c) for n = 20..NMAX, scratch outside the shared folder."""
import numpy as np, sys, os, time
sys.path.insert(0, "../knauf-field"); from knauf import levels
OUT = "/tmp/knauf_levels"; os.makedirs(OUT, exist_ok=True)
def expand(u, v): return np.concatenate([u, u+v]), np.concatenate([u+v, v])
if __name__ == "__main__":
    NMAX = int(sys.argv[1]); B = 22
    u = np.array([1], np.int64); v = np.array([1], np.int64)
    for n in range(1, 21):
        u, v = expand(u, v)
        if n <= 18: assert np.array_equal(np.bincount(u), np.bincount(np.array(levels(n)[0]).astype(np.int64).repeat(levels(n)[1].astype(int)))) if n < 6 else True
    q20, c20 = levels(20); qq = np.nonzero(np.bincount(u))[0]; assert np.array_equal(qq, q20.astype(np.int64)) and np.array_equal(np.bincount(u)[qq], c20.astype(np.int64))
    hist = {n: None for n in range(21, NMAX+1)}; t = time.time()
    blk = 2**12
    for i in range(0, u.size, blk):
        uu, vv = u[i:i+blk], v[i:i+blk]
        for n in range(21, NMAX+1):
            uu, vv = expand(uu, vv); h = np.bincount(uu)
            if hist[n] is None or hist[n].size < h.size: hh = np.zeros(h.size, np.int64); hh[:0 if hist[n] is None else hist[n].size] = 0 if hist[n] is None else hist[n]; hist[n] = hh
            hist[n][:h.size] += h
        if i % (blk*32) == 0: print("block", i//blk, "of", u.size//blk, "%.0fs" % (time.time()-t), flush=True)
    for n, h in hist.items():
        q = np.nonzero(h)[0]; np.savez(os.path.join(OUT, "levels_%d.npz" % n), q=q, c=h[q]); print(n, "states", h.sum(), "distinct", q.size)
