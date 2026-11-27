import numpy as np, random
from exactcheck import *
from engine import Tval, strip
its, L, G, rp = item_arrays(12)
# compare hashes to exact K mod p for all words n<=3, gap sum <= 7
for n in (0, 1, 2, 3):
    tot = 0
    for h, dec in enumerate_hashes(n, 7, its, L, G, rp):
        for k in range(len(h)):
            w = dec(k)
            K = strip(Tval(list(w)))
            assert int(h[k]) == ((K % P1) << 32) | (K % P2), (w,)
            tot += 1
    print('n=%d: %d words hashed correctly' % (n, tot))
counts, bad = run({0: 0, 1: 30, 2: 14, 3: 7})
