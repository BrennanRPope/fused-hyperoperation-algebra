import numpy as np, sys
sys.path.insert(0, "../knauf-field")
from knauf import chain
for n in (12, 16, 20):
    phi = chain(n).astype(np.int64); P = phi*phi[::-1]
    c = np.bincount(P); Pv = np.nonzero(c)[0]; cv = c[Pv]
    print(n, "min P", Pv[:12], "counts", cv[:12])
    for b in (3.0, 4.0, 6.0):
        print("   beta", b, "Z_even =", (cv*Pv.astype(float)**(-b/2)).sum())
