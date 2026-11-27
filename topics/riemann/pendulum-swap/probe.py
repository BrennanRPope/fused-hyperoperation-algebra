import numpy as np, sys
sys.path.insert(0, "../knauf-field"); from knauf import chain
def walsh(f):
    a = f.astype(float).copy(); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
for n in (8, 12):
    p = chain(n).astype(float); H = np.log(p); He = 0.5*(H+H[::-1]); J = walsh(He)
    sz = np.array([bin(k).count("1") for k in range(2**n)])
    tot = (J[1:]**2).sum()
    print("n=%d even-half link power by order:" % n, {k: round(float((J[sz == k]**2).sum()/tot), 4) for k in range(2, n+1, 2)})
    # pair links matrix (bit i <-> position)
    P = np.zeros((n, n))
    for k in np.nonzero(sz == 2)[0]:
        i, j = [b for b in range(n) if k >> b & 1]; P[i, j] = P[j, i] = -J[k]
    np.set_printoptions(precision=3, suppress=True, linewidth=150)
    print(P[:8, :8])
