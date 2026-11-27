import numpy as np, sys, json
sys.path.insert(0, "../knauf-field")
from knauf import chain
out = {}
for n in (8, 10, 12):
    phi = chain(n).astype(float); H = np.log(phi); He = 0.5*(H + H[::-1]); Ho = 0.5*(H - H[::-1])
    up = np.array([bin(k).count("1") for k in range(2**n)])
    res = {}
    for beta in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        coef = np.bincount(up, weights=np.exp(-beta*He), minlength=n+1)
        r = np.roots(coef[::-1]); dev = float(np.abs(np.abs(r)-1).max())
        res[beta] = dict(maxdev=dev, roots=[[float(z.real), float(z.imag)] for z in r])
        print(n, beta, "max ||w|-1| =", "%.2e" % dev)
    # how big is the built-in field: correlation of odd part with magnetisation
    m = n - 2*up
    print("   odd part vs magnetisation: corr = %.3f" % np.corrcoef(Ho, -m)[0,1])
    out[n] = res
json.dump(out, open("circle_even.json", "w"))
