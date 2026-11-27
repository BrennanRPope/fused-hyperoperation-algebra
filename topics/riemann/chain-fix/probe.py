import numpy as np, sys
sys.path.insert(0, "../knauf-field")
from knauf import chain
def walsh(f):
    """f over 2^n states -> coefficients J_A, f(s) = sum_A J_A prod_{i in A} s_i, s_i = 1-2*bit."""
    a = f.astype(float).copy(); n = int(np.log2(len(a))); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
for n in (6, 8, 10):
    phi = chain(n).astype(float); H = np.log(phi)
    J = walsh(H); A = np.arange(len(J)); size = np.array([bin(a).count("1") for a in A])
    # energy H; Boltzmann weight e^{-beta H}; ferromagnetic means e^{-beta H} favours aligned spins: coefficient of prod s_A in -H is >= 0
    neg = -J
    print(n, "terms", len(J), "| -J_A >= 0 for all nonempty A:", np.all(neg[1:] >= -1e-12), "| odd-size terms present:", np.abs(J[(size % 2) == 1]).max() > 1e-12,
          "| min -J_A even:", neg[(size % 2 == 0) & (size > 0)].min(), "| min -J_A odd:", neg[size % 2 == 1].min())
    comp = phi[::-1]   # flipping all bits: index k -> 2^n-1-k
    print("   phi(all down)=%g phi(all up)=%g ; flip symmetric? %s" % (phi[0], phi[-1], np.allclose(phi, comp)))
