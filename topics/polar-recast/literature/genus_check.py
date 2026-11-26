"""Check: the plane-dial law F(x,y)=(x+y+kxy)/(1-dxy) is Hirzebruch's two-parameter
Todd genus formal group law. Its genus of CP^n, Gr(2,4) and the flag variety Fl(3)
should be q-numbers / q-binomials / q-factorials with q = beta/alpha, where
alpha+beta = k, alpha*beta = d (roots of t^2 - k t + d).
Genus by Atiyah-Bott localization: phi(M) = sum over torus fixed points of prod 1/f(w),
f = exponential (inverse of logarithm) of F."""
import itertools, cmath, random

def f_exp(t, a, b):
    # logarithm g'(x)=1/((1+a x)(1+b x)); inverse f(t)=(e^{(a-b)t}-1)/(a-b e^{(a-b)t})
    if abs(a-b) < 1e-12:
        return t/(1-a*t)
    E = cmath.exp((a-b)*t)
    return (E-1)/(a-b*E)

def genus(fixed_points_weights, a, b):
    s = 0
    for ws in fixed_points_weights:
        p = 1
        for w in ws: p /= f_exp(w, a, b)
        s += p
    return s

def h(n, a, b):  # complete homogeneous: sum a^i b^(n-i)
    return sum(a**i * b**(n-i) for i in range(n+1))

random.seed(1)
t = [random.uniform(-1,1) for _ in range(6)]
out = []
for (k, d, label) in [(1.0, 0.0, "core law d=0"), (0.0, 1.0, "angle k^2<4d"),
                      (3.0, 1.0, "radius k^2>4d"), (2.0, 1.0, "straight k^2=4d"), (1.0, 2.0, "angle")]:
    disc = cmath.sqrt(k*k-4*d)
    a, b = (k+disc)/2, (k-disc)/2   # x^2 coefficient: (1+ax)(1+bx)=1+kx+dx^2
    q = b/a if abs(a) > 1e-12 else None
    # CP^n fixed points e_i, weights t_j - t_i
    res = [label, f"k={k} d={d} |q|={abs(q):.4f}" if q is not None else "q=0"]
    for n in (2, 3):
        fp = [[t[j]-t[i] for j in range(n+1) if j != i] for i in range(n+1)]
        g = genus(fp, a, b); pred = (-1)**n * h(n, a, b)
        res.append(f"CP{n}: {g.real:+.6f} vs {pred.real:+.6f}")
    # Gr(2,4): fixed points = 2-subsets S, weights t_j - t_i, i in S, j not in S
    fp = []
    for S in itertools.combinations(range(4), 2):
        fp.append([t[j]-t[i] for i in S for j in range(4) if j not in S])
    g = genus(fp, a, b)
    # homogeneous q-binomial [4 choose 2]: sum over 2-subsets of a^(4-area) b^(area)
    pred = sum(a**(4-ar) * b**ar for ar in [0,1,2,2,3,4])
    res.append(f"Gr(2,4): {g.real:+.6f} vs {pred.real:+.6f}")
    # Fl(3): fixed points = permutations, weights t_{s(j)} - t_{s(i)} for i<j
    fp = []
    for s in itertools.permutations(range(3)):
        fp.append([t[s[j]]-t[s[i]] for i in range(3) for j in range(i+1,3)])
    g = genus(fp, a, b)
    inv = lambda s: sum(1 for i in range(3) for j in range(i+1,3) if s[i] > s[j])
    pred = (-1)**3 * sum(a**(3-inv(s)) * b**inv(s) for s in itertools.permutations(range(3)))
    res.append(f"Fl(3): {g.real:+.6f} vs {pred.real:+.6f}")
    out.append(" | ".join(res))
print("\n".join(out))
