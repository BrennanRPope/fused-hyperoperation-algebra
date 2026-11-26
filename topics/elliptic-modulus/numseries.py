"""Numeric power-series tools (mpmath, high precision). A series is a list c with c[n] = coefficient of u^n."""
import mpmath as mp
mp.mp.dps = 60
def mul(a, b, N):
    c = [mp.mpf(0)]*N
    for i, ai in enumerate(a[:N]):
        if ai == 0: continue
        for j in range(0, N-i):
            c[i+j] += ai*b[j]
    return c
def integ(a, N): return [mp.mpf(0)] + [a[n]/(n+1) for n in range(N-1)]
def inv(a, N):  # 1/a, a[0] != 0
    b = [mp.mpf(0)]*N; b[0] = 1/a[0]
    for n in range(1, N):
        b[n] = -sum(a[k]*b[n-k] for k in range(1, n+1))/a[0]
    return b
def jacobi(m, N):
    m = mp.mpf(m)
    sn = [0, 1] + [0]*(N-2); cn = [1] + [0]*(N-1); dn = [1] + [0]*(N-1)
    sn = [mp.mpf(v) for v in sn]; cn = [mp.mpf(v) for v in cn]; dn = [mp.mpf(v) for v in dn]
    for it in range(N):
        sn2 = integ(mul(cn, dn, N), N)
        cn2 = [1 - v if i == 0 else -v for i, v in enumerate(integ(mul(sn, dn, N), N))]
        dn2 = [1 - v if i == 0 else -v for i, v in enumerate(integ([m*t for t in mul(sn, cn, N)], N))]
        sn, cn, dn = sn2, cn2, dn2
    return sn, cn, dn
def slope(m, N):
    sn, cn, dn = jacobi(m, N)
    return mul(mul(sn, dn, N), inv(cn, N), N)
def deriv(a, N): return [a[n+1]*(n+1) for n in range(N-1)] + [mp.mpf(0)]
def compose(f, g, N):  # f(g(u)), g[0]=0
    res = [mp.mpf(0)]*N; p = [mp.mpf(1)] + [mp.mpf(0)]*(N-1)
    for n in range(N):
        if n > 0: p = mul(p, g, N)
        if f[n] != 0:
            res = [r + f[n]*q for r, q in zip(res, p)]
    return res
def julia(f, N):
    """Iterative logarithm V: exp(V d/du) u = f, solved order by order (V[n] enters the flow at order n with weight 1)."""
    V = [mp.mpf(0)]*N
    for n in range(3, N):
        V[n] = 0
        V[n] = f[n] - flow(V, 1, n+1)[n]
    return V
def flow(V, t, N):
    term = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)]*(N-2); total = term[:]
    for n in range(1, N):
        term = [t*c/n for c in mul(V, deriv(term, N), N)]
        total = [a + b for a, b in zip(total, term)]
        if all(abs(c) < mp.mpf(10)**(-mp.mp.dps+5) for c in term): break
    return total
