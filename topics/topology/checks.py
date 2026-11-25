"""Worked checks: topology tools that fit mult-add  a(+)b = a+b+ab."""
import itertools, math
import numpy as np
import sympy as sp
from math import comb
x, t = sp.symbols('x t')

def trunc(e, n):  # work in Z[x]/(x^(n+1))
    p = sp.Poly(sp.expand(e), x)
    return sp.expand(sum(c*x**m for (m,), c in p.terms() if m <= n))

print("== 1. K-theory: Adams operations are mult-add powers ==")
n = 6
psi = lambda k, e: trunc(e.subs(x, (1+x)**k - 1), n)   # psi^k(L)=L^k, x=L-1
for j, k in [(2,3),(3,5),(2,7)]:
    lhs = psi(j, psi(k, x)); rhs = psi(j*k, x)
    print(f" psi^{j} psi^{k} = psi^{j*k} on K(CP^{n}):", sp.simplify(lhs-rhs) == 0)
for p in [2,3,5,7]:
    d = sp.Poly(trunc((1+x)**p - 1 - x**p, n), x)
    print(f" psi^{p}(x) - x^{p} has every coefficient divisible by {p}:", all(c % p == 0 for c in d.coeffs()))
# Chern character: ch(x)=e^t-1, and ch takes mult-add to + via log
u, v = sp.symbols('u v', real=True); X, Y = sp.exp(u)-1, sp.exp(v)-1     # ch of two line-bundle classes
print(" log(1 + ch(x) (+) ch(y)) = u + v :", sp.simplify(sp.log(1 + X + Y + X*Y) - (u+v)) == 0)

print("\n== Hopf invariant one (Adams-Atiyah 1966) ==")
# CP^2: a=x (2-cell), b=x^2 (4-cell). psi^k(a)=k a + mu_k b
for k in [2,3]:
    print(f" CP^2: psi^{k}(x) =", trunc((1+x)**k-1, 2))
hits = [m for m in range(1, 5001) if (3**m - 1) % (2**m) == 0]
print(" n <= 5000 with 2^n | 3^n-1 :", hits, " -> spheres S^1,S^3,S^7 (dims 2,4,8 algebras)")

print("\n== 2. Todd genus and Riemann-Roch ==")
lg = sp.series(sp.log(1+x), x, 0, 9).removeO()
print(" log(1+x) =", lg)
print(" coefficient of x^(m+1) times (m+1), sign removed = Todd genus of CP^m:",
      [abs(lg.coeff(x, m+1)*(m+1)) for m in range(8)])
ok = True
for N in range(1, 7):
    td = (t/(1-sp.exp(-t)))**(N+1)
    for k in range(-N-2, 8):
        val = sp.series(sp.exp(k*t)*td, t, 0, N+1).removeO().coeff(t, N)
        ok &= (val == sp.binomial(N+k, N))
print(" chi(O(k)) on CP^N = C(N+k, N) for N<=6, k in [-N-2,7]:", ok)

print("\n== 3. Torus homology = the Pascal cube ==")
def torus_betti(dim, s=3):
    # cubical complex of T^dim, s cells per side; cell = (base vertex, set of directions)
    cells = {k: [] for k in range(dim+1)}
    for v in itertools.product(range(s), repeat=dim):
        for S in itertools.product([0,1], repeat=dim):
            cells[sum(S)].append((v, S))
    idx = {k: {c: i for i, c in enumerate(cells[k])} for k in cells}
    ranks = {}
    for k in range(1, dim+1):
        M = np.zeros((len(cells[k-1]), len(cells[k])))
        for j, (v, S) in enumerate(cells[k]):
            dirs = [i for i in range(dim) if S[i]]
            for pos, i in enumerate(dirs):
                S2 = list(S); S2[i] = 0; S2 = tuple(S2)
                v2 = list(v); v2[i] = (v2[i]+1) % s; v2 = tuple(v2)
                sign = (-1)**pos
                M[idx[k-1][(v2, S2)], j] += sign
                M[idx[k-1][(v, S2)], j] -= sign
        ranks[k] = np.linalg.matrix_rank(M)
    b = []
    for k in range(dim+1):
        b.append(len(cells[k]) - ranks.get(k, 0) - ranks.get(k+1, 0))
    return b, [len(cells[k]) for k in range(dim+1)]
for d in [1,2,3,4]:
    b, c = torus_betti(d)
    print(f" T^{d}: cells {c}, Betti {b}, Pascal {[comb(d,k) for k in range(d+1)]}, Euler {sum((-1)**k*bb for k,bb in enumerate(b))}")
a = sp.symbols('a')
for name, cnt in [("P^1", a+2), ("P^2", a**2+3*a+3), ("P^1xP^1", (a+2)**2), ("torus T^2", a**2), ("P^3", sp.expand(((1+a)**4-1)/a))]:
    print(f" {name}: count in mult-add coords {sp.expand(cnt)}  -> Euler characteristic = value at a=0 = {cnt.subs(a,0)}")

print("\n== 4. Winding and Chern numbers by mult-add walks ==")
def qwz(kx, ky, m):
    d = np.array([np.sin(kx), np.sin(ky), m + np.cos(kx) + np.cos(ky)])
    H = np.array([[d[2], d[0]-1j*d[1]], [d[0]+1j*d[1], -d[2]]])
    w, v = np.linalg.eigh(H); return v[:, 0]
def chern(m, N=40):
    ks = np.linspace(0, 2*np.pi, N, endpoint=False)
    u = np.array([[qwz(a_, b_, m) for b_ in ks] for a_ in ks])
    C = 0.0; F = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            u1, u2, u3, u4 = u[i,j], u[(i+1)%N,j], u[(i+1)%N,(j+1)%N], u[i,(j+1)%N]
            U = np.vdot(u1,u2)*np.vdot(u2,u3)*np.vdot(u3,u4)*np.vdot(u4,u1)  # multiply around the plaquette
            F[i,j] = np.angle(U); C += F[i,j]                               # then take the log
    return C/(2*np.pi), F
for m in [-3, -1.5, -0.5, 0.5, 1.5, 3]:
    print(f" QWZ model m={m:+}: Chern number {chern(m)[0]:+.6f}")
# winding of a loop as a mult-add walk
z = np.exp(1j*np.linspace(0, 2*np.pi, 400))
f = (z**3)*(1+0.4*z**-5)        # degree 3 map with a wiggle
steps = f[1:]/f[:-1] - 1          # mult-add steps a_k with f_{k+1} = f_k (1+a_k)
print(" winding of f = sum log(1+a_k) / 2 pi i =", round((np.sum(np.log1p(steps))/(2j*np.pi)).real, 6))

print("\n== 5. Torus knots ==")
def alex(p, q):
    T = sp.symbols('T'); A = sp.symbols('A')
    mp = lambda k: (1+A)**k - 1                      # mult-add power, A = T-1
    D = sp.cancel(mp(p*q)*mp(1)/(mp(p)*mp(q)))
    return sp.expand(D.subs(A, T-1)), D
for p, q in [(2,3),(2,5),(3,4),(3,5)]:
    DT, DA = alex(p, q)
    print(f" T({p},{q}): Alexander = {DT}; degree {sp.degree(DT)} = (p-1)(q-1) = {(p-1)*(q-1)}; value at T=1: {DT.subs(sp.Symbol('T'),1)}")
def gauss_link(c1, c2):
    d1 = np.roll(c1, -1, 0) - c1; d2 = np.roll(c2, -1, 0) - c2
    m1 = c1 + d1/2; m2 = c2 + d2/2; L = 0.0
    for i in range(len(c1)):
        r = m1[i] - m2; L += np.sum(np.einsum('j,kj->k', d1[i], np.cross(d2, r)) / np.linalg.norm(r, axis=1)**3)
    return -L/(4*np.pi)
s = np.linspace(0, 2*np.pi, 900, endpoint=False)
def tor(th, ph, R=2, r=0.7): return np.stack([(R+r*np.cos(th))*np.cos(ph), (R+r*np.cos(th))*np.sin(ph), r*np.sin(th)], 1)
for k in [1, 2, 3]:
    c1 = tor(k*s, s); c2 = tor(k*s + np.pi, s)   # two strands of the (2,2k) torus link
    print(f" two helices wound {k} times around each other on the torus: Gauss linking number {gauss_link(c1, c2):+.4f}")

print("\n== 6. Lefschetz: fixed points of the k-th mult-add power ==")
for k in [2,3,5]:
    for dim in [1,2,3]:
        fixed = (k-1)**dim
        print(f" z->z^{k} on T^{dim}: fixed points {fixed}, Lefschetz number (1-k)^{dim} = {(1-k)**dim}")
