"""Checks for 'Linear algebra from the circle operation'. Run: python3 checks.py
Every identity the page and the Lean file state is checked here symbolically
(sympy) or, where noted, numerically (numpy)."""
import itertools, math
import numpy as np
import sympy as sp

results = []
def check(name, ok):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)

def Z(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.expand(v) == 0 for v in e)
    return sp.expand(e) == 0

k = sp.Symbol('k')
def mat(n, s):
    return sp.Matrix(n, n, lambda i, j: sp.Symbol(f'{s}{i+1}{j+1}'))
def madd(A, B, w=1):
    return A + B + w * A * B
I2, I3 = sp.eye(2), sp.eye(3)
A2, B2 = mat(2, 'a'), mat(2, 'b')
A3, B3 = mat(3, 'a'), mat(3, 'b')

def esym(A):
    n = A.shape[0]
    return [sum(A.extract(list(S), list(S)).det() for S in itertools.combinations(range(n), j))
            if j else 1 for j in range(n + 1)]

# 1. matrix lift and closure
check("lift 2x2: I+(A∘₁B) = (I+A)(I+B)", Z(I2 + madd(A2, B2) - (I2 + A2) * (I2 + B2)))
check("lift 3x3", Z(I3 + madd(A3, B3) - (I3 + A3) * (I3 + B3)))
check("weighted lift 3x3: I+k(A∘ₖB) = (I+kA)(I+kB)", Z(I3 + k * madd(A3, B3, k) - (I3 + k * A3) * (I3 + k * B3)))
C3 = mat(3, 'c')
check("matrix ∘₁ associative 3x3", Z(madd(madd(A3, B3), C3) - madd(A3, madd(B3, C3))))
check("matrix ∘₁ not commutative in general", sp.expand(madd(A2, B2) - madd(B2, A2)) != sp.zeros(2))

# 2. det(I+A) = 1 + e1 + e2 + e3 (principal minors; counts 1,3,3,1)
e = esym(A3)
check("det(I+A) 3x3 = 1+e1+e2+e3", Z((I3 + A3).det() - sum(e)))
check("det(I+kA) 3x3 = 1+k e1+k^2 e2+k^3 e3", Z((I3 + k * A3).det() - sum(k**j * e[j] for j in range(4))))
check("minor counts are 1,3,3,1", [math.comb(3, j) for j in range(4)] == [1, 3, 3, 1])
e2 = esym(A2)
check("det(I+kA) 2x2 = 1 + k tr + k^2 det", Z((I2 + k * A2).det() - (1 + k * A2.trace() + k**2 * A2.det())))
a, b, c = sp.symbols('a b c')
check("diagonal: det(I+diag)-1 = a∘₁b∘₁c", Z((I3 + sp.diag(a, b, c)).det() - 1 - (lambda x, y: x + y + x * y)((lambda x, y: x + y + x * y)(a, b), c)))

# 3. the dial homomorphism D_k(A ∘ₖ B) = D_k(A) ∘ₖ D_k(B)
def D2(A, w): return A.trace() + w * A.det()
def D3(A, w):
    ee = esym(A); return ee[1] + w * ee[2] + w**2 * ee[3]
wadd = lambda w, x, y: x + y + w * x * y
check("D_k homomorphism 2x2 (all k)", Z(D2(madd(A2, B2, k), k) - wadd(k, D2(A2, k), D2(B2, k))))
check("D_k homomorphism 3x3 (all k)", Z(D3(madd(A3, B3, k), k) - wadd(k, D3(A3, k), D3(B3, k))))
check("k=0 case is trace additivity", Z(D3(A3 + B3, 0) - A3.trace() - B3.trace()))
check("det multiplicative 3x3", Z((A3 * B3).det() - A3.det() * B3.det()))

# 4. 2x2: off-diagonal coupling sets the dial weight
m11, m12, m21, m22, w = sp.symbols('m11 m12 m21 m22 w')
M = sp.Matrix([[m11, m12], [m21, m22]])
check("det(I+M)-1 = m11 ∘_w m22 when m12 m21 = (1-w) m11 m22",
      Z(((I2 + M).det() - 1 - wadd(w, m11, m22)).subs(m12, (1 - w) * m11 * m22 / m21)))

# 5. orthogonal group as Xᵀ ∘ₖ X = 0, circle picture
x, y = sp.symbols('x y')
X = sp.Matrix([[x, -y], [y, x]])
cond = sp.expand((X.T + X + k * X.T * X)[0, 0])
check("2x2 rotation-type: Xᵀ∘ₖX=0 ⇔ 2x + k(x²+y²) = 0", Z(cond - (2 * x + k * (x**2 + y**2))) and Z((X.T + X + k * X.T * X)[0, 1]))
check("circle: (x+1/k)² + y² = 1/k² same equation", Z(sp.expand(k * ((x + 1 / k)**2 + y**2 - 1 / k**2)) - (2 * x + k * (x**2 + y**2))))
Q = (I3 + k * A3)
check("I+kX orthogonal ⇔ Xᵀ ∘ₖ X = 0 (identity)", Z(Q.T * Q - I3 - k * (A3.T + A3 + k * A3.T * A3)))

# 6. Cayley transform = Aᵀ ⊟ A  (numeric, skew A)
rng = np.random.default_rng(1)
S = rng.normal(size=(4, 4)); S = S - S.T
In = np.eye(4)
cay = np.linalg.solve(In + S, In - S)          # (I+S)^{-1}(I-S)
msub = (S.T - S) @ np.linalg.inv(In + S)         # (Aᵀ - A)(I+A)^{-1}, right-division form
# Circle-operation subtraction X ⊟ A with X ⊟ A ∘₁ A = X:  X ⊟ A = (X - A)(I+A)^{-1}
check("Cayley (I+S)^-1(I-S) = I + (Sᵀ ⊟ S)", np.allclose(cay, In + msub))
check("Cayley is orthogonal", np.allclose(cay.T @ cay, In))
check("⊟ undoes ∘₁: (X⊟A)∘₁A = X", np.allclose((msub) + S + msub @ S, S.T))

# 7. eigenvalues: det(I+A)-1 = ∘₁ of eigenvalues ; log det = tr log
G = rng.normal(size=(5, 5)) * 0.12
lam = np.linalg.eigvals(G)
fold = 0
for l in lam: fold = fold + l + fold * l
check("det(I+A)-1 = λ1∘₁…∘₁λn", np.isclose(np.linalg.det(np.eye(5) + G) - 1, fold))
from numpy.linalg import det
def logm_series(A, n=200):
    out, P = np.zeros_like(A), np.eye(len(A))
    for m in range(1, n):
        P = P @ A; out += (-1)**(m + 1) * P / m
    return out
check("log det(I+A) = tr log(I+A)", np.isclose(np.log(det(np.eye(5) + G)), np.trace(logm_series(G))))
def expm(A, n=60):
    out, P = np.eye(len(A)), np.eye(len(A))
    for m in range(1, n):
        P = P @ A / m; out += P
    return out
check("det(exp Y) = exp(tr Y)", np.isclose(det(expm(G)), np.exp(np.trace(G))))

# 8. exterior / symmetric powers and the matrix dial det(I+kA)^(1/k)
H = rng.normal(size=(4, 4)) * 0.2; H = (H + H.T) / 2
mu = np.linalg.eigvalsh(H)
ext = sum(np.prod(c) for j in range(5) for c in itertools.combinations(mu, j))
check("Σ tr Λ^j A = det(I+A)", np.isclose(ext, det(np.eye(4) + H)))
symsum = np.prod(1 / (1 - mu))
check("Σ tr S^j A = 1/det(I-A)", np.isclose(symsum, 1 / det(np.eye(4) - H)))
dial = lambda kk: det(np.eye(4) + kk * H) ** (1 / kk)
check("dial k→0 gives exp(tr A)", np.isclose(dial(1e-7), np.exp(np.trace(H)), rtol=1e-5))
check("dial k=-1 gives bosons 1/det(I-A)", np.isclose(dial(-1), symsum))
# k = -2: real Gaussian  E[exp(xᵀ A x)] = det(I-2A)^(-1/2) — Monte Carlo
xs = rng.normal(size=(400000, 4))
mc = np.mean(np.exp(np.einsum('ni,ij,nj->n', xs, H, xs)))
check("dial k=-2 = real Gaussian E[exp(xᵀAx)] (Monte Carlo, 1%)", abs(mc / dial(-2) - 1) < 0.01)
# alpha-permanent: Σ_S per_α(A_S) = det(I-αA)^(-1/α), α = -k  (exact on 3x3)
def per_alpha(A, al):
    n = A.shape[0]; tot = 0
    for p in itertools.permutations(range(n)):
        seen, cyc = set(), 0
        for i in range(n):
            if i not in seen:
                cyc += 1; j = i
                while j not in seen: seen.add(j); j = p[j]
        tot += al**(n - cyc) * math.prod(A[i, p[i]] for i in range(n))
    return tot
# Vere-Jones: det(I-αA)^(-1/α) = Σ_n (1/n!) Σ_{i in [m]^n} per_α(A[i,i])  (indices may repeat)
H2 = H[:2, :2] * 0.5
for al in (-1.0, 1.0, 2.0, 0.5):
    tot = 1.0
    for n in range(1, 8):
        tot += sum(per_alpha(H2[np.ix_(ix, ix)], al) for ix in itertools.product(range(2), repeat=n)) / math.factorial(n)
    check(f"α-permanent generating function, α={al} (= dial k={-al})", np.isclose(tot, det(np.eye(2) - al * H2) ** (-1 / al), rtol=1e-6))

# 9. KdV: Hirota's weighted cube is det(I+M) for a Cauchy-type matrix
ks = np.array([0.7, 1.3, 2.1]); eta = np.array([0.4, -0.9, 1.7])
Mk = np.array([[2 * math.sqrt(ks[i] * ks[j]) / (ks[i] + ks[j]) * math.exp((eta[i] + eta[j]) / 2)
                for j in range(3)] for i in range(3)])
Aij = lambda i, j: ((ks[i] - ks[j]) / (ks[i] + ks[j]))**2
hirota = sum(math.prod(math.exp(eta[i]) for i in S) * math.prod(Aij(i, j) for i, j in itertools.combinations(S, 2))
             for r in range(4) for S in itertools.combinations(range(3), r))
check("KdV 3-soliton: Hirota τ = det(I+M), M Cauchy-type", np.isclose(hirota, det(np.eye(3) + Mk)))
# each principal minor factors into pair weights (Cauchy determinant)
ok = all(np.isclose(det(Mk[np.ix_(S, S)]),
                    math.prod(math.exp(eta[i]) for i in S) * math.prod(Aij(i, j) for i, j in itertools.combinations(S, 2)))
         for r in range(1, 4) for S in itertools.combinations(range(3), r))
check("Cauchy minors factor into pair weights", ok)
# a generic symmetric matrix does NOT: the triple term is not the product of pair weights
Gs = rng.normal(size=(3, 3)); Gs = Gs @ Gs.T + np.eye(3)
pair = lambda i, j: det(Gs[np.ix_([i, j], [i, j])]) / (Gs[i, i] * Gs[j, j])
trip = det(Gs) / np.prod(np.diag(Gs))
check("generic symmetric 3x3: triple weight ≠ product of pair weights", not np.isclose(trip, pair(0, 1) * pair(0, 2) * pair(1, 2)))
# u = 2 (log τ)_xx solves KdV u_t + 6 u u_x + u_xxx = 0 (exact derivatives, evaluated at sample points)
X_, T_ = sp.symbols('x t')
kk = [sp.Rational(7, 10), sp.Rational(13, 10), sp.Rational(21, 10)]; ee = [sp.Rational(2, 5), sp.Rational(-9, 10), sp.Rational(17, 10)]
th = [kk[i] * X_ - kk[i]**3 * T_ + ee[i] for i in range(3)]
Ms = sp.Matrix(3, 3, lambda i, j: 2 * sp.sqrt(kk[i] * kk[j]) / (kk[i] + kk[j]) * sp.exp((th[i] + th[j]) / 2))
tau_s = (sp.eye(3) + Ms).det()
uu = 2 * sp.diff(sp.log(tau_s), X_, 2)
res = sp.diff(uu, T_) + 6 * uu * sp.diff(uu, X_) + sp.diff(uu, X_, 3)
vals = [abs(sp.N(res.subs({X_: xv, T_: tv}), 30)) for xv, tv in ((sp.Rational(3, 10), sp.Rational(1, 10)), (sp.Rational(-6, 5), sp.Rational(2, 5)), (2, sp.Rational(-3, 10)))]
check("u = 2(log det(I+M))_xx solves KdV (3 solitons, residual < 1e-20)", max(vals) < 1e-20)

print()
n_pass = sum(ok for n, ok in results if not n.startswith("(info)"))
n_all = sum(1 for n, ok in results if not n.startswith("(info)"))
for n, ok in results:
    if n.startswith("(info)"): print(n)
print(f"{n_pass}/{n_all} checks pass")
