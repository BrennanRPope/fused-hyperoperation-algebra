"""Dial sweep for linear algebra from mult-add. Run: python3 checks_dial.py"""
import itertools, math, cmath
import numpy as np
import sympy as sp

res = []
def check(n, ok):
    res.append((n, bool(ok))); print(("PASS " if ok else "FAIL ") + n)
Z = lambda e: all(sp.expand(v) == 0 for v in (e if isinstance(e, sp.MatrixBase) else [e]))
rng = np.random.default_rng(7)
k, a, b, q, kap, t = sp.symbols('k a b q kappa t')
A2 = sp.Matrix(2, 2, sp.symbols('a11 a12 a21 a22')); B2 = sp.Matrix(2, 2, sp.symbols('b11 b12 b21 b22'))
A3 = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'a{i+1}{j+1}')); I2, I3 = sp.eye(2), sp.eye(3)
wm = lambda X, Y, w: X + Y + w * X * Y

# 1. real dial: one law in other units for every k != 0
check("k(A ⊞k B) = (kA) ⊞ (kB)", Z(k * wm(A2, B2, k) - wm(k * A2, k * B2, 1)))
# det(I+kA) is the reversed characteristic polynomial; walls at k = -1/λ
p3 = (t * I3 - A3).det()
check("det(I+kA) = (-k)^3 p_A(-1/k)", Z((I3 + k * A3).det() - (-k)**3 * p3.subs(t, -1 / k)))
check("k = -1: det(I-A) = p_A(1), Fredholm λ = -k", Z((I3 - A3).det() - p3.subs(t, 1)))
lam = np.array([2.0, 0.5, -1.0]); P = rng.normal(size=(3, 3)); A = P @ np.diag(lam) @ np.linalg.inv(P)
zeros = [-1 / l for l in lam]
check("walls: det(I+kA)=0 exactly at k = -1/λ", all(abs(np.linalg.det(np.eye(3) + z * A)) < 1e-9 for z in zeros))
# k -> infinity: det(I+kA)/k^n -> det A ; k -> 0: (det(I+kA)-1)/k -> tr A
check("k→∞: det(I+kA)/k^3 → det A", np.isclose(np.linalg.det(np.eye(3) + 1e6 * A) / 1e18, np.linalg.det(A), rtol=1e-4))
check("k→0: (det(I+kA)-1)/k → tr A", np.isclose((np.linalg.det(np.eye(3) + 1e-7 * A) - 1) / 1e-7, np.trace(A), rtol=1e-5))
# shares of slabs/beams/corner in D_k for positive eigenvalues
lp = np.array([0.8, 0.6, 0.5])
def shares(kk):
    e1, e2, e3 = lp.sum(), lp[0]*lp[1]+lp[0]*lp[2]+lp[1]*lp[2], lp.prod()
    v = np.array([e1, kk * e2, kk * kk * e3]); return v / v.sum()
check("slabs dominate at small k, corner at large k", shares(0.01)[0] > .95 and shares(1e3)[2] > .95)

# 2. rotations and Cayley at weight k
x, y = sp.symbols('x y', real=True)
kr, ki = sp.symbols('kr ki', real=True); m2 = kr**2 + ki**2
z = x + sp.I * y; kc = kr + sp.I * ki
lhs = sp.expand((1 + kc * z) * sp.conjugate(1 + kc * z) - 1)
rhs = sp.expand(m2 * ((x + kr / m2)**2 + (y - ki / m2)**2 - 1 / m2))
check("|1 + k z| = 1 ⇔ circle through 0 centred at -1/k (complex k)", sp.simplify(lhs - rhs) == 0)
s = sp.Symbol('s'); S = sp.Matrix([[0, -s], [s, 0]])
Ck = (I2 + k * S).inv() * (I2 - k * S)
check("weighted Cayley = Cayley of k·s (angle 2·atan(ks))", Z(sp.simplify(Ck - ((I2 + S.subs(s, k * s)).inv() * (I2 - S.subs(s, k * s))))))
H = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); H = (H + H.conj().T) / 2
Ui = np.linalg.solve(np.eye(4) + 1j * H, np.eye(4) - 1j * H)
check("k = i: Hermitian H gives unitary (I+iH)^-1(I-iH)", np.allclose(Ui.conj().T @ Ui, np.eye(4)))
check("k = i: wall -1/i = i is never an eigenvalue of Hermitian H", np.min(np.abs(np.linalg.eigvalsh(H) - 1j)) > 0.99)
X = rng.normal(size=(3, 3)) * .3
for kk in (-1.0, 0.5, 3.0):
    inv = -X @ np.linalg.inv(np.eye(3) + kk * X)
    ok = np.allclose(X + inv + kk * X @ inv, 0)
    check(f"weighted inverse at k={kk}: -X(I+kX)^-1", ok)
check("k = -1: I - X orthogonal ⇔ X + Xᵀ - XᵀX = 0 (circle centred +1)",
      Z((I3 - A3).T * (I3 - A3) - I3 - (-(A3.T + A3) + A3.T * A3)))

# 3. a matrix on the dial: U with U² = kU - qI gives the whole landscape plane
U = sp.Matrix([[0, -q], [1, k]])
check("companion U² = kU - qI", Z(U * U - (k * U - q * I2)))
check("det(I + aU) = 1 + k a + q a²  (the landscape wall polynomial)", Z((I2 + a * U).det() - (1 + k * a + q * a**2)))
check("(I+aU)(I+bU) = (1 - q ab) I + (a + b + k ab) U", Z((I2 + a * U) * (I2 + b * U) - ((1 - q * a * b) * I2 + (a + b + k * a * b) * U)))
l1, l2 = sp.symbols('l1 l2')
check("walls = -1/eigenvalues of U: det(I+aU) = (1+a l1)(1+a l2)", Z(((I2 + a * U).det() - (1 + a * l1) * (1 + a * l2)).subs({k: l1 + l2, q: l1 * l2})))
check("countries = discriminant of U: k² - 4q = (l1 - l2)²", Z((k**2 - 4 * q - (l1 - l2)**2).subs({k: l1 + l2, q: l1 * l2})))
check("or line k + q = -1 ⇔ wall at a = 1 ⇔ U has eigenvalue -1", Z((I2 + U).det() - (1 + k + q)))
Uk = sp.Matrix([[0, kap], [1, 0]])
check("pure unit U² = κI: κ = -1 rotation, 0 shear, +1 boost", Z(Uk * Uk - kap * I2))
check("projected law at κ: (a+b)/(1+κab)", Z((I2 + a * Uk) * (I2 + b * Uk) - ((1 + kap * a * b) * I2 + (a + b) * Uk)))
th1, th2 = 0.4, 0.9
check("κ = -1 gives tangent addition", np.isclose((math.tan(th1) + math.tan(th2)) / (1 - math.tan(th1) * math.tan(th2)), math.tan(th1 + th2)))
check("κ = +1 gives velocity addition (tanh), walls at ±1 = light speed", np.isclose((math.tanh(th1) + math.tanh(th2)) / (1 + math.tanh(th1) * math.tanh(th2)), math.tanh(th1 + th2)))
check("κ = +1: det(I + aU) = (1-a)(1+a)", Z((I2 + a * Uk).det().subs(kap, 1) - (1 - a) * (1 + a)))
# shape number s = k²/(k²-4q) invariant under units a -> c a (k -> c k, q -> c² q)
c = sp.Symbol('c')
check("shape s = k²/(k²-4q) unchanged by units", Z(sp.simplify((c*k)**2 / ((c*k)**2 - 4 * c**2 * q) - k**2 / (k**2 - 4 * q))))

# 4. matrix dial K inside the law: A ⊞_K B = A + B + AKB ; Sylvester ; rank-one = addition
K2 = sp.Matrix(2, 2, sp.symbols('k11 k12 k21 k22'))
check("(I + AK)(I + BK) = I + (A + B + AKB)K", Z((I2 + A2 * K2) * (I2 + B2 * K2) - (I2 + (A2 + B2 + A2 * K2 * B2) * K2)))
check("Sylvester 2x2: det(I + AK) = det(I + KA)", Z((I2 + A2 * K2).det() - (I2 + K2 * A2).det()))
Ar = rng.normal(size=(3, 5)); Kr = rng.normal(size=(5, 3))
check("Sylvester rectangular 3x5 / 5x3: cube in dim 3 = cube in dim 5", np.isclose(np.linalg.det(np.eye(3) + Ar @ Kr), np.linalg.det(np.eye(5) + Kr @ Ar)))
u = sp.Matrix(sp.symbols('u1 u2 u3')); v = sp.Matrix(sp.symbols('v1 v2 v3'))
check("rank-one dial: det(I + A u vᵀ) = 1 + vᵀ A u (no corner terms)", Z((I3 + A3 * u * v.T).det() - (1 + (v.T * A3 * u)[0])))

print(); n = sum(o for _, o in res if _ != "info"); m = sum(1 for _, o in res if _ != "info")
print(f"{n}/{m} checks pass")
