"""Asymptotic Fatou (Abel) coordinate of g(x) = e^x - 1 at its parabolic point x = 0.

phi(x) = -2/x + c*log(-x) + sum_{k>=1} a_k x^k   with   phi(g(x)) = phi(x) + 1.
Exact rational arithmetic on truncated power series (known method, Ecalle).
The series diverges (Gevrey-1), so it is only used very close to 0.
"""
from fractions import Fraction as F
from math import factorial

N = 30  # truncation order


def mul(p, q):
    r = [F(0)] * N
    for i, pi in enumerate(p):
        if pi:
            for j in range(N - i):
                r[i + j] += pi * q[j]
    return r


def inv(p):
    r = [F(0)] * N
    r[0] = 1 / p[0]
    for n in range(1, N):
        r[n] = -sum(p[k] * r[n - k] for k in range(1, n + 1)) / p[0]
    return r


def log_series(p):  # log p, p[0] = 1
    dp = [(i + 1) * p[i + 1] for i in range(N - 1)] + [F(0)]
    q = mul(dp, inv(p))
    return [F(0)] + [q[i - 1] / i for i in range(1, N)]


r = [F(1, factorial(n + 1)) for n in range(N)]          # g(x)/x
ir = inv(r)
term_pole = [-2 * ir[i + 1] for i in range(N - 1)] + [F(0)]  # (2/x)(1 - 1/r)
term_log = log_series(r)
powers = []  # x^k (r^k - 1), k = 1..N-2
rk = [F(1)] + [F(0)] * (N - 1)
for k in range(1, N - 1):
    rk = mul(rk, r)
    s = [F(0)] * N
    for i in range(N - k):
        s[i + k] = rk[i] - (1 if i == 0 else 0)
    powers.append(s)

assert term_pole[0] == 1
c = -term_pole[1] / term_log[1]
a = []
for n in range(2, N - 1):
    known = term_pole[n] + c * term_log[n] + sum(a[k] * powers[k][n] for k in range(len(a)))
    a.append(-known / powers[n - 2][n])

if __name__ == "__main__":
    print("c =", c)
    for i, v in enumerate(a[:10], 1):
        print(f"a{i} =", v)
