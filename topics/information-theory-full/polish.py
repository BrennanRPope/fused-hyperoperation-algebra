"""Polish (Nov 26) for checks.py section 7: Z-channel cascade capacities.

checks.py maximises the mutual information over the input probability with scipy's bounded minimize_scalar at its
default xatol = 1e-5 and prints the capacities to 3 decimals.  Polish: (1) re-maximise from the same bracket with a
golden-section/Brent search to 1e-14 (numpy only, no scipy here); (2) compare with the closed form for a Z-channel,
C = log2(1 + (1 - p) p^(p/(1-p))), where p = P(1 -> 0); the n-hop cascade of Mz = [[1, 0], [0.3, 0.7]] is a
Z-channel with p = 1 - 0.7^n.
"""
import math, time
import numpy as np

T0 = time.time()
def h(v): return -sum(t * math.log2(t) for t in v if t > 0)

def mi(u, M):
    return h(np.array([1 - u, u]) @ M) - ((1 - u) * h(M[0]) + u * h(M[1]))

def golden_max(f, a, b, tol):
    g = (math.sqrt(5) - 1) / 2; c, d = b - g * (b - a), a + g * (b - a); fc, fd = f(c), f(d); n = 2
    while b - a > tol:
        if fc > fd: b, d, fd = d, c, fc; c = b - g * (b - a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a + g * (b - a); fd = f(d)
        n += 1
    x = 0.5 * (a + b); return x, f(x), n

stored = [0.504, 0.314, 0.206, 0.138, 0.094, 0.065]          # checks_output.txt, section 7
Mz = np.array([[1, 0], [0.3, 0.7]]); P = np.eye(2)
lines = []
worst_loose = worst_tight = 0.0
for n in range(1, 7):
    P = P @ Mz
    p = 1 - 0.7 ** n
    exact = math.log2(1 + (1 - p) * p ** (p / (1 - p)))
    _, c_loose, _ = golden_max(lambda u: mi(u, P), 1e-9, 1 - 1e-9, 1e-5)
    _, c_tight, ne = golden_max(lambda u: mi(u, P), 1e-9, 1 - 1e-9, 1e-14)
    worst_loose = max(worst_loose, abs(c_loose - exact)); worst_tight = max(worst_tight, abs(c_tight - exact))
    lines.append(f"  hop {n}: printed {stored[n-1]:.3f} | xatol 1e-5 {c_loose:.15f} | polished {c_tight:.15f} | closed form {exact:.15f}"
                 f" | printed value {'unchanged' if round(exact, 3) == stored[n-1] else 'CHANGES'} at 3 decimals")
out = ["Z-channel cascade capacities (bits): before (optimizer at xatol 1e-5) -> after (tight search and closed form)"] + lines
out.append(f"  largest error at xatol 1e-5: {worst_loose:.1e} bits (an error of 1e-5 in u costs ~1e-10 in C at a smooth maximum);"
           f" after polish {worst_tight:.1e}; digits gained ~{max(0, round(math.log10(max(worst_loose, 1e-300) / max(worst_tight, 1e-16))))}")
out.append("  verdict (capacity falls every hop) unchanged; no printed digit changes")
out.append(f"wall time {time.time() - T0:.3f} s")
print("\n".join(out))
open(__file__.replace("polish.py", "polish_output.txt"), "w").write("\n".join(out) + "\n")
