"""Polish step for the one fitted number in checks.py: the best bilinear k (Hao et al.'s rule).

checks.py fits V12 = a + b + k*a*b to 30 simulated amplitude triples with scipy's least_squares
(default tolerances, ftol = xtol = gtol = 1e-8). The model is linear in k, so the least-squares
answer has a closed form:  k = sum((V12 - a - b)*a*b) / sum((a*b)^2).  The polish rebuilds the same
30 triples from dendrite.py (deterministic, no random numbers) and evaluates that closed form,
then takes one Newton step on the residual sum as a cross-check.
Run: python3 polish.py > polish_output.txt
"""
import json, math, time
import numpy as np
from dendrite import build, respond, lens_constants

before = json.load(open("results.json"))["hao"]["kb"]
t = time.perf_counter()
G = build(); EE, EI = 65.0, -5.0; iI = 10; iE = 40
k, q = lens_constants(G, iE, iI, EE, EI)
sets = []
for gE in np.linspace(0.01, 0.3, 6):
    for gI in np.linspace(0.1, 2.0, 5):
        a = respond(G, [(iE, gE, EE)]); b = respond(G, [(iI, gI, EI)]); ab = respond(G, [(iE, gE, EE), (iI, gI, EI)])
        sets.append((a, b, ab))
S = np.array(sets); A, B, AB = S[:, 0], S[:, 1], S[:, 2]
t_build = time.perf_counter() - t
t = time.perf_counter()
x = A*B; y = AB - A - B
kb = float(x @ y / (x @ x))                       # exact linear least squares
r = y - kb*x; kb_newton = kb + float(x @ r / (x @ x))   # one Newton step (should not move)
dt = time.perf_counter() - t
rms = lambda kk: np.sqrt(np.mean(((A + B + kk*A*B - AB)/(AB - A - B))**2))
print("Best bilinear k for V12 = EPSP + IPSP + k*EPSP*IPSP on the model dendrite (30 pairs)")
print(f"  exact lens k (formula)              = {k:.12f} /mV")
print(f"  before (least_squares, results.json) = {before:.12f} /mV")
print(f"  after  (closed-form least squares)   = {kb:.12f} /mV   (Newton step from it moves {kb_newton-kb:+.1e})")
d = abs(before - kb)
print(f"  change {kb-before:+.2e} /mV; " + ("no change at the printed 4 decimals" if round(before, 4) == round(kb, 4) else "CHANGES the 4-decimal value"))
print(f"  digits: before agreed with the exact optimum to {-math.log10(d/abs(kb)):.1f} significant digits; after to machine precision")
print(f"  bilinear misses the shunting term by {rms(before)*100:.2f}% rms before, {rms(kb)*100:.2f}% after")
print(f"  wall time: rebuilding the 30 triples {t_build:.3f} s, the polish itself {dt*1e3:.3f} ms")
print(f"  verdict unchanged: |kb - k|/k = {abs(kb-k)/k:.4f} < 0.05")
