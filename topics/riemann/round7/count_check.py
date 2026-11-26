"""Total zero count in a window by the argument principle (S(T) from horizontal continuation),
compared with the count of sign changes on the line, to confirm the off-line deficit."""
import math, numpy as np, json
from rmt import fam, zeros_on_line

def S_of_T(f, T, sig0=4.0, n=4000):
    sig = np.linspace(sig0, 0.5, n)
    v = f(sig + 1j*T)
    ph = np.unwrap(np.angle(v))
    return (ph[-1] - 0.0)/math.pi   # arg f(sig0+iT) ~ 0 (principal, f ~ 1 there)

res = {}
for name, wins in [('x2+5y2', [(500, 520), (900, 920), (1300, 1320)]), ('DH', [(500, 540), (1500, 1540)]),
                   ('x2+y2', [(500, 520)])]:
    th, f = fam(name)
    for (T1, T2) in wins:
        tot = (th(np.array([T2]))[0]-th(np.array([T1]))[0])/math.pi + S_of_T(f, T2) - S_of_T(f, T1)
        g, _ = zeros_on_line(name, T1, T2, step=0.002)
        print(f"{name} [{T1},{T2}]: total by argument principle {tot:.3f}, on the line (fine scan) {len(g)}, off line {tot-len(g):.2f}")
        res[f"{name} {T1}-{T2}"] = [tot, len(g)]
json.dump(res, open('count_check.json', 'w'))
