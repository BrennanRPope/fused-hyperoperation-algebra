import numpy as np
from rotate import setup, F, start_zeros, newton
for n in (8, 10, 12):
    S = setup(n); z = start_zeros(S, 36.0)
    ths = np.radians(np.linspace(0, 6, 601))
    for t in ths[1:]: z = newton(S, z, t)
    for deg in (0.5, 1, 2, 4, 6): pass
    S2 = setup(n); z2 = start_zeros(S2, 36.0); rows = []
    for k, t in enumerate(ths[1:], 1):
        z2 = newton(S2, z2, t)
        if k in (10, 50, 100, 200, 400, 600):
            f, _ = F(S2, z2, t); rows.append("%.1fdeg max|Re z|=%.2e (resid %.0e)" % (np.degrees(t), np.abs(z2.real).max(), np.abs(f).max()/S2['w'].sum()))
    print(n, " | ".join(rows))
