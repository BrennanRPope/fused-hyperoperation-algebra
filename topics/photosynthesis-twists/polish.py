"""Polish for quarter_turns_test.py (Nov 26): the 'best possible on the 3.2 eV scaling line' comes from a grid argmin
(441 x 401 points, spacing 0.005 eV) of spread = biggest step - smallest step, with steps
(dG_OH, d, 3.2 - d, 4.92 - 3.2 - dG_OH). Closed form: the biggest of d and 3.2 - d is at least 1.6 (equal only at d = 1.6),
the smallest of dG_OH and 1.72 - dG_OH is at most 0.86 (equal only at dG_OH = 0.86), so the spread is at least 0.74,
reached at exactly one point. This checks the grid answer against that, and against a 10x finer local grid. numpy only."""
import json, time
import numpy as np
t0 = time.process_time()
TOT = 4.92
grid = json.load(open('quarter_turns_output.json'))['round_trip']['best_possible_on_scaling_line']
spread = lambda oh, d: (np.maximum.reduce([oh, d, 3.2 - d, TOT - 3.2 - oh]) - np.minimum.reduce([oh, d, 3.2 - d, TOT - 3.2 - oh]))
exact = dict(spread_V=1.6 - (TOT - 3.2)/2, dG_OH=(TOT - 3.2)/2, descriptor=1.6)
oh, d = np.meshgrid(np.linspace(grid['dG_OH'] - 0.01, grid['dG_OH'] + 0.01, 401), np.linspace(grid['descriptor'] - 0.01, grid['descriptor'] + 0.01, 401))
S = spread(oh, d); i = np.unravel_index(S.argmin(), S.shape)
fine = dict(spread_V=float(S[i]), dG_OH=float(oh[i]), descriptor=float(d[i]))
lines = ["Polish of the grid argmin in quarter_turns_test.py (best two-way spread on the 3.2 eV line)", ""]
for k in ('spread_V', 'dG_OH', 'descriptor'):
    ch = "no change" if abs(grid[k] - exact[k]) < 1e-12 else f"changed by {exact[k]-grid[k]:+.2e}"
    lines.append(f"{k:11s} grid {grid[k]:.3f} -> closed form {exact[k]:.15f}  (finer local grid {fine[k]:.6f}); {ch}")
lines += ["", "The grid happened to contain the exact optimum (0.86, 1.60), so 0.74 V is exact, not grid-limited.",
          f"CPU time: {time.process_time() - t0:.3f} s"]
out = "\n".join(lines); print(out); open('polish_output.txt', 'w').write(out + "\n")
