"""Round 9: Gram's law, Gram blocks and Rosser's rule on the first 5500 zeros.

Facts used (all standard): theta(g_m) = m pi, so S(g_m) = N(g_m) - m - 1 is an integer;
g_m is 'good' iff (-1)^m Z(g_m) > 0 iff S(g_m) is even (all zeros here are simple and on
the line, so Z changes sign at each zero and Z(0) = zeta(1/2) < 0).
The number of zeros in [g_m, g_{m+1}) is 1 + S(g_{m+1}) - S(g_m).
"""
import json
import numpy as np
from stringlib import theta, gram, true_S, String, mean_spacing
from math import pi

gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
OUT = {}

# Gram points covering the zeros
M = int(np.floor(theta(gam[-1]) / pi)) - 1      # last Gram index fully below the last zero
ms = np.arange(-1, M + 1)
g = gram(ms)
N_at = np.searchsorted(gam, g, side='right')
S_at = N_at - ms - 1
assert np.all(np.abs(S_at - np.round(S_at)) < 1e-9)
good = (S_at % 2 == 0)
counts = np.diff(N_at)                 # zeros in [g_m, g_{m+1}) for m = -1..M-1
print(f"Gram points g_-1..g_{M}: {len(ms)}; zeros covered {N_at[-1]}")
hist = {int(c): int(np.sum(counts == c)) for c in np.unique(counts)}
nint = len(counts)
print("Gram intervals by zero count:", hist, " total", nint)
law = float(np.mean(counts == 1))
print(f"Gram's law (exactly one zero) holds in {law:.4f} of intervals; fails {1-law:.4f}")
OUT['hist'] = hist; OUT['n_intervals'] = nint; OUT['law_rate'] = law

# first failures
first_fail = [int(ms[i]) for i in np.nonzero(counts != 1)[0][:6]]
print("first intervals [g_m, g_m+1) breaking the law: m =", first_fail)
OUT['first_fail'] = first_fail
bad_idx = ms[~good]
print("first bad Gram points:", bad_idx[:8].tolist(), " number bad:", int((~good).sum()), "of", len(ms))
OUT['bad_first'] = bad_idx[:10].tolist(); OUT['n_bad'] = int((~good).sum()); OUT['n_gram'] = int(len(ms))

# law rate vs height (bins of 500 intervals)
bins = []
for i in range(0, nint, 500):
    c = counts[i:i + 500]
    bins.append(dict(m_from=int(ms[i]), t_from=float(g[i]), rate=float(np.mean(c == 1)),
                     zero=float(np.mean(c == 0)), two=float(np.mean(c == 2))))
OUT['law_by_height'] = bins
print("law rate by blocks of 500 intervals:", [round(b['rate'], 3) for b in bins])

# S(t) statistics: why it fails
Sg = S_at
print("S at Gram points: values", {int(v): int(np.sum(Sg == v)) for v in np.unique(Sg)})
OUT['S_gram_values'] = {int(v): int(np.sum(Sg == v)) for v in np.unique(Sg)}
tt = np.linspace(15, gam[-1] - 1, 200000)
St = true_S(tt, gam)
print(f"S(t) on [15, {gam[-1]:.0f}]: min {St.min():.3f}, max {St.max():.3f}, rms {St.std():.3f}")
OUT['S_range'] = [float(St.min()), float(St.max()), float(St.std())]
# variance growth vs Selberg (1/2pi^2) log log t, in windows
win = []
for a, b in [(100, 600), (600, 1500), (1500, 3000), (3000, 5900)]:
    sel = (tt >= a) & (tt < b)
    win.append(dict(a=a, b=b, var=float(St[sel].var()),
                    selberg=float(np.log(np.log(np.sqrt(a * b))) / (2 * pi ** 2))))
OUT['S_var'] = win
print("Var S vs (1/2pi^2) log log t:", [(w['a'], round(w['var'], 4), round(w['selberg'], 4)) for w in win])

# Gram blocks: between consecutive good Gram points
gi = np.nonzero(good)[0]
blocks = []
for a, b in zip(gi[:-1], gi[1:]):
    L = b - a
    z = int(N_at[b] - N_at[a])
    pattern = tuple(int(x) for x in counts[a:b])
    blocks.append((int(ms[a]), L, z, pattern))
Ls = {}
for _, L, z, pat in blocks:
    Ls[L] = Ls.get(L, 0) + 1
ross_fail = [bk for bk in blocks if bk[1] != bk[2]]
print("Gram blocks by length:", dict(sorted(Ls.items())), " Rosser failures:", len(ross_fail))
OUT['blocks_by_length'] = {int(k): v for k, v in sorted(Ls.items())}
OUT['rosser_failures'] = len(ross_fail)
types = {}
for _, L, z, pat in blocks:
    if L >= 2:
        types[str(pat)] = types.get(str(pat), 0) + 1
OUT['block_types'] = dict(sorted(types.items(), key=lambda kv: -kv[1])[:12])
print("commonest block types (zeros per interval):", OUT['block_types'])
first_long = {}
for m0, L, z, pat in blocks:
    if L not in first_long:
        first_long[L] = (m0, float(gram(m0)))
OUT['first_block_of_length'] = {int(k): v for k, v in sorted(first_long.items())}
print("first block of each length (m, t):", OUT['first_block_of_length'])

# pattern: are failures clustered? P(next interval fails | this fails)
f = (counts != 1)
pff = float(np.mean(f[1:][f[:-1]]))
print(f"P(fail | previous fails) = {pff:.3f} vs P(fail) = {f.mean():.3f}")
OUT['p_fail_given_prev'] = pff
# zero-interval followed by a two-interval?
z0 = np.nonzero(counts[:-1] == 0)[0]
OUT['zero_then_two'] = float(np.mean(counts[z0 + 1] == 2))
OUT['zero_then_two_or_prev_two'] = float(np.mean((counts[z0 + 1] == 2) | (counts[np.maximum(z0 - 1, 0)] == 2)))
print(f"an empty Gram interval is next to a two-zero interval in {OUT['zero_then_two_or_prev_two']:.3f} of cases")

# Can the prime string tell good from bad Gram points?  good iff S even; predict
# with the nearest even/odd integer to S_model(g_m).
res = {}
for X in [10, 30, 100, 300, 1000, 3000]:
    st = String(X=float(X))
    Sm = st.S(g)
    pred_good = (np.round(Sm) % 2 == 0)
    acc = float(np.mean(pred_good == good))
    # balanced: recall on bad points
    rec_bad = float(np.mean(~pred_good[~good]))
    res[X] = dict(acc=acc, bad_recall=rec_bad, exactS=float(np.mean(np.round(Sm) == Sg)))
    print(f"X={X}: good/bad accuracy {acc:.3f} (always-good baseline {good.mean():.3f}), "
          f"bad recall {rec_bad:.3f}, exact S(g_m) {res[X]['exactS']:.3f}")
OUT['classifier'] = res
OUT['always_good'] = float(good.mean())

json.dump(dict(OUT, S_gram=Sg.tolist(), counts=counts.tolist(), g=g.tolist(), ms=ms.tolist()),
          open('gram.json', 'w'))
