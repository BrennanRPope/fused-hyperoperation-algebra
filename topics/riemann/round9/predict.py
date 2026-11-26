"""Round 9: predicting the next zero on held-out zeros (501-5500).

Point predictors (error in mean spacings):
  G   Gram midpoint g_{n-3/2} (the smooth law; Gram's law in point form)
  L   last zero plus the mean gap: gamma_{n-1} + 2 pi / log(gamma_{n-1}/2pi)
  P   the fitted prime string (settings fitted on zeros 1-500)
  PL  prime string anchored on the last true zero: F(E) = F(gamma_{n-1}) + 1
Interval predictors:
  Gram's law   gamma_n in [g_{n-2}, g_{n-1})
  Rosser       gamma_n in the Gram block that contains g_{n-3/2}
"""
import json
import numpy as np
from math import pi
from stringlib import String, theta, gram, mean_spacing, score

gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
TE = np.arange(500, NZ)
best = dict(taper='gauss', X=200.0)   # the regulator used for the sweeps (probe.py)
OUT = {'setting': best}
n = np.arange(1, NZ + 1)


def stats(pred, idx=TE):
    e = (pred[idx] - gam[idx]) / mean_spacing(gam[idx])
    return dict(rms=float(np.sqrt(np.nanmean(e ** 2))), mae=float(np.nanmean(np.abs(e))),
                near=float(np.nanmean(np.abs(e) < 0.5)),
                gram=float(np.nanmean(np.floor(theta(pred[idx]) / pi) == np.floor(theta(gam[idx]) / pi))))


G = gram(n - 1.5)
Lp = np.concatenate([[np.nan], gam[:-1] + mean_spacing(gam[:-1])])
st = String(X=best['X'], taper=best['taper'])
P = st.levels(NZ)
# anchored
lo = gam[:-1].copy(); hi = gam[:-1] + 3 * mean_spacing(gam[:-1])
target = st.F(gam[:-1]) + 1.0
for _ in range(48):
    mid = 0.5 * (lo + hi)
    neg = st.F(mid) < target
    lo = np.where(neg, mid, lo); hi = np.where(neg, hi, mid)
PL = np.concatenate([[np.nan], 0.5 * (lo + hi)])
# larger cutoff string, same taper (more orbits)
st2 = String(X=3000.0, taper='sharp')
P2 = st2.levels(NZ)

for lab, pr in [('G', G), ('L', Lp), ('P', P), ('PL', PL), ('P_X3000', P2)]:
    OUT[lab] = stats(pr)
    print(lab, OUT[lab])

# interval rules
gi = np.floor(theta(gam) / pi)            # Gram interval index of each zero
OUT['gram_law_index_form'] = float(np.mean(gi[TE] == (n[TE] - 2)))
# Rosser blocks
g_json = json.load(open('gram.json'))
ms = np.array(g_json['ms']); S_gram = np.array(g_json['S_gram'])
good = (S_gram % 2 == 0)
goodm = ms[good]
blk_of = lambda m: np.searchsorted(goodm, m, side='right') - 1   # block index of interval [g_m, g_m+1)
blk_true = blk_of(gi[TE])
blk_pred = blk_of(np.floor(n[TE] - 1.5))
OUT['rosser_block_hit'] = float(np.mean(blk_true == blk_pred))
for lab, pr in [('P', P), ('PL', PL), ('P_X3000', P2)]:
    OUT[f'{lab}_block_hit'] = float(np.mean(blk_of(np.floor(theta(pr[TE]) / pi)) == blk_true))
print("Gram's law (index form):", OUT['gram_law_index_form'], " Rosser block hit:", OUT['rosser_block_hit'])
print("string block hits:", {k: OUT[k] for k in OUT if k.endswith('_block_hit')})

# where does each predictor fail: in Gram intervals that break the law?
counts = np.array(g_json['counts'])
mlist = ms[:-1]
cnt_of_int = dict(zip(mlist.tolist(), counts.tolist()))
in_bad = np.array([cnt_of_int.get(int(m), 1) != 1 for m in gi[TE]])
for lab, pr in [('G', G), ('P', P), ('PL', PL)]:
    e = np.abs(pr[TE] - gam[TE]) / mean_spacing(gam[TE])
    OUT[f'{lab}_err_lawful'] = float(np.nanmean(e[~in_bad]))
    OUT[f'{lab}_err_lawless'] = float(np.nanmean(e[in_bad]))
print({k: round(v, 4) for k, v in OUT.items() if k.endswith('_lawful') or k.endswith('_lawless')})

json.dump(dict(OUT, G=G.tolist(), P=P.tolist(), PL=PL.tolist(), L=Lp.tolist()), open('predict.json', 'w'),
          default=float)
