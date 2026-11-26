"""Round 9: how exactly must the orbit lengths be log p?  Score against the ratio dial r
on a fine scale, split by height: the allowed error shrinks like 1/(E log X)."""
import json
import numpy as np
from stringlib import String, score
gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
bands = {'1-500': np.arange(0, 500), '501-1500': np.arange(500, 1500), '1501-5500': np.arange(1500, NZ)}
OUT = {}
for dr in [0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]:
    st = String(taper='gauss', X=200.0, r=1.0 + dr)
    E = st.levels(NZ)
    OUT[str(dr)] = {b: score(E, gam, idx)['rms'] for b, idx in bands.items()}
    print(dr, {b: round(v, 4) for b, v in OUT[str(dr)].items()})
json.dump(OUT, open('rwidth.json', 'w'))
