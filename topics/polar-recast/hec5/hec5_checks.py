"""Sanity checks for hec5_signs.py: the 372 facets hold on a perfect-tensor state and GHZ breaks some."""
import itertools, importlib.util, sys
spec = importlib.util.spec_from_file_location("h", "hec5_signs.py"); h = importlib.util.module_from_spec(spec)
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(h)
from collections import Counter
def all_instances():
    out = []
    for k, (L, R) in h.REPS.items():
        base = h.parse(L, 1); base.update(h.parse(R, -1)); seen = set()
        for p in itertools.permutations(h.SYM):
            mp = dict(zip(h.SYM, p)); inst = h.canon(Counter({frozenset(mp[x] for x in I): c for I, c in base.items()}))
            if inst not in seen: seen.add(inst); out.append((k, inst))
    return out
inst = all_instances()
def value(inst, S): return sum(c*S(I) for I, c in inst)
pt = lambda I: min(len(I), 6 - len(I))          # 6-qubit perfect tensor (AME(6,2)), regions A..E + purifier
ghz = lambda I: 1                                # 6-qubit GHZ: every nonempty proper subset has entropy 1
v_pt = [value(i, pt) for _, i in inst]; v_ghz = [value(i, ghz) for _, i in inst]
print(f"facets: {len(inst)}")
print(f"PASS perfect tensor satisfies all 372: min value {min(v_pt)}" if min(v_pt) >= 0 else "FAIL perfect tensor")
bad = Counter(h.NAME[k] for (k, _), v in zip(inst, v_ghz) if v < 0)
print(f"PASS GHZ breaks some holographic rules: {dict(bad)}" if bad else "FAIL GHZ breaks none")
# corner weights of the two states (y = Mobius of S on regions A..E)
def ycorners(S):
    y = {}
    for r in range(1, 6):
        for J in itertools.combinations("ABCDE", r):
            y["".join(J)] = sum((-1)**(r - len(K))*S(frozenset(K)) for q in range(r + 1) for K in itertools.combinations(J, q) if K)
    return y
for nm, S in [("perfect tensor", pt), ("GHZ", ghz)]:
    y = ycorners(S); by = {g: sorted({v for J, v in y.items() if len(J) == g}) for g in range(1, 6)}
    print(f"corner weights by grade, {nm}: {by}")
