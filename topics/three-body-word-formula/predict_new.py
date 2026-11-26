"""DRAFT. Generative test: words that are NOT in Li & Liao's list. Every catalog word read from its start is
s = x + mirror(x) with mirror(x) = reverse of x with a<->b swapped. Make new half-words x by splicing halves of two
catalog words, keep runs of length 2-4 and length <= 40, drop anything equivalent to a catalog word, then ask the
word -> start pipeline (letters + run-law period + closest return + polish) for an orbit with exactly that word."""
import json, os, sys, time, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import predict2 as P
from tbw import polish, raw
inv = lambda w: w[::-1].swapcase(); mir = lambda x: x[::-1].translate(str.maketrans('abAB', 'baBA'))
def reduced(w): return all(w[i] != w[(i + 1) % len(w)].swapcase() for i in range(len(w)))
def canon(w):
    ab = lambda x: x.translate(str.maketrans('abAB', 'baBA')); c = set()
    for v in (w, inv(w), w.swapcase(), inv(w.swapcase())):
        for u in (v, ab(v)):
            for k in range(len(u)): c.add(u[k:] + u[:k])
    return min(c)
CAT = {canon(o['word']) for o in P.d}
halves = [inv(o['word'])[:o['Lf'] // 2] for o in P.d]
coef = np.linalg.lstsq(P.RF, P.TS, rcond=None)[0]
rng = np.random.default_rng(2026); words = []; seen = set()
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 30
while len(words) < NW:
    x1, x2 = halves[rng.integers(len(halves))], halves[rng.integers(len(halves))]
    i, j = rng.integers(1, len(x1) + 1), rng.integers(0, len(x2))
    x = x1[:i] + x2[j:]; s = x + mir(x)
    if not (8 <= len(s) <= 40) or not reduced(s) or not set(P.runs(s)) <= {2, 3, 4}: continue
    c = canon(s)
    if c in CAT or c in seen: continue
    seen.add(c); words.append(s)
outp = os.path.join(H, 'predict_new.json'); R = json.load(open(outp)) if os.path.exists(outp) else {}
for s in words:
    if s in R: continue
    t0 = time.time(); Ts = float(np.array(P.rfeat(inv(s))) @ coef)
    cands, k, nreg = P.predict(s, Ts); fins = []
    for (best, c1, c2, tg, _, _) in cands[:4]:
        q, r = polish(c1, c2, tg)
        wq = raw(q[0], q[1], q[2] * 0.999999, maxev=4 * len(s) + 8)[0] if r < 1e-6 else ''
        fins.append(dict(guess=[float(c1), float(c2), float(tg)], guess_return=float(best), v1=float(q[0]), v2=float(q[1]), T=float(q[2]), res=r, word=wq, exact=(wq == s and r < 1e-8)))
        if fins[-1]['exact']: break
    found = [f for f in fins if f['exact']]
    R[s] = dict(Lf=len(s), Tstar_pred=Ts, k=k, region_cells=nreg, finals=fins, found=bool(found),
                Tstar_found=(found[0]['T'] * abs(-2.5 + 3 * (found[0]['v1']**2 + found[0]['v2']**2))**1.5 if found else None), secs=time.time() - t0)
    json.dump(R, open(outp, 'w'))
    print(len(s), s, 'FOUND' if found else 'none', 'k', k, 'cells', nreg, 'T* pred %.3f found %s' % (Ts, ('%.3f' % R[s]['Tstar_found']) if found else '-'), flush=True)
