"""DRAFT. Second pass on the orbits the first pass missed: more candidate starts (12), a finer zoom (6 levels) and a
polish from each. Same rules: no catalog start used, period law refitted without the target."""
import json, os, sys, time, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import predict2 as P
from tbw import polish, raw
R = json.load(open(os.path.join(H, 'predict2_sample.json'))); outp = os.path.join(H, 'retry.json')
Q = json.load(open(outp)) if os.path.exists(outp) else {}
HIT = lambda f, o: f['res'] < 1e-7 and f['dist'] < 1e-6
for cls, rec in R.items():
    if cls in Q: continue
    if any(HIT(f, None) for r in rec['rotations'] for f in r['finals']): continue
    i = [j for j, o in enumerate(P.d) if o['cls'] == cls][0]; o = P.d[i]; t0 = time.time()
    s = P.inverse(o['word']); Ts = float(np.array(P.rfeat(o['word'])) @ P.period_law_without(i))
    cands, k, n = P.predict(s, Ts, nbest=8, levels=5, shrink=3); fins = []
    for (best, c1, c2, tg, _, _) in cands:
        q, r = polish(c1, c2, tg); dist = float(np.hypot(q[0] - o['v1'], q[1] - o['v2']))
        wq = raw(q[0], q[1], q[2] * 0.999999, maxev=4 * len(s) + 8)[0] if r < 1e-6 else ''
        fins.append(dict(guess_dist=float(np.hypot(c1 - o['v1'], c2 - o['v2'])), res=r, dist=dist, word_ok=(wq == s), v1=float(q[0]), v2=float(q[1]), T=float(q[2])))
        if HIT(fins[-1], o): break
    Q[cls] = dict(Lf=o['Lf'], finals=fins, hit=any(HIT(f, o) for f in fins), twin=any(f['res'] < 1e-7 and f['word_ok'] and f['dist'] >= 1e-6 for f in fins), secs=time.time() - t0)
    json.dump(Q, open(outp, 'w'))
    print(cls, o['Lf'], 'HIT' if Q[cls]['hit'] else ('twin' if Q[cls]['twin'] else 'miss'), 'best guess %.1e' % min(f['guess_dist'] for f in fins), 'secs %.0f' % Q[cls]['secs'], flush=True)
