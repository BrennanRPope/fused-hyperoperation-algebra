"""Resolution run with checkpoints (resumes after a restart): pair alone, nu=0.001, cube, grid N.
usage: WORKERS=4 GAMMA=8 DELTA=0.5 python3 res256.py N T out.json ckpt.npz"""
import sys, os, json, time
import numpy as np
import tubes
from flow3d import Box, NS3
N, T, out, ck = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
nu, every = 0.001, 0.1
box = Box(N, 1.0, Nz=N); ns = NS3(box, nu, 'T3')
if os.path.exists(ck):
    d = np.load(ck, allow_pickle=True); U = d['U']; t = float(d['t']); nxt = float(d['nxt'])
    rec = json.loads(str(d['rec'])); meta = json.loads(str(d['meta'])); print('resumed at t=%.3f' % t, flush=True)
else:
    U = tubes.initial(box, 'none', 0, 0, 0.25, 2); t = 0.0; nxt = 0.0
    rec = dict(t=[], wmax=[], Z=[], E=[], tail=[])
    meta = dict(E0=ns.energy(U), Z0=ns.enstrophy(U), W0=ns.wmax(U), secs=0.0)
t0 = time.time(); last_ck = time.time()
def save():
    m = dict(meta); m['secs'] = meta['secs'] + time.time() - t0
    np.savez(ck + '.tmp.npz', U=U, t=t, nxt=nxt, rec=json.dumps(rec), meta=json.dumps(m))
    os.replace(ck + '.tmp.npz', ck)
while t < T - 1e-12:
    if t >= nxt - 1e-12:
        rec['t'].append(t); rec['wmax'].append(ns.wmax(U)); rec['Z'].append(ns.enstrophy(U))
        rec['E'].append(ns.energy(U)); rec['tail'].append(tubes.tail(ns, U)); nxt += every
        print('t=%.2f wmax=%.1f Z/Z0=%.3f tail=%.1e' % (t, rec['wmax'][-1], rec['Z'][-1]/meta['Z0'], rec['tail'][-1]), flush=True)
        save(); last_ck = time.time()
    dt = min(ns.dt_cfl(U, 0.8), T - t, 0.05)
    U = ns.step(U, dt); t += dt
    if time.time() - last_ck > 600:
        save(); last_ck = time.time()
rec['t'].append(t); rec['wmax'].append(ns.wmax(U)); rec['Z'].append(ns.enstrophy(U))
rec['E'].append(ns.energy(U)); rec['tail'].append(tubes.tail(ns, U))
i = int(np.argmax(rec['wmax']))
r = dict(N=N, nu=nu, T=T, cls='none', seed=0, eps=0, a=0.25, kbend=2, E0=meta['E0'], Z0=meta['Z0'], W0=meta['W0'],
         wpeak=rec['wmax'][i], tpeak=rec['t'][i], Zpeak=max(rec['Z']), tailmax=max(rec['tail']),
         secs=meta['secs'] + time.time() - t0, **rec)
json.dump(r, open(out, 'w')); save()
print('done wpeak %.1f at t=%.2f tail %.1e' % (r['wpeak'], r['tpeak'], r['tailmax']), flush=True)
