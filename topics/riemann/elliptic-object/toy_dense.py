"""Detection toy at every n (no sampling gaps): split a double zero at height g0 off the line by delta,
and find the first n with lambda_n < 0. Background = zero sum over 5,000 zeros + density tail
(tail interpolated from dial.json, it is smooth in n). Rewrites dial.json['toy']."""
import json, math, numpy as np
d = json.load(open('dial.json'))
gam = np.array([r['gamma'] for r in d['C']]); th = 2 * np.arctan(1 / (2 * gam))
nD = np.array([r['n'] for r in d['D']]); tD = np.array([r['tail'] for r in d['D']])
N = np.arange(1, 20001)
tail = np.interp(N, nD, tD)
zs = np.zeros(len(N))
for i in range(0, len(N), 500):
    nn = N[i:i + 500, None]; zs[i:i + 500] = (2 * (1 - np.cos(nn * th[None, :]))).sum(1)
base = zs + tail
def quartet(g0, dl):
    tot = np.zeros(len(N))
    for b in (0.5 + dl, 0.5 - dl):
        for sg in (1, -1):
            z = 1 - 1 / complex(b, sg * g0)
            tot += (1 - np.exp(N * np.log(z))).real
    return tot
toy = []
for g0, dl in ((5.0, 0.25), (5.0, 0.45), (14.134725, 0.45), (14.134725, 0.25)):
    th0 = 2 * math.atan(1 / (2 * g0))
    moved = base - 4 * (1 - np.cos(N * th0)) + quartet(g0, dl)
    neg = np.nonzero(moved < 0)[0]
    fn = int(N[neg[0]]) if len(neg) else None
    stop = min(len(N), int(1.15 * fn) if fn else len(N))
    eps = 0.5 * math.log1p(2 * dl / ((0.5 - dl) ** 2 + g0 ** 2))
    toy.append({"gamma0": g0, "delta": dl, "eps": eps, "first_negative_n": fn,
                "curve": [[int(N[i]), float(moved[i])] for i in range(0, stop)]})
    print(g0, dl, 'first negative n', fn, ' 1/eps', 1 / eps)
d['toy'] = toy
d['base_dense_check'] = [[int(n), float(base[n - 1])] for n in (1, 10, 100, 1000)]
json.dump(d, open('dial.json', 'w'))
