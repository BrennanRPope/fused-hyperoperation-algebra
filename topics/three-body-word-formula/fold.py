"""DRAFT. Why can one word carry two orbits? Follow both orbits of each same-word pair as the third mass m3 moves
away from 1 (masses 1, 1, m3). If the two branches run into each other and vanish together, they were born as a pair
(a fold, or saddle-node). Writes fold.json."""
import json, os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import tbm
def follow(p, direction, step=float(os.environ.get("STEP", 0.004)), mmax=float(os.environ.get("MMAX", 0.25)), ok=float(os.environ.get("OK", 1e-9))):
    path = [(1.0, *p)]; m = 1.0; h = step
    while abs(m - 1) < mmax and h > step / 64:
        mn = m + direction * h
        if len(path) >= 2:
            a, b = np.array(path[-2][1:]), np.array(path[-1][1:]); g = b + (b - a) * h / abs(path[-1][0] - path[-2][0])
        else: g = np.array(path[-1][1:])
        tbm.masses(mn); q, r = tbm.polish(*g)
        if r < ok and np.hypot(*(q[:2] - path[-1][1:3])) < 0.05 * max(h / step, 1):
            path.append((mn, *map(float, q))); m = mn; h = min(h * 1.5, step)
        else: h /= 2
    tbm.masses(1.0); return path
def run(name, a, b):
    out = {}
    for sgn, lab in ((1, 'up'), (-1, 'down')):
        A = follow(a, sgn); B = follow(b, sgn)
        out[lab] = dict(A=A, B=B)
        mA, mB = A[-1][0], B[-1][0]
        # distance between branches at common masses
        common = sorted(set(round(x[0], 6) for x in A) & set(round(x[0], 6) for x in B))
        dA = {round(x[0], 6): x for x in A}; dB = {round(x[0], 6): x for x in B}
        out[lab]['dist'] = [(m, float(np.hypot(dA[m][1] - dB[m][1], dA[m][2] - dB[m][2]))) for m in common]
        print(name, lab, 'A ends at m3=%.4f, B ends at m3=%.4f' % (mA, mB), 'last common dist %.2e' % out[lab]['dist'][-1][1] if out[lab]['dist'] else '', flush=True)
    return out
if __name__ == '__main__':
    T = json.load(open(os.path.join(H, 'twins.json'))); which = sys.argv[1:] or list(T)
    fp = os.path.join(H, os.environ.get('FOLDOUT', 'fold.json')); R = json.load(open(fp)) if os.path.exists(fp) else {}
    for c in which:
        R[c] = run(c, tuple(T[c]['cat']), tuple(T[c]['twin'])); json.dump(R, open(fp, 'w'))
