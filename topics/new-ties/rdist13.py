"""Round 13 (DRAFT): two-sided estimate of dist(rho, R) by branch and bound on cylinder pairs.

R = { x/y : x coded by an infinite reduced word starting with N, y by one starting with N or K }
(first letter outermost; after m only m or N may follow). A pair of cylinders u(I), v(I) holds ratios in
[min u / max v, max u / min v]; its distance from rho is a LOWER bound for the pair, and u(0)/v(0) is a
genuine point of R (0 is the fixed point of m, which may follow any letter), giving an UPPER bound.
Best-first refinement until LB >= UB / 2 (or a node cap). Returns (LB, UB, x-word, y-word, nodes).
"""
import heapq, mpmath as mp
mp.mp.dps = 120
LAM, MU = mp.mpf(19) / 200, mp.mpf(1) / 20
LET = {"l": (LAM, mp.mpf(0)), "m": (MU, mp.mpf(0)), "N": (mp.mpf(1) / 19, mp.mpf(1) / 19), "K": (mp.mpf(1) / 10, mp.mpf(1) / 10)}
AFTER = {"l": "lmNK", "m": "mN", "N": "lmNK", "K": "lmNK"}
NINTH = mp.mpf(1) / 9

def ends(w):   # affine u(x) = a x + c ; cylinder u(I) = [c, a/9 + c]
    a, c = w[1], w[2]; return c, a * NINTH + c

def ext(w, x):
    s, t = LET[x]; return (w[0] + x, w[1] * s, w[1] * t + w[2])

def dist_R(rho, ratio=2, cap=200000):
    rho = mp.mpf(rho)
    start = [(("N", *LET["N"][::-1][::-1]))]
    X0 = [("N", LET["N"][0], LET["N"][1])]; Y0 = [(L, LET[L][0], LET[L][1]) for L in "NK"]
    best = (mp.inf, None, None); heap = []; cnt = 0
    def push(u, v):
        nonlocal best, cnt
        x0, x1 = ends(u); y0, y1 = ends(v)
        lo, hi = x0 / y1, x1 / y0
        lb = lo - rho if rho < lo else (rho - hi if rho > hi else mp.mpf(0))
        ub = abs(rho - u[2] / v[2])
        if ub < best[0]: best = (ub, u[0], v[0])
        cnt += 1; heapq.heappush(heap, (lb, cnt, u, v))
    for u in X0:
        for v in Y0: push(u, v)
    nodes = 0
    while heap and nodes < cap:
        lb, _, u, v = heapq.heappop(heap)
        if lb * ratio >= best[0]: return (lb, best[0], best[1], best[2], nodes, True)
        nodes += 1
        x0, x1 = ends(u); y0, y1 = ends(v)
        if (x1 - x0) / x0 >= (y1 - y0) / y0:
            for x in AFTER[u[0][-1]]: push(ext(u, x), v)
        else:
            for y in AFTER[v[0][-1]]: push(u, ext(v, y))
    lb = heap[0][0] if heap else best[0]
    return (lb, best[0], best[1], best[2], nodes, False)

if __name__ == "__main__":
    import sys
    LL, L20 = mp.log(LAM), mp.log(20)
    for s in sys.argv[1:]:
        a, b = map(int, s.split(","))
        r = mp.exp(a * LL + b * L20)
        lb, ub, xw, yw, n, done = dist_R(r)
        print(a, b, mp.nstr(r, 12), "log2 LB %.2f UB %.2f" % (float(mp.log(lb, 2)) if lb > 0 else -999, float(mp.log(ub, 2))), n, done, xw[:40], yw[:40])
