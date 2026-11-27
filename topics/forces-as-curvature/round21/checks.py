"""DRAFT round 21 checks: why Lorentzian? What in the law, if anything, makes time's sign differ from space's.

Round 20 built de Sitter spacetime from the law with one dial q, but only by giving its dot product one minus sign. Here the
law is asked what that sign does. Every result below uses only q.eta(v, v), the dial times the dot product along a direction,
because that is all the law sees: flipping both signs together leaves the law unchanged.
Run: python3 checks.py  (numpy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools, time
import numpy as np

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

def oplus(a, b, qq, eta):         # the law of rounds 12-20, word for word, dot product eta
    d = lambda u, v: (u * v * eta).sum(-1, keepdims=True)
    ab, aa, bb = d(a, b), d(a, a), d(b, b)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
rg = np.random.default_rng(21)

# ================================================================== W. What the sign does
say("== W. Why Lorentzian?")
# W1: along one direction the law is a Moebius map with matrix [[1, a], [-c a, 1]], c = q.eta(v, v).  Its kind decides everything.
def orbit(c, a, n):
    x, out = 0.0, [0.0]
    for _ in range(n):
        den = 1 - c * a * x
        if abs(den) < 1e-12:
            out.append(np.inf); x = np.inf; continue
        x = (x + a) / den if np.isfinite(x) else -1 / (c * a); out.append(x)
    return np.array(out)
N = 12; ang = math.pi / N
per = orbit(1.0, math.tan(ang / 2), 2 * N)                     # c > 0: 2N kicks of angle 2*pi/(2N)... returns to start
hyp = orbit(-1.0, math.tanh(0.3 / 2), 40)                      # c < 0: rapidity steps 0.3
par = orbit(0.0, 0.1, 40)
ret_per = abs(per[-1]) < 1e-9 and not np.isfinite(per[N])   # reaches the far point half way round, then comes back exactly
mono = np.all(np.diff(hyp) > 0) and hyp.max() < 1.0
rap = 2 * np.arctanh(hyp)
e_rap = np.abs(np.diff(rap) - 0.3).max()
e_par = np.abs(np.diff(par) - 0.1).max()
kinds = {}
for c in (1.0, 0.0, -1.0):
    a = 0.4; M = np.array([[1, a], [-c * a, 1]]); kinds[c] = float(np.trace(M)**2 / np.linalg.det(M))
check("W1 (exact) along any one direction the law is a Moebius map with matrix [[1, a], [-c a, 1]], c = q.eta(v, v). If c > 0 it is a rotation: repeated kicks go round and come back exactly (space closing, rounds 12-13). If c < 0 it is a boost: repeated kicks march one way for ever, in equal rapidity steps, never reaching the edge 1/sqrt(-c) (an ordered direction with a horizon). If c = 0 it is plain addition. Which kind is fixed by the sign of c alone",
      ret_per and mono and e_rap < 1e-10 and e_par < 1e-12 and kinds[1.0] < 4 and kinds[-1.0] > 4 and abs(kinds[0.0] - 4) < 1e-12,
      f"c = +1: {2 * N} kicks pass the far point half way and return to the start ({abs(per[-1]):.0e}); c = -1: 40 kicks strictly increasing, largest {hyp.max():.6f} < 1, rapidity steps equal to {e_rap:.0e}; "
      f"c = 0: equal steps to {e_par:.0e}; trace^2/det = {kinds[1.0]:.3f} (< 4, rotation), {kinds[0.0]:.0f} (= 4, plain), {kinds[-1.0]:.3f} (> 4, boost)")
results["kinds"] = dict(periodic=per[np.isfinite(per)].tolist(), boost=hyp.tolist(), plain=par.tolist(), trace2_over_det=kinds)

# W2: light needs an indefinite dot.  A direction where the law ignores the dial exists iff eta(v, v) = 0 for some v != 0.
eu = np.ones(4); lo = np.array([-1.0, 1, 1, 1])
V = rg.normal(size=(5000, 4)); s1, s2 = rg.uniform(-1, 1, (5000, 1)), rg.uniform(-1, 1, (5000, 1))
one_d = lambda a_, b_, c: (a_ + b_) / (1 - c * a_ * b_)        # the 1D law with dial c
e_line = 0.0
for eta in (eu, lo):
    c = 0.5 * (V * V * eta).sum(1, keepdims=True)              # q.eta(v, v) with q = 0.5
    ok_ = np.abs(1 - c * s1 * s2)[:, 0] > 0.05                  # away from the far point, where both blow up
    e_line = max(e_line, (np.abs(oplus(s1 * V, s2 * V, 0.5, eta) - one_d(s1, s2, c) * V)[ok_] / np.abs(one_d(s1, s2, c) * V)[ok_].max(1, keepdims=True)).max())
c_eu = 0.5 * (V * V).sum(1)
nul = np.concatenate([np.ones((5000, 1)), V[:, 1:] / np.linalg.norm(V[:, 1:], axis=1)[:, None]], 1)
dev_null = np.abs(oplus(s1 * nul, s2 * nul, 0.5, lo) - (s1 + s2) * nul).max()
etaV = (V * V * lo).sum(1)
check("W2 (exact) light needs a minus sign. Along any direction v the law is the 1D law with dial c = q.eta(v, v), so it ignores the dial (plain addition, a signal whose speed no q changes) exactly when eta(v, v) = 0. With an all-plus dot product eta(v, v) = |v|^2 > 0 for every v != 0, so no direction escapes the dial. With one minus sign there is a whole cone of them",
      e_line < 1e-12 and c_eu.min() > 0 and dev_null < 1e-12,
      f"5000 random directions, both dot products: the law along v vs the 1D law with c = q.eta(v, v) {e_line:.0e} (relative, away from the far point); all-plus dot: smallest c {c_eu.min():.3f} > 0; "
      f"one minus sign: along the light cone, law vs plain addition {dev_null:.0e}")
results["light"] = dict(line_is_1d_law=float(e_line), min_c_euclid=float(c_eu.min()), null=float(dev_null))

# W3: time that loops.  Count the directions with c < 0 (ordered, 'time-like').  With two or more, a closed curve can stay time-like
# throughout (it turns round inside the time plane); with exactly one, every closed curve must stop moving in time somewhere.
def timelike_loop_exists(p, trials=4000):
    eta = np.array([-1.0] * p + [1.0] * (4 - p))
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    if p >= 2:
        loop = np.stack([np.cos(th), np.sin(th)] + [0 * th] * 2, 1)
        vel = np.stack([-np.sin(th), np.cos(th)] + [0 * th] * 2, 1)
        explicit = float(((vel * vel * eta).sum(1)).max())         # < 0 means time-like all the way round
    else:
        explicit = None
    worst = -np.inf; found = False
    for _ in range(trials):                                     # random closed loops (a few Fourier modes)
        C = rg.normal(size=(4, 4)) / (1 + np.arange(4))[None, :, None].squeeze(-1)
        S = rg.normal(size=(4, 4)) / (1 + np.arange(4))
        k = np.arange(1, 5)
        vel = (np.cos(np.outer(th, k)) @ (S * k).T - np.sin(np.outer(th, k)) @ (C * k).T)
        m = (vel * vel * eta).sum(1).max()                      # the loop's least time-like point
        if m < 0:
            found = True
        worst = max(worst, -m)
    return explicit, found
rows3 = {}
for p in range(5):
    ex, fd = timelike_loop_exists(p)
    rows3[p] = dict(explicit_max_eta=ex, random_found=fd)
proof_ok = all(not rows3[p]["random_found"] for p in (0, 1)) and all(rows3[p]["explicit_max_eta"] < 0 for p in (2, 3, 4))
check("W3 (exact for the explicit loops; random search for the rest) time can loop unless exactly one direction is ordered. With one minus sign, any closed curve must stop moving in time somewhere, so it can't stay time-like all the way round: time is ordered. With two or more, a circle in the time plane is time-like all the way round, a closed time-like loop",
      proof_ok,
      "; ".join(f"{p} minus sign{'s' if p != 1 else ''}: " + (f"circle in the time plane, largest eta(v, v) {rows3[p]['explicit_max_eta']:.2f} < 0" if p >= 2
                 else f"random closed loops time-like throughout: {'found' if rows3[p]['random_found'] else 'none'}") for p in range(5)))
results["loops"] = {str(k): v for k, v in rows3.items()}

# W4: the fused shift (round 15) is a forward run and its inverse run.  On an ordered direction they meet once (the bounce);
# on a closed direction they meet again and again, so a 'time' built that way has no single start.
def meetings(c, a, n=400):
    f = orbit(c, a, n); b = orbit(c, -a, n)
    return int(np.sum(np.isclose(f[1:], b[1:], atol=1e-9) | (~np.isfinite(f[1:]) & ~np.isfinite(b[1:]))))
m_time = meetings(-1.0, math.tanh(0.05)); m_space = meetings(1.0, math.tan(math.pi / 40))
check("W4 (exact) round 15 read time as the fused shift, a forward run of the law plus its inverse run. Along an ordered direction (c < 0) the two runs share only their start, the bounce. Along a closed direction (c > 0) they meet again every half turn, so time built that way would have no single bounce. So round 15's reading itself needs time to be the ordered kind",
      m_time == 0 and m_space >= 5,
      f"400 steps each way: ordered direction, meetings after the start {m_time}; closed direction (40 steps per half turn) {m_space}")
results["fused_shift"] = dict(ordered_meetings=m_time, closed_meetings=m_space)

# W5: the table.  For each number p of minus signs (q > 0 throughout; q < 0 swaps p and 4 - p), which needs are met?
table = []
for p in range(5):
    eta = np.array([-1.0] * p + [1.0] * (4 - p))
    c_signs = set(np.sign(eta))
    table.append(dict(minus_signs=p, space_closes=bool(1.0 in c_signs), time_ordered=bool(-1.0 in c_signs),
                      light=bool(len(c_signs) == 2), no_time_loops=bool(p <= 1),
                      dims_space=int((eta > 0).sum()), dims_time=int((eta < 0).sum())))
for row in table:
    row["all"] = row["space_closes"] and row["time_ordered"] and row["light"] and row["no_time_loops"]
winners = [r["minus_signs"] for r in table if r["all"]]
sym = np.abs(oplus(V[:9] * 0.3, V[9:18] * 0.3, 0.7, lo) - oplus(V[:9] * 0.3, V[9:18] * 0.3, -0.7, -lo)).max()
check("W5 (exact) only one minus sign meets all four needs: a space that closes (rounds 12-13), a direction that is ordered (time), light that ignores the dial (W2), and no time loops (W3). Flipping the dial's sign and the dot product's together leaves the law unchanged, so q < 0 with three minus signs is the same answer. So the law's own behaviour, plus these four needs, picks one time and three space directions",
      winners == [1] and sym < 1e-12,
      "minus signs 0-4: " + "; ".join(f"{r['minus_signs']}: closes {'y' if r['space_closes'] else 'n'}, ordered {'y' if r['time_ordered'] else 'n'}, light {'y' if r['light'] else 'n'}, no loops {'y' if r['no_time_loops'] else 'n'}" for r in table)
      + f"; law unchanged when both signs flip {sym:.0e}")
results["table"] = table

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
