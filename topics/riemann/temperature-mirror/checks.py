"""Checks for temperature-mirror. Run: python3 checks.py. DRAFT."""
import numpy as np, json, sys
import mirror
sys.path.insert(0, "../pendulum-swap"); from linecount import Fv
from oddrewind_dual import setup_frame
from oddpath import setup_odd
from frames import to_bonds
from knauf import chain
ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
def load(fr, n): return dict(json.loads(l) for l in open("mz_%s_%d_compass.jsonl" % (fr, n)))
z = np.array([0.3+5j, -0.2+11.3j, 0.05+27.1j])
for fr, n, Sref in (("spin", 10, setup_odd(10)), ("dual", 10, setup_frame(to_bonds(np.log(chain(10).astype(float)), 10), 10))):
    S = mirror.frame(fr, n); good = True
    for d in (0, 12, 50, 90):
        w, X = mirror.terms(S, np.radians(d), "compass"); G = lambda q: mirror.ev(w, X, q)
        good &= np.allclose(G(z) + G(-z), Fv(Sref, z, np.radians(d)), rtol=1e-10)
    check("mirror built in: F = G(z) + G(-z) exactly, %s frame n=10, 4 angles" % fr, good)
for n in (8, 10, 12):
    check("spin frame n=%d: G one-sided (no wrong-side zeros) at every tested angle" % n, all(v == 0 for v in load("spin", n).values()))
check("spin frame n=14: G one-sided at 13 angles 0..90", all(v == 0 for v in load("spin", 14).values()) and len(load("spin", 14)) == 13)
lc = json.load(open("../pendulum-swap/linecount.json"))
for n in (10, 12):
    first_fail = min(k for k, v in load("dual", n).items() if v > 0); first_leave = min(r[0] for r in lc["dual frame n=%d" % n])
    check("dual frame n=%d: one-sidedness fails (%d deg) before zeros leave (%d deg)" % (n, first_fail, first_leave), first_fail < first_leave)
    check("dual frame n=%d: wherever F leaves the line, G has a wrong-side zero (forced by the theorem)" % n,
          all(load("dual", n).get(r[0], 1) > 0 for r in lc["dual frame n=%d" % n]))
check("dual frame n=14: one-sidedness fails from 2 deg", load("dual", 14)[2] > 0 and load("dual", 14)[0] == 0)
# Hermite-Biehler inequality and interlacing, spin frame n=10 at 45 deg
S = mirror.frame("spin", 10); w, X = mirror.terms(S, np.radians(45), "compass")
xs = np.linspace(0.02, 3, 60); ys = np.linspace(0.3, 36, 700); Z = (xs[None, :] + 1j*ys[:, None]).ravel()
r = np.abs(mirror.ev(w, X, Z))/np.abs(mirror.ev(w, X, -Z))
check("spin frame n=10, 45 deg: |G(z)| > |G(-z)| on the whole right-half grid", r.min() > 1)
y = np.arange(0.3, 36, 0.0005); g = mirror.ev(w, X, 1j*y)
fz = y[1:][np.sign(g.real[1:]) != np.sign(g.real[:-1])]; oz = y[1:][np.sign(g.imag[1:]) != np.sign(g.imag[:-1])]
allz = sorted([(v, 0) for v in fz] + [(v, 1) for v in oz]); kinds = [k for _, k in allz]
check("zeros of F and of G(z) - G(-z) interlace on the line (spin, n=10, 45 deg)", all(kinds[i] != kinds[i+1] for i in range(len(kinds)-1)))
check("their count matches the box winding (all F zeros accounted for on the line)", len(fz) == mirror.winding(np.concatenate([w, w]), np.concatenate([X, -X]), -1, 1, 0.3, 36.0))
fe = json.load(open("farend.json")); pos = json.load(open("farend_pos.json"))
check("far end: the full chain has zeros right of Re beta = 3/2 at every length 6..16", all(v["Z_zeros_right_of_line"] > 0 for v in fe.values()))
check("far end: mirrored full chain keeps under half its zeros on the line, n=6..16", all(v["mirrored_on_line"] < 0.5*v["mirrored_total"] for v in fe.values()))
zt = (14.1347, 21.0220, 25.0109, 30.4249, 32.9351)
copies = [p for n in pos for p in pos[n]["right"] if min(abs(p[1]-t) for t in zt) < 0.15]
check("far end: right-side zeros near zeta's heights sit within 0.04 of the line", len(copies) >= 5 and all(p[0] < 0.04 for p in copies))
print("%d/%d passed" % (sum(ok), len(ok)))
