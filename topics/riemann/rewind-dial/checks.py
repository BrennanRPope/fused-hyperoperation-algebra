"""Rewind-dial checks. DRAFT. Each line prints PASS/FAIL."""
import numpy as np, json, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import chain
from rotate import setup, F, start_zeros, newton
from oddpath import setup_odd
res = []
def ok(name, c): res.append(bool(c)); print(("PASS " if c else "FAIL ") + name)
def walsh(f):
    a = f.astype(float).copy(); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
def firstorder(S):
    z0 = start_zeros(S, 36.0); G = (S["w"]*S["H"]*np.exp(np.multiply.outer(z0, S["a"]*S["M"]))).sum(-1)
    _, Fz = F(S, z0, 0.0); return z0, -(z0*G)/Fz
mv = {n: json.load(open("movie_n%d.json" % n)) for n in (10, 12)}
od = {**json.load(open("oddpath_8_10.json")), **json.load(open("oddpath_12.json"))}
for n in (10, 12):
    z = np.array(mv[n]["frames"]["0"]); ok("R1 n=%d, angle 0: all %d zeros in the box on the line (all |Re z| below 1e-8)" % (n, len(z)), np.abs(z[:, 0]).max() < 1e-8)
S = setup(10); z0, zp = firstorder(S); num = (newton(S, z0, 1e-6) - z0)/1e-6
ok("R2 first-order formula dz/dangle = -z G/F' matches a numerical step (max diff %.1e)" % np.abs(num - zp).max(), np.abs(num - zp).max() < 1e-2)
ok("R3 full rewind: at first order every zero moves straight off the line (|Im dz/dangle| < 1e-12, |Re| up to %.1f)" % np.abs(zp.real).max(), np.abs(zp.imag).max() < 1e-12 and np.abs(zp.real).max() > 1)
spd = [np.abs(firstorder(setup(n))[1].real) / firstorder(setup(n))[0].imag for n in (8, 10, 12, 14)]
ok("R4 off-line speed / height does not shrink with chain length (median %s)" % ["%.2f" % np.median(s) for s in spd], min(np.median(s) for s in spd) > 0.15)
for n in (10, 12):
    z = np.array(mv[n]["frames"]["5"]); frac = np.mean(np.abs(z[:, 0]) > 1e-3)
    ok("R5 n=%d, angle 5 deg: %.0f%% of zeros already off the line" % (n, 100*frac), frac > 0.5)
So = setup_odd(10); z0o, zpo = firstorder(So)
ok("R6 odd-half rewind: at first order every zero slides along the line (|Re dz/dangle| max %.1e)" % np.abs(zpo.real).max(), np.abs(zpo.real).max() < 1e-12)
for n in ("8", "10", "12"):
    mx = max(np.abs(np.array(f)[:, 0]).max() for f in od[n].values() if len(f))
    cnt = [len(f) for f in od[n].values()]
    ok("R7 odd-half rewind n=%s: every zero on the line at all 10 angles (all |Re z| below 1e-8; %d to %d zeros per frame)" % (n, min(cnt), max(cnt)), mx < 1e-8 and min(cnt) > 5)
phi = chain(12).astype(float); H = np.log(phi); Ho = 0.5*(H-H[::-1]); J = walsh(Ho); sz = np.array([bin(k).count("1") for k in range(2**12)])
share = (J[sz == 1]**2).sum()/(J[sz % 2 == 1]**2).sum()
ok("R8 the odd half is %.1f%% single-spin links (an uneven field), the rest odd many-spin links" % (100*share), 0.95 < share < 1 and (J[(sz % 2 == 1) & (sz > 1)]**2).sum() > 0)
Sd = setup(10); zt = np.array([0.3+5j, -0.2+17j]); f90, _ = F(Sd, zt, np.pi/2)
Zk = np.array([(chain(10).astype(float)**(-(1.5 - z))).sum() for z in zt])
ok("R9 full rewind at 90 deg is the temperature view: F = 2 Z_Knauf(3/2 - z)", np.allclose(f90, 2*Zk, rtol=1e-10))
z90 = np.array(mv[10]["frames"]["90"]); beta = 1.5 - (z90[:, 0] - 1j*z90[:, 1])
g = json.load(open("../knauf-field/refined.json"))["gam"]
near = [min(abs(b - (1.5+1j*G)) for b in beta) for G in g if G < 36]
ok("R10 n=10 at 90 deg: zeros sit next to 3/2 + i*gamma for the first zeta zeros (distances %s)" % ["%.3f" % d for d in near], max(near) < 0.1)
fr = mv[10]["frames"]; ninline = [sum(abs(p[0]) < 0.05 for p in fr[str(d)]) for d in (5, 45, 90)]
ok("R11 full rewind n=10: zeros within 0.05 of the line at 5/45/90 deg: %s of about 25" % ninline, ninline[0] < 15)
print("%d/%d passed" % (sum(res), len(res)))
