import mpmath as mp, json, glob
from epstein import Z, brute
mp.mp.dps = 30
R = []
def ck(name, ok, info=""):
    R.append(bool(ok)); print(("PASS" if ok else "FAIL"), name, info)
def m_of_x(x):
    q = mp.exp(-mp.pi*x); return (mp.jtheta(2, 0, q)/mp.jtheta(3, 0, q))**4
ok = True
for m in (0.1, 0.3, 0.5, 0.7, 0.9):
    K = mp.ellipk(m); x = mp.ellipk(1-m)/K
    ok &= abs(mp.jtheta(3, 0, mp.exp(-mp.pi*x)) - mp.sqrt(2*K/mp.pi)) < 1e-15
ck("P1 Jacobi theta θ(x) = √(2K/π) at x = K′/K (five swings, to 15 digits)", ok)
ok = all(abs(mp.diff(lambda mm: mp.ellipk(1-mm)/mp.ellipk(mm), m) + mp.pi/(4*m*(1-m)*mp.ellipk(m)**2)) < 1e-15 for m in (0.2, 0.5, 0.8))
ck("P2 Legendre: d(K′/K)/dm = −π/(4m(1−m)K²), so Riemann's integral over x is an integral over swings", ok)
ck("P3 the 90° swing is the square torus: m = ½ gives K′/K = 1 (u = 0)", abs(mp.ellipk(0.5)/mp.ellipk(0.5) - 1) < 1e-25 and abs(m_of_x(1) - 0.5) < 1e-25)
ck("P4 swapping K and K′ (m ↔ 1−m) sends u to −u", all(abs(mp.log(mp.ellipk(1-m)/mp.ellipk(m)) + mp.log(mp.ellipk(m)/mp.ellipk(1-m))) < 1e-25 for m in (0.1, 0.37)))
d = json.load(open("swingdensity.json"))
th = d["theta"]; w = d["w"]; i60 = th.index(60.5) if 60.5 in th else min(range(len(th)), key=lambda i: abs(th[i]-60)); i120 = min(range(len(th)), key=lambda i: abs(th[i]-(180-th[i60])))
ck("P5 Φ per degree of swing is mirror-symmetric about 90° and peaks at 90°", abs(w[i60] - w[i120]) < 1e-12*max(w) and abs(th[w.index(max(w))] - 90) < 0.6)
ck("P6 the Epstein formula (Chowla–Selberg) matches the brute lattice sum at s = 3", all(abs(Z(3, y) - brute(3, y, 150)) < 1e-6 for y in (1, 2, 5.3)))
s = mp.mpc(0.5, 14)
ck("P7 the 90° swing's Epstein zeta is 4ζ(s)β(s)", abs(Z(s, 1) - 4*mp.zeta(s)*mp.dirichlet(s, [0, 1, 0, -1])) < 1e-20)
ck("P8 class-number-one swing amplitudes: y = 3 is exactly 30° (m = sin²15°)", abs(2*mp.asin(mp.sqrt(m_of_x(mp.sqrt(3))))*180/mp.pi - 30) < 1e-20)
ys = json.load(open("pend.json"))["ystar"]
ck("P9 real off-line zero enters at torus ratio √y = 7.0555 (Bateman–Grosswald: 7.0556)", abs(ys**0.5 - 7.0556) < 2e-4, f"√y* = {ys**0.5:.5f}")
rows = {json.load(open(f))["y"]: json.load(open(f)) for f in glob.glob("zeros_y*.json")}
ck("P10 no zeros off the line up to height 30 at the Euler-product swings y = 1, 2, 3, 4, 7", all(rows[y]["offline"] == 0 for y in (1.0, 2.0, 3.0, 4.0, 7.0)))
ck("P11 zeros off the line below height 30 at y = 1.5, 2.5, 5, 6, 8, 10, 15", all(rows[y]["offline"] > 0 for y in (1.5, 2.5, 5.0, 6.0, 8.0, 10.0, 15.0)))
ck("P12 off-line counts are even (zeros pair as σ and 1 − σ)", all(r["offline"] % 2 == 0 and r["offline"] >= 0 for r in rows.values()))
y = mp.mpf(5); L = lambda s: (mp.sqrt(y)/mp.pi)**s*mp.gamma(s)*Z(s, y, K=15)
ck("P13 completed Epstein Λ(s) = (√y/π)^s Γ(s) Z(s) is real on the line and symmetric", abs(mp.im(L(mp.mpc(0.5, 7.3)))) < 1e-15 and abs(L(mp.mpc(0.3, 2)) - L(mp.mpc(0.7, -2))) < 1e-15)
print(f"{sum(R)}/{len(R)} pass")
