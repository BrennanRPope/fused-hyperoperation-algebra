"""Data center energy through Fused Hyperoperation Algebra.
All inputs are published figures (see SOURCES.md). Mult-add: a(+)b = a+b+ab = (1+a)(1+b)-1.
A stage with efficiency eta has overhead odds r = 1/eta - 1 (input watts per output watt, minus 1)."""
import numpy as np, itertools, math
madd = lambda a, b: a + b + a*b
def fold(rs):
    t = 0.0
    for r in rs: t = madd(t, r)
    return t
odds = lambda eta: 1/eta - 1
out = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

p("=== Lead 1: the power chain composes by mult-add ===")
# mid-range losses: MV/LV distribution 2-3%, online UPS 2-4% (DCD: 94-97% eff), PSU 4-6%, VRM 8-12% (arXiv 2610.09060)
chain = [("Distribution and transformers", 0.975), ("UPS, double conversion", 0.96), ("Server PSU", 0.95), ("Board VRM", 0.90)]
rs = [odds(e) for _, e in chain]
eta = np.prod([e for _, e in chain])
p("stage odds", [round(r, 4) for r in rs])
p("chain efficiency", round(eta, 4), " fused odds", round(fold(rs), 4), " check 1/eta-1", round(1/eta-1, 4))
p("naive sum of fractional losses", round(sum(1-e for _, e in chain), 4), " true lost fraction", round(1-eta, 4))
inside = madd(rs[2], rs[3]); outside = madd(rs[0], rs[1])
p("inside the IT meter (PSU (+) VRM) odds", round(inside, 4), "-> 'IT' power is", round(1+inside, 4), "x chip power")
p("outside (dist (+) UPS) odds", round(outside, 4))

p("\n=== Lead 2: PUE overhead = electrical (+) cooling; the corner is cooling the power chain ===")
for name, PUE, PLF in [("Uptime 2024 average, PLF from chain above", 1.56, outside),
                       ("Hamilton 2008 (IT 59, dist 8, mech 33)", 100/59, 8/59)]:
    c = PUE/(1+PLF) - 1          # cooling+other watts per watt of heat (IT + electrical losses)
    corner = PLF*c
    CLF_add = PUE - 1 - PLF
    p(f"{name}: PUE {PUE:.3f}  PLF {PLF:.4f}  fused cooling odds c {c:.4f}  additive CLF {CLF_add:.4f}")
    p(f"   corner PLF*c = {corner:.4f} = {100*corner/(PUE-1):.1f}% of overhead; 1 W of electrical loss costs {1+c:.3f} W")
    if PUE < 1.6:
        p(f"   per chip watt: {PUE:.3f} * {1+inside:.4f} = {PUE*(1+inside):.3f} W at the meter; fused odds {madd(PUE-1, inside):.4f}")
US = 176.0  # TWh 2023, LBNL
c = 1.56/(1+outside)-1
p(f"US scale: overhead {US*0.56/1.56:.1f} TWh/yr; corner (cooling the power chain) {US*outside*c/1.56:.1f} TWh/yr (Uptime-average model)")

p("\n--- order-free blame: nepers (log of 1+odds) vs waterfall watts ---")
stages = chain + [("Cooling and other overhead", 1/(1+c))]
L = np.array([-math.log(e) for _, e in stages])
p("nepers", dict((n, round(l, 4)) for (n, _), l in zip(stages, L)), "total", round(L.sum(), 4), "=ln", round(math.exp(L.sum()), 4))
p("shares", dict((n, f"{100*l/L.sum():.1f}%") for (n, _), l in zip(stages, L)))
# position independence: improving any one stage's efficiency by the same ratio cuts the facility by the same fraction
p("upgrades (facility saving = 1 - eta_old/eta_new, independent of where the stage sits):")
for n, a, b in [("UPS double conversion -> eco/high-efficiency mode", 0.96, 0.985),
                ("PSU 80 PLUS Platinum 94% -> Titanium 96% (50% load)", 0.94, 0.96),
                ("VRM 88% -> 92%", 0.88, 0.92)]:
    s = 1 - a/b
    p(f"   {n}: {100*s:.2f}% of facility energy (US-scale {US*s:.1f} TWh/yr if applied everywhere)")

p("\n=== Lead 3: part load. Loss = a + b x + c x^2 (fixed, proportional, square). ===")
x = np.array([0.1, 0.2, 0.5, 1.0]); e = np.array([0.90, 0.94, 0.96, 0.91])  # 80 PLUS Titanium 230V internal redundant minimums
loss = x*(1/e - 1)
A = np.vstack([np.ones_like(x), x, x**2]).T
(a, b, cq), *_ = np.linalg.lstsq(A, loss, rcond=None)
p(f"Titanium-minimum fit: a={a:.5f} b={b:.5f} c={cq:.5f}; residuals {np.round(A@[a,b,cq]-loss,5)}")
xs = math.sqrt(a/cq)
p(f"peak efficiency at x* = sqrt(a/c) = {xs:.3f} of rating, eta* = {xs/(xs+a+b*xs+cq*xs**2):.4f}")
p(f"share between 2 units only above sqrt(2)*x* = {math.sqrt(2)*xs:.3f} of one unit's rating; k -> k+1 units above x*·sqrt(k(k+1))")
for X in [0.2, 0.3, 0.5, 0.8]:
    one = a + b*X + cq*X**2; two = 2*(a + b*X/2 + cq*(X/2)**2)
    p(f"   load {X:.1f}: one unit loses {one:.4f}, two sharing lose {two:.4f} -> {'one' if one < two else 'two'} better; input saving {(two-one)/(X+two):+.2%}")
p("Thomas-Krenn worked example (2x1000 W Platinum, 500 W load): 556 W shared vs 532 W one active ->", f"{1-532/556:.2%} less")

p("\n=== Lead 4: fans and the cube ===")
for d in [-0.1, -0.2, -0.3]:
    p(f"speed change {d:+.0%}: power x{(1+d)**3:.3f}; fused cube (1+d)^3-1 = 3d+3d^2+d^3 = {3*d:+.3f}{3*d*d:+.3f}{d**3:+.4f} = {(1+d)**3-1:+.4f}")
for N in [2, 4, 8]:
    for m, lab in [(3, "cube law"), (2.75, "Blasius turbulent duct (dp ~ Q^1.75)"), (2.0, "laminar (dp ~ Q)")]:
        r = (N/(N+1))**(m-1)
        p(f"   N={N}: run N+1 units at N/(N+1) speed instead of N at full, {lab}: power x{r:.3f} ({1-r:.1%} saving)")
p("excess airflow removed by slowing fans (cube) vs switching units off (linear), 30% less air:", f"{1-0.7**3:.1%} vs 30.0% fan saving")

p("\n=== Lead 5: setpoint raise compounds per degree (rung 3 power) ===")
for s in [0.04, 0.05]:
    for n in [3, 5, 9]:
        p(f"   {s:.0%}/°F, {n}°F: compounded {1-(1-s)**n:.1%} vs naive {s*n:.0%}")

p("\n=== Stacking independent savings composes by mult-add on the kept fraction ===")
sv = [0.025, 0.03, 0.03]
p("UPS eco 2.5%, cold-redundant PSUs 3%, fan sharing 3% ->", f"{1-np.prod([1-s for s in sv]):.2%}", "(naive", f"{sum(sv):.1%})", f"US-scale {US*(1-np.prod([1-s for s in sv])):.1f} TWh/yr")
open("/mnt/project-files/math-theory/datacenter-energy/checks_output.txt", "w").write("\n".join(out) + "\n")
