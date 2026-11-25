"""
Rulecraft x Fused Hyperoperation Algebra: every number on the page.

Run:  python3 checks.py > checks_output.txt
Mult-add:  a (+) b = a + b + ab = (1+a)(1+b) - 1
Weighted:  a (+)_A b = a + b + A*ab      (A = coupling weight of the corner)
Inverse:   (-)a = -a/(1+a)               (undefined only at a = -1)
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
import sympy as sp

def ma(a, b): return a + b + a*b
def ma_list(xs):
    out = 0
    for x in xs: out = ma(out, x)
    return out
def inv(a): return -a/(1+a)

print("="*70)
print("0. Algebra identities (symbolic check with sympy)")
a, b, c, A, d, s, e1, e2, i1, i2, k = sp.symbols('a b c A d s e1 e2 i1 i2 k')
M = lambda x, y: x + y + x*y
checks = {
  "product form       a(+)b = (1+a)(1+b)-1": sp.expand(M(a,b) - ((1+a)*(1+b)-1)),
  "commutative        a(+)b = b(+)a": sp.expand(M(a,b) - M(b,a)),
  "associative": sp.expand(M(M(a,b),c) - M(a,M(b,c))),
  "-1 absorbs         (-1)(+)b = -1": sp.expand(M(-1,b) + 1),
  "inverse            a(+)(-a/(1+a)) = 0": sp.simplify(M(a, -a/(1+a))),
  "thermostat residue d(+)(-d) = -d^2": sp.expand(M(d,-d) + d**2),
  "weighted conjugacy A*(a(+)_A b) = (Aa)(+)(Ab)": sp.expand(A*(a+b+A*a*b) - M(A*a, A*b)),
  "balance            (s+d)(s-d) = s^2 - d^2": sp.expand((s+d)*(s-d) - (s**2-d**2)),
  "panel grid         (e1+e2)(i1+i2) = own + cross": sp.expand((e1+e2)*(i1+i2) - (e1*i1+e2*i2+e1*i2+e2*i1)),
  "fused voice price  (1+k)^2-1 = k^2+2k": sp.expand((1+k)**2-1-(k**2+2*k)),
}
for name, val in checks.items():
    print(f"  {'OK ' if val == 0 else 'FAIL'} {name}")

print()
print("Boundary theorem: mult-add over S++T = X_S (+) X_T, cross term = straddling corners")
xs = sp.symbols('x1:6')
for split in (1, 2, 3):
    S, T = xs[:split], xs[split:]
    XS = sp.expand(prod(1+x for x in S) - 1); XT = sp.expand(prod(1+x for x in T) - 1)
    total = sp.expand(prod(1+x for x in xs) - 1)
    straddle = sum(prod(sub) for r in range(1, 6) for sub in combinations(xs, r)
                   if any(x in S for x in sub) and any(x in T for x in sub))
    print(f"  |S|={split}, |T|={5-split}: total = XS+XT+XS*XT? {sp.expand(total-(XS+XT+XS*XT))==0};"
          f"  XS*XT = sum of straddling corners? {sp.expand(XS*XT - straddle)==0};"
          f"  straddling corners = {(2**split-1)*(2**(5-split)-1)} of {2**5-1}")

print()
print("="*70)
print("1. SCOPE MATCHES SPILLOVER: three towns on a river")
# Each town's treatment plant multiplies water quality at its own point and
# everywhere downstream by (1+q). Town j's gain = VALUE * (prod_{i<=j}(1+q_i) - 1).
q, VALUE, COST = F(2,5), 10, 6
towns = ["Upper", "Middle", "Lower"]
def town_gain(j, built):          # built: set of town indices with plants
    return VALUE * (prod((1+q) if i in built else 1 for i in range(j+1)) - 1)
def basin_gain(built): return sum(town_gain(j, built) for j in range(3))
print(f"  q = {q} (each plant lifts quality 40% at and below it), value {VALUE} per unit, plant cost {COST}")
print("  Own-gain test (each town counts only its own water):")
for j in range(3):
    for up in ([], [0], [0,1]) if j == 2 else ([], [0]) if j == 1 else ([],):
        up = [u for u in up if u < j]
        g = town_gain(j, set(up)|{j}) - town_gain(j, set(up))
        print(f"    {towns[j]:6s} plant, upstream built {[towns[u] for u in up]}: own gain {float(g):5.2f} vs cost {COST} -> {'build' if g>COST else 'skip'}")
alone = set()   # Upper skips (4<6), Middle skips (4<6 with no Upper), Lower skips (4<6 with none)
print(f"  Towns deciding alone: nothing gets built. Basin net = 0")
allb = {0,1,2}
print(f"  Basin body (all three): gain {float(basin_gain(allb)):.2f} - cost {3*COST} = net {float(basin_gain(allb)-3*COST):+.2f}")
for j in range(3):
    print(f"    {towns[j]} gets {float(town_gain(j, allb)):.2f}")
# Lower town's gain as a cube (layers by how many plants a corner multiplies)
lay = [comb(3,r)*q**r for r in range(1,4)]
print(f"  Lower town's gain on the cube (1.4^3 - 1 = {float((1+q)**3-1)}):"
      f" layer1 {float(VALUE*lay[0]):.2f}, layer2 {float(VALUE*lay[1]):.2f}, layer3 {float(VALUE*lay[2]):.2f}")
own_D = VALUE*(q + 2*q*q + q**3)     # corners containing Lower's own plant
print(f"    corners that contain Lower's own plant: {float(own_D):.2f};"
      f" corners owed to upstream plants only: {float(VALUE*(2*q+q*q)):.2f}")
# Upper plant's full value to the basin vs what Upper sees
up_full = basin_gain(allb) - basin_gain({1,2})
print(f"  Upper plant, everyone built: worth {float(up_full):.2f} to the basin; Upper sees {float(town_gain(0,{0})):.2f}."
      f" Spillover share {float(1 - town_gain(0,{0})/up_full):.0%}")
# boundary {Upper} | {Middle, Lower}: cross term
XS, XT = q, ma(q, q)
print(f"  Boundary Upper | Middle+Lower at Lower's tap: XS={float(XS)}, XT={float(XT)}, cross XS*XT={float(XS*XT):.3f}"
      f" = {float(XS*XT/ma(XS,XT)):.0%} of the effect")

print()
print("  Taxi ride vs a carbon-heavy shipment (Harsanyi corners)")
taxi = {("R",): 0, ("D",): 0, ("R","D"): 12}     # rider values ride at 30, driver's cost 18
print(f"   Taxi: corners {taxi}; every corner inside the pair -> no spillover; a cap removing 1000 rides destroys {12*1000}")
for surplus in (30, 15):
    damage = 18
    pair_sees, world_net = surplus, surplus - damage
    print(f"   Shipment, surplus {surplus}, world corner {{F,C,World}} = -{damage}: pair sees {pair_sees} -> trades;"
          f" world net {world_net:+d} -> {'good trade' if world_net>0 else 'bad trade made anyway'}")
print("   A price or cap of 18 per shipment puts the -18 corner inside the pair's view: trade only if surplus > 18")

print()
print("  Counting corners: a body with s of n parties sees 2^s-1 corners;"
      " (2^s-1)(2^(n-s)-1) straddle its edge")
for n_, s_ in ((3,1),(3,2),(10,3),(50,10)):
    print(f"    n={n_:3d}, s={s_:3d}: inside {2**s_-1}, straddling {(2**s_-1)*(2**(n_-s_)-1)}, outside-only {2**(n_-s_)-1}")

print()
print("="*70)
print("2. WISDOM AS AN AREA, AND MIXED PANELS")
def fused(e, i): return e + i + e*i
for name, (e_, i_) in {"expert": (1,9), "resident": (9,1), "balanced": (5,5), "no experience": (0,10)}.items():
    print(f"  {name:13s} e={e_}, i={i_}: strips e+i = {e_+i_}, wisdom corner e*i = {e_*i_}, fused total = {fused(e_,i_)}")
panels = {"expert+resident": [(1,9),(9,1)], "two experts": [(1,9),(1,9)], "two residents": [(9,1),(9,1)]}
for name, mem in panels.items():
    E = sum(m[0] for m in mem); I = sum(m[1] for m in mem)
    own = sum(m[0]*m[1] for m in mem); full = E*I
    print(f"  {name:16s}: E={E}, I={I}, area {full} = own corners {own} + cross corners {full-own}"
          f"   (balance: 10^2 - {abs(E-I)//2}^2 = {100-(abs(E-I)//2)**2})")
print("  Joining weight A (how well deliberation joins members): wisdom = own + A*cross")
for A_ in (0, 0.1, 0.5, 0.7, 1):
    print(f"    A={A_:<4}: expert+resident {18+82*A_:6.1f}   two experts {18+18*A_:5.1f}")

print()
print("="*70)
print("3. THE EVENT HORIZON AT -1, THE THERMOSTAT AND COMPOUNDING")
for loss in (F(1,5), F(1,2), F(9,10), F(99,100), F(999,1000)):
    print(f"  lose {float(loss):6.1%} -> exact way back {float(inv(-loss)):9.1%}")
print("  -1 absorbs: (-1)(+)b =", ma(-1, F(7,3)), "for any b; no way back exists")
print("  Additive thermostat (undo +d with -d) leaves -d^2:")
for dd in (F(1,10), F(1,4), F(1,2)):
    print(f"    drift {float(dd):.0%}: additive fix leaves {float(ma(dd,-dd)):+.2%}; exact fix {float(inv(dd)):+.2%} leaves {float(ma(dd,inv(dd))):+.0%}")
r, g = F(5,100), F(3,100)
rg = ma(r, inv(g))
print(f"  Debt trap: r=5%, g=3%: naive r-g = 2.000%, exact r(-)g = (r-g)/(1+g) = {float(rg):.4%}")
print(f"    debt 100% of GDP needs a primary surplus of {float(rg):.3%} of GDP to hold steady")
n = 30
tot = (1+rg)**n - 1
layers = [comb(n,j)*rg**j for j in range(1, n+1)]
print(f"    30 years, no surplus: debt ratio x{float(1+tot):.4f} -> {float(100*(1+tot)):.1f}% of GDP")
print("    the 30-fold mult-add on the cube, by layer: " +
      ", ".join(f"L{j+1} {float(100*layers[j]):.2f}" for j in range(5)) +
      f", rest {float(100*sum(layers[5:])):.3f}")
print(f"    linear layer {float(100*layers[0]):.2f} points; compounding layers (2 and up) {float(100*sum(layers[1:])):.2f} points")
nstar = 1 + 2/rg
print(f"    layer 2 overtakes layer 1 after n = 1 + 2/x = {float(nstar):.0f} years (L2/L1 = (n-1)x/2)")

print()
print("="*70)
print("4. QUICK MATCHES")
print("  Order of reforms: mult-add commutes, so percent reforms on one quantity are order-free:",
      ma(F(1,10), F(1,5)) == ma(F(1,5), F(1,10)))
print("  Fused exponent rung (1+a)^(1+b)-1 does not commute: a=1,b=2 ->", 2**3-1, "; a=2,b=1 ->", 3**2-1)
print(f"    with a=10%, b=20%: {1.1**1.2-1:.4%} vs {1.2**1.1-1:.4%}")
print("  Voice pricing with a free floor vote (value u, price 0.5 per credit):")
for u in (1, 2, 4):
    k_plain = u            # plain k^2 on extras: marginal 2k*0.5 = k = u
    V_fused = max(1, u)    # fused (1+k)^2-1 on total V=1+k: marginal 2V*0.5 = V = u
    print(f"    value {u}: plain -> total voice {1+k_plain} (cost {k_plain**2});"
          f" fused -> total voice {V_fused} (cost {V_fused**2-1})")
print("    plain pricing gives voice 2:3:5 for values 1:2:4; fused pricing gives 1:2:4")
print("  Saved influence that earned any return would compound: 10 periods at 5% ->",
      f"{1.05**10:.3f}x, layers 1..3: {10*0.05:.3f}, {45*0.0025:.4f}, {120*0.000125:.5f}")
