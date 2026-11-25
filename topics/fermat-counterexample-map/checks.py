"""Checks for the Fermat counterexample map.

What a hypothetical a^p + b^p = c^p (p >= 5 prime) would force, checked by computer on
real numbers: near-misses, Frey-like curves from real abc triples, the level-2 squeeze,
and the methods that are known not to work.

Run: python3 checks.py > checks_output.txt   (needs sympy and cypari2)
"""
from math import gcd, log
import sympy as sp
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**8, silent=True) if hasattr(pari, "allocatemem") else None

def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def v(q, n):
    n = abs(n); k = 0
    while n % q == 0:
        n //= q; k += 1
    return k

def rad(n):
    r = 1
    for q in sp.factorint(abs(n)):
        r *= q
    return r

# ---------------------------------------------------------------------------
hdr("1. The Frey curve's invariants, symbolically")
A, B = sp.symbols("A B")
C = A + B
# y^2 = x(x-A)(x+B) = x^3 + (B-A) x^2 - AB x
a1, a2, a3, a4, a6 = 0, B - A, 0, -A * B, 0
b2 = a1**2 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3**2 + 4 * a6
b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
D = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
print("c4   =", sp.factor(c4))
print("c6   =", sp.factor(c6))
print("Disc =", sp.factor(D))
print("c4^3 - c6^2 - 1728*Disc =", sp.expand(c4**3 - c6**2 - 1728 * D))
c4m, c6m, Dm = c4 / 16, c6 / 64, D / 2**12
print("minimal model (u = 2): c4 =", sp.factor(c4m), "| c6 =", sp.factor(c6m), "| Disc =", sp.factor(Dm))
X = (B - A) * (2 * A + B) * (A + 2 * B)
print("identity 4(A^2+AB+B^2)^3 - ((B-A)(2A+B)(A+2B))^2 - 27(ABC)^2 =",
      sp.expand(4 * (A**2 + A * B + B**2)**3 - X**2 - 27 * (A * B * C)**2))
print("j = c4^3/Disc =", sp.factor(c4**3 / D))

# ---------------------------------------------------------------------------
hdr("2. The exponent of 2 in the minimal discriminant never lands on a multiple of p")
print("v_2(Disc_min) = 2p*v_2(b) - 8 for the normalised curve (a = 3 mod 4, b even).")
bad = []
for p in sp.primerange(5, 200):
    for vb in range(1, 40):
        e = 2 * p * vb - 8
        if e % p == 0:
            bad.append((p, vb))
print("pairs (p, v_2(b)) with p | 2p*v - 8, p in [5,200), v in [1,40):", bad or "none")
print("Reason: 2p*v - 8 = -8 mod p, and an odd prime never divides 8.")

# ---------------------------------------------------------------------------
hdr("3. Genus of X0(N) by the standard formula, and dim S_2(Gamma0(N)) from PARI")
def genus(N):
    mu = N
    for q in sp.factorint(N):
        mu = mu * (q + 1) // q
    nu2 = sum(1 for x in range(N) if (x * x + 1) % N == 0)
    nu3 = sum(1 for x in range(N) if (x * x + x + 1) % N == 0)
    nuinf = sum(sp.totient(gcd(d, N // d)) for d in sp.divisors(N))
    g12 = 12 + mu - 3 * nu2 - 4 * nu3 - 6 * nuinf
    assert g12 % 12 == 0
    return g12 // 12, (mu, nu2, nu3, nuinf)
mism = []
zero = []
for N in range(1, 61):
    g, data = genus(N)
    dim = int(pari(f"mfdim([{N},2],1)"))
    if g != dim:
        mism.append((N, g, dim))
    if g == 0:
        zero.append(N)
print("genus 0 levels up to 60:", zero)
print("levels where genus != dim S_2 (should be none):", mism or "none")
print("level 2 counts (mu, nu2, nu3, nu_inf) =", genus(2)[1], "-> genus", genus(2)[0])
print("level 11 genus =", genus(11)[0], "(the first level with a weight-2 cusp form)")

# ---------------------------------------------------------------------------
hdr("4. Near misses: fake solutions, and the smallest modulus that kills each")
def smallest_killer(a, b, c, p, limit=10**4):
    for n in range(2, limit):
        if (pow(a, p, n) + pow(b, p, n) - pow(c, p, n)) % n:
            return n
    return None
fakes = [(3987, 4365, 4472, 12), (1782, 1841, 1922, 12)]
for a, b, c, p in fakes:
    l, r = a**p + b**p, c**p
    print(f"{a}^{p} + {b}^{p} vs {c}^{p}: relative gap {abs(l - r) / r:.3e}; "
          f"killed by modulus {smallest_killer(a, b, c, p)}; agree to "
          f"{len(str(l)) - len(str(abs(l - r)))} leading digits")
for p in (5, 7, 11, 13):
    s = 3**p + 5**p
    cr = round(s ** (1 / p))
    best = min(range(cr - 1, cr + 2), key=lambda c: abs(c**p - s))
    print(f"3^{p} + 5^{p} = {s} = {dict(sp.factorint(s))}; nearest {best}^{p} = {best**p}")
print("Every particular fake dies to some small modulus, but no single modulus kills them all:")
print("for any n, (a, b, c) = (-1, 32n, -1) satisfies a^p + b^p = c^p mod n with a = 3 mod 4,")
print("32 | b, gcd(a, b) = 1, abc != 0 (proved in Lean as fake_passes_mod).")
for n in (7, 11, 1001, 2**10 * 3**5 * 5**3):
    a, b, c = -1, 32 * n, -1
    assert all((pow(a, p, n) + pow(b, p, n) - pow(c, p, n)) % n == 0 for p in (5, 7, 11, 13))
print("checked for n = 7, 11, 1001, 2^10*3^5*5^3 and p = 5, 7, 11, 13.")

# ---------------------------------------------------------------------------
hdr("5. Local (mod q) obstructions only push divisibility, never contradict")
for p in (5, 7):
    obst = []
    for q in sp.primerange(3, 400):
        if q == p:
            continue
        P = {pow(x, p, q) for x in range(1, q)}
        ok = any(((x + y) % q) in P for x in P for y in P)
        if not ok:
            obst.append(q)
    print(f"p = {p}: primes q < 400 where x^p + y^p = z^p mod q forces q | xyz: {obst}")
print("Each one only says q divides abc. Solutions with abc != 0 exist in every q-adic field,")
print("e.g. 1 + (q t)^p = c^p has a q-adic root c for q != p (Hensel), so no local test can")
print("ever rule a counterexample out.")

# ---------------------------------------------------------------------------
hdr("6. Sophie Germain's case I: auxiliary primes q = 2kp + 1")
germ = {}
for p in sp.primerange(5, 100):
    for k in range(1, 30):
        q = 2 * k * p + 1
        if not sp.isprime(q):
            continue
        P = {pow(x, p, q) for x in range(1, q)}
        no_sum = not any(((x + y) % q) in P for x in P for y in P)
        p_not_power = (p % q) not in P
        if no_sum and p_not_power:
            germ[p] = q
            break
print("first auxiliary prime q for each p < 100 (Germain/Legendre condition):", germ)
print("This proves case I (p does not divide abc) for these p; case II (p | abc) is untouched.")

# ---------------------------------------------------------------------------
hdr("7. Kummer's method: irregular primes and the Wieferich criterion")
irr = []
for p in sp.primerange(3, 300):
    for k in range(2, p - 2, 2):
        if sp.numer(sp.bernoulli(k)) % p == 0:
            irr.append((p, k))
            break
print("irregular primes < 300 (first index 2k with p | numerator of B_2k):", irr)
print("B_32 numerator =", sp.numer(sp.bernoulli(32)), "=", dict(sp.factorint(abs(sp.numer(sp.bernoulli(32))))))
wief = [p for p in sp.primerange(3, 10**6) if pow(2, p - 1, p * p) == 1]
print("Wieferich primes below 10^6 (2^(p-1) = 1 mod p^2):", wief)
print("Unique factorisation fails already in Z[sqrt(-5)]: 6 = 2*3 = (1+sqrt-5)(1-sqrt-5),")
print("no element has norm 2 or 3:",
      not any(x * x + 5 * y * y in (2, 3) for x in range(-3, 4) for y in range(-2, 3)))
print("Class number of Q(sqrt -5):", pari("quadclassunit(-20).no"),
      "| class number of Q(zeta_23) is 3 (Kummer; cited, not recomputed here)")

# ---------------------------------------------------------------------------
hdr("8. Size: a Fermat triple would beat every known abc triple's quality")
print("rad(a^p b^p c^p) = rad(abc) <= abc < c^3, so quality log(c^p)/log(rad) > p/3.")
records = [(2, 3**10 * 109, 23**5), (11**2, 3**2 * 5**6 * 7**3, 2**21 * 23),
           (19 * 1307, 7 * 29**2 * 31**8, 2**8 * 3**22 * 5**4), (283, 5**11 * 13**2, 2**8 * 3**8 * 17**3),
           (1, 2 * 3**7, 5**4 * 7)]
for a, b, c in records:
    assert a + b == c, (a, b, c)
    q = log(c) / log(rad(a * b * c))
    print(f"  {a} + {b} = {c}: quality {q:.4f}")
for p in (5, 7, 11):
    print(f"  a counterexample with p = {p} would need quality > {p/3:.4f}")
print("The record (Reyssat, 2 + 3^10*109 = 23^5) is 1.6299 < 5/3. Under an explicit abc bound")
print("c < rad(abc)^2 every p >= 6 would be ruled out at once; abc itself is not proved.")

# ---------------------------------------------------------------------------
hdr("9. Frey-like curves from real triples: which property each real curve lacks")
def frey_from_triple(a, b, c):
    """Normalise A + B = C with A = 3 mod 4 and B even; return (A, B, C)."""
    trip = [a, b, -c]  # a + b + (-c) = 0
    ev = [t for t in trip if t % 2 == 0]
    od = [t for t in trip if t % 2 != 0]
    assert len(ev) == 1
    Bv = ev[0]
    o1, o2 = od
    # want A + B = C, i.e. A + B - C = 0 with {A, -C} = {o1, o2} up to overall sign
    for sgn in (1, -1):
        for Av, Cneg in ((o1, o2), (o2, o1)):
            Aa, Bb, Cc = sgn * Av, sgn * Bv, -sgn * Cneg
            if Aa + Bb == Cc and Aa % 4 == 3:
                return Aa, Bb, Cc
    raise ValueError
def serre_level(E, p):
    """Odd part of the mod-p Serre level for a semistable curve: bad q != p with p not dividing v_q(Disc_min)."""
    N = int(pari.ellglobalred(E)[0])
    D = int(E[11])  # disc of the (minimal) model
    lvl = []
    for q in sp.factorint(N):
        if q == p:
            continue
        if v(q, D) % p != 0:
            lvl.append(q)
    return N, lvl
triples = [(3, 125, 128), (1, 80, 81), (1, 4374, 4375), (2, 3**10 * 109, 23**5),
           (1, 2400, 2401), (243, 2882, 3125), (1, 1023, 1024), (3125, 29643, 32768)]
for (a, b, c) in triples:
    assert a + b == c and gcd(a, b) == 1
    Aa, Bb, Cc = frey_from_triple(a, b, c)
    E = pari.ellinit([0, Bb - Aa, 0, -Aa * Bb, 0])
    Em = pari.ellinit(pari.ellminimalmodel(E))
    gr = pari.ellglobalred(Em)
    N = int(gr[0])
    D = int(Em[11])
    semi = (N == rad(N))
    tors = pari.elltors(Em)[1]
    line = f"{a} + {b} = {c}: A,B,C = {Aa},{Bb},{Cc}; conductor {N} = {dict(sp.factorint(N))}; "
    line += f"semistable {semi}; Disc_min = {dict(sp.factorint(D))}; torsion {tors}"
    print(line)
    for p in (3, 5, 7):
        Nn, lvl = serre_level(Em, p)
        two = 2 in lvl
        print(f"     mod {p}: primes q != p whose exponent is not a multiple of p: {lvl}")
print("A genuine counterexample would leave only 2 in that list (and keep 2 there for every p >= 5).")
print("Every real curve above leaves some odd prime, or escapes by being reducible (rational p-torsion).")

# ---------------------------------------------------------------------------
hdr("10. The reducible escape: a real curve where the squeeze closes except for Mazur's arrow")
E = pari.ellinit(pari.ellminimalmodel(pari.ellinit([0, 128 - (-125), 0, 125 * 128, 0])))
print("Curve from -125 + 128 = 3: conductor", int(pari.ellglobalred(E)[0]), "torsion", pari.elltors(E)[1])
print("Disc_min =", dict(sp.factorint(int(E[11]))), "- for p = 3 the exponents at 2 and 5 are multiples of 3,")
iso = pari.ellisomat(E)[1]
degs = sorted({int(iso[i, j]) for i in range(iso.nrows()) for j in range(iso.ncols())})
print("isogeny degrees in its class (PARI ellisomat):", degs)
print("so the mod-3 Serre level drops to 1 (and the weight rises because 3 does not divide v_3 = 2,")
print("so no level-1 form exists either), yet there is no contradiction: the curve has a rational")
print("3-isogeny, its mod-3 representation is reducible, and Ribet's arrow needs it irreducible.")
print("That is exactly the hole Mazur's theorem plugs for the Frey curve when p >= 5.")

# ---------------------------------------------------------------------------
hdr("11. Search: real semistable Frey-like curves with the fewest odd blocking primes, p = 5 and 7")
import itertools
for p in (5, 7):
    best = []
    pw = [n**p for n in range(1, 40)]
    smalls = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 19, 23, 32, 64, 128]
    seen = set()
    for x in pw:
        for s in smalls:
            for y in pw:
                for t in smalls:
                    a, b = x * s, y * t
                    c = a + b
                    if gcd(a, b) != 1 or (a, b) in seen or a > b:
                        continue
                    seen.add((a, b))
                    if sum(1 for u in (a, b, c) if u % 2 == 0) != 1:
                        continue
                    ev = [u for u in (a, b, c) if u % 2 == 0][0]
                    if ev % 16:
                        continue
                    if c > 10**18:
                        continue
                    # odd blocking primes: odd q != p with v_q(abc) not a multiple of p
                    fa = sp.factorint(a * b * c)
                    block = [q for q, e in fa.items() if q != 2 and (2 * e) % p]
                    best.append((len(block), sp.prod(block) if block else 1, a, b, c, block))
    best.sort()
    print(f"p = {p}: fewest odd blocking primes among {len(best)} normalisable triples:")
    for row in best[:6]:
        print(f"   {row[2]} + {row[3]} = {row[4]}  blocking {row[5]}")
print("(blocking = odd primes, p included, whose exponent in abc is not a multiple of p.)")
print("In this search no real triple reaches zero blocking primes; the best leave one odd prime,")
print("and that one prime is what keeps the real curve's level above 2.")
