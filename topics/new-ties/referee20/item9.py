"""Item 9: W = N K l K K m^M N (five letters, items K, lK, K, m^M N)."""
from fractions import Fraction as F
def T(items):
    Tv, A = 1, 0
    for (L, Mm, y) in items:
        A += L
        Tv = 19 ** A + (19 if y == 'N' else 10) * 200 ** L * 20 ** Mm * Tv
    return Tv
def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v
# independent check by composing the affine maps
MAPS = {'l': (F(19, 200), F(0)), 'm': (F(1, 20), F(0)), 'N': (F(1, 19), F(1, 19)), 'K': (F(1, 10), F(1, 10))}
def psi(word):
    s, t = F(1), F(0)
    for ch in word:
        sx, tx = MAPS[ch]; s, t = s * sx, s * tx + t
    return t / s
T3 = T([(0, 0, 'K'), (1, 0, 'K'), (0, 0, 'K')])
assert psi('NKlKK') == F(T3, 19)
print('T(NKlKK) =', T3, ' mod 361 =', T3 % 361, '(-1 mod 361 is 360)')
for Mm in range(0, 6):
    w = 'NKlKK' + 'm' * Mm + 'N'
    Tw = T([(0, 0, 'K'), (1, 0, 'K'), (0, 0, 'K'), (0, Mm, 'N')])
    assert Tw == 19 + 19 * 20 ** Mm * T3
    assert psi(w) == F(Tw, 19)
# 19-adic mu: 20^mu = -1/T3 in Z_19, computed digit by digit mod 19^J
J = 12
mod = 19 ** J
c = (-pow(T3, -1, mod)) % mod
mu = 0
for j in range(1, J):
    # want 20^mu = c mod 19^(j+1); 20^x mod 19^(j+1) depends on x mod 19^j
    while pow(20, mu, 19 ** (j + 1)) != c % 19 ** (j + 1):
        mu += 19 ** (j - 1)
print('mu mod 19^%d =' % (J - 1), mu, ' v19(mu) =', vp(mu, 19))
first = {}
for Mm in range(0, 200000):
    Tw = 19 + 19 * 20 ** Mm * T3 if Mm < 50 else None
    if Mm < 50:
        v = vp(Tw, 19)
    else:
        v = vp((19 + 19 * pow(20, Mm, 19 ** 12) * T3) % 19 ** 12, 19)
    assert v == 2 + min(vp(Mm - mu, 19) if Mm != mu else 99, 10), Mm
    if v not in first:
        first[v] = Mm
print('v19(T) = 2 + v19(M - mu) verified for M < 200000; first M with v19 = 3,4,5,6:', {k: first.get(k) for k in (3, 4, 5, 6)})
