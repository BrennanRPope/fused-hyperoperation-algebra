"""Brute force directly from the affine maps.
1. map facts: ml=lm, mK=lN; leading l/m keep psi; trailing l/m scale psi by 2,5,19 units;
   first letter K ~ N; psi(NN)=psi(Nm)=20.
2. psi from map composition == T/19^L and T3/19^(L1+L2) closed forms.
3. small searches: three-letter fingerprint 1 (16.1); two-letter collisions (16.2);
   three-vs-two collisions (17.1).
4. raw-word search: all words over {l,m,N,K} up to length 8, reduced by the map facts,
   fingerprints compared by number of N/K letters of the reduced word.
"""
import itertools, sys
from fractions import Fraction as F
from common import *

# ---- 1. facts
assert compose_word('ml') == compose_word('lm')
assert compose_word('mK') == compose_word('lN')
assert psi_word('NN') == psi_word('Nm') == 20
assert psi_word('K') == psi_word('N') == 1
for w in ['NlN', 'NmK', 'NllmN', 'NKmN']:
    for pre in ['l', 'm', 'lm', 'mml']:
        assert psi_word(pre + w) == psi_word(w)
    for suf in ['l', 'm', 'lm', 'mml']:
        q = psi_word(w + suf) / psi_word(w)
        assert kappa_of_fraction(q) == 1
print('map facts ok')

# ---- 2. closed forms
for (L, M, x) in two_letter_words(8):
    assert psi_word(word_string2(L, M, x)) == F(T2(L, M, x), 19 ** L)
cnt = 0
for w in three_letter_words(5, 5):
    assert psi_word(word_string3(w)) == F(T3(*w), 19 ** (w[0] + w[3]))
    cnt += 1
print('closed forms ok (%d three-letter words)' % cnt)

# also check g = v19(T) <= 1 for two-letter words, and g = [x=N and L>=1]
for (L, M, x) in two_letter_words(60):
    k, g, T = kappa2(L, M, x)
    assert g == (1 if (x == 'N' and L >= 1) else 0)
print('two-letter g = [x=N, L>=1] verified up to H=60')

# ---- 3. small searches on reduced words (closed forms, verified above)
H2max = int(sys.argv[1]) if len(sys.argv) > 1 else 40
H3max = int(sys.argv[2]) if len(sys.argv) > 2 else 14
fp2 = {}
for (L, M, x) in two_letter_words(H2max):
    k, g, T = kappa2(L, M, x)
    assert k != 1
    if k in fp2:
        print('16.2 COLLISION', fp2[k], (L, M, x))
    fp2[k] = (L, M, x)
print('16.2: %d two-letter words with H<=%d, fingerprints distinct' % (len(fp2), H2max))
n3 = 0
for w in three_letter_words(H3max, H3max):
    k = strip(T3(*w))
    n3 += 1
    if k == 1:
        print('16.1 COUNTEREXAMPLE', w)
    if k in fp2:
        print('17.1 COLLISION', w, fp2[k])
print('16.1/17.1: %d three-letter words with H1,H2<=%d: none with kappa 1, none collide with two-letter words H<=%d'
      % (n3, H3max, H2max))

# ---- 4. raw words up to length 8 (independent of the reduced-word bookkeeping)
def reduce_raw(word):
    """drop leading l/m, trailing l/m; first K->N; rewrite mK->lN, ml->lm; NN->Nm.
    returns (canonical string, number of N/K letters)"""
    w = word.lstrip('lm').rstrip('lm')
    if not w:
        return '', 0
    if w[0] == 'K':
        w = 'N' + w[1:]
    changed = True
    while changed:
        changed = False
        for a, b in (('mK', 'lN'), ('ml', 'lm')):
            if a in w:
                w = w.replace(a, b, 1)
                changed = True
        if w.startswith('NN'):
            w = 'Nm' + w[2:]
            changed = True
        w = w.rstrip('lm')
    return w, sum(ch in 'NK' for ch in w)

byfp = {}
Lmax = 8
for length in range(1, Lmax + 1):
    for tup in itertools.product('lmNK', repeat=length):
        word = ''.join(tup)
        if 'N' not in word and 'K' not in word:
            continue
        can, k = reduce_raw(word)
        kp = kappa_of_fraction(psi_word(word))
        assert kp == kappa_of_fraction(psi_word(can)), (word, can)
        byfp.setdefault(kp, {}).setdefault(k, set()).add(can)
bad = 0
for kp, d in byfp.items():
    if kp == 1:
        others = {k: v for k, v in d.items() if k >= 3}
        if others:
            print('raw: kappa 1 with >=3 letters', others); bad += 1
    if 2 in d and len(d[2]) > 1:
        print('raw: two-letter collision', d[2]); bad += 1
    if 2 in d and 3 in d and kp != 1:
        print('raw: 3 vs 2 collision', d); bad += 1
print('raw words up to length %d: %d fingerprints, problems: %d' % (Lmax, len(byfp), bad))
