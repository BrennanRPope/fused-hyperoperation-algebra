"""Round 15 (DRAFT): the fingerprint kappa(psi(w)) of every reduced word up to n letters, found exhaustively in exact
integer arithmetic.

psi(w) = t_w/s_w (Round 14). Write psi(w) = num / 19^e with num an integer; the fingerprint kappa(w) is num with every
factor 2, 5 and 19 removed. A path off the line with words P and Q comes back to the line exactly when kappa(P) = kappa(Q)
(Lemma 15.1). Front relations (Proposition 15.2): psi(K w) = psi(N w); psi(NN w) = psi(KN w) = psi(Nm w) = psi(Km w);
mK = lN and ml = lm as maps; leading l, m do not change psi; trailing l, m multiply psi by a unit. canon() applies them.
Conjecture 15.3: two reduced words have the same fingerprint only if they have the same canonical form.
Counts are cumulative over lengths 1..L and cover words ending in N or K (trailing l, m only multiply psi by a unit).
Conjecture 15.3 holds to n letters exactly when distinct canonical forms = distinct fingerprints at every length.
Usage: python3 kappa15.py [n=13]  (n = 14: 29 million words, about a minute, a few GB of memory). Writes kappa15.json.
"""
import sys, re, json
from collections import defaultdict, Counter
n = int(sys.argv[1]) if len(sys.argv) > 1 else 13
AFTER = {'l': 'lmNK', 'm': 'mN', 'N': 'lmNK', 'K': 'lmNK'}
def kappa_int(a):
    for p in (2, 5, 19):
        while a % p == 0: a //= p
    return a
def canon(w):
    w = w.lstrip('lm')
    w = 'N' + w[1:]
    if len(w) > 1 and w[1] == 'N': w = 'Nm' + w[2:]
    while True:
        v = w.replace('mK', 'lN').replace('ml', 'lm')
        if v == w: break
        w = v
    return w.rstrip('lm')
def shape(w): return re.sub(r'[lm]+', '.', w)
if __name__ == "__main__":
    first = {}            # fingerprint -> canonical form seen first
    clash = []; words = 0; per_len = []; unit_words = Counter(); classes_by_len = []; raw = Counter(); nk_by_len = []; canon_by_len = []
    canons = set(); nk = 0
    ndk = Counter()   # words ending in N or K whose psi is a unit, by number of N/K letters
    lev = [('N', 1, 0), ('K', 1, 0)]       # (word, num, e) with psi = num / 19^e
    for L in range(1, n + 1):
        for w, num, e in lev:
            words += 1
            if w[-1] not in 'NK': continue
            k = kappa_int(num); c = canon(w); raw[k] += 1; nk += 1; canons.add(c)
            if k == 1: unit_words[shape(w)] += 1; ndk[sum(1 for ch in w if ch in 'NK')] += 1
            if k in first:
                if first[k] != c: clash.append((first[k], c, k))
            else: first[k] = c
        per_len.append(len(lev)); classes_by_len.append(len(first)); nk_by_len.append(nk); canon_by_len.append(len(canons))
        if L == n: break
        p19 = None; new = []
        for w, num, e in lev:
            z = 19 ** e
            for x in AFTER[w[-1]]:
                if x == 'l': new.append((w + x, num * 200, e + 1))
                elif x == 'm': new.append((w + x, num * 20, e))
                elif x == 'N': new.append((w + x, num * 19 + z, e))
                else: new.append((w + x, num * 10 + z, e))
        lev = new
        print("length", L + 1, "words", len(lev), flush=True)
    out = dict(n=n, words=words, words_per_length=per_len, distinct_fingerprints_by_length=classes_by_len,
               words_ending_in_N_or_K_by_length=nk_by_len, distinct_canonical_forms_by_length=canon_by_len,
               most_words_sharing_one_fingerprint=max(raw.values()), fingerprints_shared_by_several_words=sum(1 for v in raw.values() if v > 1), clashes_after_canonical_form=len(clash), clash_examples=clash[:20], unit_psi_shapes=dict(unit_words), unit_psi_by_NK_letters=dict(ndk))
    json.dump(out, open("kappa15.json", "w"), indent=1)
    print("n", n, "words", words, "distinct fingerprints", len(first), "clashes after canonical form", len(clash))
    print("words ending in N or K whose psi is a {2,5,19}-unit, by shape (l, m runs shown as '.'):", dict(unit_words))
