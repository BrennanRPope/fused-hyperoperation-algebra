"""
Consensus through Fused Hyperoperation Algebra: numeric checks.

Every claim on the page that has a number comes from here.  Run:
    python3 checks.py > checks_output.txt
Seeds are fixed, so the output is repeatable.

Notation
    the circle operation        a (+) b = a + b + a*b,       1 + (a (+) b) = (1+a)(1+b)
    weighted        a (+)k b = a + b + k*a*b     k = -1 is "either one happens"
    circle power  a^(+n) = (1+a)^n - 1
"""
import math
import numpy as np

rng = np.random.default_rng(20261125)


def madd(a, b, k=1.0):
    return a + b + k * a * b


def mpow(a, n, k=1.0):
    """n-fold circle operation with dial k of a with itself."""
    x = 0.0
    for _ in range(n):
        x = madd(x, a, k)
    return x


def hdr(s):
    print()
    print("=" * 72)
    print(s)
    print("=" * 72)


# ----------------------------------------------------------------------------
hdr("C1  Proof of work: catching up z blocks is a circle power")
# Gambler's ruin: attacker with share q, honest p = 1-q, starting z behind,
# ever catches up with probability (q/p)^z when q < p, else 1.
# Write r = q/p and a = r - 1.  Then r^z = 1 + a^(+z): the depth is a
# Circle power of a, and q = 1/2 makes a = 0, the circle-operation identity.
print("q      r=q/p   a=r-1    r^10 (gambler)   1 + mpow(a,10)")
for q in [0.10, 0.25, 0.40, 0.45, 0.49, 0.50]:
    p = 1 - q
    r = q / p
    a = r - 1
    print(f"{q:.2f}  {r:7.4f}  {a:+.4f}   {r**10:.6e}     {1 + mpow(a, 10):.6e}")


def nakamoto(q, z):
    """Satoshi's whitepaper formula (Poisson approximation)."""
    p = 1 - q
    lam = z * q / p
    s = 1.0
    for k in range(z + 1):
        poi = math.exp(-lam) * lam**k / math.factorial(k)
        s -= poi * (1 - (q / p) ** (z - k))
    return s


def rosenfeld(q, z):
    """Exact double-spend probability (Rosenfeld 2014), negative binomial race."""
    p = 1 - q
    if q >= p:
        return 1.0
    s = 0.0
    for m in range(z + 1):
        s += math.comb(m + z - 1, m) * (p**z * q**m - p**m * q**z)
    return 1 - s


def sim_double_spend(q, z, trials=200_000, cap=400):
    """Monte Carlo: honest network mines z blocks; attacker mines in secret
    meanwhile, then races. Success once the attacker draws level (Rosenfeld)."""
    wins = 0
    # number of attacker blocks while honest mines z: negative binomial
    m = rng.negative_binomial(z, 1 - q, size=trials)
    deficit = z - m  # attacker needs to overtake: deficit+1 more than honest
    for d in deficit:
        if d < 0:
            wins += 1
            continue
        gap = d  # Rosenfeld's convention: a tie is enough (attacker pre-mined one block)
        steps = 0
        while 0 < gap < cap and steps < 20_000:
            gap += -1 if rng.random() < q else 1
            steps += 1
        if gap <= 0:
            wins += 1
    return wins / trials


print("\nDouble-spend success: Nakamoto, exact (Rosenfeld), Monte Carlo")
print("q     z   Nakamoto     exact       simulated")
for q in [0.10, 0.25, 0.40]:
    for z in [1, 3, 6]:
        print(f"{q:.2f}  {z}   {nakamoto(q, z):.5f}     {rosenfeld(q, z):.5f}     {sim_double_spend(q, z, 40_000):.5f}")
print("At q = 0.5 the exact value is 1 at every depth:", [rosenfeld(0.5, z) for z in (1, 6, 60)])

# ----------------------------------------------------------------------------
hdr("C2  BFT voting: why 1/3, as worst-case overlap of two quorums")
# Two quorums of size Q in n validators overlap in at least 2Q - n.
# Safety needs the overlap to contain an honest validator: 2Q - n > f.
# Liveness needs a quorum without the f faulty: Q <= n - f.
# Both at once  <=>  n > 3f.  As fractions: f < 1/3.
for n in [4, 7, 10, 100]:
    fmax = max(f for f in range(n) if 2 * (n - f) - n > f)
    print(f"n={n:4d}: largest tolerable f = {fmax}  (n > 3f gives {(n - 1) // 3})")
# k quorums all sharing an honest validator: overlap >= kQ - (k-1)n > f
print("k quorums must share an honest node: tolerable fraction 1/(k+1)")
for k in [1, 2, 3, 4]:
    print(f"  k={k}: f < 1/{k + 1} = {1 / (k + 1):.4f}")
# Random quorums: overlap of two random fraction-q subsets is q*q (product),
# their union is q (+)_{-1} q = 2q - q^2: the k = -1 circle operation.
for q in [0.5, 0.67, 0.8]:
    print(f"random quorums of fraction {q}: union {madd(q, q, -1):.4f}, overlap {q * q:.4f}, worst-case overlap {max(0, 2 * q - 1):.4f}")

# ----------------------------------------------------------------------------
hdr("C3  Committee sampling: failure is the mass of the top layers of a cube")
# A committee of c validators drawn from a population with faulty fraction f.
# Each member is one axis of a c-cube with sides (1-f, f).  Layer j holds the
# corners with j faulty members: C(c,j) f^j (1-f)^(c-j) (Pascal-weighted).
# The committee fails when j >= c/3.
from math import comb


def committee_fail(c, f, frac=1 / 3):
    from scipy.stats import binom
    jmin = math.ceil(frac * c)
    return float(binom.sf(jmin - 1, c, f))  # sum of layers j >= jmin


print("c      f=0.20          f=0.25          f=0.30   (exact layer mass)")
for c in [10, 30, 100, 300, 1000, 2000]:
    print(f"{c:5d}  {committee_fail(c, .2):.3e}      {committee_fail(c, .25):.3e}      {committee_fail(c, .3):.3e}")
c, f = 100, 0.25
sims = rng.binomial(c, f, size=1_000_000)
print(f"Monte Carlo c=100 f=0.25: {np.mean(sims >= math.ceil(c / 3)):.4e} vs exact {committee_fail(c, f):.4e}")
# majority of 3: layers 1,3,3,1 -> 3p^2 - 2p^3, fixed point 1/2
maj3 = lambda p: p**3 + 3 * p**2 * (1 - p)
print("majority-of-3 map (layers 1,3,3,1):", {p: round(maj3(p), 4) for p in (0.2, 0.4, 0.5, 0.6)})
# hierarchical delegation: iterate the majority map (two-level delegated stake)
for p in (0.30, 0.45, 0.55):
    seq = [p]
    for _ in range(4):
        seq.append(maj3(seq[-1]))
    print(f"  random spread, start {p}: levels -> {[round(x, 4) for x in seq]}")
print("  worst case (attacker places stake): 2 levels of bare majority need only just over (1/2)^2 = 1/4")

# ----------------------------------------------------------------------------
hdr("C4  Selfish mining (Eyal & Sirer 2014): reproduce the threshold")


def selfish_mining(alpha, gamma, blocks=600_000):
    """State machine from Eyal & Sirer. Returns attacker revenue share."""
    lead = 0
    branch_race = False
    a_rev = h_rev = 0
    r = rng.random(blocks)
    g = rng.random(blocks)
    for i in range(blocks):
        if r[i] < alpha:  # attacker finds a block
            if branch_race:  # state 0': attacker wins the race
                a_rev += 2
                branch_race = False
                lead = 0
            else:
                lead += 1
        else:  # honest finds a block
            if branch_race:  # 0': honest resolve race
                if g[i] < gamma:
                    a_rev += 1
                    h_rev += 1
                else:
                    h_rev += 2
                branch_race = False
                lead = 0
            elif lead == 0:
                h_rev += 1
            elif lead == 1:
                branch_race = True  # publish, race
                lead = 0
            elif lead == 2:
                a_rev += 2
                lead = 0
            else:
                a_rev += 1
                lead -= 1
    return a_rev / (a_rev + h_rev)


def es_formula(a, g):
    num = a * (1 - a) ** 2 * (4 * a + g * (1 - 2 * a)) - a**3
    den = 1 - a * (1 + (2 - a) * a)
    return num / den


print("alpha  gamma  simulated  Eyal-Sirer formula  (honest share = alpha)")
for g in [0.0, 0.5]:
    for a in [0.20, 0.25, 0.30, 0.33, 0.35, 0.40]:
        print(f"{a:.2f}   {g:.1f}    {selfish_mining(a, g, 300_000):.4f}     {es_formula(a, g):.4f}")
for g in [0.0, 0.5, 1.0]:
    print(f"threshold (1-g)/(3-2g) at gamma={g}: {(1 - g) / (3 - 2 * g):.4f}")

# ----------------------------------------------------------------------------
hdr("C5  Candidate A: fuse stake with time (coin age), weight (1+s)^(1+t) - 1")
# The obvious fused-rung idea: a validator's weight is the fused exponent of
# its stake share s and the time t since it last produced a block.  Honest
# validators reset t when they produce; an attacker can stay silent and let t
# grow.  This is Peercoin's coin age (2012) with a different curve.


def coin_age_sim(att_share, weight, epochs=230, n_honest=40, n_att=20, window=200):
    """Each epoch one leader is picked with probability proportional to
    weight(stake, age).  Honest validators always produce when picked, and
    producing resets that account's age.  The attacker splits its stake over
    n_att accounts (as Peercoin attackers did), stays silent for `window`
    epochs while age builds up, then produces every time it is picked.
    Returns the attacker's weight share at switch-on and its share of the
    next 20 blocks."""
    n = n_honest + n_att
    stakes = np.array([(1 - att_share) / n_honest] * n_honest + [att_share / n_att] * n_att)
    is_att = np.arange(n) >= n_honest
    age = np.zeros(n)
    att_blocks = 0
    share_at_switch = None
    for e in range(epochs):
        if e < window:
            # attacker silent: honest blocks only, attacker ages
            w = weight(stakes, age) * ~is_att
        else:
            w = weight(stakes, age)
        probs = w / w.sum()
        if e == window:
            share_at_switch = probs[is_att].sum()
        leader = rng.choice(n, p=probs)
        age += 1
        age[leader] = 0
        if is_att[leader] and window <= e < window + 20:
            att_blocks += 1
    return share_at_switch, att_blocks / 20


lin = lambda s, t: s
fused = lambda s, t: (1 + s) ** (1 + t) - 1
coinage = lambda s, t: s * (1 + t)
print("attacker stake   weight rule            weight share at switch-on   share of next 20 blocks (avg of 200 runs)")
for att in [0.10, 0.20, 0.30]:
    for name, wf in [("linear stake", lin), ("coin age s*t", coinage), ("fused (1+s)^(1+t)-1", fused)]:
        shares, blocks = [], []
        for _ in range(200):
            sh, bl = coin_age_sim(att, wf)
            shares.append(sh)
            blocks.append(bl)
        print(f"{att:.2f}             {name:22s} {np.mean(shares):.3f}                       {np.mean(blocks):.3f}")

# ----------------------------------------------------------------------------
hdr("C6  Candidate B: weigh stake by log(1+s), additive under the circle operation")
# w(a (+) b) = w(a) + w(b) for w = log(1+s).  But real stake adds with +,
# so splitting s = s1 + s2 into two accounts changes the weight:
#   (1+s1)(1+s2) = 1 + s1 + s2 + s1*s2  ->  the gain is exactly the corner s1*s2.
for s, m in [(0.10, 2), (0.10, 10), (0.30, 10), (0.30, 1000)]:
    w_one = math.log1p(s)
    w_split = m * math.log1p(s / m)
    print(f"stake {s:.2f} split into {m:4d}: log weight {w_one:.4f} -> {w_split:.4f}  (gain {100 * (w_split / w_one - 1):.1f}%)")
s1, s2 = 0.12, 0.05
print("corner check: (1+s1)(1+s2) - (1+s1+s2) =", round((1 + s1) * (1 + s2) - (1 + s1 + s2), 6), "= s1*s2 =", s1 * s2)
print("other rules, stake 0.30 split 10 ways:")
for name, w in [("linear s", lambda s: s), ("sqrt s", math.sqrt), ("s^2", lambda s: s * s), ("log(1+s)", math.log1p)]:
    print(f"  {name:10s} weight ratio split/whole = {10 * w(0.03) / w(0.30):.3f}")

# ----------------------------------------------------------------------------
hdr("C7  Nothing at stake and slashing: detection is a k = -1 circle power")
# m independent watchers each catch a double vote with probability d.
# Caught = d (+)_{-1} d (+)_{-1} ... = 1 - (1-d)^m.
for d in [0.05, 0.2]:
    for m in [1, 10, 100]:
        print(f"d={d:.2f} m={m:4d}: caught {mpow(d, m, -1):.4f}  (check 1-(1-d)^m = {1 - (1 - d) ** m:.4f})")


def fork_sim(slash, n=100, rounds=60, reward=1.0, penalty=32.0, d=0.2, m=50):
    """Two forks of equal weight. Each validator votes A, B or both.
    Voting both pays the reward on whichever fork wins; it costs
    penalty * P(caught) if slashing exists. A fork wins once it holds
    2/3 of votes cast that round. Returns rounds until resolution."""
    caught = 1 - (1 - d) ** m
    for t in range(rounds):
        prefer_a = rng.random(n) < 0.5
        both_value = reward - (penalty * caught if slash else 0.0)
        single_value = 0.5 * reward  # expected, fork is a coin flip
        votes_both = both_value > single_value
        if votes_both:
            a = b = n
        else:
            a = int(prefer_a.sum())
            # once a side leads, rational validators join it next round
            if t > 0:
                a = n if a >= n / 2 else 0
            b = n - a
        if max(a, b) >= 2 / 3 * (a + b) and a != b:
            return t + 1
    return None


print("rounds to settle a fork, no slashing:", fork_sim(False), "(None = never: everyone votes both)")
print("rounds to settle a fork, with slashing:", fork_sim(True))

# ----------------------------------------------------------------------------
hdr("C8  Long-range attacks: exited stake accumulates, weak subjectivity")
# If a fraction c of the original validator set leaves each epoch, the part of
# the old set that has left after t epochs is c (+)_{-1} c ... = 1 - (1-c)^t.
# Old keys that have left cannot be slashed, so an attacker can buy them.
# Once they hold 1/3 of an old set, a fork from that point can't be told
# apart by safety alone: new nodes need a recent checkpoint.
staked = 36.0e6  # ETH staked, January 2026 (Dune via news reports)
cap = 256.0  # ETH per epoch, Electra MAX_PER_EPOCH_ACTIVATION_EXIT_CHURN_LIMIT
epoch_min = 6.4
c = cap / staked
t_lin = (staked / 3) / cap
t_mult = math.log(2 / 3) / math.log(1 - c)
print(f"exit churn cap {cap:.0f} ETH/epoch on {staked / 1e6:.0f}M ETH: c = {c:.3e} per epoch")
print(f"epochs until 1/3 of an old set could have left:")
print(f"  fixed cap (linear):            {t_lin:,.0f} epochs = {t_lin * epoch_min / 1440:.0f} days")
print(f"  proportional churn (k=-1 power): {t_mult:,.0f} epochs = {t_mult * epoch_min / 1440:.0f} days")
print("  (Ethereum's own weak-subjectivity rule is stricter; this is the upper bound of the idea, not their number)")

# ----------------------------------------------------------------------------
hdr("C9  Verifiable delay functions: repeated circle doubling")
# Wesolowski / Pietrzak VDFs compute y = x^(2^T): T squarings in a row.
# In circle-operation coordinates squaring is a (+) a, so T doublings climb a cube of
# dimension 2^T.  Pietrzak's proof splits T in half: x^(2^T) = (x^(2^(T/2)))^(2^(T/2)).
a = 0.001
x = a
for _ in range(10):
    x = madd(x, x)
print("10 doublings of a=0.001:", x, " vs (1+a)^(2^10)-1 =", (1 + a) ** 1024 - 1)
N = 1000003 * 999983  # toy modulus, not secure
x0, T = 5, 20
y = x0
for _ in range(T):
    y = y * y % N
half = x0
for _ in range(T // 2):
    half = half * half % N
y2 = half
for _ in range(T // 2):
    y2 = y2 * y2 % N
print(f"toy VDF: x^(2^{T}) mod N = {y}, via two halves = {y2}, equal: {y == y2}")

# ----------------------------------------------------------------------------
hdr("C10  Energy: proof of work vs proof of stake (published figures)")
btc_twh = 138.0  # CCAF Digital Mining Industry Report, 28 Apr 2025
eth_pow_twh = 23.0  # CCRI for ConsenSys, Sept 2022: ~23 million MWh/yr before the Merge
eth_pos_twh = 0.0026  # CCRI: ~2,600 MWh/yr after the Merge
print(f"Bitcoin (PoW):           {btc_twh:8.1f} TWh/yr")
print(f"Ethereum before Merge:   {eth_pow_twh:8.1f} TWh/yr")
print(f"Ethereum after Merge:    {eth_pos_twh:8.4f} TWh/yr")
print(f"Ethereum reduction:      {100 * (1 - eth_pos_twh / eth_pow_twh):.3f}%")
print(f"Bitcoin / post-Merge Ethereum: {btc_twh / eth_pos_twh:,.0f}x")
# Security as a flow vs a stock: PoW must keep paying for hash rate every block;
# PoS locks a stake once and burns slashed stake only on misbehaviour.
