"""Does the law's union rule plus 'information adds' force the 2 sqrt 2 circle? (DRAFT, Nov 28 2026)
Run:  python3 checks.py > checks_output.txt   (writes results.json)

  F1  in Pawlowski's tree, a measure that grows like |E|^a near a fair coin allows exactly the l^a ball
      |E1|^a + |E2|^a <= 1; its CHSH ceiling is 2^(2 - 1/a): a = 1 gives 2, a = 2 gives 2 sqrt 2, a -> inf gives 4
  F2  the law's bit flip is composing with the certain bit (p o_-2 1 = 1 - p), which turns the lift E into -E;
      a smooth measure symmetric under it starts at E^2 (a = 2); Shannon, Renyi-2 and chi^2 do; the guessing
      measure (Renyi-inf) and total variation are not smooth at E = 0 and start at |E| (a = 1)
  F3  the union rule (independent parts) makes every Renyi entropy add; only Shannon also obeys the chain rule
      H(X, Y) = H(X) + H(Y | X), which is what Pawlowski's proof uses
"""
import json
import numpy as np
import sympy as sp

rng = np.random.default_rng(112826)
results, passed, total = {}, 0, 0


def check(name, ok, detail, **extra):
    global passed, total
    total += 1; passed += bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    results[name] = {"pass": bool(ok), "detail": detail, **extra}


# F1: sum over the 2^n branches of |product of biases|^a = (|E1|^a + |E2|^a)^n: bounded iff the l^a norm is <= 1
a, E1, E2, n = sp.symbols("a E1 E2 n", positive=True)
branch_sum = sp.summation(sp.binomial(n, sp.Symbol("j")) * E1 ** (a * (n - sp.Symbol("j"))) * E2 ** (a * sp.Symbol("j")),
                          (sp.Symbol("j"), 0, n))
ok_sum = all(abs(float(branch_sum.subs({a: av, E1: 0.6, E2: 0.5, n: nv})) - (0.6 ** av + 0.5 ** av) ** nv) < 1e-9
             for av in (1, 2, 3.5) for nv in (1, 4, 9))
ceil = {}
for av in [1, 1.5, 2, 3, 4, 8]:
    r = 2 ** (-1 / av)                                   # symmetric point on the l^a ball, where E1 + E2 is largest
    ceil[av] = 4 * r
ok_ceil = abs(ceil[1] - 2) < 1e-12 and abs(ceil[2] - 2 * np.sqrt(2)) < 1e-12 and ceil[8] < 4
# confirm the symmetric point is the maximum (Lagrange, a = 3 as a test)
x, y, l = sp.symbols("x y l", positive=True)
sol = sp.solve([1 - l * 3 * x ** 2, 1 - l * 3 * y ** 2, x ** 3 + y ** 3 - 1], [x, y, l], dict=True)
ok_lag = abs(float(2 * (sol[0][x] + sol[0][y])) - ceil[3]) < 1e-12
check("F1", ok_sum and ok_ceil and ok_lag,
      "Pawlowski's tree multiplies biases along each branch (the law at k = -2), so a measure that grows like |E|^a "
      "near a fair coin sums over the 2^n branches to (|E1|^a + |E2|^a)^n exactly. One sent bit can only cover that "
      "for every n if |E1|^a + |E2|^a <= 1. The CHSH ceiling is then 2^(2 - 1/a): "
      + ", ".join(f"a = {k}: {v:.4f}" for k, v in ceil.items())
      + ". a = 1 gives the local 2, a = 2 gives 2 sqrt 2, large a heads to the PR box's 4. The ceiling is set by "
        "one number: how fast information grows with the bias", ceilings={str(k): v for k, v in ceil.items()})

# F2: the law's bit flip, and the small-bias exponent of standard measures
p = sp.symbols("p")
law = lambda u, v, kk: u + v + kk * u * v
flip = sp.simplify(law(p, 1, -2))
lift_flip = sp.simplify((1 - 2 * flip) + (1 - 2 * p))
measures = {
    "Shannon (1 - h)": lambda E: ((1 + E) * np.log2(1 + E) + (1 - E) * np.log2(1 - E)) / 2,
    "Renyi-2 (collision)": lambda E: np.log2(1 + E ** 2),
    "chi-squared": lambda E: E ** 2,
    "Renyi-inf (guessing)": lambda E: np.log2(1 + abs(E)),
    "total variation": lambda E: abs(E),
}
expo, sym = {}, {}
for name, f in measures.items():
    e1, e2 = 1e-3, 2e-3
    expo[name] = float(np.log(f(e2) / f(e1)) / np.log(e2 / e1))
    sym[name] = abs(f(0.37) - f(-0.37)) < 1e-12
want = {"Shannon (1 - h)": 2, "Renyi-2 (collision)": 2, "chi-squared": 2, "Renyi-inf (guessing)": 1, "total variation": 1}
ok2 = flip == 1 - p and lift_flip == 0 and all(abs(expo[k] - want[k]) < 1e-3 for k in want) and all(sym.values())
check("F2", ok2,
      f"flipping a bit is composing it with the certain bit, p o_(-2) 1 = {flip}, which sends the lift E = 1 - 2p to "
      f"-E. Every measure tested respects that flip; the smooth ones then must start at E^2: "
      + "; ".join(f"{k}: a = {v:.3f}" for k, v in expo.items())
      + ". Smooth + flip-symmetric gives a = 2 and the circle; the guessing measure and total variation have a "
        "corner at the fair coin and give a = 1, the local bound 2",
      exponents=expo)


# F3: union rule (independence) vs chain rule
def renyi(P, alpha):
    P = P[P > 0]
    if alpha == 1: return float(-(P * np.log2(P)).sum())
    if alpha == np.inf: return float(-np.log2(P.max()))
    return float(np.log2((P ** alpha).sum()) / (1 - alpha))


add_err, chain_err = {}, {}
for alpha in (1, 2, np.inf):
    worst_add = worst_chain = 0
    for _ in range(200):
        px = rng.dirichlet(np.ones(3)); py = rng.dirichlet(np.ones(4))
        joint_ind = np.outer(px, py)
        worst_add = max(worst_add, abs(renyi(joint_ind.ravel(), alpha) - renyi(px, alpha) - renyi(py, alpha)))
        J = rng.dirichlet(np.ones(12)).reshape(3, 4)          # correlated pair
        mx = J.sum(1)
        cond = sum(mx[i] * renyi(J[i] / mx[i], alpha) for i in range(3))
        worst_chain = max(worst_chain, abs(renyi(J.ravel(), alpha) - renyi(mx, alpha) - cond))
    add_err[str(alpha)], chain_err[str(alpha)] = worst_add, worst_chain
ok3 = all(v < 1e-12 for v in add_err.values()) and chain_err["1"] < 1e-12 and chain_err["2"] > 1e-3 and chain_err["inf"] > 1e-3
check("F3", ok3,
      f"for independent parts (the union rule's 'share nothing') every Renyi entropy adds: worst error "
      f"{max(add_err.values()):.0e} over 200 random pairs for orders 1, 2 and infinity. For correlated parts only "
      f"Shannon obeys the chain rule H(X, Y) = H(X) + H(Y | X): errors {chain_err['1']:.0e} (Shannon), "
      f"{chain_err['2']:.2f} (Renyi-2), {chain_err['inf']:.2f} (Renyi-inf), with the averaged conditional. So the union "
      f"rule cannot pick the measure; the chain rule picks Shannon (Faddeev 1956), and Shannon is smooth, so a = 2",
      additivity=add_err, chain=chain_err)

print(f"\n{passed}/{total} pass")
results["summary"] = f"{passed}/{total} pass"
json.dump(results, open("results.json", "w"), indent=1, default=float)
