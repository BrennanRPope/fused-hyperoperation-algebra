"""Does the law force Shannon's chain rule? (DRAFT, Nov 28 2026)
Run:  python3 checks.py > checks_output.txt   (writes results.json)

  C1  the union rule applied step by step is exact for dependent events: P(A or B) = P(A) o_-1 P(B | not A);
      its lift 1 - p multiplies. The law has a chain rule for probabilities.
  C2  for entropies, every Renyi order obeys a chain rule of the law's own shape:
      H(X,Y) = H(X) + [lift, average, unlift] of the H(Y|x), with the law's lift at dial k = (1 - alpha) ln 2
      and weights p_x^alpha. Shannon is the dial k = 0, the only one whose weights are the plain probabilities.
  C3  every finite dial gives information ~ E^2 near a fair coin, hence the same circle in Pawlowski's tree;
      only the dial at -infinity (guessing) is not smooth and gives the local bound
"""
import json
import numpy as np

rng = np.random.default_rng(112827)
results, passed, total = {}, 0, 0


def check(name, ok, detail, **extra):
    global passed, total
    total += 1; passed += bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    results[name] = {"pass": bool(ok), "detail": detail, **extra}


law = lambda a, b, k: a + b + k * a * b

# C1
worst = 0
for _ in range(500):
    J = rng.dirichlet(np.ones(4)).reshape(2, 2)      # J[a, b], dependent events A = {a=1}, B = {b=1}
    pA = J[1].sum()
    pB_notA = J[0, 1] / J[0].sum()
    worst = max(worst, abs((1 - J[0, 0]) - law(pA, pB_notA, -1)))
check("C1", worst < 1e-14,
      f"for 500 random dependent pairs, P(A or B) = P(A) o_-1 P(B | not A) to {worst:.0e}: the union rule, applied one "
      f"step at a time with the conditional chance, is exact. Its lift 1 - p multiplies (1 - P(A or B) = "
      f"(1 - P(A))(1 - P(B | not A))). This is the ordinary rule of conditional probability, in the law's form")


# C2
def H(P, al):
    P = P[P > 0]
    if al == 1: return float(-(P * np.log2(P)).sum())
    return float(np.log2((P ** al).sum()) / (1 - al))


def law_mean(vals, w, k):
    vals, w = np.asarray(vals), np.asarray(w)
    if k == 0: return float((w * vals).sum())
    Ek = (np.exp(k * vals) - 1) / k                  # exponential map
    return float(np.log(1 + k * (w * Ek).sum()) / k)  # lift of the weighted average


alphas = [0.5, 1.0, 2.0, 3.0]
err_law, err_plain, wgap = {}, {}, {}
for al in alphas:
    k = (1 - al) * np.log(2)
    e_law = e_plain = g = 0
    for _ in range(300):
        J = rng.dirichlet(np.ones(12)).reshape(3, 4)
        px = J.sum(1)
        cond = [H(J[i] / px[i], al) for i in range(3)]
        w_esc = px ** al / (px ** al).sum()
        e_law = max(e_law, abs(H(J.ravel(), al) - H(px, al) - law_mean(cond, w_esc, k)))
        e_plain = max(e_plain, abs(H(J.ravel(), al) - H(px, al) - float((px * np.array(cond)).sum())))
        g = max(g, np.abs(w_esc - px).max())
    err_law[str(al)], err_plain[str(al)], wgap[str(al)] = e_law, e_plain, g
ok2 = (all(v < 1e-12 for v in err_law.values()) and err_plain["1.0"] < 1e-12 and wgap["1.0"] < 1e-15
       and all(err_plain[str(a)] > 1e-3 and wgap[str(a)] > 1e-3 for a in alphas if a != 1.0))
check("C2", ok2,
      "every Renyi order alpha obeys a chain rule shaped like the law: H(X,Y) = H(X) + (lift, weighted average, "
      "unlift) of the conditional entropies, using the law's lift at dial k = (1 - alpha) ln 2 and weights "
      "p_x^alpha. Worst error " + ", ".join(f"alpha = {a}: {err_law[str(a)]:.0e}" for a in alphas)
      + ". Only k = 0 (Shannon) uses the plain probabilities as weights (gap "
      + ", ".join(f"{wgap[str(a)]:.2f}" for a in alphas)
      + "); with plain weights the other orders miss by " + ", ".join(f"{err_plain[str(a)]:.2f}" for a in alphas if a != 1.0)
      + ". So the law does not force Shannon: it holds the whole family, one chain rule per dial",
      law_chain_error=err_law, plain_chain_error=err_plain, weight_gap=wgap)


# C3
def info(E, al):     # 1 - H_alpha of a coin with bias E
    P = np.array([(1 + E) / 2, (1 - E) / 2])
    if al == np.inf: return float(1 + np.log2(P.max()))
    if al == 1:
        return float(((1 + E) * np.log1p(E) + (1 - E) * np.log1p(-E)) / (2 * np.log(2)))
    return float(1 - np.log2((P ** al).sum()) / (1 - al))


expo = {}
for al in [0.5, 1, 2, 3, 10, np.inf]:
    e1, e2 = 1e-4, 2e-4
    expo[str(al)] = float(np.log(info(e2, al) / info(e1, al)) / np.log(2))
ok3 = all(abs(v - 2) < 1e-3 for k_, v in expo.items() if k_ != "inf") and abs(expo["inf"] - 1) < 1e-3
check("C3", ok3,
      "near a fair coin, information grows like bias^a with a = " + ", ".join(f"{v:.3f} (alpha {k_})" for k_, v in expo.items())
      + ". Every finite dial is smooth and gives a = 2, so in Pawlowski's tree every one of them would draw the same "
        "circle and the same 2 sqrt 2 (forces-circle/, F1; the constant in front drops out). Only alpha = infinity "
        "(guessing, dial k -> -infinity) gives a = 1 and the local 2. What differs between dials is not the shape but "
        "the proof that one sent bit carries at most one bit: Pawlowski's proof uses Shannon's (k = 0) chain rule",
      exponents=expo)

print(f"\n{passed}/{total} pass")
results["summary"] = f"{passed}/{total} pass"
json.dump(results, open("results.json", "w"), indent=1, default=float)
