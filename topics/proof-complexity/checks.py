"""Checks for the proof-complexity write-up. Run: python3 checks.py
Reads yard.json (written by lean/ProofComplexity/Measure.lean) and asserts every number-based
claim made in README.md, printing the table it rests on."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "yard.json")))
by = {(r["result"], r["route"]): r for r in rows}
results = sorted({r["result"] for r in rows})

def cost(res, route, vs):
    """Route-specific cost: own declaration lines + lines of every theorem in this route's
    dependency closure that the other route's closure does not contain."""
    r = by[(res, route)]
    return r["own_lines"] + r[f"unshared_lines_vs_{vs}"]

n_checks = 0
def check(cond, msg):
    global n_checks
    n_checks += 1
    print(("PASS " if cond else "FAIL ") + msg)
    assert cond, msg

print(f"{'result':26} {'std':>6} {'op':>6} {'ctl':>6}   (route-specific lines, each vs standard; std vs op)")
for res in results:
    s = cost(res, "standard", "operation"); o = cost(res, "operation", "standard")
    c = cost(res, "control", "standard") if (res, "control") in by else None
    print(f"{res:26} {s:6} {o:6} {c if c is not None else '-':>6}")

# 1. Surface length: the standard proof is never longer than the operation proof.
for res in results:
    check(by[(res, "standard")]["own_lines"] <= by[(res, "operation")]["own_lines"],
          f"{res}: standard own lines {by[(res,'standard')]['own_lines']} <= operation {by[(res,'operation')]['own_lines']}")
    check(by[(res, "standard")]["term"] < by[(res, "operation")]["term"],
          f"{res}: standard proof term {by[(res,'standard')]['term']} < operation {by[(res,'operation')]['term']}")

# 2. Route-specific cost, operation vs standard: lower for R1-R3, higher for R4-R6.
for res in results:
    s, o = cost(res, "standard", "operation"), cost(res, "operation", "standard")
    if res[:2] in ("R1", "R2", "R3"):
        check(o < s, f"{res}: operation route cost {o} < standard {s}")
    else:
        check(o > s, f"{res}: operation route cost {o} > standard {s}")

# 3. The control (same plan, no operation) is never costlier than the operation route.
for res in results:
    if (res, "control") not in by: continue
    o_vs_c = cost(res, "operation", "control"); c_vs_o = cost(res, "control", "operation")
    check(c_vs_o < o_vs_c, f"{res}: control {c_vs_o} < operation {o_vs_c} (head to head)")
    check(cost(res, "control", "standard") <= cost(res, "operation", "standard"),
          f"{res}: vs standard, control {cost(res,'control','standard')} <= operation {cost(res,'operation','standard')}")

# 4. Wherever the operation route beats the standard one, the control beats it too.
for res in results:
    if (res, "control") in by and cost(res, "operation", "standard") < cost(res, "standard", "operation"):
        check(cost(res, "control", "standard") < cost(res, "standard", "control"),
              f"{res}: control also beats standard ({cost(res,'control','standard')} < {cost(res,'standard','control')})")

# 5. Depth of the dependency chain differs by at most 4 between routes.
for res in results:
    ds = [by[(res, r)]["depth"] for r in ("standard", "operation", "control") if (res, r) in by]
    check(max(ds) - min(ds) <= 4, f"{res}: depths {ds} within 4")

print(f"\n{n_checks} checks passed.")
