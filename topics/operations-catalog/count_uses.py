"""Count how often each catalogued move appears across math-theory (md, py, lean, txt; not build output, html or this folder).
Usage: python3 -I count_uses.py /mnt/project-files/math-theory"""
import os, re, sys, json
root = sys.argv[1]
C = {
 "two-dial law (lens)": r"\blens|Coleman|landscape law|two-dial law|∘ₖ,d",
 "lift 1+kx": r"\blift|1 \+ k·?x\b|1\+kx|1 \+ k \* x|1 \+ k\*x|x ↦ 1 \+ x",
 "way back (quasi-inverse)": r"way back|quasi-?inverse|IsQuasiregular",
 "wall -1/k": r"\bwall",
 "circle power": r"\bmpow|circle power|mult-add power|\^\(✳n\)|\^\(✳",
 "n-fold combine / cube corners": r"foldW|n-fold|elementary symmetric|corner weight|\bD_k\b|\bDk\b|det\(I \+ kA\)",
 "Möbius inversion": r"Möbius inversion|Mobius inversion|\bRota\b",
 "flip x/(1+kx) (odds/probability)": r"\bodds\b|\bprob k\b|weight_prob_flip|additive duality|probability side",
 "mirror x -> -x": r"\bmirror",
 "logarithm L_k (bend / ruler / flat coordinate)": r"\bbend(s|ing)?\b|flat coordinate|\bruler|log\(1 \+ kx\)|log\(1\+kx\)|log1p",
 "k-log / k-exp (Box-Cox, Tsallis)": r"k-log|q-log|k-exp|q-exp|Box.Cox|q-logarithm|q-exponential",
 "slide y/(1+cy) (chart change)": r"shift_conj|chart change|y/\(1 \+ cy\)|x/\(1 − cx\)|coordinate change",
 "rescale x -> λx (units, Maslov)": r"rescal|units_change_weight|units_num|wAdd_scale|dequantiz|Maslov|zoom out",
 "classes by discriminant (countries)": r"\bcountr(y|ies)\b|k² − 4q|k\^2 - 4q|p² − 4q|discriminant",
 "turning a dial: deformation / contraction": r"contraction|Inönü|deformation",
 "height flow E^t (regular iteration)": r"E\^s|E\^t|E\^\{|E\^\(|regular iteration|fractional iteration|half.rung|rung flow|Abel function|Fatou",
 "height-s operation B_s": r"\bB_s\b|B_\{s|\bB_-?\d|B_\(",
 "climb a rung / conjugation": r"\bclimb|conjugat|transport",
 "fuse (shift by 1)": r"\bfused? rung|\bsadd\b|\bfexp\b|\bftet\b|\bfuse\b",
 "translation T_a (motion)": r"\bT_a\b|T_\(a|motion rule|translation by a|app_comp|T_comp|stretch-and-shift",
 "exponential tilting": r"\btilt|Esscher",
 "zoom (Moran exponent)": r"\bzoom|Moran",
 "order gap (wedge)": r"order_gap|order gap|\bwedge",
}
C = {k: re.compile(v) for k, v in C.items()}
counts = {k: [0, 0] for k in C}
nfiles = 0
for dp, dn, fn in os.walk(root):
    dn[:] = [d for d in dn if d not in ("build", ".lake", "operations-catalog", "__pycache__", ".git")]
    for f in fn:
        if not f.endswith((".md", ".py", ".lean", ".txt")): continue
        p = os.path.join(dp, f)
        try:
            if os.path.getsize(p) > 2_000_000: continue
            t = open(p, encoding="utf-8", errors="ignore").read()
        except OSError: continue
        nfiles += 1
        for k, rx in C.items():
            n = len(rx.findall(t))
            if n: counts[k][0] += 1; counts[k][1] += n
print(f"files scanned: {nfiles}\n")
print(f"{'move':48s} {'files':>6s} {'uses':>6s}")
for k, (f, u) in sorted(counts.items(), key=lambda kv: -kv[1][1]):
    print(f"{k:48s} {f:6d} {u:6d}")
json.dump({k: {"files": v[0], "uses": v[1]} for k, v in counts.items()}, open("counts.json", "w"), indent=1, ensure_ascii=False)
