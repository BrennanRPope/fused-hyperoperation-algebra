"""DRAFT: pictures for the polish section (numpy + svgkit + PIL; no matplotlib). Run: python3 figures_polish.py <dir with tbpy+libtb.so>"""
import sys, os, re, json, math
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "..")
HERE = os.path.dirname(os.path.abspath(__file__)); TB = os.path.dirname(HERE)
sys.path.insert(1, "/mnt/project-files/math-theory/newton-g")
import numpy as np, tbpy
from svgkit import Fig, Axes, BLUE, ORANGE, GRAY, RED, DARK, GREEN, PURPLE
FIG = os.path.join(TB, "figures")
P = {r["id"]: r for r in json.load(open(os.path.join(HERE, "polished.json")))}
A = {o["id"]: o for o in json.load(open(os.path.join(HERE, "words_all.json")))}
txt = open(os.path.join(HERE, "polish_output.txt")).read()
hist = {}
for blk in txt.split("ORBIT ")[1:]:
    i = int(blk.split()[0])
    s0 = max(abs(float(x)) for x in re.search(r"start_halfres (\S+) (\S+) (\S+)", blk).groups())
    dq = [float(a) for a, b in re.findall(r"\[(\S+) (\S+)\]", re.search(r"double rc .*", blk).group(0))]
    q = [float(a) for a, b in re.findall(r"\[(\S+) (\S+)\]", re.search(r"quad0 rc .*", blk).group(0))]
    hist[i] = [s0, min(dq)] + q[1:]
listed = lambda i: bool(A[i]["new_match"])

# ---------------- Figure 12: polishing every orbit
F = Fig(1200, 580); F.title("Figure 12. Polishing every orbit: Newton steps on half a lap, in 34-digit arithmetic")
ax = Axes(F, 90, 80, 470, 380, (-0.3, 4.3), (-34, -6))
ax.frame(xticks=[(0, "search"), (1, "double"), (2, "quad 1"), (3, "quad 2"), (4, "quad 3")],
         yticks=[(-k, f"1e-{k}") for k in (8, 12, 16, 20, 24, 28, 32)], grid=True, ylabel="how far the half-lap mirror test misses")
for i, h in hist.items():
    ys = [math.log10(max(v, 1e-34)) for v in h][:5]
    ax.plot(range(len(ys)), ys, color=ORANGE if listed(i) else BLUE, width=1.3, opacity=0.55)
ax.hline(-16, color=GRAY, dash="5,4", width=1); ax.text(3.0, -15.2, "ordinary double precision", size=12, color="#555")
F.lines(90, 515, "A. Each line is one orbit. Each quad step roughly doubles the correct digits,\nuntil the 34-digit arithmetic runs out.", size=13)
bx = Axes(F, 690, 80, 470, 380, (0, 44), (0, 34))
bx.frame(xticks=[0, 10, 20, 30, 40], yticks=[0, 5, 10, 15, 20, 25, 30], grid=True, xlabel="period T", ylabel="correct digits")
for i, r in P.items():
    T = float(r["T"]); dd = min(16, -math.log10(max(r["double_close"], 1e-17))); qd = r["digits"]
    bx.f.line(bx.X(T), bx.Y(dd), bx.X(T), bx.Y(qd), color="#bbb", width=1)
    bx.dot(T, dd, r=3, color=GRAY); bx.dot(T, qd, r=4.5, color=ORANGE if listed(i) else BLUE)
bx.dot(2, 4, r=3, color=GRAY); bx.text(3, 3.3, "double precision", size=12)
bx.dot(2, 1.5, r=4.5, color=ORANGE); bx.text(3, 0.8, "polished, in the 2017 list", size=12)
bx.dot(24, 1.5, r=4.5, color=BLUE); bx.text(25, 0.8, "polished, not in it", size=12)
F.lines(690, 515, "B. Grey: best double precision can do (8 to 14 digits). Coloured: polished,\nall 66 orbits close to 26 to 32 digits.", size=13)
F.save(os.path.join(FIG, "fig12-polishing-every-orbit.svg"))

# ---------------- Figure 13: loop types re-read
C = {o["id"]: o for o in json.load(open(os.path.join(HERE, "words_converge.json")))}
from words_converge import word_at
changed = [i for i in A if A[i]["word_changed"]]
lens = {}
for i in changed:
    if i in C: lens[i] = [len(C[i]["words"][k]) for k in ("3000", "12000", "48000", "192000")]
    else:
        p1, p2, T = float(P[i]["p1"]), float(P[i]["p2"]), float(P[i]["T"])
        lens[i] = [len(word_at(p1, p2, T, k)) for k in (3000, 12000, 48000, 192000)]
    lens[i] = [len(P[i]["word"])] + lens[i][1:]
F = Fig(1200, 620); F.title("Figure 13. Re-reading the loop types: 11 words were cut short by coarse sampling")
ax = Axes(F, 90, 80, 470, 400, (-0.2, 3.2), (0, 50))
ax.frame(xticks=[(0, "3,000 (old)"), (1, "12,000"), (2, "48,000"), (3, "192,000")], yticks=[0, 10, 20, 30, 40, 50],
         grid=True, xlabel="samples per unit time", ylabel="letters in the loop word")
for i in changed:
    col = BLUE if not listed(i) else (GREEN if not A[i]["old_match"] else ORANGE)
    ax.plot(range(4), lens[i], color=col, width=2, opacity=0.8)
F.lines(90, 530, "A. Each line is one orbit. The word grows, then stops changing by 48,000\nsamples (the three hardest checked to 768,000). Green: now in the 2017 list.\nOrange: still listed, but a different family. Blue: still not in the list.", size=13)
# B: why: the tight pair spins between samples
i = 35; p1, p2, T = float(P[i]["p1"]), float(P[i]["p2"]), float(P[i]["T"])
t, X = tbpy.traj(tbpy.suvakov(p1, p2), T, int(T * 3000), 1e-13)
r = X[:, :6].reshape(-1, 3, 2); pairs = [(0, 1), (0, 2), (1, 2)]
d = np.array([np.linalg.norm(r[:, a] - r[:, b], axis=1) for a, b in pairs]); k = np.unravel_index(np.argmin(d), d.shape)
pa, kk = pairs[k[0]], k[1]; tc = t[kk]
w = 0.004
t2, X2 = tbpy.traj(X[max(kk - 40, 0)], (min(kk + 40, len(t) - 1) - max(kk - 40, 0)) * (t[1] - t[0]), 20000, 1e-13)
t2 = t2 + t[max(kk - 40, 0)]
rel = lambda Y: Y[:, 2*pa[1]:2*pa[1]+2] - Y[:, 2*pa[0]:2*pa[0]+2]
R2 = rel(X2); R1 = rel(X[max(kk - 40, 0):kk + 41])
rmin_ = np.linalg.norm(R2, axis=1).min(); s = 5 * rmin_
bx = Axes(F, 720, 90, 400, 400, (-s, s), (-s, s))
bx.frame(xticks=[], yticks=[])
m = np.abs(R2).max(axis=1) < s; bx.plot(R2[m, 0], R2[m, 1], color=BLUE, width=2)
m1 = np.abs(R1).max(axis=1) < s
for x, y in R1[m1]: bx.dot(x, y, r=5, color=RED)
bx.dot(0, 0, r=3, color=DARK); bx.text(-0.15*s, -0.12*s, f"body {pa[0]+1}", size=12)
F.lines(720, 530, f"B. Orbit T = 27.8: body {pa[1]+1} as seen from body {pa[0]+1}, box {2*s:.3f} wide. Blue: the real\npath. Red: the only look the original 3,000-per-unit sampling took inside the\nbox. The pair whips half a turn round each other between looks, so the\nline-ups with the third body during the turn were never counted.", size=13)
F.save(os.path.join(FIG, "fig13-loop-types-re-read.svg"))

# ---------------- Figure 10b: the 6 orbits still not in the 2017 list
un = sorted([i for i in A if not listed(i)], key=lambda i: float(P[i]["T"]))
F = Fig(1200, 820); F.title("Figure 10b. The 6 polished orbits whose loop types are not in the 2017 list")
cols = [ORANGE, BLUE, GREEN]
for n, i in enumerate(un):
    p1, p2, T = float(P[i]["p1"]), float(P[i]["p2"]), float(P[i]["T"])
    t, X = tbpy.traj(tbpy.suvakov(p1, p2), T, int(T * 4000), 1e-13)
    R = X[:, :6].reshape(-1, 3, 2); L = 1.1 * np.abs(R).max()
    cx, cy = 40 + (n % 3) * 390, 60 + (n // 3) * 380
    ax = Axes(F, cx, cy, 340, 300, (-L * 340 / 300, L * 340 / 300), (-L, L))
    ax.frame(xticks=[], yticks=[])
    for b in range(3): ax.plot(R[::3, b, 0], R[::3, b, 1], color=cols[b], width=1.1, opacity=0.85)
    F.text(cx + 170, cy + 325, f"T = {T:.2f}, {len(A[i]['word'])}-letter word, {P[i]['digits']} digits", size=13, anchor="middle")
F.save(os.path.join(FIG, "fig10b-unlisted-polished.svg"))
print("ok")
