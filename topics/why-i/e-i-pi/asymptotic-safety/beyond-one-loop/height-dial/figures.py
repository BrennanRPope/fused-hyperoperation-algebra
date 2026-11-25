"""Figures for the height-dial test. Run: python3 figures.py (numpy + PIL via ../../../svgkit.py)."""
from hd import *
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../../..")))
from svgkit import *
_lines = Fig.lines
Fig.lines = lambda self, x, y, s, **kw: _lines(self, x, y, "\n".join(s) if isinstance(s, list) else s, **kw)
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True); P = lambda n: os.path.join(OUT, n)
def save(f, n): f.save(P(n + ".svg")); f.png(P(n + ".png"))

v = newton(eh_beta(), [0.7, 0.19]); th = np.linalg.eigvals(-jac(eh_beta(), v))
L_lit = 23/22; s_star = s_for(L_lit)
nums = [("Litim, λ = 0: θ = 23/11", L_lit, GREEN),
        ("Bonanno–Reuter's own β", 1.18, ORANGE),
        ("two loops, 2 + ε at ε = 2", ak_ratio(2.0), RED),
        ("Litim, full EH (real part)", th[0].real/2, PURPLE),
        ("large N (exact)", 1.0, BLUE)]

# ---------- h1: one height per number ----------
f = Fig(1180, 700); f.title("Each published exponent needs its own height: the dial is re-set every time", size=17)
ax = Axes(f, 90, 80, 520, 420, (0, 2), (1.9, 4.5))
ax.frame(xticks=[0, 0.5, 1, 1.5, 2], yticks=[2, 2.5, 3, 3.5, 4, 4.5], xlabel="rung height s", ylabel="e at height s,  e_s = 1 + E^s(1)", grid=True)
S = np.linspace(0, 1.78, 72); ax.plot(S, [e_s(s) for s in S], color=DARK, width=3)
for name, L, col in nums:
    s = s_for(L); y = e_s(s)
    ax.plot([0, s], [y, y], color=col, width=1.5, dash="5,4"); ax.plot([s, s], [1.9, y], color=col, width=1.5, dash="5,4")
    ax.dot(s, y, r=6, color=col)
ax.text(1.02, 2.62, "e_1 = e", size=12, color=BLUE)
ax.text(0.04, 3.0, f"s* = {s_star:.3f}, e_s* = e^(23/22) = {e_s(s_star):.3f}", size=12, color=GREEN)
for i, (name, L, col) in enumerate(nums):
    f.line(100, 558 + i*22, 128, 558 + i*22, color=col, width=3); f.text(136, 563 + i*22, f"{name}:  s = {s_for(L):.3f}", size=12)
ax2 = Axes(f, 700, 80, 430, 420, (2, 6), (0.95, 1.75))
ax2.frame(xticks=[2, 3, 4, 5, 6], yticks=[1.0, 1.2, 1.4, 1.6], xlabel="spacetime dimension d", ylabel="height needed to match θ", grid=True)
D = np.linspace(2.02, 6, 60); ax2.plot(D, [s_for(litim0_ratio(d)) for d in D], color=GREEN, width=2.5)
D2 = np.linspace(2.02, 4, 40); ax2.plot(D2, [s_for(ak_ratio(d - 2)) for d in D2], color=RED, width=2.5)
ax2.hline(s_star, color=GREEN, dash="3,3", width=1); ax2.dot(4, s_star, r=5, color=GREEN)
ax2.text(4.6, 1.15, "Litim, λ = 0", size=12, color=GREEN); ax2.text(3.35, 1.45, "two loops", size=12, color=RED)
ax2.text(2.1, 0.965, "every curve starts at s = 1: one loop is exact near d = 2", size=12, color=GRAY)
f.lines(700, 560, ["Left: black curve is ours (e-every-height). Each dot is the",
                   "height that turns one published exponent into θ = (d − 2) ln e_s.",
                   "Five numbers, five heights.",
                   "Right: the height that works moves with d, so it is not a",
                   "fixed property of e.",
                   "Known: the five exponents (README). Ours: e_s and the heights.",
                   "Checks H2–H6."], size=12, color=GRAY)
save(f, "h1-one-height-per-number")

# ---------- h2: the low-energy end ----------
f = Fig(1180, 600); f.title("Only height 1 lets G stop running at low energy", size=17)
ax = Axes(f, 100, 80, 640, 400, (0, 5), (0, 6))
ax.frame(xticks=range(0, 6), yticks=range(0, 7), xlabel="log10 a   (a = 1/g; large a = low energy, G nearly constant)", ylabel="rate at that point,  V_s(a) / a", grid=True)
A = np.logspace(0, 5, 40); cols = {0.5: GRAY, 1: BLUE, 1.1: GREEN, 1.5: ORANGE, 2: RED}
for s, col in cols.items():
    R = np.array([V(s, a)/a for a in A]); m = R <= 6
    ax.plot(np.log10(A)[m], R[m], color=col, width=3 if s == 1 else 2)
    ax.text(5.05, min(5.9, V(s, 1e5)/1e5) - 0.05, f"s = {s}", size=12, color=col)
ax.text(2.2, 0.62, "s = 1:  V(a) = 1 + a, rate → 1 exactly", size=12, color=BLUE)
f.lines(800, 110, ["What this shows",
                   "",
                   "At low energy, G is a constant and",
                   "g = G k² grows exactly like k².",
                   "That needs a fixed rate at a = ∞.",
                   "",
                   "s = 1 settles to 1.",
                   "s = 1.1 (the fitted height) keeps",
                   "climbing, slowly, like a log.",
                   "s < 1 sinks to 0, s > 1 climbs.",
                   "",
                   "At the other end, the wall E^s(−∞),",
                   "the rate is 1 for every s (check H8).",
                   "So the height can't make the two",
                   "ends' rates differ by 23/22.",
                   "It only spoils the low-energy end."], size=13)
f.lines(100, 540, ["Curves are ours: V_s(a) = (E^s)′(E^−s(a)) is the rate of the height-s law B_s(a, b) = E^s(E^−s a + E^−s b), computed by regular iteration.",
                   "Rescaling to (1/k) B_s(ka, kb) multiplies both end rates by k, so the ratio is unchanged (check H10). Checks H7–H10."], size=12, color=GRAY)
save(f, "h2-low-energy-end")

# ---------- h3: scorecard ----------
f = Fig(1180, 560); f.title("Fix s* with 23/11, then test it: no second number comes out right", size=17)
ax = Axes(f, 120, 80, 640, 360, (-0.5, 4.5), (0, 3.4))
ax.frame(yticks=[0, 1, 2, 3], ylabel="critical exponent θ", grid=True)
pred = 2*L_lit
for i, (name, L, col) in enumerate(nums):
    ax.plot([i - 0.3, i - 0.3, i - 0.02, i - 0.02], [0, 2*L, 2*L, 0], color=col, fill=col, opacity=0.85, close=True)
    ax.plot([i + 0.02, i + 0.02, i + 0.3, i + 0.3], [0, pred, pred, 0], color=GRAY, fill=GRAY, opacity=0.35, close=True)
    off = pred/(2*L) - 1
    ax.text(i, max(2*L, pred) + 0.1, "fitted" if i == 0 else f"{100*off:+.0f}%", size=13, anchor="middle", weight="bold", color=DARK if i else GREEN)
    f.lines(ax.X(i), 465, ["Litim\nλ = 0", "Bonanno–Reuter\nown β", "two loops\nε = 2", "Litim full EH\n(real part)", "large N\n(exact)"][i].split("\n"), size=11, anchor="middle")
ax.text(2.7, 1.62, "+ 3.04 i in the real one", size=11, color=PURPLE, anchor="middle")
f.lines(800, 110, ["colored: published θ",
                   "gray: what s* = 1.101 says",
                   "",
                   "Large N is exact at s = 1.",
                   "s* would make it 2.09 and 4.18",
                   "in place of 2 and 4.",
                   "",
                   "The full Einstein–Hilbert",
                   "exponent is complex. A real",
                   "height gives a real number, and",
                   "a complex height can't keep G",
                   "real (see README §3).",
                   "",
                   "One match from one free",
                   "parameter, zero after that."], size=13)
f.lines(120, 520, ["Known: published exponents. Ours: the s* prediction θ = 2 ln e_s* = 23/11 applied to each. Checks H4, H5, H12."], size=12, color=GRAY)
save(f, "h3-scorecard")
print("ok")
