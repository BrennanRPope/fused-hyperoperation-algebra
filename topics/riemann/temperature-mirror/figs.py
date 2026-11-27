"""Temperature-mirror figures. DRAFT."""
import numpy as np, json, sys, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import mirror
from farend_pos import zeros
plt.rcParams.update({"font.size": 10})
def load(fr, n): return dict(json.loads(l) for l in open("mz_%s_%d_compass.jsonl" % (fr, n)))
def line_split(w, X, y):
    g = mirror.ev(w, X, 1j*y); return y[1:][np.sign(g.real[1:]) != np.sign(g.real[:-1])], y[1:][np.sign(g.imag[1:]) != np.sign(g.imag[:-1])]

def figA():
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2)); S = mirror.frame("spin", 10); th = np.radians(45)
    w, X = mirror.terms(S, th, "compass")
    zg = zeros(w, X, -2.5, 2.5, 0.3, 20, h=0.01)
    y = np.arange(0.3, 20, 0.0005); fz, oz = line_split(w, X, y)
    a = ax[0]; a.axvspan(-2.5, 0, color="tab:blue", alpha=0.07); a.axvspan(0, 2.5, color="tab:red", alpha=0.05)
    a.plot(zg.real, zg.imag, "o", color="tab:blue", ms=6, label="zeros of G (compass-up half): all on the left")
    a.plot(-zg.real, zg.imag, "o", mfc="none", color="tab:red", ms=6, label="zeros of G(−z), the mirror image: all on the right")
    a.plot(np.zeros_like(fz), fz, "D", color="k", ms=6, label="zeros of F = G(z) + G(−z): forced onto the line")
    a.plot(np.full_like(oz, 0.18), oz, "_", color="tab:green", ms=12, mew=2.5, label="zeros of G(z) − G(−z) (drawn just right): interlace with F's")
    a.axvline(0, color="k", lw=0.8); a.set_xlim(-2.5, 2.5); a.set_ylim(0, 26); a.set_xlabel("Re z  (0 = the line)"); a.set_ylabel("Im z")
    a.legend(fontsize=7.5, loc="upper center", ncol=1); a.set_title("A1. Mirror + one-sided half ⇒ every zero on the line\n(spin frame, n = 10, rewind 45°)", fontsize=10)
    b = ax[1]; b.axis("off")
    txt = ("The compass spin already builds the mirror in:\n"
           "    F(z) = G(z) + G(−z)  exactly, at every rewind angle.\n\n"
           "A mirror alone does NOT keep zeros on the line:\n"
           "    zeros may still leave, as mirror pairs (z, −z̄).\n\n"
           "Mirror + ONE-SIDEDNESS does (Hermite–Biehler):\n"
           "    if G has no zeros on the right half, then\n"
           "    |G(z)| > |G(−z)| on the right, |G(z)| < |G(−z)| on the left,\n"
           "    and on the line they are equal. So F = 0 only on the line,\n"
           "    every zero is simple, and they interlace with G(z) − G(−z).\n\n"
           "At 0° this is exactly the chain-fix field condition\n"
           "(P zero-free in the disk). The test: does it survive\n"
           "every angle of the rewind?")
    b.text(0.02, 0.98, txt, va="top", fontsize=10, family="DejaVu Sans")
    b.set_title("A2. What a temperature mirror would force", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figA_mirror_forces.png", dpi=150); plt.close(fig)

def figB():
    lc = json.load(open("../pendulum-swap/linecount.json")); fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    a = ax[0]
    for n, ls in ((8, ":"), (10, "--"), (12, "-")):
        for fr, c in (("spin", "tab:blue"), ("dual", "tab:orange")):
            d = load(fr, n); ks = sorted(d); a.plot(ks, [d[k] + (0.05 if fr == "dual" else 0) for k in ks], ls, color=c, marker="o", ms=2.5, label="%s frame, n = %d" % (fr, n))
    for n, col in ((10, "tab:orange"), (12, "tab:orange")):
        off = [r[0] for r in lc["dual frame n=%d" % n]]; a.axvspan(min(off), max(off), color=col, alpha=0.08 if n == 10 else 0.05)
    a.text(6, 3.6, "shaded: where the dual frame's zeros\nactually leave the line (n = 10, 12)", fontsize=8.5, color="tab:orange")
    a.set_xlabel("rewind angle θ (degrees), safe half only"); a.set_ylabel("zeros of G on the wrong side (Im ≤ 36)")
    a.set_ylim(-0.2, 4.4); a.legend(fontsize=8, loc="center right")
    a.set_title("B1. Spin frame: one-sided at every angle (0 wrong-side zeros).\nDual frame: fails from 1–3°, just before zeros leave", fontsize=10)
    b = ax[1]; S = mirror.frame("dual", 10); th = np.radians(10); w, X = mirror.terms(S, th, "compass")
    zg = zeros(w, X, -2.5, 2.5, 0.3, 20, h=0.01); W = np.concatenate([w, w]); XX = np.concatenate([X, -X])
    zf = zeros(W, XX, -1.0, 1.0, 0.3, 20, h=0.01)
    good = zg[zg.real < 0]; bad = zg[zg.real > 0]
    b.axvspan(0, 2.5, color="tab:red", alpha=0.05)
    b.plot(good.real, good.imag, "o", color="tab:blue", ms=6, label="zeros of G on the allowed side")
    b.plot(bad.real, bad.imag, "o", color="tab:red", ms=9, label="zeros of G on the WRONG side")
    on = zf[np.abs(zf.real) < 1e-6]; off = zf[np.abs(zf.real) >= 1e-6]
    b.plot(on.real, on.imag, "D", color="k", ms=5, label="zeros of F on the line")
    b.plot(off.real, off.imag, "*", color="tab:purple", ms=14, label="zeros of F OFF the line (a mirror pair)")
    b.axvline(0, color="k", lw=0.8); b.set_xlim(-2.5, 2.5); b.set_ylim(0, 20); b.set_xlabel("Re z"); b.legend(fontsize=8, loc="upper right")
    b.set_title("B2. Dual frame, n = 10, 10°: a G zero crosses to the wrong side;\nnear it F loses its line zeros as a mirror pair", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figB_one_sided_test.png", dpi=150); plt.close(fig)

def figC():
    pos = json.load(open("farend_pos.json")); fe = json.load(open("farend.json")); fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    a = ax[0]
    for n, c in (("10", "tab:green"), ("12", "tab:orange"), ("16", "tab:purple")):
        L = np.array(pos[n]["left"]); R = np.array(pos[n]["right"])
        a.plot(L[:, 0], L[:, 1], "o", mfc="none", color=c, ms=5); a.plot(R[:, 0], R[:, 1], "o", color=c, ms=7, label="n = %s (filled: right of the line)" % n)
    for t in (14.1347, 21.0220, 25.0109, 30.4249, 32.9351): a.axhline(t, color="0.75", ls=":", lw=0.8)
    a.axvline(0, color="k"); a.axvspan(0, 1, color="tab:red", alpha=0.05); a.set_xlim(-1, 0.5); a.set_ylim(0, 36)
    a.set_xlabel("Re β − 3/2"); a.set_ylabel("Im β"); a.legend(fontsize=8, loc="lower left")
    a.text(0.03, 34.5, "dotted: heights of ζ's zeros", fontsize=8, color="0.4")
    a.set_title("C1. Full chain at 90°: its copies of ζ's zeros wobble just\nright of the line, so no side is empty", fontsize=10)
    b = ax[1]; ns = sorted(int(k) for k in fe)
    b.bar([n-0.4 for n in ns], [fe[str(n)]["mirrored_total"] for n in ns], width=0.8, color="0.8", label="all zeros of the mirrored chain (box)")
    b.bar([n-0.4 for n in ns], [fe[str(n)]["mirrored_on_line"] for n in ns], width=0.8, color="tab:green", label="those on the line")
    b.plot(ns, [fe[str(n)]["Z_zeros_right_of_line"] for n in ns], "o-", color="tab:red", label="chain zeros right of the line (must be 0)")
    b.set_xlabel("chain length n"); b.legend(fontsize=8); b.set_xticks(ns)
    b.set_title("C2. Mirroring the full chain, Z(3/2+z) + Z(3/2−z):\nmost zeros stay off the line", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figC_far_end.png", dpi=150); plt.close(fig)

if __name__ == "__main__":
    for w in (sys.argv[1:] or ["A", "B", "C"]): globals()["fig" + w](); print("fig", w, "done", flush=True)
