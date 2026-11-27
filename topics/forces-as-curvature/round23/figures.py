"""DRAFT picture for round 23. Run after checks.py. Blue = the law's own geometry, orange = what the mass adds, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()
q = 1.0; MN = 1 / math.sqrt(27 * q)
def hz(M):
    rr = np.roots([-q, 0, 1, -2 * M]); rr = np.sort(rr[np.abs(rr.imag) < 1e-12].real); return rr[rr > 0]

# A. f(r) for growing mass
a = ax[0]; r = np.linspace(0.01, 1.25, 800)
for frac, col, lw in ((0, BL, 2.6), (0.3, OR, 1.6), (0.9, OR, 1.6), (1.0, OR, 2.2), (1.5, OR, 1.0)):
    M = frac * MN; f = 1 - 2 * M / r - q * r**2
    lab = "M = 0: the law's de Sitter, 1 − q r²" if frac == 0 else f"M = {frac:g} × Nariai mass"
    a.plot(r, f, color=col, lw=lw, ls="-" if frac <= 1 else ":", label=lab)
    for h in (hz(M) if frac < 1 else []):
        a.plot(h, 0, "o", color=col, ms=6)
a.plot(1 / math.sqrt(3 * q), 0, "s", color=OR, ms=8, mfc="none")
a.annotate("Nariai: the two horizons merge\nat r = 1/√(3q)", (1 / math.sqrt(3 * q), 0), (0.62, -0.98), fontsize=9, color=OR, arrowprops=dict(arrowstyle="->", color=OR))
a.annotate("the law's cosmic horizon\nr = 1/√q", (1, 0), (0.95, 0.5), fontsize=9, color=BL, arrowprops=dict(arrowstyle="->", color=BL))
a.axhline(0, color=GR, lw=0.8); a.set_ylim(-1.2, 1.05); a.set_xlim(0, 1.25)
a.set_xlabel("areal radius r (units of 1/√q)"); a.set_ylabel("f(r) = 1 − 2M/r − q r²   (horizons where f = 0)")
a.legend(fontsize=8.5, loc="lower left"); a.set_title("A. A point mass in the law's de Sitter: black-hole and cosmic horizons (E1–E3)")

# B. embedding of the slice at rest: a tighter law sphere glued smoothly into the outside
b = ax[1]
cs = R["dust_ball"]["cases"][0]; a0, chi0, M = cs["a0"], cs["chi0"], cs["M"]
chi = np.linspace(0, math.pi / 2, 300); b.plot(np.sin(chi), -np.cos(chi), color=BL, lw=1.4, ls="--", label="no matter: the law's 3-sphere at q (equatorial slice)")
ch = np.linspace(0, chi0, 200); Rin = a0 * np.sin(ch); zin = -a0 * np.cos(ch)
b.plot(Rin, zin, color=BL, lw=2.8, label=f"inside the ball: the law's 3-sphere at a larger dial, q_in = {cs['q_inside']:.2f}")
Re = a0 * math.sin(chi0); rc = hz(M)[1]
u = np.linspace(math.sqrt(rc - Re), 0, 4000); Rs = rc - u**2                # R = rc - u^2 removes the 1/sqrt(f) end point
fR = 1 - 2 * M / Rs - q * Rs**2; fp = abs(2 * M / rc**2 - 2 * q * rc)
with np.errstate(all="ignore"):
    g = np.where(u > 1e-6, 2 * u * np.sqrt(np.clip((1 - fR) / fR, 0, None)), 2 * np.sqrt((1 - fR) / fp))
z = zin[-1] + np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * -np.diff(u))])
b.plot(Rs, z, color=OR, lw=2.4, label=f"outside: the slice of Schwarzschild–de Sitter, M = {M:.3f}")
b.plot(Re, zin[-1], "o", color="k", ms=6); b.annotate(f"edge of the ball (angle χ0 = {chi0}):\nslopes match because f = cos²χ0 there", (Re, zin[-1]), (0.36, -0.75), fontsize=9, arrowprops=dict(arrowstyle="->"))
b.plot(rc, z[-1], "s", color=OR, ms=7); b.text(rc - 0.03, z[-1] - 0.06, "cosmic horizon", ha="right", color=OR, fontsize=9)
b.set_aspect("equal"); b.set_xlim(0, 1.15); b.set_ylim(-1.1, max(0.3, z[-1] + 0.35)); b.set_xlabel("areal radius R"); b.set_ylabel("height in the embedding (only shape matters)")
b.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 1.0)); b.set_title("B. The ball's slice at rest, drawn as a surface (E4)")

# C. sectional curvatures on that slice: inside equal (the law), outside split by the tide
c = ax[2]
rr = np.linspace(0.02, rc, 600); inside = rr < Re
across = np.where(inside, cs["q_inside"], q + 2 * M / rr**3); along = np.where(inside, cs["q_inside"], q - M / rr**3)
c.plot(rr, across, color=OR, lw=2.2, label="outside, bending across the sphere: q + 2M/r³")
c.plot(rr, along, color=OR, lw=2.2, ls="--", label="outside, bending along the radius: q − M/r³")
c.plot(rr, (across + 2 * along) / 3, color=BL, lw=2.2, label="average (scalar curvature / 6): q outside, q_in inside")
c.plot(rr[inside], across[inside], color=BL, lw=4, alpha=0.5, label="inside: every plane bends alike (the law at q_in)")
c.axhline(q, color=GR, lw=0.8, ls=":"); c.axvline(Re, color="k", lw=0.8, ls=":"); c.text(Re + 0.01, 5.6, "edge", fontsize=9)
c.set_ylim(-2.5, 6.2); c.set_xlabel("areal radius r"); c.set_ylabel("sectional curvature of the slice at rest")
c.legend(fontsize=8.5, loc="lower right"); c.set_title("C. Where the law stops: the tidal split M/r³ outside the lump (E5)")

# D. the edge falls freely
d = ax[3]
for cs_, col in zip(R["dust_ball"]["cases"], (OR, BL)):
    s = np.array(cs_["s"]); Rv = np.array(cs_["R"])
    d.plot(s, Rv, "o", color=col, ms=4, label=f"{cs_['case']}: the edge, inside (Friedmann) = Kottler geodesic to {cs_['match']:.0e}")
    for h, ls in zip(cs_["horizons"], (":", "--")):
        d.axhline(h, color=col, lw=0.9, ls=ls)
d.text(1.15, R["dust_ball"]["cases"][0]["horizons"][0] + 0.03, "black-hole horizon (dense ball)", color=OR, fontsize=8.5)
d.text(0.05, R["dust_ball"]["cases"][1]["horizons"][1] + 0.03, "cosmic horizons (both balls)", color=GR, fontsize=8.5)
d.set_ylim(0, 2.5); d.set_xlim(-0.05, 2.5); d.set_xlabel("proper time of the ball's dust"); d.set_ylabel("areal radius of the ball's edge")
d.legend(fontsize=8.5, loc="upper left"); d.set_title("D. A dense ball collapses to a black hole; a light one is carried off")

fig.suptitle("DRAFT round 23: a lump of matter in the law's de Sitter spacetime", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("figures/f30_matter_curves_the_law.png", bbox_inches="tight")
