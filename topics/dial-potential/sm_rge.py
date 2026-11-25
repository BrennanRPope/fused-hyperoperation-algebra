"""Standard Model running of the Higgs quartic, for the "flat dial at the top scale" guesses.

Conventions (Buttazzo et al. 2013, arXiv:1307.3536):  V = lambda (H^dag H - v^2/2)^2,  m_h^2 = 2 lambda v^2,
gY hypercharge coupling (not GUT normalised).

Inputs at mu = M_t are the NNLO matching fits of Buttazzo et al. (their section 6), linear in the measured
pole masses and alpha_s.  Running: full two-loop gauge, one-loop plus the leading two-loop (top, QCD and
lambda) terms for y_t and lambda.  Electroweak two-loop terms in y_t and lambda are left out; they shift
lambda(M_Pl) by a few 10^-3, so this code is a cross-check.  The precise number used for scoring is
Buttazzo's own NNLO fit formula for lambda(M_Pl), quoted in `lambda_mpl_buttazzo`.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp

K = 1.0 / (16 * math.pi ** 2)
MPL_FULL = 1.22089e19
MPL_RED = 2.43536e18


def inputs_at_mt(Mh, Mt, as_mz, MW=80.3692):
    """MS-bar couplings at mu = Mt (Buttazzo et al. 2013, eqs. 58-62, reference Mt=173.34, Mh=125.15)."""
    a3 = (as_mz - 0.1184) / 0.0007
    lam = 0.12604 + 0.00206 * (Mh - 125.15) - 0.00004 * (Mt - 173.34)
    yt = 0.93690 + 0.00556 * (Mt - 173.34) - 0.00042 * a3
    g3 = 1.1666 + 0.00314 * a3
    g2 = 0.64779 + 0.00004 * (Mt - 173.34) + 0.00011 * (MW - 80.384) / 0.014
    gY = 0.35830 + 0.00011 * (Mt - 173.34) - 0.00020 * (MW - 80.384) / 0.014
    return np.array([gY, g2, g3, yt, lam])


def betas(t, c, two_loop=True):
    gY, g2, g3, yt, lam = c
    Y2, W2, G2, T2 = gY * gY, g2 * g2, g3 * g3, yt * yt
    # one loop
    bY = 41 / 6 * gY ** 3
    bW = -19 / 6 * g2 ** 3
    bG = -7 * g3 ** 3
    bt = yt * (4.5 * T2 - 8 * G2 - 2.25 * W2 - 17 / 12 * Y2)
    bl = (24 * lam ** 2 - 6 * T2 ** 2 + 3 / 8 * (2 * W2 ** 2 + (W2 + Y2) ** 2)
          + lam * (12 * T2 - 9 * W2 - 3 * Y2))
    out = K * np.array([bY, bW, bG, bt, bl])
    if two_loop:
        # two-loop gauge (gY normalisation): rows (Y, W, G), columns (Y, W, G), top column
        bY2 = gY ** 3 * (199 / 18 * Y2 + 9 / 2 * W2 + 44 / 3 * G2 - 17 / 6 * T2)
        bW2 = g2 ** 3 * (3 / 2 * Y2 + 35 / 6 * W2 + 12 * G2 - 3 / 2 * T2)
        bG2 = g3 ** 3 * (11 / 6 * Y2 + 9 / 2 * W2 - 26 * G2 - 2 * T2)
        bt2 = yt * (-12 * T2 ** 2 + 6 * lam ** 2 - 12 * lam * T2 + 36 * G2 * T2 - 108 * G2 ** 2)
        bl2 = (-312 * lam ** 3 - 144 * lam ** 2 * T2 - 3 * lam * T2 ** 2 + 30 * T2 ** 3
               + 80 * lam * G2 * T2 - 32 * G2 * T2 ** 2)
        out = out + K * K * np.array([bY2, bW2, bG2, bt2, bl2])
    return out


def run(Mh, Mt, as_mz, mu_end, two_loop=True, n=400):
    c0 = inputs_at_mt(Mh, Mt, as_mz)
    t_end = math.log(mu_end / Mt)
    ts = np.linspace(0, t_end, n)
    sol = solve_ivp(betas, (0, t_end), c0, t_eval=ts, rtol=1e-10, atol=1e-12, args=(two_loop,))
    return Mt * np.exp(sol.t), sol.y


def lambda_at(Mh, Mt, as_mz, mu, two_loop=True):
    mus, ys = run(Mh, Mt, as_mz, mu, two_loop, n=2)
    return ys[4, -1], betas(0, ys[:, -1], two_loop)[4]


def lambda_mpl_buttazzo(Mh, Mt, as_mz):
    """Buttazzo et al. 2013, NNLO: lambda at M_Pl = 1.22e19 GeV, linear fit around the reference point."""
    return (-0.0143 - 0.0066 * (Mt - 173.34) + 0.0018 * (as_mz - 0.1184) / 0.0007
            + 0.0029 * (Mh - 125.15))


def mh_for_flat(Mt, as_mz, mu=MPL_FULL, use_buttazzo=True, offset=0.0):
    """Higgs mass that makes lambda(mu) = 0."""
    if use_buttazzo and mu == MPL_FULL:
        return 125.15 + (0.0143 + 0.0066 * (Mt - 173.34) - 0.0018 * (as_mz - 0.1184) / 0.0007) / 0.0029
    from scipy.optimize import brentq
    return brentq(lambda m: lambda_at(m, Mt, as_mz, mu)[0] + offset, 100, 160)


if __name__ == "__main__":
    # calibration: Buttazzo reference point gives lambda(M_Pl) = -0.0143 at NNLO
    l1, b1 = lambda_at(125.15, 173.34, 0.1184, MPL_FULL, two_loop=False)
    l2, b2 = lambda_at(125.15, 173.34, 0.1184, MPL_FULL, two_loop=True)
    print(f"reference point: lambda(MPl) one-loop {l1:.4f}, our two-loop {l2:.4f}, Buttazzo NNLO -0.0143")
    print(f"beta_lambda(MPl): one-loop {b1:.5f}, two-loop {b2:.5f}")
