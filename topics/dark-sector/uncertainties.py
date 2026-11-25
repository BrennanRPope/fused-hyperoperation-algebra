"""1-sigma errors from the curvature of chi2 (numerical Hessian) for the one-parameter LCDM fit and the free layer weights."""
import json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import expansion_fits as X
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "expansion_results.json")))
def joint(model, p): return X.sn_chi2(model, p) + X.bao_chi2(model, p)[0]
def hess(f, p, h=2e-3):
    p = np.array(p, float); n = len(p); H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            e_i = np.eye(n)[i] * h; e_j = np.eye(n)[j] * h
            H[i, j] = (f(p + e_i + e_j) - f(p + e_i - e_j) - f(p - e_i + e_j) + f(p - e_i - e_j)) / (4 * h * h)
    return H
out = {}
for key, which in (("Flat LCDM (cosmological constant)", "sn"), ("Flat LCDM (cosmological constant)", "sn+bao"),
                   ("Mult-add cube, free weight per layer", "sn+bao"), ("Flat wCDM (dark energy with w free)", "sn+bao")):
    model = X.MODELS[key][0]; p = R["fits"][key][which]["params"]
    f = (lambda q: X.sn_chi2(model, q)) if which == "sn" else (lambda q: joint(model, q))
    C = np.linalg.inv(hess(f, p) / 2)
    out[f"{key} [{which}]"] = dict(params=p, sigma=np.sqrt(np.diag(C)).tolist(), cov=C.tolist())
    print(key, which, np.round(p, 3), "+/-", np.round(np.sqrt(np.diag(C)), 3))
# test equal layer weights (LCDM) inside the cube model
key = "Mult-add cube, free weight per layer"; p = np.array(out[f"{key} [sn+bao]"]["params"]); C = np.array(out[f"{key} [sn+bao]"]["cov"])
d = np.array([p[0] - p[1], p[1] - p[2]]); J = np.array([[1, -1, 0], [0, 1, -1]]); Cd = J @ C @ J.T
out["equal_weights_test"] = dict(diffs=d.tolist(), chi2=float(d @ np.linalg.inv(Cd) @ d), dof=2)
print("equal-weights test chi2 =", out["equal_weights_test"]["chi2"], "for 2 dof")
json.dump(out, open(os.path.join(HERE, "uncertainties.json"), "w"), indent=1)
