"""Riemann's Phi(u) as a weight over pendulum swings: u = 1/2 ln(K'/K), m = sin^2(theta0/2)."""
import mpmath as mp, numpy as np, json
mp.mp.dps = 30
def Phi(u):
    u = mp.mpf(u); return sum((2*mp.pi**2*n**4*mp.e**(9*u) - 3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u)) for n in range(1, 12))
th = np.linspace(0.5, 179.5, 359)
u = []; w = []
for t in th:
    m = mp.sin(mp.radians(t)/2)**2; K = mp.ellipk(m); Kp = mp.ellipk(1-m)
    uu = 0.5*mp.log(Kp/K)
    # du/dtheta0 = (1/2) (dx/dm)/x * dm/dtheta0 ; dx/dm = -pi/(4 m (1-m) K^2); dm/dtheta = sin(theta0)/2
    dudt = 0.5*(-mp.pi/(4*m*(1-m)*K**2))/(Kp/K)*mp.sin(mp.radians(t))/2*mp.pi/180
    u.append(float(uu)); w.append(float(Phi(uu)*abs(dudt)))
w = np.array(w); u = np.array(u)
print("u at 90deg:", u[np.argmin(abs(th-90))], " weight peak at", th[np.argmax(w)], "deg; half-max band", th[w > w.max()/2].min(), th[w > w.max()/2].max())
print("mirror check: weight(90-d) vs weight(90+d)", w[np.argmin(abs(th-60))], w[np.argmin(abs(th-120))])
json.dump(dict(theta=th.tolist(), u=u.tolist(), w=w.tolist()), open("swingdensity.json", "w"))
