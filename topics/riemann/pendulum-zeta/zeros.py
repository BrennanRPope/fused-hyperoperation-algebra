"""Zeros of the completed Epstein zeta Lam(s)=(2pi/sqrt y)^-s Gamma(s) Z(s,y) for 0.5<=Im s<=T:
total in the strip by argument principle, on the line by sign changes."""
import mpmath as mp, sys, json
from epstein import Z
mp.mp.dps = 20
y = mp.mpf(sys.argv[1]); T = float(sys.argv[2]) if len(sys.argv) > 2 else 30
def Lam(s): return (mp.sqrt(y)/mp.pi)**s*mp.gamma(s)*Z(s, y, K=int(T/(2*3.14159*float(mp.sqrt(y))))+8)
# on-line: Lam(1/2+it) is real
ts = [0.5 + i*0.05 for i in range(int((T-0.5)/0.05)+1)]
vals = [mp.re(Lam(mp.mpc(0.5, t))) for t in ts]
online = [(ts[i]+ts[i+1])/2 for i in range(len(ts)-1) if vals[i]*vals[i+1] < 0]
# argument principle on rectangle [-1,2] x [0.5,T]
def seg(a, b, n):
    return [a + (b-a)*j/n for j in range(n)]
path = seg(mp.mpc(2, 0.5), mp.mpc(2, T), int(T*6)) + seg(mp.mpc(2, T), mp.mpc(-1, T), 60) + seg(mp.mpc(-1, T), mp.mpc(-1, 0.5), int(T*6)) + seg(mp.mpc(-1, 0.5), mp.mpc(2, 0.5), 60)
path.append(path[0])
ph = [Lam(s) for s in path]
wind = 0.0
for i in range(len(ph)-1):
    wind += float(mp.im(mp.log(ph[i+1]/ph[i])))
total = round(wind/(2*3.141592653589793))
out = dict(y=float(y), T=T, total=total, online=len(online), online_t=online, offline=total-len(online))
print(json.dumps(out))
json.dump(out, open(f"zeros_y{float(y):g}.json", "w"))
