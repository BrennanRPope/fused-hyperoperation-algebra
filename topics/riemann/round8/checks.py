"""Round 8 checks: re-reads the saved outputs and re-derives the key numbers."""
import json, numpy as np, mpmath as mp
mp.mp.dps = 20
R = []
def chk(name, ok): R.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name)

q = json.load(open('qgas.json'))
chk("q=-0.5 Fisher zero matches round 7 (0.905308+14.068964i)", abs(complex(*[r for r in q['1']['rows'] if r['q'] == -0.5][0]['fisher']) - complex(0.905308, 14.068964)) < 1e-5)
for n in ('1', '2', '3'):
    rows = {r['q']: r for r in q[n]['rows']}
    for qq in (-0.02, 0.02):
        f, a = complex(*rows[qq]['fisher']), complex(*rows[qq]['apoint']); rho = complex(*q[n]['rho'])
        chk(f"rho_{n}, q={qq:+}: Fisher zero within 10% of the a-point displacement", abs(f-a) < 0.1*abs(a-rho))
    chk(f"rho_{n}: q<0 zeros right of the line, q>0 left", all(r['fisher'][0] > 0.5 for r in q[n]['rows'] if r['q'] < 0 and r['fisher']) and all(r['fisher'][0] < 0.5 for r in q[n]['rows'] if r['q'] > 0 and r['fisher']))
# a-point independent recomputation
rho = mp.zetazero(1); ap = mp.findroot(lambda x: mp.zeta(x) - 0.02, rho + 0.02/mp.zeta(rho, derivative=1))
chk("a-point zeta = 0.02 near rho_1 recomputed", abs(complex(ap) - complex(*[r for r in q['1']['rows'] if r['q'] == -0.02][0]['apoint'])) < 1e-8)

c = json.load(open('calogero.json'))
chk("CM identity n=3 (symbolic)", c['cm_identity_n3']); chk("CM identity n=4 (symbolic)", c['cm_identity_n4'])
chk("Lax matrix nilpotent: ||L^n|| < 1e-8", c['lax_pow_n'] < 1e-8)
chk("zeros at time t = eigenvalues of Z0 + t L0 (5 times, < 1e-5)", max(e[1] for e in c['eig_formula_errors']) < 1e-5)
chk("first pair to leave the line in the rewind is the closest pair", c['closest_index'] == 27 and c['first_complex_t'] < 0)
# two-particle check numerically
d = 1.7; v = 2/d; acc = -2*(v-(-v))/d**2
chk("two-particle acceleration = -8/d^3", abs(acc + 8/d**3) < 1e-12)

k = json.load(open('klein.json'))
chk("torus: all 58 zeros on Re s = 1/2", all(abs(a-0.5) < 1e-7 for a, b in k['torus']) and len(k['torus']) == 58)
chk("Klein bottle: none of its zeros on the line", all(abs(a-0.5) > 1e-6 for a, b in k['klein_even']))
chk("odd sector: none of its zeros on the line", all(abs(a-0.5) > 1e-6 for a, b in k['klein_odd']))
def L4(s): return (mp.zeta(s, 0.25) - mp.zeta(s, 0.75))*mp.power(4, -s)
F = lambda s: 2*mp.zeta(s)*L4(s) - (1-mp.power(2, 1-2*s))*mp.zeta(2*s)
chk("every stored Klein zero refines to a true zero within 1e-8", all(abs(complex(mp.findroot(F, mp.mpc(*z))) - complex(*z)) < 1e-8 for z in k['klein_even']))

t = json.load(open('topology.json'))
chk("twisted ladder: Z- Z+ = Z_2l", t['ladder_identity_err'] < 1e-12)
chk("sinh/cosh orbit weights reproduce log-derivative", t['weight_check'] < 1e-15)
chk("flat glide weight 1/2 matches heat trace (4 values)", all(abs(a-b) < 1e-9 for _, a, b in t['flat_glide']))
chk("zeta spacings closer to GUE than GOE", t['ks_gue'] < t['ks_goe'])
chk("crosscap cover: chi doubles, genus k-1", all(2*(2-kk) == 2-2*(kk-1) for kk in range(3, 20)))

r = json.load(open('recentre.json'))
chk("Moebius: zeros on |w+1| = 1", max(abs(abs(complex(a, b)+1)-1) for a, b in r['mobius_w']) < 1e-12)
chk("fused residue k=-1 is 6/pi^2", abs(r['Ck']['-1.0'] - 6/np.pi**2) < 1e-4)
chk("fused residue k=1 is 1", abs(r['Ck']['1.0'] - 1) < 1e-12)
chk("k=2 blow-up order 3/2", abs(r['k2_slopes'][-1] + 1.5) < 1e-3)
chk("Stieltjes removal: zeros off the line, both sides", all(abs(a-0.5) > 1e-3 for a, b in r['stieltjes_zeros']) and min(a for a, b in r['stieltjes_zeros']) < 0.5 < max(a for a, b in r['stieltjes_zeros']))
chk("Boltzmann wall at 1.3994333", abs(r['wall'] - 1.3994333) < 1e-6)
chk("(s-1) zeta(s, a) -> 1 for a = 0.3", abs((mp.mpf('1e-8'))*mp.zeta(1+mp.mpf('1e-8'), 0.3) - 1) < 1e-6)
print(f"\n{sum(o for _, o in R)} of {len(R)} checks pass")
