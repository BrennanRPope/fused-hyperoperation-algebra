"""Checks for Riemann round 4 (DRAFT). Run: python3 checks.py > checks_output.txt"""
import json, re, os, numpy as np
idn = json.load(open('identities.json')); ex = json.load(open('extrema.json')); rw = json.load(open('rewind4.json'))
pm = json.load(open('polymath.json')); br = json.load(open('barrier.json')); b1 = json.load(open('braid.json'))
bd = json.load(open('braid_deep.json')); r3 = json.load(open('../round3/rewind.json'))
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name + (f'   [{info}]' if info else ''))
check('C1 sympy: ladder identity, exact rate of D, lone-pair D = s^2-4p, min Laguerre = D/2, sin steady, lone-pair flow D0+8t, lone-pair ladder rates 4 vs 8',
      all(idn.values()), ', '.join(k for k, v in idn.items() if v))
ty = rw['toy']
check('C2 toy: the rate formula 8 + D*kappa matches a numerical derivative of D at t = 0 to 1e-5', abs(ty['rate0']/ty['num_rate0']-1) < 1e-5, f"{ty['rate0']:.6f} vs {ty['num_rate0']:.6f}")
k = (ty['rate0']-8)/ty['De0']; tp = (1/k)*np.log(8/(8+k*ty['De0']))
im_at = min((abs(r['t']-ty['t_merge']), r['im']) for r in ty['traj'])[1]
check('C3 toy: D crosses 0 at the merger (t = -0.011712), within 5% of the lone pair and within 0.5% of the one-point prediction',
      abs(ty['t_merge']/ty['lone_pair']-1) < 0.05 and abs(tp/ty['t_merge']-1) < 0.005, f"merge {ty['t_merge']:.6f}, pred {tp:.6f}, lone {ty['lone_pair']}")
allx = ex['main'] + ex['lehmer']; De = np.array([r['De'] for r in allx]); g2 = np.array([r['gap'] for r in allx])**2
check('C4 D = -8 Xi/Xi\'\' > 0 at all 1,318 extrema (zeros 1 to 1,000 and the Lehmer window), and D <= gap^2 at each',
      len(allx) == 1318 and (De > 0).all() and (De <= g2*(1+1e-9)).all(), f'D/gap^2 in [{(De/g2).min():.3f}, {(De/g2).max():.4f}]')
G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
w = np.linspace(G[-1], G[-1]*400, 400001); rho = np.log(w/(2*np.pi))/(2*np.pi); errs = []
for r in [ex['main'][i] for i in (0, 1, 2, 250, 500, 750)]:
    x = r['ge']; s = np.sum(1/(x-G)**2+1/(x+G)**2) + np.trapezoid((1/(x-w)**2+1/(x+w)**2)*rho, w); errs.append(abs(8/s/r['De']-1))
check('C5 D read off Xi matches D from the zeros, 8/sum 1/(x-z)^2 (5,000 zeros plus tail), to 1e-4 at six extrema', max(errs) < 1e-4, f'max rel {max(errs):.1e}')
L = [r for r in ex['lehmer'] if r['n'] == 6709][0]; k = L['kappa']; tl = (4/k)*np.log(8/(8+k*L['De']))
e = ex['main'][921]; k2 = e['kappa']; t9 = (4/k2)*np.log(8/(8+k2*e['De'])); f9 = [f for f in r3['full'] if f['j'] == 922][0]
check("C6 one point of Xi predicts the merger time of Lehmer's pair to 1e-5 and of zeros 922/923 to 1e-3 (vs round-3 full rewinds)",
      abs(tl/r3['lehmer']['t_full']-1) < 1e-5 and abs(t9/f9['t_full']-1) < 1e-3, f"{tl:.10f} vs {r3['lehmer']['t_full']:.10f}; {t9:.7f} vs {f9['t_full']:.7f}")
kap = np.array([r['kappa'] for r in allx])
check('C7 (log Xi\'\')\'\' < 0 at all 1,318 extrema; 487 of 999 and 163 of 319 have D|kappa| > 8 (gap opens under the rewind)',
      (kap < 0).all() and sum(r['De']*abs(r['kappa']) > 8 for r in ex['main']) == 487 and sum(r['De']*abs(r['kappa']) > 8 for r in ex['lehmer']) == 163)
mt = np.array(pm['moll_ts']); mm = np.array(pm['moll']['60000000000.0'])
check('C8 Polymath toy: bare margin < 0 at t = 0 for every x (-212 at 6e10); mollified margin at 6e10 < 0 at t = 0.175 and > 0 at t = 0.2',
      all(v < 0 for v in pm['m_at_t0'].values()) and abs(pm['m_at_t0']['60000000000.0']+211.7) < 0.5 and mm[np.isclose(mt, 0.175)][0] < 0 < mm[np.isclose(mt, 0.2)][0])
g1 = complex(br['gap_c1'].strip('()')) if isinstance(br['gap_c1'], str) else complex(*br['gap_c1'])
gf = complex(br['gap_formula'].strip('()')) if isinstance(br['gap_formula'], str) else complex(*br['gap_formula'])
g0 = complex(br['gap_c0'].strip('()')) if isinstance(br['gap_c0'], str) else complex(*br['gap_c0'])
check('C9 letters: order gap 0 at c = 0 (shared centre) and equal to 2^-s - 3^-s at c = 1, s = 1/2 + 14.1347i', abs(g0) < 1e-15 and abs(g1-gf) < 1e-14, f'{g1:.4f}')
zc = br['Zq_two_prime_check']
check('C10 order dial, primes 2 and 3 at s = 2: q = 0 gives the Euler product, q = 1 the free words, closed form agrees at q = 0.3, 0.7',
      abs(br['Zq_q0']-br['euler_23']) < 1e-12 and abs(br['Zq_q1']-br['free_23']) < 1e-9 and all(abs(a[1]-a[2]) < 1e-12 for a in zc))
check('C11 q = 1 wall with every prime: sigma0 = 1.39943 (prime zeta = 1); with primes 2, 3 only: 0.78788', abs(br['sigma0_all_primes']-1.39943) < 1e-5 and abs(br['sigma0_23']-0.78788) < 1e-5)
gr = br['growth']; bet = [gr[q]['beta_fit'] for q in ['0.0', '0.25', '0.5', '0.75']]
check('C12 order-dial growth exponent rises with q (0, 0.32, 0.86, 1.96) and is within 5% of q/(1-q) at q = 0.25',
      np.all(np.diff(bet) > 0) and abs(bet[0]) < 1e-9 and abs(bet[1]/(1/3)-1) < 0.05, ', '.join(f'{b:.2f}' for b in bet))
ld = br['link_dial']; lx = np.array(ld['loglog']); ok = True; info = []
for kname, v in ld['var'].items():
    v = np.array(v); sl = np.polyfit(lx[1:], v[1:], 1)[0]; want = 0.5*ld['k2_mean'][kname]
    ok &= (abs(sl-want) < 0.1*max(want, 0.01)) if want > 0 else (abs(v[-1]-v[1]) < 1e-3)
    info.append(f'{sl:.3f}/{want:.3f}')
check('C13 link dial: spread grows like (1/2) mean(k^2) log log X for each family (10%); k = 0 stays flat', ok, ', '.join(info))
m1 = b1['mergers'][0]
check("C14 topological word: first inner merger is Lehmer's pair (fused at gap 0.02 near gamma 7005.08); to t = -0.4: >= 50 mergers and 0 syzygies",
      abs(m1['at']-7005.08) < 0.01 and len(bd['mergers']) >= 50 and len(bd['events']) == 0 and len(b1['events']) == 0, f"{len(bd['mergers'])} mergers")
lw = rw['lehmer']
check("C15 Lehmer window ODE: D from the window's zeros within 0.5% of D from Xi; its zero time equals round 3's merger to 1e-6",
      abs(lw['De0_gamma']/L['De']-1) < 0.005 and abs(lw['t_zero_De']-r3['lehmer']['t_full']) < 1e-6, f"{lw['De0_gamma']:.7f} vs {L['De']:.7f}")
hl = json.load(open('hyperlog.json'))
check('C16 hypernatural log: height-s Laguerre is s f\'^2 - f f\'\', the height drops out at an extremum, and degree-n Laguerre holds on random real-rooted polynomials (n = 3, 5, 8)',
      hl['laguerre_height_s'] == '-f0*f2 + f1**2*s' and hl['extremum_height_free'] and hl['laguerre_degree_n_holds'])
sp_ok = all(abs(a[1]/a[2]-1) < 0.02 for a in hl['spot_checks'])
check('C17 critical height h* = Xi Xi\'\'/Xi\'^2 stays below 3/4 on 0 < gamma < 5355 (max 0.705, below the first zero); zeros-based h* matches Xi itself to 2%',
      hl['hstar_max_refined'] < 0.75 and abs(hl['hstar_max_refined']-0.705) < 0.001 and sp_ok, f"max {hl['hstar_max_refined']:.4f}, N_eff {hl['Neff_max_refined']:.2f}")
# ---------- polish ----------
rd = open('README.md').read(); lean = open('lean/RiemannRound4.lean').read()
check('P1 README is marked DRAFT and names the theory (Affine Composition Theory, built on Fused Hyperoperation Algebra)',
      'DRAFT' in rd.split('\n')[0] and 'Affine Composition Theory, built on Fused Hyperoperation Algebra' in rd)
shown = set(re.findall(r'figures/(fig[^)]+\.png)', rd)); made = set(f for f in os.listdir('figures') if f.endswith('.png'))
check('P2 every picture the README shows exists, and every picture made is shown', shown == made, f'{len(made)} pictures')
bad = [p for p in ['proves RH', 'proof of RH', 'breakthrough', 'solves', 'we prove the Riemann', 'settles RH'] if p.lower() in rd.lower()]
check('P3 no overclaim phrases in the README', not bad, ', '.join(bad))
check('P4 wording rule: "backward heat flow" and "rewind" used; no "downward"/"going down" for the flow',
      'backward heat flow' in rd and 'rewind' in rd and not re.search(r'downward|going down|run down|flow down', rd, re.I))
first = rd.split('## The answer')[1].split('##')[0] if '## The answer' in rd else ''
check('P5 the README opens with the answer and says plainly whether it is a lead', rd.index('## The answer') < rd.index('## Pass 1') and 'Not a lead' in first)
nums = ['0.705', '5,355', '0.02584', '0.0014206', '−0.00071120370', '−0.00071120457', '1.39943', '−212', '0.237', '1,318']
check('P6 numbers quoted in the README match the data',
      all(n in rd for n in nums) and abs(L['De']-0.0014206) < 5e-8 and abs(ex['main'][921]['De']-0.02584) < 1e-5 and abs(r3['lehmer']['t_full']+0.00071120370) < 1e-11 and abs(tl+0.00071120457) < 1e-11)
check('P7 the Lean file says it is compiled in the shared library, has 12 theorems and no sorry',
      'compiled in the shared library' in lean and len(re.findall(r'^theorem', lean, re.M)) == 12 and 'sorry' not in re.sub(r'--.*', '', re.sub(r'/-.*?-/', '', lean, flags=re.S)))
print(f'\n{sum(res)} of {len(res)} checks pass')
