"""DRAFT. Checks for ns-two-planes, from algebra_results.json, algebra3_results.json, runs/feed, twist_sectors.json and runs/*.json.
Run: python3 checks.py  -> checks_output.txt, checks_results.json (a few seconds; run algebra.py, batches first)."""
import json, os, numpy as np
L = lambda f: json.load(open(f))
res = []
def check(name, ok, note):
    res.append(dict(name=name, ok=bool(ok), note=note))
alg = L('algebra_results.json')
check('C1 exact algebra (algebra.py)', all(v['pass'] for v in alg.values()), f"{sum(v['pass'] for v in alg.values())} of {len(alg)} pass")
runs = lambda nu, m, k='1': [L(f'runs/n64_nu{nu}_{m}_k{k}_s{s}.json')['rows'] for s in range(3)]
flat = [r for nu in ['0.01', '0.005'] for m in ['A', 'B'] for r in runs(nu, m)]
grow = max(max(x['Z'] for x in r)/r[0]['Z'] for r in flat)
check('C2 a flat flow alone never grows its total spin', grow <= 1 + 1e-9 and all(r[-1]['off_frac'] < 1e-9 for r in flat),
      f'largest peak/start over 12 flat runs = {grow:.6f}; stays exactly flat')
two = {(nu, m): runs(nu, m) for nu in ['0.01', '0.005'] for m in ['twisted', 'untwisted']}
ratios = {k: [max(x['Z'] for x in r)/r[0]['Z'] for r in v] for k, v in two.items()}
check('C3 two crossing flat flows grow total spin in all 12 runs', all(x > 1.0 for v in ratios.values() for x in v),
      '; '.join(f"{k[0]} {k[1]}: {', '.join(f'{x:.2f}' for x in v)}" for k, v in ratios.items()))
offs = [r[-1]['off_frac'] for v in two.values() for r in v]
check('C4 most energy leaves the two planes by t = 8', min(offs) > 0.5, f'off-plane share at t=8: {min(offs):.2f} to {max(offs):.2f}')
def peak_to(f, T=4.0):
    r = [x for x in L(f)['rows'] if x['t'] <= T + 1e-9]; return max(x['Z'] for x in r)/r[0]['Z']
p9 = {m: [peak_to(f'runs/n64_nu0.005_{m}_k1_s{s}.json') for s in range(9)] for m in ['twisted', 'untwisted']}
diff = [a - b for a, b in zip(p9['twisted'], p9['untwisted'])]
check('C5 the twist does not tame the crossing (9 seeds, viscosity 0.005, peak by t = 4)',
      min(diff) < 0 < max(diff) and abs(np.mean(diff)) < 0.1*np.mean(p9['untwisted']),
      f'twisted lower in {sum(d < 0 for d in diff)} of 9; differences {", ".join(f"{d:+.2f}" for d in diff)}; '
      f'means {np.mean(p9["twisted"]):.2f} vs {np.mean(p9["untwisted"]):.2f}')
kd = {}
for m in ['twisted', 'untwisted']:
    kd[m] = []
    for k in ['0', '0.1', '0.25', '0.5', '1']:
        r = L(f'runs/n64_nu0.005_{m}_k{k}_s0.json')['rows']; kd[m].append((max(x['Z'] for x in r)/r[0]['Z'], r[-1]['off_frac']))
mono = all(all(b[0] >= a[0] - 1e-3 for a, b in zip(v, v[1:])) for v in kd.values())
check('C6 k setting on the crossing: peak spin rises with k (k = 0, 0.1, 0.25, 0.5, 1)', mono,
      '; '.join(f"{m}: {', '.join(f'{a:.2f}' for a, _ in v)}" for m, v in kd.items()))
k0 = [kd[m][0] for m in kd]
check('C7 k = 0: planes stay nearly flat and spin does not grow', all(z < 1.02 and o < 0.15 for z, o in k0),
      f"peak/start {k0[0][0]:.3f}, {k0[1][0]:.3f}; off-plane share {k0[0][1]:.3f}, {k0[1][1]:.3f} (the leak is the out-of-plane part fed through the shared line)")
a, b = L('runs/n64_nu0.005_twisted_k1_s0.json')['rows'], L('runs/n96_nu0.005_twisted_k1_s0.json')['rows']
za, zb = max(x['Z'] for x in a), max(x['Z'] for x in b)
ta = np.array([x['t'] for x in a]); tb = np.array([x['t'] for x in b])
gap = max(abs(np.interp(t, ta, [x['Z'] for x in a]) - z)/z for t, z in zip(tb, [x['Z'] for x in b]))
check('C8 grid check: 64³ and 96³ agree on total spin within 3% at every time', gap < 0.03,
      f'peak {za:.3f} vs {zb:.3f}; largest gap {100*gap:.1f}%; peak spin {max(x["wmax"] for x in a):.1f} vs {max(x["wmax"] for x in b):.1f}')
ts = L('twist_sectors.json')
sh = ts['random field: energy share by sector (1, i, -1, -i)']
cl = ts['nonlinear term of sector 0: share landing in sectors 0..3']
check('C9 twisted part: about a quarter of a random flow, closed under the nonlinear term', 0.2 < sh[0] < 0.35 and cl[0] > 1 - 1e-6,
      f'share {sh[0]:.3f}; sector -1 self-interaction lands in twisted part: {ts["nonlinear term of sector 2: share landing in sectors 0..3"][0]:.3f}')
if os.path.exists('runs/klein_planes.json'):
    kr = L('runs/klein_planes.json')['rows']
    fff = max(abs(x['P']['FFF']) for x in kr); tot = max(abs(sum(x['P'].values())) for x in kr)
    check('C10 Klein pair: the flat part alone does no stretching; all of it involves the bend', fff < 1e-6*max(tot, 1),
          f'max |P(flat,flat,flat)| = {fff:.1e} vs max total {tot:.2f}; flat share of spin {kr[0]["plane"]:.2f} -> {min(x["plane"] for x in kr):.2f}')

# ---- Round 2: three perpendicular planes (algebra3.py, threeplanes.py, batch4.sh) ----
if os.path.exists('algebra3_results.json'):
    a3 = L('algebra3_results.json')
    check('C11 three-plane algebra (algebra3.py): exact split, spins add up, all stretching on crossing lines, no cyclic twist on the cube',
          all(v['pass'] for v in a3.values()), f"{sum(v['pass'] for v in a3.values())} of {len(a3)} pass")
if os.path.exists('runs/p3_two-un_s0.json'):
    def lines_of(T):
        T = np.array(T); return np.array([T[1, :, 2].sum()+T[2, :, 1].sum(), T[0, :, 2].sum()+T[2, :, 0].sum(), T[0, :, 1].sum()+T[1, :, 0].sum()])
    P3 = {m: [L(f'runs/p3_{m}_s{s}.json')['rows'] for s in range(3)] for m in ['two-tw', 'two-un', 'three-tw', 'three-un']}
    def stats(r):
        t = np.array([x['t'] for x in r]); Z = np.array([x['Z'] for x in r]); ip = int(np.argmax(Z))
        comp = np.array([np.array(x['T']).sum(axis=(1, 2)) for x in r]); ln = np.array([lines_of(x['T']) for x in r])
        g = np.trapezoid(comp[:ip+1], t[:ip+1], axis=0); G = np.trapezoid(ln[:ip+1], t[:ip+1], axis=0)
        self_ = max(abs(np.array(x['T'])[i, :, i].sum()) for x in r for i in range(3)); tot = max(abs(np.array(x['T']).sum()) for x in r)
        return dict(peak=Z[ip]/Z[0], gain=g/g.sum(), line=G/G.sum(), early=ln[2]/ln[2].sum(), self_=self_/tot,
                    compgrow=np.array(r[ip]['comp_Z'])/np.maximum(np.array(r[0]['comp_Z']), 1e-12), endshare=np.array(r[-1]['comp_Z'])/r[-1]['Z'])
    S = {m: [stats(r) for r in v] for m, v in P3.items()}
    allS = [s for v in S.values() for s in v]
    check('C12 runs: no plane stretches its own spin; all stretching sits on the three crossing lines (12 runs)',
          max(s['self_'] for s in allS) < 1e-6, f"largest own-plane term / total = {max(s['self_'] for s in allS):.1e}")
    two = S['two-tw'] + S['two-un']
    check('C13 two-plane starts: the third plane (y, no spin at start) takes at least ~0.3 of the spin gain and ends with ~1/3 of the spin',
          all(s['gain'][1] > 0.3 and s['endshare'][1] > 0.3 for s in two),
          'gain share of plane y: ' + ', '.join(f"{s['gain'][1]:.2f}" for s in two) + '; spin share at t=4: ' + ', '.join(f"{s['endshare'][1]:.2f}" for s in two))
    check('C14 two-plane starts: stretching starts on the shared line y (>0.7 at t=0.4), then spreads to all three lines',
          all(s['early'][1] > 0.7 and min(s['line']) > 0.14 for s in two),
          'line y early: ' + ', '.join(f"{s['early'][1]:.2f}" for s in two) + '; smallest line share up to the peak: ' + ', '.join(f"{min(s['line']):.2f}" for s in two))
    pk = {m: [s['peak'] for s in v] for m, v in S.items()}
    check('C15 three planes grow more than two at the same energy (every three-plane seed above every two-plane seed)',
          min(pk['three-tw'] + pk['three-un']) > max(pk['two-tw'] + pk['two-un']),
          '; '.join(f"{m}: {', '.join(f'{x:.2f}' for x in v)}" for m, v in pk.items()))
    d3 = [a - b for a, b in zip(pk['three-tw'], pk['three-un'])]
    check('C16 the twist does not lower growth with three planes counted (3 seeds)', np.mean(d3) > -0.05 and sum(d < 0 for d in d3) < 3,
          f"twisted minus untwisted: {', '.join(f'{d:+.2f}' for d in d3)}; means {np.mean(pk['three-tw']):.2f} vs {np.mean(pk['three-un']):.2f}")
    three = S['three-tw'] + S['three-un']
    check('C17 three-plane starts: every plane\'s own spin grows (no plane stays quiet), each plane ~1/3 of the gain',
          all(min(s['compgrow']) > 1.5 and max(abs(s['gain'] - 1/3)) < 0.06 for s in three),
          'smallest per-plane growth at the peak: ' + ', '.join(f"{min(s['compgrow']):.1f}x" for s in three)
          + '; gain shares within ' + f"{max(max(abs(s['gain'] - 1/3)) for s in three):.2f} of 1/3")

# ---- Round 3: how fast can two planes feed each other? (feed.py, batch5.sh, batch6.sh) ----
if os.path.exists('runs/feed'):
    import glob, feed
    gc = {}
    for c in ['any', 'two', 'twod', 'one']:
        n, a = feed.gradcheck(feed.Prob(16, 0.1, c)); gc[c] = abs(n - a)/abs(n)
    check('C18 fastest-growth search: its gradient matches finite differences in every class', max(gc.values()) < 1e-6,
          ', '.join(f'{c} {v:.0e}' for c, v in gc.items()))
    FR = [L(f) for f in glob.glob('runs/feed/f_*.json')]
    B = {}
    for r in FR:
        k = (r['cls'], r['N'], r['nu'])
        if k not in B or r['R'] > B[k]['R']: B[k] = r
    offy = max(abs(np.array(r['T'])).sum() - abs(np.array(r['T'])[0, :, 2]).sum() - abs(np.array(r['T'])[2, :, 0]).sum() for r in FR if r['cls'] == 'two')
    check('C19 in a two-plane flow all stretching sits on the shared line (every two-plane optimum)', offy < 1e-9,
          f'largest stretching off the shared line = {offy:.1e}')
    hv = [(N, nu) for N, nu in [(32, 0.1), (32, 0.07), (64, 0.07)]]
    check('C20 at high viscosity no two-plane flow (shared line on a box axis) can grow its spin, while some 3D flow can',
          all(B[('two', N, nu)]['R'] < 0 < B[('any', N, nu)]['R'] for N, nu in hv),
          '; '.join(f"{N}³ visc {nu}: two {B[('two', N, nu)]['R']:+.3f}, any {B[('any', N, nu)]['R']:+.3f}" for N, nu in hv))
    rat = [(nu, B[('two', 64, nu)]['R']/B[('any', 64, nu)]['R']) for nu in [0.05, 0.035, 0.025]]
    check('C21 at lower viscosity two planes on a box axis come close to the fastest 3D growth, closer as viscosity drops (64³)',
          all(b >= a - 1e-3 for (_, a), (_, b) in zip(rat, rat[1:])) and rat[-1][1] > 0.9,
          ', '.join(f'visc {nu}: {q:.2f}' for nu, q in rat))
    lo = [(f"{r['N']}³ visc {r['nu']} seed {r['seed']}", r['M_eig'][0]/sum(r['M_eig'])) for r in FR
          if r['cls'] == 'any' and r.get('M_eig') and r['nu'] <= 0.05 and r['N'] >= 64]
    dg = [(N, nu, B[('twod', N, nu)]['R']/B[('any', N, nu)]['R']) for (c, N, nu) in B if c == 'twod']
    check('C22 the fastest-growing 3D flow has no spin along one fixed direction (a two-plane flow), and two planes sharing the box diagonal reach its growth within 0.2%',
          all(v < 1e-3 for _, v in lo) and all(abs(q - 1) < 2e-3 for *_, q in dg),
          'smallest spin-direction share: ' + ', '.join(f'{n} {v:.0e}' for n, v in lo) + '; diagonal ÷ any: ' + ', '.join(f'{N}³ visc {nu} {q:.4f}' for N, nu, q in sorted(dg)))
    if ('two', 96, 0.025) in B:
        q64 = B[('two', 64, 0.025)]['R']/B[('any', 64, 0.025)]['R']; q96 = B[('two', 96, 0.025)]['R']/B[('any', 96, 0.025)]['R']
        a64, a96 = B[('any', 64, 0.025)]['R'], B[('any', 96, 0.025)]['R']
        check('C23 grid check at viscosity 0.025: the two-plane ÷ any comparison holds at 96³ (the rate itself is not settled)', abs(q96 - q64) < 0.05,
              f'two ÷ any {q64:.3f} (64³) vs {q96:.3f} (96³); fastest rate itself {a64:.3f} vs {a96:.3f}, so 64³ is too coarse for the rate at this viscosity')
    ones = [r for r in FR if r['cls'] == 'one']
    check('C24 one layered flat flow alone can stretch (a strict one cannot, C2); best found growth is far below two planes',
          all(r['P'] > 0 for r in ones) and max(r['R'] for r in ones if r['nu'] == 0.025) < 0.2*B[('two', 64, 0.025)]['R'],
          f"stretching > 0 in all {len(ones)} one-plane optima; best found at visc 0.025: {max(r['R'] for r in ones if r['nu'] == 0.025):.3f} vs two planes {B[('two', 64, 0.025)]['R']:.3f} (a search result, not a bound)")
# ---- Round 4: a 3D view of a 4D flow (feedd.py, batch7.sh, drain.py) ----
if os.path.exists('runs/feedd'):
    import glob, feedd, feed, drain
    g3 = feedd.gradcheck(feedd.Prob(3, 16, 0.1)); g4 = feedd.gradcheck(feedd.Prob(4, 8, 0.1))
    e3, e4 = abs(g3[0] - g3[1])/abs(g3[0]), abs(g4[0] - g4[1])/abs(g4[0])
    p3 = feedd.Prob(3, 16, 0.05); uh = feedd.init(p3, 0, 1.0); va = p3.value(uh)[0]; vb = feed.Prob(16, 0.05, 'any').value(uh)[0]
    check('C25 any-dimension search (feedd.py): gradient matches finite differences in 3D and 4D, and in 3D it gives the same growth as feed.py on the same flow',
          max(e3, e4) < 1e-6 and abs(va - vb) < 1e-10, f'gradient error 3D {e3:.0e}, 4D {e4:.0e}; growth feedd {va:.10f} vs feed {vb:.10f}')
    u4 = np.load('runs/feedd/pic_d4_n16_nu0.05_s0_u.npz')['u'].astype(float); r4 = drain.rate4(u4, (16,)*4, 0.05); p4 = feedd.Prob(4, 16, 0.05); Rp = p4.value(p4.fft(u4))[0]
    check('C26 in 4D the search formula equals the growth computed straight from 4D Navier–Stokes (Leray-projected advection) on the 4D winner',
          abs(r4['total'] - Rp) < 1e-8*max(1, abs(Rp)), f"formula {Rp:.10f} vs direct {r4['total']:.10f} (saved field)")
    GR = [L(f) for f in glob.glob('runs/feedd/g_*.json')]; B = {}
    for r in GR:
        k = (r['d'], r['N'], r['nu'])
        if k not in B or r['R'] > B[k]['R']: B[k] = r
    s3 = [B[(3, N, 0.1)]['R'] for N in (16, 24, 32, 48)]; t3 = [B[(3, N, 0.05)]['R'] for N in (16, 24, 32, 48)]
    q4 = [B[(4, 24, nu)]['R']/B[(4, 16, nu)]['R'] for nu in (0.1, 0.05)]; q32 = B[(4, 32, 0.05)]['R']/B[(4, 24, 0.05)]['R']
    check('C27 the 3D worst case settles as the grid refines, while the 4D one keeps rising (at least 1.5x per refinement)',
          max(s3) - min(s3) < 1e-4 and abs(t3[-1] - t3[-2]) < 0.05*t3[-1] and min(q4 + [q32]) > 1.5,
          f"3D visc 0.1: {', '.join(f'{v:.5f}' for v in s3)}; 3D visc 0.05: {', '.join(f'{v:.3f}' for v in t3)}; 4D 16→24: x{q4[0]:.2f} (0.1), x{q4[1]:.2f} (0.05); 24→32 (visc 0.05, climb stopped at 150 steps, so a lower bound): x{q32:.2f}")
    w4 = [r for r in GR if r['d'] == 4]
    check('C28 every 4D winner is genuinely 4D: it varies along every direction (no direction carries over half its spin, and none is nearly constant)',
          all(max(r['flat_share']) < 0.5 and min(r['vary_eig']) > 0.05 for r in w4),
          f"largest flat share {max(max(r['flat_share']) for r in w4):.2f}; smallest share of variation along a direction {min(min(r['vary_eig']) for r in w4):.2f} ({len(w4)} runs)")
    Dr = L('runs/drain/drain_rings.json'); sc = Dr['scan']
    pred = [Dr['R3_view_alone']*(1 - q['frac_of_cancel']) for q in sc]
    lam = Dr['lambda_kw1.0']; c_pred = Dr['R3_view_alone']/lam
    check('C29 hidden 4th-direction motion lowers the 3D view\'s growth exactly in proportion to its spin, at the rate the eigenproblem gives',
          max(abs(q['view'] - p) for q, p in zip(sc, pred)) < 1e-8 and abs(Dr['hidden_spin_to_cancel'] - c_pred) < 1e-6,
          f"view growth {', '.join(f'{q['view']:+.4f}' for q in sc)}; largest gap from the straight line {max(abs(q['view'] - p) for q, p in zip(sc, pred)):.0e}; strongest drain λ = {lam:.3f}, hidden spin to cancel {Dr['hidden_spin_to_cancel']:.4f} = 0.2022/λ")
    tot = [q['total'] for q in sc]
    check('C30 the hidden motion does not limit anything: total 4D growth rises as the 3D view\'s growth falls',
          all(b > a for a, b in zip(tot, tot[1:])), f"total growth {', '.join(f'{v:.4f}' for v in tot)}")
# ---- Blow-up hunt: the push-back dial (hunt.py, hunt2.py, batch10.sh) ----
if os.path.exists('runs/hunt2'):
    import hunt_an as HA
    HR = HA.load_all('runs/hunt2')
    a = L('runs/hunt2/h_tg_n64_nu0_th1.json'); b = L('runs/hunt/h_tg_n64_nu0_th1.json')
    za = min(a['rows'], key=lambda r: abs(r['t'] - 2))['Z']; zb = min(b['rows'], key=lambda r: abs(r['t'] - 2))['Z']
    check('C31 the faster solver (hunt2.py, single precision) matches the first solver (hunt.py) on Taylor–Green, 64³',
          abs(HA.t_trust(a) - HA.t_trust(b)) < 1e-3 and abs(za - zb) < 1e-4*zb,
          f"how far the grid follows it: {HA.t_trust(a):.4f} vs {HA.t_trust(b):.4f}; spin at t = 2: {za:.6f} vs {zb:.6f}; time {a['secs']:.0f}s vs {b['secs']:.0f}s")
    eu = [r for k, r in HR.items() if k[2] == 0 and k[3] == 1.0]
    dE = max(max(abs(q['E'] - r['rows'][0]['E']) for q in r['rows'])/r['rows'][0]['E'] for r in eu)
    check('C32 with full push-back and no viscosity (Euler) energy is conserved in every run, while the grid follows the flow',
          dE < 1e-5, f'largest energy drift {dE:.1e} over {len(eu)} runs')
    g0 = [HA.grid_close(HR, 'tg', 0.0, c=c) for c in (1, 2)]
    check('C33 pressureless dust (no push-back): the three grids close in on the exact fold time 1 within 4%',
          all(abs(g['T'] - 1) < 0.04 and g['ratio'] < 0.8 for g in g0),
          '; '.join(f"threshold {c}: stops at {', '.join(f'{x:.3f}' for x in g['t'])} on 64³/128³/256³, step ratio {g['ratio']:.2f}, closes in on {g['T']:.3f}" for c, g in zip((1, 2), g0)))
    g5 = [HA.grid_close(HR, 'tg', 0.5, c=c) for c in (1, 2)]
    check('C34 half push-back: the grids still close in on one blow-up time (step ratio below 0.7)',
          all(g['ratio'] < 0.7 for g in g5), '; '.join(f"threshold {c}: step ratio {g['ratio']:.2f}, blow-up at {g['T']:.2f}" for c, g in zip((1, 2), g5)))
    fd = [L(f'runs/hunt2/h_rings_n{N}_nu0_th0.json') for N in (64, 128)]
    check('C35 colliding rings as dust: both grids stop just before the exact fold time 0.1105 (within 15%)',
          all(0.85*0.1105 < r['t_end'] <= 0.1105*1.02 for r in fd), f"stops at {fd[0]['t_end']:.3f} (64³), {fd[1]['t_end']:.3f} (128³)")
    ns = [HR[(ic, 128, 0.005, 1.0)] for ic in ('tg', 'rings', 'klein')]
    mn = [min(q['delta'] for q in r['rows'] if q.get('nfit', 99) >= 0.5*(128/3 - 3) and np.isfinite(q['delta']))*128/3 for r in ns]
    check('C36 Navier–Stokes (viscosity 0.005), 128³: all three candidates stay resolved to the end, peak spin rises and falls back',
          all(r['stop'] == '' for r in ns) and min(mn) > 2, '; '.join(f"{ic}: to t = {r['t_end']:.0f}, peak spin {max(q['wmax'] for q in r['rows']):.1f}, closest singularity {m:.1f} grid spacings"
                                                                      for ic, r, m in zip(('Taylor–Green', 'rings', 'tube pair'), ns, mn)))
    g1 = [HA.grid_close(HR, 'tg', 1.0, c=c) for c in (2, 4)]
    check('C37 real flow (full push-back), Taylor–Green, Euler: the grids do not close in on a blow-up time (step ratio at least 0.95)',
          all(g['ratio'] >= 0.95 for g in g1), '; '.join(f"threshold {c}: stops at {', '.join(f'{x:.2f}' for x in g['t'])} on 64³/128³/256³, step ratio {g['ratio']:.2f}" for c, g in zip((2, 4), g1)))
    TH = sorted(th for (i, n, nu, th) in HR if i == 'tg' and n == 256 and nu == 0)
    rs = [HA.grid_close(HR, 'tg', th, c=2)['ratio'] for th in TH]
    check('C38 the closing-in fades as push-back rises: on Taylor–Green the step ratio (threshold 2) rises from dust to real flow',
          rs[0] < 0.6 and rs[-1] > 0.95 and all(b > a - 0.05 for a, b in zip(rs, rs[1:])), ', '.join(f'θ {th:g}: {x:.2f}' for th, x in zip(TH, rs)))
    gk = HA.grid_close(HR, 'klein', 1.0, c=2); rk = HR[('klein', 256, 0.0, 1.0)]
    rows = [q for q in rk['rows'] if q.get('nfit', 0) >= 0.5*(256/3 - 3) and np.isfinite(q['delta']) and q['delta'] > 0]
    tk = np.array([q['t'] for q in rows]); dk = np.log([q['delta'] for q in rows]); m = len(tk)//3
    r0 = -(dk[m] - dk[0])/(tk[m] - tk[0]); r1 = -(dk[-1] - dk[-m])/(tk[-1] - tk[-m])
    check('C39 real flow, bent tube pair, Euler: 256³ follows it much further than 128³ (no closing in), and the singularity distance shrinks more slowly over time',
          gk['t'][2] - gk['t'][1] > 0.2 and r1 < 0.6*r0, f"threshold 2: stops at {', '.join(f'{x:.2f}' for x in gk['t'])} on 64³/128³/256³; δ shrinks at rate {r0:.2f} early and {r1:.2f} late on 256³")
lines = [f"{'PASS' if r['ok'] else 'FAIL'} {r['name']}: {r['note']}" for r in res]
lines.append(f"{sum(r['ok'] for r in res)} of {len(res)} pass")
open('checks_output.txt', 'w').write('\n'.join(lines) + '\n'); json.dump(res, open('checks_results.json', 'w'), indent=1)
print('\n'.join(lines))
