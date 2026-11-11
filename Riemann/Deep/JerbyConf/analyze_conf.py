# analyze_conf.py -- jerbyconf: confined-pair statistics per window (Jerby 2610.21952 l. 653-663 liminf law, l. 2039 lower scale)
# Inputs: JerbyFals window logs (CAND lines complete for mu<5 or mu/theta'<0.62), JerbyConf run_M.log (21 windows rerun, complete for
# mu/theta'<0.80), resonance logs R1-R3 (CAND complete for mu<5 or mu/theta'<0.62 inside the checked regions).
# Confined := zeros_left + zeros_right == 2 (one interval holds both; parity forces it); zero-free := 0.  s*mu = pi*mu/theta'.
import re, math, sys, collections
F = '../JerbyFals/logs/'
def parse(fn):
    C, W = [], {}
    for line in open(fn):
        if line.startswith('CAND'):
            d = dict(re.findall(r'(\w+)=([-\d.e+]+)', line)); w = line.split()[1]
            zs = [float(x) for x in re.search(r'Z\[-2..2\]=([^ ]+)', line).group(1).split(',')]
            C.append(dict(w=w, n=int(d['n']), t=float(d['t']), mu=float(d['mu']), nm=float(d['norm']), Z=float(d['Z']), Zp=float(d['Zp']),
                          Zs=zs, zl=int(d['zeros_left']), zr=int(d['zeros_right'])))
        elif line.startswith('WIN'):
            d = dict(re.findall(r'(\w+)=([-\d.e+\[\],]+)', line)); w = line.split()[1]
            t0, t1 = [float(x) for x in re.search(r't=\[([^\]]+)\]', line).group(1).split(',')]
            W[w] = dict(t0=t0, t1=t1, cnt=int(d['cnt']), iso=int(d['iso']), thp=float(d['thetap']))
        elif line.startswith('RESO'):
            d = dict(re.findall(r'(\w+)=([-\d.e+]+)', line)); W[line.split()[1]] = dict(T=float(d['T']), H=float(d['H']), iso=int(d['iso']), gp=int(d['gram_points']))
    return C, W
CL, WL = [], {}
for fn in ['run_L1.log', 'run_L2a.log', 'run_L2b.log']:
    c, w = parse(F + fn); CL += c; WL.update(w)
CM, WM = parse('logs/run_M.log')
rerun = set(WM)
rows = []
for j in range(64):
    w = 'L%02d' % j
    src = 'M' if w in rerun else 'L'
    cands = [x for x in (CM if src == 'M' else CL) if x['w'] == w]
    thr = 0.80 if src == 'M' else 0.62
    conf = [x for x in cands if x['zl'] + x['zr'] == 2]
    zf = [x for x in cands if x['zl'] + x['zr'] == 0]
    oth = [x for x in cands if x['zl'] + x['zr'] not in (0, 2)]
    # counts below fixed normalised levels (complete for nm < 0.62 in every window)
    lev = [1.80, 1.90, 1.948]
    cnt = {L: sum(1 for x in conf if math.pi*x['nm'] < L) for L in lev}
    cntzf = {L: sum(1 for x in zf if math.pi*x['nm'] < L) for L in lev}
    m = min(conf, key=lambda x: x['nm']) if conf else None
    mm = min(conf, key=lambda x: x['mu']) if conf else None
    rows.append(dict(w=w, src=src, thr=thr, t=WL[w]['t0'], iso=WL[w]['iso'], cnt=WL[w]['cnt'], nconf=len(conf), nzf=len(zf), noth=len(oth),
                     min_smu=(math.pi*m['nm'] if m else None), at=m, min_mu_conf=(mm['mu'] if mm else None), cnt_lev=cnt, cntzf_lev=cntzf,
                     min_smu_zf=(math.pi*min(x['nm'] for x in zf) if zf else None)))
print('| win | t (start) | Gram pts | isolated bad | complete below s*mu | confined: min s*mu | its mu | its n | confined: min mu | # conf s*mu<1.948 | # zero-free s*mu<1.948 | zero-free: min s*mu |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    a = r['at']
    print(f"| {r['w']} | {r['t']:.4e} | {r['cnt']:,} | {r['iso']:,} | {math.pi*r['thr']:.3f} | " +
          (f"{r['min_smu']:.4f} | {a['mu']:.4f} | {a['n']} | {r['min_mu_conf']:.4f} |" if a else f"> {math.pi*r['thr']:.3f} | – | – | – |") +
          f" {r['cnt_lev'][1.948]} | {r['cntzf_lev'][1.948]} | " + (f"{r['min_smu_zf']:.4f} |" if r['min_smu_zf'] else '– |'))
# trend: per-window min s*mu (windows where determined) vs log10 t
xs = [math.log10(r['t']) for r in rows if r['min_smu']]; ys = [r['min_smu'] for r in rows if r['min_smu']]
def fit(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys)); syy = sum((y-my)**2 for y in ys)
    b = sxy/sxx; return my - b*mx, b, sxy/math.sqrt(sxx*syy)
a, b, rr = fit(xs, ys)
print(f'\nTREND per-window min s*mu over confined ({len(xs)} windows with a confined point below the completeness level): s*mu_min = {a:.3f} + {b:.4f} log10 t (r = {rr:.2f}); max of minima {max(ys):.4f}, min {min(ys):.4f}')
xs2 = [math.log10(r['t']) for r in rows]; ys2 = [r['min_mu_conf'] for r in rows]
# binned fractions per isolated bad point
print('\n| height bin (log10 t) | windows | isolated bad | confined s*mu<1.80 | <1.90 | <1.948 | per 10^5 iso (<1.948) | zero-free s*mu<1.948 | per 10^5 |')
print('|---|---|---|---|---|---|---|---|---|')
bins = collections.OrderedDict()
for r in rows:
    k = math.floor(2*math.log10(r['t']))/2
    bins.setdefault(k, []).append(r)
for k, rs in bins.items():
    iso = sum(r['iso'] for r in rs); c = {L: sum(r['cnt_lev'][L] for r in rs) for L in (1.80, 1.90, 1.948)}; z = sum(r['cntzf_lev'][1.948] for r in rs)
    print(f"| {k:.1f} | {len(rs)} | {iso:,} | {c[1.80]} | {c[1.90]} | {c[1.948]} | {1e5*c[1.948]/iso:.2f} | {z} | {1e5*z/iso:.2f} |")
# resonance scans
print('\nRESONANCE scans (confined = 2 zeros in the block):')
for fn in ['reso_R1.log', 'reso_R2a.log', 'reso_R2b.log', 'reso_R3.log']:
    c, w = parse(F + fn)
    conf = [x for x in c if x['zl'] + x['zr'] == 2]; zf = [x for x in c if x['zl'] + x['zr'] == 0]
    m = min(conf, key=lambda x: x['nm']); mm = min(conf, key=lambda x: x['mu'])
    print(f"{fn}: iso {list(w.values())[0]['iso']:,}; confined CAND {len(conf)}, min s*mu {math.pi*m['nm']:.4f} (n={m['n']}, t={m['t']:.4e}, mu={m['mu']:.4f}); min mu {mm['mu']:.4f} (n={mm['n']}, s*mu={math.pi*mm['nm']:.4f}); # conf s*mu<1.948: {sum(1 for x in conf if math.pi*x['nm']<1.948)}; zero-free CAND {len(zf)}")
# list of lowest confined s*mu overall
allc = [x for x in CL + CM if x['zl'] + x['zr'] == 2]
for fn in ['reso_R1.log', 'reso_R2a.log', 'reso_R2b.log', 'reso_R3.log']:
    allc += [x for x in parse(F + fn)[0] if x['zl'] + x['zr'] == 2]
seen = set(); low = []
for x in sorted(allc, key=lambda x: x['nm']):
    if x['n'] in seen: continue
    seen.add(x['n']); low.append(x)
print('\nLOWEST confined s*mu (all samples):')
for x in low[:15]:
    print(f"  {x['w']} n={x['n']} t={x['t']:.6f} mu={x['mu']:.6f} s*mu={math.pi*x['nm']:.5f} Z={x['Z']:.4f} Zp={x['Zp']:.4f} Z[-1,+1]={x['Zs'][1]:.4f},{x['Zs'][3]:.4f} zl={x['zl']} zr={x['zr']}")
lowmu = sorted({x['n']: x for x in allc}.values(), key=lambda x: x['mu'])[:8]
print('LOWEST confined mu (all samples):')
for x in lowmu:
    print(f"  {x['w']} n={x['n']} t={x['t']:.6f} mu={x['mu']:.6f} s*mu={math.pi*x['nm']:.5f} zl={x['zl']} zr={x['zr']}")
