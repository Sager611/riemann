# trunk.py -- jerbyconf: Jerby's trunk-slope law (2610.21952 l. 1826-1851): |Z'(g_n)| ~ A_act/(g_{n+1}-g_n), A_act = max(|Z(g_{n-1})|, |Z(g_{n+1})|).
# Ratio rho = |Z'(g_n)| (g_{n+1}-g_n) / A_act with g_{n+1}-g_n = pi/theta'(g_n) + O(1/t), on the low-viscosity isolated bad points of the scans
# (CAND lines, complete for mu/theta' < 0.62).  'coherent' (l. 1824-1826: one neighbour small, the other much larger): min/max neighbour < 0.2.
import re, math, collections
F = '../JerbyFals/logs/'
rows = []
for fn in [F + x for x in ['run_L1.log', 'run_L2a.log', 'run_L2b.log', 'reso_R1.log', 'reso_R2a.log', 'reso_R2b.log', 'reso_R3.log']] + ['logs/run_M.log']:
    for line in open(fn):
        if not line.startswith('CAND'): continue
        d = dict(re.findall(r'(\w+)=([-\d.e+]+)', line)); zs = [float(x) for x in re.search(r'Z\[-2..2\]=([^ ]+)', line).group(1).split(',')]
        t = float(d['t']); mu = float(d['mu']); nm = float(d['norm']); thp = mu/nm
        if nm >= 0.62: continue
        A = max(abs(zs[1]), abs(zs[3])); rho = abs(float(d['Zp']))*(math.pi/thp)/A
        coh = min(abs(zs[1]), abs(zs[3]))/A < 0.2
        typ = {0: 'zero-free', 2: 'confined'}.get(int(d['zeros_left']) + int(d['zeros_right']), 'other')
        rows.append(dict(n=int(d['n']), t=t, mu=mu, rho=rho, coh=coh, typ=typ, Z=abs(float(d['Z'])), A=A))
uniq = {r['n']: r for r in rows}.values()
bins = collections.OrderedDict()
for r in sorted(uniq, key=lambda r: r['t']):
    bins.setdefault(math.floor(2*math.log10(r['t']))/2, []).append(r)
print('| log10 t bin | points (mu/theta\' < 0.62) | 40 lowest mu: rho min / median / max | coherent among them | coherent: rho median | |Z(g_n)|/A_act median |')
print('|---|---|---|---|---|---|')
for k, rs in bins.items():
    lo = sorted(rs, key=lambda r: r['mu'])[:40]; rr = sorted(r['rho'] for r in lo); co = sorted(r['rho'] for r in lo if r['coh'])
    za = sorted(r['Z']/r['A'] for r in lo)
    print(f"| {k:.1f} | {len(rs)} | {rr[0]:.2f} / {rr[len(rr)//2]:.2f} / {rr[-1]:.1f} | {len(co)} | " + (f"{co[len(co)//2]:.2f}" if co else '–') + f" | {za[len(za)//2]:.2f} |")
for lab, ns in [('jerbyfals counterexamples', [1565734083, 2565007238, 3315629581])]:
    print(lab, ', '.join(f"n={n}: rho={[r for r in uniq if r['n'] == n][0]['rho']:.2f}" for n in ns))
allr = sorted(r['rho'] for r in uniq); print(f'all {len(allr)} points: rho min {allr[0]:.3f}, median {allr[len(allr)//2]:.2f}, max {allr[-1]:.0f}')
