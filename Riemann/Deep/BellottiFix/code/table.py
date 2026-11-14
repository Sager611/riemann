
import json, math
o = json.load(open('logs/orig.json')); f = json.load(open('logs/fix.json')); ind = json.load(open('logs/indep_all.json'))
def up(c):   # Bellotti's convention, directed upward: 2 decimals below 10^6, else 3 significant figures
    if c < 1e6: return '%.2f' % (math.ceil(c * 100 - 1e-9) / 100)
    e = int(math.floor(math.log10(c))); m = math.ceil(c / 10**e * 100 - 1e-9) / 100
    if m >= 10: m /= 10; e += 1
    return '%.2f·10^%d' % (m, e)
mx_code = max(abs(float(ind['%d_code' % r['row']]['C']) / r['code']['Constant C'] - 1) for r in o)
mx_fix = max(abs(float(ind['%d_flaws2' % r['row']]['C']) / x['code']['Constant C'] - 1) for r, x in zip(o, f))
print('indep(code) vs his code: max rel %.2e; indep(flaws2) vs corrected code: max rel %.2e' % (mx_code, mx_fix))
share = [float(ind['%d_code' % r['row']]['c3']) / sum(float(ind['%d_code' % r['row']][k]) for k in ('c2', 'c3', 'c4')) for r in o]
print('c3 share of c2+c3+c4: %.3f-%.3f' % (min(share), max(share)))
lT = lambda s: math.log(3e12) if s.startswith('3') else float(s[5:-1])
lw = [float(r['uxvw'][3]) * lT(r['T0']) for r in o]
print('log W range %.2f-%.2f' % (min(lw), max(lw)))
print('F4need max %.4f; H<=Hb all: %s; min cond4 all-mode %.4f' % (max(float(ind['%d_all' % r['row']]['F4need']) for r in o),
      all(float(ind['%d_all' % r['row']]['H']) <= float(ind['%d_all' % r['row']]['Hb']) for r in o), min(float(ind['%d_all' % r['row']]['cond4']) for r in o)))
lines = ['| # | T0 | T1 | alpha0 | B | C printed | C his code | C both flaws fixed | ratio | C all fixes | ratio |', '|---|---|---|---|---|---|---|---|---|---|---|']
rf = []; ra = []
for r in o:
    n = r['row']; P = r['paper']['Constant C']; cf = float(ind['%d_flaws2' % n]['C']); ca = float(ind['%d_all' % n]['C'])
    rf.append(cf / P); ra.append(ca / P)
    t0 = r['T0'].replace('\\cdot 10^{12}', '·10^12').replace('\\exp', 'e^').replace('(', '').replace(')', '')
    t1 = r['T1'].replace('\\exp', 'e^').replace('(', '').replace(')', '')
    lines.append('| %d | %s | %s | %s | %.3f | %s | %s | %s | %.3f | %s | %.3f |' % (n, t0, t1, r['alpha0'], r['paper']['Constant B'],
                 up(P) if P < 1e6 else up(P), up(r['code']['Constant C']), up(cf), cf / P, up(ca), ca / P))
open('logs/table.md', 'w').write('\n'.join(lines) + '\n')
print('ratio both-flaws/printed %.3f-%.3f; all/printed %.3f-%.3f' % (min(rf), max(rf), min(ra), max(ra)))
print('\n'.join(lines[:4] + lines[-2:]))
